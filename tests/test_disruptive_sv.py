"""Single-gene SVs remain valid; same-terminus pairs are unsupported."""

import pytest

from fusion_function import annotate_fusion_domains


SLC45A2 = ("ENST00000296589", "5:33973126", "N")
AMACR = ("ENST00000335606", "5:34006836", "C")


def utr_breakpoint(transcript, utr):
    """Choose a base immediately before or after the complete CDS on either strand."""
    if utr == "5prime":
        block = min(transcript["cds_blocks"], key=lambda b: b["cds_start"])
        assert block["premrna_start"] > 1, "Test transcript needs a 5' UTR"
        position = (
            block["genomic_start"] - 1 if transcript["strand"] == 1 else block["genomic_end"] + 1
        )
    else:
        block = max(transcript["cds_blocks"], key=lambda b: b["cds_end"])
        assert block["premrna_end"] < transcript["premrna_length"], "Test transcript needs a 3' UTR"
        position = (
            block["genomic_end"] + 1 if transcript["strand"] == 1 else block["genomic_start"] - 1
        )
    return f"{transcript['chromosome']}:{position}"


def coding_breakpoint(transcript, cds_position):
    """Locate a chosen original CDS base using the annotated exon blocks."""
    block = next(
        b for b in transcript["cds_blocks"] if b["cds_start"] <= cds_position <= b["cds_end"]
    )
    offset = cds_position - block["cds_start"]
    position = (
        block["genomic_start"] + offset
        if transcript["strand"] == 1
        else block["genomic_end"] - offset
    )
    return f"{transcript['chromosome']}:{position}"


def single_partner_result(reference, transcript_id, breakpoint, terminus):
    """Explicitly supply None for the unknown transcript, as the batch caller does."""
    slot = 1 if terminus == "N" else 2
    missing_slot = 2 if slot == 1 else 1
    result = annotate_fusion_domains(
        **{
            f"transcript{slot}_id": transcript_id,
            f"breakpoint{slot}": breakpoint,
            f"gene{slot}_terminus": terminus,
            f"transcript{missing_slot}_id": None,
            f"breakpoint{missing_slot}": None,
            f"gene{missing_slot}_terminus": None,
            "reference": reference,
        }
    )
    assert "error" not in result, result
    assert result["assumed_disruptive"] is True
    assert result["domains"]
    assert {d["transcript_id"] for d in result["domains"]} == {transcript_id}
    for domain in result["domains"]:
        if domain["breakpoint_based_status"] == "excluded":
            assert domain["post_splicing_status"] is None
            assert domain["post_translation_status"] is None
        else:
            assert domain["post_splicing_status"] is not None
            if "preserved" in domain["post_splicing_status"].split("/"):
                assert domain["post_translation_status"] is not None
    return result


@pytest.mark.integration
@pytest.mark.parametrize(
    "transcript_id", ["ENST00000305877", "ENST00000296589"], ids=["BCR-plus", "SLC45A2-minus"]
)
@pytest.mark.parametrize(
    "terminus,utr,status,percent",
    [
        pytest.param("C", "5prime", "included", 100.0, id="C-5UTR"),
        pytest.param("C", "3prime", "excluded", 0.0, id="C-3UTR"),
        pytest.param("N", "5prime", "excluded", 0.0, id="N-5UTR"),
        pytest.param("N", "3prime", "included", 100.0, id="N-3UTR"),
    ],
)
def test_single_gene_utr_breakpoints(reference_db, transcript_id, terminus, utr, status, percent):
    """UTR cuts either retain all original protein features or exclude them all."""
    transcript = reference_db.get_transcript(transcript_id)
    assert "error" not in transcript, transcript
    result = single_partner_result(
        reference_db, transcript_id, utr_breakpoint(transcript, utr), terminus
    )
    assert all(d["breakpoint_based_status"] == status for d in result["domains"])
    assert all(d["breakpoint_retained_percent"] == percent for d in result["domains"])
    if status == "included":
        assert result["frame_status"] == "in_frame"
        assert "alternative_start_used" not in result
        assert all(d["post_splicing_status"] == "preserved" for d in result["domains"])
        assert all(d["post_translation_status"] == "preserved" for d in result["domains"])
    else:
        assert result["frame_status"] is None


@pytest.mark.integration
@pytest.mark.parametrize(
    "transcript_id", ["ENST00000305877", "ENST00000296589"], ids=["BCR-plus", "SLC45A2-minus"]
)
def test_c_terminal_start_codon_loss_retains_downstream_domains(reference_db, transcript_id):
    transcript = reference_db.get_transcript(transcript_id)
    assert "error" not in transcript, transcript
    first = min(transcript["cds_blocks"], key=lambda b: b["cds_start"])
    assert first["cds_start"] == 1
    assert first["premrna_end"] - first["premrna_start"] >= 3
    assert (
        transcript["premrna_sequence"][first["premrna_start"] - 1 : first["premrna_start"] + 2]
        == "ATG"
    )
    # Retaining from CDS base 4 removes the entire native start codon. A domain
    # beginning at residue 2 or later can nevertheless remain completely intact
    # at the breakpoint level. Translation is checked from the next retained ATG.
    result = single_partner_result(
        reference_db, transcript_id, coding_breakpoint(transcript, 4), "C"
    )
    assert "alternative_start_found" in result["translation_start"][transcript_id].split("/")
    downstream = [d for d in result["domains"] if d["start"] >= 2]
    assert downstream
    assert all(d["breakpoint_based_status"] == "included" for d in downstream)
    assert all(d["breakpoint_retained_percent"] == 100.0 for d in downstream)


@pytest.mark.integration
@pytest.mark.parametrize(
    "transcript_id", ["ENST00000305877", "ENST00000296589"], ids=["BCR-plus", "SLC45A2-minus"]
)
def test_c_terminal_breakpoint_cuts_domain(reference_db, transcript_id):
    transcript = reference_db.get_transcript(transcript_id)
    assert "error" not in transcript, transcript
    full = single_partner_result(
        reference_db, transcript_id, utr_breakpoint(transcript, "5prime"), "C"
    )
    target = max(full["domains"], key=lambda d: d["end"] - d["start"])
    assert target["end"] > target["start"] + 1
    cds_position = 3 * (target["start"] - 1) + 1 + 3 * (target["end"] - target["start"] + 1) // 2
    result = single_partner_result(
        reference_db, transcript_id, coding_breakpoint(transcript, cds_position), "C"
    )
    (cut,) = [
        d
        for d in result["domains"]
        if (d["interpro_id"], d["name"], d["start"], d["end"])
        == (target["interpro_id"], target["name"], target["start"], target["end"])
    ]
    assert "disrupted" in cut["breakpoint_based_status"].split("/")
    assert 0.0 < cut["breakpoint_retained_percent"] < 100.0


@pytest.mark.integration
@pytest.mark.parametrize(
    "transcript_id", ["ENST00000305877", "ENST00000296589"], ids=["BCR-plus", "SLC45A2-minus"]
)
def test_late_n_terminal_breakpoint_retains_all_domains(reference_db, transcript_id):
    """The disruptive event assumption must not override complete domain retention."""
    transcript = reference_db.get_transcript(transcript_id)
    assert "error" not in transcript, transcript
    last_cds = max(transcript["cds_blocks"], key=lambda block: block["cds_end"])
    # The last CDS base is a late breakpoint on either genomic strand. It retains
    # the entire original protein feature span, despite the unknown partner.
    position = last_cds["genomic_end"] if transcript["strand"] == 1 else last_cds["genomic_start"]
    result = annotate_fusion_domains(
        transcript1_id=transcript_id,
        breakpoint1=f"{transcript['chromosome']}:{position}",
        gene1_terminus="N",
        reference=reference_db,
    )
    assert "error" not in result, result
    assert result["assumed_disruptive"] is True
    assert result["domains"]
    assert all(domain["breakpoint_based_status"] == "included" for domain in result["domains"])
    assert all(domain["breakpoint_retained_percent"] == 100.0 for domain in result["domains"])
    assert result["frame_status"] == "in_frame"
    assert all(d["post_splicing_status"] == "preserved" for d in result["domains"])
    assert all(d["post_translation_status"] == "preserved" for d in result["domains"])


@pytest.mark.integration
@pytest.mark.parametrize("known_slot", [1, 2])
@pytest.mark.parametrize("partner", [SLC45A2, AMACR], ids=["SLC45A2", "AMACR"])
def test_single_known_partner_is_assumed_disruptive(reference_db, known_slot, partner):
    transcript, breakpoint, terminus = partner
    result = annotate_fusion_domains(
        **{
            f"transcript{known_slot}_id": transcript,
            f"breakpoint{known_slot}": breakpoint,
            f"gene{known_slot}_terminus": terminus,
            "reference": reference_db,
        }
    )
    assert "error" not in result, result
    assert result["assumed_disruptive"] is True
    assert result["frame_status"] in ("in_frame", "out_of_frame", "in_frame/out_of_frame")
    assert result["domains"]
    assert {d["transcript_id"] for d in result["domains"]} == {transcript}
    assert any(d["breakpoint_based_status"] == "disrupted" for d in result["domains"])
    assert any(0 < d["breakpoint_retained_percent"] < 100 for d in result["domains"])
    retained = [d for d in result["domains"] if d["breakpoint_based_status"] != "excluded"]
    assert all(d["post_splicing_status"] is not None for d in retained)
    assert any(d["post_translation_status"] is not None for d in retained)


@pytest.mark.parametrize("terminus", ["N", "C"])
@pytest.mark.parametrize("same_transcript", [False, True])
@pytest.mark.parametrize("with_breakpoints", [False, True])
def test_same_terminus_pair_raises_before_database_lookup(
    monkeypatch, terminus, same_transcript, with_breakpoints
):
    def unexpected_lookup(*, reference=None, database=None, release=None):
        pytest.fail("Unsupported pair must be rejected before loading a reference")

    monkeypatch.setattr("fusion_function.fusion.get_reference", unexpected_lookup)
    with pytest.raises(NotImplementedError, match="N/N and C/C"):
        annotate_fusion_domains(
            SLC45A2[0],
            SLC45A2[0] if same_transcript else AMACR[0],
            SLC45A2[1] if with_breakpoints else None,
            AMACR[1] if with_breakpoints else None,
            terminus,
            terminus,
        )


@pytest.mark.integration
def test_c_n_order_matches_n_c_product(reference_db):
    forward = annotate_fusion_domains(
        SLC45A2[0], AMACR[0], SLC45A2[1], AMACR[1], "N", "C", reference=reference_db
    )
    reverse = annotate_fusion_domains(
        AMACR[0], SLC45A2[0], AMACR[1], SLC45A2[1], "C", "N", reference=reference_db
    )
    assert "error" not in forward, forward
    assert reverse == forward
    assert "assumed_disruptive" not in reverse
    assert reverse["frame_status"] is not None


@pytest.mark.parametrize(
    "arguments",
    [
        {},
        {"transcript1_id": "", "transcript2_id": None},
        {"transcript1_id": SLC45A2[0], "breakpoint1": SLC45A2[1]},
        {"transcript2_id": AMACR[0], "gene2_terminus": "C"},
        {"transcript1_id": SLC45A2[0], "breakpoint1": SLC45A2[1], "gene1_terminus": "n"},
    ],
)
def test_missing_required_input_returns_error_before_database_lookup(monkeypatch, arguments):
    def unexpected_lookup(*, reference=None, database=None, release=None):
        pytest.fail("Incomplete event should be rejected before loading a reference")

    monkeypatch.setattr("fusion_function.fusion.get_reference", unexpected_lookup)
    result = annotate_fusion_domains(**arguments)
    assert set(result) == {"error"}
    assert result["error"]


@pytest.mark.parametrize("breakpoint", ["1:33973126", "5:1", "5:33973126-33973127"])
def test_invalid_single_partner_breakpoint_returns_error(reference_db, breakpoint):
    result = annotate_fusion_domains(
        transcript1_id=SLC45A2[0],
        breakpoint1=breakpoint,
        gene1_terminus="N",
        reference=reference_db,
    )
    assert set(result) == {"error"}
    assert SLC45A2[0] in result["error"]


def test_single_partner_does_not_hide_transcript_lookup_failure(reference_db):
    result = annotate_fusion_domains(
        transcript1_id="ENST00000000000",
        breakpoint1="5:33973126",
        gene1_terminus="N",
        reference=reference_db,
    )
    assert set(result) == {"error"}
