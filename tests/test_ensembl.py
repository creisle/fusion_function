import pytest

from fusion_function.ensembl import get_protein_domains

pytestmark = pytest.mark.integration

DOMAIN_CASES = [
    # gene, transcript, expected feature/domain text
    ("TP53", "ENST00000269305", "DNA-binding"),
    ("BRCA1", "ENST00000357654", "BRCT"),
    ("EGFR", "ENST00000275493", "protein kinase"),
    ("BRAF", "ENST00000646891", "protein kinase"),
    ("AMACR", "ENST00000335606", "racemase"),
]


def _feature_names(result: dict) -> list[str]:
    """Collect names exposed by Ensembl plus PANTHER subfamily enrichment."""
    names: list[str] = []
    for feature in result["protein_features"]:
        if feature["description"]:
            names.append(feature["description"])
        if feature["panther_subfamily_description"]:
            names.append(feature["panther_subfamily_description"])
    return names


@pytest.mark.parametrize(
    "gene,transcript_id,expected_domain", DOMAIN_CASES, ids=[case[0] for case in DOMAIN_CASES]
)
def test_get_protein_domains_contains_expected_domain(
    gene: str, transcript_id: str, expected_domain: str, reference_db
):
    result = get_protein_domains(transcript_id, reference=reference_db)
    assert "error" not in result, f"{gene} ({transcript_id}) returned an error: {result}"
    assert result["protein_features"], f"{gene} ({transcript_id}) returned no protein features"
    names = _feature_names(result)
    assert names, f"{gene} ({transcript_id}) returned protein features but no names"
    assert any(expected_domain.casefold() in name.casefold() for name in names), (
        f"{gene} ({transcript_id}) did not contain expected domain "
        f"{expected_domain!r}.\n"
        "Returned feature names:\n" + "\n".join(f"  - {name}" for name in names)
    )


def test_amacr_panther_subfamily_is_alpha_methylacyl_coa_racemase(reference_db):
    result = get_protein_domains("ENST00000335606", reference=reference_db)
    assert "error" not in result
    panther_features = [
        feature
        for feature in result["protein_features"]
        if (feature["source"] or "").casefold() == "panther"
    ]
    assert panther_features, "AMACR returned no PANTHER protein feature"
    assert any(
        feature["panther_subfamily_id"] == "PTHR48228:SF5" for feature in panther_features
    ), (
        "AMACR did not map to expected PANTHER subfamily PTHR48228:SF5.\n"
        f"Returned PANTHER features: {panther_features!r}"
    )
    assert any(
        "racemase" in (feature["panther_subfamily_description"] or "").casefold()
        for feature in panther_features
    ), (
        "AMACR PANTHER subfamily description did not contain 'racemase'.\n"
        f"Returned PANTHER features: {panther_features!r}"
    )
