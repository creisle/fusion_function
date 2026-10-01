"""Check indexed CDS mapping independently across split codons."""

import pytest

from fusion_function import data


@pytest.mark.parametrize("strand", [1, -1])
@pytest.mark.parametrize("start,end", [(1, 1), (23, 24), (34, 58), (68, 68), (1, 68)])
def test_indexed_feature_mapping_across_split_codons(strand, start, end):
    blocks = []
    cds_position = 1
    for offset, length in enumerate([70, 101, 33]):
        genomic_start = 100 + 300 * (offset if strand == 1 else 2 - offset)
        blocks.append(
            {
                "cds_start": cds_position,
                "cds_end": cds_position + length - 1,
                "genomic_start": genomic_start,
                "genomic_end": genomic_start + length - 1,
                "chromosome": "1",
                "strand": strand,
                "assembly_name": "GRCh38",
            }
        )
        cds_position += length
    segments = data.protein_segments(
        start, end, blocks, block_ends=[block["cds_end"] for block in blocks]
    )
    reference_bases = []
    for block in blocks:
        positions = list(range(block["genomic_start"], block["genomic_end"] + 1))
        reference_bases.extend(positions if strand == 1 else reversed(positions))
    mapped_bases = []
    for segment in segments:
        positions = list(range(segment["start"], segment["end"] + 1))
        mapped_bases.extend(positions if strand == 1 else reversed(positions))
    assert mapped_bases == reference_bases[(start - 1) * 3 : end * 3]
    assert segments == data.protein_segments(start, end, blocks)
