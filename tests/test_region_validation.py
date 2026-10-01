"""Cross-assembly names must not create false FASTA/core mismatches."""

import sqlite3

import pytest

from fusion_function import data


def database():
    db = sqlite3.connect(":memory:")
    db.executescript("""
        CREATE TABLE ensembl_coord_system (coord_system_id INTEGER PRIMARY KEY, version TEXT);
        CREATE TABLE ensembl_seq_region (seq_region_id INTEGER, name TEXT, coord_system_id INTEGER, length INTEGER);
        CREATE TABLE sequences (kind TEXT, sequence_id TEXT, length INTEGER);
        INSERT INTO ensembl_coord_system VALUES (3,'GRCh38'),(12,'GRCh37');
        INSERT INTO ensembl_seq_region VALUES
            (131744,'HSCHR10_1_CTG2',3,309802),
            (1000321812,'HSCHR10_1_CTG2',12,135582047);
        INSERT INTO sequences VALUES ('dna','HSCHR10_1_CTG2',309802);
    """)
    return db


def test_real_release_116_cross_assembly_collision():
    with database() as db:
        # This is the erroneous join that failed on the user's real reference.
        assert db.execute(
            "SELECT r.length FROM ensembl_seq_region r JOIN sequences s "
            "ON r.name=s.sequence_id WHERE r.length != s.length"
        ).fetchone() == (135582047,)
        data.validate_dna_regions(db)


def test_rejects_actual_grch38_length_mismatch():
    with database() as db:
        db.execute("UPDATE sequences SET length=309803")
        with pytest.raises(
            ValueError, match=r"FASTA length 309803; GRCh38 core lengths \[309802\]"
        ):
            data.validate_dna_regions(db)


def test_wrong_assembly_exact_length_does_not_hide_mismatch():
    with database() as db:
        db.execute("UPDATE sequences SET length=135582047")
        with pytest.raises(ValueError, match="GRCh38 core lengths"):
            data.validate_dna_regions(db)


def test_requires_matching_region_in_supported_assembly():
    with database() as db:
        db.execute("DELETE FROM ensembl_seq_region WHERE coord_system_id=3")
        with pytest.raises(ValueError, match="no matching region"):
            data.validate_dna_regions(db)


def test_multiple_coordinate_systems_accept_correct_assembly_match():
    with database() as db:
        db.execute("INSERT INTO ensembl_coord_system VALUES (4,'GRCh38')")
        db.execute("INSERT INTO ensembl_seq_region VALUES (2,'HSCHR10_1_CTG2',4,135582047)")
        data.validate_dna_regions(db)
