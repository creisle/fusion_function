import gzip
import sqlite3
from pathlib import Path
from fusion_function import data as e


def fixture(path):
    db = sqlite3.connect(path)
    definitions = {
        "transcript": "transcript_id INTEGER, stable_id TEXT, version INTEGER, seq_region_id INTEGER, seq_region_start INTEGER, seq_region_end INTEGER, seq_region_strand INTEGER, canonical_translation_id INTEGER, is_current INTEGER",
        "translation": "translation_id INTEGER, transcript_id INTEGER, stable_id TEXT, version INTEGER, seq_start INTEGER, start_exon_id INTEGER, seq_end INTEGER, end_exon_id INTEGER",
        "exon": "exon_id INTEGER, seq_region_id INTEGER, seq_region_start INTEGER, seq_region_end INTEGER, seq_region_strand INTEGER, phase INTEGER, end_phase INTEGER",
        "exon_transcript": "exon_id INTEGER, transcript_id INTEGER, rank INTEGER",
        "seq_region": "seq_region_id INTEGER, name TEXT, coord_system_id INTEGER, length INTEGER",
        "coord_system": "coord_system_id INTEGER, version TEXT",
        "analysis": "analysis_id INTEGER, logic_name TEXT",
        "protein_feature": "protein_feature_id INTEGER, translation_id INTEGER, seq_start INTEGER, seq_end INTEGER, hit_name TEXT, hit_description TEXT, analysis_id INTEGER",
        "interpro": "interpro_ac TEXT, id TEXT",
        "xref": "dbprimary_acc TEXT, description TEXT, display_label TEXT",
    }
    for name, cols in definitions.items():
        db.execute(f"CREATE TABLE ensembl_{name} ({cols})")
    db.execute("INSERT INTO ensembl_seq_region VALUES (1,'1',1,1000)")
    db.execute("INSERT INTO ensembl_coord_system VALUES (1,'GRCh38')")
    for tid, strand in [(1, 1), (2, -1), (4, 1)]:
        db.execute(
            "INSERT INTO ensembl_transcript VALUES (?,?,?,?,?,?,?,?,?)",
            (tid, f"ENST{tid:011}", 3, 1, 100, 159, strand, tid, 1),
        )
        starts = [(100, 111), (148, 159)] if strand == 1 else [(148, 159), (100, 111)]
        for rank, (a, b) in enumerate(starts, 1):
            eid = tid * 10 + rank
            db.execute(
                "INSERT INTO ensembl_exon VALUES (?,?,?,?,?,?,?)", (eid, 1, a, b, strand, -1, -1)
            )
            db.execute("INSERT INTO ensembl_exon_transcript VALUES (?,?,?)", (eid, tid, rank))
        db.execute(
            "INSERT INTO ensembl_translation VALUES (?,?,?,?,?,?,?,?)",
            (tid, tid, f"ENSP{tid:011}", 2, 4, tid * 10 + 1, 9, tid * 10 + 2),
        )
        for fid, a, b, hit in [(tid * 10 + 1, 2, 5, "PF1"), (tid * 10 + 2, 3, 5, "PF2")]:
            db.execute(
                "INSERT INTO ensembl_protein_feature VALUES (?,?,?,?,?,?,?)",
                (fid, tid, a, b, hit, "member name", 1),
            )
    db.execute("INSERT INTO ensembl_transcript VALUES (3,'ENST00000000003',1,1,200,211,1,NULL,1)")
    db.execute("INSERT INTO ensembl_exon VALUES (31,1,200,211,1,-1,-1)")
    db.execute("INSERT INTO ensembl_exon_transcript VALUES (31,3,1)")
    db.execute("INSERT INTO ensembl_analysis VALUES (1,'Pfam')")
    db.executemany("INSERT INTO ensembl_interpro VALUES (?,?)", [("IPR1", "PF1"), ("IPR1", "PF2")])
    db.execute("INSERT INTO ensembl_xref VALUES ('IPR1','Domain from core','IPR1')")
    e.create_sequence_tables(db)
    fa = path.with_suffix(".fa.gz")
    with gzip.open(fa, "wb") as f:
        f.write(b">1\n" + b"ACGT" * 250 + b"\n")
    e.import_fasta(db, "dna", fa)
    with gzip.open(fa, "wb") as f:
        f.write(
            b">ENSP00000000001.2\nMKLFIN\n>ENSP00000000002.2\nMKLFIN\n>ENSP00000000004.2\nMKLFINS\n"
        )
    e.import_fasta(db, "pep", fa)
    db.execute("CREATE TABLE build_metadata (key TEXT PRIMARY KEY,value TEXT NOT NULL)")
    db.executemany(
        "INSERT INTO build_metadata VALUES (?,?)",
        [
            ("format_version", "1"),
            ("species", "homo_sapiens"),
            ("release", "116"),
            ("assembly", "GRCh38"),
            ("sequence_chunk_size", str(e.CHUNK_SIZE)),
            ("sequence_codec", "zlib"),
        ],
    )
    db.commit()
    return db


def uniprot_file(tmp_path: Path) -> Path:
    """Small reviewed human XML with no mapped features for offline CLI tests."""
    path = tmp_path / "uniprot.xml"
    path.write_text(
        '<uniprot xmlns="http://uniprot.org/uniprot">'
        '<entry dataset="Swiss-Prot"><accession>PTEST</accession>'
        '<organism><dbReference type="NCBI Taxonomy" id="9606"/></organism>'
        "<sequence>MKLFIN</sequence></entry></uniprot>"
    )
    return path
