# Reference configuration

Build the database with `fusion-function prepare-data`. See the README for release selection, caching, offline analysis, and regenerating derived annotations.

Reference selection precedence is: explicit `reference=`, explicit `database=`, `FUSION_FUNCTION_DB`, then the cache path. A supplied reader cannot be combined with database/release arguments. `FUSION_FUNCTION_CACHEDIR` changes the default cache; `--cache-dir` applies to preparation. If a custom cache directory is passed only at preparation time, also set the environment variable for subsequent analysis.

Analysis chooses the newest installed release numerically, or an explicit `release=`. Setup chooses the latest numbered FTP release when no release is supplied. Missing or unprocessed references fail clearly; annotation never downloads data implicitly.

Readers open SQLite in read-only mode and validate the reference format, human species, GRCh38 assembly, preprocessing version, and requested release. Explicit transcript versions must match the stored version.

`ReferenceDatabase(..., cached_chunks=64)` controls the bounded LRU genome cache. Use `cached_chunks=0` to disable it. Context managers close explicit readers; `close_default_reference()` closes the current thread's automatic reader. Database replacement or configuration changes cause automatic readers to reopen on the next call.

Preparation uses atomic replacement for a new database and transactions for derived-table regeneration. Source-file checksums and preprocessing metadata remain in the database. The provided test fixture has sparse genomic chunks and is exclusively for regression tests.

During preprocessing, SQLite's `secure_delete` setting is temporarily disabled.
The data are public reference annotations, so obsolete table contents need no
secure erasure. This avoids rewriting large tables and their rollback journals
just to zero deleted content. Transactional rollback remains enabled, and the
previous setting is restored on success or failure. Each derived-table deletion
has its own elapsed-time display.

Prepared transcript models use lossless zlib-compressed JSON at level 1 in an
ordinary rowid table. This reduces payload size and overflow-page traversal;
large models are not stored in a `WITHOUT ROWID` primary-key tree. Errors remain
plain JSON for SQL inspection. The reader returns the same model dictionaries
and also reads existing plain-JSON references. For direct SQL model inspection,
decode `payload` with `fusion_function.data.decode_transcript_payload`.
The database records `transcript_payload_codec=zlib-json-v1`.

This storage change does not accelerate deletion of an already existing large,
uncompressed table. Its first reprocessing still has that one-time cleanup cost.
Preprocessing logs the effective secure-delete, journal, synchronous, auto-vacuum
and page-size settings. Transactional rollback and durability settings are preserved.

Full builds save progress in `<output>.building`, beside the final database.
Failures and interrupts retain this checkpoint. Repeat the same preparation
command to resume: completed tables and FASTA files are skipped, and an
interrupted FASTA import reuses committed complete records from the same
checksum-verified input. Compressed FASTA must still be reread to reach the next
record, but completed sequences are not recompressed or inserted again.
An unfinished table is reimported; unfinished transcript preprocessing is
restarted as an atomic transaction. Completed preprocessing is reused after a
later failure unless its implementation or metadata input changes. Integrity
checking always runs before publication. Existing final references remain in
place until a completed checkpoint atomically replaces them.

Checkpoint settings pin the FTP root, release, assembly, schema checksum and
storage format. Incompatible checkpoints fail clearly; use another output path
or remove the `.building` database to start again. `--force` permits replacing
the final reference; it does not discard a resumable checkpoint. The persistent
`<output>.prepare.lock` file prevents concurrent builders from writing the same
checkpoint; its OS lock is released when the process exits. Both paths are
printed at startup.

Preparation always loads InterPro entry metadata, including `--preprocess-only`. By default it fetches the entry list or reuses a checksum-verified cached download. Use `--interpro-entries FILE` to supply local metadata for offline preparation. Missing entry types on usable transcripts fail preparation. Runtime annotation still never accesses the network.

Ensembl can retain accessions omitted from the current InterPro list. Default fetching searches archived entry lists for missing types before transcript preprocessing, using the newest available historical metadata for each accession. Current types are never overwritten. The database records `interpro_historical_entries` (accession to source release) and `interpro_historical_sources` (release, URL, SHA-256). Archived files and partial downloads are stored under `<cache>/metadata/interpro/releases/<release>/`. Explicit local lists remain offline. Reprocessing preserves previously recovered historical metadata unless the supplied list overrides it.

Protein sources use `analysis.db` (for example, `PANTHER`), with `logic_name` only as a fallback. The pipeline name `hmmpanther` is not a distinct annotation source. Preparation fetches PANTHER's human classification file for the `analysis.db_version` recorded by Ensembl and caches it under `<cache>/metadata/panther/<version>/`. `--panther-classifications FILE` supplies a local TSV instead. The database stores the UniProt/subfamily lookup in `ff_panther` and its URL or local path, version when fetched, and SHA-256 in `panther_classifications_source`.

AlphaFold and SIFTS structure mappings are excluded from functional features,
even when their intervals fit the peptide. They are recognized by pipeline or
database name before checking feature bounds, so an oversized structure mapping
does not reject an otherwise usable transcript. Functional domain/site intervals
still require valid peptide coordinates; they are never silently clipped.
References prepared with this filter record `feature_annotation_version=3`.

The downloader accepts both `PTHR<version>_human` and the historical
`PTHR<version>_human_` filename (used by PANTHER 14.1), always within the exact
requested release. An existing verified cache takes precedence; a 404 allows
the alternate filename, while other errors propagate. Automatic PANTHER
metadata is checked before genome/peptide imports, so unavailable metadata does
not first consume the full sequence-import time.

The parser detects both human classification layouts: PANTHER 14.1 stores the
subfamily ID and name in columns 3 and 5, while newer files use columns 4 and 6.
This also applies to local TSV overrides. UniProt accessions come from the
first-column identifier; family descriptions and ontology terms are not used
as subfamily names. Records without a named subfamily are omitted.

Subfamilies enrich existing PANTHER feature intervals only when exact protein cross-references and family IDs agree. Ambiguous assignments are omitted; gene names do not propagate assignments across proteins. No PANTHER API calls or new domain coordinates are needed. Programmatic reprocessing without a new classification file can reuse `ff_panther`; fully offline command-line preparation requires cached downloads or local InterPro, PANTHER and UniProt files.

To regenerate derived annotations from imported Ensembl data:

```bash
fusion-function prepare-data --preprocess-only /path/to/ensembl.sqlite
```

This regenerates derived records from the stored Ensembl tables and sequences. It does not repeat the genome import or the full-build SQLite integrity check.

Long SQLite phases report their names, elapsed time and completion. The terminal timer refreshes independently of SQLite progress callbacks, so table replacement, query preparation, indexing, commits and full integrity checks remain visible. Redirected logs receive a heartbeat every 30 seconds. These phases have no reliable completion percentage or ETA; the elapsed timer indicates a pending operation, not how much work has finished.

## Reviewed human UniProt features

Preparation fetches reviewed human XML with UniProt's bulk download endpoint. It is cached as `<cache>/metadata/uniprot/stream`, with a SHA-256 sidecar; partial downloads are beside it. `--uniprot-features FILE` accepts local XML or gzip regardless of its extension. UniProt metadata releases are independent of the selected Ensembl release. The input URL/path and checksum are saved in `uniprot_features_source`; `uniprot_counts` records import and sequence-match counts.

Only reviewed human domains, active sites, binding sites, functional sites and motifs are imported. Coordinates must be exact. Ensembl protein IDs in UniProt, or translation-level UniProt cross-references in Ensembl, identify candidates. A complete peptide sequence match is mandatory, including when mapping an older Ensembl release. Gene names and sequence length alone do not establish a mapping.

Explicitly described isoforms are reconstructed from nonoverlapping UniProt splice-variant records. Canonical annotations transfer only through unchanged residues; coordinates shift for edits before the feature. Features touching changed/deleted residues are omitted. Explicit isoform-specific annotations use that isoform's coordinates. External/undescribed isoforms and sequence mismatches are omitted, rather than approximately aligned.

`ff_uniprot_features` stores verified lookups by Ensembl peptide sequence ID. The lookup, transcript records and provenance are replaced in one transaction; a failure preserves the previous reference. Programmatic `preprocess_reference` calls can reuse this persisted lookup without downloading. The CLI defaults to fetching/reusing the bulk XML, including during reprocessing.

All feature annotations are stored locally and evaluated without runtime requests.

Repeated preparation reuses sequence-verified UniProt lookups when the XML
checksum, mapping implementation, recorded Ensembl source-file checksums and
reference release/storage metadata are unchanged. Changes invalidate reuse.
References without recorded source checksums are reverified. Transcript models
and their genomic feature mappings are still regenerated.

## Transcript errors

For a 5'-incomplete CDS, Ensembl's starting exon phase identifies one or two
missing leading codon bases. Prepared CDS coordinates preserve that peptide
offset while mapping only bases present in the genome. Length validation includes
the missing bases and accepts the optional terminal stop codon; it does not use
a general mismatch tolerance. Nonzero offsets are stored as `cds_start_phase`.
Features crossing the incomplete first codon have only their available bases
mapped and retained. That codon cannot supply a native start; initiation uses
the existing alternative-ATG heuristic. Other length mismatches remain errors.
These references record `cds_mapping_version=2`.

Noncoding transcripts have their own status and are not included in the error
count. Errors identify models the package cannot currently use, including
unsupported coordinate systems, missing peptides and CDS/peptide mismatches.
Counts alone do not establish whether the causes are expected source limitations
or preprocessing defects.

Preparation logs grouped reasons and up to three example transcript IDs per
reason, and saves the summary in `preprocessing_errors`. Full messages remain
in the transcript payloads. Inspect a completed reference without rebuilding:

```python
from fusion_function import ReferenceDatabase

with ReferenceDatabase(release=116) as reference:
    for reason, details in reference.transcript_error_summary().items():
        print(details["count"], reason, details["example_transcripts"])
```

If no saved summary exists, the reader derives it from the stored error records.
Failed transcript annotation returns an error; an unsupported model is not
evidence that a fusion lacks functional features.
