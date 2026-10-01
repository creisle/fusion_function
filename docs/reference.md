# Reference configuration

Build the database with `fusion-function prepare-data`. See the README for release selection, caching, offline analysis, and upgrading an existing database.

Reference selection precedence is: explicit `reference=`, explicit `database=`, `FUSION_FUNCTION_DB`, then the cache path. A supplied reader cannot be combined with database/release arguments. `FUSION_FUNCTION_CACHEDIR` changes the default cache; `--cache-dir` applies to preparation. If a custom cache directory is passed only at preparation time, also set the environment variable for subsequent analysis.

Analysis chooses the newest installed release numerically, or an explicit `release=`. Setup chooses the latest numbered FTP release when no release is supplied. Missing or unprocessed references fail clearly; annotation never downloads data implicitly.

Readers open SQLite in read-only mode and validate the reference format, human species, GRCh38 assembly, preprocessing version, and requested release. Explicit transcript versions must match the stored version. Prior standalone databases without an assembly metadata key are accepted as GRCh38 because that preprocessor only emitted GRCh38 models.

`ReferenceDatabase(..., cached_chunks=64)` controls the bounded LRU genome cache. Use `cached_chunks=0` to disable it. Context managers close explicit readers; `close_default_reference()` closes the current thread's automatic reader. Database replacement or configuration changes cause automatic readers to reopen on the next call.

Preparation uses atomic replacement for a new database and transactions for derived-table regeneration. Source-file checksums and preprocessing metadata remain in the database. The provided test fixture has sparse genomic chunks and is exclusively for regression tests.

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
printed at startup. A failed build made by older versions cannot recover raw
imports if that version already deleted its temporary SQLite database.

Preparation always loads InterPro entry metadata, including `--preprocess-only`. By default it fetches the entry list or reuses a checksum-verified cached download. Use `--interpro-entries FILE` to supply local metadata for offline preparation. Missing entry types on usable transcripts fail preparation. Runtime annotation still never accesses the network.

Ensembl can retain accessions omitted from the current InterPro list. Default fetching searches archived entry lists for missing types before transcript preprocessing, using the newest available historical metadata for each accession. Current types are never overwritten. The database records `interpro_historical_entries` (accession to source release) and `interpro_historical_sources` (release, URL, SHA-256). Archived files and partial downloads are stored under `<cache>/metadata/interpro/releases/<release>/`. Explicit local lists remain offline. Reprocessing preserves previously recovered historical metadata unless the supplied list overrides it.

Protein sources use `analysis.db` (for example, `PANTHER`), with `logic_name` only as a fallback. The pipeline name `hmmpanther` is not a distinct annotation source. Preparation fetches PANTHER's human classification file for the `analysis.db_version` recorded by Ensembl and caches it under `<cache>/metadata/panther/<version>/`. `--panther-classifications FILE` supplies a local TSV instead. The database stores the UniProt/subfamily lookup in `ff_panther` and its URL or local path, version when fetched, and SHA-256 in `panther_classifications_source`.

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

To update feature types, add subfamily annotations and import UniProt features into an existing reference, install the updated package and run:

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

The added annotations increase output coverage without adding runtime requests or dependencies. Previously built references remain readable and receive the corrected filtering/boundary behavior immediately. They need `--preprocess-only` once to acquire UniProt features.
