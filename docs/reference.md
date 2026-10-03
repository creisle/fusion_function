# Reference configuration

For basic installation and annotation, see the [quick start](../README.md#quick-start). This page covers reference configuration, preparation, updates and debugging.

## Selecting a reference

Reference selection precedence is: explicit `reference=`, explicit `database=`, `FUSION_FUNCTION_DB`, then the cache path. A supplied reader cannot be combined with database/release arguments. `FUSION_FUNCTION_CACHEDIR` changes the default cache; `--cache-dir` applies to preparation. If a custom cache directory is passed only at preparation time, also set the environment variable for subsequent analysis.

Analysis chooses the newest installed release numerically, or an explicit `release=`. Setup chooses the latest numbered FTP release when no release is supplied. Missing or unprocessed references fail clearly; annotation never downloads data implicitly.

Readers open SQLite in read-only mode and validate the reference format, human species, GRCh38 assembly, preprocessing version, and requested release. Explicit transcript versions must match the stored version.

`ReferenceDatabase(..., cached_chunks=64)` controls the bounded LRU genome cache. Use `cached_chunks=0` to disable it. Context managers close explicit readers; `close_default_reference()` closes the current thread's automatic reader. Database replacement or configuration changes cause automatic readers to reopen on the next call.

## Installing the reference

```bash
# Latest numbered Ensembl release
fusion-function prepare-data

# Select a release
fusion-function prepare-data --release 116
```

`prepare-data` first looks for a compatible prebuilt human GRCh38 reference for the selected Ensembl release. Without `--release`, it resolves the latest numbered Ensembl FTP release; it does not substitute an older published reference. If no compatible artifact is registered, it logs that decision and builds from Ensembl source files. The initial catalog is empty until the first reference is published.

Prebuilt installation has four numbered steps: download and verify the archive, expand the database, validate its format and release, and install it. Download and expansion progress bars show byte totals in a terminal. The selected URL, build revision, destination, archive size, expanded size and partial-file locations are logged. A failed artifact download or checksum check fails clearly rather than silently starting a source build.

Interrupted downloads resume from `<cache>/prebuilt/<archive-sha256>/ensembl.sqlite.gz.part`; if the server ignores the resume request, the download restarts. Verified archives remain beside that file. Expansion uses `<output>.prebuilt.part` beside the destination and removes that temporary file on failure. Both compressed and expanded files require disk space; replacing an existing database temporarily also requires space for the old copy. SHA-256 and byte sizes are checked for both files, and the existing reference stays in place until successful atomic replacement. The exporter performs full SQLite integrity checking before publication; clients verify the exact artifact and runtime compatibility without repeating that long scan.

An existing compatible database is reused unless `--force` is supplied. Published build revisions distinguish updated annotations for the same Ensembl release. Use `fusion-function prepare-data --release 116 --force` to install its newest compatible published revision.

To build from source explicitly:

```bash
fusion-function prepare-data --from-source --release 116
```

Source builds download Ensembl core tables, genomic and peptide FASTA, and annotation metadata, then import and preprocess them locally. They need several GB of downloads and working space and can take an hour or longer. Explicit `--base-url`, `--interpro-entries`, `--panther-classifications` or `--uniprot-features` overrides also select source preparation, so custom inputs are not ignored.

The default catalog is `https://raw.githubusercontent.com/creisle/fusion_function/main/src/fusion_function/reference_catalog.json`, with a bundled copy used if that URL is unavailable. It can be updated independently of package releases. Set `FUSION_FUNCTION_REFERENCE_CATALOG` or pass `--reference-catalog FILE_OR_HTTPS_URL` to select another catalog; explicit catalogs must load successfully. Catalog entries pin the species, assembly, Ensembl release, build revision, storage/preprocessing versions, download URL, sizes and SHA-256 checksums. Only compatible entries are selected. A catalog cannot be combined with explicit source-build options.

The default database path is `<cache>/homo_sapiens/GRCh38/release-<release>/ensembl.sqlite`. The cache defaults to `~/.cache/fusion_function` on Linux, `~/Library/Caches/fusion_function` on macOS, or `%LOCALAPPDATA%/fusion_function` on Windows. Set `FUSION_FUNCTION_CACHEDIR` to change it for both preparation and analysis. `--cache-dir` overrides it for preparation only.

Preparation uses atomic replacement for a new database and transactions for derived-table regeneration. Source-file checksums and preprocessing metadata remain in the database. The provided test fixture has sparse genomic chunks and is exclusively for regression tests.

Prepared transcript models use lossless zlib-compressed JSON at level 1 in an ordinary rowid table. This reduces payload size and overflow-page traversal; large models are not stored in a `WITHOUT ROWID` primary-key tree. Errors remain plain JSON for SQL inspection. The reader returns the same model dictionaries and also reads existing plain-JSON references. For direct SQL model inspection, decode `payload` with `fusion_function.data.decode_transcript_payload`. The database records `transcript_payload_codec=zlib-json-v1`.

## Resuming a full build

Source builds save progress in `<output>.building`, beside the final database. Failures and interrupts retain this checkpoint. Repeat the same preparation command to resume: completed tables and FASTA files are skipped, and an interrupted FASTA import reuses committed complete records from the same checksum-verified input. Compressed FASTA must still be reread to reach the next record, but completed sequences are not recompressed or inserted again. An unfinished table is reimported; unfinished transcript preprocessing is restarted as an atomic transaction. Completed preprocessing is reused after a later failure unless its implementation or metadata input changes. Integrity checking always runs before publication. Existing final references remain in place until a completed checkpoint atomically replaces them.

Checkpoint settings pin the FTP root, release, assembly, schema checksum and storage format. Incompatible checkpoints fail clearly; use another output path or remove the `.building` database to start again. `--force` permits replacing the final reference; it does not discard a resumable checkpoint. The persistent `<output>.prepare.lock` file prevents concurrent builders from writing the same checkpoint; its OS lock is released when the process exits. Both paths are printed at startup.

## Annotation inputs

Source preparation always loads InterPro entry metadata; prebuilt references already contain it. By default it fetches the entry list or reuses a checksum-verified cached download. Use `--interpro-entries FILE` to supply local metadata for offline preparation. Missing entry types on usable transcripts fail preparation. Runtime annotation still never accesses the network.

Ensembl can retain accessions omitted from the current InterPro list. Default fetching searches archived entry lists for missing types before transcript preprocessing, using the newest available historical metadata for each accession. Current types are never overwritten. The database records `interpro_historical_entries` (accession to source release) and `interpro_historical_sources` (release, URL, SHA-256). Archived files and partial downloads are stored under `<cache>/metadata/interpro/releases/<release>/`. Explicit local lists remain offline.

Protein sources use `analysis.db` (for example, `PANTHER`), with `logic_name` only as a fallback. The pipeline name `hmmpanther` is not a distinct annotation source. Preparation fetches PANTHER's human classification file for the `analysis.db_version` recorded by Ensembl and caches it under `<cache>/metadata/panther/<version>/`. `--panther-classifications FILE` supplies a local TSV instead. The database stores the UniProt/subfamily lookup in `ff_panther` and its URL or local path, version when fetched, and SHA-256 in `panther_classifications_source`.

AlphaFold and SIFTS structure mappings are excluded from functional features, even when their intervals fit the peptide. They are recognized by pipeline or database name before checking feature bounds, so an oversized structure mapping does not reject an otherwise usable transcript. Functional domain/site intervals still require valid peptide coordinates; they are never silently clipped. References prepared with this filter record `feature_annotation_version=3`.

The downloader accepts both `PTHR<version>_human` and the historical `PTHR<version>_human_` filename (used by PANTHER 14.1), always within the exact requested release. An existing verified cache takes precedence; a 404 allows the alternate filename, while other errors propagate. Automatic PANTHER metadata is checked before genome/peptide imports, so unavailable metadata does not first consume the full sequence-import time.

The parser detects both human classification layouts: PANTHER 14.1 stores the subfamily ID and name in columns 3 and 5, while newer files use columns 4 and 6. This also applies to local TSV overrides. UniProt accessions come from the first-column identifier; family descriptions and ontology terms are not used as subfamily names. Records without a named subfamily are omitted.

Subfamilies enrich existing PANTHER feature intervals only when exact protein cross-references and family IDs agree. Ambiguous assignments are omitted; gene names do not propagate assignments across proteins. No PANTHER API calls or new domain coordinates are needed. Fully offline preparation requires cached downloads or local InterPro, PANTHER and UniProt files.

SQLite phases log their start and outcome, with no background timer, periodic elapsed-time output or heartbeat. Operations that can exceed ten minutes print an advance notice; transcript preprocessing and removal of old transcript models can take an hour or longer on some systems. These notices describe possible durations, not ETAs. Downloads and record imports retain rate-limited progress bars in a terminal; redirected logs suppress those redraws.

## Reviewed human UniProt features

Preparation fetches reviewed human XML with UniProt's bulk download endpoint. It is cached as `<cache>/metadata/uniprot/stream`, with a SHA-256 sidecar; partial downloads are beside it. `--uniprot-features FILE` accepts local XML or gzip regardless of its extension. UniProt metadata releases are independent of the selected Ensembl release. The input URL/path and checksum are saved in `uniprot_features_source`; `uniprot_counts` records import and sequence-match counts.

Only reviewed human domains, active sites, binding sites, functional sites and motifs are imported. Coordinates must be exact. Ensembl protein IDs in UniProt, or translation-level UniProt cross-references in Ensembl, identify candidates. A complete peptide sequence match is mandatory, including when mapping an older Ensembl release. Gene names and sequence length alone do not establish a mapping.

Explicitly described isoforms are reconstructed from nonoverlapping UniProt splice-variant records. Canonical annotations transfer only through unchanged residues; coordinates shift for edits before the feature. Features touching changed/deleted residues are omitted. Explicit isoform-specific annotations use that isoform's coordinates. External/undescribed isoforms and sequence mismatches are omitted, rather than approximately aligned.

`ff_uniprot_features` stores verified lookups by Ensembl peptide sequence ID. The source and sequence-match counts are recorded in the reference; see [Updating an existing reference](#updating-an-existing-reference) for lookup reuse during updates.

All feature annotations are stored locally and evaluated without runtime requests.

## Updating an existing reference

### Update the package

```bash
pip install --upgrade fusion-function
```

Use the updated package for both preparation and analysis. New references can use a storage format that older readers cannot decode.

### Regenerate annotations for the same Ensembl release

For a full source-built database, regenerate transcript models, protein features and splice-site annotations from the Ensembl tables and sequences already stored in it:

```bash
fusion-function prepare-data --preprocess-only /path/to/ensembl.sqlite
```

Compact prebuilt references omit the raw Ensembl tables and cannot use `--preprocess-only`. Install a newer prebuilt with `--release RELEASE --force`, or create a full source reference with `--from-source --release RELEASE --force` before reprocessing. Exporting a compact reference leaves the original full database intact.

Reprocessing a full source-built database does not repeat the genome import or the full-build SQLite integrity check. It keeps the database's Ensembl release and assembly; `--release`, `--output`, `--force`, `--from-source` and `--reference-catalog` cannot be combined with `--preprocess-only`.

Use this after changes to preprocessing logic or annotation inputs. It is unnecessary for documentation, logging or other changes that do not affect the prepared records.

The default database path is `<cache>/homo_sapiens/GRCh38/release-<release>/ensembl.sqlite`. See [Selecting a reference](#selecting-a-reference) for cache settings and database selection.

### Annotation inputs and reuse

Reprocessing loads InterPro metadata, PANTHER classifications for the version recorded by Ensembl, and reviewed human UniProt features. It reuses checksum-verified cached files when available. Rerunning alone does not necessarily fetch a newer metadata snapshot.

Ensembl release selection does not pin the independently released InterPro or UniProt metadata. Source paths/URLs and checksums are recorded in the database. Supply `--interpro-entries FILE`, `--panther-classifications FILE` and `--uniprot-features FILE` to select local inputs or prepare offline.

Reprocessing preserves previously recovered historical metadata unless the supplied list overrides it.

Repeated preparation reuses sequence-verified UniProt lookups when the XML checksum, mapping implementation, recorded Ensembl source-file checksums and reference release/storage metadata are unchanged. Changes invalidate reuse. References without recorded source checksums are reverified. Transcript models and their genomic feature mappings are still regenerated.

Programmatic reprocessing without a new classification file can reuse `ff_panther`. Programmatic `preprocess_reference` calls can also reuse stored UniProt lookups without downloading.

### Transactions and long-running work

Derived transcript records, lookup tables and provenance are replaced in one transaction. A failure or interrupt rolls back to the previous reference. If interrupted before commit, rerunning repeats preprocessing; it does not resume partway through that transaction.

Transcript preprocessing and removal of old models can take an hour or longer on some systems. Advance notices describe possible durations rather than completion estimates.

During preprocessing, SQLite's `secure_delete` setting is temporarily disabled. The data are public reference annotations, so obsolete table contents need no secure erasure. This avoids rewriting large tables and their rollback journals just to zero deleted content. Transactional rollback remains enabled, and the previous setting is restored on success or failure. Each derived-table deletion logs its start and outcome.

This storage change does not accelerate deletion of an already existing large, uncompressed table. Its first reprocessing still has that one-time cleanup cost. Preprocessing logs the effective secure-delete, journal, synchronous, auto-vacuum and page-size settings. Transactional rollback and durability settings are preserved.

Compressed models can reuse freed pages, but an in-place update does not necessarily shrink the database file.

### Change or rebuild an Ensembl release

Install a different Ensembl release using `fusion-function prepare-data --release RELEASE`; a compatible prebuilt is downloaded when available, otherwise a source build runs. It has its own default database path. Analysis selects the newest prepared local release unless you pin it with `release=` or select a database explicitly.

To rebuild an existing release from source files, add `--from-source --force`. Verified cached downloads are reused. The previous final database stays in place until a completed build passes integrity checking and is published. `--force` does not discard an existing `.building` checkpoint; see [Resuming a full build](#resuming-a-full-build).

## Debugging transcript errors

For a 5'-incomplete CDS, Ensembl's starting exon phase identifies one or two missing leading codon bases. Prepared CDS coordinates preserve that peptide offset while mapping only bases present in the genome. Length validation includes the missing bases and accepts the optional terminal stop codon; it does not use a general mismatch tolerance. Nonzero offsets are stored as `cds_start_phase`. Features crossing the incomplete first codon have only their available bases mapped and retained. That codon cannot supply a native start; initiation uses the existing alternative-ATG heuristic. Other length mismatches remain errors. These references record `cds_mapping_version=2`.

Noncoding transcripts have their own status and are not included in the error count. Errors identify models the package cannot currently use, including unsupported coordinate systems, missing peptides and CDS/peptide mismatches. Counts alone do not establish whether the causes are expected source limitations or preprocessing defects.

Preparation logs grouped reasons and up to three example transcript IDs per reason, and saves the summary in `preprocessing_errors`. Full messages remain in the transcript payloads. Inspect a completed reference without rebuilding:

```python
from fusion_function import ReferenceDatabase

with ReferenceDatabase(release=116) as reference:
    for reason, details in reference.transcript_error_summary().items():
        print(details["count"], reason, details["example_transcripts"])
```

If no saved summary exists, the reader derives it from the stored error records. Failed transcript annotation returns an error; an unsupported model is not evidence that a fusion lacks functional features.

## Publishing a prebuilt reference

Keep the full source-built database for maintenance and reprocessing. Run its built-reference tests before exporting:

```bash
pytest -m human_reference --human-reference-db=/path/to/ensembl.sqlite -v
fusion-function export-reference --output reference-exports --build-revision 1
```

With no database path, export uses `FUSION_FUNCTION_DB` if set, otherwise the newest installed local release in the default cache, including a cache selected by `FUSION_FUNCTION_CACHEDIR`. Use `--release 116` to select an installed Ensembl release, or supply a database path explicitly. An explicit path takes precedence over `FUSION_FUNCTION_DB`; `--release` is checked against either selected database. Export never downloads or builds a missing reference. The selected source path is logged before export.

```bash
fusion-function export-reference --release 116 --output reference-exports --build-revision 1
```

`--release` selects the Ensembl release being exported; `--build-revision` labels this exported build of that release.

The export creates a fresh SQLite snapshot containing transcript models, InterPro metadata, genomic and peptide sequence chunks, source checksums when present, and build provenance. Raw Ensembl tables and preprocessing lookup tables are omitted. The source is opened read-only and remains unchanged. The fresh database is integrity-checked, compressed with gzip, and accompanied by a JSON manifest containing exact compatibility metadata, sizes, hashes and source notices. Temporary paths are printed. Copying, integrity checking and compression can each take over ten minutes for a large reference.

To publish from your `creisle` Zenodo account:

1. Create a Zenodo upload and reserve its record ID. The account username identifies the owner; downloads require an exact public record/file URL.
2. Export with `--zenodo-record RECORD_ID` to populate that URL in the manifest. Without this option, fill the manifest's `url` after publication; entries with a null URL are not downloadable.
3. Upload the `.sqlite.gz` and `.json` files, include source attribution and applicable notices, and publish the record. Resolve the data-license questions below before choosing a license for the combined artifact.
4. Append the manifest's reference entry to `src/fusion_function/reference_catalog.json` under `references`, keeping `catalog_version` equal to 1, and commit it to the repository's `main` branch. Test using `prepare-data --release RELEASE --reference-catalog MANIFEST.json --output /path/to/test.sqlite` with the published manifest.

Keep published artifacts immutable. Increment `--build-revision` for new preprocessing or metadata of the same Ensembl release; the exporter refuses to overwrite an existing export. Use a new Zenodo version/record and its exact file URL for that build. The catalog selects the highest compatible build revision for each release. Do not add a SQLite file to the Python wheel or Git history. Record Zenodo's DOI in the release documentation for citation and reproducibility.

## Data attribution and reuse

These source terms are separate from the package's software license. The export embeds source-credit and modification notices in `build_metadata` and its manifest, along with input provenance. These summaries are not a blanket license for the combined database or a substitute for any required upstream copyright notices. This review covers the sources currently imported; it does not cover RMC, MTR or AlphaMissense, which are not imported.

| Source | Published terms and implication for redistribution |
|---|---|
| [Ensembl](https://www.ensembl.org/info/about/legal/disclaimer.html) core tables and DNA/peptide FASTA | Ensembl imposes no restrictions on project-generated data, but explicitly preserves third-party constraints. Its Apache 2.0 software license does not license all database content. Retain release and source provenance. |
| [UniProt](https://www.uniprot.org/help/license) reviewed human features and sequences | Copyrightable database content is [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). Redistribution and adaptation are allowed with credit, a license link and identification of changes. Credit the UniProt Consortium and describe the derived feature mappings/compression. |
| [InterPro](https://interpro-documentation.readthedocs.io/en/latest/license.html) entry names and types | Current downloadable InterPro data are CC0 1.0. Historical entry lists should retain their source-release notices; current InterPro terms do not automatically license every member database's signature collections. |
| [PANTHER human classifications](https://data.pantherdb.org/ftp/sequence_classifications/) | PANTHER data and data products are reported as [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) by the [Reusable Data Project](https://reusabledata.org/panther.html), citing [PANTHER’s Terms of Use](https://pantherdb.org/tou.jsp). Credit PANTHER and its contributors, link to the license and identify modifications. Historical 14.1 and 17.0 FTP READMEs contain GPL-2.0-or-later notices; retain and reconcile those notices for those specific inputs rather than treating them as a blanket license for all PANTHER data. |
| Member annotations carried through Ensembl | The exact sources vary by release. [PROSITE](https://prosite.expasy.org/prosite_license.html) database terms are CC BY-NC-ND 4.0 with commercial licensing, and [SMART](https://smart.embl.de/about.cgi) models/alignments/thresholds require a license. The package imports derived match annotations rather than their models; confirm how those terms apply to redistribution of the matches. Do not assume Ensembl or InterPro removes these conditions. |

Before publishing, inventory member sources in the full database and verify the terms of the exact imported releases:

```sql
SELECT a.logic_name, a.db, a.db_version, COUNT(*) AS matches
FROM ensembl_protein_feature AS pf
JOIN ensembl_analysis AS a USING (analysis_id)
GROUP BY a.logic_name, a.db, a.db_version;
```

Retain the applicable license files, copyright notices, source URLs and checksums with the Zenodo deposit. Member-annotation terms still require separate review before assigning a single permissive license to the combined database.
