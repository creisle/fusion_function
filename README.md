# fusion_function

Predict fusion reading frames and retained, disrupted, or excluded functional protein features from human GRCh38 transcript breakpoints. Annotation uses a preprocessed local SQLite reference and makes no network requests. Preparation uses tqdm progress bars.

## Install and prepare the reference

```bash
pip install .

# Latest numbered Ensembl FTP release; default cache location
fusion-function prepare-data

# Pin a release for reproducibility
fusion-function prepare-data --release 116
```

The setup command downloads Ensembl core table dumps and genomic/peptide FASTA, imports them into SQLite, and precomputes exon/CDS coordinates, protein-feature mappings, canonical splice-site windows, and feature-group membership. Preparation commands also download InterPro's small entry list by default, including `--preprocess-only` for canonical names and entry types. Progress is logged throughout. The complete human reference needs several GB of downloads and working space.

The default database location is:

```text
<cache>/homo_sapiens/GRCh38/release-<release>/ensembl.sqlite
```

`<cache>` is the OS user cache (`~/.cache/fusion_function` on Linux, `~/Library/Caches/fusion_function` on macOS, `%LOCALAPPDATA%/fusion_function` on Windows). Set `FUSION_FUNCTION_CACHEDIR` or pass `--cache-dir` to change it. Analysis selects the newest prepared local release; it does not contact FTP or download data automatically. Supply `release=` to pin analysis to an installed release.

If preparation fails or is interrupted, repeat the same command to resume from
`ensembl.sqlite.building` beside the output. Completed table/FASTA imports and
completed preprocessing are reused. Interrupted FASTA imports reuse committed
records from unchanged input; an unfinished preprocessing transaction restarts.
The final reference is published only after integrity checking succeeds.

### Upgrade an existing database

A database created by the earlier standalone script remains compatible. Add InterPro metadata and regenerate the derived tables without downloading Ensembl again:

```bash
fusion-function prepare-data --preprocess-only /path/to/ensembl.sqlite
```

Preparation always loads InterPro metadata, including `--preprocess-only`. By default it downloads the entry list or reuses a checksum-verified cached copy. Use `--interpro-entries FILE` to supply local metadata for offline preparation. Ensembl release selection does not pin InterPro's independently released entry list: the imported file's checksum is recorded. Default fetching supplements missing types from archived InterPro entry lists, newest first, before processing transcripts. Current metadata takes precedence; historical source releases and checksums are recorded in `build_metadata`. Archives are cached under `<cache>/metadata/interpro/releases/<release>/`. Explicit local lists do not trigger archive downloads. Missing entry types still fail preparation while preserving any existing reference.

## Annotate a fusion

```python
from fusion_function import annotate_fusion_domains

result = annotate_fusion_domains(
    transcript1_id="ENST00000305877",
    transcript2_id="ENST00000318560",
    breakpoint1="22:23290413",
    breakpoint2="9:130854064",
    gene1_terminus="N",
    gene2_terminus="C",
    release=116,  # optional; otherwise newest prepared local release
)
```

Use `database="/path/to/ensembl.sqlite"` or `FUSION_FUNCTION_DB` for a database outside the default cache. An explicit `database=` takes precedence over that environment variable. Supplied transcript versions are checked against the prepared release.

For batches, reuse a reader:

```python
from fusion_function import ReferenceDatabase, annotate_fusion_domains

with ReferenceDatabase(release=116) as reference:
    result = annotate_fusion_domains(
        transcript1_id="ENST00000305877", transcript2_id="ENST00000318560",
        breakpoint1="22:23290413", breakpoint2="9:130854064",
        gene1_terminus="N", gene2_terminus="C", reference=reference,
    )
```

Default calls also reuse one SQLite connection per thread and cache at most 64 decompressed 1 MiB genome chunks. Use separate readers per worker thread/process. `close_default_reference()` releases the current thread's automatic reader. Full transcript sequences are not duplicated in the database or retained in this cache.

Results contain `frame_status`, `domains`, and a `translation_start` summary mapping initiating transcript IDs to status strings. Alternative product outcomes are slash-separated. Domain fields cover retention, splicing, frame, and premature termination. The former HTTP `session=` parameter has been removed. `inserted_sequence` remains supported; reference options are keyword-only.

## Reference limitations

Only human GRCh38 is supported. Breakpoint-dependent clipping, splice-product construction, frame comparison, and stop-codon detection remain runtime calculations. Unsupported transcript models, missing sequences, and CDS/peptide length mismatches are recorded in `ff_transcripts` with `status='error'`; noncoding transcripts return an error through the annotation API. Sequence-edited models are not silently corrected.

Preparation downloads PANTHER's human classification file for the version recorded in Ensembl's analysis metadata. Exact UniProt protein cross-references and matching PANTHER family IDs supply the subfamily names used by the existing family fallback logic. No gene-wide propagation or new domain intervals are inferred. The bulk lookup and source checksum are saved in the database; annotation makes no PANTHER API calls. Use `--panther-classifications FILE` for a local classification file. Member features without a resolved InterPro entry retain their source annotation; entry types are not guessed from the source database.

## Tests and documentation

```bash
pip install --group dev
pytest
```

Tests use a small, offline SQLite fixture converted from the supplied response cache, complete legacy fusion outputs, synthetic preprocessing fixtures, and network blocking. The sparse regression fixture is not a production reference and does not represent the whole genome or a verified Ensembl release. A full human FTP-to-database build still needs validation on production data.

The original integration cases run against both references, with `[bundled]` and `[built]` test IDs. This includes SLC45A2::AMACR and the PANTHER enrichment assertions; their expectations are unchanged. The bundled cases need no reference build. Built cases and the additional full-reference checks run automatically when a completed human GRCh38 database is available in the default cache, including a cache set with `FUSION_FUNCTION_CACHEDIR`. The newest installed release is selected. `FUSION_FUNCTION_DB` selects an existing database, and `--human-reference-db` overrides that selection:

```bash
pytest -m human_reference --human-reference-db=/path/to/ensembl.sqlite -v
```

Without a local build, the full-reference tests are skipped. To run only the lightweight tests even when a build is installed, use `pytest -m 'not human_reference'`.

- [API](docs/api.md)
- [Reference configuration](docs/reference.md)
- [Methodology](docs/methodology.md)
- [Testing](docs/testing.md)

## Citation

This package ports and extends selected fusion-annotation logic from [MAVIS](https://github.com/bcgsc/mavis). Please cite:

Reisle C, Mungall KL, Choo C, et al. MAVIS: merging, annotation, validation, and illustration of structural variants. *Bioinformatics*. 2019;35(3):515–517. PMID:30016509.

Predicted structural consequences do not establish fusion expression, oncogenicity, pathogenicity, or clinical actionability.
