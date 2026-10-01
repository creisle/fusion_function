# Testing

Install development dependencies and run `pytest`. The original tests still use their bundled fixtures and need no reference build. Network access is blocked for all tests. No test downloads reference files or builds a replacement database.

## Test a fully built human GRCh38 reference locally

From the package root, run `pytest`. The full-reference tests automatically use the newest installed release in the default cache (`homo_sapiens/GRCh38/release-N/ensembl.sqlite`). `FUSION_FUNCTION_CACHEDIR` changes that cache location, and `FUSION_FUNCTION_DB` selects a specific existing database. To run just these tests with automatic discovery:

```bash
pytest -m human_reference -v
```

The original integration tests in `test_fusion.py` and `test_ensembl.py` each run twice: `[bundled]` uses the original offline fixture and `[built]` uses the completed reference. Assertions are identical for both variants, including SLC45A2::AMACR racemase preservation and AMACR PANTHER subfamily metadata. Missing features in the built reference produce test failures. The byte-for-byte legacy output comparisons remain tied to the bundled reference.

To run SLC45A2::AMACR against both sources, or only the build:

```bash
pytest tests/test_fusion.py -k slc45a2 -v
pytest tests/test_fusion.py -k slc45a2 -m human_reference -v
```

To override the discovered database:

```bash
pytest -m human_reference --human-reference-db=/path/to/ensembl.sqlite -v
```

The selected database must be a completed, preprocessed reference. These tests open it in read-only mode and reject the small bundled regression fixture. The tests cover:

- Source and preprocessing counts, species/assembly/release metadata, and presence of primary chromosomes.
- Known GRCh38 chromosome lengths and assembly-aware matching for every DNA record.
- The real `HSCHR10_1_CTG2` cross-assembly name collision that failed the original build.
- Sequence chunk boundaries and reverse-complement retrieval.
- Canonical BCR, ABL1, and TP53 models, including both strands, source exon ranks, pre-mRNA, CDS/peptide lengths, and versioned IDs.
- InterPro metadata for their protein features.
- BCR–ABL1 fusion annotation using the supplied full reference.

The 55 built-reference cases (39 integration cases plus 16 additional checks) are skipped when no local reference is available. An explicitly configured missing database or an invalid selected database fails the tests instead of being silently skipped. No transcript data is mocked in the built variants. They exercise representative transcript and fusion paths; sequence-region and count checks cover the entire selected database. They do not claim exhaustive biological validation of every possible fusion.

To run all tests against a specific reference:

```bash
pytest --human-reference-db=/path/to/ensembl.sqlite
```

## Original and lightweight tests

`pytest -m 'not human_reference'` runs only the lightweight tests, regardless of whether a production reference is installed. `pytest -m 'integration and not human_reference'` selects just the original fusion regression scenarios. The fixture `tests/data/reference.sqlite` was converted from the user's cached API responses. It contains 15 transcript models and sparse genome chunks; its release label is for configuration tests, not a verified source release. It must not be used for production annotation or as the full-reference test input. `fusion_golden.json` contains 15 complete legacy outputs, including a same-transcript fusion, compared field-for-field with the local workflow.

Other tests cover reference configuration, invalid CLI options, both-strand preprocessing, rollback, progress totals and cleanup, the cross-assembly validation bug, and PANTHER source normalization and bulk subfamily matching. PANTHER checks include a missing-feature reproduction, family/protein mismatches, ambiguous assignments, offline reuse, provenance and rollback. A run that skips the full-reference tests does not validate a production database.

Build tests inject failures after DNA import, before preprocessing, and during
integrity checking. Retrying skips completed imports and preprocessing while
preserving the previous final reference until publication. Controlled FASTA
interrupts verify committed-record reuse and orphan-chunk cleanup. Additional
checks reject mismatched checkpoint settings and concurrent writers. PANTHER
14.1 tests cover the trailing-underscore filename and older column layout,
exact-release provenance, verified offline cache reuse, and propagation of
errors other than 404. Column-format checks cover unassigned and unnamed
records, malformed IDs, and keeping family descriptions and GO terms out of
subfamily names.

Progress tests also verify elapsed-time updates during a blocked SQLite call without progress callbacks, redirected-log heartbeats, named preprocessing phases, and timer cleanup on failures. Database calls remain on the caller's thread.

Single-gene integration cases verify retained-portion splicing and translation predictions alongside the explicit disruptive assumption. A C/N case verifies identical results to the corresponding N/C reconstruction. Lightweight checks cover missing inputs, invalid breakpoint errors, and `NotImplementedError` for N/N and C/C pairs before database access. Same-transcript N/C fixtures verify one initiation summary for the final product, including a UTR-only N-terminal portion followed by a retained C-terminal native start.

Late N-terminal breakpoint cases on both strands verify that all domains can remain included at 100% retention while the event is labelled assumed disruptive.

Single-gene UTR cases cover C/5' UTR, C/3' UTR, N/5' UTR, and N/3' UTR on both strands with the other transcript explicitly None. Start-loss and domain-cut cases also run against the installed reference. Controlled sequence fixtures verify native-start retention, next-ATG initiation after splicing, out-of-frame first starts, no remaining ATG, domains before or overlapping the selected start, and downstream premature termination.
