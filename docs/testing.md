# Testing

Install development dependencies and run `pytest`. The original tests still use their bundled fixtures and need no reference build. Network access is blocked for all tests. No test downloads reference files or builds a replacement database.

## Test a fully built human GRCh38 reference locally

From the package root, run `pytest`. The full-reference tests automatically use the newest installed release in the default cache (`homo_sapiens/GRCh38/release-N/ensembl.sqlite`). `FUSION_FUNCTION_CACHEDIR` changes that cache location, and `FUSION_FUNCTION_DB` selects a specific existing database. To run just these tests with automatic discovery:

```bash
pytest -m human_reference -v
```

The original integration tests in `test_fusion.py` and `test_ensembl.py` each run twice: `[bundled]` uses the original offline fixture and `[built]` uses the completed reference. Assertions are identical for both variants, including retention of the AMACR family match in SLC45A2::AMACR and AMACR PANTHER subfamily metadata. Missing features in the built reference produce test failures. Complete output snapshots remain tied to the bundled reference; they were updated for the expanded feature types, exact source boundaries and provenance fields. Existing frame and translation-start predictions are unchanged.

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

The 63 built-reference cases are skipped when no local reference is available. An explicitly configured missing database or an invalid selected database fails the tests instead of being silently skipped. No transcript data is mocked in the built variants. They exercise representative transcript and fusion paths; sequence-region and count checks cover the entire selected database. They do not claim exhaustive biological validation of every possible fusion.

To run all tests against a specific reference:

```bash
pytest --human-reference-db=/path/to/ensembl.sqlite
```

## Original and lightweight tests

`pytest -m 'not human_reference'` runs only the lightweight tests, regardless of whether a production reference is installed. `pytest -m 'integration and not human_reference'` selects just the original fusion regression scenarios. The fixture `tests/data/reference.sqlite` was converted from the user's cached API responses. It contains 15 transcript models and sparse genome chunks; its release label is for configuration tests, not a verified source release. It must not be used for production annotation or as the full-reference test input. `fusion_golden.json` contains 15 complete fusion regression outputs, including a same-transcript fusion, compared field-for-field with the local workflow.

Other tests cover reference configuration, invalid CLI options, both-strand preprocessing, rollback, progress totals and cleanup, the cross-assembly validation bug, and PANTHER source normalization and bulk subfamily matching. PANTHER checks include a missing-feature reproduction, family/protein mismatches, ambiguous assignments, offline reuse, provenance and rollback. A run that skips the full-reference tests does not validate a production database.

Build tests inject failures after DNA import, before preprocessing, and during integrity checking. Retrying skips completed imports and preprocessing while preserving the previous final reference until publication. Controlled FASTA interrupts verify committed-record reuse and orphan-chunk cleanup. Additional checks reject mismatched checkpoint settings and concurrent writers. PANTHER 14.1 tests cover the trailing-underscore filename and older column layout, exact-release provenance, verified offline cache reuse, and propagation of errors other than 404. Column-format checks cover unassigned and unnamed records, malformed IDs, and keeping family descriptions and GO terms out of subfamily names.

Progress tests also verify start/outcome logs without SQLite timer bars, advance notices for potentially long operations, named preprocessing phases, and progress-handler cleanup on failures. Database calls remain on the caller's thread.

Single-gene integration cases verify retained-portion splicing and translation predictions alongside the explicit disruptive assumption. A C/N case verifies identical results to the corresponding N/C reconstruction. Lightweight checks cover missing inputs, invalid breakpoint errors, and `NotImplementedError` for N/N and C/C pairs before database access. Same-transcript N/C fixtures verify one initiation summary for the final product, including a UTR-only N-terminal portion followed by a retained C-terminal native start.

Late N-terminal breakpoint cases on both strands verify that all domains can remain included at 100% retention while the event is labelled assumed disruptive.

Single-gene UTR cases cover C/5' UTR, C/3' UTR, N/5' UTR, and N/3' UTR on both strands with the other transcript explicitly None. Start-loss and domain-cut cases also run against the installed reference. Controlled sequence fixtures verify native-start retention, next-ATG initiation after splicing, out-of-frame first starts, no remaining ATG, domains before or overlapping the selected start, and downstream premature termination.

UniProt tests run offline using captured human AMACR XML and independently downloaded canonical/isoform protein sequences. They verify sequence identity, described isoform reconstruction, omission of changed/fuzzy intervals, evidence, accession mapping, genomic coordinates, retention of catalytic versus substrate-binding sites, offline reuse and rollback. Synthetic tests exercise the CLI's default download and local-file override. The two built-reference UniProt checks verify actual AMACR annotations and their status in the fusion product. A selected database missing these annotations fails; these checks skip only when no built reference is available. Lightweight regression tests verify that missing annotations fail instead of silently skipping. No full human reference is constructed or downloaded by pytest.

The AMACR fixtures were captured from `https://rest.uniprot.org/uniprotkb/Q9UHK6.xml`, its canonical sequence, and `https://rest.uniprot.org/uniprotkb/Q9UHK6-5.fasta` on 2026-10-01. They are small source captures, not a built human reference.

Repeated-preparation tests verify UniProt lookup reuse and invalidation when the XML, mapping implementation, Ensembl source checksums, release or lookup table changes. References without source checksums are reverified. Error-summary tests check grouping, sample limits and read-only inspection of existing references without a saved summary.

Partial-codon regressions cover both strands, missing one or two leading bases, presence/absence of a terminal stop codon, and single N- or C-terminal partners. They independently check per-base genomic feature mapping, incomplete first residues, alternative initiation and preserved downstream features. Built-reference tests check the three reported release-100/116 transcripts; older prepared models that still store their length errors fail these checks.

Structure-filter tests verify that AlphaFold/SIFTS mappings are omitted whether their intervals fit the peptide or exceed it, including differing pipeline and database labels. Actual functional features remain unchanged, and invalid domain coordinates still reject the transcript. Three built-reference checks cover the reported AlphaFold failures; older prepared models storing these errors fail.

Compact-storage tests compare all 15 complete fusion regression outputs against compressed transcript models, check ordinary row storage and readable JSON errors, and reject unknown codecs. Existing plain-JSON fixtures remain readable. Rollback checks exercise the new storage while preserving the prior reference and the caller's secure-delete, journal and synchronous settings.
