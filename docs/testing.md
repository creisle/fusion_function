# Testing

From the repository root:

```bash
pip install -e . --group dev
pytest
```

Tests run offline using bundled fixtures and synthetic data. They never download or build a full reference.

If a prepared human GRCh38 database is installed, `pytest` also runs the reference integration tests against it. The newest local release is selected automatically; `FUSION_FUNCTION_CACHEDIR` changes the cache location and `FUSION_FUNCTION_DB` selects a database. These tests skip when no reference is available. An explicitly selected missing or invalid database fails.

```bash
# Lightweight tests only
pytest -m 'not human_reference'

# Reference integration tests only
pytest -m human_reference -v

# Run all tests against a specific reference
pytest --human-reference-db=/path/to/ensembl.sqlite
```

Reference tests accept full source-built databases and compact prebuilt references, and open them read-only. Bundled cases remain available without either.

The suite covers fusion frames and protein features, translation initiation, preprocessing, configuration, rollback, and prebuilt export/download behavior. Passing lightweight tests does not validate a production reference.
