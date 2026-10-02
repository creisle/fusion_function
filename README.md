# fusion_function

![Coverage](https://raw.githubusercontent.com/creisle/fusion_function/badges/badges/main/coverage.svg)

Predict fusion reading frames and retained, disrupted, or excluded functional protein features from human GRCh38 transcript breakpoints. Predictions include effects of splicing, translation initiation, and premature termination. Annotation uses a preprocessed local SQLite reference.

## Quick start

```bash
pip install fusion-function
fusion-function prepare-data
```

Preparation downloads and builds the latest Ensembl reference in the default user cache. It can take an hour or longer and needs several GB of disk space.

```python
from fusion_function import annotate_fusion_domains

result = annotate_fusion_domains(
    transcript1_id="ENST00000305877",
    transcript2_id="ENST00000318560",
    breakpoint1="22:23290413",
    breakpoint2="9:130854064",
    gene1_terminus="N",
    gene2_terminus="C",
)
```

Analysis uses the newest prepared local release by default.

## Documentation

- [API](docs/api.md)
- [Reference configuration and preparation](docs/reference.md)
- [Methodology](docs/methodology.md)
- [Testing](docs/testing.md)

## Citation

This package ports and extends selected fusion-annotation logic from [MAVIS](https://github.com/bcgsc/mavis). Please cite:

Reisle C, Mungall KL, Choo C, et al. MAVIS: merging, annotation, validation, and illustration of structural variants. *Bioinformatics*. 2019;35(3):515–517. PMID:30016509.

Predicted structural consequences do not establish fusion expression, oncogenicity, pathogenicity, or clinical actionability.
