# Fusion-function

![Coverage](https://raw.githubusercontent.com/creisle/fusion_function/badges/badges/main/coverage.svg)

Predict fusion reading frames and retained, disrupted, or excluded functional protein features from human GRCh38 transcript breakpoints. Predictions include effects of splicing, translation initiation, and premature termination. Annotation uses a preprocessed local SQLite reference.

## Quick start

```bash
pip install fusion-function
fusion-function prepare-data
```

Preparation downloads and builds the latest Ensembl reference in the default user cache. It can take an hour or longer and needs several GB of disk space.

Use the annotate_fusion_domains function to get information about the status of various domains in the expected fusion product

```python
from fusion_function import ReferenceDatabase, annotate_fusion_domains

with ReferenceDatabase() as ref:
    # ex. BCR::ABL1
    result = annotate_fusion_domains(
        transcript1_id="ENST00000305877", # BCR
        transcript2_id="ENST00000318560", # ABL1
        breakpoint1="22:23290413",
        breakpoint2="9:130854064",
        gene1_terminus="N",
        gene2_terminus="C",
        reference=ref
    )
```

Thisw will return an object with the following shape. See the [api](./docs/api.md) for details.

```json
{
    "frame_status": "in_frame",
    "domains": [
        {
            "transcript_id": "ENST00000305877",
            "interpro_id": "IPR036481",
            "name": "Bcr-Abl oncoprotein oligomerisation domain superfamily",
            "domain_type": "homologous_superfamily",
            "start": 1,
            "end": 67,
            "sources": ["SuperFamily"],
            "feature_ids": ["SSF69036"],
            "breakpoint_based_status": "included",
            "breakpoint_retained_percent": 100.0,
            "post_splicing_status": "preserved",
            "post_translation_status": "preserved"
        },
        ...
    ],
    "translation_start": {"ENST00000305877": "native_start_retained"}
}
```

Analysis uses the newest prepared local release by default.

## Documentation

- [API](docs/api.md)
- [Reference configuration and preparation](docs/reference.md)
- [Methodology](docs/methodology.md)
- [Testing](docs/testing.md)

## Important Limitations

- This package uses the splicing model defined by [MAVIS](https://github.com/bcgsc/mavis), this is non-exhaustive. It assumes splice sites to be disrupted based on a breakpoint being within 2bp but there are many other ways to disrupt splicing that are difficult to predict computationally (ex. deep intronic). This package only covers the standard scenarios
- Currently we only support Hg38
- Only exact breakpoints are supported
- Predicted structural consequences do not establish fusion expression, oncogenicity, pathogenicity, or clinical actionability.

## Citation

This package ports and extends selected fusion-annotation logic from [MAVIS](https://github.com/bcgsc/mavis). Please cite:

Reisle C, Mungall KL, Choo C, et al. MAVIS: merging, annotation, validation, and illustration of structural variants. *Bioinformatics*. 2019;35(3):515–517. PMID:30016509.
