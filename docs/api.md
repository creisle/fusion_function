# API Reference

## `annotate_fusion_domains`

```python
from fusion_function.fusion import annotate_fusion_domains
```

```python
annotate_fusion_domains(
    transcript1_id: str | None = None,
    transcript2_id: str | None = None,
    breakpoint1: str | None = None,
    breakpoint2: str | None = None,
    gene1_terminus: FusionTerminus | None = None,
    gene2_terminus: FusionTerminus | None = None,
    inserted_sequence: str = "",
    *,
    reference: ReferenceDatabase | None = None,
    database: str | Path | None = None,
    release: int | None = None,
) -> FusionDomainResult | EnsemblError
```

Annotates retained functional features for one or two Ensembl transcripts. A product is reconstructed when two partners provide one N-terminal and one C-terminal portion.

## Parameters

### `transcript1_id`, `transcript2_id`

Ensembl transcript IDs for the two fusion partners.

For a single known gene, omit the other transcript or pass `None`. The known
partner still requires its breakpoint and retained terminus. Either partner
slot can hold the known transcript.

Example:

```python
"ENST00000305877"
```

### `breakpoint1`, `breakpoint2`

Exact genomic breakpoints in `chromosome:position` format.

Example:

```python
"22:23290413"
```

Breakpoint intervals are not supported.

Coordinates must use the same genome assembly as the Ensembl transcript annotations used by the package. The current workflow assumes GRCh38.

### `gene1_terminus`, `gene2_terminus`

The retained terminus of each fusion partner:

```python
"N"
"C"
```

An N/C or C/N pair supports product reconstruction. N/N and C/C pairs raise
`NotImplementedError`. A single known N- or C-terminal partner remains valid
and is assumed disruptive. Its retained portion is spliced and checked.

Both single-partner events and N/C products use a surviving native ATG.
If no native start survives and some annotated coding sequence remains, the
first ATG in the complete spliced sequence is used, in any phase. This includes
ATGs created across a reconstructed fusion junction. Features are checked for
start-related loss, compatibility with their original coding frame, and
premature termination.

The `translation_start` summary reports when an alternative start is selected.
This is an ATG-based initiation heuristic applied to the known sequence. For a
single-partner event, the unknown partner and its junction are not reconstructed.

For example, to annotate one known partner:

```python
result = annotate_fusion_domains(
    transcript1_id="ENST00000296589",
    breakpoint1="5:33973126",
    gene1_terminus="N",
    reference=reference,
)
```

### `reference`, `database`, and `release`

Reference data are read locally. Reuse a `ReferenceDatabase` through `reference=`, or choose a file through `database=` and optionally validate its `release=`. A supplied reader cannot be combined with the other reference options. Without these arguments, the newest prepared local release is discovered under the default cache; `FUSION_FUNCTION_DB` overrides discovery. The removed `session=` parameter is no longer accepted.

See [`reference.md`](reference.md).

### `inserted_sequence`

Optional nucleotide sequence inserted between the two fusion breakpoints.

Default:

```python
""
```

Inserted sequence is included during fusion reconstruction and can affect reading frame and premature-stop predictions.

## Return value

A successful result contains:

```python
{
    "frame_status": "in_frame",
    "domains": [...]
}
```

### `frame_status`

Possible values:

| Value | Meaning |
|---|---|
| `in_frame` | All predicted splice products preserve source coding phase |
| `out_of_frame` | All predicted splice products disrupt source coding phase |
| `in_frame/out_of_frame` | Different predicted splice products produce different frame states |
| `None` | No translated original CDS remains, or no start was found in an isolated retained portion |

Single-partner results also contain `assumed_disruptive: True`.
This labels the event assumption independently of domain predictions. Their
frame summary describes the known retained portion. Domains
still report the retained, partially retained, and excluded features of each
known partner, together with their splicing and sequence checks.

### `translation_start`

Results include one initiation summary keyed by the transcript supplying the
selected native start, or the intended initiating transcript if its start is
lost or an alternative start is used. Initiation is assessed in the complete
spliced fusion product. Using the same transcript for both partners does not
produce duplicate statuses:

```python
"translation_start": {
    "ENST...": "alternative_start_not_found",
}
```

Each value uses the following controlled vocabulary:

| Value | Meaning |
|---|---|
| `native_start_retained` | A native ATG survives and supplies initiation |
| `alternative_start_found` | The native start is lost; an ATG search selected an alternative start |
| `alternative_start_not_found` | The native start is lost; an ATG search found no alternative start |
| `native_start_lost` | The native start is lost and no annotated coding sequence remains, so no alternative ATG search was performed |

Distinct outcomes from alternative splice products are deduplicated and joined
with `/`, in the table order. For example:

```python
"translation_start": {
    "ENST...": "native_start_retained/alternative_start_found",
}
```

Separate transcript keys can occur if alternative products initiate from
different transcripts; each product has only one selected initiation site.

These outcomes describe the retained, spliced sequence. They do not predict
initiation supplied by an unknown partner. A native ATG cut in half is lost.
For an N/C product whose N-terminal portion is UTR-only, a surviving native
start in the C-terminal portion supplies initiation.

### `domains`

Each functional feature contains:

```python
{
    "transcript_id": str,
    "interpro_id": str | None,
    "name": str,
    "domain_type": str,
    "start": int,
    "end": int,
    "breakpoint_based_status": str,
    "breakpoint_retained_percent": float,
    "splicing_based_status": str | None,
    "fusion_sequence_status": str | None,
}
```

#### `domain_type`

Returned functional feature types include:

- `domain`
- `binding_site`
- `active_site`
- `conserved_site`

#### `breakpoint_based_status`

| Value | Meaning |
|---|---|
| `included` | Entire feature is retained by the genomic breakpoint |
| `disrupted` | Breakpoint passes through the feature |
| `excluded` | Feature is completely outside the retained portion |

Overlapping annotations representing the same functional feature may produce a slash-delimited combined status.

#### `breakpoint_retained_percent`

Percentage of the annotated feature retained directly by the breakpoint.

This reflects genomic retention before considering fusion splicing or translation.

#### `splicing_based_status`

| Value | Meaning |
|---|---|
| `preserved` | The portion retained by the breakpoint survives predicted splicing |
| `lost` | Predicted splicing removes part or all of that retained portion |
| `None` | Downstream splicing status is not applicable or not evaluated |

Multiple possible predictions may be slash-delimited.

These checks apply to both `included` and `disrupted` features. For a disrupted
feature, downstream `preserved` describes only the retained portion; the original
domain remains disrupted. Completely excluded features have no downstream status.

#### `fusion_sequence_status`

| Value | Meaning |
|---|---|
| `preserved` | Retained portion remains in the expected coding frame without premature termination |
| `frame_disrupted` | Feature is retained but translated in an incompatible frame |
| `translation_start_disrupted` | The selected start falls within the retained feature, omitting its beginning |
| `translation_start_excluded` | The feature lies before the selected start, or no ATG is found by the C-terminal heuristic |
| `premature_termination_disrupted` | Translation terminates within the retained portion |
| `premature_termination_excluded` | Translation terminates before the retained portion |
| `None` | Sequence-level status is not applicable or not evaluated |

Multiple possible predictions may be slash-delimited.

## Errors

Transcript retrieval failures and invalid partner inputs return an
`EnsemblError` dictionary containing an `error` message. This includes:

- no known transcript
- a known partner without a breakpoint or N/C terminus
- breakpoint intervals
- breakpoint chromosome/transcript chromosome mismatch
- breakpoints outside the transcript span

Reference configuration failures and missing required preprocessing metadata
still raise exceptions. Two known partners with N/N or C/C termini raise
`NotImplementedError`. N/None and C/None remain valid and are annotated as
assumed disruptive; metadata on the missing partner is ignored.
