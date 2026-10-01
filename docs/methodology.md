# Methodology

`fusion_function` predicts how genomic fusion breakpoints affect transcript structure, coding frame, and annotated functional protein features.

The analysis is transcript-specific and depends on the supplied Ensembl transcript IDs.

## 1. Transcript annotation

For each fusion partner, transcript structure and protein-feature annotations are read from preprocessed Ensembl FTP data in a local SQLite database.

The annotation includes information such as:

- genomic transcript coordinates
- strand
- exon coordinates
- coding-sequence blocks
- transcript-oriented pre-mRNA sequence
- protein features

Reference preprocessing resolves available InterPro links and metadata before runtime.

## 2. Breakpoint orientation

The supplied fusion terminus determines which side of each transcript is retained.

For each partner:

- `N` retains the transcript-oriented sequence from the transcript start through the breakpoint
- `C` retains the transcript-oriented sequence from the breakpoint through the transcript end

Transcript strand is accounted for when genomic coordinates are converted into transcript-oriented coordinates.

A valid fusion requires one N-terminal and one C-terminal partner.

## 3. Fusion layout

The retained N-terminal partner is placed before the retained C-terminal partner.

Any `inserted_sequence` is placed between them.

The result is an unspliced fusion pre-mRNA representation containing retained sequence, coding segments, and splice sites from both partners.

## 4. Splice-site prediction

Canonical donor and acceptor positions and their disruption windows are precomputed for every transcript. Runtime retains the sites on the selected fusion side and checks whether the breakpoint disrupts them.

A splice site is considered unavailable if the breakpoint directly disrupts the genomic region surrounding that site.

Plausible donor/acceptor combinations are used to generate candidate fusion splice products.

This is a structural prediction and does not model tissue-specific splice abundance or experimentally measured fusion RNA.

## 5. Coding-frame prediction

For each predicted splice product, retained CDS fragments are mapped into fusion-product coordinates.

The first retained native ATG in the spliced product supplies initiation. This
normally comes from the N-terminal partner; a surviving C-terminal native start
can supply initiation when the N-terminal portion contains only UTR sequence.
Each product has one selected initiation site, including when both breakpoints
use the same transcript. Distinct initiation statuses from alternative splice
products are deduplicated and joined with `/`, as for other status fields.

Coding phase is compared with the original source CDS to determine whether each product is:

- `in_frame`
- `out_of_frame`

If no native ATG survives and the product retains annotated coding sequence,
translation is modelled from its first ATG. This applies to N/C fusions and to
single N- or C-terminal portions. Searching the complete spliced sequence also
includes ATGs formed across the fusion junction. The search does not skip an
earlier out-of-frame ATG in favour of a later in-frame one.

The overall fusion result summarizes all predicted products.

## 6. Premature termination

The reconstructed spliced nucleotide sequence is scanned in-frame from the predicted translation start.

The first in-frame stop codon is recorded and used when classifying downstream functional features.

## 7. Functional feature annotation

Protein features are evaluated at three levels.

### Breakpoint-based status

This considers only whether the original feature coordinates fall within the retained portion of the fusion partner:

- `included`
- `disrupted`
- `excluded`

`breakpoint_retained_percent` reports how much of the feature remains.

### Splicing-based status

For both fully and partially retained features, the portion retained by the
breakpoint is mapped through each predicted splice product. Splicing removes a
feature if it removes any of that retained portion. Bases already excluded by
the breakpoint do not count as a further splicing loss.

A feature can be:

- `preserved`
- `lost`

### Fusion-sequence status

For features preserved through splicing, the reconstructed translation is assessed.

Frame comparisons use the original CDS position of the retained portion,
including when a breakpoint cuts within a codon. A downstream `preserved` status
does not restore a disrupted domain: it describes the surviving portion only.

Possible classifications include:

- `preserved`
- `frame_disrupted`
- `translation_start_disrupted`
- `translation_start_excluded`
- `premature_termination_disrupted`
- `premature_termination_excluded`

For a single known partner, the retained portion is modelled separately.
N/None and C/None are supported; two known partners with N/N or C/C termini
raise `NotImplementedError`. The normal splice-site rules apply to the known sequence.
Both N- and C-terminal portions use a native ATG when present after splicing;
otherwise, if coding sequence remains, translation is modelled from the first
remaining ATG under the same rules as a two-partner fusion.
Features before or overlapping the selected start receive the corresponding
translation-start exclusion or disruption status. The event assumption remains
`assumed_disruptive: True`, independent of the individual feature predictions.

## 8. InterPro annotation

Member-database features are linked to InterPro entries during reference preprocessing. No runtime annotation API calls are made.

Canonical InterPro names are preferred for real functional features such as domains, binding sites, active sites, and conserved sites.

Overlapping annotations representing the same functional feature are collapsed to reduce redundant output.

Protein- or family-specific annotations may be retained as fallback domain annotations when no corresponding functional InterPro feature represents that region.

## Interpretation

The output describes predicted structural consequences of a proposed fusion.

It does not establish that the fusion:

- is expressed
- is oncogenic
- is pathogenic
- is clinically actionable
- occurs in a specific tumour clone
- produces the predicted splice product in vivo

Interpretation should therefore consider transcript choice, breakpoint accuracy, RNA evidence, tumour context, and experimental or clinical evidence separately.
