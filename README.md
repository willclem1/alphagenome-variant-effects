# Can an AI read a DNA typo? Validating AlphaGenome on variants we already understand

AlphaGenome (Google DeepMind) is a deep learning model that predicts what a
single-letter DNA change does: does it break a splice site, change how much
RNA a gene makes, or open up chromatin? The claims are impressive. This repo
asks the obvious question: **does it get the easy ones right?**

We scored five famous, publicly documented variants, each chosen because its
real biology is already settled, and checked whether the model's predictions
match. On the cases with known answers, it does.

All coordinates are GRCh38, verified against NCBI dbSNP. No private or
personal data is used anywhere in this repo.

## The idea in 30 seconds

Give the model a variant and a window of DNA. It predicts the difference
between the reference and alternate allele across molecular readouts: splice
site strength and usage, splice junction counts, RNA expression,
polyadenylation, and chromatin accessibility, across hundreds of human cell
types and tissues. We query the hosted AlphaGenome API and summarize the
results.

## The panel: five variants with known answers

| Variant | Gene | Type | Known biology | What a correct prediction looks like |
|---|---|---|---|---|
| rs4244285 (chr10:94,781,859 G>A) | CYP2C19 | splice defect (*2 allele) | Creates an aberrant splice site; knocks out enzyme activity, alters drug metabolism | Strong splicing signal |
| rs4988235 (chr2:135,851,076 G>A) | MCM6 / LCT | intronic, regulatory | Lactase persistence: a noncoding enhancer variant keeping LCT expressed into adulthood | Strong chromatin signal, modest expression shift |
| rs334 (chr11:5,227,002 T>A) | HBB | missense (Glu6Val) | Sickle cell: a protein-level change with a huge phenotype | Quiet on splicing/expression |
| rs429358 (chr19:44,908,684 T>C) | APOE | missense | Defines the APOE-e4 haplotype, the strongest common risk factor for late-onset Alzheimer disease | Quiet (protein-level effect) |
| rs7412 (chr19:44,908,822 C>T) | APOE | missense | Defines the protective APOE-e2 haplotype | Quiet (protein-level effect) |

The APOE pair is the control experiment: two missense variants, one gene,
opposite risk directions, both acting at the protein level.

## What the model predicted

Strongest predicted effect per variant and modality (log scale):

![Effect sizes across modalities](figures/effects_by_modality.png)

| Variant | Splicing | Expression | Chromatin |
|---|---|---|---|
| rs4244285 (CYP2C19) | **4.75** (splice junctions) | not tested | not tested |
| rs4988235 (MCM6) | not tested | 0.10 (RNA-seq) | **1.69** (ATAC-seq) |
| rs334 (HBB) | 0.05 | 0.09 | not tested |
| rs429358 (APOE) | 0.03 | 0.04 | not tested |
| rs7412 (APOE) | 0.02 | 0.03 | not tested |

## Scorecard: did it get the easy ones right?

- **CYP2C19 \*2: yes.** A new acceptor site (score 0.84) with junction log
  fold change 4.75, strongest in liver tracks, matching the established
  splice-defect mechanism. The model found the broken splice site.
- **Lactase persistence: yes.** The strongest chromatin effects in the panel
  (ATAC log2 fold change up to 1.69, top tissues including small intestine),
  exactly what an enhancer variant should show. Expression effects are
  modest, consistent with a subtle regulatory tweak rather than a knockout.
- **The three missense controls: correctly null.** rs334, rs429358, and
  rs7412 are near-zero across these modalities, which is the honest answer:
  their effects play out at the protein level, which splicing and expression
  scorers are not designed to capture. A null prediction here is
  informative, not a failure, and a useful reminder that no single modality
  tells the whole story.

## Run it yourself

```bash
pip install -r requirements.txt
python score_variants.py    # queries the API, writes results/variant_scores.json
python make_figure.py       # rebuilds figures/effects_by_modality.png
```

`score_variants.py` needs an AlphaGenome API key. It reads the key from a
local credential store at runtime (see the top of the script); the key is
never hard-coded or committed. Raw per-track scores for every scorer live in
`results/variant_scores.json`.

## Go further

The same script scales: point the variant panel at any candidate list and
rank by strongest predicted effect. That is how this code is used in
practice, as a triage step that decides which variants deserve experimental
follow-up.

## Files

- `score_variants.py`: variant panel, API queries, JSON summaries
- `make_figure.py`: rebuilds the summary figure
- `results/variant_scores.json`: full per-track scores
- `figures/effects_by_modality.png`: summary plot
