"""Bar chart of max predicted effect per variant and scoring modality."""
import json
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

with open("results/variant_scores.json") as f:
    results = json.load(f)


def scorer_name(s):
    """Extract a clean modality name from a scorer repr string."""
    m = re.search(r"requested_output=([A-Z_]+)", s)
    if m:
        return m.group(1)
    m = re.search(r"(\w+Scorer)\(", s)
    name = m.group(1) if m else s
    return {"SpliceJunctionScorer": "SPLICE_JUNCTIONS"}.get(name, name)

ORDER = ["rs334", "rs4988235", "rs429358", "rs7412", "rs4244285"]
LABELS = {
    "rs334": "rs334\nHBB missense",
    "rs4988235": "rs4988235\nMCM6 regulatory",
    "rs429358": "rs429358\nAPOE missense",
    "rs7412": "rs7412\nAPOE missense",
    "rs4244285": "rs4244285\nCYP2C19 splice",
}
SCORER_LABELS = {
    "SPLICE_SITES": "Splice sites",
    "SPLICE_SITE_USAGE": "Splice site usage",
    "SPLICE_JUNCTIONS": "Splice junctions",
    "RNA_SEQ": "Expression (RNA-seq)",
    "DNASE": "Chromatin (DNase)",
    "ATAC": "Chromatin (ATAC)",
}

scorers = []
for rs in ORDER:
    for s in results[rs]["scores"]:
        name = scorer_name(s["scorer"])
        if name not in scorers:
            scorers.append(name)

vals = {rs: {} for rs in ORDER}
for rs in ORDER:
    for s in results[rs]["scores"]:
        name = scorer_name(s["scorer"])
        vals[rs][name] = s.get("max_abs_effect", 0.0)

x = np.arange(len(ORDER))
width = 0.8 / len(scorers)
fig, ax = plt.subplots(figsize=(11, 6))
for i, sc in enumerate(scorers):
    v = [vals[rs].get(sc, 0.0) for rs in ORDER]
    ax.bar(x + (i - (len(scorers) - 1) / 2) * width, v, width,
           label=SCORER_LABELS.get(sc, sc))

ax.set_yscale("log")
ax.set_xticks(x)
ax.set_xticklabels([LABELS[rs] for rs in ORDER])
ax.set_ylabel("Max predicted effect size (log scale)")
ax.set_title("AlphaGenome: strongest predicted effect per variant and modality")
ax.legend(frameon=False, fontsize=9)
fig.tight_layout()
fig.savefig("figures/effects_by_modality.png", dpi=150)
print("Wrote figures/effects_by_modality.png")
