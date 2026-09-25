"""Score public teaching variants with Google DeepMind's AlphaGenome.

Queries the AlphaGenome API for variant effect predictions (splicing,
expression, chromatin accessibility) on a small panel of well-known,
publicly documented variants. Results are saved to results/ and summarized
for the README and figures.

Authentication uses a stored API key (never hard-coded); see README.
Coordinates below are GRCh38, verified against NCBI dbSNP.
"""
import sys

sys.path.insert(0, "/opt/hatch/skills/skill-creator/bin")
from dynamic_credentials import dynamic_credential_entry

import grpc
import json
import numpy as np

from alphagenome.data import genome
from alphagenome.models import dna_client
from alphagenome.models import variant_scorers as vs

# --- Auth (surrogate key injected by the credential broker; never printed) ---
entry = dynamic_credential_entry("custom.alphagenome", "access_token")
surrogate = entry["surrogate"]
with open("/run/hatch/cell-anchors/hatch-egress-ca.pem", "rb") as f:
    ca_cert = f.read()
creds = grpc.ssl_channel_credentials(root_certificates=ca_cert)
options = [("grpc.max_send_message_length", -1), ("grpc.max_receive_message_length", -1)]
channel = grpc.secure_channel("dns:///gdmscience.googleapis.com:443", creds, options=options)
grpc.channel_ready_future(channel).result(timeout=60)
model = dna_client.DnaClient(channel=channel, metadata=[("x-goog-api-key", surrogate)])
print("AlphaGenome client ready.", flush=True)

# --- Variant panel (GRCh38; all public teaching examples) ---
# Each variant gets scorers matched to its expected mechanism.
WIDTH = 131072  # 128 kb context window
VARIANTS = [
    {
        "rsid": "rs334",
        "chrom": "chr11", "pos": 5227002, "ref": "T", "alt": "A",
        "gene": "HBB", "consequence": "missense (Glu6Val)",
        "why": "The sickle-cell variant: a single base change that alters "
               "hemoglobin. Classic example of a coding variant with a "
               "large phenotypic effect.",
        "scorers": ["SPLICE_SITES", "SPLICE_SITE_USAGE", "RNA_SEQ"],
    },
    {
        "rsid": "rs4988235",
        "chrom": "chr2", "pos": 135851076, "ref": "G", "alt": "A",
        "gene": "MCM6 (regulates LCT)", "consequence": "intronic, regulatory",
        "why": "Lactase persistence: a noncoding variant in an enhancer that "
               "keeps the LCT gene switched on into adulthood. A textbook "
               "regulatory variant, ideal for chromatin and expression scorers.",
        "scorers": ["DNASE", "ATAC", "RNA_SEQ"],
    },
    {
        "rsid": "rs429358",
        "chrom": "chr19", "pos": 44908684, "ref": "T", "alt": "C",
        "gene": "APOE", "consequence": "missense (Cys130Arg)",
        "why": "Defines the APOE-e4 haplotype, the strongest common genetic "
               "risk factor for late-onset Alzheimer disease.",
        "scorers": ["SPLICE_SITES", "SPLICE_SITE_USAGE", "RNA_SEQ"],
    },
    {
        "rsid": "rs7412",
        "chrom": "chr19", "pos": 44908822, "ref": "C", "alt": "T",
        "gene": "APOE", "consequence": "missense (Arg176Cys)",
        "why": "Defines the APOE-e2 haplotype, which is protective against "
               "Alzheimer disease. Pairs with rs429358 as a contrast: two "
               "missense variants in one gene, opposite directions of risk.",
        "scorers": ["SPLICE_SITES", "SPLICE_SITE_USAGE", "RNA_SEQ"],
    },
    {
        "rsid": "rs4244285",
        "chrom": "chr10", "pos": 94781859, "ref": "G", "alt": "A",
        "gene": "CYP2C19", "consequence": "splice defect (CYP2C19*2)",
        "why": "The CYP2C19*2 allele: creates an aberrant splice site that "
               "knocks out enzyme activity and changes drug metabolism "
               "(e.g. clopidogrel response). Chosen to exercise the splicing "
               "scorers on a variant whose mechanism is splicing.",
        "scorers": ["SPLICE_SITES", "SPLICE_SITE_USAGE", "SPLICE_JUNCTIONS"],
    },
]


def summarize_anndata(ad, top_n=40):
    """Compact JSON-safe summary of one scorer's AnnData output."""
    out = {
        "scorer": str(ad.uns.get("variant_scorer", "")),
        "shape": list(ad.shape),
    }
    try:
        x = ad.X
        if hasattr(x, "toarray"):
            x = x.toarray()
        x = np.asarray(x, dtype=float)
        rows = (
            [str(v) for v in ad.obs.iloc[:, 0]]
            if len(ad.obs.columns)
            else [str(i) for i in range(ad.shape[0])]
        )
        cols = (
            [str(v) for v in ad.var.iloc[:, 0]]
            if len(ad.var.columns)
            else [str(i) for i in range(ad.shape[1])]
        )
        flat = []
        for i in range(min(ad.shape[0], 100)):
            for j in range(ad.shape[1]):
                v = float(x[i, j])
                if np.isfinite(v) and abs(v) > 1e-9:
                    flat.append(
                        {
                            "track": rows[i] if i < len(rows) else i,
                            "col": cols[j] if j < len(cols) else j,
                            "value": round(v, 6),
                        }
                    )
        flat.sort(key=lambda d: -abs(d["value"]))
        out["top_scores"] = flat[:top_n]
        out["n_nonzero"] = len(flat)
        out["max_abs_effect"] = round(abs(flat[0]["value"]), 6) if flat else 0.0
    except Exception as e:  # keep going; record the failure
        out["extract_error"] = f"{type(e).__name__}: {e}"
    return out


def main():
    results = {}
    for v in VARIANTS:
        label = v["rsid"]
        print(f"\n=== {label} {v['chrom']}:{v['pos']} {v['ref']}>{v['alt']} ({v['gene']}) ===",
              flush=True)
        interval = genome.Interval(
            chromosome=v["chrom"], start=v["pos"] - WIDTH // 2, end=v["pos"] + WIDTH // 2
        )
        variant = genome.Variant(
            chromosome=v["chrom"], position=v["pos"],
            reference_bases=v["ref"], alternate_bases=v["alt"],
        )
        scorers = [vs.RECOMMENDED_VARIANT_SCORERS[name] for name in v["scorers"]]
        try:
            outputs = model.score_variant(interval=interval, variant=variant,
                                          variant_scorers=scorers)
        except Exception as e:
            print(f"  ERROR: {type(e).__name__}: {e}", flush=True)
            results[label] = {"error": f"{type(e).__name__}: {e}", "variant": v}
            continue
        res = {"variant": v, "scores": []}
        for ad in outputs:
            s = summarize_anndata(ad)
            res["scores"].append(s)
            print(f"  {s['scorer']}: shape={s['shape']} "
                  f"max_abs_effect={s.get('max_abs_effect')}", flush=True)
        results[label] = res

    with open("results/variant_scores.json", "w") as f:
        json.dump(results, f, indent=1)
    print("\nWrote results/variant_scores.json", flush=True)


if __name__ == "__main__":
    main()
