"""Check that the inputs are present, and say exactly where each one comes from if it is not.

The models are large (about 50MB of SBML) and freely redistributable from BiGG, so they are fetched
rather than committed. The three experimental tables are publisher supplementary files behind pages that
do not serve them to scripts, so they are downloaded by hand and committed, which is also what makes the
analysis reproducible without a live internet connection.

Usage:  PYTHONPATH=src python scripts/01_fetch.py [--download]
"""

from __future__ import annotations

import argparse
import sys
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MODELS = ROOT / "data" / "external" / "models"
RAW = ROOT / "data" / "raw"

# Every genome-scale reconstruction of E. coli K-12 MG1655 in the BiGG series, oldest to newest, plus
# the published correction of the newest one.
MODEL_URLS = {
    "iJR904.xml": "http://bigg.ucsd.edu/static/models/iJR904.xml",
    "iAF1260.xml": "http://bigg.ucsd.edu/static/models/iAF1260.xml",
    "iJO1366.xml": "http://bigg.ucsd.edu/static/models/iJO1366.xml",
    "iML1515.xml": "http://bigg.ucsd.edu/static/models/iML1515.xml",
}

# Not on BiGG: the corrected model is a supplementary file of the paper that made the corrections.
BY_HAND = {
    MODELS / "iML1515_model_adjusted_all_corrections.xml": (
        "Bernstein et al. (2023), 'Evaluating E. coli genome-scale metabolic model accuracy with "
        "high-throughput mutant fitness data', Molecular Systems Biology 19:e11566.\n"
        "      Supplementary model file (Dataset EV / model adjusted, all corrections)."
    ),
    RAW / "keio_baba2006_table3.xls": (
        "Baba et al. (2006), 'Construction of Escherichia coli K-12 in-frame, single-gene knockout "
        "mutants: the Keio collection', Molecular Systems Biology 2:2006.0008.\n"
        "      Supplementary Table 3 (44320_2006_BFMSB4100050_MOESM5_ESM.xls)."
    ),
    RAW / "tradis_goodall2018_tableS4.xlsx": (
        "Goodall et al. (2018), 'The Essential Genome of Escherichia coli K-12', mBio 9:e02096-17.\n"
        "      Supplementary Table S4 (mbo001183726st4.xlsx)."
    ),
    RAW / "pec_genes.txt.rtf": (
        "Profiling of the E. coli Chromosome (PEC) database, gene table export.\n"
        "      https://shigen.nig.ac.jp/ecoli/pec/"
    ),
}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--download", action="store_true",
                        help="fetch the four BiGG models (about 37MB)")
    args = parser.parse_args()

    MODELS.mkdir(parents=True, exist_ok=True)
    missing = []

    for name, url in MODEL_URLS.items():
        target = MODELS / name
        if target.exists():
            print(f"  present   {name}")
        elif args.download:
            print(f"  fetching  {name} from {url}", flush=True)
            urllib.request.urlopen(url)
            urllib.request.urlretrieve(url, target)
            print(f"  fetched   {name}")
        else:
            missing.append((target, f"BiGG Models: {url}\n      or re-run this script with --download"))

    for target, source in BY_HAND.items():
        if target.exists():
            print(f"  present   {target.name}")
        else:
            missing.append((target, source))

    if missing:
        print(f"\n{len(missing)} input(s) missing. Each has to be downloaded to the path shown:\n")
        for target, source in missing:
            print(f"  {target.relative_to(ROOT)}")
            print(f"      {source}\n")
        return 1
    print("\nall inputs present")
    return 0


if __name__ == "__main__":
    sys.exit(main())
