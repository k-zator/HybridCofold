import argparse
import csv
import json
import re
from pathlib import Path

# Regular expressions to parse GNINA .dat files and DynamicBind filenames
AFF_RE = re.compile(r"^Affinity:\s+(-?\d+(?:\.\d+)?)")
CNN_AFF_RE = re.compile(r"^CNNaffinity:\s+(-?\d+(?:\.\d+)?)")
LDDT_RE = re.compile(r"lddt(?P<lddt>\d+(?:\.\d+)?)")
DBSCORE_RE = re.compile(r"affinity(?P<aff>-?\d+(?:\.\d+)?)")

def parse_gnina_dat(dat_path: Path) -> dict:
    """Parse GNINA .dat file with regular expression to extract Affinity and CNNaffinity."""
    affinity = None
    cnn_affinity = None

    for line in dat_path.read_text().splitlines():
        m = AFF_RE.match(line.strip())
        if m:
            affinity = float(m.group(1))
            continue
        m = CNN_AFF_RE.match(line.strip())
        if m:
            cnn_affinity = float(m.group(1))
            continue

    if affinity is None or cnn_affinity is None:
        raise ValueError(f"Could not parse Affinity/CNNaffinity from: {dat_path}")

    return {"gnina_affinity": affinity, "gnina_cnnaffinity": cnn_affinity}

def parse_dynamicbind_from_meta_or_filename(meta_path: Path | None, search_dir: Path) -> dict:
    """Using JSON metadata from select_best_dynamicbind_outputs.py, read in the DynamicBind score and LDDT.
       If no metadata is provided, infer from best filenames inside search_dir."""

    if meta_path and meta_path.exists():
        meta = json.loads(meta_path.read_text())
        return {"dynamicbind_score": float(meta["dynamicbind_score"]), "lddt": float(meta["lddt"])}

    # If no JSON, infer from best filenames inside search_dir
    best_lddt = None
    best_aff = None
    for p in search_dir.rglob("*rank*_complex_lddt*_affinity*_relaxed.pdb"):
        m1 = LDDT_RE.search(p.name)
        m2 = DBSCORE_RE.search(p.name)
        if not (m1 and m2):
            continue
        lddt = float(m1.group("lddt"))
        aff = float(m2.group("aff"))
        key = (lddt, aff)
        if best_lddt is None or key > (best_lddt, best_aff):
            best_lddt, best_aff = lddt, aff

    if best_lddt is None:
        raise FileNotFoundError(f"Could not infer DynamicBind scores from filenames under: {search_dir}")

    return {"dynamicbind_score": best_aff, "lddt": best_lddt}

def main() -> int:
    """"Create an output CSV score file combining DynamicBind and GNINA scores for a given ligand."""
    ap = argparse.ArgumentParser()
    ap.add_argument("--name", required=True, help="Ligand name (e.g., {LIGAND})")
    ap.add_argument("--search-dir", required=True, help="DynamicBind output dir searched by select_best (e.g., results/{LIGAND})")
    ap.add_argument("--gnina-dat", required=True, help="GNINA rescore .dat path")
    ap.add_argument("--db-meta", default=None, help="Optional JSON from select_best_dynamicbind_outputs.py")
    ap.add_argument("--out", required=True, help="Output CSV path (e.g., results/{LIGAND}_scores.csv)")
    args = ap.parse_args()

    # make sure paths contain the files
    search_dir = Path(args.search_dir).resolve() 
    gnina_dat = Path(args.gnina_dat).resolve()
    out = Path(args.out)

    if not search_dir.exists():
        raise FileNotFoundError(search_dir)
    if not gnina_dat.exists():
        raise FileNotFoundError(gnina_dat)

    meta_path = Path(args.db_meta).resolve() if args.db_meta else None

    # now extract the values
    db = parse_dynamicbind_from_meta_or_filename(meta_path, search_dir)
    gn = parse_gnina_dat(gnina_dat)

    if gn["gnina_affinity"] is None or gn["gnina_affinity"] > 0:
        raise ValueError(f"Invalid GNINA affinity from {gnina_dat}: {gn['gnina_affinity']}")
    
    # Format to match example output
    # name,dynamicbind_score,gnina_affinity,gnina_cnnaffinity,lddt
    # ligand_1,4.74,-6.53057,6.37874,0.7

    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as f:
        w = csv.DictWriter(
            f,
            fieldnames=["name", "dynamicbind_score", "gnina_affinity", "gnina_cnnaffinity", "lddt"],
        )
        w.writeheader()
        w.writerow(
            {
                "name": args.name,
                "dynamicbind_score": db["dynamicbind_score"],
                "gnina_affinity": gn["gnina_affinity"],
                "gnina_cnnaffinity": gn["gnina_cnnaffinity"],
                "lddt": db["lddt"],
            }
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())