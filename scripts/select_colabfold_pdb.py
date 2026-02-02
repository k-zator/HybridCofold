import argparse
import shutil
from pathlib import Path

def pick_pdb(colabfold_out: Path) -> Path:
    """Pick a PDB file from ColabFold output directory.
       Prefer "relaxed_rank_001" if present, otherwise take the first .pdb found."""

    preferred = sorted(colabfold_out.rglob("*relaxed_rank_001*.pdb")) # default colabfold naming
    if preferred:
        return preferred[0]

    any_pdb = sorted(colabfold_out.rglob("*.pdb"))
    if any_pdb:
        return any_pdb[0]

    raise FileNotFoundError(f"No .pdb files found under {colabfold_out}")

def main():
    """Select a PDB file from ColabFold output directory and copy it to the specified output path."""
    ap = argparse.ArgumentParser()
    ap.add_argument("--colabfold_out", required=True, help="ColabFold output directory (e.g., results/colabfold)")
    ap.add_argument("--out", required=True, help="Output protein pdb path")
    args = ap.parse_args()

    colabfold_out = Path(args.colabfold_out).resolve()
    if not colabfold_out.exists():
        raise FileNotFoundError(colabfold_out)

    src = pick_pdb(colabfold_out)
    dst = Path(args.out)
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(src, dst)
    return 0

if __name__ == "__main__":
    raise SystemExit(main())