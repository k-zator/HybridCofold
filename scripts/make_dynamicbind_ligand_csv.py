import argparse
from pathlib import Path

def main():
    """Create a ligand.csv file for DynamicBind - formatted to match input spec."""
    ap = argparse.ArgumentParser()
    ap.add_argument("--ligand-sdf", required=True, help="Path to ligand PDB (e.g., results/ligands/{LIGAND}.sdf)")
    ap.add_argument("--out", required=True, help="Output ligand.csv path")
    args = ap.parse_args()

    ligand_sdf = Path(args.ligand_sdf).resolve()
    if not ligand_sdf.exists():
        raise FileNotFoundError(f"ligand sdf not found: {ligand_sdf}")

    # File to look like:
    # ligand
    # /path/to/ligand.sdf

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("ligand\n" + str(ligand_sdf) + "\n")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())