import argparse
import os
from pathlib import Path

from rdkit_utils import generate_3D_conformers_from_smiles

def main():
    """Generate 3D conformers from SMILES (ligand.smi) input file."""
    ap = argparse.ArgumentParser()
    ap.add_argument("--smiles-file", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--num-conformers", type=int, default=100)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--outname", default="ligand_1.sdf")
    args = ap.parse_args()

    smiles = Path(args.smiles_file).read_text().strip()
    if not smiles:
        raise ValueError("Empty SMILES")
    
    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    outpath = outdir / args.outname

    # Run the conformer generation with RDKit
    generate_3D_conformers_from_smiles(
        smiles,
        num_conformers=args.num_conformers,
        random_seed=args.seed,
        out_sdf=str(outpath),
    )

    if not outpath.exists():
        raise FileNotFoundError(f"Expected {outpath} to be created")

    return 0

if __name__ == "__main__":
    raise SystemExit(main())