import argparse
from rdkit import Chem   # type: ignore

def atom_count_sdf(path: str) -> int:
    suppl = Chem.SDMolSupplier(path, removeHs=False)
    mol = next((m for m in suppl if m is not None), None)
    if mol is None:
        raise ValueError(f"Could not read SDF: {path}")
    return mol.GetNumAtoms()

def atom_count_pdb(path: str) -> int:
    mol = Chem.MolFromPDBFile(path, removeHs=False)
    if mol is None:
        raise ValueError(f"Could not read PDB: {path}")
    return mol.GetNumAtoms()

def main() -> int:
    """Check for atom count consistency between two SDF or PDB files."""
    ap = argparse.ArgumentParser()
    ap.add_argument("--before", required=True)
    ap.add_argument("--after", required=True)
    ap.add_argument("--allow-delta", type=int, default=0, help="Allowed difference in atom count")
    args = ap.parse_args()

    # for proteins (PDB extension or ligand SDF)
    if args.before.lower().endswith(".sdf"):
        b = atom_count_sdf(args.before)
        a = atom_count_sdf(args.after)
    elif args.before.lower().endswith(".pdb"):
        b = atom_count_pdb(args.before)
        a = atom_count_pdb(args.after)
    else:
        raise SystemExit("Unsupported file extension. Use .sdf or .pdb")

    if abs(a - b) > args.allow_delta:
        raise SystemExit(
            f"Atom count mismatch: before={b} after={a} "
            f"(allowed delta={args.allow_delta}).\n"
            f"before={args.before}\nafter={args.after}"
        )
    return 0

if __name__ == "__main__":
    raise SystemExit(main())