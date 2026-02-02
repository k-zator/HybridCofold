import argparse
from rdkit_utils import protonate_molecule

def main():
    """Wrapper for protonating molecules using RDKit."""
    ap = argparse.ArgumentParser()
    ap.add_argument("--infile", required=True, help="Input molecule file (PDB or SDF)")
    ap.add_argument("--out", required=True, help="Output protonated molecule file (same format as input)")
    args = ap.parse_args()

    # protonation with stanard pH 7.4 to assure correct Gly/His/Ser/Thr/Tyr states
    protonate_molecule(args.infile, args.out, pH=7.4)
    return 0

if __name__ == "__main__":
    raise SystemExit(main())