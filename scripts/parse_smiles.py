import argparse
from pathlib import Path

def main():
    """Parse a .smi file and write the first SMILES to an output text file.
       Assure correct formatting for downstream processing."""
    ap = argparse.ArgumentParser()
    ap.add_argument("--smi-file", required=True, help="Input .smi file (first token on first non-empty line is used)")
    ap.add_argument("--out", required=True, help="Output text file containing exactly one SMILES")
    args = ap.parse_args()

    smi_path = Path(args.smi_file)
    if not smi_path.exists():
        raise FileNotFoundError(smi_path)

    smiles = None
    for line in smi_path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"): # ignore comments and empty lines
            continue
        smiles = line.split()[0]  # first whitespace-separated token
        break

    if not smiles:
        raise ValueError(f"No SMILES found in {smi_path}")

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(smiles + "\n")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())