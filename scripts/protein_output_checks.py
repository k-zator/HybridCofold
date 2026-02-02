import argparse
import math
from pathlib import Path

# headless-safe
import matplotlib # type: ignore
matplotlib.use("Agg")
import matplotlib.pyplot as plt # type: ignore

from Bio.PDB import PDBParser # type: ignore
from Bio.PDB.Polypeptide import PPBuilder # type: ignore


def plot_ramachandran(pdb_path: str, out_png: str, title: str | None = None) -> None:
    """Generate a Ramachandran plot from a protein PDB file to visually check protein prediction quality."""

    pdb_path = str(pdb_path)
    out_png = Path(out_png)
    out_png.parent.mkdir(parents=True, exist_ok=True)

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("prot", pdb_path)

    ppb = PPBuilder()
    phis: list[float] = []
    psis: list[float] = []

    for model in structure:
        for chain in model:
            for poly in ppb.build_peptides(chain):
                for phi, psi in poly.get_phi_psi_list():
                    if phi is None or psi is None:
                        continue
                    phis.append(math.degrees(phi))
                    psis.append(math.degrees(psi))

    if not phis:
        raise ValueError(f"No phi/psi angles found (is this a protein PDB?): {pdb_path}")

    fig = plt.figure(figsize=(6, 6), dpi=160)
    ax = fig.add_subplot(111)
    ax.scatter(phis, psis, s=8, alpha=0.6)

    ax.set_xlim(-180, 180)
    ax.set_ylim(-180, 180)
    ax.set_xlabel("Phi (deg)")
    ax.set_ylabel("Psi (deg)")
    ax.set_title(title or "Ramachandran plot")
    ax.grid(True, linewidth=0.3, alpha=0.4)

    fig.tight_layout()
    fig.savefig(out_png)
    plt.close(fig)


def main() -> int:
    """Run sanity check for generated protein structure by creating a Ramachandran plot and with quality."""
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb", required=True, help="Input protein PDB")
    ap.add_argument("--out-png", required=True, help="Output PNG path")
    ap.add_argument("--title", default=None)
    args = ap.parse_args()

    plot_ramachandran(args.pdb, args.out_png, title=args.title)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())