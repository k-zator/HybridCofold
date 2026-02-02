import argparse
import json
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

# DynamicBind "complex" filename convention
COMPLEX_RE = re.compile(
    r"^rank\d+_complex_lddt(?P<lddt>\d+(?:\.\d+)?)_affinity(?P<aff>-?\d+(?:\.\d+)?)_relaxed\.pdb$"
)


def _count_atoms_in_pdb(p: Path) -> int:
    n = 0
    for ln in p.read_text().splitlines():
        if ln.startswith("ATOM") or ln.startswith("HETATM"):
            n += 1
    return n


def _split_complex_with_obabel(complex_pdb: Path, out_protein_pdb: Path, out_ligand_sdf: Path) -> None:
    obabel = shutil.which("obabel") or shutil.which("babel")
    if not obabel:
        raise FileNotFoundError(
            "OpenBabel not found (expected `obabel` or `babel` on PATH). "
            "Install it in the active environment (e.g. posebusters env)."
        )

    with tempfile.TemporaryDirectory(prefix="split_complex_") as td:
        tmpdir = Path(td)
        # OpenBabel will create part1.pdb, part2.pdb, ... in tmpdir
        prefix = tmpdir / "part.pdb"
        subprocess.run(
            [obabel, str(complex_pdb), "-m", "--separate", "-O", str(prefix)],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )

        parts = sorted(tmpdir.glob("part*.pdb"))
        if len(parts) < 2:
            raise RuntimeError(
                f"OpenBabel did not produce multiple molecules from {complex_pdb}. "
                f"Expected at least 2 parts, got {len(parts)}."
            )

        parts_by_atoms = sorted(((p, _count_atoms_in_pdb(p)) for p in parts), key=lambda x: x[1])
        ligand_pdb, ligand_atoms = parts_by_atoms[0]
        protein_pdb, protein_atoms = parts_by_atoms[-1]
        if ligand_atoms <= 0 or protein_atoms <= 0:
            raise RuntimeError(f"Bad split from {complex_pdb}: ligand_atoms={ligand_atoms}, protein_atoms={protein_atoms}")

        out_protein_pdb.parent.mkdir(parents=True, exist_ok=True)
        out_ligand_sdf.parent.mkdir(parents=True, exist_ok=True)

        shutil.copyfile(protein_pdb, out_protein_pdb)

        # Convert ligand PDB -> SDF
        subprocess.run(
            [obabel, str(ligand_pdb), "-O", str(out_ligand_sdf)],
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )

def main() -> int:
    """Select best complex pose from DynamicBind outputs and split into protein+ligand."""
    ap = argparse.ArgumentParser()
    ap.add_argument("--search-dir", required=True, help="Directory to search (e.g., results/{LIGAND})")
    ap.add_argument("--out-protein", required=True, help="Where to copy best receptor pdb")
    ap.add_argument("--out-ligand", required=True, help="Where to copy best ligand sdf")
    ap.add_argument("--out-metadata", default=None,
        help="Optional JSON output path with selected lddt/affinity and source paths",
    )
    args = ap.parse_args()

    search_dir = Path(args.search_dir)
    if not search_dir.is_dir():
        raise FileNotFoundError(f"DynamicBind results directory not found: {search_dir}")

    best_lddt = float("-inf")
    best_aff = float("-inf")
    best_complex = None

    for complex_path in search_dir.rglob("rank*_complex*_relaxed.pdb"):
        m = COMPLEX_RE.match(complex_path.name)
        if not m:
            continue
        lddt = float(m.group("lddt"))
        aff = float(m.group("aff"))
        if lddt > best_lddt or (lddt == best_lddt and aff > best_aff):
            best_lddt = lddt
            best_aff = aff
            best_complex = complex_path

    if best_complex is None:
        raise FileNotFoundError(
            f"No relaxed complex PDB found under {search_dir} (expected rank*_complex_*_relaxed.pdb)."
        )

    out_protein = Path(args.out_protein)
    out_ligand = Path(args.out_ligand)
    out_protein.parent.mkdir(parents=True, exist_ok=True)
    out_ligand.parent.mkdir(parents=True, exist_ok=True)

    _split_complex_with_obabel(best_complex, out_protein, out_ligand)

    # write JSON to rely the binding information
    if args.out_metadata:
        meta_path = Path(args.out_metadata)
        meta_path.parent.mkdir(parents=True, exist_ok=True)
        meta = {
            "lddt": best_lddt,
            "dynamicbind_score": best_aff,
        }
        meta_path.write_text(json.dumps(meta, indent=2) + "\n")

    print(f"Selected best pair: lddt={best_lddt} affinity={best_aff}")
    print(f"Complex: {best_complex}")
    print(f"Protein: {out_protein}")
    print(f"Ligand:  {out_ligand}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())