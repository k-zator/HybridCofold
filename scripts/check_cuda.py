import argparse
import os
import shutil
import subprocess
from pathlib import Path

def run(cmd: list[str]) -> tuple[int, str]:
    """Run a command to test CUDA availability."""
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, check=False)
        out = (p.stdout or "") + (p.stderr or "")
        return p.returncode, out.strip()
    except Exception as e:
        return 1, str(e)

def check_nvidia_smi() -> dict:
    """Check if nvidia-smi is available and return its status."""
    if not shutil.which("nvidia-smi"):
        return {"available": False, "detail": "nvidia-smi not found"}
    rc, out = run(["nvidia-smi", "-L"])
    return {"available": rc == 0, "detail": out}

def check_torch() -> dict:
    """Check CUDA availability via PyTorch for GNINA rescoring."""
    try:
        import torch  # type: ignore

        ok = bool(torch.cuda.is_available())
        n = int(torch.cuda.device_count()) if ok else 0
        names = []
        if ok:
            for i in range(n):
                try:
                    names.append(torch.cuda.get_device_name(i))
                except Exception:
                    names.append(f"cuda:{i}")
        return {"available": ok, "devices": n, "names": names, "version": getattr(torch, "__version__", "unknown")}
    except Exception as e:
        return {"available": False, "detail": f"torch import failed: {e}"}

def main() -> int:
    """Check CUDA availability and write a report to a file. 
       The pipeline is critically dependent on CUDA availability for GPU-accelerated docking and rescoring."""
    ap = argparse.ArgumentParser()
    ap.add_argument("--framework", choices=["auto", "torch"], default="auto")
    ap.add_argument("--out", required=True, help="Output stamp file")
    args = ap.parse_args()

    report = {
        "env": os.environ.get("CONDA_PREFIX") or os.environ.get("VIRTUAL_ENV") or "",
        "nvidia_smi": check_nvidia_smi(),
    }

    if args.framework in ("auto", "torch"):
        report["torch"] = check_torch()

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(str(report) + "\n")

    # Non-zero exit if requested framework CUDA is unavailable
    if args.framework == "torch" and not report.get("torch", {}).get("available", False):
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())