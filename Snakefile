#Making a pipeline for the co-folding task using snakemake
import os
from pathlib import Path

configfile: "snakemake_config.yaml"

RANDOMSEED = config["inputs"]["random_seed"] # Set for reproducibility
OUTDIR     = config["paths"]["outdir"]
FASTA      = config["inputs"]["protein_fasta"]
SMI        = config["inputs"]["ligands_smi"]
LIGAND     = Path(SMI).stem # Basename for folder/file naming

rule all:
    input:
        f"{OUTDIR}/colabfold/done.stamp",
        f"{OUTDIR}/ligands/{LIGAND}.sdf",
        f"{OUTDIR}/{LIGAND}_scores.csv",

#############################################
# 1) Setup localcolabfold env via pixi
#############################################

rule setup_localcolabfold:
    # The pipeline requires three environments, this is No.1 with explicit custom pixi installation
    output:
        stamp = f"{OUTDIR}/setup/localcolabfold.stamp"
    params:
        pixi = lambda wc: config["paths"].get("pixi", "pixi"),
        repo = lambda wc: config["paths"].get("localcolabfold_repo", "https://github.com/yoshitakamo/localcolabfold.git"),
        lc_dir = lambda wc: config["paths"]["localcolabfold_dir"],
    shell:
        r"""
        set -euo pipefail
        mkdir -p {OUTDIR}/setup

        # Ensure pixi exists for PATH persistence
        if ! command -v {params.pixi} >/dev/null 2>&1; then
          echo "ERROR: pixi not found. Install pixi or set paths.pixi to the pixi binary in snakemake_config.yaml" >&2
          exit 2
        fi

        # Clone if needed
        if [ ! -d "{params.lc_dir}/.git" ]; then
          git clone "{params.repo}" "{params.lc_dir}"
        fi

        echo "Setting up localcolabfold env"

        cd "{params.lc_dir}"
        {params.pixi} install
        {params.pixi} run setup

        touch "{workflow.basedir}/{output.stamp}"
        """

############################################
# 2) CUDA checks (per-env)
############################################

rule cuda_check_colabfold:
    # Check CUDA connection as calculations run primarily on GPUs 
    input:
        setup = f"{OUTDIR}/setup/localcolabfold.stamp"
    output:
        stamp = f"{OUTDIR}/setup/cuda_colabfold.stamp"
    params:
        lc_dir = lambda wc: config["paths"]["localcolabfold_dir"],
        pixi = lambda wc: config["paths"].get("pixi", "pixi")
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{OUTDIR}/setup"
        cd "{params.lc_dir}"

        "{params.pixi}" run python "{workflow.basedir}/scripts/check_cuda.py" \
          --framework auto --out "{workflow.basedir}/{output.stamp}"
        """

rule cuda_check_dynamicbind:
    # Check CUDA connection as calculations run primarily on GPUs 
    output:
        stamp = f"{OUTDIR}/setup/cuda_dynamicbind.stamp"
    conda:
        "envs/dynamicbind.yml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{OUTDIR}/setup"
        python "{workflow.basedir}/scripts/check_cuda.py" \
          --framework torch --out "{workflow.basedir}/{output.stamp}"
        """

#############################################
# 3) ColabFold run + Protein structure checks
#############################################

rule run_colabfold:
    # Full protein structure prediction run from FASTA to PDB with localcolabfold via pixi 
    input:
        setup = f"{OUTDIR}/setup/localcolabfold.stamp",
        fasta = FASTA
    output:
        stamp = f"{OUTDIR}/colabfold/done.stamp"
    params:
        num_recycle = lambda wc: config["colabfold"]["num_recycle"],
        num_model = lambda wc: config["colabfold"]["num_model"],
        pixi = lambda wc: config["paths"].get("pixi", "pixi"),
        lc_dir = lambda wc: config["paths"]["localcolabfold_dir"],
        script = lambda wc: config["paths"]["colabfold_script"],
        outdir = f"{OUTDIR}/colabfold"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.outdir}
        cd {params.lc_dir}

        echo "Running protein structure prediction"
        # run colab fold
        {params.pixi} run bash "{workflow.basedir}/scripts/{params.script}" \
          "{workflow.basedir}/{input.fasta}" \
          "{workflow.basedir}/{params.outdir}" \
          "{params.num_recycle}" \
          "{params.num_model}" \
          "{RANDOMSEED}"
        touch {workflow.basedir}/{output.stamp}
        """

rule protein_output_checks:
    # Sanity checks for protein using a Ramachandran, more checks after docking
    input:
        done = f"{OUTDIR}/colabfold/done.stamp"
    output:
        png = f"{OUTDIR}/colabfold/ramachandran.png",
        stamp = f"{OUTDIR}/colabfold/checks.stamp"
    conda:
        "envs/posebusters.yml"
    params:
        selected_pdb = f"{OUTDIR}/colabfold/selected_protein.pdb"
    shell:
        r"""
        set -euo pipefail

        echo "Checking generated protein structure"

        # find the structure in output folder, based on custom name
        python "{workflow.basedir}/scripts/select_colabfold_pdb.py" \
          --colabfold_out "{OUTDIR}/colabfold" \
          --out "{params.selected_pdb}"

        # and run code for Ramachandran plot
        python "{workflow.basedir}/scripts/protein_output_checks.py" \
                    --pdb "{params.selected_pdb}" \
                    --out-png "{output.png}" \
                    --title "Ramachandran plot for generated protein"

        touch "{output.stamp}"
        """

#############################################
# 4) Generate 3D ligand conformers + checks
#############################################

rule parse_smiles:
    # Read in ligand.smi and check formatting
    input:
        smi = SMI
    output:
        smiles = f"{OUTDIR}/ligands/smiles.txt"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {OUTDIR}/ligands
        python {workflow.basedir}/scripts/parse_smiles.py \
          --smi-file "{input.smi}" --out "{output.smiles}"
        """

rule ligand_3d:
    # Generate 3D conformations of ligand and select the most stable one
    input:
        smiles = f"{OUTDIR}/ligands/smiles.txt"
    output:
        pdb = f"{OUTDIR}/ligands/{LIGAND}.sdf"
    conda:
        "envs/dynamicbind.yml"
    params:
        outdir = f"{OUTDIR}/ligands"
    shell:
        r"""
        set -euo pipefail

        echo "Generating 3D conformations for ligand"

        PYTHONPATH="{workflow.basedir}" \
          python "{workflow.basedir}/scripts/smiles_to_3d.py" \
          --smiles-file "{input.smiles}" \
          --outdir "{params.outdir}" \
          --seed "{RANDOMSEED}" \
          --outname "{LIGAND}.sdf"
        """

#############################################
# 5) DynamicBind Install, Score & Relax
#############################################

rule setup_dynamicbind:
    # Move model parameters directory to DynamicBind installation site
    output:
        stamp = f"{OUTDIR}/setup/dynamicbind_workdir.stamp"
    params:
        repo = lambda wc: config["paths"]["dynamicbind_repo"],
        ddir = lambda wc: config["paths"]["dynamicbind_dir"],
        workdir_src = lambda wc: config["dynamicbind"]["workdir_src"],
    shell:
        r"""
        set -euo pipefail

        echo "Setting up dynamicbind env"

        if [ ! -d "{params.ddir}/.git" ]; then
          git clone "{params.repo}" "{params.ddir}"
        fi

        SRC="{workflow.basedir}/{params.workdir_src}"
        DST="{params.ddir}/workdir"

        if [ -e "$DST" ] && [ ! -L "$DST" ]; then
          rm -rf "$DST"
        fi
    
        cp -r "$SRC" "$DST"

        # Sanity check: confirm the workdir has been moved correctly 
        # as DynamicBind hardcodes paths to the model weights
        if ! find "$DST" -maxdepth 4 -name "*.pt" | head -n 1 | grep -q ".pt"; then
          echo "ERROR: No .pt model parameter files found under $DST" >&2
          exit 4
        fi

        touch "{workflow.basedir}/{output.stamp}"
        """

rule run_dynamicbind:
    # A single-command but two-stage process of dynamic docking with DynamicBind (+ relax)
    input:
        protein_done = f"{OUTDIR}/colabfold/done.stamp",
        ligand_done = f"{OUTDIR}/ligands/{LIGAND}.sdf",
        weights = f"{OUTDIR}/setup/dynamicbind_workdir.stamp"
    output:
        stamp = f"{OUTDIR}/binding/done.stamp"
    conda:
        "envs/dynamicbind.yml"
    params:
        outdir = lambda wc: f"{OUTDIR}/binding",
        ddir = lambda wc: config["paths"]["dynamicbind_dir"],
        relax_env = lambda wc: config["env_names"]["relax"], # third env implicitly read
        savings = lambda wc: config["dynamicbind"]["savings_per_complex"], # run variable
        steps = lambda wc: config["dynamicbind"]["inference_steps"], #run variable
        header = LIGAND
    shell:
        r"""
        set -euo pipefail
        mkdir -p {params.outdir}

        echo "Running co-folded docking with DynamicBind"

        # Create custom ligand.csv file with {LIGAND}.sdf path - required formatting
        python "{workflow.basedir}/scripts/make_dynamicbind_ligand_csv.py" \
          --ligand-sdf "{input.ligand_done}" \
          --out "{params.outdir}/ligand.csv"

        # Detect and copy the correct protein PDB from colabfold results (the relaxed version)
        python "{workflow.basedir}/scripts/select_colabfold_pdb.py" \
          --colabfold_out "{OUTDIR}/colabfold" \
          --out "{params.outdir}/protein.pdb"

        # Set second env for MM relaxation of structure
        RELAX_PY="$(conda run -n "{params.relax_env}" python -c 'import sys; print(sys.executable)')"
        # mkl export necessary for compatibility with libgomp.so.1
        export MKL_THREADING_LAYER=GNU

        python "{params.ddir}/run_single_protein_inference.py" -p \
          "{params.outdir}/protein.pdb" \
          "{params.outdir}/ligand.csv" \
          --savings_per_complex "{params.savings}" \
          --inference_steps "{params.steps}" \
          --header "{params.header}" \
          --python "$(which python)" \
          --relax_python "$RELAX_PY"

        touch "{output.stamp}"
        """

#############################################
# 6) GNINA 1.3.2 Rescore
#############################################

rule setup_gnina:
    # GNINA download and activation
    output:
        stamp = f"{OUTDIR}/setup/gnina.stamp"
    params:
        bin = lambda wc: config["paths"]["gnina_bin"],
        url = lambda wc: config["paths"]["gnina_repo"]
    shell:
        r"""
        set -euo pipefail
        mkdir -p tools/gnina "{OUTDIR}/setup"

        if [ ! -x "{params.bin}" ]; then
          wget -O "{params.bin}" "{params.url}"
          chmod +x "{params.bin}"
        fi

        touch "{workflow.basedir}/{output.stamp}"
        """

rule obtain_docked_protein_ligand:
    # Select best DynamicBind complex pose and split into protein + ligand
    # DynamicBind outputs are under results/{ligand}/index0_idx_0 and will be copied to rescoring folder
    input:
        bind_done = f"{OUTDIR}/binding/done.stamp",
    output:
        proteinH = f"{OUTDIR}/rescoring/{LIGAND}/proteinH.pdb",
        ligandH  = f"{OUTDIR}/rescoring/{LIGAND}/ligandH.sdf",
        db_meta  = f"{OUTDIR}/rescoring/{LIGAND}/dynamicbind_selected.json"
    conda:
        "envs/posebusters.yml"
    params:
        search_dir = f"{OUTDIR}/{LIGAND}",
                db_meta = f"{OUTDIR}/rescoring/{LIGAND}/dynamicbind_selected.json"
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{OUTDIR}/rescoring/{LIGAND}"

        python "{workflow.basedir}/scripts/select_best_dynamicbind_outputs.py" \
          --search-dir "{params.search_dir}" \
                    --out-protein "{output.proteinH}" \
                    --out-ligand "{output.ligandH}" \
          --out-metadata "{params.db_meta}"
        """

rule run_rescore_gnina:
    # Run rescore only on top pose from DynamicBind
    input:
        setup = f"{OUTDIR}/setup/gnina.stamp",
        proteinH = f"{OUTDIR}/rescoring/{LIGAND}/proteinH.pdb",
        ligandH  = f"{OUTDIR}/rescoring/{LIGAND}/ligandH.sdf"
    output:
        dat = f"{OUTDIR}/rescoring/{LIGAND}/{LIGAND}_rescore.dat",
        stamp = f"{OUTDIR}/rescoring/{LIGAND}/rescore.stamp",
    params:
        gnina = lambda wc: config["paths"]["gnina_bin"],
        cudnn_lib = lambda wc: os.path.expanduser(config["paths"]["cudnn_lib"])
    shell:
        r"""
        set -euo pipefail
        mkdir -p "{OUTDIR}/rescoring/{LIGAND}"

        export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:{params.cudnn_lib}"
    
        echo "Rescoring top pose with GNINA"

        "{params.gnina}" --score_only --cnn_scoring rescore \
          -r "{input.proteinH}" \
          -l "{input.ligandH}" > "{output.dat}"

        touch "{workflow.basedir}/{output.stamp}"
        """

rule protein_pose_checks:
    # Geometrical and chemical sanity checks re structure quality 
    input:
        proteinH = f"{OUTDIR}/rescoring/{LIGAND}/proteinH.pdb",
        ligandH  = f"{OUTDIR}/rescoring/{LIGAND}/ligandH.sdf",
        ligand_before  = f"{OUTDIR}/ligands/{LIGAND}.sdf",
        protein_before = f"{OUTDIR}/binding/protein.pdb",
        done = f"{OUTDIR}/binding/done.stamp"
    output:
        report_ligand = f"{OUTDIR}/rescoring/{LIGAND}/posebusters_report_ligand.txt",
        report_ligand2ligand = f"{OUTDIR}/rescoring/{LIGAND}/posebusters_report_ligand2ligand.txt",
        report_proteinligand = f"{OUTDIR}/rescoring/{LIGAND}/posebusters_report_proteinligand.txt",
        stamp  = f"{OUTDIR}/rescoring/{LIGAND}/posebusters.stamp"
    conda:
        "envs/posebusters.yml"
    shell:
        r"""
        set -euo pipefail

        echo "PoseBusters sanity checks"
        
        # ligand conformation check
        bust --outfmt long "{input.ligandH}" > "{output.report_ligand}" 2>&1
        # ligand conformational change check
        bust --outfmt long "{input.ligandH}" -l "{input.ligand_before}" > "{output.report_ligand2ligand}" 2>&1 
        # protein-ligand interaction check
        bust --outfmt long "{input.ligandH}" -p "{input.proteinH}" > "{output.report_proteinligand}" 2>&1

        # and atom count check before and after docking/rescoring
        python "{workflow.basedir}/scripts/check_atom_counts.py" \
          --before "{input.ligand_before}" \
          --after "{input.ligandH}" \
          --allow-delta 0

        python "{workflow.basedir}/scripts/check_atom_counts.py" \
          --before "{input.protein_before}" \
          --after "{input.proteinH}" \
          --allow-delta 0

        touch "{output.stamp}"
        """

rule summarise_scores:
    # Write out a csv file with collated binding scores and lddt
    input:
        gnina = f"{OUTDIR}/rescoring/{LIGAND}/{LIGAND}_rescore.dat",
        db_meta = f"{OUTDIR}/rescoring/{LIGAND}/dynamicbind_selected.json",
        bind_done = f"{OUTDIR}/binding/done.stamp"
    output:
        csv = f"{OUTDIR}/{LIGAND}_scores.csv"
    conda:
        "envs/dynamicbind.yml"
    params:
        search_dir = f"{OUTDIR}/{LIGAND}"
    shell:
        r"""
        set -euo pipefail

        echo "Summarising calculated scores {LIGAND}" # to explicitly print in CLI

        python "{workflow.basedir}/scripts/make_score_table.py" \
          --name "{LIGAND}" \
          --search-dir "{params.search_dir}" \
          --gnina-dat "{input.gnina}" \
          --db-meta "{input.db_meta}" \
          --out "{output.csv}"
        """