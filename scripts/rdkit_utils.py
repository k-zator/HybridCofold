#!/usr/bin/env python3
from rdkit import Chem # type: ignore
from rdkit.Chem import AllChem # type: ignore
from pathlib import Path

def generate_3D_conformers_from_smiles(smiles_list, num_conformers=100, random_seed=42, out_sdf="ligand.sdf"):
    """Generate 3D conformers from a SMILES strings.
       Area for further improvement: ensemble docking with varied generated conformers."""

    if isinstance(smiles_list, str):
        smiles_list = [smiles_list]
    if len(smiles_list) != 1:
        raise ValueError("This pipeline currently expects exactly one SMILES (ligand) per run.")
    
    smi = smiles_list[0].strip()
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smi}")
    
    mol = Chem.AddHs(mol, addCoords=True)

    params = AllChem.ETKDG()
    params.randomSeed = random_seed # reproducibility
    params.enforceChirality = True # respect stereochemistry
    params.useExpTorsionAnglePrefs = True # use experimental torsion angle for better conformers
    params.useBasicKnowledge = True # use basic chemical knowledge for better conformers

    cids = AllChem.EmbedMultipleConfs(mol, numConfs=num_conformers, params=params) # Generate conformers
    results = AllChem.MMFFOptimizeMoleculeConfs(mol, numThreads=0) # Optimize conformers
    if not results:
        raise RuntimeError(f"MMFF optimization failed for SMILES: {smi}")
    
    min_energy_idx = min(range(len(results)), key=lambda i: results[i][1]) # Best conformer
    best_conf_id = cids[min_energy_idx]

    out_sdf = Path(out_sdf)
    out_sdf.parent.mkdir(parents=True, exist_ok=True)
    w = Chem.SDWriter(str(out_sdf))
    w.write(mol, confId=int(best_conf_id))
    w.close()

def protonate_molecule(mol_path, out_path, pH=7.4):
    """Protonate a molecule from a file at a given pH and save as PDB or SDF. 
       To be used for both proteins and ligands."""

    mol_path = str(mol_path)
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    # check extension to read correctly
    if mol_path.endswith('.pdb'):
        mol = Chem.MolFromPDBFile(mol_path, removeHs=False)
        if mol is None:
            raise ValueError(f"Failed to read PDB: {mol_path}")
        mol = Chem.AddHs(mol, addCoords=True)
        Chem.MolToPDBFile(mol, str(out_path))

    elif mol_path.endswith('.sdf'):
        suppl = Chem.SDMolSupplier(mol_path, removeHs=False)
        mol = next((m for m in suppl if m is not None), None)
        if mol is None:
            raise ValueError(f"Failed to read SDF: {mol_path}")
        mol = Chem.AddHs(mol, addCoords=True)
        w = Chem.SDWriter(str(out_path))
        w.write(mol)
        w.close()   
    else:
        raise ValueError("Unsupported file format. Please provide a PDB or SDF file.")
    


