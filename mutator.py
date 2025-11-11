import warnings
import argparse
import MDAnalysis as mda
import numpy as np

warnings.filterwarnings("ignore",message="Unit cell dimensions not found")
warnings.filterwarnings("ignore",message="Found no information for attr: 'formalcharges' Using default value of '0'")

parser = argparse.ArgumentParser()
parser.add_argument('-i',dest='in_file')
parser.add_argument('-o',dest='out_file')
args = parser.parse_args()

def main():
    u = mda.Universe(args.in_file)
    chain_id = set(u.atoms.chainIDs)
    chains = sorted(chain_id)

    selection = " or ".join([f"chainID {cid}" for cid in chains])
    target_atoms = u.select_atoms(f"(resname LYS or resname LYM or resname ARG) and ({selection})")
    atoms_to_remove = set()

    for res in target_atoms.residues:
        orig = res.resname

        if orig == 'ARG':
            for atom in res.atoms:
                if atom.name in {'CE','CZ','NH1','NH2','HZ1','HZ2','HZ3','HE1','HE2','1HH1','2HH1','1HH2','2HH2','HE'}:
                    atoms_to_remove.add(atom.index)
                if atom.name == 'CD':
                    atom.name = 'ND'
                    atom.element = 'N'
                if atom.name == 'NE':
                    atom.name = 'HD3'
                    atom.element = 'H'
            res.resname = 'DAB'

        elif orig in {'LYS','LYM'}:
            for atom in res.atoms:
                if atom.name in {'CE','NZ','NE','CZ','NH1','NH2','HZ1','HZ2','HZ3','HE1','HE2','1HH1','2HH1','1HH2','2HH2','HE','HG3','HD1','HD2'}:
                    atoms_to_remove.add(atom.index)
                if atom.name == 'CG':
                    atom.name = 'NG'
                    atom.element = 'N'
                if atom.name == 'CD':
                    atom.name = 'HG3'
                    atom.element = 'H'
            res.resname = 'DPR'
    mask = ~np.isin(u.atoms.indices, list(atoms_to_remove))
    kept_atoms = u.atoms[mask]
    kept_atoms.write(args.out_file)
    print(f"\nMutated structure save to: {args.out_file}")

main()
