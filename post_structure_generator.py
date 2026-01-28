from pymatgen.io.ase import AseAtomsAdaptor as aseatoms
from ase.visualize import view
from pymatgen.core import Structure, Lattice, Molecule
from pymatgen.core.surface import SlabGenerator
from pymatgen.io.vasp import Poscar
from pymatgen.analysis.adsorption import AdsorbateSiteFinder
import numpy as np

base_structure = Structure.from_file("POSCAR_Cu211_Base")
co_molecule = Molecule(["C", "O"], [[0,0,0], [0,0,1.13]]) ## what is the meaning of coordinates?

asf = AdsorbateSiteFinder(base_structure)
sites_dict = asf.find_adsorption_sites(distance=1.0)



for label, coords_list in sites_dict.items():
    if label == 'all': continue  # Skip the summary list
    
    for i, top_coord in enumerate(coords_list):
        struct = base_structure.copy()
        mol_top = co_molecule.copy()
        mol_top.translate_sites(indices=range(len(mol_top)), vector=top_coord)

        for atom in mol_top:
            struct.append(atom.specie, atom.coords, coords_are_cartesian=True, properties={"selective_dynamics": [True, True, True]})
            frac_top = struct.lattice.get_fractional_coords(top_coord)

            frac_bottom = 1.0 - frac_top

            bottom_coord = struct.lattice.get_cartesian_coords(frac_bottom)
            mol_bottom = co_molecule.copy()
            mol_bottom.rotate_sites(theta=180, axis=[1,0,0])
            mol_bottom.translate_sites(indices=range(len(mol_bottom)), vector = bottom_coord)

            for atom in mol_bottom:
                struct.append(atom.specie, atom.coords, coords_are_cartesian=True, properties={"selective_dynamics": [True, True, True]})

            
        filename = f"POSCAR_{label}_{i}.vasp"
        Poscar(struct).write_file(filename)
        
        print(f"Saved: {filename}")


