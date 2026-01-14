from pymatgen.io.ase import AseAtomsAdaptor as aseatoms
from ase.visualize import view
from pymatgen.core import Structure, Lattice, Molecule
from pymatgen.core.surface import SlabGenerator
from pymatgen.io.vasp import Poscar
from pymatgen.analysis.adsorption import AdsorbateSiteFinder


cu_bulk = Structure.from_spacegroup("Fm-3m", Lattice.cubic(3.61), ["Cu"], [[0,0,0]]) ## meaning of coordinates?
co_molecule = Molecule(["C", "O"], [[0,0,0], [0,0,1.13]]) ## what is the meaning of coordinates?
gen = SlabGenerator(cu_bulk, 
              miller_index=(2, 1, 1), 
             min_slab_size=9.0, 
              min_vacuum_size=15.0, 
              center_slab=True, lll_reduce=True)

slabs = gen.get_slabs()
cu_211_slab = slabs[0]

orthogonal_slab = cu_211_slab.get_orthogonal_c_slab() ## this corrects tilted c axis, but it would alter inherent symmetry of structure.

print("before_conversion(alpha, beta, gamma):", cu_211_slab.lattice.angles)
print("after_conversion(alpha, beta, gamma):", orthogonal_slab.lattice.angles)
# 4. 3x4 Supercell
# surface unit cell
orthogonal_slab.make_supercell([[3, 0, 0], [0, 1, 0], [0, 0, 1]])


asf = AdsorbateSiteFinder(orthogonal_slab)
sites_dict = asf.find_adsorption_sites(distance=1.0)
adsorbed_structures = asf.generate_adsorption_structures(
    co_molecule, 
    repeat=[1, 1, 1], 
    min_lw=10.0,         # Minimum vacuum size
    find_args={"distance": 2.0}  # Height of C atom above surface
)
print(f"Total distinct sites found: {len(adsorbed_structures)}")

for label, coords_list in sites_dict.items():
    if label == 'all': continue  # Skip the summary list
    
    for i, coord in enumerate(coords_list):
        struct = orthogonal_slab.copy()
        mol_copy = co_molecule.copy()
        mol_copy.translate_sites(indices=range(len(mol_copy)),vector=coord)
        # 3. Save immediately with the correct name
        for atom in mol_copy:
            struct.append(atom.specie, atom.coords, coords_are_cartesian=True)

        filename = f"POSCAR_{label}_{i}.vasp"
        Poscar(struct).write_file(filename)
        
        print(f"Saved: {filename}")





poscar = Poscar(orthogonal_slab)
poscar.write_file("POSCAR_Cu211_Base")
view(aseatoms.get_atoms(orthogonal_slab))

