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

adsorbed_structures = asf.generate_adsorption_structures(
    co_molecule, 
    repeat=[1, 1, 1], 
    min_lw=10.0,         # Minimum vacuum size
    find_args={"distance": 2.0}  # Height of C atom above surface
)
print(adsorbed_structures)
print(f"Total distinct sites found: {len(adsorbed_structures)}")

for i, structure in enumerate(adsorbed_structures):
    # Save each structure to a separate POSCAR file
    filename = f"POSCAR_site_{i}.vasp"
    Poscar(structure).write_file(filename)
    print(f"Saved: {filename}")

# --- Extra: Visualizing Site Coordinates ---
# If you just want to see where the points are without generating structures:
sites = asf.find_adsorption_sites()
print("\nCoordinates of found sites:")
for site_type, coords in sites.items():
    print(f"{site_type}: {len(coords)} sites")

poscar = Poscar(orthogonal_slab)
poscar.write_file("POSCAR_Cu211_Base")
view(aseatoms.get_atoms(orthogonal_slab))

