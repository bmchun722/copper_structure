from pymatgen.io.ase import AseAtomsAdaptor as aseatoms
from ase.visualize import view
from pymatgen.core import Structure, Lattice, Molecule
from pymatgen.core.surface import SlabGenerator
from pymatgen.io.vasp import Poscar
from pymatgen.analysis.adsorption import AdsorbateSiteFinder
import numpy as np

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
    find_args={"distance": 1.0}  # Height of C atom above surface
)
sites_dict = asf.find_adsorption_sites(distance=0.0)

target_neighbors = {
    'ontop': 1,
    'bridge': 2,
    'hollow': 3
}
for label, coords_list in sites_dict.items():
    if label == 'all': continue
    
    # Determine how many atoms support this type (e.g., bridge uses 2)
    n_target = target_neighbors.get(label, 1) 
    
    for i, crude_coord in enumerate(coords_list):
        
        # --- [1] Geometric Refinement of Coordinates ---
        
        # Find all Cu atoms within a sphere (radius 2.5A) around the coordinate.
        # (Considers Periodic Boundary Conditions)
        neighbors = orthogonal_slab.get_sites_in_sphere(crude_coord, r=2.5)
        
        # Sort by distance and pick only the closest n_target atoms.
        # e.g., Top 2 for Bridge, Top 3 for Hollow.
        sorted_neighbors = sorted(neighbors, key=lambda x: x[1]) # x[1] is distance
        closest_sites = [n[0] for n in sorted_neighbors[:n_target]]
        
        if not closest_sites: continue # Skip if no neighbors are found
        
        # Extract coordinates of the found Cu atoms
        cu_coords = [site.coords for site in closest_sites]
        cu_coords = np.array(cu_coords)
        
        # A. X, Y coordinates: Mean of Cu atoms (Perfect symmetry point)
        center_x = np.mean(cu_coords[:, 0])
        center_y = np.mean(cu_coords[:, 1])
        
        # B. Z coordinate: Based on the 'highest' atom among neighbors
        # This prevents the molecule from being buried in the step edge.
        max_z = np.max(cu_coords[:, 2])
        safe_z = max_z + 1.85  # Cu-C bond length approx 1.85 A
        
        # Finalize C position
        exact_c_pos = np.array([center_x, center_y, safe_z])
        
        
        # --- [2] Generate Structure ---
        struct = orthogonal_slab.copy()
        
        # Add C (at refined position)
        struct.append("C", exact_c_pos, coords_are_cartesian=True)
        
        # Add O (Vertical, directly above C)
        # Let VASP find the optimal tilt during relaxation.
        exact_o_pos = exact_c_pos + np.array([0, 0, 1.13])
        struct.append("O", exact_o_pos, coords_are_cartesian=True)
        
        
        # --- [3] Save ---
        filename = f"POSCAR_{label}_{i}_exact.vasp"
        Poscar(struct).write_file(filename)
        print(f"Saved: {filename} (Refined using {n_target} Cu atoms)")




poscar = Poscar(orthogonal_slab)
poscar.write_file("POSCAR_Cu211_Base")
view(aseatoms.get_atoms(orthogonal_slab))

