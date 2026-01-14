from ase.visualize import view
from ase.build import bulk
from catkit.gen.surface import SlabGenerator as Slab

a = bulk('Cu')
b = Slab(a, (2,1,1), 6, 20, 0, 'sym',standardize_bulk=True).get_slab().repeat([1,1,1])
view(b)
