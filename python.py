from rdkit import Chem
from rdkit.Chem.rdMolTransforms import ComputeCentroid
from vina import Vina
v = Vina(sf_name='vina')
lig = Chem.SDMolSupplier('ligand1.sdf', removeHs=False)[0]
 
#RDKit can't read pdbqt format so I use SDF for centroid calculation
centroid = ComputeCentroid(lig.GetConformer())
 
v.set_receptor('protein.pdbqt')
v.set_ligand_from_file('ligand.pdbqt')
 
print(centroid.x, centroid.y, centroid.z)

 
v.compute_vina_maps(center=[centroid.x, centroid.y, centroid.z], box_size=[20, 20, 20])
 
energy_minimized = v.optimize()
print('Score after minimization : %.3f (kcal/mol)' % energy_minimized[0])
v.write_pose('minimized.pdbqt', overwrite=True)
 
# Dock the ligand
v.dock(exhaustiveness=32, n_poses=20)
v.write_poses('ligand_vina_out.pdbqt', n_poses=10, overwrite=True)
