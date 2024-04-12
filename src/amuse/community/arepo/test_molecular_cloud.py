"""
Test Arepo with molecular cloud ICs from AMUSE
"""

""" load libraries """
import sys    # system specific calls
import numpy as np    # scientific computing package
# import h5py    # hdf5 format
from amuse.units import nbody_system, units
from amuse.datamodel import Particles
from amuse.io import write_set_to_file

from amuse.ic.molecular_cloud import new_molecular_cloud

from amuse.community.arepo import Arepo

np.random.seed(123)

converter = nbody_system.nbody_to_si(1000.0 | units.MSun, 2.0 | units.RSun)
# instance = Arepo(converter, redirection="none")
instance = Arepo(redirection="none")
instance.parameters.periodic_box_size = 10 | nbody_system.length

cloud = new_molecular_cloud(targetN=20000)  #, convert_nbody=converter)

# move particles to the center of the box (not yet sure if this box is periodic)
cloud.x += (0.5 * instance.parameters.periodic_box_size)
cloud.y += (0.5 * instance.parameters.periodic_box_size)
cloud.z += (0.5 * instance.parameters.periodic_box_size)

instance.gas_particles.add_particles(cloud)

write_set_to_file(instance.gas_particles, "cloudorig.amuse", overwrite_file=True)

print(len(cloud))
print(len(instance.gas_particles))

# print(instance.parameters)

# Modified from arepo/examples/noh3d/create.py
NumPart = np.array([NumberOfCells, 0, 0, 0, 0, 0], dtype=IntType)
params_to_set = {
    # "NumPart_ThisFile": NumPart,
    # "NumPart_Total": NumPart,
    # "NumPart_Total_HighWord": np.zeros(6, dtype=IntType),
    # "MassTable": np.zeros(6, dtype=IntType),
    # "Time": 0.0,
    # "Redshift": 0.0,
    # "BoxSize": 10*Boxsize,
    # "NumFilesPerSnapshot": 1,
    # "Omega0": 0.0,
    # "OmegaB": 0.0,
    # "OmegaLambda": 0.0,
    # "HubbleParam": 1.0,
    # "Flag_Sfr": 0,
    # "Flag_Cooling": 0,
    # "Flag_StellarAge": 0,
    # "Flag_Metals": 0,
    # "Flag_Feedback": 0,
}
# if Pos.dtype == np.float64:
#     params_to_set["Flag_DoublePrecision"] = 1
# else:
#     params_to_set["Flag_DoublePrecision"] = 0

# for k, v in params_to_set.items():
#     try:
#         instance.parameters.__setattr__(k, v)
#     except:
#         print(f"Could not set {k} to {v}")

# instance.parameters.minimum_timestep = 0.1 * instance.parameters.minimum_timestep
# instance.commit_particles()
print(instance.parameters)
print(instance.get_number_of_particles())
instance.evolve_model(1.9 | nbody_system.time)

write_set_to_file(instance.gas_particles, "cloud_end.amuse", overwrite_file=True)
