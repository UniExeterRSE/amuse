""" @package ./examples/Noh_3d/create.py
Code that creates 3d Noh test problem initial conditions

created by Rainer Weinberger, last modified 24.02.2019
modified for AMUSE by Steven Rieder, 17 November 2023
"""

""" load libraries """
import sys    # system specific calls
import numpy as np    # scientific computing package
# import h5py    # hdf5 format
from amuse.units import nbody_system
from amuse.datamodel import Particles
from amuse.io import write_set_to_file

from amuse.community.arepo import Arepo

simulation_directory = str(sys.argv[1])
print("examples/Noh_3d/create.py: creating ICs in directory " + simulation_directory)

""" initial condition parameters """
FilePath = f"{simulation_directory}/IC.amuse"

FloatType = np.float64  # double precision: np.float64, for single use np.float32
IntType = np.int32

Boxsize = FloatType(6.0)
CellsPerDimension = IntType(30)
NumberOfCells = CellsPerDimension * CellsPerDimension * CellsPerDimension

## initial state
density_0 = 1.0
velocity_radial_0 = -1.0 ## radial inflow velocity
pressure_0 = 1.0e-4
gamma = 5./3.  ## note: this has to be consistent with the parameter settings for Arepo!
utherm_0 = pressure_0 / ( gamma - 1.0 ) / density_0

""" set up grid: cartesian 3d grid """
## spacing
dx = Boxsize / FloatType(CellsPerDimension)
## position of first and last cell
pos_first, pos_last = 0.5 * dx, Boxsize - 0.5 * dx

## set up grid
Grid1d = np.linspace(
    pos_first,
    pos_last,
    CellsPerDimension,
    dtype=FloatType,
)
xx, yy, zz = np.meshgrid(Grid1d, Grid1d, Grid1d)
Pos = np.zeros([NumberOfCells, 3], dtype=FloatType)
Pos[:,0] = xx.reshape(NumberOfCells)
Pos[:,1] = yy.reshape(NumberOfCells)
Pos[:,2] = zz.reshape(NumberOfCells)
## calculate distance from center
xPosFromCenter = (Pos[:,0] - 0.5 * Boxsize)
yPosFromCenter = (Pos[:,1] - 0.5 * Boxsize)
zPosFromCenter = (Pos[:,2] - 0.5 * Boxsize)
Radius = np.sqrt( xPosFromCenter**2 + yPosFromCenter**2 + zPosFromCenter**2 )

""" set up hydrodynamical quantities """
## mass instead of density
Mass = np.full(NumberOfCells, density_0*dx*dx*dx, dtype=FloatType)
## velocity
Velocity = np.zeros([NumberOfCells,3], dtype=FloatType)
Velocity[:,0] = velocity_radial_0 * xPosFromCenter / Radius
Velocity[:,1] = velocity_radial_0 * yPosFromCenter / Radius
Velocity[:,2] = velocity_radial_0 * zPosFromCenter / Radius
## specific internal energy
Uthermal = np.full(NumberOfCells, utherm_0, dtype=FloatType)

p = Particles(len(Pos))
p.position = Pos[:,:] | nbody_system.length
p.velocity = Velocity[:,:] | nbody_system.speed
p.mass = Mass | nbody_system.mass
p.u = Uthermal | nbody_system.speed**2

write_set_to_file(p, f'{simulation_directory}/IC.amuse')

instance = Arepo(redirection="none")
#instance.parameters
