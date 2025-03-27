"""
Noh3d example, taken from arepo examples.

Prior to running this script, perform the following steps:

1. Run `configToArepoflags.sh` against `configs/noh_config.sh` to generate `Arepoflags.mk`.
2. Run `make clean && make`.
3. Move file `params/noh3d_params.txt` to `params.txt`.
"""

import random
import sys

from matplotlib import pyplot as plt
from matplotlib import cm
import numpy as np

from amuse.community.arepo import Arepo
from amuse.datamodel import Particles
from amuse.units import nbody_system


N_PLOT_PARTICLES = 3000

def plot_particles_3d(
    plot_count, instance, tracked_ids,
    vmin=0.1, vmax=100,
):
    time = instance.get_time()
    x, y, z = instance.get_position(tracked_ids)
    dens = instance.get_density(tracked_ids)

    max_density = max(dens)
    min_density = min(dens)
    print(max_density, min_density)
    fig, ax = plt.subplots(subplot_kw=dict(projection='3d',), figsize=(10, 7))
    cmap = cm.plasma
    plot = ax.scatter(
        x.number, y.number, z.number,
        c=dens,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        edgecolors='none',
    )
    plt.colorbar(plot, ax=ax)

    ax.set_xlim((0, 6))
    ax.set_ylim((0, 6))
    ax.set_zlim((0, 6))
    fig.suptitle(f"Densities at t = {time.number:04.2f}")
    fig.savefig(f"noh_positions_{plot_count:04d}.png")
    plt.close(fig)
    return plot_count + 1

def noh_3d_data():
    """
    Noh 3d example data.
    
    Adapted from Arepo example script.
    """

    FloatType = np.float64
    IntType = np.int32

    Boxsize = 6.0
    CellsPerDimension = 30
    NumberOfCells = CellsPerDimension * CellsPerDimension * CellsPerDimension

    ## initial state
    density_0 = 1.0
    velocity_radial_0 = -1.0    ## radial inflow velocity
    pressure_0 = 1.0e-4
    gamma = 5./3.  ## note: this has to be consistent with the parameter settings for Arepo!
    utherm_0 = pressure_0 / ( gamma - 1.0 ) / density_0


    """ set up grid: cartesian 3d grid """
    ## spacing
    dx = Boxsize / FloatType(CellsPerDimension)
    ## position of first and last cell
    pos_first, pos_last = 0.5 * dx, Boxsize - 0.5 * dx


    ## set up grid
    Grid1d = np.linspace(pos_first, pos_last, CellsPerDimension, dtype=FloatType)
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

    """ set up hydrodynamical quantitites """
    ## mass insetad of density
    Mass = np.full(NumberOfCells, density_0*dx*dx*dx, dtype=FloatType)
    ## velocity
    Velocity = np.zeros([NumberOfCells,3], dtype=FloatType)
    Velocity[:,0] = velocity_radial_0 * xPosFromCenter / Radius
    Velocity[:,1] = velocity_radial_0 * yPosFromCenter / Radius
    Velocity[:,2] = velocity_radial_0 * zPosFromCenter / Radius
    ## specific internal energy
    Uthermal = np.full(NumberOfCells, utherm_0, dtype=FloatType)

    ParticleIDs = np.arange(1, NumberOfCells+1)

    return {"ParticleIDs": ParticleIDs,
            "Coordinates": Pos,
            "Masses": Mass,
            "Velocities": Velocity,
            "InternalEnergy": Uthermal,
            }


def evolve(instance, target_time):
    instance.evolve_model(target_time)
    time = instance.get_time()
    print(f"Evolving to time = {time} (offset from requested: {time - target_time})")


def main():
    # Check code runs without errors
    # instance = Arepo(redirection="none")
    instance = Arepo(redirection="none")
    instance.initialize_code()
    print(instance.get_box_size())

    d = noh_3d_data()

    p = Particles(len(d["Masses"]))
    p.mass = d["Masses"] | nbody_system.mass
    p.x = d["Coordinates"][:,0] | nbody_system.length
    p.y = d["Coordinates"][:,1] | nbody_system.length
    p.z = d["Coordinates"][:,2] | nbody_system.length
    p.vx = d["Velocities"][:,0] | nbody_system.speed
    p.vy = d["Velocities"][:,1] | nbody_system.speed
    p.vz = d["Velocities"][:,2] | nbody_system.speed
    p.u = d["InternalEnergy"] | nbody_system.specific_energy

    print("AMUSE: BoxSize:")
    print(instance.get_box_size())
    instance.gas_particles.add_particles(p)

    print(instance.gas_particles)
    
    n_particles_total = instance.get_number_of_particles()
    print(f'AMUSE: number of particles: {n_particles_total}')
    random.seed(123)
    tracked_ids = random.sample(range(n_particles_total), k=N_PLOT_PARTICLES)
    
    plot_count = 0
    for t in np.linspace(0.1, 2, 190):
        evolve(instance, t | nbody_system.time)
        plot_count = plot_particles_3d(plot_count, instance, tracked_ids)
    sys.exit()


if __name__ == "__main__":
    main()
