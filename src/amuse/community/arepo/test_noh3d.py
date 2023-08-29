import sys
import random
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm
from amuse.community.arepo import Arepo
from amuse.datamodel import Particles, Particle
# import unit seconds as s

N_PLOT_PARTICLES = 3000

def plot_particles_3d(
    plot_count, instance, tracked_ids,
    vmin=0.1, vmax=100,
):
    time = instance.get_time()
    x, y, z = instance.get_position(tracked_ids)
    dens = instance.get_density(tracked_ids)

    # for i in range(len(x)):
    #     print(x[i], y[i], z[i], dens[i])
    # sys.exit()
    # positions = {
    #     index: [instance.get_position(index)] for index in tracked_ids
    # }
    # densities = {
    #     index: [instance.get_density(index)] for index in tracked_ids
    # }
    max_density = max(dens)  # densities.values())[0]
    min_density = min(dens)  # densities.values())[0]
    print(max_density, min_density)
    fig, ax = plt.subplots(subplot_kw=dict(projection='3d',), figsize=(10, 7))
    cmap = cm.plasma
    plot = ax.scatter(
        x, y, z,
        c=dens,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        edgecolors='none',
        s=5,
    )
    plt.colorbar(plot, ax=ax)
    # for index, pos in positions.items():
    #     ax.plot(
    #         *np.array(pos).T,
    #         color=cmap(
    #             (densities[index] - min_density)
    #             / (max_density - min_density)
    #         )
    #     )
    ax.set_xlim((0, 6))
    ax.set_ylim((0, 6))
    ax.set_zlim((0, 6))
    fig.suptitle(f"Densities at t = {time:04.2f}")
    fig.savefig(f"noh_positions_{plot_count:04d}.png")
    plt.close(fig)
    return plot_count + 1


def evolve(instance, target_time):
    instance.evolve_model(target_time)
    time = instance.get_time()
    print(f"Evolving to time = {time:.4f} (offset from requested: {time - target_time:.4f})")


def main():
    # Check code runs without errors
    # instance = Arepo(redirection="none")
    instance = Arepo(redirection="null")
    instance.initialize_code()
    
    n_particles_total = instance.get_number_of_particles()
    print(f'AMUSE: number of particles: {n_particles_total}')
    random.seed(123)
    tracked_ids = random.sample(range(n_particles_total), k=N_PLOT_PARTICLES)
    
    plot_count = 0
    for t in np.linspace(0.1, 2, 240):
        evolve(instance, t)
        plot_count = plot_particles_3d(plot_count, instance, tracked_ids)
    sys.exit()


if __name__ == "__main__":
    main()
