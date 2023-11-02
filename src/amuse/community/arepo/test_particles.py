from amuse.ic.plummer import new_plummer_model
from amuse.community.gadget2 import Gadget2
from amuse.community.arepo import Arepo
from amuse.units import nbody_system

particles = new_plummer_model(1000)

instance_gg = Gadget2(redirection="none")
instance_gg.dm_particles.add_particles(particles)
instance_gg.evolve_model(0.01 | nbody_system.time)
print(f"gadget time: {instance_gg.model_time}")

instance_ar = Arepo(redirection="none")
instance_ar.dm_particles.add_particles(particles)
instance_ar.evolve_model(0.01 | nbody_system.time)
print(f"arepo time: {instance_ar.model_time}")

# then do some comparison between instance_gg.dm_particles and instance_ar.dm_particles
