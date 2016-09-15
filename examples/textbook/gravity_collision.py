import numpy
from amuse.lab import *
from amuse.io import store
from optparse import OptionParser
from amuse.community.seba.interface import SeBa
from amuse.ext.LagrangianRadii import LagrangianRadii

def merge_two_stars(bodies, particles_in_encounter, tcoll):
    com_pos = particles_in_encounter.center_of_mass()
    com_vel = particles_in_encounter.center_of_mass_velocity()
    new_particle=Particles(1)
    new_particle.mass = particles_in_encounter.total_mass()
    new_particle.collision_time = tcoll
    new_particle.position = com_pos
    new_particle.velocity = com_vel
    new_particle.radius = particles_in_encounter.radius.sum()/len(particles_in_encounter.radius)
    bodies.add_particles(new_particle)
    bodies.remove_particles(particles_in_encounter)
    return new_particle
    
def main(N=10, W0=7.0, t_end=10, nsteps=10, filename="gravity_stellar.hdf5", Mmax= 100, Qvir=0.5):
    t_end = t_end | nbody_system.time
    dt = t_end/float(nsteps)

    bodies = new_king_model(N, W0)
    masses = new_salpeter_mass_distribution(N, mass_min=0.1|units.MSun, mass_max=10|units.MSun)
#    masses = masses.value_in(units.MSun)
    masses /= masses.sum()
    bodies.mass = masses | nbody_system.mass

    """
    alpha1 = 2.35 + 1
    factor = (pow(Mmax, alpha1) - 1.0 )
    masses = (pow(1.0 + (factor * numpy.random.random(N)), 1.0 / alpha1))
    masses = masses / masses.sum()
    Mtot_init = masses.sum()
    """
    bodies.radius = 0.001 | nbody_system.length
    Mtot_init = 1 | nbody_system.mass

    gravity = ph4(number_of_workers=4)
#    gravity.parameters.timestep_parameter = 0.01
    gravity.particles.add_particles(bodies)
    bodies.scale_to_standard(virial_ratio=Qvir)

    channel_from_gd_to_framework = gravity.particles.new_channel_to(bodies)

    stopping_condition = gravity.stopping_conditions.collision_detection
    stopping_condition.enable()
    
    Etot_init = gravity.kinetic_energy + gravity.potential_energy
    Etot_prev = Etot_init

    time = 0.0 | t_end.unit
    Nenc = 0
    dEk_enc = zero    
    dEp_enc = zero
    while time < t_end:

        RL = LagrangianRadii(gravity.particles)
        pos,coreradius,coredens=gravity.particles.densitycentre_coreradius_coredens()
        print "Cluster at time=", time, "core_radius=", coreradius, "L_radius=",
        for rl in RL:
            print rl.number, " ",
        print "length"

        gravity.evolve_model(time+dt)
        time = gravity.model_time

        if stopping_condition.is_set():
            Ek_enc = gravity.kinetic_energy 
            Ep_enc = gravity.potential_energy
            print "At time=", time, "number of encounters=", len(stopping_condition.particles(0))
            for ci in range(len(stopping_condition.particles(0))): 
                particles_in_encounter = Particles(particles=[stopping_condition.particles(0)[ci], stopping_condition.particles(1)[ci]])
                particles_in_encounter = particles_in_encounter.get_intersecting_subset_in(bodies)

                new_particle = merge_two_stars(bodies, particles_in_encounter, gravity.model_time)
                bodies.synchronize_to(gravity.particles)
                Nenc+=1
                RL = LagrangianRadii(gravity.particles)

                print "Resolve encounter Number:", Nenc
                pos,coreradius,coredens=gravity.particles.densitycentre_coreradius_coredens()
                print "Collision at time=", time, new_particle.mass.sum(), new_particle.position.length(), "Nstars= ", len(bodies), "Ncoll=", Nenc, "core_radius=", coreradius, "L_radius=", 
                for rl in RL:
                    print rl.number, " ",
                print "length"

            dEk_enc += Ek_enc - gravity.kinetic_energy 
            dEp_enc += Ep_enc - gravity.potential_energy

        channel_from_gd_to_framework.copy()

#        write_set_to_file(bodies.savepoint(time), filename, 'hdf5')

        Ekin = gravity.kinetic_energy 
        Epot = gravity.potential_energy
        Etot = Ekin + Epot
        dE = Etot_prev-Etot
        Mtot = bodies.mass.sum()
        print "T=", time, 
        print "M=", Mtot, "(dM[SE]=", Mtot/Mtot_init, ")",
        print "E= ", Etot, "Q= ", Ekin/Epot,
        print "dE=", (Etot_init-Etot)/Etot, "ddE=", (Etot_prev-Etot)/Etot, 
        print "dE(enc)=", dEk_enc, dEp_enc
        Etot_init -= dE
        Etot_prev = Etot


    gravity.stop()
    
def new_option_parser():
    result = OptionParser()
    result.add_option("-f", dest="filename", default = "gravity_stellar.hdf5",
                      help="output filename [gravity_stellar.hdf5]")
    result.add_option("-N", dest="N", type="int",default = 1000,
                      help="number of stars [100]")
    result.add_option("-n", dest="nsteps", type="int",default = 100,
                      help="number of output steps [10]")
    result.add_option("-M", dest="Mmax", type="float",default = 100,
                      help="maximal stellar mass [100]")
    result.add_option("-t", dest="t_end", type="float", default = 100,
                      help="end time of the simulation [1]")
    result.add_option("-Q", dest="Qvir", type="float", default = 0.5,
                      help="initial virial radio [%default]")
    result.add_option("-W", dest="W0", type="float", default = 7.0,
                      help="Dimension-less depth of the King potential (W0) [7.0]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    main(**o.__dict__)

