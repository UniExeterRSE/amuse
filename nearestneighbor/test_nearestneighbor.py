from amuse.datamodel import Particles
from amuse.test.amusetest import TestWithMPI
from amuse.units import generic_unit_system

from .interface import NearestNeighborInterface
from .interface import NearestNeighbor

class NearestNeighborInterfaceTests(TestWithMPI):
    
    def test1(self):
        instance = NearestNeighborInterface()
        result,error = instance.new_particle(1.0, 1.0, 2.0)
        self.assertEquals(error, 0)
        self.assertEquals(result, 0)
        instance.stop()
        
    def test2(self):
        instance = NearestNeighborInterface()
        result,error = instance.new_particle(1.0, 1.0, 2.0)
        self.assertEquals(error, 0)
        self.assertEquals(result, 0)
        result,error = instance.new_particle(2.0, 3.0, 2.0)
        self.assertEquals(error, 0)
        self.assertEquals(result, 1)
        #result,error = instance.new_particle(2.0, 3.0, 2.0)
        #self.assertEquals(error, -1)
        error = instance.delete_particle(1)
        self.assertEquals(error, 0)
        result,error = instance.new_particle(2.0, 3.0, 2.0)
        self.assertEquals(error, 0)
        self.assertEquals(result, 1)
        instance.stop()
        
    def test3(self):
        instance = NearestNeighbor()
        #instance.set_maximum_number_of_particles(2)
        #instance.commit_parameters()
        result = instance.new_particle(
            1.0 | generic_unit_system.length,
            2.0 | generic_unit_system.length,
            3.0 | generic_unit_system.length
        )
        self.assertEquals(result, 0)
        result = instance.new_particle(
            1.0 | generic_unit_system.length,
            1.0 | generic_unit_system.length,
            2.0 | generic_unit_system.length
        )
        self.assertEquals(result, 1)
        x,y,z = instance.get_position(0)
        self.assertEquals(1.0 | generic_unit_system.length, x)
        self.assertEquals(2.0 | generic_unit_system.length, y)
        self.assertEquals(3.0 | generic_unit_system.length, z)
        instance.stop()
        
    def test4(self):
        instance = NearestNeighbor()
        #instance.set_maximum_number_of_particles(100)
        #instance.commit_parameters()

        particles = Particles(4)
        particles.x = [0.0, 1.0, 4.0, 7.5] | generic_unit_system.length
        particles.y = 0.0 | generic_unit_system.length
        particles.z = 0.0 | generic_unit_system.length

        instance.particles.add_particles(particles)
        instance.run()

        self.assertEqual(instance.particles[0].neighbor0, instance.particles[1])
        self.assertEqual(instance.particles[1].neighbor0, instance.particles[0])
        self.assertEqual(instance.particles[2].neighbor0, instance.particles[1])
        self.assertEqual(instance.particles[3].neighbor0, instance.particles[2])

        instance.stop()

