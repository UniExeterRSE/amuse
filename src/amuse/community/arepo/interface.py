from amuse.community import CodeInterface
from amuse.community import LegacyFunctionSpecification
from amuse.community import legacy_function
from amuse.community import LiteratureReferencesMixIn

from amuse.community.interface.gd import GravitationalDynamicsInterface
from amuse.community.interface.gd import GravitationalDynamics

from amuse.units import generic_unit_system


class ArepoInterface(
    CodeInterface,
    GravitationalDynamicsInterface,
    LiteratureReferencesMixIn
):
    """
    Arepo is a cosmological magnetohydrodynamical moving-mesh simulation code,
    descended from GADGET.

    References:
        .. [#] Springel, V., 2010, MNRAS, 401, 791 (Arepo) [2010MNRAS.401..791S]
        .. [#] Pakmor, R., Bauer, A., Springel, V., 2011, MNRAS, 418, 1392 (Magnetohydrodynamics Module) [2011MNRAS.418.1392P]
        .. [#] Pakmor, R. et al., 2016, MNRAS, 455, 1134 (Gradient Estimation) [2016MNRAS.455.1134P]
        .. [#] Weinberger, R., Springel, V., Pakmor, R., 2020, ApJS, 248, 32 (Public Code Release) [2020ApJS..248...32W]
    """

    include_headers = ["worker_code.h", "interface.h"]

    def __init__(self, **keyword_arguments):
        CodeInterface.__init__(
            self, name_of_the_worker="arepo_worker", **keyword_arguments
        )
        LiteratureReferencesMixIn.__init__(self)
        # TODO: Determine whether need to inherit from CodeWithDataDirectories.

    @legacy_function
    def get_minimum_time_step():
        """
        Retrieve the model timestep.
        """
        function = LegacyFunctionSpecification()
        function.addParameter(
            'minimum_time_step', dtype='float64', direction=function.OUT,
            description="The minimum timestep"
        )
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_minimum_time_step():
        """
        Retrieve the model timestep.
        """
        function = LegacyFunctionSpecification()
        function.addParameter(
            'minimum_time_step', dtype='float64', direction=function.IN,
            description="The minimum timestep"
        )
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_pressure():
        function = LegacyFunctionSpecification()
        function.addParameter(
            "index_of_the_particle", dtype="int32", direction=function.IN
        )
        function.addParameter("p", dtype="float64", direction=function.OUT)
        function.result_type = "int32"
        function.can_handle_array = True
        return function

    @legacy_function
    def get_density():
        function = LegacyFunctionSpecification()
        function.addParameter(
            "index_of_the_particle", dtype="int32", direction=function.IN
        )
        function.addParameter("rho", dtype="float64", direction=function.OUT)
        function.result_type = "int32"
        function.can_handle_array = True
        return function

    @legacy_function
    def get_internal_energy():
        function = LegacyFunctionSpecification()
        function.addParameter(
            "index_of_the_particle", dtype="int32", direction=function.IN
        )
        function.addParameter("u", dtype="float64", direction=function.OUT)
        function.result_type = "int32"
        function.can_handle_array = True
        return function
    
    @legacy_function
    def set_internal_energy():
        function = LegacyFunctionSpecification()
        function.addParameter(
            "index_of_the_particle", dtype="int32", direction=function.IN)
        function.addParameter("u", dtype="float64", direction=function.IN)
        function.result_type = "int32"
        function.can_handle_array = True        
        return function

    @legacy_function
    def new_dm_particle():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            "index_of_the_particle", dtype="int32", direction=function.OUT
        )
        for x in ["mass", "x", "y", "z", "vx", "vy", "vz"]:
            function.addParameter(x, dtype="float64", direction=function.IN)
        function.result_type = "int32"
        return function

    @legacy_function
    def new_gas_particle():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            "index_of_the_particle", dtype="int32", direction=function.OUT
        )
        for x in ["mass", "x", "y", "z", "vx", "vy", "vz", "u"]:
            function.addParameter(x, dtype="float64", direction=function.IN)
        function.result_type = "int32"
        return function

    def new_particle(self, mass, x, y, z, vx, vy, vz):
        return self.new_dm_particle(mass, x, y, z, vx, vy, vz)

    @legacy_function
    def get_state_gas():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_particle', dtype='int32', direction=function.IN,
            description=(
                "Index of the particle to get the state from. "
                "This index must have been returned by an earlier call to "
                ":meth:`new_particle`"
            )
        )
        for x in ['mass','x','y','z','vx','vy','vz','u']:
            function.addParameter(x, dtype='float64', direction=function.OUT)
        # function.addParameter('length', 'int32', function.LENGTH)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_state_gas():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter(
            'index_of_the_particle', dtype='int32', direction=function.IN,
            description=(
                "Index of the particle for which the state is to be updated. "
                "This index must have been returned by an earlier call to "
                ":meth:`new_particle`"
            )
        )
        for x in ['mass','x','y','z','vx','vy','vz','u']:
            function.addParameter(x, dtype='float64', direction=function.IN)
        # function.addParameter('length', 'int32', function.LENGTH)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_box_size():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function
        
    @legacy_function
    def get_box_size():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function
    
    @legacy_function
    def get_omega_zero():
        function = LegacyFunctionSpecification()
        function.addParameter('omega_zero', dtype='d', direction=function.OUT,
            description = "Cosmological matter density parameter in units of the critical density at z=0.")
        function.result_type = 'i'
        return function
    @legacy_function
    def set_omega_zero():
        function = LegacyFunctionSpecification()
        function.addParameter('omega_zero', dtype='d', direction=function.IN,
            description = "Cosmological matter density parameter in units of the critical density at z=0.")
        function.result_type = 'i'
        return function
    
    @legacy_function
    def get_omega_lambda():
        function = LegacyFunctionSpecification()
        function.addParameter('omega_lambda', dtype='d', direction=function.OUT,
            description = "Cosmological vacuum energy density parameter in units of the critical density at z=0.")
        function.result_type = 'i'
        return function
    @legacy_function
    def set_omega_lambda():
        function = LegacyFunctionSpecification()
        function.addParameter('omega_lambda', dtype='d', direction=function.IN,
            description = "Cosmological vacuum energy density parameter in units of the critical density at z=0.")
        function.result_type = 'i'
        return function
    
    @legacy_function
    def get_omega_baryon():
        function = LegacyFunctionSpecification()
        function.addParameter('omega_baryon', dtype='d', direction=function.OUT,
            description = "Cosmological baryonic density parameter in units of the critical density at z=0.")
        function.result_type = 'i'
        return function
    @legacy_function
    def set_omega_baryon():
        function = LegacyFunctionSpecification()
        function.addParameter('omega_baryon', dtype='d', direction=function.IN,
            description = "Cosmological baryonic density parameter in units of the critical density at z=0.")
        function.result_type = 'i'
        return function
    
    @legacy_function
    def get_hubble_param():
        function = LegacyFunctionSpecification()
        function.addParameter('hubble_param', dtype='d', direction=function.OUT,
            description = "The cosmological Hubble parameter.")
        function.result_type = 'i'
        return function
    @legacy_function
    def set_hubble_param():
        function = LegacyFunctionSpecification()
        function.addParameter('hubble_param', dtype='d', direction=function.IN,
            description = "The cosmological Hubble parameter.")
        function.result_type = 'i'
        return function

    # This function has been kept as a basic template for future functions.
    # @legacy_function
    # def set_parameters():
    #     function = LegacyFunctionSpecification()
    #     function.addParameter("param_file", dtype="string", direction=function.IN)
    #     function.result_type = "int32"
    #     return function


class Arepo(GravitationalDynamics):
    def __init__(self, unit_converter=None, **options):
        GravitationalDynamics.__init__(
            self, ArepoInterface(**options), unit_converter, **options
        )

    def initialize_code(self):
        result = self.overridden().initialize_code()

        return result

    def define_methods(self, handler):
        GravitationalDynamics.define_methods(self, handler)
        # TODO: Determine how to link this to Arepo's run() - the main simulation loop.
        handler.add_method("run_sim", (), (handler.ERROR_CODE))
        # When simulation is finished, shutdown HDF5 & MPI, and exit(0)
        handler.add_method("cleanup_code", (), (handler.ERROR_CODE))
        handler.add_method(
            "new_dm_particle",
            (
                generic_unit_system.mass,
                generic_unit_system.length,
                generic_unit_system.length,
                generic_unit_system.length,
                generic_unit_system.speed,
                generic_unit_system.speed,
                generic_unit_system.speed,
            ),
            (
                handler.INDEX,
                handler.ERROR_CODE,
            ),
        )
        handler.add_method(
            "new_gas_particle",
            (
                generic_unit_system.mass,
                generic_unit_system.length,
                generic_unit_system.length,
                generic_unit_system.length,
                generic_unit_system.speed,
                generic_unit_system.speed,
                generic_unit_system.speed,
                generic_unit_system.specific_energy,
            ),
            (
                handler.INDEX,
                handler.ERROR_CODE,
            )
        )
        handler.add_method(
            "get_state_gas",
            (
                handler.INDEX,
            ),
            (
                generic_unit_system.mass,
                generic_unit_system.length,
                generic_unit_system.length,
                generic_unit_system.length,
                generic_unit_system.speed,
                generic_unit_system.speed,
                generic_unit_system.speed,
                generic_unit_system.specific_energy,
                handler.ERROR_CODE
            )
        )
        handler.add_method(
            "set_state_gas",
            (
                handler.INDEX,
                generic_unit_system.mass,
                generic_unit_system.length,
                generic_unit_system.length,
                generic_unit_system.length,
                generic_unit_system.speed,
                generic_unit_system.speed,
                generic_unit_system.speed,
                generic_unit_system.specific_energy,
            ),
            (
                handler.ERROR_CODE,
            )
        )

        handler.add_method(
            "get_box_size",
            (),
            (generic_unit_system.length, handler.ERROR_CODE,)
        )
        
        handler.add_method(
            "set_box_size",
            (generic_unit_system.length, ),
            (handler.ERROR_CODE,)
        )

        handler.add_method(
            "get_minimum_time_step",
            (),
            (generic_unit_system.time, handler.ERROR_CODE,)
        )

        handler.add_method(
            "set_minimum_time_step",
            (generic_unit_system.time),
            (handler.ERROR_CODE,)
        )
        

    def define_particle_sets(self, handler):
        handler.define_super_set('particles', ['dm_particles','gas_particles'], 
            index_to_default_set = 0)
        handler.define_set('dm_particles', 'index_of_the_particle')
        handler.set_new('dm_particles', 'new_dm_particle')
        handler.set_delete('dm_particles', 'delete_particle')
        handler.add_setter('dm_particles', 'set_state')
        handler.add_getter('dm_particles', 'get_state')
        handler.add_setter('dm_particles', 'set_mass')
        handler.add_getter('dm_particles', 'get_mass', names=('mass',))
        handler.add_setter('dm_particles', 'set_position')
        handler.add_getter('dm_particles', 'get_position')
        handler.add_setter('dm_particles', 'set_velocity')
        handler.add_getter('dm_particles', 'get_velocity')
        # handler.add_getter('dm_particles', 'get_acceleration')
        # handler.add_getter('dm_particles', 'get_epsilon_dm_part', names = ('radius',))
        # handler.add_getter('dm_particles', 'get_epsilon_dm_part', names = ('epsilon',))

        handler.define_set('gas_particles', 'index_of_the_particle')
        handler.set_new('gas_particles', 'new_gas_particle')
        handler.set_delete('gas_particles', 'delete_particle')
        handler.add_setter('gas_particles', 'set_state_gas')
        handler.add_getter('gas_particles', 'get_state_gas')
        handler.add_setter('gas_particles', 'set_mass')
        handler.add_getter('gas_particles', 'get_mass', names=('mass',))
        handler.add_setter('gas_particles', 'set_position')
        handler.add_getter('gas_particles', 'get_position')
        handler.add_setter('gas_particles', 'set_velocity')
        handler.add_getter('gas_particles', 'get_velocity')
        # handler.add_getter('gas_particles', 'get_acceleration')
        handler.add_setter('gas_particles', 'set_internal_energy')
        handler.add_getter('gas_particles', 'get_internal_energy')
        # handler.add_getter('gas_particles', 'get_smoothing_length')
        # handler.add_getter('gas_particles', 'get_density', names=('rho',))
        # handler.add_getter('gas_particles', 'get_density', names=('density',))
        # handler.add_getter('gas_particles', 'get_pressure')
        # handler.add_getter('gas_particles', 'get_d_internal_energy_dt')
        # handler.add_getter('gas_particles', 'get_n_neighbours')
        # handler.add_getter('gas_particles', 'get_epsilon_gas_part', names = ('radius',))
        # handler.add_getter('gas_particles', 'get_epsilon_gas_part', names = ('epsilon',))


    def define_parameters(self, handler):
        handler.add_method_parameter(
            "get_time_step", 
            None,
            "timestep", 
            "timestep for the system.", 
            default_value=1.0 | generic_unit_system.time
        )

        handler.add_method_parameter(
            "get_minimum_time_step", 
            "set_minimum_time_step",
            "minimum_timestep", 
            "minimum timestep for the system.", 
            default_value=1.0e-5 | generic_unit_system.time
        )

        handler.add_method_parameter(
            "get_omega_zero", 
            "set_omega_zero",
            "omega_zero", 
            "Cosmological matter density parameter in units of the critical density at z=0.", 
            default_value = 0.0
        )
        
        handler.add_method_parameter(
            "get_omega_lambda", 
            "set_omega_lambda",
            "omega_lambda", 
            "Cosmological vacuum energy density parameter in units of the critical density at z=0.", 
            default_value = 0.0
        )
        
        handler.add_method_parameter(
            "get_omega_baryon", 
            "set_omega_baryon",
            "omega_baryon", 
            "Cosmological baryonic density parameter in units of the critical density at z=0.", 
            default_value = 0.0
        )
        
        handler.add_method_parameter(
            "get_hubble_param", 
            "set_hubble_param",
            "hubble_parameter", 
            "The cosmological Hubble parameter, value of Hubble constant in units of 100 km/s / Mpc.", 
            default_value = 0.7
        )
        handler.add_method_parameter(
            "get_box_size", 
            "set_box_size",
            "periodic_box_size", 
            "The size of the box in case of periodic boundary conditions.", 
            default_value = 1.0 | generic_unit_system.length
        )

    def define_state(self, handler):
        GravitationalDynamics.define_state(self, handler)

        handler.add_method("EDIT", "new_dm_particle")
        handler.add_method("EDIT", "new_gas_particle")
        handler.add_method("UPDATE", "new_dm_particle")
        handler.add_method("UPDATE", "new_gas_particle")
        handler.add_transition("RUN", "UPDATE", "new_dm_particle", False)
        handler.add_transition("RUN", "UPDATE", "new_gas_particle", False)

        handler.add_method('RUN', 'get_state_dm')
        handler.add_method('RUN', 'get_state_gas')

