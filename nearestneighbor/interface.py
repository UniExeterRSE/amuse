from amuse.community import *

class NearestNeighborInterface(CodeInterface):
    
    include_headers = ['worker_code.h']
    
    def __init__(self, **keyword_arguments):
        CodeInterface.__init__(self, name_of_the_worker="nearestneighbor_worker", **keyword_arguments)
    
    @legacy_function
    def new_particle():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.OUT)
        function.addParameter('x', dtype='float64', direction=function.IN)
        function.addParameter('y', dtype='float64', direction=function.IN)
        function.addParameter('z', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def delete_particle():
        function = LegacyFunctionSpecification()  
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_state():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True 
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN)
        function.addParameter('x', dtype='float64', direction=function.OUT)
        function.addParameter('y', dtype='float64', direction=function.OUT)
        function.addParameter('z', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def set_state():
        function = LegacyFunctionSpecification()  
        function.can_handle_array = True 
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN)
        function.addParameter('x', dtype='float64', direction=function.IN)
        function.addParameter('y', dtype='float64', direction=function.IN)
        function.addParameter('z', dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function    

    @legacy_function
    def find_nearest_neighbors():
        function = LegacyFunctionSpecification()
        function.result_type = 'int32'
        return function    
    
    @legacy_function
    def get_close_neighbors():
        function = LegacyFunctionSpecification()  
        function.can_handle_array = True 
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN)
        function.addParameter('index_of_first_neighbor', dtype='float64', direction=function.OUT)
        function.addParameter('index_of_second_neighbor', dtype='float64', direction=function.OUT)
        function.addParameter('index_of_third_neighbor', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function    
    
    @legacy_function
    def get_nearest_neighbor():
        function = LegacyFunctionSpecification()  
        function.can_handle_array = True 
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN)
        function.addParameter('index_of_the_neighbor', dtype='float64', direction=function.OUT)
        function.addParameter('distance', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def get_number_of_particles():
        function = LegacyFunctionSpecification()  
        function.can_handle_array = True 
        function.addParameter('value', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function    
    
    
class NearestNeighbor(InCodeComponentImplementation):

    def __init__(self, **options):
        InCodeComponentImplementation.__init__(self,  NearestNeighborInterface(), **options)


    def define_methods(self, builder):

        builder.add_method(
            "new_particle",
            (generic_unit_system.length, generic_unit_system.length, generic_unit_system.length,),
            (builder.INDEX, builder.ERROR_CODE)
        )

        builder.add_method(
            "delete_particle",
            (builder.INDEX,),
            (builder.ERROR_CODE)
        )

        builder.add_method(
            "get_state",
            (builder.INDEX,),
            (generic_unit_system.length, generic_unit_system.length, generic_unit_system.length, builder.ERROR_CODE),
            public_name = "get_position"
        )

        builder.add_method(
            "set_state",
            (builder.INDEX, generic_unit_system.length, generic_unit_system.length, generic_unit_system.length,),
            (builder.ERROR_CODE),
            public_name = "set_position"
        )

        builder.add_method(
            "find_nearest_neighbors",
            (),
            (builder.ERROR_CODE),
            public_name = "run"
        )

        builder.add_method(
            "get_close_neighbors",
            (builder.INDEX,),
            (builder.LINK('particles'), builder.LINK('particles'), builder.LINK('particles'), builder.ERROR_CODE),
        )

        builder.add_method(
            "get_nearest_neighbor",
            (builder.INDEX,),
            (builder.LINK('particles'), generic_unit_system.length, builder.ERROR_CODE),
        )

        builder.add_method(
            "get_number_of_particles",
            (),
            (builder.NO_UNIT, builder.ERROR_CODE),
        )
        
    def define_particle_sets(self, builder):
        builder.define_set('particles', 'index_of_the_particle')
        builder.set_new('particles', 'new_particle')
        builder.set_delete('particles', 'delete_particle')
        builder.add_setter('particles', 'set_position')
        builder.add_getter('particles', 'get_position')
        builder.add_getter('particles', 'get_close_neighbors', names=('neighbor0', 'neighbor1', 'neighbor2') )
