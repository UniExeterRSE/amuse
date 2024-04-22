#!/bin/bash        	# this line only there to enable syntax highlighting in this file

#########################################################
#  Enable/Disable compile-time options as needed        #
#  examples/galaxy_merger_star_formation_3d/Config.sh   #
#########################################################


#--------------------------------------- Mesh motion and regularization
AREPOFLAGS += -DREGULARIZE_MESH_CM_DRIFT                 # Mesh regularization; Move mesh generating point towards center of mass to make cells rounder.
AREPOFLAGS += -DREGULARIZE_MESH_CM_DRIFT_USE_SOUNDSPEED  # Limit mesh regularization speed by local sound speed
AREPOFLAGS += -DREGULARIZE_MESH_FACE_ANGLE               # Use maximum face angle as roundness criterion in mesh regularization


#--------------------------------------- Refinement and derefinement
AREPOFLAGS += -DREFINEMENT_SPLIT_CELLS                   # Refinement
AREPOFLAGS += -DREFINEMENT_MERGE_CELLS                   # Derefinement
AREPOFLAGS += -DREFINEMENT_VOLUME_LIMIT                  # Limit the volume of cells and the maximum volume difference between neighboring cels
AREPOFLAGS += -DNODEREFINE_BACKGROUND_GRID               # Do not de-refine low-res gas cells in zoom simulations


#--------------------------------------- Time integration options
AREPOFLAGS += -DTREE_BASED_TIMESTEPS                     # non-local timestep criterion (take 'signal speed' into account)


#--------------------------------------- Gravity treatment
AREPOFLAGS += -DSELFGRAVITY                              # gravitational intraction between simulation particles/cells 	 
AREPOFLAGS += -DHIERARCHICAL_GRAVITY                     # use hierarchical splitting of the time integration of the gravity
AREPOFLAGS += -DCELL_CENTER_GRAVITY                      # uses geometric centers to calculate gravity of cells, only possible with HIERARCHICAL_GRAVITY
AREPOFLAGS += -DALLOW_DIRECT_SUMMATION                   # Performed direct summation instead of tree-based gravity if number of active particles < DIRECT_SUMMATION_THRESHOLD (= 3000 unless specified differently here)
AREPOFLAGS += -DDIRECT_SUMMATION_THRESHOLD=500           # Overrides maximum number of active particles for which direct summation is performed instead of tree based calculation
AREPOFLAGS += -DGRAVITY_NOT_PERIODIC                     # gravity is not treated periodically


#--------------------------------------- Gravity softening
AREPOFLAGS += -DNSOFTTYPES=2                             # Number of different softening values to which particle types can be mapped.
AREPOFLAGS += -DMULTIPLE_NODE_SOFTENING                  # If a tree node is to be used which is softened, this is done with the softenings of its different mass components
AREPOFLAGS += -DINDIVIDUAL_GRAVITY_SOFTENING=32          # bitmask with particle types where the softenig type should be chosen with that of parttype 1 as a reference type
AREPOFLAGS += -DADAPTIVE_HYDRO_SOFTENING                 # Adaptive softening of gas cells depending on their size


#--------------------------------------- Single/Double Precision
AREPOFLAGS += -DDOUBLEPRECISION=1                        # Mode of double precision: not defined: single; 1: full double precision 2: mixed, 3: mixed, fewer single precisions; unless short of memory, use 1.
AREPOFLAGS += -DNGB_TREE_DOUBLEPRECISION                 # if this is enabled, double precision is used for the neighbor node extension


#-------------------------------------------- Things for special behaviour
AREPOFLAGS += -DPROCESS_TIMES_OF_OUTPUTLIST              # goes through times of output list prior to starting the simulaiton to ensure that outputs are written as close to the desired time as possible (as opposed to at next possible time if this flag is not active)
AREPOFLAGS += -DOVERRIDE_PEANOGRID_WARNING               # don't stop if peanogrid is not fine enough

#--------------------------------------- Output/Input options
AREPOFLAGS += -DHAVE_HDF5                                # needed when HDF5 I/O support is desired (recommended)


#--------------------------------------- Testing and Debugging options
#AREPOFLAGS += -DDEBUG                                    # enables core-dumps


#--------------------------------------- non-standard phyiscs
AREPOFLAGS += -DENFORCE_JEANS_STABILITY_OF_CELLS         # this imposes an adaptive floor for the temperature
AREPOFLAGS += -DCOOLING                                  # Simple primordial cooling
AREPOFLAGS += -DUSE_SFR                                  # Star formation model, turning dense gas into collisionless partices
