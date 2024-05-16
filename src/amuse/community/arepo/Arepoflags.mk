#!/bin/bash            # this line only there to enable syntax highlighting in this file

## Modified from examples/Noh_3d/Config.sh
## config file for 3d Noh problem, modified with galaxy merger flags

#--------------------------------------- Basic operation mode of code
AREPOFLAGS += -DREFLECTIVE_X=2                           # in-/outflow boundary conditions in x direction
AREPOFLAGS += -DREFLECTIVE_Y=2                           # in-/outflow boundary conditions in y direction
AREPOFLAGS += -DREFLECTIVE_Z=2                           # in-/outflow boundary conditions in z direction

#--------------------------------------- Refinement and derefinement
AREPOFLAGS += -DREFINEMENT_SPLIT_CELLS                   # Refinement
AREPOFLAGS += -DREFINEMENT_MERGE_CELLS                   # Derefinement
AREPOFLAGS += -DREFINEMENT_VOLUME_LIMIT                  # Limit the volume of cells and the maximum volume difference between neighboring cels
AREPOFLAGS += -DNODEREFINE_BACKGROUND_GRID               # Do not de-refine low-res gas cells in zoom simulations

#--------------------------------------- Mesh motion and regularization
AREPOFLAGS += -DREGULARIZE_MESH_CM_DRIFT                 # Mesh regularization; Move mesh generating point towards center of mass to make cells rounder.
AREPOFLAGS += -DREGULARIZE_MESH_CM_DRIFT_USE_SOUNDSPEED  # Limit mesh regularization speed by local sound speed
AREPOFLAGS += -DREGULARIZE_MESH_FACE_ANGLE               # Use maximum face angle as roundness criterion in mesh regularization

#--------------------------------------- Time integration options
AREPOFLAGS += -DTREE_BASED_TIMESTEPS                     # non-local timestep criterion (take 'signal speed' into account)

#---------------------------------------- Single/Double Precision
AREPOFLAGS += -DDOUBLEPRECISION=1                        # Mode of double precision: not defined: single; 1: full double precision 2: mixed, 3: mixed, fewer single precisions; unless short of memory, use 1.
AREPOFLAGS += -DINPUT_IN_DOUBLEPRECISION                 # initial conditions are in double precision
AREPOFLAGS += -DOUTPUT_CENTER_OF_MASS                    # output centers of cells

#--------------------------------------- Gravity softening
AREPOFLAGS += -DNSOFTTYPES=2                             # Number of different softening values to which particle types can be mapped.
AREPOFLAGS += -DMULTIPLE_NODE_SOFTENING                  # If a tree node is to be used which is softened, this is done with the softenings of its different mass components
AREPOFLAGS += -DINDIVIDUAL_GRAVITY_SOFTENING=32          # bitmask with particle types where the softenig type should be chosen with that of parttype 1 as a reference type
AREPOFLAGS += -DADAPTIVE_HYDRO_SOFTENING                 # Adaptive softening of gas cells depending on their size

#--------------------------------------- Output/Input options
AREPOFLAGS += -DHAVE_HDF5                                # needed when HDF5 I/O support is desired; should this be standard?

#--------------------------------------- Testing and Debugging options
AREPOFLAGS += -DDEBUG                                    # enables core-dumps, should this be standard?
