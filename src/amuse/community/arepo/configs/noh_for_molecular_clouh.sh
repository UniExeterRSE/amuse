#!/bin/bash            # this line only there to enable syntax highlighting in this file

## Modified from examples/Noh_3d/Config.sh
## config file for 3d Noh problem, modified with galaxy merger flags

#--------------------------------------- Basic operation mode of code
REFLECTIVE_X=2                           # in-/outflow boundary conditions in x direction
REFLECTIVE_Y=2                           # in-/outflow boundary conditions in y direction
REFLECTIVE_Z=2                           # in-/outflow boundary conditions in z direction

#--------------------------------------- Refinement and derefinement
REFINEMENT_SPLIT_CELLS                   # Refinement
REFINEMENT_MERGE_CELLS                   # Derefinement
REFINEMENT_VOLUME_LIMIT                  # Limit the volume of cells and the maximum volume difference between neighboring cels
NODEREFINE_BACKGROUND_GRID               # Do not de-refine low-res gas cells in zoom simulations

#--------------------------------------- Mesh motion and regularization
REGULARIZE_MESH_CM_DRIFT                 # Mesh regularization; Move mesh generating point towards center of mass to make cells rounder.
REGULARIZE_MESH_CM_DRIFT_USE_SOUNDSPEED  # Limit mesh regularization speed by local sound speed
REGULARIZE_MESH_FACE_ANGLE               # Use maximum face angle as roundness criterion in mesh regularization

#--------------------------------------- Time integration options
TREE_BASED_TIMESTEPS                     # non-local timestep criterion (take 'signal speed' into account)

#---------------------------------------- Single/Double Precision
DOUBLEPRECISION=1                        # Mode of double precision: not defined: single; 1: full double precision 2: mixed, 3: mixed, fewer single precisions; unless short of memory, use 1.
INPUT_IN_DOUBLEPRECISION                 # initial conditions are in double precision
OUTPUT_CENTER_OF_MASS                    # output centers of cells

#--------------------------------------- Gravity softening
NSOFTTYPES=2                             # Number of different softening values to which particle types can be mapped.
MULTIPLE_NODE_SOFTENING                  # If a tree node is to be used which is softened, this is done with the softenings of its different mass components
INDIVIDUAL_GRAVITY_SOFTENING=32          # bitmask with particle types where the softenig type should be chosen with that of parttype 1 as a reference type
ADAPTIVE_HYDRO_SOFTENING                 # Adaptive softening of gas cells depending on their size

#--------------------------------------- Output/Input options
HAVE_HDF5                                # needed when HDF5 I/O support is desired; should this be standard?

#--------------------------------------- Testing and Debugging options
DEBUG                                    # enables core-dumps, should this be standard?
