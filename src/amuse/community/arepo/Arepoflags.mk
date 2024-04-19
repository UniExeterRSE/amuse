## examples/Noh_3d/Config.sh
## config file for 3d Noh probelm

#--------------------------------------- Basic operation mode of code
AREPOFLAGS += -DREFLECTIVE_X=2                           # in-/outflow boundary conditions in x direction
AREPOFLAGS += -DREFLECTIVE_Y=2                           # in-/outflow boundary conditions in y direction
AREPOFLAGS += -DREFLECTIVE_Z=2                           # in-/outflow boundary conditions in z direction

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

#--------------------------------------- Output/Input options
AREPOFLAGS += -DHAVE_HDF5                                # needed when HDF5 I/O support is desired; should this be standard?

#--------------------------------------- Testing and Debugging options
#AREPOFLAGS += -DDEBUG                                    # enables core-dumps, should this be standard?