#ifdef __cplusplus
extern "C" {
#define ___cplusplus
#undef __cplusplus
#endif

#include "src/main/allvars.h"
#include "src/main/proto.h"

typedef struct {
    double mass;                                        /// mass
    double x, y, z;                                     /// position
    double vx, vy, vz;                                  /// velocity
} dynamics_state;

typedef struct {
    double mass;                                        /// mass
    double x, y, z;                                     /// position
    double vx, vy, vz;                                  /// velocity
    double u;                                           /// entropy
} gas_state;

// query_struct.h
void query_density(const int npoints, const double x[], const double y[], const double z[], double density[]);

#ifdef ___cplusplus
}
#define __cplusplus
#endif
