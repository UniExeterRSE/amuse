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

#ifdef ___cplusplus
}
#define __cplusplus
#endif
