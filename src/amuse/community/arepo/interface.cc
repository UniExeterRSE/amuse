#include <cstdio>
#include <cstring>
#include <map>

#ifndef NOMPI
#include <mpi.h>
#endif

#include "worker_code.h"
#include "interface.h"


// general interface functions:

using namespace std;


// Global ID_RLOOKUP will be initalized by create_ID_reverse_lookup in
// initialize_code()
map<MyIDType, size_t> ID_RLOOKUP;
static void create_ID_reverse_lookup();

long long dm_particles_in_buffer = 0;
long long gas_particles_in_buffer = 0;
long long particle_id_counter = 0;  // dm and sph use same counter
map<long long, dynamics_state> dm_states;
map<long long, gas_state> gas_states;

void set_default_parameters(){
  // Relevant files
  strcpy(All.InitCondFile, "./ICs");
  strcpy(All.OutputDir,   "./output");
  strcpy(All.SnapshotFileBase, "snap");
  strcpy(All.OutputListFilename, "./output_list.txt");

  // File formats
  All.ICFormat = 1;
  All.SnapFormat = 1;

  // CPU-time LimitUBelowThisDensity
  All.TimeLimitCPU = 93000;
  All.CpuTimeBetRestartFile = 12000;
  All.ResubmitOn = 0;
  strcpy(All.ResubmitCommand, "my-scriptfile");

  // Memory allocation
  All.MaxMemSize = 2500;

  // Characteristics of run
  All.TimeBegin = 0.0;
  All.TimeMax = 2.7;
  All.TimeStep = 0.00314159;

  // Basic code options that set simulation type
  All.ComovingIntegrationOn = 0;
  All.PeriodicBoundariesOn = 0;
  All.CoolingOn = 0;
  All.StarformationOn = 0;

  // Cosmological parameters
  All.Omega0 = 0.0;
  All.OmegaLambda = 0.0;
  All.OmegaBaryon = 0.0;
  All.HubbleParam = 1.0;
  All.BoxSize = 100000.0;

  // Output frequency and output parameters
  All.OutputListOn = 1;
  All.TimeBetSnapshot = 0.1;
  All.TimeOfFirstSnapshot = 0.0;
  All.TimeBetStatistics = 0.01;
  All.NumFilesPerSnapshot = 1;
  All.NumFilesWrittenInParallel = 1;

  // Integration timing accuracy
  All.TypeOfTimestepCriterion = 0;
  All.ErrTolIntAccuracy = 0.012;
  All.CourantFac = 0.3;
  All.MaxSizeTimestep = 0.05;
  All.MinSizeTimestep = 2.0e-9;

  // Treatment of empty space and temp limits
  All.InitGasTemp = 244.8095;
  All.MinGasTemp = 5.0;
  All.MinimumDensityOnStartUp = 1.0e-20;
  All.LimitUBelowThisDensity = 0.0;
  All.LimitUBelowCertainDensityToThisValue = 0.0;
  All.MinEgySpec = 0.0;

  // Tree algorithm, force accuracy, domain update frequency
  All.TypeOfOpeningCriterion = 1;
  All.ErrTolTheta = 0.7;
  All.ErrTolForceAcc = 0.0025;
  All.MultipleDomains = 8;
  All.TopNodeFactor = 2.5;
  All.ActivePartFracForNewDomainDecomp = 0.01;

  // Initial density estimates
  All.DesNumNgb = 64;
  All.MaxNumNgbDeviation = 4;

  // System of Units
  All.UnitLength_in_cm = 3.085678e21;
  All.UnitMass_in_g = 1.989e43;
  All.UnitVelocity_in_cm_per_s = 1e5;

  // Gravitational softening lengths
  All.SofteningComoving[0] = 1.0;
  All.SofteningComoving[1] = 1.0;
  All.SofteningMaxPhys[0] = 1.0;
  All.SofteningMaxPhys[1] = 1.0;
  All.GasSoftFactor = 2.5;


  All.SofteningTypeOfPartType[0] = 0;
  All.SofteningTypeOfPartType[1] = 1;
  All.SofteningTypeOfPartType[2] = 1;
  All.SofteningTypeOfPartType[3] = 1;
  All.SofteningTypeOfPartType[4] = 1;
  All.SofteningTypeOfPartType[5] = 1;
  #ifdef ADAPTIVE_HYDRO_SOFTENING
    All.MinimumComovingHydroSoftening = 1.0;
    All.AdaptiveHydroSofteningSpacing = 1.2;
  #endif

  // Mesh regularization options
  All.CellShapingSpeed = 0.5;
  #ifndef REGULARIZE_MESH_FACE_ANGLE  // Compiler error if flag defined
    All.CellShapingFactor = 1.0;
  #endif

  // parameters that are fixed for AMUSE:
  All.TreeAllocFactor = 0.8; // Memory allocation parameter
  All.ResubmitOn = 0;              // Keep this turned off!
  All.OutputListOn = 0;            // Keep this turned off
  All.GravityConstantInternal = 0; // Keep this turned off
}

void set_noh_3d_parameters(){
  // Relevant files
  strcpy(All.InitCondFile, "./IC");
  strcpy(All.OutputDir,   "./output");
  strcpy(All.SnapshotFileBase, "snap");
  strcpy(All.OutputListFilename, "./output_list.txt");

  All.ICFormat = 3;

  All.SnapFormat = 3;
  All.NumFilesPerSnapshot = 1;
  All.NumFilesWrittenInParallel = 1;

  All.ResubmitOn = 0;
  strcpy(All.ResubmitCommand, "my-scriptfile");
  All.OutputListOn = 0;

  All.CoolingOn = 0;
  All.StarformationOn = 0;

  All.Omega0 = 0.0;
  All.OmegaBaryon = 0.0;
  All.OmegaLambda = 0.0;
  All.HubbleParam = 1.0;

  All.BoxSize = 6.0;
  All.PeriodicBoundariesOn = 1;
  All.ComovingIntegrationOn = 0;

  All.MaxMemSize = 2500;

  All.TimeOfFirstSnapshot = 10.0;
  All.CpuTimeBetRestartFile = 9000;
  All.TimeLimitCPU = 90000;

  All.TimeBetStatistics = 0.005;
  All.TimeBegin = 0.0;
  All.TimeMax = 2.0;
  All.TimeBetSnapshot = 0.5;

  All.UnitVelocity_in_cm_per_s = 1.0;
  All.UnitLength_in_cm = 1.0;
  All.UnitMass_in_g = 1.0;
  All.GravityConstantInternal = 0.0;

  All.ErrTolIntAccuracy = 0.1;
  All.ErrTolTheta = 0.1;
  All.ErrTolForceAcc = 0.1;

  All.MaxSizeTimestep = 0.5;
  All.MinSizeTimestep = 1e-5;
  All.CourantFac = 0.3;

  All.LimitUBelowThisDensity = 0.0;
  All.LimitUBelowCertainDensityToThisValue = 0.0;
  All.DesNumNgb = 64;
  All.MaxNumNgbDeviation = 2;

  All.MultipleDomains = 2;
  All.TopNodeFactor = 4;
  All.ActivePartFracForNewDomainDecomp = 0.1;

  All.TypeOfTimestepCriterion = 0;
  All.TypeOfOpeningCriterion = 1;
  All.GasSoftFactor = 0.01;

  All.SofteningComoving[0] = 0.1;
  All.SofteningComoving[1] = 0.1;
  All.SofteningComoving[2] = 0.1;
  All.SofteningComoving[3] = 0.1;
  All.SofteningComoving[4] = 0.1;
  All.SofteningComoving[5] = 0.1;

  All.SofteningMaxPhys[0] =  0.1;
  All.SofteningMaxPhys[1] =  0.1;
  All.SofteningMaxPhys[2] =  0.1;
  All.SofteningMaxPhys[3] =  0.1;
  All.SofteningMaxPhys[4] =  0.1;
  All.SofteningMaxPhys[5] =  0.1;

  All.SofteningTypeOfPartType[0] = 0;
  All.SofteningTypeOfPartType[1] = 1;
  All.SofteningTypeOfPartType[2] = 1;
  All.SofteningTypeOfPartType[3] = 1;
  All.SofteningTypeOfPartType[4] = 1;
  All.SofteningTypeOfPartType[5] = 1;

  All.InitGasTemp = 0.0;
  All.MinGasTemp = 0.0;
  All.MinEgySpec = 0.0;
  All.MinimumDensityOnStartUp = 0.0;

  All.CellShapingSpeed = 0.5;
  #ifdef REGULARIZE_MESH_FACE_ANGLE
    All.CellMaxAngleFactor = 2.25;
  #endif

}

int initialize_code(){

#ifndef NOMPI
    MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
    MPI_Comm_size(MPI_COMM_WORLD, &NTask);
#else
    ThisTask = 0;
    NTask = 1;
#endif

  /* output a welcome message */
  hello();

  /* initialize CPU-time/Wallclock-time measurement */
  init_cpu_log();

  determine_compute_nodes();
  // Needed to check available memory
  mpi_report_committable_memory();

  // set_default_parameters();
  set_noh_3d_parameters();


  RestartFlag = 0;

  // May not need to do this (we want AMUSE to manage this)
  // MPI_Bcast(&All, sizeof(struct global_data_all_processes), MPI_BYTE, 0, MPI_COMM_WORLD);
  begrun1(); /* set-up run  */
  // this needs to be called by "commit_parameters" in AMUSE, transitioning EDIT->RUN

  char fname[MAXLEN_PATH];
  strcpy(fname, All.InitCondFile);

  /* now we can load the file */

// #ifdef READ_DM_AS_GAS
//       read_ic(fname, (RestartFlag == 14) ? 0x02 : LOAD_TYPES);
// #else  /* #ifdef READ_DM_AS_GAS */
//       read_ic(fname, (RestartFlag == 14) ? 0x01 : LOAD_TYPES);
// #endif /* #ifdef READ_DM_AS_GAS #else */



  return 0;
}

int run_sim() {
  /* This run command is for the Arepo simulation */
  run();
  return 0;
}

int cleanup_code(){
  printf("Code run for %f seconds!\n", timediff(StartOfRun, second()));
  printf("endrun called, calling MPI_Finalize()\nbye!\n\n");
  fflush(stdout);

#ifdef HAVE_HDF5
  /*The hdf5 library will sometimes register an atexit() handler that calls its
   * error handler. In AREPO this is set to my_hdf_error_handler, which calls
   * MPI_Abort. Calling MPI_Abort after MPI_Finalize is not allowed.
   * Hence unset the HDF error handler here
   */
  H5Eset_auto(NULL, NULL);
#endif /* #ifdef HAVE_HDF5 */

  MPI_Finalize();
  exit(0);
  return 0;
}

static void create_ID_reverse_lookup() {
  map<MyIDType, size_t> id_rlookup_local;
  for (size_t i = 0; i < NumPart; i++) {
    MyIDType id = P[i].ID;
    id_rlookup_local[id] = i;
  }
  ID_RLOOKUP = id_rlookup_local;
}

static int find_particle_with_ID(int particle_id) {
  for (int j = 0; j < 2; ++j) {
    // If there's a problem, rebuild ID_RLOOKUP and retry once
    // These cautious checks may add unnecessary overhead

    auto it = ID_RLOOKUP.find(particle_id);

    if (it == ID_RLOOKUP.end()) {
      // particle_id wasn't in the map - rebuild ID_RLOOKUP and try again
      printf("AMUSE: Rebuilding particle_ID lookup (ID not found).\n");
      create_ID_reverse_lookup();
      continue;
    }

    size_t particle_pos = (*it).second;

    if (P[particle_pos].ID != particle_id) {
      // particle_id had the wrong value - rebuild ID_RLOOKUP and try again
      printf("AMUSE: Rebuilding particle ID lookup table (ID index changed).\n");
      create_ID_reverse_lookup();
      continue;
    }

    return particle_pos;
  }
  return -1;
}

int get_mass(int index_of_the_particle, double * mass){
  int p = find_particle_with_ID(index_of_the_particle);
  if (p >= 0) {
    *mass = P[p].Mass;
    return 0;
  }
  return -3;
}

int commit_particles(){
  // Refer to gadget2 interface.cc, also read_ic.c
  double t0, t1;
  int i, j;
  All.TotNumPart = gas_particles_in_buffer + dm_particles_in_buffer;
  All.TotNumGas = gas_particles_in_buffer;
  // All.MaxPart = All.TreeAllocFactor * (All.TotNumPart / NTask); // TODO: Check TreeAllocFactor is right
  // All.MaxPartSph = 0; // TODO:

  printf("Initialising %lli particles, of which %lli gas\n", All.TotNumPart, All.TotNumGas);

  double a;
  if (All.ComovingIntegrationOn) {
    a = All.Time;
  } else {
    a = 1.0;
  }
  double a_inv = 1.0 / a;

  NumPart = dm_states.size() + gas_states.size();
  NumGas = gas_states.size();

  // TODO: Should be max of NumParts if working across multiple processors
  All.MaxPart    = NumPart / (1.0 - 2 * ALLOC_TOLERANCE);
  All.MaxPartSph = NumGas / (1.0 - 2 * ALLOC_TOLERANCE);
  allocate_memory();
  dump_memory_table();

  // Initialize gas particles; From gadget2 interface.cc
  i = 0;
  for (map<long long, gas_state>::iterator state_iter = gas_states.begin();
          state_iter != gas_states.end(); state_iter++, i++){
      P[i].ID = (*state_iter).first;
      P[i].Mass = (*state_iter).second.mass;
      P[i].Pos[0] = (*state_iter).second.x * a_inv;
      P[i].Pos[1] = (*state_iter).second.y * a_inv;
      P[i].Pos[2] = (*state_iter).second.z * a_inv;
      P[i].Vel[0] = (*state_iter).second.vx * a;
      P[i].Vel[1] = (*state_iter).second.vy * a;
      P[i].Vel[2] = (*state_iter).second.vz * a;
      P[i].Type = 0; // SPH particles (dark matter particles have type 1)
      SphP[i].Utherm = (*state_iter).second.u;
      SphP[i].Density = -1;
      SphP[i].Hsml = 0;
#ifdef MORRIS97VISC
      SphP[i].Alpha = (*state_iter).second.alpha;
      SphP[i].DAlphaDt = (*state_iter).second.dalphadt;
#endif
  }
  gas_states.clear();

  // initialize dark matter particles
  for (map<long long, dynamics_state>::iterator state_iter = dm_states.begin();
          state_iter != dm_states.end(); state_iter++, i++){
      P[i].ID = (*state_iter).first;
      P[i].Mass = (*state_iter).second.mass;
      P[i].Pos[0] = (*state_iter).second.x * a_inv;
      P[i].Pos[1] = (*state_iter).second.y * a_inv;
      P[i].Pos[2] = (*state_iter).second.z * a_inv;
      P[i].Vel[0] = (*state_iter).second.vx * a;
      P[i].Vel[1] = (*state_iter).second.vy * a;
      P[i].Vel[2] = (*state_iter).second.vz * a;
      P[i].Type = 1; // 1=halo type (dm)
  }
  dm_states.clear();
  All.TimeBegin += All.Ti_Current * All.Timebase_interval; // TODO: Why?
  All.Ti_Current = 0;
  All.Time = All.TimeBegin;
  set_softenings(); // TODO: Why?
  for (i = 0; i < NumPart; i++){ // Start-up initialization
      P[i].GravAccel[0] = 0;
      P[i].GravAccel[1] = 0;
      P[i].GravAccel[2] = 0;
      // TODO: PMGRID?
      // TODO: Ti_endstep? Ti_begstep?
      P[i].OldAcc = 0;
      // P[i].GravCost = 1;
      // TODO: Potential?
  }

  /* init returns a status code, where a value of >=0 means that endrun() should be called. */
  int status = init();

  if(status >= 0)
    {
      if(status > 0)
        printf("init() returned with %d\n", status);

      cleanup_code();
    }

  begrun2();

  return 0; // Here
}

int get_time(double * time){
  // Return error code if calling from non-zero task
  if (ThisTask) {return 0;}
  *time = All.Time;
  return 0;
}

int set_mass(int index_of_the_particle, double mass){
  return 0;
}

int get_index_of_first_particle(int * index_of_the_particle){
  return 0;
}

int get_total_radius(double * radius){
  return 0;
}

int new_dm_particle(int * index_of_the_particle, double mass, double x,
  double y, double z, double vx, double vy, double vz){
  particle_id_counter++;
  if (ThisTask == 0)
      *index_of_the_particle = particle_id_counter; // id undefined?

  // Divide the particles equally over all Tasks, arepo will redistribute them later.
  if (ThisTask == (dm_particles_in_buffer % NTask)){
      dynamics_state state;
      state.mass = mass;
      state.x = x;
      state.y = y;
      state.z = z;
      state.vx = vx;
      state.vy = vy;
      state.vz = vz;
      dm_states.insert(std::pair<long long, dynamics_state>(particle_id_counter, state));
  }
  dm_particles_in_buffer++;
  return 0;

}

int new_gas_particle(int * index_of_the_particle, double mass, double x,
  double y, double z, double vx, double vy, double vz, double u){
  particle_id_counter++;
  if (ThisTask == 0)
      *index_of_the_particle = particle_id_counter; // id undefined?

  // Divide the particles equally over all Tasks, arepo will redistribute them later.
  if (ThisTask == (gas_particles_in_buffer % NTask)){
      gas_state state;  // TODO: Check these properties match up
      state.mass = mass;
      state.x = x;
      state.y = y;
      state.z = z;
      state.vx = vx;
      state.vy = vy;
      state.vz = vz;
      state.u = u;
      gas_states.insert(std::pair<long long, gas_state>(particle_id_counter, state));
  }
  gas_particles_in_buffer++;  // TODO: initialised?
  return 0;

}

int get_total_mass(double * mass){
  return 0;
}

int evolve_model(double time){
  printf("AMUSE interface: setting TimeMax from %g to %g\n", All.TimeMax, time);
  All.TimeMax = time;
  //All.TimeStep = time - All.Time;
  run();
  return 0;
}

int set_eps2(double epsilon_squared){
  // This looks bizarre
  if (ThisTask) {return 0;}
  return -2;
}

int get_begin_time(double * time){
  if (ThisTask) {return 0;}
  *time = All.TimeBegin;
  return 0;
}

int get_eps2(double * epsilon_squared){
  if (ThisTask) {return 0;}
  return -2;
}

int get_index_of_next_particle(int index_of_the_particle,
  int * index_of_the_next_particle){
  return 0;
}

int delete_particle(int index_of_the_particle){
  return 0;
}

int get_potential(int index_of_the_particle, double * potential){
  return 0;
}

int synchronize_model(){
  return 0;
}

int set_state(int index_of_the_particle, double mass, double x, double y,
  double z, double vx, double vy, double vz, double radius){
  return 0;
}

int get_state(int index_of_the_particle, double * mass, double * x,
  double * y, double * z, double * vx, double * vy, double * vz,
  double * radius){
    // Arepo has Delaunay cell radii.
  return 0;
}

int set_state_gas(int index_of_the_particle, double mass, double x, double y,
  double z, double vx, double vy, double vz, double u){
  return -1;
}

int get_state_gas(int index_of_the_particle, double * mass, double * x,
  double * y, double * z, double * vx, double * vy, double * vz,
  double * u){
  int p_idx = find_particle_with_ID(index_of_the_particle);
  if (p_idx < 0) {
    printf("AREPO: Particle with ID %d not found in P", index_of_the_particle);
    return p_idx;
  }

  if (P[p_idx].Type > 0){
    printf("AREPO: Particle with index %d not gas", index_of_the_particle);
    return -2;
  }
  *mass = P[p_idx].Mass;
  *x = P[p_idx].Pos[0];
  *y = P[p_idx].Pos[1];
  *z = P[p_idx].Pos[2];
  *vx = P[p_idx].Vel[0];
  *vy = P[p_idx].Vel[1];
  *vz = P[p_idx].Vel[2];
  *u = SphP[p_idx].Utherm;
  return 0;
}

int get_minimum_time_step(double * minimum_time_step){
  if (ThisTask) {return 0;}
  *minimum_time_step = All.MinSizeTimestep;
  return 0;
}

int set_minimum_time_step(double minimum_time_step){
  if (ThisTask) {return 0;}
  All.MinSizeTimestep = minimum_time_step;
  return 0;
}

int get_time_step(double * time_step){
  if (ThisTask) {return 0;}
  *time_step = All.TimeStep;
  return 0;
}

int recommit_particles(){
  return 0;
}

int get_kinetic_energy(double * kinetic_energy){
  return 0;
}

int get_number_of_particles(int * number_of_particles){
  if (ThisTask) {return 0;}
  *number_of_particles = NumPart;
  return 0;
}

int get_number_of_gas_particles(int * number_of_particles){
  if (ThisTask) {return 0;}
  *number_of_particles = All.TotNumPart;
  return 0;
}

int set_acceleration(int index_of_the_particle, double ax, double ay,
  double az){
  return -2;
}

int get_center_of_mass_position(double * x, double * y, double * z){
  return 0;
}

int get_center_of_mass_velocity(double * vx, double * vy, double * vz){
  return 0;
}

int get_radius(int index_of_the_particle, double * radius){
  return -2;
}

int set_begin_time(double time){
  All.TimeBegin = time;
  return 0;
}

int set_radius(int index_of_the_particle, double radius){
  return -2;
}

int recommit_parameters(){
  return 0;
}

int get_potential_energy(double * potential_energy){
  return 0;
}

int get_velocity(int index_of_the_particle, double * vx, double * vy,
  double * vz){
  int p = find_particle_with_ID(index_of_the_particle);
  if (p >= 0) {
    *vx = P[p].Vel[0];
    *vy = P[p].Vel[1];
    *vz = P[p].Vel[2];
    return 0;
  }
  return -3;
}

int get_position(int index_of_the_particle, double * x, double * y,
  double * z){
    int p = find_particle_with_ID(index_of_the_particle);
    if (p >= 0) {
      *x = P[p].Pos[0];
      *y = P[p].Pos[1];
      *z = P[p].Pos[2];
      return 0;
    }
    return -3;
}

int set_position(int index_of_the_particle, double x, double y, double z){
  return 0;
}

int get_acceleration(int index_of_the_particle, double * ax, double * ay,
  double * az){
    int p = find_particle_with_ID(index_of_the_particle);
    if (p < 0) {
      return p;
    }
    *ax = P[p].GravAccel[0];
    *ay = P[p].GravAccel[1];
    *az = P[p].GravAccel[2];
    return 0;
}

int commit_parameters(){
  return 0;
}

int set_parameters(char * param_file){
  return 0;
}

int set_velocity(int index_of_the_particle, double vx, double vy,
  double vz){
  return 0;
}

int get_pressure(int index_of_the_particle, double * p){
  int p_idx = find_particle_with_ID(index_of_the_particle);
  if (p_idx < 0) {
    printf("AREPO: Particle with ID %d not found in P", index_of_the_particle);
    return p_idx;
  }

  if (P[p_idx].Type > 0){
    printf("AREPO: Particle with index %d not gas", index_of_the_particle);
    return -2;
  }

  *p = SphP[p_idx].Pressure;
  return 0;
}

int get_density(int index_of_the_particle, double * rho){
  int p_idx = find_particle_with_ID(index_of_the_particle);
  if (p_idx < 0) {
    printf("AREPO: Particle with ID %d not found in P", index_of_the_particle);
    return p_idx;
  }

  if (P[p_idx].Type > 0){
    printf("AREPO: Particle with index %d not gas", index_of_the_particle);
    return -2;
  }

  *rho = SphP[p_idx].Density;
  return 0;
}

int get_internal_energy(int index_of_the_particle, double * u){
  int p_idx = find_particle_with_ID(index_of_the_particle);
  if (p_idx < 0) {
    printf("AREPO: Particle with ID %d not found in P", index_of_the_particle);
    return p_idx;
  }

  if (P[p_idx].Type > 0){
    printf("AREPO: Particle with index %d not gas", index_of_the_particle);
    return -2;
  }

  *u = SphP[p_idx].Utherm;
  return 0;
}

int set_internal_energy(int index_of_the_particle, double u){
  int p_idx = find_particle_with_ID(index_of_the_particle);
  if (p_idx < 0) {
    printf("AREPO: Particle with ID %d not found in P", index_of_the_particle);
    return p_idx;
  }

  if (P[p_idx].Type > 0){
    printf("AREPO: Particle with index %d not gas", index_of_the_particle);
    return -2;
  }

  SphP[p_idx].Utherm = u;
  return 0;
}


int get_box_size(double *value)
{
    if (ThisTask) {return 0;}
    *value = All.BoxSize;
    return 0;
}

int set_box_size(double value)
{
    All.BoxSize = value;
    return 0;
}

int get_omega_zero(double *omega_zero){
    if (ThisTask) {return 0;}
    *omega_zero = All.Omega0;
    return 0;
}
int set_omega_zero(double omega_zero){
    All.Omega0 = omega_zero;
    return 0;
}
int get_omega_lambda(double *omega_lambda){
    if (ThisTask) {return 0;}
    *omega_lambda = All.OmegaLambda;
    return 0;
}
int set_omega_lambda(double omega_lambda){
    All.OmegaLambda = omega_lambda;
    return 0;
}
int get_omega_baryon(double *omega_baryon){
    if (ThisTask) {return 0;}
    *omega_baryon = All.OmegaBaryon;
    return 0;
}
int set_omega_baryon(double omega_baryon){
    All.OmegaBaryon = omega_baryon;
    return 0;
}
int get_hubble_param(double *hubble_param){
    if (ThisTask) {return 0;}
    *hubble_param = All.HubbleParam;
    return 0;
}
int set_hubble_param(double hubble_param){
    All.HubbleParam = hubble_param;
    return 0;
}
