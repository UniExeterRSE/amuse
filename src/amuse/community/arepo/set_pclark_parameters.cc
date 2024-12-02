#include "set_pclark_parameters.h"


void set_pclark_parameters(){

  // INITIAL CONDITIONS                                                
  strcpy(All.InitCondFile, "./cDobbsSetup");
  All.ICFormat = 3;  // 3 for HDF5 Format

  // OUTPUT CONFIG
  strcpy(All.OutputDir, "./");
  All.OutputListOn = 0;
  strcpy(All.OutputListFilename, "outputs.txt");

  // SNAPSHOTS
  strcpy(All.SnapshotFileBase, "snapshot");
  All.SnapFormat = 3;
  All.NumFilesPerSnapshot = 1;
  All.NumFilesWrittenInParallel = 1;  // Set to number of cores, speeds up I/O

  All.TimeBetSnapshot = 0.1;  // ~9,000 years
  All.TimeOfFirstSnapshot = 0.;
  All.TimeBetStatistics = 0.1;

  // RUN CHARACTERISTICS
  All.TimeBegin = 0;
  All.TimeMax = 0.5;

  All.BoxSize = 620.;  // Always want to be slighly larger than IC boundaires
  All.PeriodicBoundariesOn = 0;

  // COSMOLOGICAL PARAMETERS
  All.ComovingIntegrationOn = 0;
  All.Omega0 = 0.;  // Ignored if ComovingIntegrationOn is 0
  All.OmegaLambda = 0.;  // ^
  All.OmegaBaryon = 0.;  // ^
  All.HubbleParam = 1.;  // ^

  // MODEL PARAMETERS
  All.CoolingOn = 0;
  All.StarformationOn = 0;

  // SOFTENING LENGTHS
  All.GasSoftFactor = 2.;

  All.MinimumComovingHydroSoftening = 0.0002;  // Always want to be smaller than Sink Accretion Radius
  All.AdaptiveHydroSofteningSpacing = 1.2;

  All.SofteningComoving[0] = 0.0002;  // ^
  All.SofteningComoving[1] = 0.;  //
  All.SofteningComoving[2] = 0.;  //
  All.SofteningComoving[3] = 0.;  //
  All.SofteningComoving[4] = 0.;  //
  All.SofteningComoving[5] = 0.0002;  // ^

  All.SofteningMaxPhys[0] = 0.0002;  // ^
  All.SofteningMaxPhys[1] = 0.;  //
  All.SofteningMaxPhys[2] = 0.;  //
  All.SofteningMaxPhys[3] = 0.;  //
  All.SofteningMaxPhys[4] = 0.;  //
  All.SofteningMaxPhys[5] = 0.0002;  // ^

  All.SofteningTypeOfPartType[0] = 0;
  All.SofteningTypeOfPartType[1] = 1;  
  All.SofteningTypeOfPartType[2] = 2;  
  All.SofteningTypeOfPartType[3] = 3;  
  All.SofteningTypeOfPartType[4] = 4;  
  All.SofteningTypeOfPartType[5] = 5;  

  // INITIAL VALUES
  All.DesNumNgb = 33;  // Number of particles to work out initial density
  All.MaxNumNgbDeviation = 0;

  All.InitGasTemp = 20.;

  // PESCRIBED MINIMUMS
  All.MinGasTemp = 2.73;
  All.MinimumDensityOnStartUp = 1e-10;
  All.LimitUBelowThisDensity = 0;
  All.LimitUBelowCertainDensityToThisValue = 0;
  All.MinEgySpec = 0;

  // UNIT SYSTEM
  All.UnitLength_in_cm = 1e17;
  All.UnitMass_in_g = 1.991e33;
  All.UnitVelocity_in_cm_per_s = 36447.268200;
  All.GravityConstantInternal = 0;

  // DOMAIN PARAMETERS
  All.MultipleDomains = 8;
  All.TopNodeFactor = 4;
  All.ActivePartFracForNewDomainDecomp = 0.005;

  // CELL REFINEMENT AND SHAPING 
  All.CellMaxAngleFactor = 2.;
  All.CellShapingSpeed = 0.5;

  All.ReferenceGasPartMass = 0.01;  // Should be TotMass / nParticles
  All.TargetGasMassFactor = 1.;

  All.MaxVolumeDiff = 8;
  All.MaxVolume = 1500.;
  All.MinVolume = 4.46e-11;  // Should be small enough for 16 cells to fit in sink accretion radius

  All.RefinementCriterion = 2;
  All.DerefinementCriterion = 2;

  // TIME INTEGRATION ACCURACY
  All.TypeOfTimestepCriterion = 0;
  All.ErrTolIntAccuracy = 0.05;
  All.MaxSizeTimestep = 0.1;  // Should be the same size as the snapshot dump rate
  All.MinSizeTimestep = 0.;
  All.CourantFac = 0.4;

  // TREE AND FORCE ACCURACY
  All.TypeOfOpeningCriterion = 1;
  All.ErrTolTheta = 0.4;
  All.ErrTolForceAcc = 0.005;

  // CPU Stuff
  All.TimeLimitCPU = 259200;  // Seconds
  All.ResubmitOn = 0;
  strcpy(All.ResubmitCommand, "xyz");
  All.MaxMemSize = 3000;  // MB
  All.CpuTimeBetRestartFile = 3600;  // Seconds

}