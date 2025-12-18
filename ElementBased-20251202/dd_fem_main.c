//=========================================================================================
// Element-based Parametric Dislocation Dynamics (PDD) Code
//-----------------------------------------------------------------------------------------
// Description:
// - This code was developed with the finite element method-like data structure.
// - The dislocation stress calculation is conducted with an analytical
// solution of the stress field of a straight dislocation segment.
// - Positions within dislocation elements are interpolated with a linear
// function.
//-----------------------------------------------------------------------------------------
// Visualization:
// - ParaView is recommended for the visualization of the simulation results.
// - The output file format is VTK.
//-----------------------------------------------------------------------------------------
// Additional features:
// - Superposition principle can be used to deal with the free surfaces.
// - Grid data of stress can be used to implement initial stress field into the
// simulation.
// - A simplified dislocation core model to account for the influence of core
// energy will be supported soon.
//-----------------------------------------------------------------------------------------
// Developer: A. Takahashi (Tokyo University of Science)
//=========================================================================================
#include <stdio.h>
#include <stdlib.h>

#include "dd_defs.h"
#include "dd_direct_interaction.h"
#include "dd_fem_correction.h"
#include "dd_fem_mechanical_behavior.h"
#include "dd_fem_rk_integral.h"
#include "dd_fem_trapezoidal_integral.h"
#include "dd_fem_write.h"
#include "dd_file.h"
#include "dd_plot.h"
#include "dd_read.h"
#include "dd_rearrange.h"
#include "dd_restart.h"
#include "dd_struct.h"
#include "dd_terminate.h"
#include "dd_write.h"
#include "fem_read.h"
#include "fem_solve.h"
#include "fem_stiffness.h"
#include "fem_struct.h"
#include "fem_write.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#ifndef _NO_FEM_
#ifndef _PARDISO_
#include <mpi.h>
#endif
#endif

int main(int argc, char** argv) {
  int iStep = 0;
  DD_t dd;
  FEM_t fem;
  void (*time_integral)(DD_t* dd, FEM_t* fem);

#ifndef _NO_FEM_
#ifndef _PARDISO_
  MPI_Init(&argc, &argv);
#endif
#endif
#ifdef _OPENMP
  omp_set_num_threads(8);
#endif

  // If the command line info is not correct, the usage of this code is shown
  if (argc != 2) {
    fprintf(stderr, "Usage: %s <working_directory>\n", argv[0]);
    return 1;
  }

  // Read input data
  DD_ReadInput(&dd, argv[1]);
  FEM_ReadInput(&fem, argv[1]);
  char *output_dir = "out";

  // Open log files
  DD_OpenLogFiles(&dd, argv[1]);

  // Set up time integral method
  switch (dd.step.time_integral_method) {
    case TIME_INTEGRAL_RK1:
      time_integral = DD_FEM_RungeKutta1stIntegral;
      break;
    case TIME_INTEGRAL_RK2:
      time_integral = DD_FEM_RungeKutta2ndIntegral;
      break;
    case TIME_INTEGRAL_RK3:
      time_integral = DD_FEM_RungeKutta3rdIntegral;
      break;
    case TIME_INTEGRAL_RK4:
      time_integral = DD_FEM_RungeKutta4thIntegral;
      break;
    default:
      time_integral = DD_FEM_TrapezoidalIntegral;
      break;
  }

  // Write simulation volume & mesh
  DD_WriteSimulationVolume(&dd, argv[1]);
  FEM_WriteMesh(&fem, argv[1]);

  // Read all data for restart
  DD_Restart(&iStep, &dd, argv[1], READ_RESTART_DATA);

  // Elastic constants
  fem.material.e = dd.material.e;
  fem.material.g = dd.material.g;
  fem.material.v = dd.material.v;

  // For correction stress
  FEM_StiffnessMatrix(&fem);
  FEM_Factorization(&fem);
  DD_FEM_SurfaceTractionGridStress(&dd, &fem);

  // Rearrange dislocation elements
  DD_RearrangeDislocations(&dd);

  char debug_log_path[256];
  sprintf(debug_log_path, "%s/debug.log", output_dir);
  
  FILE *fp_debug = fopen(debug_log_path, "w");
  if (fp_debug == NULL) {
      fprintf(stderr, "Cannot open %s\n", debug_log_path);
      return 1;
  }
  fprintf(fp_debug, "Step, Time, dt\n");

  
  // Dislocation dynamics
  for (; iStep < dd.step.n; iStep++) {

    // ログファイルへの書き込み (毎ステップ)
    fprintf(fp_debug, "%d, %e, %e\n", iStep, dd.step.t, dd.step.dt);
    fflush(fp_debug); // 強制書き込み

    // 標準ログ出力 (intervalごと)
    if (iStep % LOG_INTERVAL == 0) {
      DD_WriteLog(iStep, &dd);
    }

    // Write all data for restart
    DD_Restart(&iStep, &dd, argv[1], WRITE_RESTART_DATA);

    // Update correction stress
    if (iStep % CORRECTION_FIELD_UPDATE_INTERVAL == 0) {
      DD_FEM_CorrectionField(&dd, &fem);
    }

    // Update dislocation positions (時間積分)
    time_integral(&dd, &fem);

    // Strain rate test, check direct interactions, and rearrange dislocation
    DD_FEM_MechanicalBehavior(&dd, &fem);
    DD_DirectInteractions(&dd);
    DD_RearrangeDislocations(&dd);

    // Write dislocations and mechanical behavior (結果出力)
    if (iStep % dd.output.interval == 0) {
      DD_WriteTime(&dd, output_dir);
      DD_WriteDislocations(&dd, output_dir);
      DD_WriteMechanicalBehavior(&dd, output_dir);
      DD_FEM_WriteDislocationDensity(&dd, &fem, output_dir);
      // DD_PlotDislocations(&dd);
      dd.output.id += 1;
    }

    // Terminate simulation, if prescribed condition is satisfied
    DD_TerminateSimulation(&dd);
  }

  // 3. ループが終わったらファイルを閉じる
  fclose(fp_debug);
  // --- 【ここまで修正】 ---

#ifndef _NO_FEM_
#ifndef _PARDISO_
  MPI_Finalize();
#endif
#endif
  return 0;
}