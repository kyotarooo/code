//=========================================================================================
// Element-based Parametric Dislocation Dynamics (PDD) Code
//-----------------------------------------------------------------------------------------
// Description:
// - This code was developed with the finite element method-like data structure.
// - The dislocation stress calculation is conducted with an analytical
// soluation of the stress field of a straight dislocation segment.
// - Positions within dislocation elements are interpolated with a linear
// function.
//-----------------------------------------------------------------------------------------
// Visualization:
// - ParaView is recommended for the visualization of the simulation results.
// - The output file format is VTK.
//-----------------------------------------------------------------------------------------
// Additional features:
// - Superposition principle can be used to deal with the free surface.
// - Grid data of stress can be used to implement initial stress field into the
// simulation.
// - A simplified dislocation core model to account for the influence of core
// energy will be supported soon.
//-----------------------------------------------------------------------------------------
// Developer: A. Takahashi (Tokyo University of Science)
//=========================================================================================
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "dd_defs.h"
#include "dd_direct_interaction.h"
#include "dd_fem_correction.h"
#include "dd_fem_mechanical_behavior.h"
#include "dd_fem_rk_integral.h"
#include "dd_fem_trapezoidal_integral.h"
#include "dd_fem_write.h"
#include "dd_plot.h"
#include "dd_read.h"
#include "dd_rearrange.h"
#include "dd_restart.h"
#include "dd_struct.h"
#include "dd_write.h"
#include "fem_read.h"
#include "fem_solve.h"
#include "fem_stiffness.h"
#include "fem_struct.h"
#include "dd_terminate.h"
#include "fem_write.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#ifndef _NOFEM_
#ifndef _PARDISO_
#include <mpi.h>
#endif
#endif

void DD_ClearOldVTK(char *directory) {
    static int already_cleared = 0;
    char cmd[512];
    char output_dir[256];

    if (already_cleared) return;
    already_cleared = 1;

    // 出力ディレクトリの存在を確認し、なければ作成
    sprintf(output_dir, "%s/out", directory);
    if (access(output_dir, F_OK) != 0) {
        sprintf(cmd, "mkdir -p %s", output_dir);
        system(cmd);
    }
    
    // .vtkと.vtuファイルを削除
    sprintf(cmd, "rm -f %s/out/*.vtk %s/out/*.vtu", directory, directory);
    int ret = system(cmd);
    if (ret != 0) {
        fprintf(stderr, "[Warning] Failed to clear old vtk/vtu files in %s/out\n", directory);
    } else {
        fprintf(stdout, "[LOG] Cleared old vtk/vtu files\n");
    }
}

int main(int argc, char **argv) {
  int iStep = 0;
  DD_t dd;
  FEM_t fem;

#ifndef _NOFEM_
#ifndef _PARDISO_
  MPI_Init(&argc, &argv);
#endif
#endif
#ifdef _OPENMP
  omp_set_num_threads(32);
#endif

  // If the command line info is not correct, the usage of this code is shown
  if (argc != 2) {
    fprintf(stderr, "Usage: %s <working_directory>\n", argv[0]);
    return 1;
  }

  DD_ClearOldVTK(argv[1]);

  // Read input data
  DD_ReadInput(&dd, argv[1]);
  FEM_ReadInput(&fem, argv[1]);

  // Write simulation volume & mesh
  DD_WriteSimulationVolume(&dd, argv[1]);
  FEM_WriteMesh(&fem, argv[1]);

  // Read all data for restart
  //DD_Restart(&iStep, &dd, argv[1], READ_RESTART_DATA);

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

  // Dislocation dynamics
  for (; iStep < dd.step.n; iStep++) {
    if (iStep % LOG_INTERVAL == 0) {
      DD_WriteLog(iStep, &dd);
    }

    // Write all data for restart
    //DD_Restart(&iStep, &dd, argv[1], WRITE_RESTART_DATA);

    // Update correction stress
    if (iStep % CORRECTION_FIELD_UPDATE_INTERVAL == 0) {
      DD_FEM_CorrectionField(&dd, &fem);
    }

    // Update dislocation position using Runge-Kutta method
    DD_FEM_RungeKutta4thIntegral(&dd, &fem);
    // Update dislocation position using Trapezoidal integral
    // DD_FEM_TrapezoidalIntegral(&dd, &fem);

    // Strain rate test, check direct interactions, and rearrange dislocation
    // discretization
    DD_FEM_MechanicalBehavior(&dd, &fem);
    DD_DirectInteractions(&dd);
    DD_RearrangeDislocations(&dd);

    // Write dislocations and mechanical behavior
    if (iStep % dd.output.interval == 0) {
      DD_WriteTime(&dd, argv[1]);
      DD_WriteDislocations(&dd, argv[1]);
      DD_WriteMechanicalBehavior(&dd, argv[1]);
      DD_FEM_WriteDislocationDensity(&dd, &fem, argv[1]);
      DD_PlotDislocations(&dd);
      dd.output.id += 1;
    }
    // terminate simulation
    if (DD_TerminateSimulation(iStep, &dd)) {
      printf("Completed DD simulation.\n");
      exit(0);
    }
  }

#ifndef _NOFEM_
#ifndef _PARDISO_
  MPI_Finalize();
#endif
#endif
  return 0;
}
