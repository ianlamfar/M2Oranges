//
// Created by Elijah on 03/03/2018.
//

#include "../../util/clUtils/clUtils.h"
#include "../../tests/run_tests/run_tests.h"
#include "../../util/simUtils/simUtils.h"
#include <CL/cl.h>
#include <malloc.h>
#include <string.h>

#ifndef DEMORANGES_SIMRUNNER_H
#define DEMORANGES_SIMRUNNER_H

int
runSim(particle *hparticles, cl_ulong NUMPART, cl_kernel iterate_particle, cl_float particle_diameter, cl_bool periodic,
       cl_float domain_length, char prefix[], char log_dir[], float sim_length, float timestep, bool VERBOSE,
       bool LOG_DATA, bool log_vel, float log_step, cl_device_id device, cl_context context, int coupled, int analytic,
       int fixed_Re, cl_int model, int fixed_tau, cl_float tau_scale);

#endif //DEMORANGES_SIMRUNNER_H
