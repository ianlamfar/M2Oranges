//
// Created by Ian on 17/2/2021.
//

#include "simRunner_tgv.h"
#include "../../util/particleUtils/particleUtils.h"
#include <pthread.h>
int nexttgv = 0;

int runSim_tgv(particle *hparticles, cl_ulong NUMPART, cl_kernel iterate_particle, cl_float particle_diameter, cl_bool periodic,
           cl_float domain_length, char prefix[], char log_dir[], float sim_length, float timestep, bool VERBOSE,
           bool LOG_DATA, bool log_vel, float log_step, cl_device_id device, cl_context context, int coupled, int analytic,
           int fixed_Re, int model, int fixed_tau, cl_float tau_scale, cl_float tgv_scale) {
    cl_int ret;  // Variable in which to store OpenCL return values.

    cl_mem gparticles;  // GPU array of particles.

    // Check that the logging directory exists if needed. Construct directory if it does not exist.
    if (!checkDirExists(log_dir)) {
        mkdir(log_dir);
//        fprintf(stderr, "Error: Directory (%s) does not exist or cannot be accessed.\n", log_dir);
//        return 1;
    }


    // Create command queue.
    cl_command_queue queue = getCommandQueue(context, device, VERBOSE);

    // Create gparticles buffer and copy particles into it.
    gparticles = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(particle) * NUMPART, NULL, &ret);

    ret = particlesToDevice(queue, gparticles, &hparticles, NUMPART);

    int int_periodic;  // Enumerated version of periodic boolean as booleans can't be passed to GPU by default.
    if (periodic) {
        int_periodic = 1;
    } else {
        int_periodic = 0;
    }

    // Set all kernel arguments. Note that if the cl_mem object changes they must be reset (e.g. gpp_cols).
    printf("[INIT] Setting kernel arguments\n");
    ret = clSetKernelArg(iterate_particle, 0, sizeof(cl_mem), &gparticles);
    ret = clSetKernelArg(iterate_particle, 1, sizeof(cl_float), &timestep);
    ret = clSetKernelArg(iterate_particle, 2, sizeof(cl_int), &int_periodic);
    ret = clSetKernelArg(iterate_particle, 3, sizeof(cl_float), &domain_length);
    ret = clSetKernelArg(iterate_particle, 4, sizeof(cl_int), &coupled);
    ret = clSetKernelArg(iterate_particle, 5, sizeof(cl_int), &analytic);
    ret = clSetKernelArg(iterate_particle, 6, sizeof(cl_int), &fixed_Re);
    ret = clSetKernelArg(iterate_particle, 7, sizeof(cl_int), &model);
    ret = clSetKernelArg(iterate_particle, 8, sizeof(cl_int), &fixed_tau);
    ret = clSetKernelArg(iterate_particle, 9, sizeof(cl_float), &tau_scale);
    ret = clSetKernelArg(iterate_particle, 10, sizeof(cl_float), &tgv_scale);

    // Do pre-simulation output and logging.
    printf("Running sim with %llu particles, timestep %f, and log step %f, length %f.\n", NUMPART, timestep,
           log_step, sim_length);

    if (!writeTime(prefix, log_dir, NUMPART, "Start")) {
        return 1;
    }
    if (!writeSetupData(prefix, log_dir, NUMPART, timestep, sim_length, domain_length, particle_diameter,
                        hparticles[0].density, hparticles[0].fluid_viscosity, get_tau(&(hparticles[0])))) {
        return 1;
    }
    if (LOG_DATA) {
        //printf("Logging at time: 0.000000\n");
        printf("LOG DATA = TRUE\n");
        if (!writeParticles(hparticles, 0, prefix, log_dir, NUMPART, log_vel)) {
            return 1;
        }
    }

    cl_float last_write = 0; // Variable to store when logging was last run.

    // Run simulation.
    size_t local;
    clGetDeviceInfo(device, CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(local), &local, NULL);
//    printf("Workgroup size = %i\n", local);

    for (cl_float time = timestep; time <= sim_length; time += timestep) {
        if (VERBOSE) {printf("Checking particle mass...\n");}
        bool cont = parallelMassChecktgv(hparticles, NUMPART);
        if (cont) {
            if (VERBOSE) printf("   Time = %f\n", time);

            if (ret != 0) {
                fprintf(stderr, "Error: Failed to read from collision count buffer. Error code %i\n", ret);
                return 1;
            }
            // Iterate particles.
            if (VERBOSE) printf("   Iterating particles\n");
            ret = clEnqueueNDRangeKernel(queue, iterate_particle, 1, NULL, &NUMPART, 0, 0, NULL, NULL);

            if (LOG_DATA && time - last_write >= 0.99*log_step) {
                if (VERBOSE) {printf("   Writing data...\n");}
                ret = particlesToHost(queue, gparticles, &hparticles, NUMPART);

                if (!writeParticles(hparticles, time, prefix, log_dir, NUMPART, log_vel)) {
                    return 1;
                }

                last_write = time;
            }
        } else{
            printf("ALL DROPLETS REACHED ZERO MASS, ITERATION TERMINATED (Time = %f)\n", last_write);
            break;
        }
    }
    // Release OpenCL memory and command queue to prevent stacking up of memory
    clReleaseMemObject(gparticles);
//    clReleaseCommandQueue(queue);


//    if (LOG_DATA) {
//        // printf("Logging at time: %f\n", sim_length);
//        if (!writeParticles(hparticles, sim_length, prefix, log_dir, NUMPART, log_vel)) {
//            return 1;
//        }
//    }
//
//    if (!writeTime(prefix, log_dir, NUMPART, "End")) {
//        return 1;
//    }
//    return 0;
}

struct arg_struct {
    int start_ind;
    int end_ind;
    int thread_id;
    particle *hhparticles;
};

void *massChecktgv(void *arguments) {
    struct arg_struct *args = arguments;
    int start = args->start_ind;
    int end = args->end_ind;
    int id = args->thread_id;
    particle *hparticles = args->hhparticles;
    for (int j = start; j <= end; j += 1) {
        double ratio = hparticles[j].m_d / hparticles[j].initial_mass;
        if (ratio > 0.0001) {
//            printf("id = %d, start = %d, end = %d, index = %d, mass = %e, ratio = %e\n", id, start, end, j, hparticles[j].m_d, ratio);
            nexttgv += 1;
            break;
        }
    }
}

// Split mass check into parallel groups to reduce GPU idle time
bool parallelMassChecktgv(particle *hparticles, cl_ulong NUMPART) {
    nexttgv = 0;
    int threadcount = 1;
    bool cont;
    pthread_t id[threadcount];
    struct arg_struct args;
    args.hhparticles = hparticles;

    for (int i = 0; i < threadcount; i += 1) {
        args.start_ind = (NUMPART / threadcount) * i;
        args.end_ind = (NUMPART / threadcount) * (i + 1);
        args.thread_id = i;
        pthread_create(&id[i], NULL, massChecktgv, (void *)&args);
//        pthread_join(id[i], NULL); // serial code to crosscheck
    }

    for (int i = 0; i < threadcount; i += 1) {
//        printf("Join %d\n", i);
        pthread_join(id[i], NULL);
    }
//    printf("Non zero threads = %d\n", nexttgv);
//    pthread_exit(NULL);

    if (nexttgv != 0) {
        cont = TRUE;
    } else {
        cont = FALSE;
    }
    return cont;
}