//
// Created by Ian on 22/2/2021.
//

#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else

#include <CL/cl.h>

#endif

#include <stdio.h>
#include "../../util/clUtils/clUtils.h"
#include "../../util/particleUtils/particleUtils.h"
#include "../../tests/run_tests/run_tests.h"
#include "../../util/simUtils/simUtils.h"
#include "../../sims/simRunner/simRunner.h"
#include <malloc.h>
#include <string.h>

#define MAX_SOURCE_SIZE (0x100000)
#define VERBOSE FALSE
#define LOG_DATA TRUE

char prefix[50];
char folder[100];
char dir[200];
char datadir[200];

particle *hparticles;
cl_ulong NUMPART;

// Gas properties
cl_double C_p_G = 1141;
cl_double W_G = 28.97; // Molecular weight of carrier gas (air) (kg/(kg mole))
cl_double Y_G = 0;
cl_double T_G = 1000; // Gas temperature (K)
cl_double rho_G = 0.3529;

// Particle properties.
cl_double T_d = 315; //particle temperature (K)
cl_double density = 642;
cl_double particle_diameter = 0.002;
cl_double particle_effect_diameter;
cl_double fluid_viscosity;// = 1.04;//0.0000193
cl_double T_B = 447.7; // Boiling temperature
cl_double W_V = 142;  // Molecular weight of vapour phase (water) (kg/(kg mole))
cl_double C_L = 2520.5;  // (J/kg/K)
//cl_double Re_d = 17.0;

// Other properties
cl_double P_atm = 101325.0;  // Atmospheric pressure (Pa)
cl_double R_bar = 8314.5;
cl_double R = 287.0;  // Universal gas constant ()
cl_bool periodic = CL_TRUE;
cl_int coupled = 2;
cl_int analytic = 2;
cl_int fixed_Re = 0;
cl_int fixed_tau = 0; // Timescale simulation, 0 = unity tau, 1 = tau_m, 2 = tau_h, 3 = tau_v
cl_int model;
cl_float tau_scale = 1;
double T_R;
float tau;

float init_speed_mean = 1;
float init_speed_std_dev = 0.1;

cl_float timestep;
cl_float sim_length = 30;
cl_float log_step;

cl_float domain_length;

cl_context context;
cl_device_id *device;

double mu_G;

// time variables
int start;
int end;
struct timespec start1, end1;
char timefile[500];

int main() {
    // Initialize OpenCL.
    setContext((cl_device_id *) &device, &context, TRUE);

    // Run tests
    if (!run_all_tests((cl_device_id) device, context, FALSE)) {
        return 1;
    }

    // Build iterate_particle kernel.
    char *iterate_particle_files[] = {PROJECT_DIR "/util/kernelUtils.cl",
                                      PROJECT_DIR "/kernels/get_gravity/no_gravity.cl",
                                      PROJECT_DIR "/kernels/get_vel_fluid/tgv.cl",
                                      PROJECT_DIR "/kernels/iterate_particle.cl"};
    cl_kernel iterate_particle = getKernel((cl_device_id) device, context, iterate_particle_files, 4, "iterate_particle", TRUE);

    float min_scale = 2.0;
    float max_scale = 7.0; // check GPU VRAM, NUMPART 1e4-->5.12MB; 5e4-->25.6MB; ...; 1e7-->5.12GB; 5e7-->25.6GB
    float step = 1.0;
    int num = (max_scale - min_scale) / step;
    float scales[num];
    for (int i = 0; i <= num; i += 1) {
        float out = min_scale + i * step;
        scales[i] = out;
    }
    sprintf(timefile, PROJECT_DIR "analysis/m1m2_perf/data/times.txt");
    printf("%s\n", timefile);
    FILE *fd = fopen(timefile, "w");
    fprintf(fd, "Computational time\n");
    fclose(fd);

    bool five = FALSE;
    for (int i = 0; i <= num; i += 1) { // numpart sequence of 10, 50, 100, 500, 1000, 5000, ....
        float numpart_index = min_scale + (i * step);
        cl_ulong NUMPART;
        if (five) {
            NUMPART = 5 * (cl_ulong) pow(10, numpart_index);
            five = FALSE;
        } else {
            NUMPART = (cl_ulong) pow(10, numpart_index);
            five = TRUE;
            i -= 1;
        }
        if (NUMPART > pow(10, max_scale)) { // reached max numpart possible for GPU
            break;
        }

        printf("Number of particles = %llu\n", NUMPART);
        printf("Particle memory = %lu MB\n", sizeof(particle) * NUMPART / 1000 / 1000);
        hparticles = malloc(sizeof(particle) * NUMPART);
        if (hparticles == NULL) {
            fprintf(stderr, "Particles memory allocation failed.\n");
            return 1;
        }
//        printf("%d\n", sizeof(hparticles));
        mu_G = 6.109e-6 + 4.604e-8 * T_R - 1.051e-11 * pow(T_R, 2);
        //    mu_G = 6.109e-6 + 4.604e-8 * T_G - 1.051e-11 * pow(T_G, 2);


        // loop through each model
        for (int k = 1; k <= 2; k += 1) {
            model = k;
            printf("Model = %d\n", model);
            printf("[INIT] Creating particle positions.\n");
            particle_effect_diameter = (cl_float) (1.5 * particle_diameter);
            cl_float3 *positions = malloc(sizeof(cl_float3) * NUMPART);
            // Using particle_effect_diameter so that cohesion effects are considered at the appropriate range.
            float cube_length = createCubePositions(positions, NUMPART, particle_effect_diameter, 2,
                                                    (cl_float3) {0, 0, 0});
            domain_length = (cl_float) (2 * PI);

            if (cube_length > domain_length) {
                fprintf(stderr,
                        "Not all particles fit within the specified domain length for the given cube parameters (%.3f > %.3f).",
                        cube_length, domain_length);
                return 1;
            }

            cl_float3 *velocities = malloc(sizeof(cl_float3) * NUMPART);
            createNormalDistVelocities(velocities, NUMPART, init_speed_mean, init_speed_std_dev);

            particle_diameter = 0.002;

            // Initialize particles.
            initializeMonodisperseParticles(hparticles, NUMPART, density, mu_G, particle_diameter,
                                            particle_effect_diameter, C_p_G, C_L, P_atm, rho_G, R_bar, R, W_G, W_V, Y_G,
                                            T_d, T_B, T_G, positions, velocities);
            free(velocities);
            free(positions);

            if (!checkPositions(hparticles, NUMPART, domain_length)) {
                fprintf(stderr, "Particles outside domain limits.\n");
                return 1;
            }
            tau = get_tau(&(hparticles[0]));
            float tratio = 0.0001;
            timestep = tratio * tau;
            int lratio = 100;
            log_step = lratio * timestep;
            sim_length = 100 * tau;

            // Convert NUMPART to string for folder string
            char sc[20];
            sprintf(sc, "%d", NUMPART); // combine ones and decimals into ones_decimals

            // create dir string
            if (model == 1) {
                snprintf(folder, sizeof(folder), "%s%s%s", "[DUMMY DATA] c_heat_mass_transfer_m1_", sc, "/");
            }
            if (model == 2) {
                snprintf(folder, sizeof(folder), "%s%s%s", "[DUMMY DATA] c_heat_mass_transfer_m2_", sc, "/");
            }
            sprintf(dir, "%s%s", PROJECT_DIR "analysis/m1m2_perf/data/", folder);
            sprintf(datadir, "%s", PROJECT_DIR "analysis/m1m2_perf/data/");

            // calculate simulation duration
            start = time(NULL);
            clock_gettime(CLOCK_MONOTONIC, &start1);
            runSim(hparticles, NUMPART, iterate_particle, hparticles[0].diameter, periodic, domain_length,
                   prefix, dir, sim_length, timestep, VERBOSE, LOG_DATA, TRUE, log_step, (cl_device_id) device, context,
                   coupled,
                   analytic, fixed_Re, model, fixed_tau, tau_scale);
            clock_gettime(CLOCK_MONOTONIC, &end1);
            end = time(NULL);

//            int minutes = (end1.tv_sec-start1.tv_sec) / 60;
            double seconds = (end1.tv_sec - start1.tv_sec) + (end1.tv_nsec - start1.tv_nsec) / 1e9;
            printf("Simulation duration = %f s\n\n", seconds);

            FILE *fd = fopen(timefile, "a");
            fprintf(fd, "%d,%f,%d,%d,%f\n", model, tratio, lratio, (int) NUMPART, seconds);
            fclose(fd);


            clReleaseContext(context);
        }

    }
    free(hparticles);

}