//
// Created by Ian on 16/2/2021.
//
#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else

#include <CL/cl.h>

#endif

#include <stdio.h>
#include "../../../util/clUtils/clUtils.h"
#include "../../../util/particleUtils/particleUtils.h"
#include "../../../tests/run_tests/run_tests.h"
#include "../../../util/simUtils/simUtils.h"
#include "../../../sims/simRunner/simRunner.h"
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
cl_ulong NUMPART = 1e5;

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
cl_int fixed_tau = 2; // Timescale simulation, 0 = unity tau, 1 = tau_m, 2 = tau_h, 3 = tau_v
cl_int model = 2;
cl_float tau_scale;
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

int start;
int end;



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
    // ONLY WORKS FOR M2
    cl_kernel iterate_particle = getKernel((cl_device_id) device, context, iterate_particle_files, 4, "iterate_particle", TRUE);



    hparticles = malloc(sizeof(particle) * NUMPART);
    if (hparticles == NULL) {
        fprintf(stderr, "Particles memory allocation failed.\n");
        return 1;
    }

    mu_G = 6.109e-6 + 4.604e-8 * T_R - 1.051e-11 * pow(T_R, 2);
    //    mu_G = 6.109e-6 + 4.604e-8 * T_G - 1.051e-11 * pow(T_G, 2);

    float min_scale = 0.5;
    float max_scale = 5.0;
    float step = 0.5;
    int num = (max_scale - min_scale) / step;
    float scales[num];
    for (int i = 0; i <= num; i += 1) {
        float out = min_scale + i * step;
        scales[i] = out;
    }

    for (int i = 0; i <= num; i += 1) {
        tau_scale = min_scale + (i * step);
        printf("Tau scale = %f\n", tau_scale);
        printf("[INIT] Creating particle positions.\n");
        particle_effect_diameter = (cl_float) (1.5 * particle_diameter);
        cl_float3 *positions = malloc(sizeof(cl_float3) * NUMPART);
        // Using particle_effect_diameter so that cohesion effects are considered at the appropriate range.
        float cube_length = createCubePositions(positions, NUMPART, particle_effect_diameter, 2, (cl_float3) {0, 0, 0});
        domain_length = (cl_float) (2 * PI);

        if (cube_length > domain_length) {
            fprintf(stderr, "Not all particles fit within the specified domain length for the given cube parameters (%.3f > %.3f).", cube_length, domain_length);
            return 1;
        }

        cl_float3 *velocities = malloc(sizeof(cl_float3) * NUMPART);
        createNormalDistVelocities(velocities, NUMPART, init_speed_mean, init_speed_std_dev);

        particle_diameter = 0.002;

        // Initialize particles.
        initializeMonodisperseParticles(hparticles, NUMPART, density, mu_G, particle_diameter,
                                        particle_effect_diameter, C_p_G, C_L, P_atm, rho_G, R_bar, R, W_G, W_V, Y_G,
                                        T_d, T_B, T_G, positions, velocities);
        free(positions);

        if (!checkPositions(hparticles, NUMPART, domain_length)) {
            fprintf(stderr, "Particles outside domain limits.\n");
            return 1;
        }
        tau = get_tau(&(hparticles[0]));
        timestep = 0.0001 * tau;
        log_step = 100*timestep;
        sim_length = 100 * tau;

        // Convert decimal scale values xx.xx to filename format xx_xx
        char s0[20];
        int j = 0;
        sprintf(s0, "%.2f", scales[i]);
        char * splits = strtok(s0, ".");
        char * array[2];
        while (splits != NULL) {
            array[j++] = splits;
            splits = strtok(NULL, ".");
        }
        char sc[20];
        sprintf(sc, "%s_%s", array[0], array[1]); // combine ones and decimals into ones_decimals

        // create dir string
        snprintf(folder, sizeof(folder), "%s%s%s", "c_heat_mass_transfer_tau_h_", sc, "/");
//        printf("%s\n", folder);
        sprintf(dir, "%s%s", PROJECT_DIR "analysis/timescales/tau_h/data/", folder);
        sprintf(datadir, "%s", PROJECT_DIR "analysis/timescales/tau_h/data/");

        // calculate simulation duration
        start = time(NULL);
        runSim(hparticles, NUMPART, iterate_particle, hparticles[0].diameter, periodic, domain_length,
               prefix, dir, sim_length, timestep, VERBOSE, LOG_DATA, TRUE, log_step, (cl_device_id) device, context, coupled,
               analytic, fixed_Re, model, fixed_tau, tau_scale);
        end = time(NULL);

        int minutes = (end-start) / 60;
        int seconds = (end-start) - (minutes * 60);
        printf("Simulation duration = %i m %i s\n\n", minutes, seconds);

        clReleaseContext(context);
    }

}

