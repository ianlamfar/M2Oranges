#ifndef DEFAULT_PARTICLE_GRAVITY
float3 get_gravity(particle p, float delta_t) {
    return (float3) {0, -9.81, 0};
}
#endif

#ifndef DEFAULT_PARTICLE_FLUID_VEL
float3 get_vel_fluid(particle p, float delta_t) {
    return (float3) {0, 0, 0};
}
#endif

float3 get_accel(particle p, float delta_t, double tau_v) {
    float tau = (float) tau_v;
    float3 non_drag_a = get_gravity(p, delta_t);
    return (get_vel_fluid(p, delta_t) - p.vel + tau * non_drag_a) / (tau + delta_t);
}

float3 iterate_velocity(particle p, float delta_t, double tau) {
    float3 next_vel = p.vel + delta_t * get_accel(p, delta_t, tau);
    return next_vel;
}

float3 iterate_position(particle p, float delta_t, float3 next_vel) {
    return p.pos + delta_t * (next_vel + p.vel) / 2;
}

float get_u_s(particle p, float delta_t){
    float3 fluid_vel = get_vel_fluid(p, delta_t);
    float3 rel_vel = {fluid_vel.x-p.vel.x, fluid_vel.y-p.vel.y, fluid_vel.z-p.vel.z};
    return sqrt(rel_vel.x * rel_vel.x + rel_vel.y * rel_vel.y + rel_vel.z * rel_vel.z);
}

double get_Re(particle p, float delta_t) {
    double u_s = get_u_s(p, delta_t);
    double mu_G = p.fluid_viscosity;
    double rho_G = p.rho_G;
    return (rho_G*u_s*p.diameter)/(mu_G);
}

float get_tau(particle p) {
    return p.density * p.diameter * p.diameter / (18.0 * p.fluid_viscosity);
}

// complete M2 calculation
struct output {
    // Items to return: m_d, diameter, pos, vel, T_d, fluid_velocity, force
    double pmass;
    double pdiameter;
    float3 ppos;
    float3 pvel;
    double ptemp;
    float3 fvel;
    float3 pforce;
};

struct output m2(particle p, float delta_t, int coupled, int analytic, int fixed_Re, int fixed_tau, float tau_scale) {
    struct output result;
    float3 fvel_return;
    float3 force;

    // printf("%e, %f, %v3f", p.m_d, p.T_d, p.vel);


    // decane
    double L_V = (3.958*pow(10.0, 4)) * (pow((619-p.T_d), 0.38));

    double x_s_eq = (p.P_atm/p.P_G) * exp(((L_V)/(p.R_bar/p.W_V))*((1/p.T_B)-(1/p.T_d)));
    double Y_s_eq = (x_s_eq) / (x_s_eq+(1-x_s_eq)*p.theta_2);
    double B_m_eq = (Y_s_eq-p.Y_G) / (1-Y_s_eq);
    double H_M = log(1.0 + B_m_eq);
    double F_M = (pow((1+B_m_eq), 0.7) * (log(1+B_m_eq))) / B_m_eq;

    double Y_bar = Y_s_eq + (1.0/3.0) * (p.Y_G - Y_s_eq);
    double T_bar = p.T_G + (1.0/3.0) * (p.T_G - p.T_d);
    double rho_g_bar = (Y_s_eq/p.density) + (pow(((1-Y_s_eq)/ p.rho_G), -1));
    double T_bar_star = T_bar / 1000;
    double C_p_V_bar;
    if (T_bar_star <= 0.8) {
        C_p_V_bar = 106.6 + (5765.0*T_bar_star) - (1675.0*(pow(T_bar_star, 2.0))) + (473.1*(pow(T_bar_star, 3.0)));
    } else {
        C_p_V_bar = 411.1 + (5460.0*T_bar_star) - (2483.0*(pow(T_bar_star, 2.0))) + 422.9*(pow(T_bar_star, 3.0));
    }
    double lambda_bar = (1.214*pow(10.0, -2)) * pow((T_bar/300), 1.8);
    double G_bar = (5.46*pow(10.0, -6)) * pow((T_bar/300), 1.583);
    double mu_bar = (5.64*pow(10.0, -6)) + ((1.75*pow(10.0,-8)) * (T_bar-300));
    double C_p_g_bar = (C_p_V_bar * Y_bar) + (p.C_p_G * (1 - Y_s_eq));
    double Le_bar = lambda_bar / (rho_g_bar * G_bar * C_p_g_bar);

    double Pr_G_bar = mu_bar * C_p_g_bar / lambda_bar;
    double Sc_G_bar = mu_bar / (rho_g_bar * G_bar);
    double Re_d;
    Re_d = get_Re(p, delta_t);
    // printf("Re = %f", Re_d);

    double Sh_0 = 2 + 0.552*pow(Re_d, 0.5)*pow(Sc_G_bar, (1.0/3.0));
    double Nu_0 = 2 + 0.552*pow(Re_d, 0.5)*pow(Pr_G_bar, (1.0/3.0));
    double Sh_star = 2 + ((Sh_0-2) / F_M);

    // B_T_pi iteration, output is B_T_pi
    double BT0 = pow(10.0, -3);
    double tol = pow(10.0, -12);
    double B_T_pi = BT0 + 1;
    double F_T = (pow((1+BT0), 0.7) / BT0) * log(1+BT0);
    double Nu_star = 2 + ((Nu_0 - 2) / F_T);
    double phi = (C_p_V_bar / C_p_g_bar) * (Sh_star / Nu_star) * (1 / Le_bar);
    B_T_pi = pow((1+B_m_eq), phi) - 1;

    while (fabs(B_T_pi - BT0) > tol) {
        BT0 = B_T_pi;
        F_T = (pow((1+BT0), 0.7) / BT0) * log(1+BT0);
        Nu_star = 2 + ((Nu_0 - 2) / F_T);
        phi = (C_p_V_bar / C_p_g_bar) * (Sh_star / Nu_star) * (1 / Le_bar);
        B_T_pi = pow((1+B_m_eq), phi) - 1;
        }

    float tau = get_tau(p);
    float tau_m = tau;
    float tau_h = tau;
    float tau_v = tau;
    if (fixed_tau == 0) { // Timescale simulation, 0 = unity tau, 1 = tau_m, 2 = tau_h, 3 = tau_v
        ;
    } if (fixed_tau == 1) {
        tau_m = tau * tau_scale;
    } if (fixed_tau == 2) {
        tau_h = tau * tau_scale;
    } if (fixed_tau == 3) {
        tau_v = tau * tau_scale;
    }

    double mass_return;
    double diameter_return;
    double mass_for_pow;
    double density_for_pow;
    double power_3_2 = 3.0/2.0;
    double power_2_3 = 2.0/3.0;
    double power_1_3 = 1.0/3.0;
    // iterate mass, output is mass_return, m_d_dot
    mass_return = (p.m_d)/(1.0 + (delta_t*(((Sh_star)/(3.0*Sc_G_bar))*(H_M/tau_m))));

    double m_d_dot = -((((Sh_star)/(3*Sc_G_bar))*(p.m_d/tau_m))*H_M);

    // iterate temp, output temp_return
    double f_2 = ((-m_d_dot/(p.m_d*B_T_pi)) * (3*Pr_G_bar*tau_h/Nu_star)); // correlation to heat transfer
    double temp_return = (p.T_d + delta_t*((((f_2*Nu_star)/(3.0*Pr_G_bar)) * (p.theta_1/tau_h)*(p.T_G)) + ((L_V*m_d_dot)/(p.C_L*p.m_d))-p.H_deltaT)) / (1.0+delta_t*(((f_2*Nu_star)/(3.0*Pr_G_bar)) * (p.theta_1/tau_h)));

    // velocity and position
    float3 next_vel = iterate_velocity(p, delta_t, tau_v);
    float3 next_pos = iterate_position(p, delta_t, next_vel);

    double b = 1.0/3.0;
    diameter_return = pow((mass_return*6.0)/(p.density*M_PI), b);

    fvel_return = get_vel_fluid(p, delta_t);
    force = (float3) {0, 0, 0};


    // Items to return: m_d, diameter, pos, vel, T_d, fluid_velocity, force
    result.pmass = mass_return;
    result.pdiameter = diameter_return;
    result.ppos = next_pos;
    result.pvel = next_vel;
    result.ptemp = temp_return;
    result.fvel = fvel_return;
    result.pforce = force;

    return result;

}


/* Kernel to iterate particles. Specified for M2 with analytic, coupled and non-fixed Re. */
__kernel void iterate_particle(__global particle *particles, float delta_t, int int_periodic, float domain_length,
        int coupled, int analytic, int fixed_Re, int model, int fixed_tau, float tau_scale) {
    bool periodic = (bool) int_periodic;
    int gid = get_global_id(0);

    struct output p = m2(particles[gid], delta_t, coupled, analytic, fixed_Re, fixed_tau, tau_scale);
    particles[gid].m_d = p.pmass;
    particles[gid].diameter = p.pdiameter;
    particles[gid].pos = p.ppos;
    particles[gid].vel = p.pvel;
    particles[gid].T_d = p.ptemp;
    particles[gid].fluid_vel = p.fvel;
    particles[gid].forces = p.pforce;

    if (periodic) {
        if (particles[gid].pos.x > domain_length / 2) {
            particles[gid].pos.x -= domain_length;
        } else if (particles[gid].pos.x < - domain_length / 2) {
            particles[gid].pos.x += domain_length;
        }

        if (particles[gid].pos.y > domain_length / 2) {
            particles[gid].pos.y -= domain_length;
        } else if (particles[gid].pos.y < - domain_length / 2) {
            particles[gid].pos.y += domain_length;
        }

        if (particles[gid].pos.z > domain_length / 2) {
            particles[gid].pos.z -= domain_length;
        } else if (particles[gid].pos.z < - domain_length / 2) {
            particles[gid].pos.z += domain_length;
        }
    }
}