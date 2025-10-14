#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rebound.h"

// Todos os dados removidos de Madeira et al. 2018

const double J2planet           = 16290.543820e-6;
const double J4planet           = -936.700366e-6;
const double Mplanet            = 5.68683765495e26;   
const double Rplanet            = 60330e3;    
const double G = 6.6743e-11 ; 

double tmax;    // tempo máximo em unidades do código

void heartbeat(struct reb_simulation* r);
void force_J2(struct reb_simulation* r);

double coefficient_of_restitution_constant(const struct reb_simulation* const r, double v){
    double cr = 0.1;
        return cr;
}

int main(int argc, char* argv[]){
    struct reb_simulation* r = reb_simulation_create();

    tmax           = 10.0*365.25*60.*60*24;          // Maximum integration time

//    reb_simulation_start_server(r, 1234);
    r->integrator  = REB_INTEGRATOR_IAS15;
    r->gravity          = REB_GRAVITY_BASIC;
    r->boundary         = REB_BOUNDARY_OPEN;
    r->G = G;  // Gravitational constant
    r->dt = 0.1*60.*60*24; // Initial timestep (s)
//    r->softening        = 1e-6;         // Gravitational softening length
    r->collision        = REB_COLLISION_TREE;
    r->coefficient_of_restitution = coefficient_of_restitution_constant;
	r->collision_resolve = reb_collision_resolve_hardsphere; // altered by create_sim 

    const double boxsize = 1000.0*Rplanet; // Distancia de ejecao
    reb_simulation_configure_box(r,boxsize,1,1,1);

    // Saturno
    struct reb_particle saturno = {0};
    saturno.m  = Mplanet;
    saturno.r  = Rplanet;
    saturno.hash = 0;
    reb_simulation_add(r, saturno);

    // Mimas
    double mass, radius, a, e, inc, w, Omg, lmd;
    mass = 3.754e19;
    radius = pow(3.*mass/(4.*M_PI*1152.),1./3.);
    a = 1.8600466879e8;
    e = 1.7245219209e-2;
    inc = 1.5687571620*M_PI/180.;
    Omg = 259.15258436*M_PI/180.;
    w =  163.18023984*M_PI/180.-Omg;
    lmd = reb_M_to_f(e, 197.73278953*M_PI/180.-Omg-w);
    struct reb_particle mimas = reb_particle_from_orbit(r->G, saturno, mass, a, e, inc, Omg, w, lmd);
    mimas.hash = 1;
    mimas.r = radius;
    reb_simulation_add(r,mimas);

    // Aegaeon
    mass =  5.997e10;
    radius = pow(3.*mass/(4.*M_PI*500.),1./3.);
    a = 1.68005e8; //1.6803396819e8;
    e = 0.3121285727e-2;
    inc = 0.0014761087*M_PI/180.;
    Omg = 233.02211168*M_PI/180.;
    w =  145.76280592*M_PI/180.-Omg;
    lmd = reb_M_to_f(e, 5.3594249045*M_PI/180.-Omg-w);
    struct reb_particle aegaeon = reb_particle_from_orbit(r->G, saturno, mass, a, e, inc, Omg, w, lmd);
    aegaeon.hash = 2;
    aegaeon.r = radius;
    reb_simulation_add(r,aegaeon);

    reb_simulation_move_to_com(r);

    // Setup callback functions
    r->N_active         = r->N;            

// Se for incluir particulas, elas devem entrar aqui

    r->heartbeat        = heartbeat;
    r->additional_forces    = force_J2;

    reb_simulation_integrate(r, tmax);

    reb_simulation_free(r);
    printf("\n");

}

void force_J2(struct reb_simulation* r){
    if (J2planet==0&&J4planet==0) return;
    // Star
    int i,k;
    k=0;
    for (i=0; i<r->N; i++){
        if (r->particles[i].m>=Mplanet) k = i;
    }
    const struct reb_particle planet = r->particles[k];     // cache
    const int N = r->N;
#pragma omp parallel for
    for (i=0;i<N;i++){
        if (i==k) continue ;
        const struct reb_particle p = r->particles[i];      // cache
        const double prx = p.x-planet.x;
        const double pry = p.y-planet.y;
        const double prz = p.z-planet.z;
        const double pr2  = prx*prx + pry*pry + prz*prz;        // distance^2 relative to planet
        const double pr = sqrt(pr2);        
        const double costheta = prz/pr;
        const double costheta2 = costheta*costheta;

        const double f1 = 3.0/2.0*G*Mplanet*J2planet*Rplanet*Rplanet/pr2/pr2/pr;
        const double f2 = 5.0*costheta2 - 1.0;
        const double f3 = f2 - 2.0;

        const double f4 = 5.0/8.0*G*Mplanet*J4planet*Rplanet*Rplanet*Rplanet*Rplanet/pr2/pr2/pr2/pr;
        const double f5 = 63.0*costheta2*costheta2 - 42.0*costheta2 + 3.0;
        const double f6 = f5 - 28.0*costheta2 + 12.0;

        const double pax  = f1*f2*prx+f4*f5*prx;
        const double pay  = f1*f2*pry+f4*f5*pry;
        const double paz  = f1*f3*prz+f4*f6*prz;

        r->particles[i].ax += pax;
        r->particles[i].ay += pay;
        r->particles[i].az += paz;

        const double fac = p.m/planet.m;

        r->particles[k].ax -= fac*pax;
        r->particles[k].ay -= fac*pay;
        r->particles[k].az -= fac*paz;

    }
}

void heartbeat(struct reb_simulation* r){
    reb_simulation_update_tree(r);
    reb_simulation_move_to_com(r);
    if(reb_simulation_output_check(r, 10.*r->dt)){          
        reb_simulation_output_timing(r, tmax);
        int i,k;
        for (i=0; i<r->N; i++){
            if (r->particles[i].m>=Mplanet) k = i;
        }
        const struct reb_particle sat = r->particles[k];
        for (i=0; i<r->N; i++){
            if (k==i) continue;
            const struct reb_particle p = r->particles[i];
            struct reb_orbit o = reb_orbit_from_particle(r->G,p,sat);
            //get orbital elements
            if (p.hash == 1){
                FILE* fm = fopen("mimas.aei","ab");
                fprintf(fm,"%.15e\t%.15e\t%.15e\t%.15e\t%.15e\t%.15e\t%.15e\t%.15e\t%.15e\t%.15e\n",r->t,p.x-sat.x,p.y-sat.y,p.z-sat.z,p.vx-sat.vx,p.vy-sat.vy,p.vz-sat.vz,p.m,o.a,o.e);
                fclose(fm);
            }
            if (p.hash == 2){
                FILE* fm = fopen("aegaeon.aei","ab");
                fprintf(fm,"%.15e\t%.15e\t%.15e\t%.15e\t%.15e\t%.15e\t%.15e\t%.15e\t%.15e\t%.15e\n",r->t,p.x-sat.x,p.y-sat.y,p.z-sat.z,p.vx-sat.vx,p.vy-sat.vy,p.vz-sat.vz,p.m,o.a,o.e);
                fclose(fm);
            }
        }
    }
}

