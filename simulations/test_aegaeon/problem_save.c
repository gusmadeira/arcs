#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rebound.h"

// Todos os dados removidos de Madeira et al. 2018

const double J2planet           = 16290.543820e-6;        
const double Mplanet            = 5.68683765495e26;   
const double Rplanet            = 60330e3;    
const double ObliquityPlanet    = 0.;  

double tmax;    // tempo máximo em unidades do código

void heartbeat(struct reb_simulation* r);
void force_J2(struct reb_simulation* r);

double coefficient_of_restitution_constant(const struct reb_simulation* const r, double v){
    double cr = 0.1;
        return cr;
}

int main(int argc, char* argv[]){
    struct reb_simulation* r = reb_simulation_create();

    tmax           = 0.1*365.25*60.*60*24;          // Maximum integration time

//    reb_simulation_start_server(r, 1234);
    r->integrator  = REB_INTEGRATOR_IAS15;
    r->gravity          = REB_GRAVITY_TREE;
    r->boundary         = REB_BOUNDARY_OPEN;
    r->opening_angle2   = pow(0.5, 2);          // referente a arvore da GRAVIDADE
    r->G = 6.6743e-11;  // Gravitational constant
    r->dt = 60.*60*24; // Initial timestep (s)
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
    mass =  0.0*5.997e10;
    radius = pow(3.*mass/(4.*M_PI*500.),1./3.);
    a = 167693.73e3; //1.6803396819e8;
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
}

void force_J2(struct reb_simulation* r){
    if (J2planet==0) return;
    // Star
    const struct reb_particle planet = r->particles[0];     // cache
    const int N = r->N;
#pragma omp parallel for
    for (int i=1;i<N;i++){
        const struct reb_particle p = r->particles[i];      // cache
        const double sprx = p.x-planet.x;
        const double spry = p.y-planet.y;
        const double sprz = p.z-planet.z;
        const double prx  = sprx*cos(-ObliquityPlanet) + sprz*sin(-ObliquityPlanet);
        const double pry  = spry;
        const double prz  =-sprx*sin(-ObliquityPlanet) + sprz*cos(-ObliquityPlanet);
        const double pr2  = prx*prx + pry*pry + prz*prz;        // distance^2 relative to planet
        const double fac  = 3.*r->G*J2planet*planet.m*Rplanet*Rplanet/2./pow(pr2,3.5);
        
        const double pax  = fac*prx*(prx*prx + pry*pry - 4.*prz*prz);
        const double pay  = fac*pry*(prx*prx + pry*pry - 4.*prz*prz);
        const double paz  = fac*prz*(3.*(prx*prx + pry*pry) - 2.*prz*prz);

        r->particles[i].ax += pax*cos(ObliquityPlanet) + paz*sin(ObliquityPlanet);
        r->particles[i].ay += pay;
        r->particles[i].az +=-pax*sin(ObliquityPlanet) + paz*cos(ObliquityPlanet);
//        r->particles[0].ax -= (p.m/planet.m)*(pax*cos(ObliquityPlanet) + paz*sin(ObliquityPlanet));
//        r->particles[0].ay -= (p.m/planet.m)*(pay);
//        r->particles[0].az -=-(p.m/planet.m)*(pax*sin(ObliquityPlanet) + paz*cos(ObliquityPlanet));

    }
}

void heartbeat(struct reb_simulation* r){
    if(reb_simulation_output_check(r, 1.*r->dt)){          
        reb_simulation_output_timing(r, tmax);
        int i,k;
        for (i=0; i<r->N; i++){
            if (r->particles[i].m>=Mplanet) k = i;
        }
        const struct reb_particle sat = r->particles[k];
        for (i=0; i<r->N; i++){
            if (k==1) continue;
            const struct reb_particle p = r->particles[i];
            struct reb_orbit o = reb_orbit_from_particle(r->G,p,sat);
            //get orbital elements
            if (p.hash == 1){
                FILE* fm = fopen("mimas.txt","ab");
                fprintf(fm,"%.15e\t%.15e\t%.15e\t%.15e\t%.15e\t%.15e\t%.15e\t%.15e\t%.15e\t%.15e\n",r->t,p.x-sat.x,p.y-sat.y,p.z-sat.z,p.vx-sat.vx,p.vy-sat.vy,p.vz-sat.vz,p.m,o.a,o.e);
                fclose(fm);
            }
            if (p.hash == 2){
                FILE* fm = fopen("aegaeon.txt","ab");
                fprintf(fm,"%.15e\t%.15e\t%.15e\t%.15e\t%.15e\t%.15e\t%.15e\t%.15e\t%.15e\t%.15e\n",r->t,p.x-sat.x,p.y-sat.y,p.z-sat.z,p.vx-sat.vx,p.vy-sat.vy,p.vz-sat.vz,p.m,o.a,o.e);
                fclose(fm);
            }
        }
    }
    reb_simulation_update_tree(r);
    reb_simulation_move_to_com(r);
}

