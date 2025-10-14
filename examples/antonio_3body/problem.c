#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rebound.h"

#define PI 3.14159265358979323846
#define MASSA_SOL 1  //massa do Sol normalizada com 1
#define GRAVITACIONAL 1 // constante G normalizada com 1 também
#define SEMI_EIXO_MAIOR 1 // distancia Terra-SOl normalizada com 1
#define MASSA_TERRA 3e-4  // a massa da terra é 0,0003% do sol
#define EXC 0.0  // a massa da terra é 0,0003% do sol

// Função heartbeat
//chama essa função periodicamente durante a integração (cada "passo de tempo") para executar tarefas de monitoramento ou gravação de dados.
void heartbeat(struct reb_simulation* r){
    //verifica se é hora de gravar a saída da simulação.

    if (reb_simulation_output_check(r, r->dt)){ //grava cada uma hora de um ano. Vai rodar 8760 vezes
        reb_simulation_output_timing(r, 0); //Grava o tempo atual da simulação no arquivo de saída padrão do REBOUND.
        struct reb_orbit o = reb_orbit_from_particle(r->G,r->particles[1],r->particles[0]);
        FILE *f = fopen("saida.dat", "a"); //abre o arquivo e escrever na proxima linha
        if(f!=NULL){
            fprintf(f, "%e\t%e\t%e\t%e\t%e\t%e\t%e\n", r->t, o.a, o.e, o.inc, o.omega, o.Omega, o.l);
            fclose(f);
        }
    }
}

double coefficient_of_restitution_constant(const struct reb_simulation* const r, double v){
    // v is the normal impact velocity.
    // Here, we just use a constant coefficient of restitution
    return 1.0;
}

int main(int argc, char* argv[]){
    //criar e inicializar a simulação do rebound
    struct reb_simulation* r = reb_simulation_create(); //cria a simulação e retorna a struct com os campos (l: 416)

    r->G            = GRAVITACIONAL; //atribuir a constante gravitacional
    r->integrator   = REB_INTEGRATOR_IAS15; //método de integração
    //r->ri_whfast.coordinates = REB_WHFAST_COORDINATES_DEMOCRATICHELIOCENTRIC; //sistema de coordenadas centrado no Sol.
    r->boundary     = REB_BOUNDARY_OPEN; //sem fronteiras artificiais
    r->dt           = 0.001; //passo de tempo, pequeno para ter mais precisão
    r->softening    = 1e-6; //parâmetro para evitar singularidade em aproximações muito pequenas.
    r->collision    = REB_COLLISION_DIRECT;
	r->collision_resolve = reb_collision_resolve_merge; // altered by create_sim 
    r->coefficient_of_restitution = coefficient_of_restitution_constant;

// Caixa de simulação em unidades absolutas
    // 3 unidades de comprimento (1 unidade = distância Terra-Sol)
    reb_simulation_configure_box(r, 10.0, 1, 1, 1); //

    reb_simulation_start_server(r, 1234);
    r->usleep = 50000;

    //definições do sol
    struct reb_particle sol = {0}; //inciar os campos com zero
    sol.m = MASSA_SOL; //massa do sol em kg
    sol.r = 0.1; //massa do sol em kg

    reb_simulation_add(r, sol); //adicionar a particula sol  r->particles[0] = sol

    //definições da terra
    struct reb_particle terra = reb_particle_from_orbit(r->G, sol, 0.0001*MASSA_TERRA, 1.0, 0.0, 0.0, 0.0, 0.0, 60.0);
    terra.m  = 0.0001*MASSA_TERRA; //massa da terra em kg
    terra.r  = 5e-2;
    reb_simulation_add(r, terra); //adicionar particula da terra r->particles[1] = terra;

    struct reb_particle teia = reb_particle_from_orbit(r->G, sol, MASSA_TERRA, 1.0/(1-EXC), EXC, 0.0, 0.0, 0.0, 0.0);
    teia.m  = MASSA_TERRA; //massa da terra em kg
    teia.r  = 5e-2;
    reb_simulation_add(r, teia); //adicionar particula da terra r->particles[1] = terra;

    r->N_active     = r->N;    // o Sol e a Terra são os únicos com massa. O resto é apenas partícula teste.

    //Move todas as partículas para que o centro de massa do sistema esteja na origem e sem momento linear total. É essencial para estabilidade numérica.
    reb_simulation_move_to_com(r);

    FILE* f = fopen("saida.dat", "w");
    if(f != NULL){
        printf("arquivo criado!");
        fclose(f);
    }

    // Conecta a função heartbeat para gravar as popsições no arquivo
    r->heartbeat = heartbeat; //a cada simulação chamar a função heartbeat


    double tempo_total = 20000.0 * PI; // Período nesta unidade absoluta

    // Rodar a simulação até o tempo final
    reb_simulation_integrate(r, tempo_total);

    return 0;
}

/*
compilar:
            gcc sol_terra_AU.c -o sol_terra_AU -I. -L. -lrebound -lm


executar:
            export LD_LIBRARY_PATH=./src:$LD_LIBRARY_PATH
            ./sol_terra_AU

python compilar:
            python3 plot_AU.py
*/