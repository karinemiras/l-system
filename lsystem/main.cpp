/*
 * Genome indirect enconding: L-System
 * Author: Karine Miras
 * Created: 02/03/17
 */

#include<iostream>
#include <random>
#include "EvolutionIndirect.h"

using namespace std;


int main(int argc, char* argv[])
{

    std::random_device rd;
    std::default_random_engine generator(rd());
    std::uniform_int_distribution<int> dist_1(1, 1000000);

    EvolutionIndirect evolve_generation = EvolutionIndirect
        ("exp_"+std::to_string(dist_1(generator)),"../../");

    /* test setup */
      evolve_generation.setupEvolution();
      int load_generation = 0;
      int ini = 1;
    /* test setup */

    for(int i=ini; i <= evolve_generation.getParams()["pop_size"]; i++)
    {
        evolve_generation.runExperiment_evolve1(i, load_generation);
        evolve_generation.runExperiment_evolve2(i,
                                                evolve_generation.getParams()["learning_iterations"]);
        load_generation = 0;
    }

    return 0;
}


