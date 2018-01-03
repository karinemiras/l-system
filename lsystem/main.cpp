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

    EvolutionIndirect evolve_generation = EvolutionIndirect("test","../../");
    evolve_generation.setupEvolution();
    evolve_generation.runExperiment_part1(1, 1);

    return 0;
}


