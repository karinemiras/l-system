//
// Created by Karine Miras on 21/03/2017.
//

#ifndef LSYSTEM_PROTO_EVOLUTION_H
#define LSYSTEM_PROTO_EVOLUTION_H

#include <map>
#include <string>
#include <vector>

#include "Aux.h"
#include "Genome.h"
#include "Measures.h"
#include "Tests.h"

//from tol.spec import get_body_spec, get_brain_spec

/**
 * Evolutionary algorithms.
 */
class  Evolution{


public:

    explicit Evolution(std::string experiment_name,
                       std::string path){

        this->experiment_name = experiment_name;
        this->path = path;

        this->measures_names.push_back("branching");
        this->measures_names.push_back("connectivity1");
        this->measures_names.push_back("connectivity2");
        this->measures_names.push_back("coverage");
        this->measures_names.push_back("effective_joints");
        this->measures_names.push_back("length_ratio");
        this->measures_names.push_back("sensors");
        this->measures_names.push_back("symmetry");
        this->measures_names.push_back("total_components");

    }

    void saveHistory(int generation);

    void readParams();
    void testGeneticString(int argc,
                           char* argv[],
                           std::string test_genome);
    void measureIndividuals(int generation,
                            std::vector<Genome>  &individuals,
                            std::string dirpath);

    void evaluateLocomotion(int generation,
                     std::vector<Genome >  &individuals);
    void savesValidity(int generation);

    int  tournament();
    void selection();
    std::vector<Genome>  getPopulation();
    std::map<std::string, double> getParams();
    double runExperiment_part1(int generation, int load_experiment);
    double runExperiment_part2(int generation);
    void exportGenerationMetrics(int generation,
                                 std::vector<int> metrics);
    void saveLocomotionFitness(std::string genome_id, double fitness);
    void saveBalanceFitness(std::string genome_id, double fitness);
    void exportPop(int generation);
    void calculateNovelty();
    void calculateFinalFitness();
    void calculateRankFitness();
    void saveParameters();
    void logsTime(std::string moment);
    void setupEvolution();
    void writesEvolutionState(int generation, int next_id);
    std::vector<std::string> readsEvolutionState();
    void loadsParams();
    void loadIndividuals(int generation, std::string type);
    std::vector<int> calculateNicheCoverage();
    void createHeader();
    void updateParameter(std::string key, double value);
    void developIndividuals(int argc, char* argv[],
                            LSystem LS,
                            int generation,
                            std::vector<Genome>  &individuals,
                            std::string dir);
    int loadExperiment();
    void initExperiment(int argc, char* argv[],
                       LSystem LS);
    void summaryNicheCoverage();
    int getGeneration_genome(std::string idgenome);
    double compareIndividual(Measures m,
                             std::string idgenome);
    double compareParents(std::string idparent1,
                          std::string idparent2);
    void  addToArchive();
    void cleanMemory(std::vector< int > index_selected);
    void cleanVertex(DecodedGeneticString::Vertex * v);

    virtual void initPopulation(LSystem LS){};
    virtual void crossover(LSystem LS){};
    virtual void mutation(LSystem LS){};



protected:

    std::vector<std::string> measures_names =
            std::vector<std::string>();

    std::map<std::string, double> params =
            std::map<std::string, double>(); // contains the list of parameters loaded from parameter file

    int next_id = 0; // id that will be given for the next genome to be created

    std::string experiment_name = ""; // name for the experiment

    std::string path = ""; // path of the lsystem

    // points in a grid representing the morphological space
    std::map<std::string, std::vector<double>>
            morphological_grid_generation =
            std::map<std::string, std::vector<double>>();

    std::map<std::string, std::vector<std::string>>
            morphological_grid_accumulated =
            std::map<std::string, std::vector<std::string>>();

     // containsgeneral auxiliar methods for the experiments
    Aux aux = Aux(this->experiment_name,this->getParams(),this->path);
    // contains methods with tests for the system
    Tests tests = Tests(this->experiment_name,
                        this->getParams(),
                        this->path);

    // contains the genomes of all the individuals of the current population
    std::vector<Genome>  population =  std::vector<Genome>();

    // contains the genomes of all the individuals the new offspring
    std::vector<Genome>  offspring =  std::vector<Genome>();

    // contains the genomes of all individuals in the archive
    std::vector<Genome>  archive = std::vector<Genome> ();


};


#endif //LSYSTEM_PROTO_EVOLUTION_H
