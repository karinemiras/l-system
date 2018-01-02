//
// Created by Karine Miras on 21/03/2017.
//

#include <algorithm>
#include <ctime>
#include <fstream>
#include <iostream>
#include <map>
#include <random>
#include <string>
#include <sstream>
#include <thread>
#include <vector>

#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>


#include <mlpack/core.hpp>
#include <mlpack/methods/neighbor_search/neighbor_search.hpp>

using namespace mlpack;
using namespace mlpack::neighbor; // NeighborSearch and NearestNeighborSort
using namespace mlpack::metric; // EuclideanDistance

#include "Aux.h"
#include "Evolution.h"
#include "Genome.h"
#include "LSystem.h"
#include "Measures.h"


/**
 * Reads parameters from file.
 **/
void Evolution::readParams()
{
  std::string line;
  std::ifstream myfile(this->path+"lsystem/configuration.txt");

  /*   pop_size - size of the population of genomes
  *    offspring_prop - proportion of the population size to dcalculate size of offspring
  *    num_initial_comp - number of initial (random) components in the production rules of the grammar
  *    show_phenotypes - flag to show the phenotype graphic
  *    export_phenotypes - if exports the phenotypes to images (1) or not (0)
  *    export_genomes - if exports the genomes to files (1) or not (0)
  *    replacement_iterations - number of replacement iterations for the l-system
  *    size_component - size of each component in pixels
  *    spacing - spacing between components in pixels
  *    num_generations - number of generations of the evolution
  *    mutation_prob - probability of adding/removing/swaping items (letters/commands) to the genetic-string in the mutation
  *    max_comps - maximum number of components allowed per phenotype
  *    prob_add_archive - probability of adding any genome to the archive
  *    grid_bins - number of bins to break morphological space
  *    logs_to_screen - if exports the logs to the screen (1) or not (0)
  *    logs_to_file - if exports logs to a file (1) or not (0)
  *    vizualize_simulation - if gazebo will be open during simulation (1) or not (0)
  *    new_experiment - if it is a new experiment (1) or a restore (0)
  */

  if (myfile.is_open())
  {
    while (getline(
        myfile,
        line))
    {
      std::vector< std::string > tokens;

      // parameters label and value separated by space
      boost::split(
          tokens,
          line,
          boost::is_any_of(" "));

      // first item is the label, second is the value
      this->params[tokens[0]] = std::stod(tokens[1]);
    }
    myfile.close();
  }
  else
  {
    this->aux.logs("Unable to open parameters file.");
  }
}

/*
 * Change value of parameter for an experiment.
 * */
void Evolution::updateParameter(
    std::string key,
    double value)
{

  this->params[key] = value;
}


/**
 * Loads parameters from saved state.
 **/
void Evolution::loadsParams()
{

  std::string line;
  std::ifstream myfile(
      this->path+"experiments/"
      + this->experiment_name +
      "/configuration.txt");

  if (myfile.is_open())
  {
    while (getline(
        myfile,
        line))
    {

      std::vector< std::string > tokens;
      // parameters label and value separated by space
      boost::split(
          tokens,
          line,
          boost::is_any_of(" "));

      // first item is the label, second is the value
      this->params[tokens[0]] = std::stod(tokens[1]);
    }
    myfile.close();
  }
  else
  {
    this->aux.logs("Unable to open parameters state.");
  }

}



/* Finds out in which generation the genome was generated.
 * @param idgenome - id of the genome for which to verify generation.
 * */
int Evolution::getGeneration_genome(std::string idgenome)
{

  int generation_genome = 0;
  int offspring_size =
      this->params["pop_size"] * this->params["offspring_prop"];

  // generation of the genome can be found by its id, considering the size of the population and the offspring
  if (this->params["offspring_prop"] == 1)
  {
    generation_genome = (int) trunc(
        std::stoi(idgenome)
        / this->params["pop_size"]) + 1;
  }
  else
  {
    generation_genome = (int) trunc((std::stoi(idgenome) - offspring_size)
                                    / offspring_size) + 1;
  }

  if (generation_genome == 0)
  { generation_genome = 1; }

  return generation_genome;

}


/**
*  Copies the phenotypes of the selection population to a separate folder.
*  @param generation - number of generation in evolution
**/
void Evolution::exportPop(int generation)
{
  std::ofstream measures_file_general;
  std::string path =
      this->path+"experiments/" + this->experiment_name + "/measures2.txt";
  measures_file_general.open(
      path,
      std::ofstream::app);

  for (int i = 0;
       i < this->population.size();
       i++)
  { // for each genome in the population


    // finds number of generation to which the genome belongs to
    int generation_genome = this->getGeneration_genome(this->population[i].getId());

    std::string filename = "/body_" + this->population[i].getId()
                           + "_p1_" +
                           this->population[i].getId_parent1() +
                           "_p2_" +
                           this->population[i].getId_parent2() +
                           ".png";

    std::string pathfrom = this->path+"experiments/"
                           + this->experiment_name + "/offspringpop" +
                           std::to_string(generation_genome);

    std::string pathto = this->path+"experiments/"
                         + this->experiment_name + "/selectedpop" +
                         std::to_string(generation);

    // copies phenotype file from offspring folder to selected population folder
    system(("exec cp " + pathfrom + filename + " " + pathto +
            filename).c_str());

    // copies values of metrics to file of selected population
    std::string line;
    std::ifstream measures(
        this->path+"experiments/" + this->experiment_name +
        "/offspringpop" +
        std::to_string(generation_genome)
        + "/measures" +
        this->population[i].getId() + ".txt");

    while (getline(
        measures,
        line))
    {

      std::vector< std::string > tokens;
      boost::split(
          tokens,
          line,
          boost::is_any_of(":"));

      measures_file_general << std::to_string(generation) << " "
                            << this->population[i].getId() << " "
                            << tokens[0] << " " << tokens[1] << std::endl;
    }
  }

  measures_file_general.close();
}


/*
 * Compare phenotype of the individual with its parent's.
 * */
double Evolution::compareIndividual(
    Measures m,
    std::string idgenome)
{


  int generation_genome = this->getGeneration_genome(idgenome);

  std::string line;
  std::ifstream measures(
      this->path+"experiments/" + this->experiment_name +
      "/offspringpop" + std::to_string(generation_genome) +
      "/measures" + idgenome + ".txt");

  double dif = 0;
  while (getline(
      measures,
      line))
  {

    std::vector< std::string > tokens;
    boost::split(
        tokens,
        line,
        boost::is_any_of(":"));

    dif += std::pow(
        m.getGenome()->getMeasures()[tokens[0]] - std::stod(tokens[1]),
        2);
  }
  dif = roundf(std::sqrt(dif) * 100) / 100;

  return dif;

}


/*
 * Compare phenotype of the parents.
 * */
double Evolution::compareParents(
    std::string idparent1,
    std::string idparent2)
{


  int generation_genome_parent1 = this->getGeneration_genome(idparent1);
  std::string line;
  std::ifstream measures(
      this->path+"experiments/" + this->experiment_name +
      "/offspringpop" +
      std::to_string(generation_genome_parent1) +
      "/measures" + idparent1 + ".txt");

  int generation_genome_parent2 = this->getGeneration_genome(idparent2);
  std::string line2;
  std::ifstream measures2(
      this->path+"experiments/" + this->experiment_name +
      "/offspringpop" +
      std::to_string(generation_genome_parent2) +
      "/measures" + idparent2 + ".txt");

  double dif = 0;
  while (getline(
      measures,
      line))
  {

    getline(
        measures2,
        line2);

    std::vector< std::string > tokens, tokens2;
    boost::split(
        tokens,
        line,
        boost::is_any_of(":"));
    boost::split(
        tokens2,
        line2,
        boost::is_any_of(":"));

    dif += std::pow(
        std::stod(tokens[1]) - std::stod(tokens2[1]),
        2);
  }
  dif = roundf(std::sqrt(dif) * 100) / 100;

  return dif;

}

/**
 * Measures all the individuals of the population for several metrics.
 *  @param argc - command line parameter
 *  @param argv[] - command line parameter
 *  @param individuals - array with genomes
 *  @param dirpath - name of the output directory
 **/
void Evolution::measureIndividuals(
    int generation,
    std::vector< Genome > &individuals,
    std::string dirpath)
{

  std::ofstream differences_file;
  std::string path =
      this->path+"experiments/" + this->experiment_name + "/differences.txt";
  differences_file.open(
      path,
      std::ofstream::app);


  // for each genome of the population
  for (int i = 0;
       i < individuals.size();
       i++)
  {

    Measures m = Measures(
        this->experiment_name,
        this->params,
        this->path
    );
    m.setGenome(individuals[i]);
    // measures phenotype
    m.measurePhenotype(
        this->params,
        dirpath,
        generation);

    // compares measures between individuals
    if (individuals[i].getId_parent1() != "N")
    {

      double dif = this->compareIndividual(
          m,
          individuals[i].getId_parent1());
      differences_file << individuals[i].getId() << " " << dif;

      dif = this->compareIndividual(
          m,
          individuals[i].getId_parent2());
      differences_file << " " << dif;

      dif = this->compareParents(
          individuals[i].getId_parent1(),
          individuals[i].getId_parent2());
      differences_file << " " << dif << std::endl;

    }
  }

  differences_file.close();
}



/*
 *
 * */




/**
 * Creates files of results containing headers.
 */

void Evolution::createHeader()
{

  std::ofstream file;

  std::string path =
      this->path+"experiments/" + this->experiment_name + "/history.txt";
  file.open(path);
  file
      << "generation idgenome fitgenome idparent1 fitparent1 idparent2 fitparent2 meandif"
      << std::endl;
  file.close();

  path = this->path+"experiments/" + this->experiment_name + "/evolution.txt";
  file.open(path);
  file
      << "generation bestgenome maxfitness meanfitness "
          "nichecoverage_generation nichecoverage_accumulated";
  for (int i = 0;
       i < this->measures_names.size();
       i++)
  {
    file << " " << this->measures_names[i];
  }
  file << std::endl;
  file.close();

  path = this->path+"experiments/" + this->experiment_name + "/measures.txt";
  file.open(path);
  file << "generation idgenome";
  for (int i = 0;
       i < this->measures_names.size();
       i++)
  {
    file << " " << this->measures_names[i];
  }
  file << std::endl;
  file.close();

  path = this->path+"experiments/" + this->experiment_name + "/measures2.txt";
  file.open(path);
  file << "generation genome measures value" << std::endl;
  file.close();

  path = this->path+"experiments/" + this->experiment_name + "/differences.txt";
  file.open(path);
  file << "idgenome difference_parent1 difference_parent2 difference_parents"
       << std::endl;
  file.close();


}

void Evolution::saveHistory(int generation)
{

  for (int i = 0;
       i < this->offspring.size();
       i++)
  {

    std::ofstream history_file;
    std::string path =
        this->path + "experiments/" + this->experiment_name + "/history.txt";
    history_file.open(
        path,
        std::ofstream::app);

    history_file << std::to_string(generation) << " "     // generation
                 << this->offspring[i].getId() << " "   // idgenome
                 << this->offspring[i].getFitness() << " "  // fitness genome
                 << this->offspring[i].getId_parent1() << " "  // id of parent1
                 << this->offspring[i].getFit_parent1() << " "  // fitness ofparent1
                 << this->offspring[i].getId_parent2() << " " // id of parent2
                 << this->offspring[i].getFit_parent2() << " " // fitness of parent2
                 << (this->offspring[i].getFitness() - this->offspring[i].getFit_parent1())
                    + (this->offspring[i].getFitness() - this->offspring[i].getFit_parent2())
                      / (float) 2// mean dif fitness from parents
                 << std::endl;

    history_file.close();
  }
}



void Evolution::savesValidity(int generation)
{
  // saves list of robots and its validits to file to be read by simulator
  std::ofstream file;
  std::string path =
      this->path + "experiments/" + this->experiment_name +
      "/offspringpop" + std::to_string(generation) + "/validity_list.txt";
  file.open(path);
  for (int i = 0;
       i <  this->offspring.size();
       i++)
  {
    file << this->offspring[i].getId() << " " <<  this->offspring[i].getValid()
         << std::endl;
  }
  file.close();
}


/**
 * Compares average distance among points inthe accumulaetd grid.
 */

void Evolution::compareIndividuals(int generation)
{

  std::ofstream file;
  std::string path = this->path+"experiments/" + this->experiment_name +
                     "/nichecoverage_distances.txt";
  file.open(
      path,
      std::ofstream::app);


  // fetches points with its dimensions
  std::vector< std::vector< double>> points;
  for (const auto &it : this->morphological_grid_accumulated)
  {

    std::vector< double > dimensions;

    std::vector< std::string > tokens;
    boost::split(
        tokens,
        it.first,
        boost::is_any_of("|"));
    std::vector< std::string > tokens2(
        tokens.begin(),
        tokens.begin() + tokens.size() - 1);

    for (int i = 0;
         i < tokens2.size();
         i++)
    {
      dimensions.push_back(
          std::stod(tokens2[i]) / this->params["grid_bins"]);
    }
    points.push_back(dimensions);
  }


  double avgdistance = 0;
  std::vector< double > avgdistance_points;
  double stddev_distance = 0;

  // for each point
  for (int i = 0;
       i < points.size();
       i++)
  {

    double avgdistance_point = 0;
    // compare to every other point
    for (int j = 0;
         j < points.size();
         j++)
    {

      // for each dimension
      double distance = 0;
      for (int d = 0;
           d < this->population[0].getMeasures().size();
           d++)
      {

        distance += std::pow(
            points[i][d] - points[j][d],
            2);
      }
      // euclidean distance
      distance = std::sqrt(distance);

      avgdistance_point += distance;
    }
    // average distance from the point to all others
    avgdistance_point /= (points.size() - 1);

    avgdistance_points.push_back(avgdistance_point);
    avgdistance += avgdistance_point;
  }

  // average distance of all points
  avgdistance /= points.size();

  for (int i = 0;
       i < avgdistance_points.size();
       i++)
  {
    stddev_distance += std::pow(
        avgdistance - avgdistance_points[i],
        2);
  }
  stddev_distance = std::sqrt(stddev_distance);


  file << generation << " " << avgdistance << " " << stddev_distance
       << std::endl;


  file.close();

}


/**
 * Selects two random genomes and compares their fitness, choosing the winner.
 * @return - the index of the winner genome
 */

int Evolution::tournament()
{

  std::random_device rd;
  std::default_random_engine generator(rd());
  std::uniform_int_distribution< int > dist_1(
      0,
      (int) this->population.size() -
      1); // size of current pop (parents+offspring)

  int genome1 = dist_1(generator); // random genome 1
  int genome2 = dist_1(generator); // random genome 2

  //random selection test
  //  return dist_1(generator);


  // return the genome with higher fitness / novelty search

  if (this->population[genome1].getFitness() >
      this->population[genome2].getFitness())
  {

    return genome1;
  }
  else
  {
    return genome2;
  }

}


/**
*  Selection of genomes in a population.
**/

void Evolution::selection()
{
  std::vector< Genome > selected = std::vector< Genome >();
  std::vector< int > index_selected = std::vector< int >();

  // selects genomes, maintaining population size
  for (int i = 0;
       i < this->params["pop_size"];
       i++)
  {
    int genome = this->tournament(); // selects one genome by tournament

    // makes sure that the same genome wont be selected more than once
    while (std::find(
        index_selected.begin(),
        index_selected.end(),
        genome) != index_selected.end())
    {
      genome = this->tournament();
    }
    selected.push_back(this->population[genome]);
    index_selected.push_back(genome);
  }

  this->cleanMemory(index_selected);

  this->population = selected; // substitutes current population for the selected subset

  // # TEST: Tests if population size remains correct.
  this->tests.testPopsize(
      this->population,
      (int) this->params["pop_size"]);
}

/*
 * Deallocate memory used by the non-selected individuals.
 * */

void Evolution::cleanMemory(std::vector< int > index_selected)
{
  // cleaning memory
  std::cout<<"START cleaning memory for non-selected individuals"<<std::endl;
  for (int i = 0;
       i < this->population.size();
       i++)
  {
    // for non-selected individuals
    if(std::find(
        index_selected.begin(),
        index_selected.end(),
        i) == index_selected.end())
    {
      //std::cout<<"developed genetic-string"<<std::endl;
      auto item = this->population[i].getGeneticString().getStart();
      while (item not_eq NULL)
      {
        auto item2 = item->next;
        delete item;
        item = item2;
      }

//      std::cout<<"grammar genetic-strings"<<std::endl;
//      for( auto &g: this->population[i].getGrammar())
//      {
//        item = g.second.getStart();
//        while (item not_eq NULL)
//        {
//          auto item2 = item->next;
//          delete item;
//          item = item2;
//        }
//      }

      //std::cout<<"decoded genetic-strings"<<std::endl;
      this->cleanVertex(this->population[i].getDgs().getRoot());

      //std::cout<<"scene"<<std::endl;
      QList<QGraphicsItem*> all = this->population[i].getScene()->items();
      for (int i = 0; i < all.size(); i++)
      {
        QGraphicsItem *gi = all[i];
        if(gi->parentItem()==NULL) delete gi;
      }
      delete this->population[i].getScene();


    }
  }
  std::cout<<"FINISH cleaning memory for non-selected individuals"<<std::endl;
}

void Evolution::cleanVertex(DecodedGeneticString::Vertex * v){

  if(v != NULL)
  {
    this->cleanVertex(v->left);
    this->cleanVertex(v->front);
    this->cleanVertex(v->right);
    if(v->item == "C")
      this->cleanVertex(v->back);
    delete v;
  }
}

/**
 *  Saves state of the generations to file.
 */

void
Evolution::exportGenerationMetrics(
    int generation,
    std::vector< int > metrics)
{


  std::ofstream evolution_file;
  std::string path =
      this->path+"experiments/" + this->experiment_name + "/evolution.txt";
  evolution_file.open(
      path,
      std::ofstream::app);

  double maximum_fitness = 0;
  std::string best_genome = "0";
  double average_fitness = 0;

  evolution_file << generation << " ";

  for (int i = 0;
       i < this->getPopulation().size();
       i++)
  {
    if (this->getPopulation()[i].getFitness() >
        maximum_fitness)
    {  // finds the maximum fitness of the population

      best_genome = this->getPopulation()[i].getId();
      maximum_fitness = this->getPopulation()[i].getFitness();
    }
    average_fitness += this->getPopulation()[i].getFitness();  //  sums all fitnesses

  }

  average_fitness /= this->getPopulation().size();  // finds the average of the fitnesses

  evolution_file << best_genome << " "<<
                    maximum_fitness << " " <<
                    average_fitness;
  for (const auto &m : metrics)
  {
    evolution_file << " " << m;
  }
  evolution_file << std::endl;
  evolution_file.close();

}


void Evolution::setupEvolution()
{
  this->aux = Aux(
      this->experiment_name,
      this->getParams(),
      this->path);

  this->readParams();

  // cleans old files and creates folders for the experiment
  aux.removeFolder(this->experiment_name);

  aux.createFolder(this->path+"experiments/"+this->experiment_name);

  // logs parameters configuration
  this->saveParameters();

}


/**
*  Develops genomes of the population: 1- grows genetic-string of the genome according to grammar, 2- decodes it, 3- constructs the phenotype
*  @param argc - default argument
*  @param argv[] - default argument
*  @param LS - Lsystem structure containing the alphabet
*  @param individuals - array with genomes
**/
void Evolution::developIndividuals(
    int argc,
    char *argv[],
    LSystem LS,
    int generation,
    std::vector< Genome > &individuals,
    std::string dir)
{
  // for each genome in the array
  for (size_t i = 0; i < individuals.size(); ++i)
  {
    // develops genome
    individuals[i].developGenomeIndirect(
        argc,
        argv,
        this->params,
        LS,
        generation,
        this->path+"experiments/"+dir);
  }
}

/**
 * Loads population of genomes from files, from previous experiment.
 * @param LS - Lsystem structure containing the alphabet.
 **/
void Evolution::loadPopulation(int generation)
{

  // deletes possible remains of unfinished generation
  std::string pathdir =
      this->path+"experiments/" + this->experiment_name + "/selectedpop" +
      std::to_string(generation + 1);
  system(("exec rm -r " + pathdir).c_str());
  pathdir = this->path+"experiments/" + this->experiment_name + "/offspringpop" +
            std::to_string(generation + 1);
  system(("exec rm -r " + pathdir).c_str());

  // generates list of files (genomes of last population)
  std::system(("ls "+this->path+"experiments/" + this->experiment_name +
               "/selectedpop" + std::to_string(generation) +
               ">"+this->path+"experiments/" + this->experiment_name +
               "/temp.txt").c_str());

  std::ifstream listgenomes((this->path+"experiments/" + this->experiment_name +
                             "/temp.txt").c_str());
  std::string linegenome;

  //for each file (genome)
  while (getline(
      listgenomes,
      linegenome))
  {
    std::vector< std::string > tokens;
    boost::split(
        tokens,
        linegenome,
        boost::is_any_of("_."));

    std::string idgenome = tokens[1];
    std::string idparent1 = tokens[3];
    std::string idparent2 = tokens[5];

    // recreates genome back to population
    // fitness of the parents is not loaded, but it doesnt matter, because at this point it has been saved to the history already
    Genome gen = Genome(
        idgenome,
        idparent1,
        idparent2,
        -1,
        -1);

    // finds number of generation to which the genome belongs to
    int generation_genome = 0;
    int offspring_size =
        this->params["pop_size"] * this->params["offspring_prop"];
    if (this->params["offspring_prop"] == 1)
    {
      generation_genome = (int) trunc(
          std::stoi(idgenome)
          / this->params["pop_size"]) + 1;
    }
    else
    {
      generation_genome =
          (int) trunc((std::stoi(idgenome) - offspring_size)
                      / offspring_size) + 1;
    }

    if (generation_genome == 0)
    { generation_genome = 1; }

    std::ifstream listalphabet(
        this->path+"experiments/" + this->experiment_name + "/offspringpop" +
        std::to_string(generation_genome) + "/genome" + idgenome +
        ".txt");
    std::string linealphabet;
    // for each letter of the alphabet
    while (getline(
        listalphabet,
        linealphabet))
    {

      // gets letter and production rule from file
      std::vector< std::string > items;
      boost::split(
          items,
          linealphabet,
          boost::is_any_of(" "));
      std::vector< std::string > items_rule(
          items.begin() + 1,
          items.begin() + items.size() -
          1);

      // build a genetic-string with the production rule for the letter
      auto lgs =  GeneticString();
      lgs = gen.build_genetic_string(
          lgs,
          items_rule);

      // adds letter and its production rule (made a genetic-string) to the grammar of the genome
      gen.addLetterGrammar(
          items[0],
          lgs);

    }

    std::ifstream listmeasures(
        this->path+"experiments/" + this->experiment_name + "/offspringpop" +
        std::to_string(generation_genome) + "/measures" + idgenome +
        ".txt");
    std::string linemeasures;
    // for each measure of the list
    while (getline(
        listmeasures,
        linemeasures))
    {

      std::vector< std::string > tokens;
      boost::split(
          tokens,
          linemeasures,
          boost::is_any_of(":"));

      gen.updateMeasure(
          tokens[0],
          std::stod(tokens[1]));
    }

    this->population.push_back(gen);  // adds genome to the population

  }


};

/**
 *  Loads state of previous experiment.
 **/
int Evolution::loadExperiment()
{

  this->logsTime("start");

  // loads state of parameters from previous experiment
  this->loadsParams();
  this->aux = Aux(
      this->experiment_name,
      this->getParams(),
      this->path);


  // loads generation number from previous  experiment
  int gi = std::stoi(this->readsEvolutionState()[0]);
  // loads next_id from previous experiment
  this->next_id = std::stoi(this->readsEvolutionState()[1]);

  // loads state of the morphological_grid_accumulated
  std::string line;
  std::ifstream myfile(
      this->path+"experiments/" + this->experiment_name +
      "/morphological_grid_accumulated.txt");
  while (getline(
      myfile,
      line))
  {
    std::vector< std::string > tokens, tokens2;
    boost::split(
        tokens,
        line,
        boost::is_any_of("-"));
    boost::split(
        tokens2,
        tokens[1],
        boost::is_any_of(" "));
    std::vector< std::string > tokens3(
        tokens2.begin(),
        tokens2.begin() + tokens2.size() -
        1);
    std::vector< std::string > points;
    for (int i = 0;
         i < tokens3.size();
         i++)
    {
      points.push_back(tokens3[i]);
    }
    this->morphological_grid_accumulated[tokens[0]] = points;
  }
  myfile.close();


  // loads experiment  population
  this->loadPopulation(gi);

  return gi;

}

/* Saves the fitness for genome after simulation.
 * */
void Evolution::saveFitness(
    int genome_index,
    double fitness)
{
  this->offspring[genome_index].updateFitness(fitness);
}



/**
*  Evolution in the search for locomotion - part 1 of the process.
**/
double Evolution::runExperiment_part1(
    int generation)
{
  int argc = 1;
  char *argv[] = { "a"};

  // loads alphabet with letters and commands
  LSystem LS;

  this->aux.logs("------------ generation " + std::to_string(generation) + " ------------");
  this->logsTime("start");

  this->aux.createFolder(this->path+"experiments/"+this->experiment_name +
                       "/offspringpop"+std::to_string(generation));

  if(generation == 1)
  {
    this->createHeader();

    // initializes population
    this->initPopulation(LS);
  }
  else{

    this->aux.createFolder(
        this->path+"experiments/"+this->experiment_name + "/selectedpop" + std::to_string(generation));

    this->offspring = std::vector< Genome >();

    // creates offspring
    this->crossover(LS);

    // mutates new individuals
    this->mutation(LS);

  }

  // develops genomes of the initial population
  this->developIndividuals(
        argc,
        argv,
        LS,
        generation,
        this->offspring,
        this->experiment_name + "/offspringpop");

  // measures phenotypes of the individuals
  this->measureIndividuals(
        generation,
        this->offspring,
        "/offspringpop");

  // updates the average measures for the population
  this->savesValidity(generation);

}


/**
*  Evolution in the search for locomotion - part 2 of the process. After
 *  fitness evaluation.
**/
double Evolution::runExperiment_part2(int generation)
{
    // saves a history with the fitnesses of new individuals
    this->saveHistory(generation);

      // adds new individuals to population
    for (int j = 0;
         j < this->offspring.size();
         j++)
    {
      this->population.push_back(this->offspring[j]);
    }

    std::vector< int > niche_measures = this->calculateNicheCoverage();

    if(generation != 1)
    {
      // selects individuals, keeping the population with a fixed size
      this->selection();

      // saves phenotypes of the selected population to a separated folder (only for visualization issues)
      this->exportPop(generation);
    }

    // saves metrics of evolution to file
    this->exportGenerationMetrics(
        generation,
        niche_measures);

    // saves the number of the last generation created/evaluated
    this->writesEvolutionState(
        generation,
        this->next_id);

    this->summaryNicheCoverage();


    this->logsTime("end");
}

void Evolution::summaryNicheCoverage()
{

  std::ofstream file;

  std::string path = this->path+"experiments/" + this->experiment_name +
                     "/morphological_grid_summary.txt";
  file.open(path);

  for (const auto &it : this->morphological_grid_accumulated)
  {

    file << it.first + " " << it.second.size() << std::endl;
  }
  file.close();

}


std::vector< Genome > Evolution::getPopulation()
{
  return this->population;
}

std::map< std::string, double > Evolution::getParams()
{
  return params;
}


/*
 * Exports the parameters of the experiment.
 **/
void Evolution::saveParameters()
{

  std::ofstream param_file;
  std::string path =
      this->path+"experiments/" + this->experiment_name + "/configuration.txt";
  param_file.open(path);

  // writes each parameter to a different line in a the file
  for (auto &it : this->getParams())
  {

    param_file << it.first << " " << it.second;
    param_file << std::endl;
  }
  param_file.close();
}


/*
 * Logs time.
 **/
void Evolution::logsTime(std::string moment)
{

  time_t sta = time(0);
  char *dtsta = ctime(&sta);
  this->aux.logs("experiment " + moment + ": " + dtsta);

}

/*
 * Logs the generation from which the recovered evolution should start from.
 * */
void Evolution::writesEvolutionState(
    int generation,
    int next_id)
{

  std::ofstream logs_file;
  std::string path = this->path+"experiments/" + this->experiment_name +
                     "/evolutionstate.txt";
  logs_file.open(path);
  logs_file << generation << " " << next_id;
  logs_file.close();
}

/*
 * Reads number of the generation from which the recovered evolution should start from.
 * */
std::vector< std::string > Evolution::readsEvolutionState()
{

  std::string line;
  std::ifstream myfile(
      this->path+"experiments/" + this->experiment_name +
      "/evolutionstate.txt");
  if (myfile.is_open())
  {

    getline(
        myfile,
        line);
    std::vector< std::string > tokens;
    // parameters label and value separated by space
    boost::split(
        tokens,
        line,
        boost::is_any_of(" "));

    return tokens;

  }
  else
  {
    this->aux.logs("Unable to open evolutionstate file.");
  }
  myfile.close();

}


/*
 * Calculates the quality metric for the novelty search: niche coverage and bins of measures
 * */
std::vector< int >
Evolution::calculateNicheCoverage()
{
  for (int i = 0;
       i < this->offspring.size();
       i++)
  {

    std::string key_point = "";
    double distance = 0;

    // for each measure (dimension)
    for (const auto &it : this->offspring[i].getMeasures())
    {
      // accounts for NC
      // for each bin
      for (int b = 1;
           b <= this->params["grid_bins"];
           b++)
      {
        // if value is zero, sets into the first bin
        if (it.second == 0 and b == 1)
        {

          key_point += std::to_string(b) + "|";
          distance +=
              -1 * (it.second - (b / this->params["grid_bins"]));
        }
        // otherwise, sets value for measure into the correct bin
        if (it.second > ((b - 1) / this->params["grid_bins"]) and
            it.second <= (b / this->params["grid_bins"]))
        {
          key_point += std::to_string(b) + "|";
          distance +=
              -1 * (it.second - (b / this->params["grid_bins"]));
        }
      }

      // accounts for occurrence of value for measure
      if (std::find(
          this->morphological_measures_accumulated[it.first].begin(),
          this->morphological_measures_accumulated[it.first].end(),
          it.second)
          == this->morphological_measures_accumulated[it.first].end())
      {

        this->morphological_measures_accumulated[it.first].push_back(it.second);
      }
    }

    // if point already exists in the array, adds an individual and the id
    if (this->morphological_grid_accumulated.count(key_point) > 0)
    {

      this->morphological_grid_accumulated[key_point].push_back
          (this->offspring[i].getId());

      // if point does not exist in the array yet, , adds new point with its first individual and the difference between them
    }
    else
    {
      std::vector< std::string > individual;
      individual.push_back(this->offspring[i].getId());
      this->morphological_grid_accumulated[key_point] = individual;
    }
  }



  // logs state of the grid
  std::ofstream myfile;
  std::string path = this->path+"experiments/" + this->experiment_name +
                     "/morphological_grid_accumulated.txt";
  myfile.open(path);
  for (const auto &it : this->morphological_grid_accumulated)
  {
    myfile << it.first + "-";
    for (int i = 0;
         i < it.second.size();
         i++)
    {
      myfile << it.second[i] << " ";
    }
    myfile << std::endl;
  }
  myfile.close();


  std::vector< int > morphological_grids;
  morphological_grids.push_back((int) this->morphological_grid_generation.size());  // pos 0
  morphological_grids.push_back((int) this->morphological_grid_accumulated.size()); // pos 1
  for (const auto &it :this->morphological_measures_accumulated)
  {
    morphological_grids.push_back((int) it.second.size());
  }

  return morphological_grids;
}

