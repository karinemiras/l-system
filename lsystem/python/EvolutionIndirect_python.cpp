
#include <boost/python.hpp>
#include "EvolutionIndirect.h"

class EvolutionIndirect_python : public EvolutionIndirect {
public:
  EvolutionIndirect_python(std::string experiment_name,
                      std::string type_experiment)
     : EvolutionIndirect(experiment_name, type_experiment)
  {}


//  double runExperiment()
//  {
//    EvolutionIndirect::runExperiment(0, nullptr);
//  }

};

BOOST_PYTHON_MODULE (lsystem_python)
{
  boost::python::class_< EvolutionIndirect_python, boost::noncopyable >(
          "EvolutionIndirect",
          boost::python::init< std::string,
                               std::string >())
      .def("setupEvolution",
           &EvolutionIndirect_python::setupEvolution)
      .def("runExperiment_evolve1",
           &EvolutionIndirect_python::runExperiment_evolve1)
      .def("runExperiment_evolve2",
           &EvolutionIndirect_python::runExperiment_evolve2)
      .def("saveLocomotionFitness",
           &EvolutionIndirect_python::saveLocomotionFitness)
      .def("runExperiment_learn1",
           &EvolutionIndirect_python::runExperiment_learn1)
      .def("runExperiment_learn2",
           &EvolutionIndirect_python::runExperiment_learn2)
      .def("writesEvolutionState",
           &EvolutionIndirect_python::writesEvolutionState)
      ;
}
