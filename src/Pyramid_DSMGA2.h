// Brian Goldman

// Implemention of the Parameter-less Population Pyramid_DSMGA2
// Full description given in our publication:
// "Fast and Efficient Black Box Optimization using the Parameter-less Population Pyramid_DSMGA2"
// by B. W. Goldman and W. F. Punch

#ifndef PYRAMID_DSMGA2_H_
#define PYRAMID_DSMGA2_H_

#include "Population_DSMGA2.h"
#include "Evaluation.h"
#include "Util.h"
#include "HillClimb.h"
#include "Configuration.h"
#include "Optimizer.h"
#include "MiddleLayer.h"

// Implements the Optimizer interface
class Pyramid_DSMGA2 : public Optimizer {
 public:
  Pyramid_DSMGA2(Random& _rand, shared_ptr<Evaluator> _evaluator,
          Configuration& _config)
      : Optimizer(_rand, _evaluator, _config),
        only_add_improvements(_config.get<int>("only_add_improvements")),
        hill_climber(_config.get<hill_climb::pointer>("hill_climber")),
        local_counter(new Middle_Layer(config, _evaluator, false)),
        cross_counter(new Middle_Layer(config, _evaluator, false)),
        restarts(0) {
  }
  // Iteratively improves the solution using the pyramid_dsmga2 of populations
  // leverages the Population class extensively
  void climb(Chromosome & chromosome, float & fitness);
  // Peforms one complete iteration of
  // * random generation
  // * hill climbing
  // * crossover with each level of the pyramid_dsmga2 (climb function)
  bool iterate() override;
  // Flatten metadata recorded about search into a string.
  string finalize() override;
  create_optimizer(Pyramid_DSMGA2);

 private:
  // Handles adding a solution to the population.
  // Returns true if it was unique and therefore added.
  bool add_unique(const Chromosome & chromosome, size_t level);

  // the pyramid_dsmga2 of populations
  vector<Population_DSMGA2> pops;
  // keeps track of the set of solutions in the pyramid_dsmga2
  std::unordered_set<vector<bool>> seen;
  // configuration options
  bool only_add_improvements;
  hill_climb::pointer hill_climber;
  // Used for recording purposes
  shared_ptr<Evaluator> local_counter, cross_counter;
  size_t restarts;

};

#endif /* PYRAMID_DSMGA2_H_ */
