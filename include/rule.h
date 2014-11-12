#ifndef RULE_H
#define RULE_H

#include <vector>
#include "generator.h"

/**
 * Rule is basically a vector of Generator objects and is responsible
 * for providing access to them and for deleting them.
 */
class Rule
{
public:
  // Constructor
  Rule()
  {
  }

  // Destructor
  virtual ~Rule()
  {
  }

  // Print the number of generators currently stored by the rule
  unsigned n_generators() const
  {
    return generators.size();
  }

  // Push a new Generator pointer onto the back of the vector.  The
  // Rule takes ownership of the pointer and is responsible for
  // deleting it.
  void push_back(const Generator & generator)
  {
    generators.push_back(generator);
  }

  // Writeable reference to the ith generator, if it exists
  Generator& operator[](unsigned i)
  {
    if (i >= generators.size())
      {
        std::cerr << "Invalid index " << i << " requested!" << std::endl;
        std::abort();
      }

    return generators[i];
  }

  void generate_points_and_weights(std::vector<Point<mpfr_class> > & generated_points,
                                   std::vector<mpfr_class> & generated_weights) const;

protected:
  std::vector<Generator> generators;
};

#endif
