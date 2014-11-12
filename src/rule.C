#include "rule.h"

unsigned Rule::n_generators() const
{
  return generators.size();
}



void Rule::push_back(const Generator & generator)
{
  generators.push_back(generator);
}



// Writeable reference to the ith generator, if it exists
Generator& Rule::operator[](unsigned i)
{
  if (i >= generators.size())
    {
      std::cerr << "Invalid index " << i << " requested!" << std::endl;
      std::abort();
    }

  return generators[i];
}



void Rule::generate_points_and_weights(std::vector<Point<mpfr_class> > & generated_points,
                                       std::vector<mpfr_class> & generated_weights) const
{
  // Vectors for individual generators
  std::vector<Point<mpfr_class> > single_generator_points;
  std::vector<mpfr_class> single_generator_weights;

  for (unsigned i=0; i<generators.size(); ++i)
    {
      generators[i].generate_points_and_weights(single_generator_points,
                                                single_generator_weights);

      // Append the points and weights to the end of the user's vectors
      generated_points.insert(/*pos=*/generated_points.end(),
                              /*first=*/single_generator_points.begin(),
                              /*last=*/single_generator_points.end());

      generated_weights.insert(/*pos=*/generated_weights.end(),
                               /*first=*/single_generator_weights.begin(),
                               /*last=*/single_generator_weights.end());
    }
}
