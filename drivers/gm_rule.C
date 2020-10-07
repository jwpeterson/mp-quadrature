// mp-quadrature includes
#include "grundmann_moller.h"
#include "exact.h"

// C++ includes
#include <iomanip>
#include <algorithm> // std::sort
#include <set>
#include <cstdlib> // std::abort
#include <getopt.h> // getopt_long()
#include <cmath> // log10

// An unoptimized pow function for mpq_class objects
mpq_class pow(const mpq_class & base, unsigned power)
{
  mpq_class result(1);
  for (unsigned i=0; i<power; ++i)
    result *= base;
  return result;
}

// Formats a rational number with decimal points and long double
// string literals for generating C++ code.
std::string format_rational(const mpq_class & number)
{
  std::string number_string = number.get_str();

  // Find the forward slash, if there is one
  size_t slash_pos = number_string.find('/');

  if (slash_pos != std::string::npos)
    {
      // Insert ".L" before the "/"
      number_string.insert(slash_pos, ".L");

      // And insert a ".L" for the denominator
      number_string += ".L";
    }

  return number_string;
}

// Function explaining command line options that gets called when the
// program is run with --help.
void usage();

// As seen on cplusplus.com, how you specify a random number generator for std::random_shuffle:
int myrandom (int i) { return random()%i; }

// Function which verifies a given quadrature rule defined by the
// points x and weights w.  The third argument is a vector defining
// the order in which to evaluate the rule.  The last two arguments
// are the maximum and average quadrature error, respectively,
// computed by the routine.
template <class T>
void verify(unsigned rule_index,
            const std::vector<Point<T> >& x,
            const std::vector<T>& w,
            const std::vector<size_t>& indices,
            T & max_error,
            T & avg_error);

// Function which prints a quadrature rule in the order implied by
// indices.  Currently only works for printing rational rules since it
// is hard-coded to use print_rational.
template <class T>
void print_ordered_rule(const std::vector<Point<T> >& x,
                        const std::vector<T>& w,
                        const std::vector<size_t>& indices);

// Note: these results are from a Mac.  The random results vary between OSX and Linux due to different orderings in random_shuffle.
// Some results for b=53
// Max error
//  s  less-than               unsorted                 abs-less-than            abs-greater-than         random                  less-than-interleaved
//  0  0                       0                        0                        0                        0                       0
//  1  6.9388939039072284e-18  3.4694469519536142e-18   3.4694469519536142e-18   6.9388939039072284e-18   2.7755575615628914e-17  6.9388939039072284e-18
//  2  5.5511151231257827e-17  3.4694469519536142e-17   4.8572257327350599e-17   8.3266726846886741e-17   2.7755575615628914e-17  1.3877787807814457e-17
//  3  1.3877787807814457e-16  5.2041704279304213e-17   3.0531133177191805e-16   1.6653345369377348e-16   1.3877787807814457e-16  1.6653345369377348e-16
//  4  1.4155343563970746e-15  1.0755285551056204e-16   2.2898349882893854e-16   1.1102230246251565e-16   3.6082248300317588e-16  1.0824674490095276e-15
//  5  3.6359804056473877e-15  2.2759572004815709e-15   6.8278716014447127e-15   6.8278716014447127e-15   1.8041124150158794e-16  2.7478019859472624e-15
//  6  7.7229889150487452e-15  1.1657341758564144e-14   9.7422070410857486e-15   1.0463852007092100e-14   6.3837823915946501e-16  3.4972025275692431e-15
//  7  4.8544501751734970e-14  2.1593837828959295e-14   4.7378767575878555e-14   5.0626169922907138e-14   5.4123372450476381e-15  9.7422070410857486e-15
//  8  2.7641777755604835e-13  3.0531133177191805e-14   7.0360384185619296e-14   6.8722805224297190e-14   2.4494295480792516e-15  1.6792123247455493e-14
//  9  8.1121220851798626e-13  1.4993561947562739e-13   1.5404344466674047e-14   1.8374191057546341e-14   8.5764728652293343e-15  8.9026008787129740e-15
// 10  2.8433644327918728e-12  3.5063618675224006e-13   4.5594084063793616e-13   4.6787573815265660e-13   2.6839641620313159e-14  3.3287261835823756e-13
//
// Avg error
// unsorted case: for ((i=0; i < 11 ; i++)) do ./drivers/gm_rule --rule-index $i --binary-digits 53 --unsorted | grep avg_error; done
// random case:   for ((i=0; i < 11 ; i++)) do ./drivers/gm_rule --rule-index $i --binary-digits 53 --random-shuffle | grep avg_error; done
// other cases:   for ((i=0; i < 11 ; i++)) do ./drivers/gm_rule --rule-index $i --binary-digits 53 | grep avg_error; done
//  s  less-than               unsorted                abs-less-than           abs-greater-than        random                  less-than-interleaved
//  0  0                       0                       0                       0                       0                       0
//  1  1.6581915579190069e-18  1.5689043201849064e-18  1.5689043201849064e-18  1.6581915579190069e-18  2.6658618123467111e-18  1.6837021972716069e-18
//  2  3.0738665838321340e-18  3.1380493067787995e-18  4.1138312620157538e-18  3.5430908731035727e-18  1.2964398619107329e-18  1.0504487086251067e-18
//  3  3.9136903941960559e-18  2.6639330983104543e-18  4.6401261206352952e-18  4.0242028552203024e-18  2.9357510388028471e-18  3.3429566984969719e-18
//  4  9.2620496069905424e-18  3.5192642969695190e-18  5.1108528666820796e-18  3.3137480489198983e-18  3.3722247932626442e-18  7.6460280889195532e-18
//  5  3.0733070788485097e-17  1.4069683811255716e-17  2.7999212009269425e-17  3.2453510115530589e-17  3.8114386857813416e-18  1.2681975992947120e-17
//  6  6.5749721438699001e-17  3.0136129406321188e-17  2.8582055703474798e-17  3.1138949086784063e-17  6.0888840309000377e-18  1.6749362669144901e-17
//  7  1.1578409343445259e-16  5.1221328895576901e-17  7.2239034144785803e-17  8.0058527055389632e-17  1.2669779380806972e-17  1.9618676330643032e-17
//  8  3.7268258211632601e-16  8.5008377664200395e-17  1.0757563952489258e-16  8.5141731208284614e-17  1.3177982821287221e-17  3.8288646032144708e-17
//  9  8.5413981279550034e-16  1.7164814304737281e-16  7.6348720727708895e-17  7.7775636199856525e-17  1.7352886783461683e-17  4.1896003236230336e-17
// 10  2.0768732594554229e-15  3.0515096041501392e-16  3.8272795924753273e-16  3.9061139591631006e-16  3.5793857834546364e-17  2.9624320337843564e-16
//
// The results are counter-intuitive with what I expected... the
// unsorted points and weights, as produced by the algorithm, have
// better round-off error characteristics than the sorted points and
// weights.
//
// The next thing I tried was sorting according to absolute value.  In
// some cases it was slightly better than the other sorted option, but
// in most cases it was even worse.
//
// The abs-greater-than sorting was not noticeably worse or better
// than the all the other methods.
//
// Randomly shuffling the points and weights seems to work
// surprisingly well as far as round-off error is concerned... it
// seems to produce the lowest of all the methods tried.  So this at
// least proves that there are some orderings that are really better
// than others.
int main(int argc, char** argv)
{
  std::cout.precision(32);
  std::cout.setf(std::ios_base::scientific);

  // By default, we sort the points and weights according to the
  // comparator used by the indirect_sort algorithm.  If false, print
  // the points and weights in the order they are computed, if true
  // print them in the weight-sorted order.  Eventually make this
  // selectable on the command line.
  bool do_sort = true;

  // If true, randomly shuffles the indices rather than using any
  // sorting algorithm (--unsorted is ignored).
  bool do_random_shuffle = false;
  unsigned n_random_trials = 1;
  unsigned random_seed = 0;

  // If true, use std::next_permutation to test a sequence of
  // point/weight permutations for the desired rule.
  bool do_permutations = false;

  // If true, attempt to find the best permutation of points and
  // weights by using an optimization algorithm.
  bool do_optimize = false;
  unsigned n_optimization_trials = 10;

  // If false, will not compute all the verification integrals.  This
  // can be nice for high-order rules which take an extremely long
  // time to verify.  Eventually make this selectable on the command
  // line.
  bool do_verification = true;

  // Possibly turn on interleaving, which might be useful with certain
  // sorting algorithms.
  bool do_interleave = false;

  // Read rule index, s, from command line.  Any integer s>=0 is valid.
  unsigned s = 0;

  // # of binary digits
  // 53 binary digits is about what you normally get with a double.
  unsigned b = 256;

  // options descriptor - this associates several "long" and "short"
  // command line options.  The last element of the longopts array has
  // to be filled with zeros.
  static struct option longopts[] =
    {
      {"no-verify",       no_argument,       NULL, 'x'},
      {"rule-index",      required_argument, NULL, 's'},
      {"unsorted",        no_argument,       NULL, 'u'},
      {"binary-digits",   required_argument, NULL, 'b'},
      {"n-random-trials", optional_argument, NULL, 'r'},
      {"random-seed",     required_argument, NULL, 'e'},
      {"interleave",      no_argument,       NULL, 'i'},
      {"optimize",        optional_argument, NULL, 'o'},
      {"permute",         no_argument,       NULL, 'p'},
      {"help",            no_argument,       NULL, 'h'},
      { NULL,             0,                 NULL,  0 }
    };

  // Parse command line options using getopt_long()
  int ch = -1;
  while ((ch = getopt_long(argc, argv, "xs:ub:r:e:io:ph", longopts, NULL)) != -1)
    {
      switch (ch)
        {
        case 'x':
          do_verification = false;
          break;

        case 's':
          s = atoi(optarg);
          break;

        case 'u':
          do_sort = false;
          break;

        case 'b':
          b = atoi(optarg);
          break;

        case 'r':
          {
          do_random_shuffle = true;
          if (optarg != NULL)
            n_random_trials = atoi(optarg);
          break;
          }

        case 'e':
          {
            if (optarg != NULL)
              random_seed = atoi(optarg);
            break;
          }

        case 'i':
          do_interleave = true;
          break;

        case 'o':
          do_optimize = true;
          if (optarg != NULL)
            n_optimization_trials = atoi(optarg);
          break;

        case 'p':
          do_permutations = true;
          break;

        case 'h':
          usage();
          return 0;

        default:
          usage();
        }
    } // end while

  std::cout << "Using random seed " << random_seed << std::endl;
  srandom(random_seed);

  if (do_random_shuffle)
    std::cout << "Performing " << n_random_trials << " random trials." << std::endl;

  // Set the default precision for real-valued numbers
  // std::cout << "Using " << b << " binary digits for mpfr_class objects." << std::endl;
  mpfr_set_default_prec(b);

  GrundmannMoller gm;
  gm.rule(s);

  // Get a reference to the points and weights arrays
  const std::vector<Point<mpq_class> >& x = gm.get_points();
  const std::vector<mpq_class>& w = gm.get_weights();

  std::cout << "Index " << s
            << " rule has degree " << 2*s+1
            << ", and " << x.size() << " point(s)." << std::endl;

  // Initialize the indices array.  This will be used to access into
  // the x and w arrays from now on...
  std::vector<size_t> indices(w.size());
  iota(indices.begin(), indices.end(), 0);

  // Execute different code depending on what the user requested
  if (do_random_shuffle)
    {
      // Convert the rational point/weight values to real values.
      std::vector<mpfr_class> w_real(w.begin(), w.end());
      std::vector<Point<mpfr_class> > x_real(x.size());
      for (unsigned i=0; i<x_real.size(); ++i)
        for (unsigned j=0; j<3; ++j)
          x_real[i](j) = x[i](j);

      mpfr_class
        max_error = 0.,
        avg_error = 0.;

      // Use 2x higher precision than the inputs to store the min error
      mpfr_class min_max_error(/*value=*/1, /*precision=*/2*b);
      mpfr_class max_max_error(/*value=*/0, /*precision=*/2*b);

      // Also keep track of the average error (maybe it has more variation than the max?)
      mpfr_class min_avg_error(/*value=*/1, /*precision=*/2*b);
      mpfr_class max_avg_error(/*value=*/0, /*precision=*/2*b);

      std::vector<size_t> best_indices;

      for (unsigned trial=0; trial<n_random_trials; ++trial)
        {
          // Print a status update every so many trials
          if (trial % 1000000 == 0)
            std::cout << "Completed trial " << trial << std::endl;

          // Choose a random order
          std::random_shuffle(indices.begin(), indices.end(), myrandom);

          // Do verification, collect error results
          verify(s, x_real, w_real, indices, max_error, avg_error);

          // Debugging:
          // std::cout << "max_error = " << max_error << std::endl;
          // std::cout << "avg_error = " << avg_error << std::endl;

          // Determination of a "best" rule based on the smallest average error found
          if (avg_error < min_avg_error)
            {
              best_indices = indices;

              std::cout << "New best quadrature error found, avg_error = " << avg_error << std::endl;
              print_ordered_rule(x, w, best_indices);

              // Print the indices
              for (unsigned i=0; i<indices.size(); ++i)
                std::cout << indices[i] << " ";
              std::cout << std::endl;
            }

          // Store the min max_error
          if (max_error < min_max_error)
            min_max_error = max_error;

          // Store the max max_error
          if (max_error > max_max_error)
            max_max_error = max_error;

          // Store the max avg_error
          if (avg_error > max_avg_error)
            max_avg_error = avg_error;

          // Store the min avg_error
          if (avg_error < min_avg_error)
            min_avg_error = avg_error;
        }

      // Report best result found at the end
      std::cout << "Smallest max error: " << min_max_error << std::endl;
      std::cout << "Largest max error: " << max_max_error << std::endl;
      std::cout << "Smallest avg error: " << min_avg_error << std::endl;
      std::cout << "Largest avg error: " << max_avg_error << std::endl;
      print_ordered_rule(x, w, best_indices);
      // Print the best_indices
      for (unsigned i=0; i<best_indices.size(); ++i)
        std::cout << best_indices[i] << " ";
      std::cout << std::endl;

      // We are done, so return from main()
      return 0;
    }

  if (do_permutations)
    {
      // Try different permutations of adding up the real-valued weights.
      // Can we pick an ordering based on the lowest error in the sum of
      // weights observed?  Note: we are aware of the different ways of
      // summing sets of numbers to reduce round-off error.  This study
      // only checks different orderings of the summation and their error.

      // Build a set of unique weight values...
      std::set<mpq_class> unique_wts(w.begin(), w.end());

      // Get counts for each unique weight
      std::vector<unsigned> counts;
      std::set<mpq_class>::iterator
        it = unique_wts.begin(),
        end = unique_wts.end();

      std::cout << "unique weights = ";
      for (; it != end; ++it)
        {
          std::cout << format_rational(*it) << " ";
          counts.push_back(std::count(w.begin(), w.end(), *it));
        }
      std::cout << std::endl;

      std::cout << "weight counts = ";
      for (unsigned i=0; i<counts.size(); ++i)
        std::cout << counts[i] << " ";
      std::cout << std::endl;


      // The permutations of the weight vector will be unique (i.e. not
      // contain repeated permutations which are identical!)
      std::vector<mpq_class> weight_copy(w);

      // Start making permutations from an initially sorted vector
      std::sort(weight_copy.begin(), weight_copy.end());

      mpq_class true_solution_rational("1/6");
      mpfr_class true_solution(true_solution_rational, /*precision=*/2*b);

      mpfr_class min_error = 1.0;
      std::vector<mpq_class> best_weights_permutation(weight_copy);
      while (std::next_permutation(weight_copy.begin(), weight_copy.end()))
        {
          // for (unsigned i=0; i<weight_copy.size(); ++i)
          //   std::cout << weight_copy[i] << " ";
          // std::cout << std::endl;

          // Compute the (real-valued) sum
          mpfr_class sum = 0.0;
          for (unsigned i=0; i<weight_copy.size(); ++i)
            sum += mpfr_class(weight_copy[i]);

          // Compute the error
          mpfr_class current_error(/*value=*/my_abs(mpfr_class(sum - true_solution)), /*precision=*/2*b);

          // Possibly update the best order
          if (current_error < min_error)
            {
              std::cout << "New minimum found=" << current_error << std::endl;
              min_error = current_error;
              best_weights_permutation = weight_copy;

              // Print the new best ordering
              for (unsigned i=0; i<best_weights_permutation.size(); ++i)
                std::cout << best_weights_permutation[i] << " ";
              std::cout << std::endl;

              // Note: we haven't actually reordered the rule, so it doesn't make sense to print it...
              // print_ordered_rule(x, w, best_indices);
            }

          // We apparently need to handle the case where current_error ==
          // 0.0.  I didn't think this would happen but apparently it
          // can...  There's no point in further searching if this happens!
          if (current_error == 0.0)
            break;
        }

      // Report the best ordering and the smallest error
      std::cout << "Smallest error = " << min_error << std::endl;
      std::cout << "Best ordering = " << std::endl;
      for (unsigned i=0; i<best_weights_permutation.size(); ++i)
        std::cout << best_weights_permutation[i] << " ";
      std::cout << std::endl;

      // We are done with permutations, so exit from main
      return 0;
    }

  if (do_optimize)
    {
      // Convert the rational point/weight values to real values.
      std::vector<mpfr_class> w_real(w.begin(), w.end());
      std::vector<Point<mpfr_class> > x_real(x.size());
      for (unsigned i=0; i<x_real.size(); ++i)
        for (unsigned j=0; j<3; ++j)
          x_real[i](j) = x[i](j);

      // Keep track of the best max error and best average error found.
      mpfr_class
        best_max_error(1., /*precision=*/2*b),
        best_avg_error(1., /*precision=*/2*b);

      // Keep track of the best ordering
      std::vector<size_t> best_indices = indices;

      for (unsigned restart=0; restart<n_optimization_trials; ++restart)
        {
          // Start the optimization routine from a random configuration
          std::random_shuffle(indices.begin(), indices.end(), myrandom);

          // std::cout << "Initial configuration:" << std::endl;
          // print_ordered_rule(x, w, indices);

          // Have the user hit a key after inspecting the initial configuration
          // std::cout << "Press enter to continue..." << std::endl;
          // std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

          // For each entry in the vector, swap it with every other entry
          // and compute the error each time, then pick the swap that most
          // reduces the error, if one exists.
          //
          // I typically see some swaps during the second sweep, but not
          // as many as during the first one.  I'm undecided on whether
          // multiple sweeps are worth doing...
          for (unsigned sweep=0; sweep<4; ++sweep)
            for (unsigned outer=0; outer<x_real.size(); ++outer)
              {
                std::vector<mpfr_class>
                  max_configuration_errors(x_real.size()),
                  avg_configuration_errors(x_real.size());

                for (unsigned inner=0; inner<x_real.size(); ++inner)
                  {
                    // When outer==inner, this is a no-op... but that's OK,
                    // we want to compute the error in the unswapped
                    // configuration too...
                    std::swap(indices[outer], indices[inner]);

                    // Compute and store the error for this set of indices
                    verify(s, x_real, w_real, indices, max_configuration_errors[inner], avg_configuration_errors[inner]);

                    // Swap back so we can compute the next error
                    std::swap(indices[outer], indices[inner]);
                  }

                // Report the error for all configurations
                // std::cout << "Max errors" << std::endl;
                // for (unsigned i=0; i<x_real.size(); ++i)
                //   std::cout << max_configuration_errors[i] << std::endl;

                // std::cout << "Avg errors" << std::endl;
                // for (unsigned i=0; i<x_real.size(); ++i)
                //   std::cout << avg_configuration_errors[i] << std::endl;

                // Get an iterator to the min avg error entry
                std::vector<mpfr_class>::iterator it =
                  std::min_element(avg_configuration_errors.begin(),
                                   avg_configuration_errors.end());

                // Get an index from this iterator
                unsigned swap_index = std::distance(avg_configuration_errors.begin(), it);
                // std::cout << "Swapping indices " << outer << " and " << swap_index << std::endl;

                // Make the swap with this entry permanent
                std::swap(indices[outer], indices[swap_index]);

                // Print the re-ordered rule at each iteration
                // std::cout << "Re-ordered rule." << std::endl;
                // print_ordered_rule(x, w, indices);
              }

          // Print the final for this iteration of the algorithm.
          mpfr_class
            max_error = 0.,
            avg_error = 0.;

          verify(s, x_real, w_real, indices, max_error, avg_error);

          // Keep track of the best_indices and print the result if it is a new best
          if (avg_error < best_avg_error)
            {
              best_avg_error = avg_error;
              best_max_error = max_error;

              best_indices = indices;

              std::cout << "\nNew optimal rule found: " << std::endl;
              std::cout << "max_error = " << max_error << std::endl;
              std::cout << "avg_error = " << avg_error << std::endl;

              print_ordered_rule(x, w, indices);

              // Print the indices
              for (unsigned i=0; i<indices.size(); ++i)
                std::cout << indices[i] << " ";
              std::cout << std::endl;
            }
        } // end for (restart)


      // Print the best overall result that we found
      std::cout << "\nBest overall rule: " << std::endl;
      std::cout << "best_max_error = " << best_max_error << std::endl;
      std::cout << "best_avg_error = " << best_avg_error << std::endl;

      print_ordered_rule(x, w, best_indices);

      // Print the best_indices
      for (unsigned i=0; i<best_indices.size(); ++i)
        std::cout << best_indices[i] << " ";
      std::cout << std::endl;

      // We are done with the optimization algorithm, so exit from main
      return 0;
    }

  // The next options are not mutually exclusive, for example, we can
  // sort and interleave, but not do verification, etc.
  if (do_sort)
    {
      // If requested, indirect sort the points and weights vectors
      indirect_sort(w.begin(), w.end(), indices);

      // Debugging: Print the sorted rule
      // std::cout << "Sorted rule: " << std::endl;
      // print_ordered_rule(x, w, indices);
    }

  if (do_interleave)
    {
      // Interleave the sorted array. This is potentially useful with
      // the less-than ordering, so that the largest negative number
      // is added to the largest positive number, the next-to-largest
      // negative number is added to the next-to-largest positive
      // number, etc.
      interleave(indices);
    }

  // Debugging: print the interleaved rule
  std::cout << "Rule order based on indices: " << std::endl;
  print_ordered_rule(x, w, indices);

  // Convert the rational point/weight values to finite precision
  // real-valued numbers with b binary digits.  The real-valued points
  // and weights can be used to simulate round-off error in a normal
  // C++ code using doubles.
  std::vector<mpfr_class> w_real(w.begin(), w.end());
  std::vector<Point<mpfr_class> > x_real(x.size());
  for (unsigned i=0; i<x_real.size(); ++i)
    for (unsigned j=0; j<3; ++j)
      x_real[i](j) = x[i](j);

  if (do_verification)
    {
      // We can run the verification for the rational values, but the
      // error is always identically zero (though it's good to verify
      // this).  The more interesting thing is what happens when we
      // convert the rational points and weights to real-valued
      // numbers with a fixed precision and see how the ordering
      // affects the maximum error.
      // std::cout << "\nVerifying rule..." << std::endl;
      // verify(s, x, w, indices);
      // std::cout << "... Verification complete!" << std::endl;

      std::cout << "\nVerifying real-valued rule..." << std::endl;
      mpfr_class
        max_error = 0.,
        avg_error = 0.;
      verify(s, x_real, w_real, indices, max_error, avg_error);
      std::cout << "max_error = " << max_error << std::endl;
      std::cout << "avg_error = " << avg_error << std::endl;
      std::cout << "... Verification complete!" << std::endl;
    }

  return 0;
}



void usage()
{
  std::cout << "\n";
  std::cout << "This program generates the points and weights for a 3D Grundmann-Moller quadrature rule.\n";
  std::cout << "\n";
  std::cout << "Valid command line options are:\n";
  std::cout << "--rule-index, -s #      = Build GM rule with index s, degree 2*s+1.\n";
  std::cout << "--no-verify, -x         = Skip verifying that the rule exactly integrates polynomials with degree<=2*s+1.\n";
  std::cout << "--unsorted, -u          = Don't sort the points and weights.\n";
  std::cout << "--binary-digits, -b #   = Default number of binary digits to use for mpfr_class variables.\n";
  std::cout << "--n-random-trials, -r # = Specify the number of random trials to attempt.  If>0, --unsorted is implied.\n";
  std::cout << "--random-seed, -e       = Specify a random seed for use with --n-random-trials>0.\n";
  std::cout << "--interleave, -i        = Use in conjunction with less-than sorting to interleave the entries of the sorted array.\n";
  std::cout << "--optimize, -o #        = Specify the number of optimization trials to attempt.  The --random-seed is used to initialize the algorithm.\n";
  std::cout << "--permute, -p           = Compute the sum of weights for all possible permutations of the weights vector.\n";
  std::cout << "--help, -h              = Print this message.\n";
  std::cout << "\n";
}



template <class T>
void verify(unsigned rule_index,
            const std::vector<Point<T> >& x,
            const std::vector<T>& w,
            const std::vector<size_t>& indices,
            T & max_error,
            T & avg_error)
{
  // When tracking round-off errors, use higher precision than the inputs.
  unsigned long double_input_precision = 2*w[0].get_prec();

  // Keep track of the maximum and average error encountered (use higher precision)
  T sum_error(/*value=*/0, /*precision=*/double_input_precision);

  // Use 2x the precision of the inputs while summing the errors.
  // This will simply be 2x the default precision if T is of type
  // mpq_class.  Note that changing the precision, invalidates these
  // (sets them to NaN) but that's OK because we set a valid value in
  // them below.
  max_error.set_prec(double_input_precision);
  avg_error.set_prec(double_input_precision);

  // Keep track of the number of calculations performed for the average.
  unsigned num_calcs = 0;

  // Mostly copied from conical_product_3D.C, the verification code
  // for tets should probably be factored out into a common
  // function...
  unsigned max_order = 2*rule_index + 1;
  for (unsigned x_power=0; x_power<max_order; ++x_power)
    for (unsigned y_power=0; y_power<max_order; ++y_power)
      for (unsigned z_power=0; z_power<max_order; ++z_power)
        {
          // Only try to integrate polynomials we can integrate exactly
          if (x_power + y_power + z_power > max_order)
            continue;

          // Compute the integral using quadrature
          T sum = 0;
          for (unsigned n_qp=0; n_qp<x.size(); ++n_qp)
            {
              unsigned mapped_index = indices[n_qp];

              sum +=
                w[mapped_index] *
                pow(x[mapped_index](0), x_power) *
                pow(x[mapped_index](1), y_power) *
                pow(x[mapped_index](2), z_power);
            }

          // // Compute the quadrature result, but save all the summands
          // // into a vector first for inspection.
          // std::vector<T> summands(x.size());
          // for (unsigned n_qp=0; n_qp<x.size(); ++n_qp)
          //   {
          //     unsigned mapped_index = indices[n_qp];
          //
          //     // Note: we store them into the summands vector via n_qp, not mappe_index, otherwise they will be stored in the original "unsorted" order!
          //     summands[n_qp] =
          //       w[mapped_index] *
          //       pow(x[mapped_index](0), x_power) *
          //       pow(x[mapped_index](1), y_power) *
          //       pow(x[mapped_index](2), z_power);
          //   }
          //
          // // If you print the summands, you see that they are not
          // // necessarily in the same ordering as the weights were
          // // placed.  And since these are the numbers that we actually
          // // add together, this is where the sorting to prevent
          // // round-off error needs to happen.
          // // for (unsigned i=0; i<x.size(); ++i)
          // //   std::cout << "summand[" << i << "]=" << summands[i] << std::endl;;
          //
          // // Note: The summands are stored in sorted order of the
          // // weights, so there's no need to map indices on the second
          // // loop.
          // T sum = 0;
          // for (unsigned i=0; i<x.size(); ++i)
          //   sum += summands[i];

          // Debugging
          // std::cout << "quadrature result = " << sum << std::endl;

          // Compute the exact solution as described above, divide
          // term-by-term to avoid huge numbers/overflows (even though
          // this is GMP).  Use higher precision when computing the
          // analytical solution.
          T analytical(/*value=*/1, /*precision=*/double_input_precision);

          // Compute the analytical integral value.
          analytical = exact_tet(x_power, y_power, z_power);

          // Debugging
          // std::cout << "analytical = " << analytical << std::endl;

          // Compute the absolute error (using higher precision):
          T abs_err(/*value=*/my_abs(T(sum - analytical)), /*precision=*/double_input_precision);

          // Debugging
          // std::cout << "abs_err = " << abs_err << std::endl;

          // Possibly update the max error
          max_error = max(abs_err, max_error);

          // Increment the sum and the calculation count, so we can
          // report the average error at the end.
          sum_error += abs_err;
          num_calcs++;

          // Debugging
          // std::cout << "sum_error = " << sum_error << std::endl;

          // Print message.  In 3D, this is just way too much output.
          // std::cout << "Computing integral of: x^" << x_power << " y^" << y_power << " z^" << z_power
          //           << ", abs_err = " << abs_err << std::endl;

          // Abort if error is too large.
          // if (abs_err > T(1.e-30))
          //   {
          //     std::cerr << "Quadrature error too large, possible problem with points and weights!" << std::endl;
          //     std::abort();
          //   }
        } // end triple-nested for loop

  // Compute the average error
  avg_error = sum_error / num_calcs;
}



template <class T>
void print_ordered_rule(const std::vector<Point<T> >& x,
                        const std::vector<T>& w,
                        const std::vector<size_t>& indices)
{
  // Generate C++ code to resize the _points and _weights vectors appropriately.
  std::cout << "_points.resize(" << x.size() << ");" << std::endl;
  std::cout << "_weights.resize(" << x.size() << ");\n" << std::endl;

  // Determine how many spaces to use to make everything line up in the printout
  unsigned print_width = 1 + static_cast<unsigned>(std::log10(static_cast<double>(x.size())));

  // Pre-compute the widest entry of each column.  We don't have to sort to do this...
  size_t
    widest_x = 0,
    widest_y = 0,
    widest_z = 0,
    widest_w = 0;
  for (unsigned i=0; i<x.size(); ++i)
    {
      widest_x = std::max(widest_x, format_rational(x[i](0)).size());
      widest_y = std::max(widest_y, format_rational(x[i](1)).size());
      widest_z = std::max(widest_z, format_rational(x[i](2)).size());
      widest_w = std::max(widest_w, format_rational(w[i]).size());
    }

  // Print the points and weights
  for (unsigned i=0; i<x.size(); ++i)
    {
      // Choose the sorted or "natural" ordering
      unsigned mapped_index = indices[i];

      // Print points in C++ code that can be used in libmesh
      std::cout << "_points[" << std::setw(print_width) << i << "](0) = " << std::setw(widest_x) << format_rational(x[mapped_index](0)) << ";  "
                << "_points[" << std::setw(print_width) << i << "](1) = " << std::setw(widest_y) << format_rational(x[mapped_index](1)) << ";  "
                << "_points[" << std::setw(print_width) << i << "](2) = " << std::setw(widest_z) << format_rational(x[mapped_index](2)) << ";  "
                << "_weights[" << std::setw(print_width) << i << "] = "   << std::setw(widest_w) << format_rational(w[mapped_index]) << ";"
                << std::endl;
    }
}
