#include <fstream>
#include "read_rule.h"

// Given an input filename and a Rule object, fills up the Rule with
// the contents of the file.
void read_rule(const std::string & filename, Rule & rule)
{
  // Create stream
  std::ifstream in(filename.c_str());

  // Read line by line
  std::string s;

  // Keep track of the number of different types of generators as we
  // encounter them in the file
  unsigned
    centroid_count  = 0,
    median_count    = 0,
    arbitrary_count = 0;

  while (true)
    {
      // Try to read something.  This may set EOF!
      std::getline(in, s);

      if (in)
        {
          // Process s, then continue

          // Delete all whitespace from the line
          s.erase(std::remove_if(s.begin(), s.end(), isspace), s.end());

          if ((s.find("Dup3(") != std::string::npos) ||
              (s.find("Perm3(") != std::string::npos))
            centroid_count++;

          else if ((s.find("Dup21(") != std::string::npos) ||
                   (s.find("Perm21(") != std::string::npos))
            median_count++;

          else if ((s.find("Dup111(") != std::string::npos) ||
                   (s.find("Perm111(") != std::string::npos))
            arbitrary_count++;

          continue;
        } // if (in)

      // If !file, check to see if EOF was set.  If so, break out
      // of while loop.
      if (in.eof())
        break;

      // If !in and !in.eof(), stream is in a bad state!
      std::cerr << "Stream is bad! Perhaps the file: "
                << filename
                << " does not exist?"
                << std::endl;
      std::abort();

    } // while

  // There should be two lines for each generator, so ensure that happened
  if ((centroid_count % 2 != 0) ||
      (median_count % 2 != 0)   ||
      (arbitrary_count % 2 != 0))
    {
      std::cerr << "Error reading the generators file." << std::endl;
      std::abort();
    }
  else
    {
      centroid_count /= 2;
      median_count /= 2;
      arbitrary_count /= 2;
    }

  std::cout << "Rule has "
            << centroid_count << " centroid generator(s), "
            << median_count << " median generator(s), and "
            << arbitrary_count << " arbitrary generator(s)."
            << std::endl;

  // We now know how many of each type of generator there are, so we can prepare the Rule object.
  for (unsigned i=0; i<centroid_count; ++i)
    rule.push_back(Generator(Generator::CENTROID));

  for (unsigned i=0; i<median_count; ++i)
    rule.push_back(Generator(Generator::MEDIAN));

  for (unsigned i=0; i<arbitrary_count; ++i)
    rule.push_back(Generator(Generator::ARBITRARY));

  // OK, now that the Rule is set up, let's read the file again and
  // extract values.  When you rewind, you also have to call clear()
  // before calling seekg because once it reaches the end of stream it
  // will have error flags set.
  in.clear();
  in.seekg(0);

  // Keep track of the index of the 'next' index in the rule of each
  // type that is "next" to receive a point or weight.  This way we
  // can handle the inputs no matter what order they are in in the
  // file.
  unsigned
    centroid_weight_index  = 0,
    median_weight_index    = centroid_count,
    median_point_index     = centroid_count,
    arbitrary_weight_index = centroid_count + median_count,
    arbitrary_point_index  = centroid_count + median_count;

  while (true)
    {
      // Try to read something.  This may set EOF!
      std::getline(in, s);

      if (in)
        {
          // Process s, then continue

          // Delete all whitespace from the line
          s.erase(std::remove_if(s.begin(), s.end(), isspace), s.end());

          if (s.find("Dup3(") != std::string::npos)
            rule[centroid_weight_index++].get_w() = extract_number(s, '(', ')');

          // Lines that start with Perm3( are ignored... they are always 1/3

          else if (s.find("Dup21(") != std::string::npos)
            rule[median_weight_index++].get_w() = extract_number(s, '(', ')');

          else if (s.find("Perm21(") != std::string::npos)
            rule[median_point_index++].get_a() = extract_number(s, '(', ')');

          else if (s.find("Dup111(") != std::string::npos)
            rule[arbitrary_weight_index++].get_w() = extract_number(s, '(', ')');

          else if (s.find("Perm111(") != std::string::npos)
            {
              // Get a reference to the next CentroidGenerator that needs a weight
              Generator & generator = rule[arbitrary_point_index];

              // Set the weight
              generator.get_a() = extract_number(s, '(', ',');
              generator.get_b() = extract_number(s, ',', ')');

              // Get ready to assign the next arbitrary generator points
              arbitrary_point_index++;
            }

          continue;
        } // if (in)

      // If !file, check to see if EOF was set.  If so, break out
      // of while loop.
      if (in.eof())
        break;

      // If !in and !in.eof(), stream is in a bad state!
      std::cerr << "Stream is bad! Perhaps the file: "
                << filename
                << " does not exist?"
                << std::endl;
      std::abort();

    } // while
}



// Extracts a number from a string (and possibly removes 'L'
// characters) from between the given 'first_delim' and 'second_delim'
// characters.
std::string extract_number(const std::string & input, char first_delim, char second_delim)
{
  // Extract the position of the opening and closing parentheses
  size_t
    beg_pos = input.find(first_delim),
    end_pos = input.find(second_delim);

  // Make sure that both delimiters were found
  if (beg_pos == std::string::npos || end_pos == std::string::npos)
    {
      std::cerr << "Cannot extract string, one of the delimiters was not found!" << std::endl;
      std::abort();
    }

  // Extract a substring
  std::string number = input.substr(beg_pos+1, end_pos - beg_pos - 1);

  // mpfr_class does not like a trailing "L" character
  // which I use in some of my formatting...
  // Note: In C++11, std::string has back() and pop_back() members
  // that make this O(1) since we know 'L' only appears at the end of
  // the string...
  number.erase(std::remove(number.begin(), number.end(), 'L'), number.end());

  return number;
}
