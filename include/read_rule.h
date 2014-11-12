#ifndef READ_RULE_H
#define READ_RULE_H

#include <string>
#include "rule.h"

// Given an input filename and a Rule object, fills up the Rule with
// the contents of the file.
void read_rule(const std::string & filename, Rule & rule);

// Extracts a number from a string (and possibly removes 'L'
// characters) from between the given 'first_delim' and 'second_delim'
// characters.
std::string extract_number(const std::string & input, char first_delim, char second_delim);

#endif

