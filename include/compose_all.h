#ifndef COMPOSE_ALL_H
#define COMPOSE_ALL_H

#include <vector>

/**
 * Helper function for computing compositions
 * This routine for computing compositions and their permutations is
 * originally due to:
 *
 * Albert Nijenhuis, Herbert Wilf,
 * Combinatorial Algorithms for Computers and Calculators,
 * Second Edition,
 * Academic Press, 1978,
 * ISBN: 0-12-519260-6,
 * LC: QA164.N54.
 *
 * I still don't understand exactly why it works, but it does.
 * s = number to be compositioned
 * p = number of partitions
 */
void compose_all(unsigned s,
                 unsigned p,
                 std::vector<std::vector<unsigned> >& result);

#endif
