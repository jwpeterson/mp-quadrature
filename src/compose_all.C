#include <assert.h>
#include "compose_all.h"

void compose_all(unsigned s,
                 unsigned p,
                 std::vector<std::vector<unsigned> >& result)
{
  // Clear out results remaining from previous calls
  result.clear();

  // Allocate storage for a workspace.  The workspace will periodically
  // be copied into the result container.
  std::vector<unsigned int> workspace(p);

  // The first result is always (s,0,...,0)
  workspace[0] = s;
  result.push_back(workspace);

  // the value of the first non-zero entry
  unsigned int head_value=s;

  // When head_index=-1, it refers to "off the front" of the array.  Therefore,
  // this needs to be a regular int rather than unsigned.  I initially tried to
  // do this with head_index unsigned and an else statement below, but then there
  // is the special case: (1,0,...,0) which does not work correctly.
  int head_index = -1;

  // At the end, all the entries will be in the final slot of workspace
  while (workspace.back() != s)
    {
      // If the previous head value is still larger than 1, reset the index
      // to "off the front" of the array
      if (head_value > 1)
        head_index = -1;

      // Either move the index onto the front of the array or on to
      // the next value.
      head_index++;

      // Get current value of the head entry
      head_value = workspace[head_index];

      // Put a zero into the head_index of the array.  If head_index==0,
      // this will be overwritten in the next line with head_value-1.
      workspace[head_index] = 0;

      // The initial entry gets the current head value, minus 1.
      // If head_value > 1, the next loop iteration will start back
      // at workspace[0] again.
      assert (head_value > 0);
      workspace[0] = head_value - 1;

      // Increment the head+1 value
      workspace[head_index+1] += 1;

      // Save this composition in the results
      result.push_back(workspace);
    }
}
