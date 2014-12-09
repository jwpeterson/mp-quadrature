#!/user/bin/env python
from sys import argv, exit, stdout
import subprocess

# This script runs the drivers/test_read_rule with all the input files
# for a given order and records the results.

# Define a function that actually runs the code
def run_code(filenames):
  ncols = len(filenames)
  all_results = []
  for i in range(0,ncols):
    # Popen is very particular about how you pass the arguments... they can't have any spaces in them at all?
    p = subprocess.Popen(['./drivers/tri_rule', '-i', filenames[i]],
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE,
                         stdin=subprocess.PIPE)

    out, err = p.communicate()

    # If err is an empty string, process out
    if not err:
      # This line does a lot: It gets stdout from the subprocess, splits
      # it on newline characters, and appends the result as a new list
      # in the all_results list-of-lists.
      all_results.append(out.splitlines())
    else:
      print err
      exit(1)

  # Print the filenames as a comment line
  print '#',
  for i in range(0, ncols):
    print filenames[i],
  stdout.write('\n')

  # Print data in columns.  Note that we're using stdout.write instead
  # of print to avoid having an extra space inserted before commas...
  nrows = len(all_results[0])
  for j in range(0, nrows):
    # Note: we are using j-1 because there is a line of garbage at the
    # top of the output.  Change to j if we remove that...
    stdout.write(str(j-1) + ', ')
    for i in range(0, ncols):
      stdout.write(str(all_results[i][j]) + ', ')
    # Print a newline
    stdout.write('\n')

# Handle command line args, We only allow 1: the degree of the rules to compute.
if (len(argv) > 2):
  print "Error! Can only pass 1 argument, which is the degree of the rule to be computed!"
  exit(1)

if (len(argv) == 2):
  n = int(argv[1])
else:
  n = 6

filenames = []

if (n==2):
  filenames = ['inputs/quad_2d_p02_bndry.out', 'inputs/quad_2d_p02_interior.out']

elif (n==3):
  filenames = ['inputs/quad_2d_p03_NI.out', 'inputs/quad_2d_p03_CP.out']

elif (n==6):
  filenames = ['inputs/quad_2d_p06_hompack1.out', 'inputs/quad_2d_p06_dunavant.out']

elif (n==7):
  filenames = ['inputs/quad_2d_p07_dunavant.out', 'inputs/quad_2d_p07_hompack1.out', 'inputs/quad_2d_p07_mine15.out',
               'inputs/quad_2d_p07_zhang_underdetermined.in', 'inputs/quad_2d_p07_gatermann.in']

elif (n==10):
  filenames = ['inputs/quad_2d_p10_dunavant.out', 'inputs/quad_2d_p10_hompack1.out', 'inputs/quad_2d_p10_hompack2.out', 'inputs/quad_2d_p10_zhang_underdetermined.in']

elif (n==11):
  filenames = ['inputs/quad_2d_p11_mine.out', 'inputs/quad_2d_p11_hompack1.out', 'inputs/quad_2d_p11_hompack2.out', 'inputs/quad_2d_p11_zhang_underdetermined.in']

elif (n==12):
  filenames = ['inputs/quad_2d_p12_mine.out', 'inputs/quad_2d_p12_dunavant.out']

elif (n==15):
  filenames = ['inputs/quad_2d_p15_mine49.out', 'inputs/quad_2d_p15_mine51.out', 'inputs/quad_2d_p15_mine52c.out',
               'inputs/quad_2d_p15_zhang_witherden.out', 'inputs/quad_2d_p15_wandzura.in']

elif (n==16):
  filenames = ['inputs/quad_2d_p16_mine55.out', 'inputs/quad_2d_p16_zhang_underdetermined.in']

elif (n==17):
  filenames = ['inputs/quad_2d_p17_mine63c.out', 'inputs/quad_2d_p17_mine63e.out', 'inputs/quad_2d_p17_dunavant_underdetermined.in', 'inputs/quad_2d_p17_zhang_xiao_underdetermined.in']

else:
  print "Error, unrecognized value n=%d" % n
  exit(1)

run_code(filenames)
