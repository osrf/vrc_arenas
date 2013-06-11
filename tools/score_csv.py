#!/usr/bin/env python

import os
import sys

def help():
    print("""Help:
\t%s test-file.xml
"""%(NAME), file=sys.stderr)
    print(sys.argv)
    sys.exit(getattr(os, 'EX_USAGE', 1))

def main():
  if len(sys.argv) < 2:
    usage()
  score_file = sys.argv[1]

  with open(score_file, 'r') as fh:
    first_line = fh.readline()
    for line in fh:
        pass

   
        
    print("Checking for test results in %s" % test_file)
    
    if not os.path.exists(test_file):
        if not os.path.exists(os.path.dirname(test_file)):
            os.makedirs(os.path.dirname(test_file))
            
        print("Cannot find results, writing failure results to", test_file)
        
        with open(test_file, 'w') as f:
            test_name = os.path.basename(test_file)
            d = {'test': test_name, 'test_file': test_file }
            f.write("""<?xml version="1.0" encoding="UTF-8"?>
<testsuite tests="1" failures="1" time="1" errors="0" name="%(test)s">
  <testcase name="test_ran" status="run" time="1" classname="Results">
    <failure message="Unable to find test results for %(test)s, test did not run.\nExpected results in %(test_file)s" type=""/>
  </testcase>
</testsuite>"""%d)

if __name__ == '__main__':
    main()
