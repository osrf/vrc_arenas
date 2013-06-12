#!/usr/bin/env python

import os
import sys

def help():
  print("""score_csv.py test-file.xml""")
  print(sys.argv)
  sys.exit(getattr(os, 'EX_USAGE', 1))

def main():
  if len(sys.argv) < 2:
    help()

  for dirname in os.listdir(sys.argv[1]):
    team_id = dirname
    finals_dir = os.path.join(sys.argv[1], dirname, "final")
    for sub_dirname in os.listdir(finals_dir):
      print "Processing: " + finals_dir + "/" + sub_dirname 
      [vrc, final, task_id_str, run_id_str] = sub_dirname.split("_")
      run_id = (int(task_id_str) - 1) * 5 + int(run_id_str)
      print "Team:" + team_id
      print "Run : %d" % run_id
     
     # tar = tarfile.open("test.tar")
     # for member in tar.getmembers():
     #   f=tar.extractfile(member)
     #   content=f.read()
     #   print "%s has %d newlines" %(member, content.count("\n"))
     #   print "%s has %d spaces" % (member,content.count(" "))
     #   print "%s has %d characters" % (member, len(content))
     #   sys.exit()
     # tar.close()

  #score_file = sys.argv[1]

  #team_id = 1
  #run_id = 1
  #start = "10/06/2013 12:25:01"
  #bits_up_total = 100
  #bits_down_total = 100
  #bits_up_score = 50
  #bits_down_score = 50

  #with open(score_file, 'r') as fh:
  #  first_line = fh.readline()
  #  for line in fh:
  #    line.strip()
  #    if line[0] == '#':
  #      continue
  #    [wall_time, sim_time, wall_time_elapsed, sim_time_elapsed, \
  #        score, falls, msg] = line.split(",")

  #print "%s, %d, %s, %f, %f, %f, %f, %d, %d, %d, %d, %d, %d"\
  #    % (team_id, run_id, start, float(wall_time), float(sim_time),\
  #    float(wall_time_elapsed), float(sim_time_elapsed), int(score),\
  #    int(falls), int(bits_up_total), int(bits_down_total),\
  #    int(bits_up_score), int(bits_down_score))
    
if __name__ == '__main__':
  main()
