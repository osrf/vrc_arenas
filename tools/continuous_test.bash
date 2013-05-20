#!/bin/bash

recipient="nate@osrfoundation.org"

while [ 1 ]; do
  task_num=`expr $RANDOM % 15 + 1`
  task="vrc_final_task$task_num"

  roslaunch vrc_finals ${task}.launch &

  # Sleep for 60 minutes
  sleep 3600
  gzlog stop
  sleep 20
  killall -INT roslaunch
  sleep 15
  num_end_tags=`grep '</gazebo_log>' /tmp/${task}/state.log | wc | awk {'print $1'}`

  if [ $num_end_tags -ne 1 ]; then
    echo "Wrong number of end tags: $num_end_tags" | mail -s "Fail Qual Task $task_num" $recipient

    echo "[FAIL] Wrong number of end tags: $num_end_tags"
  else
    echo "[PASS] Right number of end tags: $num_end_tags"
  fi

  last_line=`tail -n 1 /tmp/${task}/state.log`
  if [ $last_line != '</gazebo_log>' ]; then
    echo "Wrong last line: $last_line" | mail -s "Fail Qual Task $task_num" $recipient
    echo "[FAIL] Wrong last line: $last_line"
  else
    echo "[PASS] Correct last line: $last_line"
  fi
done
