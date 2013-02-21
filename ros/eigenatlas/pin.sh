#!/bin/bash
rostopic pub --once /atlas/mode std_msgs/String no_gravity &
sleep 0.1
rostopic pub --once /atlas/set_pose geometry_msgs/Pose '{position: {z: 1}, orientation: {w: 1}}' &
sleep 0.1
rostopic pub --once /atlas/mode std_msgs/String pinned &
