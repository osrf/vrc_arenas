#!/usr/bin/env python
#
# Software License Agreement (Apache License)
#
# Copyright 2013 Open Source Robotics Foundation
# Author: Brian Gerkey
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import roslib; roslib.load_manifest('vrc_cheats')
import rospy
import sys
from osrf_msgs.msg import JointCommands
from sensor_msgs.msg import Joy

EIGEN_POSES = {
    'tabletop_reach_r': [[0.0, 
                          0.0, 
                          0.0, 
                          0.0, 
                          0.0, 
                          0.0, 
                          0.0, 
                          0.0, 
                          0.0, 
                          0.0, 
                          0.0, 
                          0.0, 
                          0.0, 
                          0.0, 
                          0.0, 
                          0.0, 
                          0.0, 
                          0.0, 
                          0.0, 
                          0.0, 
                          0.0, 
                          0.0, 
                          0.0, 
                          0.0, 
                          0.0, 
                          0.0, 
                          0.0, 
                          0.0],
                         [0.0,
                          0.0,
                          0.0,
                          0.0,
                          0.0,
                          0.0,
                          0.0,
                          0.0,
                          0.0,
                          0.0,
                          0.0,
                          0.0,
                          0.0,
                          0.0,
                          0.0,
                          0.0,
                          0.0,
                          0.0,
                          0.0,
                          0.0,
                          0.0,
                          0.0,
-0.7668111815691914, 1.396255493523947, 1.9044557076508042, -1.1626445233703198, 0.023734655893724366, -0.22264687988712772]]
}

# Axis to move name; activate button for each axis is axis number * 2
AXIS_MAP = {
    0: 'tabletop_reach_r'
}

# Scene button comes last
SCENE_BUTTON_INDEX = -1

USAGE = 'eigen_arm.py [scene]'

class EigenArm:

    def __init__(self, argv):
        # If specified, we'll only do anything when the scene number matches
        self.scene = None
        self.arms = []
        self.scene = None

        # Prepare a message structure that we'll reuse
        self.atlasJointNames = [ 'atlas::back_lbz', 'atlas::back_mby',
                                 'atlas::back_ubx', 'atlas::neck_ay',
                                 'atlas::l_leg_uhz', 'atlas::l_leg_mhx',
                                 'atlas::l_leg_lhy', 'atlas::l_leg_kny',
                                 'atlas::l_leg_uay', 'atlas::l_leg_lax',
                                 'atlas::r_leg_uhz', 'atlas::r_leg_mhx',
                                 'atlas::r_leg_lhy', 'atlas::r_leg_kny',
                                 'atlas::r_leg_uay', 'atlas::r_leg_lax',
                                 'atlas::l_arm_usy', 'atlas::l_arm_shx',
                                 'atlas::l_arm_ely', 'atlas::l_arm_elx',
                                 'atlas::l_arm_uwy', 'atlas::l_arm_mwx',
                                 'atlas::r_arm_usy', 'atlas::r_arm_shx',
                                 'atlas::r_arm_ely', 'atlas::r_arm_elx',
                                 'atlas::r_arm_uwy', 'atlas::r_arm_mwx']
        self.joint_command = JointCommands()
        self.joint_command.name = list(self.atlasJointNames)
        n = len(self.joint_command.name)
        self.joint_command.position = [0.0] * n
        self.joint_command.velocity = [0.0] * n
        self.joint_command.effort = [0.0] * n
        self.joint_command.kp_position = [0.0] * n
        self.joint_command.ki_position = [0.0] * n
        self.joint_command.kd_position = [0.0] * n
        self.joint_command.kp_velocity = [0.0] * n
        self.joint_command.i_effort_min = [0.0] * n
        self.joint_command.i_effort_max = [0.0] * n

        if len(argv) > 2:
            self.usage()
        if len(argv) == 2:
            self.scene = int(argv[1])

        rospy.init_node('eigen_arm', anonymous=True)

        # now get gains from parameter server
        for i in xrange(len(self.joint_command.name)):
            name = self.joint_command.name[i]
            self.joint_command.kp_position[i]  = rospy.get_param('atlas_controller/gains/' + name[7::] + '/p')
            self.joint_command.ki_position[i]  = rospy.get_param('atlas_controller/gains/' + name[7::] + '/i')
            self.joint_command.kd_position[i]  = rospy.get_param('atlas_controller/gains/' + name[7::] + '/d')
            self.joint_command.i_effort_max[i] = rospy.get_param('atlas_controller/gains/' + name[7::] + '/i_clamp')
            self.joint_command.i_effort_min[i] = -self.joint_command.i_effort_max[i]

        self.joy_sub = rospy.Subscriber('joy', Joy, self.joy_cb)
        self.pub = rospy.Publisher('atlas/joint_commands', JointCommands)

    def usage(self):
        print USAGE
        sys.exit(1)

    def joy_cb(self, msg):
        if (self.scene is not None and 
            msg.buttons[SCENE_BUTTON_INDEX] != self.scene):
            return
        pose = None
        amount = None
        for k,v in AXIS_MAP.iteritems():
            if msg.buttons[k*2] == 1:
                pose = v
                amount = msg.axes[k]
                break
        if pose is None:
            return
        poses = EIGEN_POSES[pose]
        for i in xrange(0,len(poses[0])):
            self.joint_command.position[i] = poses[0][i] + amount * (poses[1][i] - poses[0][i])
        self.joint_command.header.stamp = rospy.Time.now()
        self.pub.publish(self.joint_command)

if __name__ == '__main__':
    ea = EigenArm(sys.argv)
    rospy.spin()

