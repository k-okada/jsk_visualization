#!/usr/bin/env python

import sys
import unittest
import rosgraph
import rosunit
import rospy
from geometry_msgs.msg import Point

NAME = "test_rqt_plugins"
WAIT_TIMEOUT = 10.0 # sec

class TestRqtPlugins(unittest.TestCase):

    def is_subscribed(self, topic, type, node):
        try:
            rospy.wait_for_message(topic, type, WAIT_TIMEOUT) # raise ROSException if fails
        except rospy.ROSException:
            self.assertTrue(False, 'Could not subscribe topic ({})'.format(topic))
        subs = []
        i = 0
        while subs == [] and i < 10:
            pubs, subs, _ =rosgraph.Master('/rostopic').getSystemState()
            subs = [x[1] for x in subs if x[0] == topic]
            rospy.sleep(1.0)
            i = i + 1
        if subs == []:
            self.assertTrue(False, 'No node subscribes topic ({})'.format(topic))
        if not node in subs[0]:
            self.assertTrue(False, 'Node ({}) does not subscribe topic ({}), but {} does'.format(node, topic, subs))
        return True

    def __init__(self, *args):
        unittest.TestCase.__init__(self, *args)
        rospy.init_node(NAME)

    def test_rqt_3d_plot(self):
        self.assertTrue(self.is_subscribed('/dummy_3d_data', Point, '/rqt_3d_plot'))
            
    def test_rqt_yn_btn(self):
        try:
            rospy.wait_for_service('rqt_yn_btn', WAIT_TIMEOUT) # raise ROSException if fails
            self.assertTrue(True)
        except rospy.ROSException:
            self.assertTrue(False, 'Could not find service (rqt_yn_btn)')

if __name__ == '__main__':
    rosunit.unitrun(NAME, sys.argv[0], TestRqtPlugins)
