#!/usr/bin/env python
import roslib; roslib.load_manifest('test_roslaunch')

import os
import unittest

## A sample python unit test
class TestCollada2Eus(unittest.TestCase):
    def check_euscollada(self,filename,function):
        print os.system('rosrun euscollada collada2eus test/'+filename+'.zae test/'+filename+'.l')
        callstr = 'irteusgl test/'+filename+'.l \"(objects ('+function+'))\" \"(send *viewer* :viewsurface :write-to-image-file \\"test/'+filename+'.ppm\\")\" \"(exit)\"'
        os.system(callstr)

    def test_pa10(self):
        self.check_euscollada("mitsubishi-pa10","Mitsubishi-PA10")
        self.check_euscollada("unimation-pumaarm","PUMA")
        self.check_euscollada("care-o-bot3","cob3-2")
        self.check_euscollada("darpa-arm","darpa_arm_robot")

if __name__ == '__main__':
    import rostest
    rostest.unitrun('test_collada2eus', 'test_collada2eus', TestCollada2Eus)

