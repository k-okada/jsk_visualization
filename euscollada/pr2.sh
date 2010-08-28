#!/bin/bash

rosrun collada_urdf_jsk_patch urdf_to_collada `rospack find pr2_mechanism_model`/pr2.urdf pr2.dae
if [ "$?" != 0 ] ;  then exit ; fi

rosrun euscollada collada2eus pr2.dae pr2.yaml pr2.l
if [ "$?" != 0 ] ;  then exit ; fi

rosrun euslisp irteusgl -e \
            "(progn (load \"pr2.l\") \
                    (pr2) \
                    (send *pr2* :reset-pose) \
                    (make-irtviewer)(objects (list *pr2*)) \
                    (send *pr2* :larm :inverse-kinematics (send *pr2* :rarm :end-coords) :debug-view t :rotation-axis nil :stop 100) \
                    (send *pr2* :head :look-at (send *pr2* :rarm :end-coords :worldpos) :debug-view t)
                    (send-all (send *pr2* :links) :draw-on :flush t)
\
             )"
