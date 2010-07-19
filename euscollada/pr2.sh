#!/bin/bash

rosrun collada_urdf urdf_to_collada `rospack find pr2_mechanism_model`/pr2.urdf pr2.dae
rosrun collada_eus collada2eus pr2.dae pr2.yaml pr2.l

irteusgl -e "(progn (load \"pr2.l\")(setq *pr2* (instance pr2-robot :init))(make-irtviewer)(objects (list *pr2*)))"
