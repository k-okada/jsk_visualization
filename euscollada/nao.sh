#!/bin/bash

rosrun collada_urdf_jsk_patch urdf_to_collada `rospack find nao_description`/urdf/nao_robot.xml nao.dae
if [ "$?" != 0 ] ;  then exit ; fi

rosrun collada_eus collada2eus nao.dae nao.yaml nao.l
if [ "$?" != 0 ] ;  then exit ; fi

irteusgl -e "(progn (load \"~/prog/eus/irteus/irtmodel.l\")(load \"nao.l\")(nao)(make-irtviewer)(objects (list *nao*)))"
