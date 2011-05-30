#!/bin/bash

rosrun collada_urdf_jsk_patch urdf_to_collada `rospack find nao_description`/urdf/nao_robot.xml nao.dae
#rosrun collada_urdf urdf_to_collada `rospack find nao_description`/urdf/nao_robot.xml nao.dae
if [ "$?" != 0 ] ;  then exit ; fi

rosrun euscollada collada2eus nao.dae nao.yaml nao.l
if [ "$?" != 0 ] ;  then exit ; fi

rosrun euslisp irteusgl -e "(progn (load \"nao.l\")(nao)(make-irtviewer)(objects (list *nao*)))"
