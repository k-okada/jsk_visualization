#!/usr/bin/env roseus
;;
;; Sample of Adding Object / point_pose_extractor
;;
;;   subscribe topic -> /narrow_stereo/left/ImageFeature0D
;;   publish topic   -> /narrow_stereo/left/ObjectDetection
;;

;; < 1 >
;; . make irteus models with texture
;; . object should have unique name
(defun make-sony-rimokon ()
  (let ((obj (make-cube 45 190 25)))
    (setq aface (elt (send obj :faces) 1))
    (send obj :paste-texture-to-face aface
	  :file "sony-tv-controller.png")
    (send obj :name "sony_tv_controller")
    obj))

;; < 2 >
;; . generate launch file
(defun generate-launcher ()
  (load "package://jsk_perception/euslisp/eusmodel_template_gen_utils.l")
  (gen-all-launch (list (make-sony-rimokon)) "rimokon-pose.launch"
		  :ns "/narrow_stereo/left" ))
(warn ";; (generate-launcher)~%")

;; < 3 >
;; . launch the generated file

;; < 4 >
;; . subscribe ObjectDetection topic
(defun test-detection-result ()
  (ros::roseus "test")
  (load "package://pr2eus/pr2-interface.l")
  (load "package://elevator_move_base_pr2/src/posedetectiondb-client.l")
  (pr2)
  (setq *ri* (instance pr2-interface :init))
  (setq obj (make-sony-rimokon))
  (objects (list *pr2* obj))
  ;;
  (setq client
	(instance posedetectiondb-client :init obj
		  :topic "/narrow_stereo/left/ObjectDetection"
		  :fixed-frame "/base_footprint"))
  ;;
  (do-until-key
   (ros::spin-once)
   (send *pr2* :angle-vector (send *ri* :state :potentio-vector))
   (x::window-main-one)
   (send *irtviewer* :redraw))
  ;;
  (send client :unregister)
)
(warn ";; (test-detection-result)~%")
