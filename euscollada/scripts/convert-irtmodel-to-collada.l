#!/usr/bin/env roseus

(ros::load-ros-manifest "euscollada")
(load "package://euscollada/src/eusmodel-description.l")
(load "package://euscollada/src/collada-interface.l")

(if (not (and *model-file* *output-uri*))
    (progn
      (warn "usage: convert-to-urdf.l <modelfile> <output-directory-uri>~%")
      (exit 1)))
(format t "convert ~A and output to ~A~%" *model-file* *output-uri*)

(format t "loading ~A~%" *model-file*)
(load *model-file*)

(defun estimate-class-name (model-file)
  (string-upcase (pathname-name (pathname model-file))))
(defun remove-directory-name (fname)
  (let ((p (pathname fname)))
    (format nil "~A.~A" (pathname-name p) (pathname-type p))))
(let ((model-class-name
       (estimate-class-name
        (remove-directory-name *model-file*)))
      (full-output-path (ros::resolve-ros-path *output-uri*))
      (pkg-name (send (url-pathname *output-uri*) :host)))
  (let ((model (instance (eval (intern model-class-name)) :init)))
    (if (boundp '*model-name*) (send model :name *model-name*))
    (eus2collada model *output-uri*)))

(exit)
