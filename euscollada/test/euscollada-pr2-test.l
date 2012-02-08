(load "unittest.l")
(load "package://euscollada/pr2.l")
(load "package://pr2eus/pr2-utils.l")

(init-unit-test)

(deftest pr2-test ()
  (if (not (boundp '*pr2*)) (pr2))
  (send *pr2* :move-to (make-coords) :world)
  (send *pr2* :reset-pose)
  (if (eq x::*display* 0) (setq x::*display* nil))
  (if x::*display* (objects (list *pr2*)))
  (dotimes (i 100)
   (print (list i (send *pr2* :torso :waist-z :joint-angle)))
   (setq p (make-coords :pos
                        (v+ (float-vector 400 -400 1000)
                            (float-vector
                             0
                             (* 300 (sin (* pi (/ i 25.0))))
                             (* 500 (sin (* pi (/ i 50.0))))))))
   (assert
    (send *pr2* :inverse-kinematics p
	  :look-at-target t
	  :rotation-axis :z
	  :use-base 0.1
	  )
    (format nil "check ik for ~A" p))
   (when x::*display*
     (send *irtviewer* :objects (list *pr2* p))
     (send *irtviewer* :draw-objects)
     (x::window-main-one))
   )
  (if x::*display*
    (send-all (send *pr2* :links) :draw-on :flush t))
  )

(deftest pr2-weight-test
  (let* ((strm (piped-fork "echo $(grep \"<mass\" `rospack find pr2_mechanism_model`/pr2.urdf | cut -d\\\" -f2)"))
	 (urdf-weight (reduce #'+ (read-from-string (format nil "(~A )" (read-line strm nil nil)))))
	 (euslisp-weight (* 1e-3 (send *pr2* :weight)))) ;; [g] -> [kg]
    (close strm)
    (assert (eps= urdf-weight euslisp-weight) (format nil "check weight urdf:~A euslisp:~A" urdf-weight euslisp-weight))))

(defun eusmodel-validity-check (robot)
  ;; root-link-validity-check
  (let ((root-link (car (send robot :links))))
    (assert (null (send root-link :parent-link))
	    (format nil "root link should have no parent-link!! ~A ~A~%" root-link (send root-link :parent-link)))
    )

  ;;link-joint-length-check-for-serial-link-manipulator
  (let ((joint-length (length (mapcar #'cdr (remove-if-not #'(lambda (s) (derivedp (cdr s) joint)) (send robot :slots)))))
	(link-length (length (remove-duplicates (mapcar #'cdr (remove-if-not #'(lambda (s) (derivedp (cdr s) bodyset-link)) (send robot :slots)))))))
    (assert (= (+ joint-length 1) link-length)
	    (format nil ";; link(~A) = joint(=A) + 1 <- for serial link manipulator!!~%" link-length joint-length)))

  ;; chain-validity-check
  (let ((jl (mapcar #'cdr (remove-if-not #'(lambda (s) (derivedp (cdr s) joint)) (send robot :slots)))))
    (dolist (j jl)
      ;; joint should have child-link derived from bodyset-link class
      (assert (and (send j :child-link) (derivedp (send j :child-link) bodyset-link))
	      (format nil "joint should have child-link derived from bodyset-link class ~A ~A" j (send j :child-link)))
      ;; joint should have parent-link derived from bodyset-link class
      (assert (and (send j :parent-link) (derivedp (send j :parent-link) bodyset-link))
	      (format nil "child-link should associated with parent-link ~A ~A"
		      j (send j :parent-link)))
      ;; definition of child-link and parent-link should consistent in joints and links
      (assert (and (member (send j :child-link) (send (send j :parent-link) :descendants))
		   (equal (send j :parent-link) (send (send j :child-link) :parent)))
	      (format nil "definition of child-link and parent-link should consistent in joints and links ~A ~A ~A" j (send j :parent-link) (send j :child-link)))
      ;; definition of child-link and parent-link should consistent in joints and links
      (assert (and (equal (send (send j :child-link) :parent-link) (send j :parent-link))
		   (member (send j :child-link) (send (send j :parent-link) :child-links)))
	      (format nil "definition of child-link and parent-link should consistent in joints and links ~A ~A ~A" j (send j :child-link) (send j :parent-link)))
      )))

(deftest pr2-eusmodel-validity-check
  (eusmodel-validity-check *pr2*))

(run-all-tests)
(exit)
