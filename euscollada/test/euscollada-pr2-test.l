(load "unittest.l")
(init-unit-test)
(load "package://euscollada/pr2.l")
(load "package://pr2eus/pr2-utils.l")

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
   (send *pr2* :inverse-kinematics p
         :look-at-target t
         :rotation-axis :z
         :use-base 0.1
         )
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
	 (ret
	  (eps= (* 1e-3 (send *pr2* :weight)) ;; [g] -> [kg]
		(reduce #'+ (read-from-string (format nil "(~A )" (read-line strm nil nil))))
		)))
    (close strm)
    (assert ret)))

(defun tmp-eusmodel-validity-check
  (robot)
  (labels ((root-link-validity-check
            ()
            (let ((ret (null (send (car (send robot :links)) :parent-link))))
              (unless ret (warn ";; root link should have no parent-link!!~%"))
              ret))
           (link-joint-length-check-for-serial-link-manipulator
            ()
            (let ((ret
                   (= (1+ (length (mapcar #'cdr (remove-if-not #'(lambda (s) (derivedp (cdr s) joint)) (send robot :slots)))))
                      (length (remove-duplicates (mapcar #'cdr (remove-if-not #'(lambda (s) (derivedp (cdr s) bodyset-link)) (send robot :slots))))))))
              (unless ret (warn ";; link = joint + 1 <- for serial link manipulator!!~%"))
              ret))
           (chain-validity-check
            ()
            (let ((jl (mapcar #'cdr (remove-if-not #'(lambda (s) (derivedp (cdr s) joint)) (send robot :slots)))))
              (every #'identity
                     (mapcar #'(lambda (j)
                                 (let ((ret
                                        (and
                                         ;; joint should have child-link derived from bodyset-link class
                                         (send j :child-link) (derivedp (send j :child-link) bodyset-link)
                                         ;; joint should have parent-link derived from bodyset-link class
                                         (send j :parent-link) (derivedp (send j :parent-link) bodyset-link)
                                         ;; child-link should associated with parent-link
                                         (member (send j :child-link) (send (send j :parent-link) :descendants))
                                         (equal (send j :parent-link) (send (send j :child-link) :parent))
                                         ;; definition of child-link and parent-link should consistent in joints and links
                                         (equal (send (send j :child-link) :parent-link) (send j :parent-link))
                                         (member (send j :child-link) (send (send j :parent-link) :child-links)))))
                                   (unless ret (warn ";; joint-link relationship at ~A is broken!!~%" (send j :name)))
                                   ret))
                             jl))
              )))
    (and (root-link-validity-check)
         (link-joint-length-check-for-serial-link-manipulator)
         (chain-validity-check))
    ))

(deftest pr2-eusmodel-validity-check
  (assert (tmp-eusmodel-validity-check *pr2*)))

(run-all-tests)
(exit)
