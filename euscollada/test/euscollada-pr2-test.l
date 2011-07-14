(load "unittest.l")
(init-unit-test)
(load "package://euscollada/pr2.l")
(load "package://pr2eus/pr2-utils.l")

(deftest pr2-test ()
  (if (not (boundp '*pr2*)) (pr2))
  (send *pr2* :move-to (make-coords) :world)
  (send *pr2* :reset-pose)
  (if (eq x::*display* 0) (setq x::*display* nil))
  (if (and x::*display* (objects (list *pr2*)))
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


(run-all-tests)
(exit)
