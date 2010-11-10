(load "state-machine-ros.l")
(load "state-machine-sample.l")
;;
;; smach_viewer samples from SMACH tutorials
;;

(ros::roseus "smach_sample")

;; execute state-machine and publish every step
(defun smach-exec (sm)
  (let ((insp (instance state-machine-inspector :init sm))
	(mydata (list (cons 'count 0) (cons 'hoge (list 1 2 3)))))
    (unix:sleep 1) ;; for advertize topics successfully
    (send sm :reset-state)
    (send insp :publish-structure)
    (send insp :publish-status mydata)
    (while (not (send sm :goal-reached))
      (unix:sleep 1)
      (send insp :publish-status mydata)
      ;; you can execute state-machine step by step
      (send sm :execute mydata :step -1))
    (send sm :active-state) ;; return reached goal state
    ))

(defun smach-exec-simple () (setq count 0) (smach-exec (smach-simple)))
(defun smach-exec-nested () (setq count 0) (smach-exec (smach-nested)))
(defun smach-exec-userdata () (smach-exec (smach-userdata)))

(warn ";;(smach-exec-simple)~%;;(smach-exec-nested)~%;;(smach-exec-userdata)~%(smach-exec-actionlib)~%")

(ros::roseus-add-msgs "actionlib_tutorials")
(defun func-fibonacci (userdata)
  (let* ((goal (cdr (assoc :goal userdata)))
	 (result (cdr (assoc :result userdata)))
	 (feedback (cdr (assoc :feedback userdata)))
	 (seq (send result :sequence)) (next 0))
    (when (/= (length seq) (send goal :order))
      (setq seq (make-array (send goal :order) :element-type :integer))
      (send result :sequence seq)
      (setf (elt seq 0) 1)
      (setf (elt seq 1) 1))
    (send feedback :sequence seq)
    (while t
      (cond ((eq next (send goal :order))
	     (return-from func-fibonacci t))
	    ((zerop (elt seq next))
	     (setf (elt seq next) (+ (elt seq (- next 1))
				     (elt seq (- next 2))))
	     (return-from func-fibonacci nil))
	    (t (incf next))))
    ))
(defun func-sleep1 (&rest args) (unix::sleep 1) t)
(defun func-sleep2 (&rest args) (unix::sleep 1) t)

(defun fibonacci-smach ()
  (let ((sm (instance state-machine :init)))
    ;; remap from variable names in function to names in userdata
    (send sm :add-node (instance state :init :CALC 'func-fibonacci))
    (send sm :add-node (instance state :init :SLEEP1 'func-sleep1))
    (send sm :add-node (instance state :init :SLEEP2 'func-sleep2))
    (send sm :goal-state (list :success))
    ;; set remap list later

    (send sm :start-state :CALC)
    (send sm :add-transition :CALC :SLEEP1 nil)
    (send sm :add-transition :CALC :success t)
    (send sm :add-transition :SLEEP1 :SLEEP2 t)
    (send sm :add-transition :SLEEP2 :CALC t)
    sm ))

(defun smach-exec-actionlib ()
  (warn "fibonacci server is working and press any key to stop.")
  (setq sm (fibonacci-smach))
  (setq s (instance state-machine-action-server :init "/fibonacci"
		    actionlib_tutorials::fibonacciaction sm))
  (do-until-key (ros::spin-once) (send s :worker) (unix::sleep 1)))

