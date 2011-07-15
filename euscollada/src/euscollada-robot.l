;; This euscollada-robot class is for robots converted from collada files.
;; This class provides :init-ending override.
(defclass euscollada-robot
  :super robot-model
  :slots ()
  )
(defmethod euscollada-robot
  (:init-ending
   ()
   (send-super :init-ending)
   ;; update link mass properties
   (labels ((find-parent-link-in-links
  	     (ll)
  	     (if (not (send ll :parent-link))
  		 (let ((tmp (send ll :parent)))
  		   (while (not (and (derivedp tmp bodyset-link)
  				    (find tmp (send self :links))))
  		     (setq tmp (send tmp :parent)))
  		   tmp)))
  	    (gather-all-links
  	     (&optional (ll (car (send self :links))))
  	     (append (list ll)
  		     (flatten
  		      (mapcar #'(lambda (a)
  				  (gather-all-links a))
  			      (remove-if-not #'(lambda (x) (derivedp x bodyset-link)) (send ll :descendants)))))
  	     ))
     ;; gather all links not included in (send self :links)
     (send-all (send self :links) :put :tmp-child-links nil)
     (dolist (ll (cdr (gather-all-links))) ;; but root link
       (unless (send ll :parent-link)
  	 (let ((pl (find-parent-link-in-links ll)))
  	   (send pl :put :tmp-child-links
  		 (cons ll (send pl :get :tmp-child-links)))
  	   )))
     ;; append mass properties of gathered links to their parent-links
     (dolist (ll (remove-if-not
  		  #'(lambda (x) (send x :get :tmp-child-links))
  		  (send self :links)))
       (send ll :append-mass-properties (send ll :get :tmp-child-links) :update t))
     (dolist (ll (send self :links)) (remprop ll :tmp-child-links)))

   ;; check root link validation
   (labels
       ((get-parent-link
	 (ll)
	 (or (send ll :parent-link)
	     (if (derivedp (send ll :parent) bodyset-link)
		 (send ll :parent))))
	(find-root-link
	 (ll)
	 (let ((pl (get-parent-link ll)))
	   (cond
	    ((null pl) nil)
	    ((get-parent-link pl) (find-root-link pl))
	    ((find (send ll :joint) (send self :joint-list)) pl)
	    (t ll)))))
     (let ((root-link
	    (remove-duplicates
	     (remove nil
		     (mapcar #'(lambda (l)
				 (find-root-link l))
			     (send-all (send self :joint-list) :parent-link))))))
       (when (> (length (send self :links)) 1)
	 (unless (= (length root-link) 1)
	   (error "root link definition is ambiguous!!"))
	 (unless (equal (car root-link) (car links))
	   (setq links (append root-link (cdr links))))))) ;; replace root link
   )
  )
