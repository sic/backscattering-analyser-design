! Program to generate the parameters for an analyser bank (IRIS/OSIRIS)
! Ellipse with focii at the sample and detector  
  PROGRAM ellipse
  IMPLICIT NONE
  INTEGER :: n, i   
  REAL :: xd, yd, xo, yo, ystart, sa, phie,s,t,tmpphi   
  REAL :: A,B,C,D,E,F,G,P,Q,R,qa,qb,qc,ar,br,cp,sp   
  REAL :: const_A, const_B, result, ar2   
  REAL, ALLOCATABLE :: x(:), y(:), phi(:)   
  CHARACTER(len=80) :: filename	
  
  WRITE(6,*) ' How many analyser crystals ?'	
  READ(5,*) n		ALLOCATE(x(n),y(n),phi(n))	
  WRITE(6,*) ' Output filename ?'	
  READ(5,'(A80)') filename      
  OPEN(unit=50,file=filename,status='new')	
  WRITE(6,*) ' Coordinates of detector ?'	
  READ(5,*) xd, yd	WRITE(6,*) ' Maximum y-coord of analyser bank ?'	
  READ(5,*) ystart	WRITE(6,*) ' Distance from sample to analyser at y=0 ?'	
  READ(5,*) sa! Calculate centre of ellipse.	
  
  xo = xd / 2.0	
  yo = yd / 2.0        	
  
  s = xo	
  t = xo	
  phie = atan(yd/xd)
  WRITE(6,*) ' Centre of ellipse @ ',xo,yo
  WRITE(6,*) ' Orientation angle = ',phie/0.01745329251	                                        
  
  cp = cos(phie)	
  sp = sin(phie) 

! ar,br are the radii in the x,y directions (axes local to ellipse)
! Need to find them, we know that at y=0, x=sa	
  const_A = (sa*sa*cp*cp)-(2*sa*s*cp*cp)-(2*sa*t*cp*sp)+ 
&		    (s*s*cp*cp)+(t*t*sp*sp)+(2*s*t*cp*sp)	const_B = (sa*sa*sp*sp)-(2*sa*s*sp*sp)-(2*sa*t*cp*sp)+ &		    (s*s*sp*sp)+(t*t*cp*cp)-(2*s*t*cp*sp)			br = sa/2.0	ar = sqrt((sa-(abs(xd)))**2 + (yd)**2)	WRITE(*,*) ' Secondary flight path = ',sa+ar	ar = (ar+sa-sqrt(xd**2 + yd**2))/2 +(0.5*sqrt(xd**2 + yd**2))	br = sqrt(const_B/(1-(const_A/(ar*ar))))!	br = sa-(sa/10)!      DO     ! loop over values for one of the axes of the ellipse!	br = br + 0.1! Define terms in equation Ax^2+Bx+Cy^2+Dy+Ez^2+Fz+G+Pxy+Qyz+Rxz=0	A = (cp*cp)/(ar*ar) + (sp*sp)/(br*br)	B = ((-2*s*cp*cp)/(ar*ar)) + ((-2*t*cp*sp)/(ar*ar)) + &	    ((-2*s*sp*sp)/(br*br)) + ((-2*t*cp*sp)/(br*br))	C = (sp*sp)/(ar*ar) + (cp*cp)/(br*br)	D = ((-2*t*sp*sp)/(ar*ar)) + ((-2*t*cp*cp)/(br*br)) + &	    ((-2*s*cp*sp)/(ar*ar)) + ((-2*s*cp*sp)/(br*br))	E = 0	F = 0	G = ((s*s*cp*cp)/(ar*ar)) + ((t*t*sp*sp)/(ar*ar)) + &	    ((s*s*sp*sp)/(br*br)) + ((t*t*cp*cp)/(br*br)) + &	    ((2*s*t*cp*sp)/(ar*ar)) + ((-2*s*t*cp*sp)/(br*br)) - 1	P = (2*cp*sp)/(ar*ar) + (-2*cp*sp)/(br*br)	Q = 0	R = 0	result = (A*sa**2) + (B*sa) + G	WRITE(*,*) result	!	IF (abs(result)<1E-4) EXIT!      END DO	WRITE(*,*) ' Axes of ellipse are ',ar,br! if we know Y then we have an equation of the form qa*X^2 + qb*X + qc = 0      DO i = 1, n 	  y(i) = ystart - (i-1) - 0.5	  qa = A	  qb = B + (P*y(i))	  qc = (C*y(i)**2) + (D*y(i)) + G	  x(i) = (-qb+sqrt(abs(qb**2 - (4*qa*qc))))/(2*qa)!	  x(i) = (-qb-sqrt(qb**2 - (4*qa*qc)))/(2*qa)        tmpphi = (atan((-br*br/ar*ar)*(x(i)-s)/(y(i)-t)))/0.0174532925	  phi(i) = - (90.0 - tmpphi) !	  write(6,100) x(i), y(i), phi(i)!	  write(50,100) x(i), y(i), phi(i)	  write(6,*) x(i), y(i), phi(i)	  write(50,*) x(i), y(i), phi(i)	ENDDO      100	FORMAT(1X,F7.3,1X,F7.3,1X,F7.3) 	CLOSE(50)  END PROGRAM ellipse                    
