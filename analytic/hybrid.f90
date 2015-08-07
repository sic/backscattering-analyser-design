! Program to generate the parameters for an analyser bank (IRIS/OSIRIS)
! Ellipse with focii at the sample and detector  
    PROGRAM hybrid
    IMPLICIT NONE
    INTEGER :: n, i, count, unit   
    REAL :: xd, yd, xo, yo, ystart, sa, phie,s,t, path   
    REAL :: lambda, lambdac   
    REAL :: A,B,C,D,E,F,G,P,Q,R,qa,qb,qc,ar,br,cp,sp, tmp   
    REAL :: const_A, const_B,result,pathc,ttheta,gradient   
    REAL, ALLOCATABLE :: x(:), y(:), phi(:)   
    CHARACTER(len=80) :: filename	
    
    WRITE(6,*) ' How many analyser crystals ?'	
    READ(5,*) n
    ALLOCATE(x(n),y(n),phi(n))

    WRITE(6,*) ' Output filename ?'	
    READ(5,'(A80)') filename      
    OPEN(unit=50,file=filename,status='new')      
    
    OPEN(unit=51,file='masked.ap',status='new')	
    WRITE(6,*) ' Coordinates of detector ?'	
    READ(5,*) xd, yd	
    WRITE(6,*) ' Maximum y-coord of analyser bank ?'	
    READ(5,*) ystart	
    WRITE(6,*) ' Distance from sample to analyser at y=0 ?'	
    READ(5,*) sa	
    WRITE(6,*) ' Value of constant 2theta ?'	
    READ(5,*) ttheta
    !	ttheta = ttheta*0.01745329251
    ! Calculate centre of ellipse.
    xo = xd / 2.0
    yo = yd / 2.0
    s = xo
    t = xo
    phie = atan(yd/xd)	
    
    WRITE(6,*) ' Centre of ellipse @ ',xo,yo
    WRITE(6,*) ' Orientation angle = ',phie/0.01745329251,' degrees'	      
    
    cp = cos(phie)
	sp = sin(phie)

! ar,br are the radii in the x,y directions (axes local to ellipse)
! Need to find them, we know that at y=0, x=sa
!	const_A = (sa*sa*cp*cp)-(2*sa*s*cp*cp)-(2*sa*t*cp*sp)+ &
!		    (s*s*cp*cp)+(t*t*sp*sp)+(2*s*t*cp*sp)
!	const_B = (sa*sa*sp*sp)-(2*sa*s*sp*sp)-(2*sa*t*cp*sp)+ &
!		    (s*s*sp*sp)+(t*t*cp*cp)-(2*s*t*cp*sp)
!			ar = sqrt((sa-(abs(xd)))**2 + (yd)**2)
	pathc = ar+sa
		WRITE(6,*) ' Design Secondary flight path = ',pathc
!	ar = (ar+sa-sqrt(xd**2 + yd**2))/2 +(0.5*sqrt(xd**2 + yd**2))
!	br = sqrt(const_B/(1-(const_A/(ar*ar))))
! Define terms in equation Ax^2+Bx+Cy^2+Dy+Ez^2+Fz+G+Pxy+Qyz+Rxz=0
!	A = (cp*cp)/(ar*ar) + (sp*sp)/(br*br)
!	B = ((-2*s*cp*cp)/(ar*ar)) + ((-2*t*cp*sp)/(ar*ar)) + &
!	    ((-2*s*sp*sp)/(br*br)) + ((-2*t*cp*sp)/(br*br))
!	C = (sp*sp)/(ar*ar) + (cp*cp)/(br*br)
!	D = ((-2*t*sp*sp)/(ar*ar)) + ((-2*t*cp*cp)/(br*br)) + &
!	    ((-2*s*cp*sp)/(ar*ar)) + ((-2*s*cp*sp)/(br*br))
!	E = 0
!	F = 0
!	G = ((s*s*cp*cp)/(ar*ar)) + ((t*t*sp*sp)/(ar*ar)) + &
!	    ((s*s*sp*sp)/(br*br)) + ((t*t*cp*cp)/(br*br)) + &
!	    ((2*s*t*cp*sp)/(ar*ar)) + ((-2*s*t*cp*sp)/(br*br)) - 1
!	P = (2*cp*sp)/(ar*ar) + (-2*cp*sp)/(br*br)
!	Q = 0!	R = 0
!
!	result = (A*sa**2) + (B*sa) + G
!	WRITE(*,*) ' This figure should be close to zero --> ',result	                            
!	WRITE(*,*) ' Axes of ellipse are ',ar,br	x(:) = sa-(sa/5)    

! if we know Y then we have an equation of the form qa*X^2 + qb*X + qc = 0
      DO i = 1, n
 	  y(i) = ystart - (i-1) - 0.5	  
 	  count = 0	  
 	  unit = 50
 	  !	  qa = A
 	  !	  qb = B + (P*y(i))
 	  !	  qc = (C*y(i)**2) + (D*y(i)) + G
 	  DO	   
 	      x(i) = x(i) + 0.001
 	      count = count + 1	   
 	      path = sqrt(x(i)**2 + y(i)**2) + & 		   
 	             sqrt((x(i)-xd)**2 + (y(i)-yd)**2)	  
 	      !WRITE(*,*) x(i), count, path 	   
 	      IF (abs(path-pathc)<0.001) EXIT	  
 	  END DO	  
 	  gradient = ((yd+2.5)/xd)*x(i)
 	  write(6,*) '====================================='
 	  IF (gradient > y(i)) THEN		WRITE(6,*) ' WARNING - Analyser Element &		&is being shielded by detector'	  unit = 51	  ENDIF!	  DO 
 	  !         lambdac = 6.69 * sin((3.14159265359/2.0) - phi)
 	  ! 	   IF (abs(ttheta-tthetac)<0.001) 
 	  EXIT!        END DO	  write(6,*) ' Analyser Crystal No : ',i!SIC	  phi(i) =  -atan((y(i)-yo)/(x(i)-xo))/0.01745329251		  tmp = atan(y(i)/x(i))/0.01745329251        phi(i) = -(180.0 - (90.0 - tmp) - ttheta/2.0)		  write(6,100) x(i), y(i), phi(i)	  write(6,*) ' Analysed wavelength = ',6.69*sin(tmp)	  write(6,*) ' Secondary Path Length = ',path	  write(unit,100) x(i), y(i), phi(i)	ENDDO      100	FORMAT(1X,F7.3,1X,F7.3,1X,F7.3) 	CLOSE(50)	
 	  CLOSE(51)  
 	  END PROGRAM hybrid
