!___________________________________________________________________________
!     --- MODULE for quench behavior (quench.f90)
!     --- defines Jc(B,T)
!     --- computes Tcs(B,J),TcB0(B) IF desired
!     --- distribution R1.0, Lucas Brouwer 11/15/2018
!___________________________________________________________________________


MODULE quench

IMPLICIT NONE

PRIVATE


!     SUBROUTINEs		
  PUBLIC check_quench
  PRIVATE secant,JcSCU,bisect
	
CONTAINS

!___________________________________________________________________________
!     --- Jc fit for SCU (from Diego)
!___________________________________________________________________________

	SUBROUTINE JcSCU(Temp,Bmag,Jc,TcB0)
		DOUBLE PRECISION Temp,Bmag,Jc,Bc
		DOUBLE PRECISION Jc0,Tc0,Bc0,TcB0	
		DOUBLE PRECISION t,h		


! RRP wire 0.6 mm diameter -> Diego's fit for SCU
    Jc0 = 6600;
    Tc0 = 16.7;
    Bc0 = 26.25;
		
!	Calc Jc for given B,T
    t = temp/Tc0
		Bc = Bc0*(1d0-t**1.52d0)
    h = Bmag/Bc
		TcB0 = Tc0*((1d0-Bmag/Bc0)**(1d0/1.52d0))
		IF (Bmag .lt. Bc) THEN
		  Jc = Jc0*(1d0-t**2d0)*h**(-0.5d0)*(1d0-h)**2
!		  Jc = Jc0*(1d0-(temp/Tc0)**2d0)*(Bmag/(Bc0*(1d0-(temp/Tc0)**1.52d0)))**(-0.5d0)*(1d0-(Bmag/(Bc0*(1d0-(temp/Tc0)**1.52d0))))**2
    ELSE
		  Jc = 0d0
    ENDIF	
		
		RETURN
	END SUBROUTINE JcSCU
	
	
!___________________________________________________________________________
!     --- determines quench state
!___________________________________________________________________________

	SUBROUTINE check_quench(Jcw,Temp,Bmag,Jc,qflag,IfCu,Jc0,Tc0,Bc0,alpha,p)
	
		DOUBLE PRECISION Jcw,Temp,Bmag
		DOUBLE PRECISION Jc,IfCu
		INTEGER qflag

		DOUBLE PRECISION Jc0,Tc0,Bc0,Bc,Jc0m	
		DOUBLE PRECISION t,h,alpha,p,TcB0,mm2msq


! RRP wire 0.6 mm diameter -> Diego's fit for SCU
!    Jc0 = 6600.0d0
!    Tc0 = 16.7d0
!    Bc0 = 26.25d0
!    alpha = 1.0d0	
!    p = 1.52d0		
		mm2msq = 1.0d6		
		Jc0m = Jc0*mm2msq

		
    IF ((temp .GE. Tc0) .OR. (Jcw .GE. Jc0m) .OR. (Bmag .GE. Bc0))  THEN
      qflag = 1
      IfCu = 1.0d0
      Jc = -1.0d0
		ELSE
!	    Calc Jc for given B,T
      t = temp/Tc0
		  Bc = Bc0*(1d0-t**p)
      h = Bmag/Bc
		  TcB0 = Tc0*((1d0-Bmag/Bc0)**(1d0/p))
  		IF (Bmag .lt. Bc) THEN
!        Jc = mm2msq*Jc0*(1d0-t**2d0)*h**(-0.5d0)*(1d0-h)**2d0
        Jc = mm2msq*Jc0*(1d0-t**2d0)**alpha*(1d0-t**p)**(alpha-1d0)*h**(-0.5d0)*(1d0-h)**2d0
!		    Jc = Jc0*(1d0-(temp/Tc0)**2d0)*(Bmag/(Bc0*(1d0-(temp/Tc0)**1.52d0)))**(-0.5d0)*(1d0-(Bmag/(Bc0*(1d0-(temp/Tc0)**1.52d0))))**2
        IF (Jc .ge. Jc0m) THEN
				  Jc = Jc0m
        ENDIF								
      ELSE
		    Jc = 0d0
      ENDIF	
		
      IF (temp .GE. TcB0)  THEN	
        qflag = 1
        IfCu = 1.0d0	
			ELSE IF (Jcw .LE. Jc)  THEN
        qflag = -1
        IfCu = 0.0d0			
      ELSE
        qflag = 0
        IfCu = 1.0d0-Jc/Jcw			
      ENDIF		 
		
    ENDIF	
		
		RETURN
	END SUBROUTINE check_quench


!___________________________________________________________________________
!     --- secant method for finding Tcs -> root search
!___________________________________________________________________________

	SUBROUTINE secant(Temp,Bmag,Jcw,Tcs,maxit,tol,i)
	  DOUBLE PRECISION Tcs,tol,x0,x1,y0,y1,dx,dy,x2
		DOUBLE PRECISION Temp,Bmag,Jcw,TcB0,mm2msq,Tc0
		INTEGER maxit
		INTEGER i,j
		mm2msq = 1.0d6
		
!   initial guesses -> probably dont want to overshoot b/c of multiple roots
!   bisection may be better to avoid multiple roots 
!   could also make dx very small to hopefully avoid jumping around
    i = 0
    x0 = 1.9d0
    x1 = 2.0d0
		dx = x1-x0		
		CALL JcSCU(x0,Bmag,y0,TcB0)
		y0 = y0-Jcw/mm2msq

				
		DO WHILE (ABS(dx).GT.tol)		
		  CALL JcSCU(x1,Bmag,y1,TcB0)	
		  y1 = y1-Jcw/mm2msq
			x2 = x1-y1*(x1-x0)/(y1-y0)
			x0 = x1
			x1 = x2
			dx = x1-x0
			y0 = y1

      i = i+1			
		  IF (i.eq.maxit) THEN 
		    EXIT		
      END IF				
		END DO
	
		Tcs = x0
		
		RETURN
	END SUBROUTINE secant
	
!___________________________________________________________________________
!     --- bisection method for finding Tcs -> root search
!___________________________________________________________________________

	SUBROUTINE bisect(Temp,Bmag,Jcw,Tcs,maxit,tol,i)
	  DOUBLE PRECISION Tcs,tol,x0,x1,y0,y1,dx,dy,x2
		DOUBLE PRECISION Temp,Bmag,Jcw,TcB0,mm2msq
		INTEGER maxit
		INTEGER i,j
		DOUBLE PRECISION a,b,c
				
		mm2msq = 1.0d6
		
!    search in 0.2 K ranges up to Tcs 
!      -> hopefully never dips below and back up in this size		
!      -> could be improved with safety and exit conditions
!         and by capping search range at TcB0

		DO i = 1,78,1
			x0 = 1d0+i*.2d0
			x1 = 1d0+(i+1)*.2d0		
			
			CALL JcSCU(x0,Bmag,y0,TcB0)
			y0 = y0-Jcw/mm2msq			
		  CALL JcSCU(x1,Bmag,y1,TcB0)	
		  y1 = y1-Jcw/mm2msq
			
			IF (y1*y0 .lt. 0) THEN     
				a = x0
				b = x1
				EXIT	
			END IF 		
		END DO
		
		i = 1
		DO WHILE ((b-a)/2.0d0 .GT. tol)		
			c = (a+b)/2.0d0

			CALL JcSCU(a,Bmag,y0,TcB0)
			y0 = y0-Jcw/mm2msq			
		  CALL JcSCU(c,Bmag,y1,TcB0)	
		  y1 = y1-Jcw/mm2msq

			IF (y1*y0 .lt. 0) THEN     
				b=c
			ELSE
			  a=c
			END IF 	
			
      i = i+1			
		  IF (i.eq.maxit) THEN 
		    EXIT		
      END IF	
			
		END DO
		Tcs = (a+b)/2.0d0
		
		RETURN
	END SUBROUTINE bisect
	
END MODULE quench




