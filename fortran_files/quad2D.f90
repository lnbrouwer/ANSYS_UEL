!___________________________________________________________________________
!     --- Module for 2D Quadratic Elements (8 nodes, 4 int)  (quad2D.f90)
!     --- sets shape functions and their gradients in local sys
!     --- same for all elem of this type (=> called 1x at start) 
!     --- for now ASSUMES QUADRILATERAL ELEMENTS ONLY (OR DEGEN)
!     --- distribution R1.0, Lucas Brouwer 11/15/2018
!___________________________________________________________________________

MODULE quad2D

IMPLICIT NONE

PRIVATE
!     element characteristics		
  INTEGER nnode,edim,nq1D,etype,nq	
	PARAMETER (nnode=8)      
	PARAMETER (edim=2)          ! 2D element
	PARAMETER (nq1D=2)          ! # guassian point integration (in 1D)
	PARAMETER (etype=2)         ! 1 = linear shape functions, 2 = quadratic shape functions (midside)
	PARAMETER (nq=4)            ! nq1D^edim
	
	
!     quadrature and shape functions (set 1x)		
	DOUBLE PRECISION qw(nq),dH(nnode,edim,nq)
	DOUBLE PRECISION qpt(edim,nq),HH(nq,nq),H(nnode,nq)

!     SUBROUTINEs		
	PUBLIC set_dH,set_quad,set_HH
	PUBLIC qw,dH,HH,H
	PUBLIC nnode,edim,nq
	PUBLIC INVF2
	
CONTAINS

!___________________________________________________________________________
!     --- shape function gradients w.r.t local var (at quad points)
!___________________________________________________________________________

	SUBROUTINE set_dH()
		DOUBLE PRECISION r,s,t	
		INTEGER i
		
!	set dH based on shape function gradients w.r.t local var (at quad points)

		DO i = 1,nq
			t = qpt(1,i)
			s = qpt(2,i)
			
!           2D 4 node element (linear)
! 					dh(1,1,i) = +0.25d0*(1.d0+s)
! 					dh(2,1,i) = +0.25d0*(-1.d0-s)
! 					dh(3,1,i) = +0.25d0*(-1.d0+s)
! 					dh(4,1,i) = +0.25d0*(1.d0-s)
! 					
! 					dh(1,2,i) = +0.25d0*(1.d0+r)
! 					dh(2,2,i) = +0.25d0*(1.d0-r)
! 					dh(3,2,i) = +0.25d0*(-1.d0+r)
! 					dh(4,2,i) = +0.25d0*(-1.d0-r)

!           2D 8 node element (quadratic)
			H(1,i) = +0.25d0*(+1.d0-s)*(+1.d0-t)*(-s-t-1.d0)
			H(2,i) = +0.25d0*(+1.d0+s)*(+1.d0-t)*(+s-t-1.d0)					
			H(3,i) = +0.25d0*(+1.d0+s)*(+1.d0+t)*(+s+t-1.d0)					
			H(4,i) = +0.25d0*(+1.d0-s)*(+1.d0+t)*(-s+t-1.d0)
			H(5,i) = +0.5d0*(+1.d0-s*s)*(+1.d0-t)
			H(6,i) = +0.5d0*(+1.d0+s)*(+1.d0-t*t)
			H(7,i) = +0.5d0*(+1.d0-s*s)*(+1.d0+t)
			H(8,i) = +0.5d0*(+1.d0-s)*(+1.d0-t*t)
			
!     2D 8 node element (quadratic)
			dH(1,1,i) = -0.25d0*(-1.d0+t)*(2.d0*s+t)
			dH(2,1,i) = -0.25d0*(-1.d0+t)*(2.d0*s-t)
			dH(3,1,i) = +0.25d0*(+1.d0+t)*(2.d0*s+t)	
			dH(4,1,i) = +0.25d0*(+1.d0+t)*(2.d0*s-t)
			dH(5,1,i) = s*(-1.d0+t)
			dH(6,1,i) = +0.5d0*(+1.d0-t*t)
			dH(7,1,i) = -s*(+1.d0+t)
			dH(8,1,i) = +0.5d0*(-1.d0+t*t)
			
			dH(1,2,i) = -0.25d0*(-1.d0+s)*(+2.d0*t+s)
			dH(2,2,i) = -0.25d0*(+1.d0+s)*(-2.d0*t+s)
			dH(3,2,i) = +0.25d0*(+1.d0+s)*(+2.d0*t+s)	
			dH(4,2,i) = +0.25d0*(-1.d0+s)*(-2.d0*t+s)
			dH(5,2,i) = +0.5d0*(-1.d0+s*s)
			dH(6,2,i) = -t*(+1.d0+s)
			dH(7,2,i) = +0.5d0*(+1.d0-s*s)
			dH(8,2,i) = t*(-1.d0+s)
		ENDDO						

		RETURN
	END SUBROUTINE set_dH

!___________________________________________________________________________
!     --- set quadrature points and weights
!___________________________________________________________________________

	SUBROUTINE set_quad()
		INTEGER i
		INTEGER permr(4),perms(4)	
		DOUBLE PRECISION P(nq1D),W(nq1D)			
! set up gaussian point locations and weights (1D)			
		W(1) = 1.0d0
		W(2) = 1.0d0
		P(1) = dsqrt(1.0d0/3.0d0)
		P(2) = -dsqrt(1.0d0/3.0d0)

!   permute to get
!    4   3
!    1   2
		permr = (/1,2,2,1/)    	
		perms = (/1,1,2,2/)    
	  DO i=1,nq
				qpt(1,i) = P(permr(i))
				qpt(2,i) = P(perms(i))
				qw(i) = W(permr(i))*W(perms(i))						
	  ENDDO
		
		RETURN
	END SUBROUTINE set_quad
	
!___________________________________________________________________________
!     --- extrapolation of solution var from quad points to corner nodes
!___________________________________________________________________________

	SUBROUTINE set_HH()
		DOUBLE PRECISION S3
!	set HH which extrapolates scalars from guass points to corners	(assumes 4 quad points)
		S3 = dsqrt(3.0d0)
		
		HH(1,1) = 1.0d0+0.5d0*S3
		HH(1,2) = -0.5d0			
		HH(1,3) = 1.0d0-0.5d0*S3
		HH(1,4) = -0.5d0		
		
		HH(2,1) = -0.5d0		
		HH(2,2) = 1.0d0+0.5d0*S3				
		HH(2,3) = -0.5d0					
		HH(2,4) = 1.0d0-0.5d0*S3
		
		HH(3,1) = 1.0d0-0.5d0*S3
		HH(3,2) = -0.5d0			
		HH(3,3) = 1.0d0+0.5d0*S3
		HH(3,4) = -0.5d0	
		
		HH(4,1) = -0.5d0		
		HH(4,2) = 1.0d0-0.5d0*S3				
		HH(4,3) = -0.5d0					
		HH(4,4) = 1.0d0+0.5d0*S3	
					
		RETURN
	END SUBROUTINE set_HH	
	
!____________________________________________________________________
!	 Subroutine INVF2 - from Diego
!____________________________________________________________________
  SUBROUTINE INVF2(D,F,FINV,Jac)
    INTEGER D
    DOUBLE PRECISION F(D,D),FINV(D,D),Jac
		IF (D == 2) THEN
			Jac = F(1,1)*F(2,2) - F(1,2)*F(2,1)
			FINV(1,1) = F(2,2)/Jac
			FINV(2,2) = F(1,1)/Jac
			FINV(2,1) = -F(2,1)/Jac
			FINV(1,2) = -F(1,2)/Jac
		ELSE IF (D == 3) THEN
			Jac = F(1,1)*(F(2,2)*F(3,3) - F(2,3)*F(3,2)) - F(1,2)*(F(2,1)*F(3,3) - F(3,1)*F(2,3)) +	F(1,3)*(F(2,1)*F(3,2) - F(3,1)*F(2,2))
			FINV(1,1) =  (F(2,2)*F(3,3) - F(2,3)*F(3,2))/Jac
			FINV(2,1) = -(F(2,1)*F(3,3) - F(3,1)*F(2,3))/Jac
			FINV(3,1) =  (F(2,1)*F(3,2) - F(3,1)*F(2,2))/Jac   
			FINV(1,2) = -(F(1,2)*F(3,3) - F(1,3)*F(3,2))/Jac
			FINV(2,2) =  (F(1,1)*F(3,3) - F(1,3)*F(3,1))/Jac
			FINV(3,2) = -(F(1,1)*F(3,2) - F(1,2)*F(3,1))/Jac  
			FINV(1,3) =  (F(1,2)*F(2,3) - F(1,3)*F(2,2))/Jac
			FINV(2,3) = -(F(1,1)*F(2,3) - F(2,1)*F(1,3))/Jac
			FINV(3,3) =  (F(1,1)*F(2,2) - F(2,1)*F(1,2))/Jac  
     ENDIF
     RETURN
  END SUBROUTINE INVF2
			
END MODULE quad2D




