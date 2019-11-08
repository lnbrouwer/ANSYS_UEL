!___________________________________________________________________________
!     --- material property libaray (matLIB.f90)
!     --- most fits can be found in G. Manfreda, Review of ROXIE's Material Properties Database for Quench Simulation,
!     --- CERN Internal Note: EDMS NR, vol. 1178007, 2018.
!     --- distribution R1.0, Lucas Brouwer 11/15/2018
!___________________________________________________________________________


MODULE matLIB

IMPLICIT NONE

PRIVATE

!  functions


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! copper properties 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! CUDI
	PUBLIC cvcucudi,rhocucudi,kxxcudi
	
! NIST
	PUBLIC cvcunist,rhocunist,kxxnist
	PUBLIC rhocunist_tref,kxxnist_tref
	PRIVATE kxxnist_B0
	
! MATPRO
	PUBLIC cvcumatpro,kxxmatpro,rhocumatpro
	

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Nb3Sn properties 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	

! CUDI    
	PUBLIC cvnb3sncudi

! NIST 
	PUBLIC cvnb3snnist
	

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! G10 properties 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	

! NIST    
	PUBLIC cvG10nist
	

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  SUBROUTINEs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	
!  LININT is for matpro
  PRIVATE LININT	

	    
CONTAINS




!________________________________________________________________
!                                                               
!                                                               
!                  COPPER MATERIAL PROPERTIES                   
!                                                               
!                                                               
!________________________________________________________________


!________________________________________________________________                                                         
!                           CUDI                                
!________________________________________________________________ 

FUNCTION cvcucudi(TT)

	DOUBLE PRECISION, INTENT(IN) :: TT
	DOUBLE PRECISION             :: cvcucudi

  IF(TT < 9.441D0) THEN
		cvcucudi = -0.0308 * TT**4.0D0 + 7.229 * TT**3.0D0 - 2.1286 * TT**2.0D0 + 101.89 * TT + 2.5631
	ELSE
	  IF(TT < 31.134D0) THEN
			cvcucudi = -0.3045 * TT**4.0D0 + 29.871 * TT**3.0D0 - 455.61 * TT**2.0D0 + 3469.5 * TT - 8250.3
		ELSE
			IF(TT < 123.34D0) THEN
				cvcucudi = 0.0419 * TT**4.0D0 - 14.024 * TT**3.0D0 + 1508.9 * TT**2.0D0 - 31595 * TT + 178432
			ELSE
				IF(TT < 306.12D0) THEN
					cvcucudi = -8.48E-4 * TT**4.0D0 + 0.8419 * TT**3.0D0 - 325.52 * TT**2.0D0 + 60590 * TT - 1.2851E6
				ELSE
					IF(TT < 498.15D0) THEN
						cvcucudi = -4.80E-5 * TT**4.0D0 + 0.09173 * TT**3.0D0 - 64.12 * TT**2.0D0 + 20363 * TT + 1.028E6
					ELSE
						cvcucudi = 12E-5 * TT**3.0D0 - 0.21486 * TT**2.0D0 + 1003.84 * TT + 3.1823E6
					ENDIF
				END IF
			END IF
		END IF
	END IF
	
END FUNCTION cvcucudi


	

FUNCTION kxxcudi(TT,RRR,BB)

  DOUBLE PRECISION TT, RRR, BB, kxxcudi
	kxxcudi = kxxnist_B0(TT,RRR)*rhocucudi(TT,RRR,0.0D0)/rhocucudi(TT,RRR,BB)

END FUNCTION kxxcudi



FUNCTION rhocucudi(TT,RRR,BB)

	DOUBLE PRECISION TT, RRR, BB, rhocucudi
	DOUBLE PRECISION INV_T, INV_T2
	INV_T = 1.0/TT
	INV_T2 = INV_T*INV_T
	rhocucudi=1.D-8*(1.7/RRR+1/(((2.32547E9*INV_T2 + 9.57137E5)*INV_T2 + 1.62735E2)*INV_T) + 0.005*BB)
				
END FUNCTION rhocucudi



!________________________________________________________________                                                         
!                           NIST                                
!________________________________________________________________ 



FUNCTION cvcunist(TT)

	DOUBLE PRECISION :: dc_a, dc_b, dc_c, dc_d, dc_e, dc_f, dc_g, dc_h
	DOUBLE PRECISION :: density
	DOUBLE PRECISION, intent(in) :: TT 
	DOUBLE PRECISION :: cvcunist, p

	density = 8960 

	IF (TT.LE.300) THEN

		dc_a = -1.91844
		dc_b = -0.15973 
		dc_c = 8.61013
		dc_d =-18.996 
		dc_e = 21.9661
		dc_f =-12.7328
		dc_g =  3.54322
		dc_h = -0.3797 
		
		p =   dc_a                         + dc_b * (LOG10(TT))**1.0D0 + dc_c * (LOG10(TT))**2.0D0 &
	+ dc_d * (LOG10(TT))**3.0D0 + dc_e * (LOG10(TT))**4.0D0 + dc_f * (LOG10(TT))**5.0D0 &
	+ dc_g * (LOG10(TT))**6.0D0 + dc_h * (LOG10(TT))**7.0D0
					p = 10**p

	ELSE 
					p = 361.5D0+0.093D0*TT
	ENDIF

	cvcunist = density*p	  

END FUNCTION cvcunist


FUNCTION kxxnist_B0(TT,RRR)

	DOUBLE PRECISION TT, RRR, kxxnist_B0
	DOUBLE PRECISION beta, k_0, k_i, k_i0
	DOUBLE PRECISION P1, P2, P3, P4, P5, P6, P7

	beta = 0.634/RRR
	P1 = 1.754E-8
	P2 = 2.763
	P3 = 1102
	P4 = -0.165
	P5 = 70
	P6 = 1.756
	P7 = 0.838/(beta/0.0003)**0.1661

	k_0 = beta/TT
	k_i = P1*TT**P2/(1+P1*P3*TT**(P2+P4)*exp(-(P5/TT)**P6))
	k_i0 = P7*k_i*k_0/(k_i+k_0)
	kxxnist_B0 = 1/(k_0+k_i+k_i0)

END FUNCTION kxxnist_B0


FUNCTION kxxnist(TT,RRR,BB)

  DOUBLE PRECISION TT, RRR, BB, kxxnist
	kxxnist = kxxnist_B0(TT,RRR)*rhocunist(TT,RRR,0.0D0)/rhocunist(TT,RRR,BB)

END FUNCTION kxxnist


FUNCTION kxxnist_tref(TT,RRR,BB,TREF)

  DOUBLE PRECISION TT, RRR, BB,TREF,kxxnist_tref
	kxxnist_tref = kxxnist_B0(TT,RRR)*rhocunist_tref(TT,RRR,0.0D0,TREF)/rhocunist_tref(TT,RRR,BB,tref)

END FUNCTION kxxnist_tref
			

FUNCTION rhocunist(TT,RRR,BB)
     
	DOUBLE PRECISION TT, RRR, BB, b, rhocunist
	DOUBLE PRECISION rho0, rhoi, rhoiref, rhcu, lgs, poly, corr
	b=abs(BB)
	rho0=1.553D-8/RRR
	rhoi=1.171D-17*(TT**4.49)/(1+4.48D-7*(TT**3.35)*exp(-(50/TT)**6.428))
	rhoiref=0.4531*rho0*rhoi/(rho0+rhoi)      
	rhcu=rho0+rhoi+rhoiref
	IF (b.lt.1D-1) THEN
		 rhocunist=rhcu
	ELSE
		 lgs=0.43429*log(1.553D-8*b/rhcu)
		 poly=-2.662+lgs*(0.3168+lgs*(0.6229+lgs*(-0.1839+lgs*0.01827)))
		 corr=(10**poly)
		 rhocunist=(1.+corr)*rhcu
	ENDIF
	
END FUNCTION rhocunist
			

FUNCTION rhocunist_Tref(TT,RRR,BB,TREF)


	DOUBLE PRECISION TT,RRR,BB,b,rhocunist_Tref,TREF,tref_RRR,t
	DOUBLE PRECISION c0,c0_scale,P1,P2,P3,P4,P5,P6,P7
	DOUBLE PRECISION a0,a1,a2,a3,a4
	DOUBLE PRECISION rho_0,rho_i,rho_i0,rho_n
	DOUBLE PRECISION x,log_x,f_exp,corr

	c0 = 1.553e-8  
	tref_RRR = 273

	P1 = 1.171e-17	
	P2 = 4.49	
	P3 = 3.841e10	
	P4 = 1.14	 
	P5 = 50	
	P6 = 6.428	
	P7 = 0.4531	  

	a0 = -2.662
	a1 = 0.3168
	a2 = 0.6229
	a3 = -0.1839
	a4 = 0.01827


	b=abs(BB)
	T = TT
	
	c0_scale = TREF/tref_RRR
	rho_0 = c0_scale*c0/RRR	  
	rho_i = P1*(T**P2)/(1+P1*P3*(T**(P2-P4))*exp(-((P5/T)**P6))) 
	rho_i0 = P7*rho_i*rho_0/(rho_i+rho_0)
	rho_n = rho_0+rho_i+rho_i0

	IF (b.gt.0.01d0) THEN
		x = c0 * B / rho_n
		log_x = log10(x)	  
		f_exp = a0+a1*log_x+a2*(log_x**2)+a3*(log_x**3) + a4*(log_x**4)
		corr = 10.0d0**f_exp
	ELSE
		corr=0.0d0
	ENDIF
	
	rhocunist_tref=rho_n*(1+corr);
	
END FUNCTION rhocunist_Tref
			

			
!________________________________________________________________                                                         
!                         MATPRO                                
!________________________________________________________________ 			
			

FUNCTION cvcumatpro(TT)

	DOUBLE PRECISION, INTENT(IN) :: TT
	DOUBLE PRECISION             :: cvcumatpro
	REAL TH(27), CAL(27), cp
	DOUBLE PRECISION density
	INTEGER NPOINT

	density = 8960.0

	DATA TH  / 10., 12., 15., 16., 18., 20., 25., 30., 40., 50., 60.,&
						 70., 80., 90.,100.,120.,140.,160.,180.,200.,220.,240.,&
						260.,280.,300.,500.,1000./
	DATA CAL / .86,1.42, 2.7, 3.4, 5. , 7.7, 16., 27., 60., 99.,137.,&
						173.,205.,232.,254.,288.,313.,332.,346.,356.,364.,371.,&
						376.,381.,386.,408.,408./

  IF (TT.LE.10) THEN
		cvcumatpro=(10.8E-6*TT+30.6*(TT/344.5)**3)*1.E3
	ELSE
    NPOINT=27
    CALL LININT (NPOINT,TH,CAL,REAL(TT),cp)
		cvcumatpro = DBLE(cp)
  ENDIF
	cvcumatpro=cvcumatpro*density

END FUNCTION cvcumatpro


FUNCTION kxxmatpro(TT,RRR,BB)

      DOUBLE PRECISION kxxmatpro
      DOUBLE PRECISION TT, RRR, BB
      DOUBLE PRECISION PHI, GAMM, DELTA, EPS, BETA, aa, cc

      PHI = 401.0
      GAMM = 2.5
      DELTA = 41630.0
      EPS = 0.666
      BETA = 3.0E-3
	
      aa=(1. - EXP(-1.*(TT**GAMM)/DELTA)) /PHI
      cc=EPS/TT*(BETA*(BB**.9)+(1./RRR))
      kxxmatpro=1./(aa+cc)

END FUNCTION kxxmatpro


FUNCTION rhocumatpro(TT,RRR,BB)

	DOUBLE PRECISION rhocumatpro
	DOUBLE PRECISION TT, RRR, BB
	DOUBLE PRECISION RHO273, ALFA, tau, aa, cc, rs0, xx, yy
	DATA RHO273,ALFA,tau,aa,cc/1.53E-8,.8,273,.35136808e-2,380.1501/
	rs0 = rho273*((1./RRR) + EXP(alfa*(TT-tau)/TT))
	xx  = BB*rho273 
	yy  = xx/(aa*xx+cc*rs0)
	rhocumatpro = rs0*(yy+1.)

END FUNCTION rhocumatpro



!________________________________________________________________
!                                                              
!                                                              
!                   NB3SN MATERIAL PROPERTIES                                          
!                                                              
!                                                              
!________________________________________________________________



!________________________________________________________________                                                         
!                         CUDI                                
!________________________________________________________________ 	


FUNCTION cvnb3sncudi(QFLAG,TT,BB)

	INTEGER, INTENT(IN) :: QFLAG
	DOUBLE PRECISION, INTENT(IN) :: TT			
	DOUBLE PRECISION, INTENT(IN) :: BB
	DOUBLE PRECISION cvnb3sncudi,dc1,dc2,dc3

!     quench flag = QFLAG
!       -1:  fully SC          (T<Tcs)
!        0:  current sharing   (Tcs<T<TcB0)
!       +1:  fully quenched    (T>TcB0 0)

	IF (QFLAG .lt. 0)THEN               !Superconducting
		 dc1 = 207-3.83*BB+2.86*BB**2
		 dc2 = -110*exp(-0.434*BB)
		 dc3 = 38.8-1.8*BB+0.0634*BB**2
		 cvnb3sncudi = dc3 * TT**3.0D0 + dc2 * TT**2.0D0 + dc1 * TT
	ELSE                                  !Quenched
		IF(TT < 26.113D0) THEN
			cvnb3sncudi = 7.42 * TT**3.0D0 + 1522.0D0 * TT
		ELSE
		IF(TT < 169.416D0) THEN
			cvnb3sncudi = - 61.635 * TT**2.0D0 + 19902 * TT - 305807
		ELSE
			IF(TT < 300.0D0) THEN
				cvnb3sncudi = - 7.4636 * TT**2.0D0 + 4411 * TT + 763801
			ELSE
				cvnb3sncudi = 1415377
			END IF
		END IF
	END IF
	END IF

END FUNCTION cvnb3sncudi

!________________________________________________________________                                                         
!                         NIST                                
!________________________________________________________________ 

FUNCTION cvnb3snnist(TT,BB)

	DOUBLE PRECISION, INTENT(IN) :: TT		
	DOUBLE PRECISION, INTENT(IN) :: BB
	DOUBLE PRECISION cvnb3snnist
	

	DOUBLE PRECISION a0,a1,a2,a3,a4,a5,a6,a7
	DOUBLE PRECISION beta,gamma,rho		
	DOUBLE PRECISION T1,T2,logT,fexp
	DOUBLE PRECISION Tc,Tc0,Bc20
	
	Tc0 = 17.8d0
	Bc20 = 27.012d0		
	T1 = 20d0
	T2 = 400d0

	a0 = 79.78547d0
	a1 = -247.44839d0
	a2 = 305.01434d0
	a3 = -186.90995d0
	a4 = 57.48133d0
	a5 = -6.3977d0
	a6 = -0.6827738d0
	a7 = 0.1662252d0
	beta = 1.241d0/(1000.0d0)
	gamma = 0.138d0
	rho = 8950d0 
	
	IF (BB .lt. Bc20)THEN  
		Tc = Tc0*((1.0d0-BB/Bc20)**0.59d0)        			
	ELSE 
		Tc =  0.0d0
	END IF	
	
	IF (TT .lt. Tc)THEN               !Superconducting
		cvnb3snnist = rho*((beta+3.0d0*gamma/(Tc0**2d0))*(TT**3d0) + gamma*BB/Bc20*TT)
	ELSE  IF  (  (TT .gt. Tc) .AND. (TT .le. T1) )THEN                               
		cvnb3snnist = rho*(beta*(TT**3.0d0) + gamma*TT)
	ELSE  IF  (  (TT .gt. T1) .AND. (TT .le. T2) )THEN				
		logT = LOG10(TT)
		fexp = (a0 + a1*logT + a2*(logT**2d0) + a3*(logT**3d0) + a4*(logT**4d0) + a5*(logT**5d0) + a6*(logT**6d0) + a7*(logT**7d0))
		cvnb3snnist = rho*(10d0**fexp)		
	ELSE  
		cvnb3snnist = rho*(234.89d0+0.0425d0*TT)
	END IF			
			

END FUNCTION cvnb3snnist
			
		

!________________________________________________________________
!                                                               
!                                                               
!                     G10 MATERIAL PROPERTIES                                                   
!                                                              
!                                                               
!________________________________________________________________		

!________________________________________________________________                                                         
!                         NIST                                
!________________________________________________________________ 


FUNCTION cvG10nist(TT)

	DOUBLE PRECISION :: dc_a, dc_b, dc_c, dc_d, dc_e, dc_f, dc_g, dc_h
	DOUBLE PRECISION :: density
	DOUBLE PRECISION, intent(in) :: TT 
	DOUBLE PRECISION :: cvG10nist, p

	dc_a = -2.4083
	dc_b =  7.6006 
	dc_c = -8.2982
	dc_d =  7.3301 
	dc_e = -4.2386
	dc_f =  1.4294
	dc_g = -0.24396
	dc_h =  0.015236 

	density = 1900.0 ! No value in NIST database, so equal to Cryocomp and Fermilab

	
	p = dc_a + dc_b * (LOG10(TT))**1.0D0 + dc_c * (LOG10(TT))**2.0D0 + dc_d * (LOG10(TT))**3.0D0 + dc_e * (LOG10(TT))**4.0D0 + dc_f * (LOG10(TT))**5.0D0 + dc_g * (LOG10(TT))**6.0D0 + dc_h * (LOG10(TT))**7.0D0
	cvG10nist = density * 10**p	  

END FUNCTION cvG10nist		




!________________________________________________________________
!                                                              
!                                                              
!                   Subroutines                                         
!                                                              
!                                                              
!________________________________________________________________



SUBROUTINE LININT (NPOINT,XVETT,YVETT,X,Y)		
				
	INTEGER NPOINT,KK,I
  REAL XVETT(NPOINT),YVETT(NPOINT)			
  REAL X,Y			

  KK=0
	DO i = 1,NPOINT-1		
		IF ((X.GE.XVETT(I)).AND.(X.LT.XVETT(I+1))) KK=I	
	ENDDO				
  IF (X.EQ.XVETT(NPOINT)) KK=NPOINT		
  Y= YVETT(KK)+(YVETT(KK+1)-YVETT(KK))*(X-XVETT(KK))/(XVETT(KK+1)-XVETT(KK))

END SUBROUTINE LININT
			
END MODULE matLIB


