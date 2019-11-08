c___________________________________________________________________________
c     --- USER101: 2D 8-node thermal solid (based on PLANE77)
c     --- uel101.f, element matrices, load vectors, and results
c     --- distribution R1.0, Lucas Brouwer 11/15/2018
c___________________________________________________________________________

      SUBROUTINE uel101 (elem,ielc,elmdat,eomask,nodes,locsvrL,kelreq,
     x kelfil,nr,xyz,u,kelout,zs,zass,damp,gstIF,zsc,zscnr,elvol,elmass,
     x center,elener,edindxL,lcerstL)
			
c   element MODULE SUBROUTINEs
      USE quad2DT9, ONLY : set_dH
      USE quad2DT9, ONLY : set_quad
      USE quad2DT9, ONLY : set_HH
			
c   element MODULE variables (same for all elements - local coords => only set 1x)
      USE quad2DT9, ONLY : qw,nnod,edim,nq
      USE quad2DT9, ONLY : dH,HH,H	

c   inverse			
      USE quad2DT9, ONLY : INVF2
			
c   quench MODULE - includes Jc fitting, current sharing etc.
c      USE quench, ONLY : check_quench
			
c   material property MODULE (includes roxie and others)

c     CU - CUDI,NIST,MATPRO
      USE matLIB, ONLY : cvcucudi,rhocucudi,kxxcudi	
      USE matLIB, ONLY : cvcunist,rhocunist,kxxnist
      USE matLIB, ONLY : kxxnist_tref
      USE matLIB, ONLY : cvcumatpro,kxxmatpro,rhocumatpro
			
c     Nb3Sn - CUDI,MATPRO not addeded yet			
      USE matLIB, ONLY : cvnb3sncudi
      USE matLIB, ONLY : cvnb3snnist
c      USE matLIB, ONLY : cpnb3snmatpro	
		
c     G10 - NIST			
      USE matLIB, ONLY : cvG10nist

c   MFS MODULE SUBROUTINEs					
      USE MFSdata, ONLY : initialize_MFS_TH
			
c   MFS MODULE data (for transfer)
      USE MFSdata, ONLY : MFS 
      USE MFSdata, ONLY : RCvar
			
c___________________________________________________________________________
c     --- Include Decks
c       allows for use of ANSYS definitions and data on distibution
c       see .inc at \ANSYS Inc\v171\ansys\customize\include
c___________________________________________________________________________

c   standard		
#include "impcom.inc"   
#include "echprm.inc"   
#include "elparm.inc"   
#include "ansysdef.inc" 
#include "seldcm.inc"   

c   output etc.
#include "soptcm.inc"
#include "stepcm.inc"
#include "elecom.inc"
#include "elucom.inc"
#include "outpcm.inc"
#include "usvrcm.inc"


      EXTERNAL TrackBegin, TrackEnd
      EXTERNAL erhandler,pstev1,svgidx,svrget,rvrget,plast1,
     x subrd,subwrt,xyzup3,mreuse,creep1,swell1,
     x uep101,usereo,trrot,svrput,hvrput,svpidx,mctac,
     x matxv,vzero,vmove,viinit,vamb,tmpget,fluget,tranx3,
     x vapb1,vapb,eldwrtL,eldwrnL,nlpropcheckL

      EXTERNAL propev
			
c___________________________________________________________________________
c     --- Set Parameters
c___________________________________________________________________________
			
      DOUBLE PRECISION PI,MU0
      PARAMETER (PI=4.0D0*DATAN(1.0D0)) 
      PARAMETER (MU0=4.0D0*PI/1.0D7) 
			
c___________________________________________________________________________
c     --- ANSYS Variables (in/out)
c___________________________________________________________________________

      INTEGER elem         ! element label number (index)
      INTEGER ielc(IELCSZ) ! array of element type characteristics			
      INTEGER elmdat(10)   ! array of element data			
      INTEGER eomask       ! bit pattern for element output (see outpcm)
      INTEGER nodes(nnod)  ! array of element node numbers			
      INTEGER kelreq(10)   ! matrix and load vector form requests
      INTEGER kelfil(10)   ! keys for incoming matrices/vectors			
      INTEGER nr           ! matrix and load vector size (for elem)
      INTEGER kelout(10)   ! keys indicating created matrices and vectors 
			
      LONGINT locsvrL      ! location of saved variables on .esav
      LONGINT edindxL(25)  ! element result data file indices
      LONGINT lcerstL      ! position on result file
			
      DOUBLE PRECISION xyz(6,nnod)    ! nodal coordinates (orig) and rot angle
      DOUBLE PRECISION u(nr,5)        ! element nodal solution values
      DOUBLE PRECISION zs(nr,nr)      ! stiffness/conductivity matrix    kelreq(1)
      DOUBLE PRECISION zass(nr,nr)    ! mass matrix                      kelreq(2)			
      DOUBLE PRECISION damp(nr,nr)    ! damping/specific heat matrix     kelreq(3)	
      DOUBLE PRECISION gstIF(nr,nr)   ! stress stiffness matrix          kelreq(4)	
      DOUBLE PRECISION zsc(nr)        ! applied force vector           kelreq(5)	
      DOUBLE PRECISION zscnr(nr)      ! n-r or imaginary force vector  kelreq(6)
c                                       or imaginary f vector          kelreq(7)
      DOUBLE PRECISION elvol          ! element volume	
      DOUBLE PRECISION elmass         ! element mass
      DOUBLE PRECISION center(3)      ! centroid location
      DOUBLE PRECISION elener(5)      ! element energies

 
c___________________________________________________________________________
c     --- Local Variables 
c___________________________________________________________________________
		
			
c     to calc local shape functions and deriv on first CALL only				
      LOGICAL ftime
      DATA ftime/.true./
      SAVE ftime
			
c     loop ints			
      INTEGER i,j,k,l
						
c     real constants
      DOUBLE PRECISION rvr(6)
      DOUBLE PRECISION Sc,Nc,fcond,fsc
      DOUBLE PRECISION scl,f1,f2,f3			

c     element matrices 
      DOUBLE PRECISION NtN(nnod,nnod)		
      DOUBLE PRECISION ijac(edim,edim),jac(edim,edim)	
      DOUBLE PRECISION djac
      DOUBLE PRECISION xyz2D(edim,nnod)				
      DOUBLE PRECISION HL(nnod),dHL(nnod,edim),HL2(nnod,1)
      DOUBLE PRECISION dN(nnod,edim)

c     non-linear solution (NR-OPT)
      DOUBLE PRECISION ui(nr)       

c     quench properties
      INTEGER qflag
			
c     material properties
      DOUBLE PRECISION prop(3)
      DOUBLE PRECISION dens,kxx,cv
      DOUBLE PRECISION trefRRR
      DOUBLE PRECISION kxxCu,kxxNb3Sn,kxxG10
      DOUBLE PRECISION cvCu,cvNb3Sn,cvG10
      INTEGER lp(3)	
			
c     variables for loading			
      DOUBLE PRECISION hgnev,hgn(nnod)
      DOUBLE PRECISION hgnbeg(nnod),hgnEND(nnod)
			
			
c     variables for ansys data retrive/write (mat prop)
      DOUBLE PRECISION tempev,Bev,RRR
      DOUBLE PRECISION tem(nnod)

c     saved variables 			
      DOUBLE PRECISION svindx(20)

c     saved variable vector
      DOUBLE PRECISION msvr(7)
			
c     MFS transfer loc		
      INTEGER ERpos
			
c     element output
      DOUBLE PRECISION Nmisc(24)
      DOUBLE PRECISION EGeo(11)
			
c     elmdat pointers		
      INTEGER mat,ireal,nrvr,nssvr,nnrsvr,nusvr
      INTEGER khgn
	
c     element mapping 
      INTEGER Emap
			
c			variables for tracking, logging, etc.
      INTEGER er
      DATA er/6000/	
      SAVE er
			

c________________________________________________________________________________		 
c set IELC Pointers, saved variable index, and set in/out element defaults
c________________________________________________________________________________

c     elmdat pointers defined in elparm
      mat = elmdat(EL_MAT)
      ireal = elmdat(EL_REAL)
			
c     ielc pointers defined in echprm and elccmt
      nrvr = ielc(NMTRLC)
      nssvr = ielc(NMSSVR)
      nnrsvr = ielc(NMNSVR)
      nusvr = ielc(NMUSVR)
			
c     get the svr index vector
      CALL svgidx(locsvrL,svindx(1))
			
c     element defaults
      elvol = 0.0d0
      elmass = 0.0d0
      center = 0.0d0
      elener = 0.0d0
			
c________________________________________________________________________________		 
c Get Element Real Constants -> Mapping of Magnetic to Thermal Elem
c________________________________________________________________________________
			
c --- get the element real constant data
      CALL rvrget (elem,ireal,ielc(1),nrvr,rvr(1))
      EMAP = NINT(rvr(1))
			
c________________________________________________________________________________		 
c Get Saved Variables
c________________________________________________________________________________
			
c     get saved variables (ERpos)
      CALL svrget (svindx(1),4,nnrsvr,msvr(1))
						
c________________________________________________________________________________		 
c Initialize RC either from real const or from EMAP mapping, ini. MFS too
c  -- find location in shared data of mapped EMAG elem
c  -- then set the needed parameters from there
c________________________________________________________________________________	

      IF (ielc(KYOP2) .eq. 0) THEN
			
        Bev = rvr(2)
        qflag = NINT(rvr(3))
        RRR = rvr(4)
        fcond = rvr(5)
        fsc = rvr(6)	
				
      ELSEIF (ielc(KYOP2) .eq. 1) THEN	
			
        IF (kfstps.eq.1) THEN		
          CALL initialize_MFS_TH(EMAP,ERpos)	
          msvr(1) = ERpos
        ELSE 				
          ERpos = NINT(msvr(1))			
        ENDIF			
	
        Bev = MFS(1,ERpos)
        qflag = MFS(2,ERpos)					
        RRR = RCvar(1,ERpos)
        fcond = RCvar(2,ERpos) 
        fsc = RCvar(3,ERpos) 

      ENDIF				

						
c________________________________________________________________________________		 
c Initialize FEM 
c________________________________________________________________________________	

c set up the gaussian points and then find the shape function gradient matrix in local coords (r,s,t): dN at gaussian point locations
c this is the same for all elements so ONLY NEED TO CALL ONCE PER ELEMENT TYPE (i.e. indepENDent of node location xyz)
      IF(ftime) THEN
        CALL set_quad()				
        CALL set_dH()	
        CALL set_HH()				
        ftime=.false. 
      ELSE				
      ENDIF		
					
c________________________________________________________________________________		 
c Get Temperature 
c________________________________________________________________________________

c  get nodal temps from DOF vector instead of using tmpget
      tem = u(:,1)
		 
c  set elem temp as average of nodal temps	
      scl = 0.125d0 
      tempev=scl*(tem(1)+tem(2)+tem(3)+tem(4)) 
     x     +scl*(tem(5)+tem(6)+tem(7)+tem(8))
		 
			
c________________________________________________________________________________		 
c Get Heat Generation Loads
c________________________________________________________________________________

      CALL hgnget (elem,ielc(1),nnod,nodes(1),nnod,hgnbeg(1),
     x     hgn(1),hgnEND(1),khgn)

      scl = 0.125d0 
      hgnev=scl*(hgn(1)+hgn(2)+hgn(3)+hgn(4)) 
     x     + scl*(hgn(5)+hgn(6)+hgn(7)+hgn(8))

			
c________________________________________________________________________________		 
c Get Material Properties From matLIB.f
c________________________________________________________________________________
      IF (ielc(KYOP1) .eq. 0) THEN
			
c________________________________________________________________________________		 
c Get Material Properties for stabilizer
c________________________________________________________________________________			
			
      IF (ielc(KYOP4) .eq. 0) THEN			
c     !!!!!!!!!!!!!!!!!!!!!!!!!!!
c         for copper 
c     !!!!!!!!!!!!!!!!!!!!!!!!!!!
c     PICK cv (key6)
        IF (ielc(KYOP6) .eq. 0) THEN
          cvCu = cvcunist(tempev)
        ELSEIF (ielc(KYOP6) .eq. 1) THEN				
          cvCu = cvcucudi(tempev)				
        ELSEIF (ielc(KYOP6) .eq. 2) THEN				
          cvCu = cvcumatpro(tempev)				
        ELSEIF (ielc(KYOP6) .eq. 3) THEN				
c         add userCucp.f
        ELSE
c         add warning about not selecting cvCu
        ENDIF
				
c     PICK kxx (key7)
        IF (ielc(KYOP7) .eq. 0) THEN
          kxxCu = kxxnist(tempev,RRR,Bev)
c         now include reference temp for RRR 
          trefRRR = 295.0d0					
          kxxCu = kxxnist_tref(tempev,RRR,Bev,trefRRR)
        ELSEIF (ielc(KYOP7) .eq. 1) THEN				
          kxxCu = kxxcudi(tempev,RRR,Bev)				
        ELSEIF (ielc(KYOP7) .eq. 2) THEN				
          kxxCu = kxxmatpro(tempev,RRR,Bev)			
        ELSEIF (ielc(KYOP7) .eq. 3) THEN				
c         add userCukxx.f
        ELSE
c         add warning about not selecting cvCu
        ENDIF
				
      ELSEIF (ielc(KYOP4) .eq. 1) THEN
c     !!!!!!!!!!!!!!!!!!!!!!!!!!!
c         for Ag 
c     !!!!!!!!!!!!!!!!!!!!!!!!!!!	
		
c     PICK cv (key12)
c     PICK kxx (key13)		
      ELSE
c       add error for that keyopt(4) is not set correctly
      ENDIF	

c________________________________________________________________________________		 
c Get Material Properties for superconductor
c________________________________________________________________________________

      IF (ielc(KYOP3) .eq. 0) THEN			
c     !!!!!!!!!!!!!!!!!!!!!!!!!!!
c         for NbTi 
c     !!!!!!!!!!!!!!!!!!!!!!!!!!!

      ELSEIF (ielc(KYOP3) .eq. 1) THEN
c     !!!!!!!!!!!!!!!!!!!!!!!!!!!
c         for Nb3Sn 
c     !!!!!!!!!!!!!!!!!!!!!!!!!!!				
			
c     PICK cv (key9)
        IF (ielc(KYOP9) .eq. 0) THEN
          cvNb3Sn = cvnb3snnist(tempev,Bev)	
        ELSEIF (ielc(KYOP9) .eq. 1) THEN				
          cvNb3Sn = cvnb3sncudi(qflag,tempev,Bev)					
        ELSEIF (ielc(KYOP9) .eq. 2) THEN				
c         add userNb3sncp.f			
        ELSE
c         add warning about not selecting cvG10 keyopt(11)
        ENDIF
				
				
      ELSEIF (ielc(KYOP3) .eq. 2) THEN
c     !!!!!!!!!!!!!!!!!!!!!!!!!!!
c         for Bi2212 
c     !!!!!!!!!!!!!!!!!!!!!!!!!!!			
      ELSE
c       add error for that keyopt(3) is not set correctly

      ENDIF	

c________________________________________________________________________________		 
c Get Material Properties for non-conductor
c________________________________________________________________________________

      IF (ielc(KYOP5) .eq. 0) THEN			
c     !!!!!!!!!!!!!!!!!!!!!!!!!!!
c         for G10 
c     !!!!!!!!!!!!!!!!!!!!!!!!!!!

c     PICK cv (key11)
        IF (ielc(KYOP11) .eq. 0) THEN
          cvG10 = cvG10nist(tempev)
        ELSEIF (ielc(KYOP11) .eq. 1) THEN				
c         add userG10cp.f				
        ELSE
c         add warning about not selecting cvG10 keyopt(11)
        ENDIF
				
      ELSE
c       add error for that keyopt(5) is not set correctly

      ENDIF	

c________________________________________________________________________________		 
c Homogenize Material Properties
c________________________________________________________________________________

      f1 = 1d0-fcond
      f2 = fcond*fsc
      f3 = fcond*(1d0-fsc)
			
c     note that the thermal conductivity of Nb3Sn and epoxy are so much smaller than Cu that it is neglected
      kxx = kxxCu*f3		
			
c     now just do mix using cv 
      cv =cvG10*f1+cvNb3Sn*f2+cvCu*f3		
			
      ELSEIF (ielc(KYOP1) .eq. 1) THEN
c________________________________________________________________________________		 
c Get Material Properties input Homoginized (for material it is meshed with)
c________________________________________________________________________________
			
c     get material prop from ANSYS (dens+kxx,c)
        lp(1) = 13  !dens=13
        lp(2) = 16  !kxx=16
        lp(3) = 22  !c=22
c     default also gets 17 and 26  (KYY,ENTH)			
        CALL propev(elem,mat,lp(1),tempev,prop(1),3)
						
        dens = prop(1)			
        kxx = prop(2)
        cv = prop(3)*dens					
      ELSE
c       add error for that keyopt(1) is not set correctly
      ENDIF					
				
			
c shape the full 3D node location into 2D 			
      xyz2D = xyz(1:edim,:)  

c initialize the matrices to be formed (before sum), and sEND back 0's for other matrices IF they are requested	
      IF(kelreq(1).eq.1)THEN
         zs = 0.0d0 
         kelout(1) = 1
      ENDIF
      IF(kelreq(2).eq.1)THEN
         zass = 0.0d0
c         kelout(2) = 1     removed because mreuse takes it out in plane53
         kelout(2) = 1     
      ENDIF
      IF(kelreq(3).eq.1)THEN
         damp = 0.0d0
         kelout(3) = 1
      ENDIF
      IF(kelreq(4).eq.1)THEN
         gstIF = 0.0d0
         kelout(4) = 1
      ENDIF				 		
			IF(kelreq(5).eq.1)THEN
         zsc = 0.0d0
         kelout(5) = 1
      ENDIF
      IF(kelreq(6).eq.1)THEN
         zscnr = 0.0d0
         kelout(6) = 1
      ENDIF
      IF(kelreq(7).eq.1)THEN
         zscnr = 0.0d0
         kelout(7) = 1
      ENDIF
			

c ___________________________________________________________________________________________________________	
c integration loop: do all at once to avoid duplication operations since some are shared (jac,curl-curl,etc.) 		
c ___________________________________________________________________________________________________________	

      DO 801 i = 1,nq
        HL = H(:,i)		
        HL2 = reshape(H(:,i),(/nnod,1/))
        dHL = reshape(dH(:,:,i),(/nnod,2/))
        jac = matmul(xyz2D,dHL)
        CALL INVF2(2,jac,ijac,djac)	
        dN = matmul(dhL,ijac)
        NtN = matmul(HL2,transpose(HL2))	

c       find element area
        elvol = elvol + qw(i)*djac	
					
c __________________________________________________________________________________________________________			
c IF requested, make stIFfness matrix by integrating with global coords using the nine quadrature points		
c __________________________________________________________________________________________________________		
        IF(kelreq(1).eq.1) THEN   				
          zs = zs + kxx*qw(i)*djac*matmul(dN,transpose(dN))							
        ENDIF
			
c __________________________________________________________________________________________________________			
c IF requested, make load vector by integrating with global coords using the nine quadrature points		
c __________________________________________________________________________________________________________		
        IF(kelreq(5).eq.1) THEN  
          zsc = zsc + hgnev*qw(i)*djac*HL					
        ENDIF
	
c __________________________________________________________________________________________________________			
c IF requested, make damping matrix by integrating with global coords using the nine quadrature points		
c __________________________________________________________________________________________________________		
        IF(kelreq(3).eq.1) THEN   
c          damp = damp + dens*cv*qw(i)*djac*NtN
          damp = damp + cv*qw(i)*djac*NtN
        ENDIF			
	
  801   	CONTINUE			


c __________________________________________________________________________________________________________			
c Find Restoring Force (NROPT: NL-matrices)		
c __________________________________________________________________________________________________________

c     restoring force 
      IF(kelreq(6).eq.1) THEN 
          ui = u(:,1)
          zscnr = matmul(zs,ui)
      ENDIF
			
c __________________________________________________________________________________________________________			
c output to results file	
c __________________________________________________________________________________________________________		
      IF (outkey.eq.1) THEN
				
c       SEE EOMASK Definitions in outpcm.inc
        IF (btest(eomask,PDBE)) THEN
				
c         Misc Variables
          IF (btest(eomask,13)) THEN	

						Nmisc = 0d0
			      Nmisc(1) = tempev
			      Nmisc(2) = Emap
			      Nmisc(3) = Bev
						Nmisc(4) = DBLE(qflag)				
						Nmisc(5) = RRR	
						Nmisc(6) = fcond
						Nmisc(7) = fsc								
						Nmisc(8) = hgnev												
						Nmisc(9) = cvCu
						Nmisc(10) = kxxCu	
						Nmisc(11) = cvNb3Sn
						Nmisc(12) = kxxNb3Sn	
						Nmisc(13) = cvG10
						Nmisc(14) = kxxG10							
						Nmisc(15) =	cv				
						Nmisc(16) =	kxx				
						Nmisc(17) =	dens						
				
						CALL eldwrtL (elem,13,lcerstL,edindxL(1),18,Nmisc(1))   !write NMISC to res file
					
          ENDIF		

c         FIELD FLUXES	
c         output imported B field as "thermal flux" for debugging/viewing
          IF (btest(eomask,PDBFLX)) THEN					
c						CALL eldwrtL (elem,EDEFX,lcerstL,edindxL(1),12,str(1))   !write B-flux to results file
          ENDIF		
					
c         VOLU/AREA and Energy			
          IF (btest(eomask,4)) THEN						
c						Egeo(1:5) = elener
						Egeo(1) = elvol
						CALL eldwrtL (elem,4,lcerstL,edindxL(1),11,Egeo(1))   
          ENDIF	
					
        ENDIF					
      ENDIF	
	
	
c -------------------------------------------------------------
c     Output saved variables and update pointers
c -------------------------------------------------------------			

c      write out user saved element variables
      CALL svrput (svindx(1),4,nnrsvr,msvr(1))
			
c     write out the svr index vector, which updates locsvrL pointer to next element
      CALL svpidx (locsvrL,svindx(1))
					
   
      END SUBROUTINE uel101
			


	
	