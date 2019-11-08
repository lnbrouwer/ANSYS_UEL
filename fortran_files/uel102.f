c___________________________________________________________________________
c     --- USER102: 2D 8-node magnetic solid (based on PLANE53)
c     --- uel102.f, element matrices, load vectors, and results
c     --- distribution R1.0, Lucas Brouwer 11/15/2018
c___________________________________________________________________________


      SUBROUTINE uel102 (elem,ielc,elmdat,eomask,nodes,locsvrL,kelreq,
     x kelfil,nr,xyz,u,kelout,zs,zass,damp,gstIF,zsc,zscnr,elvol,elmass,
     x center,elener,edindxL,lcerstL)

c   module subroutines
      USE quad2D, ONLY : set_dH
      USE quad2D, ONLY : set_quad
      USE quad2D, ONLY : set_HH
			
c   module variables (same for all elements - local coords => only set 1x)
      USE quad2D, ONLY : qw,nnode,edim,nq
      USE quad2D, ONLY : dH,HH,H	

c   quench module - includes Jc fitting, current sharing etc.
      USE quench, ONLY : check_quench

c   inverse			
      USE quad2D, ONLY : INVF2
			
c   material property module (includes roxie and others)

c     CU - CUDI,NIST,MATPRO
      USE matLIB, ONLY : rhocucudi
      USE matLIB, ONLY : rhocunist
      USE matLIB, ONLY : rhocunist_tref
      USE matLIB, ONLY : rhocumatpro
			
c   MFS module subroutines					
      USE MFSdata, ONLY : initialize_MFS_EM
      USE MFSdata, ONLY :	initialize_var		
      USE MFSdata, ONLY : count_ne
			
c   MFS module data (for transfer)
      USE MFSdata, ONLY : MFS 

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
     x uep102,usereo,trrot,svrput,hvrput,svpidx,mctac,
     x matxv,vzero,vmove,viinit,vamb,tmpget,fluget,tranx3,
     x vapb1,vapb,eldwrtL,eldwrnL,nlpropcheckL

      EXTERNAL erinqr,wrinqr,elLenPsvrBuf
      integer  erinqr,wrinqr,elLenPsvrBuf
      EXTERNAL         dist3d
      DOUBLE PRECISION dist3d
			
c     EXTERNAL - ADDED BY LUCAS
      EXTERNAL edcget
      EXTERNAL edciqr
      integer edciqr
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
      INTEGER nodes(nnode) ! array of element node numbers			
      INTEGER kelreq(10)   ! matrix and load vector form requests
      INTEGER kelfil(10)   ! keys for incoming matrices/vectors			
      INTEGER nr           ! matrix and load vector size (for elem)
      INTEGER kelout(10)   ! keys indicating created matrices and vectors 
			
      LONGINT locsvrL      ! location of saved variables on .esav
      LONGINT edindxL(25)  ! element result data file indices
      LONGINT lcerstL      ! position on result file
			
      DOUBLE PRECISION xyz(6,nnode)  ! nodal coordinates (orig) and rot angle
      DOUBLE PRECISION u(nr,5)       ! element nodal solution values
      DOUBLE PRECISION zs(nr,nr)     ! stiffness/conductivity matrix    kelreq(1)
      DOUBLE PRECISION zass(nr,nr)   ! mass matrix                      kelreq(2)			
      DOUBLE PRECISION damp(nr,nr)   ! damping/specific heat matrix     kelreq(3)	
      DOUBLE PRECISION gstIF(nr,nr)  ! stress stiffness matrix          kelreq(4)	
      DOUBLE PRECISION zsc(nr)       ! applied force vector           kelreq(5)	
      DOUBLE PRECISION zscnr(nr)     ! n-r or imaginary force vector  kelreq(6)
c                                      or imaginary f vector          kelreq(7)
      DOUBLE PRECISION elvol         ! element volume	
      DOUBLE PRECISION elflux        ! element linked flux	
      DOUBLE PRECISION elmass        ! element mass
      DOUBLE PRECISION center(3)     ! centroid location
      DOUBLE PRECISION elener(5)     ! element energies

 
c___________________________________________________________________________
c     --- Local Variables 
c___________________________________________________________________________

c     submatrices for stIFfness and damping (8x8 -> are THEN padded into full matrices)			
      DOUBLE PRECISION Kaa(8,8),Kai(8,8),Kei(8,8),Kee(8,8)  
      DOUBLE PRECISION Caa(8,8),Cea(8,8)

c     shape function matrix form operations (to avoid repeated generation when used for multiple matrices)			
      DOUBLE PRECISION curlcurl(8,8),NtN(8,8),Npat(8,8)
      DOUBLE PRECISION djac
			
c     non-linear solution (NR-OPT)
      DOUBLE PRECISION ui(nr)     
			
c     extracted Az from full solu vector
      DOUBLE PRECISION uB(8,1)       ! element nodal solution values
      DOUBLE PRECISION udBdt(8,1)    ! time deriv of element nodal solution values
	
c     first time -> to set local elem shape functions and their deriv	
      Logical ftime
      DATA ftime/.true./
      SAVE ftime

c     loop counters, permute, etc.	 		
      INTEGER i,j,k,l
      INTEGER perm(4)
 		
c     real constant vector
      DOUBLE PRECISION rvr(18)
      DOUBLE PRECISION Sc,Nc,fcond,fsc,dirz
      DOUBLE PRECISION Lc,Li,RRR,Lpp,feff
      DOUBLE PRECISION crat
			
c     area fractions, scaling factor		
      DOUBLE PRECISION f1,f2,f3	
      DOUBLE PRECISION scl
			
c     saved variable vector
      DOUBLE PRECISION msvr(3)
			
c     MFS field transfer position
      INTEGER ERpos
			
c     quench properties
      INTEGER qflag
      DOUBLE PRECISION Bev
      DOUBLE PRECISION I0,Je,Jcw,Jc
      DOUBLE PRECISION IfCu,scIfCu
      DOUBLE PRECISION Jc0,Tc0,Bc0,ap,p
      DOUBLE PRECISION Jstemp(3)		
			
c     IFCC properties
      DOUBLE PRECISION tau,tauMult,Ctau
      DOUBLE PRECISION Mx,My

c     element matrices (set each time)		
      DOUBLE PRECISION xyz2D(edim,nnode)		
      DOUBLE PRECISION dHL(nnode,edim),HL(nnode),dN(nnode,edim)
      DOUBLE PRECISION curldN(edim,nnode),HL2(nnode,1)
      DOUBLE PRECISION ijac(edim,edim),jac(edim,edim)	

c     element output			
      DOUBLE PRECISION Bm(edim,nq),Bout(3*nq),Bx(nq,1),By(nq,1)
      DOUBLE PRECISION Hout(3*nq)
      DOUBLE PRECISION dBdt(edim,nq)
      DOUBLE PRECISION JHtau,JHres
      DOUBLE PRECISION JHcon
      DOUBLE PRECISION JHeat(1)
      DOUBLE PRECISION EGeo(11)
      DOUBLE PRECISION Nmisc(30),Smisc(5)		
      DOUBLE PRECISION NmiscCH(18)
      DOUBLE PRECISION Bmod			
      DOUBLE PRECISION Sene
			
c     material properties	
      DOUBLE PRECISION mu,rsvx,trefRRR
      DOUBLE PRECISION prop(2)
      DOUBLE PRECISION tempev
      DOUBLE PRECISION tref,tembeg(nnode)
      DOUBLE PRECISION tem(nnode),temEND(nnode)
      INTEGER lp(2)			
      INTEGER mat,ireal,nrvr,nmmvr,nusvr,nfsvr
      INTEGER plopt,cropt,swopt,ktherm
      LONGINT keyplsL,keycrpL,keyswlL

c     je get
      DOUBLE PRECISION jbeg(4*nnode), jem(4*nnode),jend(4*nnode)		
      INTEGER iex
			
c     saved variables 			
      INTEGER svindx(20) 
			
c     variables for tracking, logging, etc.
      INTEGER er
      DATA er/600/	
      SAVE er	 
			
			
			
c________________________________________________________________________________		 
c set pointers, get index, set element defaults 
c________________________________________________________________________________				

c     elmdat pointers defined in elparm
      mat = elmdat(EL_MAT)
      ireal = elmdat(EL_REAL)
			
c     ielc pointers defined in echprm and elccmt
      nrvr = ielc(NMTRLC)
      nusvr = ielc(NMUSVR)
      nfsvr = ielc(NMFSVR)
			
c     get the svr index vector
      CALL svgidx(locsvrL,svindx(1))
			
c     default ANSYS in/out var to zero
      elvol = 0.0d0
      elflux = 0.0d0
      elmass = 0.0d0
      center = 0.0d0
      elener = 0.0d0
			
			
c________________________________________________________________________________		 
c Initialize FEM / and MFS settings if needed
c________________________________________________________________________________	

c     set up the gaussian points and then find the shape function gradient matrix in local coords (r,s,t): dN at gaussian point locations
c     this is the same for all elements so ONLY NEED TO CALL ONCE PER TYPE (i.e. indepedent of node location xyz)
c     if MFS is used with field transfer initialize the element counter etc. in the MFSdata module
      IF(ftime) THEN
        CALL set_quad()				
        CALL set_dH()	
        CALL set_HH()	
c       initialize MFS module variables				
        IF (ielc(KYOP2) .eq. 1) THEN		
          CALL initialize_var()
        ENDIF						
        ftime=.false.   
      ELSE    
      ENDIF		
			
			
c________________________________________________________________________________		 
c Get Element Real Constants
c________________________________________________________________________________
		
c     get the element real constant data
      CALL rvrget (elem,ireal,ielc(1),nrvr,rvr(1))
			Sc = rvr(1)
			Nc = rvr(2)
      fcond = rvr(3)
      fsc = rvr(4)
      dirz = rvr(5)
      Lc = rvr(6)
      Li = rvr(7)
      RRR = rvr(8)
      Lpp = rvr(9)
      feff = rvr(10)
      Ctau = rvr(11)			
      taumult = rvr(12)

c     set ratio of coil lengths
      crat = Li/Lc
			
c     Jc fit parameters			
      Jc0 = rvr(13)
      Tc0 = rvr(14)
      Bc0 = rvr(15)
      ap = rvr(16)
      p = rvr(17)			
      scIFCU = rvr(18)				

c________________________________________________________________________________		 
c Get Saved Variables
c________________________________________________________________________________
			
c     get saved variables (Bev)
      CALL svrget (svindx(1),2,nfsvr,msvr(1))
      Bev = msvr(1)
			
c________________________________________________________________________________		 
c Add element to counter for MFS module if needed
c  -- then populate mapping vector the 2nd time called (after counting)
c  -- place the location in the shared data as the 2nd saved var (ERpos)
c  -- set single transfer values RRR,fcond,fsc
c________________________________________________________________________________	

      IF (ielc(KYOP2) .eq. 1) THEN
        IF (kfstps.ne.1) THEN
          IF (msvr(3) .eq. -1) THEN	
            CALL initialize_MFS_EM(elem,ERpos,RRR,fcond,fsc)
						msvr(2) = ERpos
            msvr(3) = 1
          ELSE
            ERpos = msvr(2)
          ENDIF					
        ELSE			
          CALL count_ne()
          er = er + 1			
          msvr(3) = -1		
        ENDIF				
      ENDIF
			
c________________________________________________________________________________		 
c Compute current Bev for updated mat prop
c________________________________________________________________________________
c     compute (BeV) so the material properties can update within the time step 
c     instead of using the field from the old timestep			
      IF (kfstps.ne.1) THEN
			
        xyz2D = xyz(1:edim,:)  
        Bmod = 0.0d0			
        DO 22 i = 1,nq
          dHL = RESHAPE(dH(:,:,i),(/8,2/))
          jac = matmul(xyz2D,dhL)
          CALL INVF2(2,jac,ijac,djac)		
          dN = matmul(dhL,ijac)																
          DO 20 k = 1,nnode				
						curldN(1,k) = dN(k,2)
						curldN(2,k) = -dN(k,1)
  20   	CONTINUE			
          IF (ielc(KYOP1) .eq. 0) THEN			
            Bm(:,i) = matmul(curldN,u(:,1))	        
          ELSEIF (ielc(KYOP1) .eq. 2) THEN	
            DO 21 k = 1,8				
						  uB(k,1) = u(1+(k-1)*3,1)
  21   	  CONTINUE	
            Bm(:,i) = matmul(curldN,uB(:,1))
          ENDIF
					
c         Bmod = average of |B| at integration points f
          Bmod=Bmod+sqrt(Bm(1,i)*Bm(1,i)+Bm(2,i)*Bm(2,i))/nq				
					
  22   CONTINUE					 	
      ENDIF	

c     set element field for mat prop etc. as ave at int. points		
      Bev = Bmod

			
c________________________________________________________________________________		 
c Set a first temp, field, and current to avoid evalu prop with non-ini var
c________________________________________________________________________________
c     initial conditions and tref to this in normal ANSYS
c     in the future it would be good to adapt to this for an "IC" spec.
      IF (kfstps.eq.1) THEN
c       temp is handled with tref
        Bev = 5.0d0
        I0 = 20.0d0
      ENDIF		
			
			
c________________________________________________________________________________		 
c Get Temperature 
c________________________________________________________________________________
c     get tref
      CALL nlpropcheckL (mat,tref,keyplsL,plopt,keycrpL,cropt,
     x                  keyswlL,swopt)	
		 
c     get nodal temps
      CALL tmpget (elem,ielc(1),nnode,nodes(1),tref,nnode,tembeg(1),
     x     tem(1),temEND(1),ktherm)
	
c     set elem temp as average of nodal temps	
      scl = 0.125d0
      tempev=scl*(tem(1)+tem(2)+tem(3)+tem(4)) 
     x     +scl*(tem(5)+tem(6)+tem(7)+tem(8))			
		

c________________________________________________________________________________		 
c Get Material Properties 
c    - choice is set by Keyopt(5)
c        KYOP5 = 0   ->  ANSYS table (temp only) 
c        KYOP5 = 1   ->  NIST (temp,RRR,B)     
c        KYOP5 = 2   ->  CUDI (temp,RRR,B)  
c        KYOP5 = 3   ->  Matpro (temp,RRR,B) 
c        KYOP5 = 4   ->  userCu.f (temp,RRR,B)  [user fortran function] 
c         - need to include warnings when temp/field go out of mat prop range
c         - use err-handler for this and linear extrap prop. from there
c________________________________________________________________________________

c     default mu to mu0
        mu = MU0
				
      IF (ielc(KYOP5) .eq. 0) THEN
c       ANSYS input table
        lp(1) = 31  !murx=31
        lp(2) = 19  !rsvx = 19
        CALL propev(elem,mat,lp(1),tempev,prop(1),2)
        mu = prop(1)*MU0
        rsvx = prop(2)
				
      ELSEIF (ielc(KYOP5) .eq. 1) THEN
c       NIST
c       now define RRR from reference temp of 295 
        trefRRR = 295.0d0
        rsvx = rhocunist_tref(tempev,RRR,Bev,trefRRR)

      ELSEIF (ielc(KYOP5) .eq. 2) THEN
c       CUDI
        rsvx = rhocucudi(tempev,RRR,Bev)

      ELSEIF (ielc(KYOP5) .eq. 3) THEN
c       Matpro			
        rsvx = rhocumatpro(tempev,RRR,Bev)
				
      ELSEIF (ielc(KYOP5) .eq. 4) THEN
c       user defined fortran function -> can be compiled with dummy 

      ELSE
c       add error for wrong setting of keyopt(8)
      ENDIF	

	
			
c________________________________________________________________________________		 
c Set IFCC tau based on rho_cu(T,RRR,B), Lp, and feff
c________________________________________________________________________________
c    - choice is set by Keyopt(8)
c        KYOP8 = 0   ->  IFCC, tau calculated from mat
c        KYOP8 = 1   ->  No IFCC 
c        KYOP8 = 2   ->  IFCC with fixed tau input as real const 11

      IF (ielc(KYOP8) .eq. 0) THEN
c       use Wilson, which differs from Lorenzo by a factor of two		
        tau = MU0*Lpp*Lpp/(8d0*pi*pi*rsvx*Feff)			
      ELSEIF (ielc(KYOP8) .eq. 1) THEN
        tau = 0.0d0
      ELSEIF (ielc(KYOP8) .eq. 2) THEN
        tau = Ctau
      ELSE
c       add error for wrong setting of keyopt(8) 
      ENDIF
			
			
c_____________________________________________________________________________		 
c fractions of different regions
c  f1 = cond
c  f2 = SC
c  f3 = stabilizer
c________________________________________________________________________________
      f1 = 1d0-fcond
      f2 = fcond*fsc
      f3 = fcond*(1d0-fsc)

c________________________________________________________________________________		 
c find current densities	(Je,Jcw, etc.)
c________________________________________________________________________________

      IF (ielc(KYOP1) .eq. 0) THEN
c     Je loading k(1)=0		

c     EDCGET is replaced with CRDGET to be compatabile with ANSYS 19.1 
c			  CALL edciqr(elem,1)   !returns number of current densities stored for current element
c        IF(edciqr(elem,1).gt.0.25d0) THEN
c          CALL edcget(elem,Jstemp(1))			
c          Je = Jstemp(3)	

       CALL crdget(elem,ielc(1),8,nodes(1),8,jbeg(1),jem(1),jend(1),iex)
			 
       IF(iex.gt.0.25d0) THEN			 
       scl = 0.125d0
       Je=scl*(jem(3)+jem(7)+jem(11)+jem(15)) 
     x     +scl*(jem(19)+jem(23)+jem(27)+jem(31))

       ELSE
          Je = 0.0d0
       ENDIF
        Jcw = abs(Je/f2)	
					
      ELSEIF (ielc(KYOP1) .eq. 2) THEN
c     Circuit loading k(1)=2	
        IF (kfstps.eq.0) THEN
          I0 = u(3,1)		
        ELSE
c        here it would be better to add an user input IC
         I0 = 20.0d0
        ENDIF
        Jcw = abs(I0*Nc/(Sc*f2))					
      ELSE		
c     not conductor - or wrong selection of keyoption
        Jcw = 0.0d0
      ENDIF

c________________________________________________________________________________		 
c Find Quench State and Set Resistive Loss
c    - choice is set by Keyopt(7)
c        KYOP7 = 0   ->  check for quench with current sharing
c        KYOP7 = 1   ->  no current sharing or quenching    
c        KYOP7 = 2   ->  force full quench 
c        KYOP7 = 3   ->  check for quench with no current sharing -> abrupt transition 
c________________________________________________________________________________

      IF (ielc(KYOP7) .eq. 0) THEN
c     check for quench and use current sharing to find % in stabilizer			
c        quench flag: qflag
c         -1:  fully SC          (T<Tcs)
c          0:  current sharing   (Tcs<T<TcB0)
c         +1:  fully quenched    (T>TcB0)		
		
c       find quench state based on Jc and choice of short-sample fit
        CALL check_quench(Jcw,tempev,Bev,Jc,qflag,IfCu,Jc0,Tc0,Bc0,ap,p)
				
      ELSEIF (ielc(KYOP7) .eq. 3) THEN	
        CALL check_quench(Jcw,tempev,Bev,Jc,qflag,IfCu,Jc0,Tc0,Bc0,ap,p)	
c       force abrupt transition (no current sharing)
        IF (IfCu .gt. 0.0d0) THEN  
          IfCu = 1.0d0		
          qflag = 1					
        ENDIF 			
				
      ELSEIF (ielc(KYOP7) .eq. 1) THEN	
c     force SC		
        IfCu = 0.0d0		
        qflag = -1	
				
      ELSEIF (ielc(KYOP7) .eq. 2) THEN	
c     force quenched	
        IfCu = 1.0d0		
        qflag = 1	
				
      ELSE
c       add error for wrong setting of keyopt(7) 
			
      ENDIF

c     scale by input real constant
      IfCu = IfCu*scIFCU
				
c     set resistive losses 
      JHres = rsvx*fsc*fsc*Jcw*Jcw*IfCu*IfCu*fcond/(1.0d0-fsc)
			
c________________________________________________________________________________		 
c
c
c             BUILD ELEMENT MATRICES
c
c________________________________________________________________________________	

			
      IF (ielc(KYOP1) .eq. 0) THEN
c___________________________________________________________________________			
c       DOF = Az - Je driven cond -> keyopt(1)=0
c___________________________________________________________________________
		
c     shape the full 3D node location into 2D 			
      xyz2D = xyz(1:edim,:)  
	
c      add a boolean to check if any matrix generation is desired -> otherwise skip integration look 
c      amatreq = .false.
			
c initialize the matrices to be formed (before sum), and sEND back 0's for other matrices IF they are requested	
      IF(kelreq(1).eq.1)THEN
        Kaa = 0.0d0
        Kai = 0.0d0
        Kee = 0.0d0
        Kei = 0.0d0
        zs = 0.0d0 
        kelout(1) = 1
c        amatreq = .true.
      ENDIF
      IF(kelreq(2).eq.1)THEN
         zass = 0.0d0
c         kelout(2) = 1     removed because mreuse takes it out in plane53
c        amatreq = .true.
      ENDIF
      IF(kelreq(3).eq.1)THEN
         Caa = 0.0d0
         Cea = 0.0d0
         damp = 0.0d0
         kelout(3) = 1
c        amatreq = .true.
      ENDIF
      IF(kelreq(4).eq.1)THEN
         gstIF = 0.0d0
         kelout(4) = 1
c        amatreq = .true.
      ENDIF				 		
      IF(kelreq(5).eq.1)THEN
         zsc = 0.0d0
         kelout(5) = 1
c        amatreq = .true.
      ENDIF
      IF(kelreq(6).eq.1)THEN
         zscnr = 0.0d0
         kelout(6) = 1
c        amatreq = .true.
      ENDIF
			
c ___________________________________________________________________________________________________________	
c integration loop: do all at once to avoid duplication operations since many are shared (jac,curl-curl,etc.) 		
c ___________________________________________________________________________________________________________	
      
      DO 801 i = 1,nq
          HL = H(:,i)		
          HL2 = RESHAPE(H(:,i),(/8,1/))
          dHL = RESHAPE(dH(:,:,i),(/8,2/))
          jac = matmul(xyz2D,dHL)
          CALL INVF2(2,jac,ijac,djac)	
          dN = matmul(dhL,ijac)
					
c         find element area
          elvol = elvol + qw(i)*djac		
					
          DO 340 k = 1,nnode									
						curldN(1,k) = dN(k,2)
						curldN(2,k) = -dN(k,1)			
  340   	CONTINUE	
	
c        set shape fucntion operations/forms that will be re-used by multiple matrices
          curlcurl = matmul(transpose(curldN),curldN)
          NtN = matmul(HL2,transpose(HL2))			
          DO 345 k = 1,8		
						DO 344 l = 1,8						
							Npat(k,l) = HL2(l,1)
  344   	CONTINUE
  345   	CONTINUE
	
c __________________________________________________________________________________________________________			
c IF requested, make stIFfness matrix by integrating with global coords using the four quadrature points		
c __________________________________________________________________________________________________________		
      IF(kelreq(1).eq.1) THEN   				
        Kaa = Kaa + qw(i)*djac*curlcurl/mu					
      ENDIF					

c __________________________________________________________________________________________________________			
c IF requested, make load vector by integrating with global coords using the four quadrature points		
c __________________________________________________________________________________________________________		
      IF(kelreq(5).eq.1) THEN   				
        zsc = zsc + qw(i)*djac*Je*HL					
      ENDIF	
			
c __________________________________________________________________________________________________________			
c IF requested, make damping matrix by integrating with global coords using the four quadrature points		
c __________________________________________________________________________________________________________		
      IF(kelreq(3).eq.1) THEN   
        Caa = Caa + 2.0d0*tau*tauMult*fcond*qw(i)*djac*curlcurl/mu
      ENDIF			
						
  801   	CONTINUE			
	
      zs = Kaa
      damp = Caa
				
c __________________________________________________________________________________________________________			
c Find Restoring Force (NROPT: NL-matrices)		
c __________________________________________________________________________________________________________

c     restoring force 
      IF(kelreq(6).eq.1) THEN 
          ui = u(:,1)
          zscnr = matmul(zs,ui)
      ENDIF
			
			
      ELSEIF (ielc(KYOP1) .eq. 1) THEN
c___________________________________________________________________________
c       DOF = Az,volt - struct with eddy -> keyopt(1)=1
c___________________________________________________________________________
c      placeholder for future development 

      ELSEIF (ielc(KYOP1) .eq. 2) THEN
c___________________________________________________________________________
c       DOF = Az,curr,emf - circuit coupled cond -> keyopt(1)=2
c___________________________________________________________________________

c     shape the full 3D node location into 2D 			
      xyz2D = xyz(1:edim,:)  	

c initialize the matrices to be formed (before sum), and send back 0's for other matrices IF they are requested	
      IF(kelreq(1).eq.1)THEN
        Kaa = 0.0d0
        Kai = 0.0d0
        Kee = 0.0d0
        Kei = 0.0d0
        zs = 0.0d0 
        kelout(1) = 1
      ENDIF
      IF(kelreq(2).eq.1)THEN
         zass = 0.0d0
c         kelout(2) = 1     removed because mreuse takes it out in plane53
      ENDIF
      IF(kelreq(3).eq.1)THEN
         Caa = 0.0d0
         Cea = 0.0d0
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
				
c ___________________________________________________________________________________________________________	
c integration loop: do all at once to avoid duplication operations since many are shared (jac,curl-curl,etc.) 		
c ___________________________________________________________________________________________________________	

      DO 901 i = 1,nq
          HL = H(:,i)		
          HL2 = RESHAPE(H(:,i),(/8,1/))
          dHL = RESHAPE(dH(:,:,i),(/8,2/))
          jac = matmul(xyz2D,dHL)
          CALL INVF2(2,jac,ijac,djac)	
          dN = matmul(dhL,ijac)
					
c         find element area
          elvol = elvol + qw(i)*djac		
					
          DO 440 k = 1,nnode									
						curldN(1,k) = dN(k,2)
						curldN(2,k) = -dN(k,1)			
  440   	CONTINUE	
	
c        set shape fucntion operations/forms that will be re-used by multiple matrices
          curlcurl = matmul(transpose(curldN),curldN)
          NtN = matmul(HL2,transpose(HL2))			
          DO 445 k = 1,8		
						DO 444 l = 1,8						
							Npat(k,l) = HL2(l,1)
  444   	CONTINUE
  445   	CONTINUE
	
c __________________________________________________________________________________________________________			
c IF requested, make stIFfness matrix by integrating with global coords using the four quadrature points		
c __________________________________________________________________________________________________________		
      IF(kelreq(1).eq.1) THEN 					
        Kaa = Kaa + qw(i)*djac*curlcurl/mu
        Kai = Kai - dirz*Nc*qw(i)*djac*NtN/Sc						
c       now add in the dIFferent lengths for inductance and resistance
        Kee = Kee - qw(i)*djac*Npat/(Sc*Lc)			
        Kei = Kei + IFCU*Nc*Nc*rsvx*qw(i)*djac*Npat/(Sc*Sc*f3)		
      ENDIF					
		
c __________________________________________________________________________________________________________			
c IF requested, make damping matrix by integrating with global coords using the four quadrature points		
c __________________________________________________________________________________________________________		
      IF(kelreq(3).eq.1) THEN   				
        Caa = Caa + 2.0d0*tau*tauMult*fcond*qw(i)*djac*curlcurl/mu				
c       now add in the different lengths for inductance and resistance
        Cea = Cea + dirz*crat*Nc*qw(i)*djac*Npat/Sc
      ENDIF			
						
  901   	CONTINUE				

c __________________________________________________________________________________________________________			
c build full matrices with padding etc		
c __________________________________________________________________________________________________________
			
c     stiffness - expand damping matrix to 16x16 following plane53 
      IF(kelreq(1).eq.1) THEN 
				l = 1
        DO 449 i = 1,nr,3
				  j = 1
          DO 448 k = 1,nr,3						
						zs(i,k) = Kaa(l,j)
						zs(i,k+2) = Kai(l,j)		
						zs(i+1,k+1) = Kee(l,j)
						zs(i+1,k+2) = Kei(l,j)
						j = j+1
  448   	CONTINUE	
						l = l+1	
  449   CONTINUE	
			
      ENDIF
			
c     damping - expand damping matrix to 16x16 following plane53 
      IF(kelreq(3).eq.1) THEN 
				l = 1
        DO 459 i = 1,nr,3
				  j = 1
          DO 458 k = 1,nr,3						
						damp(i,k) = Caa(l,j)
						damp(i+1,k) = Cea(l,j)		
						j = j+1
  458   	CONTINUE	
						l = l+1	
  459   CONTINUE			
      ENDIF
			
c     the expansion can also be done using matrix multiplication
c       -> see FEM derivation reference manual
c       -> at some point should check if this improves speed

c __________________________________________________________________________________________________________			
c Find Restoring Force (NROPT: NL-matrices)		
c __________________________________________________________________________________________________________

c     restoring force 
      IF(kelreq(6).eq.1) THEN 
          ui = u(:,1)
          zscnr = matmul(zs,ui)
      ENDIF

			
      ENDIF

c __________________________________________________________________________________________________________			
c If requested, calculate field at integration points	
c __________________________________________________________________________________________________________
      elflux = 0.0d0								
      IF (outkey.eq.1) THEN
        JHtau = 0.0d0
        Bmod = 0.0d0
        Sene = 0.0d0		
        Mx = 0.0d0
        My = 0.0d0				
        DO 852 i = 1,nq
          dHL = RESHAPE(dH(:,:,i),(/8,2/))
          jac = matmul(xyz2D,dhL)
          CALL INVF2(2,jac,ijac,djac)		
          dN = matmul(dhL,ijac)	
										
c         could include Nc,Sc,and dirz here or have them be added by user in /post1 (for now assume user will take care of it to calc inductance)					
          elflux = elflux + qw(i)*djac*DOT_PRODUCT(H(:,i),u(:,1))
					
          DO 850 k = 1,nnode				
						curldN(1,k) = dN(k,2)
						curldN(2,k) = -dN(k,1)
  850   	CONTINUE			

          IF (ielc(KYOP1) .eq. 0) THEN			
            Bm(:,i) = matmul(curldN,u(:,1))	
						dBdt(:,i) = matmul(curldN,u(:,4))          
          ELSEIF (ielc(KYOP1) .eq. 2) THEN	
					
            DO 851 k = 1,8				
						  uB(k,1) = u(1+(k-1)*3,1)
						  udBdt(k,1) = u(1+(k-1)*3,4)
  851   	  CONTINUE	
            Bm(:,i) = matmul(curldN,uB(:,1))
            dBdt(:,i) = matmul(curldN,udBdt(:,1))	
						
          ENDIF
					
c         Bmod = average of |B| at integration points f
          Bmod=Bmod+ sqrt(Bm(1,i)*Bm(1,i)+Bm(2,i)*Bm(2,i))/nq	
					
c         find element stored energy integral - no M at this point
          Sene=Sene+.5d0*(Bm(1,i)*Bm(1,i)+Bm(2,i)*Bm(2,i))*qw(i)*djac/mu				

c         IFCC heat generation and magnetization
          JHcon = 2.0d0*tau*tauMult*fcond/mu	
          JHtau=JHtau+
     x     JHcon*(dBdt(1,i)*dBdt(1,i)+dBdt(2,i)*dBdt(2,i))/nq	
          Mx = Mx - 2.0d0*tau*tauMult*fcond*dBdt(1,i)/(mu*nq)
          My = My - 2.0d0*tau*tauMult*fcond*dBdt(2,i)/(mu*nq)	
					
  852   CONTINUE					 	


					
      ENDIF	
			
c ___________________________________________________________
c --- output to database
c ___________________________________________________________
					
      IF (outkey.eq.1) THEN
				
c       see EOMASK Definitions in outpcm.inc
        IF (btest(eomask,PDBE)) THEN
				

c         FIELD FLUXES	
          IF (btest(eomask,PDBFLX)) THEN		
					
c           EXTRAPOLATE FROM GAUSS POINTS TO CORNER NODES			
 						Bx = matmul(HH,RESHAPE(Bm(1,:),(/4,1/)))
 						By = matmul(HH,RESHAPE(Bm(2,:),(/4,1/)))	
						
c           reformat 2D to 3D (also matrix to a single vector)	
 						perm = (/3,2,1,4/)   !permute to get guass points aligned with nodal connect. (swap 1->3)

						k=1					
						DO 360 i = 1,nq
							Bout(k) = Bx(perm(i),1)    !Bx
							Bout(k+1) = By(perm(i),1)  !By
							Bout(k+2) = 0.0d0          !Bz	
							
c             H = B/mu-M
							Hout(k) = Bx(perm(i),1)/mu - Mx     !Hx
							Hout(k+1) = By(perm(i),1)/mu - My   !Hy
							Hout(k+2) = 0.0d0                   !Hz	
				
							k = k+3								

  360   CONTINUE	
					
						CALL eldwrtL (elem,EDEFX,lcerstL,edindxL(1),12,Bout(1))   !write B-flux to results file
c						CALL eldwrtL (elem,EDEFX,lcerstL,edindxL(1),24,Bout(1))   !write E-flux and B-flux to results file
						CALL eldwrtL (elem,5,lcerstL,edindxL(1),12,Hout(1))   !write H-flux to results file
						
          ENDIF		

c         JOULE HEATING
          IF (btest(eomask,16)) THEN		
            JHeat(1) = JHres + JHtau
						CALL eldwrtL (elem,16,lcerstL,edindxL(1),1,JHeat(1))   
          ENDIF		
					
c         NMISC Variables
          IF (btest(eomask,13)) THEN	

						Nmisc = 0d0
			      Nmisc(1) = tempev
			      Nmisc(2) = Bev
						Nmisc(3) = I0
						Nmisc(4) = Jcw						
						Nmisc(5) = Jc						
						Nmisc(6) = qflag					
						Nmisc(7) = Lc						
						Nmisc(8) = Li						
						Nmisc(9) = scIFCU
						Nmisc(10) = tau
						Nmisc(11) = JHres
						Nmisc(12) = JHtau
						Nmisc(13) = JHeat(1)
						Nmisc(14) = IfCu
						Nmisc(15) = RRR
						Nmisc(16) = fsc
						Nmisc(17) = fcond
						Nmisc(18) = rsvx
						Nmisc(19) = elem
						Nmisc(20) = elflux
						Nmisc(21) = Mx
						Nmisc(22) = My
c						Nmisc(23) = sqrt(dBdt(1,i)*dBdt(1,i)+dBdt(2,i)*dBdt(2,i))
            Nmisc(23) = 0.0d0
c           output Cern "rhoht"
c						Nmisc(24) = rsvx*IfCu*IfCu/(fcond*(1-fsc))
c           now remove IFCU
						Nmisc(24) = rsvx/(fcond*(1-fsc))
						Nmisc(25) = Jc0
						Nmisc(26) = Tc0
						Nmisc(27) = Bc0
						Nmisc(28) = ap
						Nmisc(29) = p						
						Nmisc(30) = tauMult					
				
c           also set saved variable now
						msvr(1) = Bmod
						
						CALL eldwrtL (elem,13,lcerstL,edindxL(1),30,Nmisc(1))   !write NMISC to res file
						
c -------------------------------------------------------------
c     Set shared data for thermal (Bev,qflag)
c       -- replaces binary file transfer with ~5x speed up
c -------------------------------------------------------------
            IF (ielc(KYOP2) .eq. 1) THEN		
              MFS(1,ERpos) = Bev
              MFS(2,ERpos) = qflag
				      msvr(2) = ERpos		
            ENDIF		

       ENDIF	
	
c         SMISC 			
          IF (btest(eomask,1)) THEN		
            Smisc = Sene
            Smisc(1) = Je
						Smisc(2) = elflux
						CALL eldwrtL (elem,1,lcerstL,edindxL(1),5,Smisc(1))   
          ENDIF

c         this seems to be output by default with the elvol and elener(5) variables
          elener(1) = Sene
          elener(6) = Sene	
					
c         VOLU/AREA and Energy			
          IF (btest(eomask,4)) THEN						
						Egeo(1:5) = elener
						Egeo(1) = elvol
						CALL eldwrtL (elem,4,lcerstL,edindxL(1),11,Egeo(1))   
          ENDIF	
					
        ENDIF			
				
      ENDIF
			
c -------------------------------------------------------------
c     Output saved variables and update pointers
c -------------------------------------------------------------			

c      write out user saved element variables
      CALL svrput (svindx(1),2,nfsvr,msvr(1))    
	
c     write out the svr index vector, which updates locsvrL pointer to next element
      CALL svpidx (locsvrL,svindx(1))	
			    
					
      END SUBROUTINE uel102
			
	