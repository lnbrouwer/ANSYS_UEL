c___________________________________________________________________________
c     --- USER102: 2D 8-node magnetic solid (based on PLANE53)
c     --- uec102.f, element characteristics
c     --- distribution R1.0, Lucas Brouwer 11/15/2018
c___________________________________________________________________________

      SUBROUTINE uec102 (elcdn,ielc,kerr)
			
c___________________________________________________________________________
c     --- Include Decks
c       allows for use of ANSYS definitions and data on distibution
c       see .inc at \ANSYS Inc\v171\ansys\customize\include
c___________________________________________________________________________

#include "impcom.inc"
#include "echprm.inc"

      INTEGER, INTENT(inout) :: ielc(IELCSZ)
      INTEGER, INTENT(inout) :: kerr
      CHARACTER, INTENT(inout) :: elcdn(28)

      EXTERNAL nminfo
      CALL nminfo (ielc(1),'USER102 ')
      elcdn = 'PLANE53 clone               '

c___________________________________________________________________________			
c       IndepENDent of keyopt
c___________________________________________________________________________

c     INDEPENDENT OF KEYOPT
        ielc(ISHAP) = 4
        ielc(MNODE) = 2
        ielc(KELSTO) = 1

c       add in number of real constants -> 18 total
        ielc(NMDRLC) = 18
        ielc(NMTRLC) = 18

        ielc(MATRQD) = 200				
        ielc(KMAGC) = 1
        ielc(MATRXS) = 10

        ielc(NMNDMX) = 8    
        ielc(NMNDMN) = 4    
        ielc(NMNDST) = 8    
        ielc(NMNDAC) = 8    
        ielc(NMNDNE) = 8    
        ielc(NMNDTP) = 8 

        ielc(INSS) = 1
        ielc(INLD) = 1

        ielc(NMTEMP) = 8
        ielc(NMVTDP) = 8
        ielc(NMBFFVI) = 8
				
c       surface loading info -> try to modIFy these since the element record length 
c       has issue with sfelem.f
c        ielc(NMPTSF) = 2
        ielc(NMPTSF) = 0
c        ielc(NMFGSF) = 4
        ielc(NMFGSF) = 0
c        ielc(NMFSIN) = 4
        ielc(NMFSIN) = 0
                
        ielc(NMFSVR) = 3
				
        ielc(NMNDNO) = 4
        ielc(KCONIT) = 1

c       summable misc				
        ielc(NMSMIS) = 5
				
c       non-summable misc	
        ielc(NMNMIS) = 30
				
        ielc(NMNMUP) = 24
				
        ielc(JSTPR) = 1

        ielc(NCPTM1) = 166
        ielc(NCPTM2) = 166
        ielc(KSWTTS) = 77		
				
      IF (ielc(KYOP1) .eq. 0) THEN
c___________________________________________________________________________			
c       DOF = Az - Je driven cond
c___________________________________________________________________________

        ielc(KELDOF) = 256				
        ielc(KDOFS) = 256
        ielc(NMDFPN) = 1 	
        ielc(NMCRDN) = 8
			
      ELSEIF (ielc(KYOP1) .eq. 1) THEN
c___________________________________________________________________________
c       DOF = Az,volt - struct with eddy
c___________________________________________________________________________
c     in the future an element with magnetoresisitivity for stutural eddy 
c     can be added [replicate plane53 keyopt(1)=1 with material prop fits]

      ELSEIF (ielc(KYOP1) .eq. 2) THEN
c___________________________________________________________________________
c       DOF = Az,curr,emf - circuit coupled cond
c___________________________________________________________________________

        ielc(KELDOF) = 50331904				
        ielc(KDOFS) = 50331904
        ielc(NMDFPN) = 3    
        ielc(NMCRDN) = 0
					
      ENDIF

      RETURN
      END
