c___________________________________________________________________________
c     --- USER101: 2D 8-node thermal solid (based on PLANE77)
c     --- uec101.f, element characteristics
c     --- distribution R1.0, Lucas Brouwer 11/15/2018
c___________________________________________________________________________

      SUBROUTINE uec101 (elcdn,ielc,kerr)
			
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


      CALL nminfo (ielc(1),'USER101 ')
      elcdn = 'PLANE55 clone               '
	
c     - see elccmt.inc for labels of what each ielc is
c     - see echprm.inc for relationship between wording and positions in the ielc array 
	
      ielc(ISHAP) = 4
      ielc(MNODE) = 2
			
      ielc(KELDOF) = 524288				
      ielc(KDOFS) = 524288

c     number of documented (NMDRLC) and total (NMTRLC) real constants			
      ielc(NMDRLC) = 6
      ielc(NMTRLC) = 6
			

      ielc(MATRQD) = 16
 
      ielc(MATRXS) = 10
			
      ielc(NMNDMX) = 8    
      ielc(NMNDMN) = 4    
      ielc(NMNDST) = 8    
      ielc(NMNDAC) = 8    
      ielc(NMNDNE) = 8    
			ielc(NMNDTP) = 8    
			ielc(NMDFPN) = 1    
			
      ielc(INSS) = 4
      ielc(INLD) = 4

c  surface loading info -> try to modIFy these since the element record length 
c has issue with sfelem.f
			
c      ielc(NMPTSF) = 2
      ielc(NMPTSF) = 0

c      ielc(NMCONV) = 4
      ielc(NMCONV) = 0
			
c      ielc(NMFSIN) = 4
      ielc(NMFSIN) = 0

c      ielc(NMRDSF) = 4
      ielc(NMRDSF) = 0
					
      ielc(NMHTGN) = 8
      ielc(NMBFFVI) = 8

c     saved variables -> Emap + Prev Loadstep # + nmisc-pass (Bev,etc.)
c      ielc(NMNSVR) = 2
      ielc(NMNSVR) = 7				
			
			ielc(NMNDNO) = 4
      ielc(KCONIT) = 1
      ielc(KNORM) = 1
			
			ielc(NMNMIS) = 24
c			ielc(NMNMIS) = 0

c      ielc(NMNMUP) = 24
      ielc(NMNMUP) = 0
			
      ielc(JSTPR) = 1
						
      ielc(NCPTM1) = 119
      ielc(NCPTM2) = 107
      ielc(KSWTTS) = 183		
      ielc(NOEMAT) = 1
			
					
      RETURN
      END
