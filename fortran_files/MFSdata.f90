!___________________________________________________________________________
!     --- module for shared EM/THERMAL properties and passed loads
!     --- this replaces binary file writing used in previous versions
!     --- distribution R1.0, Lucas Brouwer 11/15/2018
!___________________________________________________________________________


MODULE MFSdata

IMPLICIT NONE

!     SUBROUTINEs		
  PUBLIC initialize_MFS_EM
	PUBLIC initialize_MFS_TH
  PUBLIC initialize_var, count_ne
  PUBLIC check_nc, check_rc
	
!     shared data
  DOUBLE PRECISION, PUBLIC, ALLOCATABLE, SAVE :: MFS(:,:)	
  DOUBLE PRECISION, PUBLIC, ALLOCATABLE, SAVE :: RCvar(:,:)	
	
	
PRIVATE
	INTEGER ne,npos
	LOGICAL isDeclared
  DOUBLE PRECISION, ALLOCATABLE :: EMmap(:)
	
	
CONTAINS 

!_____________________________________________________________________________________________
!     --- intitialize module variables (set counters etc.)
!_____________________________________________________________________________________________

	SUBROUTINE initialize_var()
		ne = 0
		isDeclared = .false.
		RETURN
	END SUBROUTINE initialize_var
		
	
!___________________________________________________________________________
!     --- count number of elements with passing on first iteration (or input as real const)
!___________________________________________________________________________

	SUBROUTINE count_ne()
		ne = ne + 1
		RETURN
	END SUBROUTINE count_ne

!___________________________________________________________________________
!     --- count number of elements with passing on first iteration (or input as real const)
!___________________________________________________________________________

	SUBROUTINE check_nc(neout)
		INTEGER neout
		neout = ne
		RETURN
	END SUBROUTINE check_nc

!___________________________________________________________________________
!     --- count number of elements with passing on first iteration (or input as real const)
!___________________________________________________________________________

	SUBROUTINE check_rc(neout,nnpos)
		INTEGER neout,nnpos
		neout = RCvar(1,nnpos)
		RETURN
	END SUBROUTINE check_rc
	

	
!___________________________________________________________________________________
!     --- set up initial mapping for EM
!___________________________________________________________________________________

	SUBROUTINE initialize_MFS_EM(nEM,nposEM,RRR,fcond,fsc)
		INTEGER nEM,nposEM
    DOUBLE PRECISION RRR,fcond,fsc
		
!   initialize the arrays based on the number of elements to be mapped
!   store the order of the EM elements in the array EMmap (and RETURN to uel)
!    this will be used to set values and also searched by the thermal element 
!    to find the mapped element 

		IF (isDeclared) THEN	
		
!     set up mapping
			nposEM = npos
			EMmap(npos) = nEM
			RCvar(1,npos) = RRR
			RCvar(2,npos) = fcond
			RCvar(3,npos) = fsc
			npos = npos + 1
			
		ELSE
		
			npos=1
			nposEM = npos
			
!     Bev, qflag  (transfered during Multi-field Solve)
		  ALLOCATE(MFS(2,ne))
!     RRR, fcond, fsc (only transfered 1x)			
		  ALLOCATE(RCvar(3,ne))
!     order of EM elements in arrays (to be matched with thermal)
		  ALLOCATE(EMmap(ne))		
			isDeclared = .true.

			EMmap(npos) = nEM
			RCvar(1,npos) = RRR
			RCvar(2,npos) = fcond
			RCvar(3,npos) = fsc
			npos = npos + 1


    ENDIF	
		
		RETURN
	END SUBROUTINE initialize_MFS_EM
	
	
!___________________________________________________________________________________
!     --- set up initial mapping for TH
!___________________________________________________________________________________

	SUBROUTINE initialize_MFS_TH(nEMmap,nposTH)
		INTEGER nEMmap,nposTH,i

!   search the EM map to find out the position of the thermal element in the transfer
!   RETURN this value (nposTH) to the uel 

		DO i = 1,ne
				IF (EMmap(i) .EQ. nEMmap) THEN	
					nposTH = i
					EXIT
				ENDIF	
		ENDDO
		
		RETURN
	END SUBROUTINE initialize_MFS_TH


	
END MODULE MFSdata




