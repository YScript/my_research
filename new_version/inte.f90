! program main

! 	implicit none
! 	interface interface_1
! 			subroutine helloWorld(i,j,k)
! 				integer, intent(inout) :: i,j,k
! ! 				integer,optional::k
! 				integer,parameter:: myKind = SELECTED_INT_KIND(8)
! 				real(kind=myKind):: name
! 				
! 			end subroutine helloWorld
				
! 			end interface ! interfa		
! 	call helloWorld(1,2,3)
	print*,'helloWorld'
! end program main
