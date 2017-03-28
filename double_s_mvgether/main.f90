program main
	
!!	i make the init solvents ad thf at init_coor.f90,
!!	then change part of solvents as methanol with a ratio of 0.875 at the replica.f90;
!	
	use globalPara
	use random
	use init_coor
	use caculate_e
	use replica
!	use simAnnealing

	implicit none
	integer test
	
	iseed = -191
	e0 = 0.0
	open(unit=3,file="atoms.txt")
 	open(unit=4,file="chains.txt")
	open(unit=7,file="annealing.txt")
 	open(unit = 8 ,file = "graftedCh.txt")
 	open(unit = 9,file = "freeCh.txt")
	open(unit = 12,file = "solvents.cc1")
	open(unit = 13, file ="test_info.txt")

	call coor_ini()
	do 5 i = 1,lx
		nct = i*ly*lz
		write(13,*)nct,box(nct,1),box(nct,2),box(nct,3)
5 	continue

	call cacul_e()
! 	print*, 'the front of programma before this position is correct'
	print*,'type_change before,et',et
! 	call monteC(3,0)
! 	call cacul_e()
! 	print*,'et,e0',et,e0 
	call type_change()
	call cacul_e()
	e0  = et
	print*,'type_change before,et',et
	print*,'et,e0',et,e0 
! 	print*,'good solvent moving finished'

	call monteC(2,3,4)
 	do 2 i=1,numOftch
     	if ( i .le.numOfgCh ) then
     		write(4,*)i,(nchain(i,j),j=1,lenOfgch)
     		write(8,*)i,(nchain(i,j),j=1,lenOfgch)
     	else
      		write(4,*)i,(nchain(i,j),j=1,lenOfft)
      		write(9,*)i,(nchain(i,j),j=1,lenOfft)
      	end if
 2	 continue

 	do 3 i  = 1,nLattp
 		write(3,*) box(i,1),box(i,2),box(i,3),icha(i)
 3 	continue
	
 	nct = 0
 	do 4 i= 1,nLattp
 		if(icha(i).eq.3) then
 			nct = nct +1
 			write(12,100)'Na', nct,box(i,1),box(i,2),box(i,3),110	
 		else if(icha(i) .eq.4)then
 			nct = nct +1
 			write(12,100)'Cl',nct,box(i,1),box(i,2),box(i,3),190
 		endif
 4 continue

 100    format(A3,1x,I8,3(1x,I4),1x,I4,1x,I5)
	stop "end"

end program main
