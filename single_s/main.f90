program main

	use globalPara
	use random
	use init_coor
	use caculate_e
	use replica
!	use simAnnealing

	implicit none
	integer test,nempty3
	
	iseed = -191
	e0 = 0.0

	open(unit=3,file="atoms.txt")
	open(unit=4,file="chains.txt")
	open(unit=7,file="annealing.txt")
	open(unit = 8 ,file = "graftedCh.txt")
	open(unit = 9,file = "freeCh.txt")
	open(unit = 12,file = "solvents.cc1")

	call coor_ini()
	call cacul_e()
	write(*,*)'the front of programma before this position is correct'

	call monteC(3,0)
	call cacul_e()
	write(*,*)'et,e0',et,e0 
	nct = 0
	do 1 i = 1, nLattp
		if ( icha(i) .eq. 3 ) then
			nct = nct +1
			empty3(nct) = i
		end if

1  continue
	
	nempty3 = nct
	write(*,*)'nempty3 =',nct,nOftype3
	nct = 0
	nOftype4 =int( nOftype3 *0.875)
	do while (nct .lt. nOftype4)
		nx1 =int(1+ nOftype3*ranf())
		nx = empty3(nx1)
		icha(nx) = 4
		ichap(nx) = 4
		nct = nct +1
	enddo
	write(*,*)'nOftype4=',nct,nOftype4
	call cacul_e()
	e0  = et 
	write(*,*)'et,e0',et,e0 
	write(*,*)'good solvent moving finished'

	call monteC(4,0)
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
