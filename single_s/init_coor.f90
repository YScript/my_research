module init_coor

	use globalPara	! call the others module
	use random
	use caculate_e
	implicit none
!****************** variable  list ******************

contains

!	initialise the interaction array

subroutine coor_ini()             ! define a subroutine ,build the coordinates
	integer nn1

	do 1001 i = 1,5
 		do 1002 j = 1,5					!   initialization process ,set the each interaction be 0;
 			eab(i,j) = 0.0				! 	specific value set after that;
 1002 	continue
 1001 continue

 	eab(1,2) = 1.0
 	eab(1,3) = -1.0
	eab(1,4) =  1.0
 	eab(1,5) = 3.0
 	eab(2,3) = -1.0
	eab(2,4) = -1.0
 	eab(2,5) = -0.5

 	do 1003 i = 1,5
 		do 1004 j = 1,5
			eab(j,i) = eab(i,j)
1004	continue
1003 continue
!
!	do 1008 i = 1,5
!		do 1009 j = 1,5
!			write(*,*)i,j,eab(i,j)		
! 1009 	continue
! 1008 continue
!
!	initialise the lattice points coordinates
	nct = 0
	do 1005 i = 0,lx-1
		do 1006 j = 0,ly-1 
			do 1007 k = 0,lz-1  !  do not just select the moving atom
				nct = nct + 1
				box(nct,1) = i
				box(nct,2) = j
				box(nct,3) = k
1007 		continue
1006 	continue
1005 continue
	
	nLattp=nct
	nct = 0;
	write(*,*)'nLattp',nLattp,lxyz
	 do 1008 i =1,nLattp
		if ( (box(i,3).gt.0).and.(box(i,3).lt.(lz-1)) )then
			nct = nct +1
			icha(i) = 3
		else if ( box(i,3).eq.0 )then
			icha(i) = 5 				!! the bottum sub is 5 
		else
			icha(i) = 4					!! the top	sub is 4
		endif

1008 continue

	nAtom = nct

!	write(*,*)'now, make the near points coor'
	allocate(nearP(1:nLattp))

	do 1010 i = 1, nLattp
		nearP(i) = 0
1010 continue
!	write(*,*)'to be continue 1010'
	nct = 0
	allocate(nearP_coor(nOfNeiP,3))
! allocate(nearP_coor(nOfNeiP,3), stat=err)
! if (err /= 0) print *, "array: Allocation request denied"
! if (allocated(nearP_coor(nOfNeiP,3))) deallocate(nearP_coor(nOfNeiP,3), stat=err)
! if (err /= 0) print *, "array: Deallocation request denied"
	do 1011 i = -1,1
		do 1012 j = -1,1
			do 1013 k = -1, 1
				x= i
				y= j
				z =k
				rr= (x*x)+(y*y)+(z*z)
				if ((rr.gt.0).and.(rr.le.lenOfNeiP))then
					nct = nct + 1
!					write(*,*)nct,'test array'
					nearP_coor(nct,1) = x
					nearP_coor(nct,2) = y
					nearP_coor(nct,3) = z	
				else
				endif
1013		continue
1012	continue
1011 continue

!	write(*,*)'to be continue 1011'
	ncountNP = nct
! 	do 1014 i = 1,nct
! 		write(*,*)i,nearP_coor(i,1),nearP_coor(i,2),nearP_coor(i,3)
! 1014 continue
	allocate(idOfNearP(nLattp,ncountNP))
	do 1015 i = 1, nLattp
		do 1016 j = 1, ncountNP
		x = box(i,1) + nearP_coor(j,1)
		y = box(i,2) + nearP_coor(j,2)
		z = box(i,3) + nearP_coor(j,3)
		if ( x.gt.(lx-1) )then
			x = -lx+x
		elseif	(x.lt.0)then
			x = lx+x
		else
		endif
		if (y.gt.(ly-1) )then
			y = -ly+y
		else if(y .lt.0)then
			y = ly +y
		else
		endif
		if ( z.gt.(lz-1) )then
			z = -lz+z
		elseif(z.lt.0)then
			z = lz+z
		else
		endif
		idOfNearP(i,j) =x*ly*lz+y*lz+z+1		!the identifier is different from the front,is the box_coor identifier
												! count from (0,0,0) as the 1st atom;

1016 continue
1015 continue
!	write(*,*)'to be continue 1016 ' 
! 	do 1017 i = 1,100
! 		do 1018 j = 1,nOfNeiP
!				write(*,*)'216000,1',idOfNearP(216000,18)
! 1018	continue		
! 1017 continue
	nct = 0
	allocate(nc(lx*ly)) !a layer has lx*ly atoms;
	do 1019 i = 1,lxyz
		if ( box(i,3).eq.1)then
			nct = nct + 1
			nc(nct) = i
		endif
1019 continue
!	write(*,*)'nct=',nct
!	write(*,*)'to be continue 1019'
	allocate(nchain(numoftch,lenOfft))   !!build a [numoftch *lenOfft] array
	ii = 0
	do 1020 i = 1, numOfgCh
		if ( mod(ii,lx).eq.0 )then
			ii=ii+(ndx-1)*ly+ndx
		else
			ii=ii+ndx
		endif
		nx = nc(ii) 
!		write(*,*)'i=',i,'ii',ii,'nx=',nx,nc(60),box(205142,1),box(205142,2)
		if(icha(nx).eq.3)then
		 	icha(nx)=1
	     	nct=1
	     	nchain(i,nct)=nx
		else
		end if
		do while(nct.ne.lenOfgch)
	     	nx2=nx+1
		    if(icha(nx2).eq.3)then
	     		icha(nx2)=1
	     		nct=nct+1
	     		nchain(i,nct)=nx2
	   	 		nx=nx2
         	else
	     	end if
		end do
1020 continue
!	write(*,*)'to be continue 1020'
	nx = 0
	nct = 0
!!!!    arrange the free AB copolymer chains
	do 1021 i = numOfgCh+1,numoftch
		nn1 = 1
		do while(nn1.ne.3)
			
			nx = int(1.0+nLattp*ranf()) !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			write(*,*)'nx=',nx
			nn1 = icha(nx)
			if(nn1.eq.3)then
! 				if(atom(nx,3).ge.lenOfgch)then
	   				icha(nx)=1
	   				nct=1
					nchain(i,nct)=nx
! 				else
! 					nn1=1
! 				endif
			else
		 	endif
		enddo
		do while(nct.ne.lenOffas)		 
	      	nx1=INT(1.0+nOfNeiP*ranf())
	  		nx2= idOfNearP(nx,nx1)
	  		if(icha(nx2).eq.3)then
	  			icha(nx2)=1
		  		nct=nct+1
		 		nchain(i,nct)=nx2
				nx=nx2
			else
			endif
		enddo
		do while(nct.ne.lenOfft)		
		  	nx1=INT(1.0+nOfNeiP*ranf())
	  	    nx2= idOfNearP(nx,nx1)
			if(icha(nx2).eq.3)then
			   	icha(nx2)=2
				nct=nct+1
				nchain(i,nct)=nx2		 
				nx=nx2				
			else
			endif
		enddo
1021 continue
!	write(*,*)'to be continue 1021бо
	nct = 0 
	do 1022 i = 1,nLattp
	    if(icha(i).eq.3)then
	            nct = nct+1 
		else
	  	endif
1022  continue
     nOftype3=nct
    allocate(empty3(nOftype3))
    do 1023 i=1,nLattp    !				if the Range of values of i1 is not equal to nBox ,than the value of ichap()  
		ichap(i)=icha(i)	!would be 0 ,and array bounds exceed; when the var i of ichap(i) ge.natom;	
 
1023     continue



!!!!					test sector
! nct=0
! 	e0=0
	     
! 	do 2001 i=1,nLattp
! 		do 2002 j=1,nOfNeiP
! 			if ( (box(i,3).gt.0).and.(box(i,3).lt.(lz-1)) ) then
! 		        ii=icha(idOfNearP(i,j))
! 		        if(icha(ii).le.ii)then
! 		        	e0=e0+eab(icha(i),ii)
! 		        else
! 		        endif
! 		    else
! 		    endif
! 2002	continue
! 2001 continue
! et = e0
! write(*,*)'1e0=',e0,'et = ',et
! call cacul_e()
!!!!!					end test
end subroutine coor_ini

!subroutine coor_test use to test the correct of the function reference;

subroutine coor_test(testa)

	integer,intent(in)::testa
	allocate(a(1))
	a(1)= testa+10


	do 1024 i = 1, numoftch
		if ( i.le.numOfgCh )then
			write(4,*)i,(nchain(i,j),j=1,lenOfgch)
		else
			write(4,*)i,(nchain(i,k),k=1,lenOfft)
		endif
1024 continue

end subroutine coor_test





end module init_coor

