module caculate_e
	use globalPara
	implicit none

	contains
	subroutine cacul_e()  		!caculate the configuration energy;

		et=0     
		do 2001 i=1,nLattp
			do 2002 j=1,nOfNeiP
	!			if ( (box(i,3).gt.0).and.(box(i,3).lt.(lz-1)) ) then
			        ii=icha(idOfNearP(i,j))
			        if(icha(i).le.ii)then
			        	et=et+eab(icha(i),ii)
			        else
			        endif

	2002	continue
	2001 continue

	end subroutine cacul_e

end module caculate_e
