module attchange_snake
	use globalPara
	use random
 	use caculate_e
	use change_part
	implicit none
contains
	subroutine attChange()		
	double precision pp
	integer ichanx,ichanx2,nearPofa ,nearPofb, ichaNPa,ichaNPb,tmpichaNPa,tmpichaNPb
			do  3101 ii = 1, nOfNeiP
			 		nearPofa = idOfNearP(nx,ii)
			 		nearPofb = idOfNearP(nx2,ii)
			 		ichanx = icha(nx)
			 		ichanx2 = icha(nx2)
			 	
			 		if ( nearPofa .ne.nx2 )then
			 			ichaNPa =icha(nearPofa)
			 			tmpichaNPa = ichap(nearPofa)
			 			de =  de +eab(ichanx2,tmpichaNPa)-eab(ichanx,ichaNPa) 
			 		else
			 		endif

			 		if ( nearPofb .ne.nx ) then
						ichaNPb = icha(nearPofb)
			 			tmpichaNPb = ichap(nearPofb) 			
			 			de = de + eab(ichanx,tmpichaNPb)-eab(ichanx2,ichaNPb)
			 		else
			 		endif
3101 		continue
			
			ichap(nx2) = icha(nx)
			ichap(nx) = icha(nx2)

			if ( de .le.0 ) then
				call change1()
			
			else
				pp = exp(-de*simAnnFactor)
				if ( pp .ge.ranf() ) then
				call change1()

				else
					ichap(nx2)=icha(nx2)
	                ichap(nx)=icha(nx)
				endif
			endif

	end subroutine attChange



	subroutine snake()
		integer jcc,jcc1,nccpp,njcc,njcc1
		integer ichanx,ichanx1,ichanx2,ichanccp,tmpichanccp
		integer nearPofa ,nearPofb, ichaNPa,ichaNPb,tmpichaNPa,tmpichaNPb
		
		double precision pp
	!	integer mmp
		ncc1  = ncc -1 
		ncc2 = ncc +1
!		write(*,*)'ncc=',ncc,'nx=',nx
		nx1 = nx
		ichanx1 = ichap(nx1)
		ichapp(1) = nx2

		do 3201 jcc = 2,ncc2
			ichapp(jcc) = nx1
			if ( jcc .ne. ncc2 ) then
				mmp = nxx+(jcc-1)*mm
!				write(*,*)'nnx=',nnx,'mmp=',mmp
				nx1 = nchain(nnx,mmp)
			else
			endif
3201 	continue
		do 3202 jcc = 1,ncc
			njcc = ichapp(jcc)
			njcc1 = ichapp(jcc+1)
			ichap(njcc) = icha(njcc1)
3202 	continue
		ichap(ichapp(ncc2)) = icha(nx2)

		do 3203 jcc = 1,ncc
			njcc=ichapp(jcc)
			njcc1=ichapp(jcc+1)
			de=de+eab(icha(njcc),icha(njcc1))-eab(ichap(njcc),ichap(njcc1)) 
3203 	continue

		do 3204 jcc = 1,ncc2
			nccp = ichapp(jcc)    	!!!!!the pair of atom coordinates
			ichanccp = icha(nccp)
			tmpichanccp = ichap(nccp)

			do 3205 ii = 1,nOfNeiP
				nearPofa  = idOfNearP(nccp,ii)
				ichaNPa = icha(nearPofa)
				tmpichaNPa = ichap(nearPofa)

				de = de +eab(tmpichanccp,tmpichaNPa)-eab(ichanccp,ichaNPa)
3205		continue

3204 	continue

		do 3206 jcc = 1,ncc
			nccp = ichapp(jcc)
			do 3207 ii = 1,nOfNeiP
				nearPofa =idOfNearP(nccp,ii)
				do 3208 kk = jcc+2,ncc2
					nccpp = ichapp(kk)
					if ( nearPofa .eq. nccpp ) then
						de = de+eab(icha(nccp),icha(nccpp))-eab(ichap(nccp),ichap(nccpp))
					else
					endif
3208			continue
3207		continue
3206 	continue
		
		if (de .le.0) then
			
			call change2()


		else
			pp = exp(-de*simAnnFactor)
			if ( pp .gt.ranf() ) then
			call change2()
					
			else
				do  3209 jcc=1,ncc2
					njcc=ichapp(jcc)
					ichap(njcc)=icha(njcc)
3209 			continue
				write(*,*)'refute the attChange moving'
			endif
		endif



end subroutine snake



subroutine check_snake()
	integer,parameter::hflx = lx/2
	integer,parameter::hfly = ly/2

	call  cacul_e()

	if ( e0 -et .gt.1.0 .or.et - e0 .gt. 1.0 ) then
		write(*,*)'error','e0',e0,'et',et
	end if



	do 3301 i = 1,numOftCh
		if ( i .le. numOfgCh ) then
			lenOfMoveC = lenOfgch
		else
			lenOfMoveC = lenOfft			
		end if

		do 3302 j = 1,(lenOfMoveC-1)
			x= box(nchain(i,j),1) - box(nchain(i,j+1),1)
			y= box(nchain(i,j),2) - box(nchain(i,j+1),2)
			z= box(nchain(i,j),3) - box(nchain(i,j+1),3)
			if(x.gt.hflx) then
	          	x=-lx+x
	        else if(x.lt.-hflx) then
	          	x=lx+x
			else
			end if

			if(y.gt.hfly) then
				y=-ly+y
			else if(y.lt.-hfly) then
				y=ly+y
			else
			end if

			rr = x*x+y*y+z*z
			if(rr.gt.lenOfNeiP)then
				write(*,*) 'error, check_snake'
!				pause
			else
			endif			
3302 	continue
3301 continue
	
end subroutine check_snake



end module attchange_snake
