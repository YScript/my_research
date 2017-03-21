 

 
	   integer, allocatable, dimension(:,:) ::  nna, nchain
	   REAL, allocatable, dimension(:,:) :: atom
	   integer, allocatable, dimension(:) :: nn, icha, mc, mw,ichap, nempty, mne, nnd0

!/////////////////////复制粘贴区域///////////////////////////////////////////////////////////////////



	nnd_A=6
	nnd1=nnd_A
    nnd_B=6
	nnd_di=nnd_A+nnd_B
	nnd=nnd_A+nnd_B
	
	ntotc1=900
	ntotc=1900
	ntotc2=ntotc-ntotc1

 
	nnn=26
	 
	  LXxx=60
	  LYyy=LXxx
 	  LZzz=LXxx
	  lx=lxxx
	  ly=lyyy
	  lz=lzzz
	 XLXxx=FLOAT(LXxx)
     XLYyy=FLOAT(LYyy)
     XLZzz=FLOAT(LZzz)
     XLXxxH=XLXxx/2
     XLYyyH=XLYYY/2
     XLZzzH=XLZZZ/2
	  	  ntot=LXxx*LYyy*LZzz
 
  
 
      
	   allocate(atom(ntot,3),nna(ntot,26),nchain(ntotc,nnd),nn(ntot),icha(ntot),mc(ntot),mw(ntot),&
	            ichap(ntot),nempty(ntot),mne(ntot),nnd0(ntotc))

!/////////////////////复制粘贴区域///////////////////////////////////////////////////////////////////
 
 
 
	do i=1,ntotc
	 
	    	nnd0(i)=nnd_di
		 
	enddo

	open(1,file='d.dat',status='old')
		do  kkk=1,ntotc
			nnd=nnd0(kkk)
			read(1,*)i,(nchain(kkk,j),j=1,NND)
		enddo
	close(1)



 
		atom(1,1)=0.
		atom(1,2)=0.
	  	atom(1,3)=0.0
	 	
		   kkk=0
		 do  i1=1,LXxx

    	 do  i2=1,LYyy 
	 
	   	 do  i3=1,LZzz
			  kkk=kkk+1
	      	  atom(kkk,1)=atom(1,1)+(i1-1)			
	          atom(kkk,2)=atom(1,2)+(i2-1)
			  atom(kkk,3)=atom(1,3)+(i3-1)
		  enddo
		  enddo
		  enddo
  

	!//////////////////////  Rg/////////////////////////////////////////////////////////////

     !////////////// diblock  Rg   接枝链的Rg  ///////////

			RcmAx_di = 0
			RcmAY_di = 0
			RcmAZ_di = 0
			RcmBx_di = 0
			RcmBY_di = 0
			RcmBZ_di = 0
			RcmABx_di = 0
			RcmABY_di = 0
			RcmABZ_di = 0
		Rg2A_di = 0
		Rg2B_di = 0
		Rg2AB_di = 0

	do i=1,ntotc1
 			llx=lxxx
			lly=lyyy
			llz=lzzz

		kkk=1
		do while(kkk>0)
		kx=0
		ky=0
		kz=0
		do j=2,nnd_di
 		x1=atom( nchain(i,j) ,  1)
		y1=atom( nchain(i,j) ,  2)
		z1=atom( nchain(i,j) ,  3)
		x0=atom( nchain(i,j-1) ,  1) 
		y0=atom( nchain(i,j-1) ,  2) 
		z0=atom( nchain(i,j-1) ,  3) 
		if(x1 > llx ) x1=x1-lxxx
		if(y1 > lly ) y1=y1-lyyy
		if(z1 > llz ) z1=z1-lzzz
		if(x0 > llx ) x0=x0-lxxx
		if(y0 > lly ) y0=y0-lyyy
		if(z0 > llz ) z0=z0-lzzz
			rrx=x1-x0
			rry=y1-y0
			rrz=z1-z0
			
			if(abs(rrx) > XLXxxH) kx=kx+1 !数有多少个断点（即周期性过去的断点）
			if(abs(rry) > XLXxxH) ky=ky+1 !数有多少个断点（即周期性过去的断点）
			if(abs(rrz) > XLXxxH) kz=kz+1 !数有多少个断点（即周期性过去的断点）

		enddo

		if(kx /=0 ) llx=llx-1
		if(ky /=0 ) lly=lly-1
		if(kz /=0 ) llz=llz-1

		
		if(kx /=0 .or. ky /=0 .or. kz /=0 ) kkk=1
		if(kx ==0 .and. ky ==0 .and. kz ==0 ) kkk=0
		enddo
		write(*,*)i,llx,lly,llz
		if(llz<lzzz) then
			write(*,*)'llz<lzzz'
			pause
		endif




			RcmAx_di=0
			RcmAY_di=0
			RcmAZ_di=0	
			RcmBx_di=0
			RcmBY_di=0
			RcmBZ_di=0
			RcmABx_di=0
			RcmABY_di=0
			RcmABZ_di=0
		do j=1,nnd_A  
			x0=atom( nchain(i,j) ,  1) 
			y0=atom( nchain(i,j) ,  2) 
			z0=atom( nchain(i,j) ,  3) 
			if(x0 > llx ) x0=x0-lxxx
			if(y0 > lly ) y0=y0-lyyy
			if(z0 > llz ) z0=z0-lzzz
			RcmAx_di = RcmAx_di + x0 ! 将原点（0,0,0）设为不动点
			RcmAY_di = RcmAY_di + y0
			RcmAZ_di = RcmAZ_di + z0
		enddo
		do j=nnd_A+1,nnd_di
			x0=atom( nchain(i,j) ,  1) 
			y0=atom( nchain(i,j) ,  2) 
			z0=atom( nchain(i,j) ,  3) 
			if(x0 > llx ) x0=x0-lxxx
			if(y0 > lly ) y0=y0-lyyy
			if(z0 > llz ) z0=z0-lzzz
			RcmBx_di = RcmBx_di + x0 ! 将原点（0,0,0）设为不动点
			RcmBY_di = RcmBY_di + y0
			RcmBZ_di = RcmBZ_di + z0
		enddo
		do j=1,nnd_di
			x0=atom( nchain(i,j) ,  1) 
			y0=atom( nchain(i,j) ,  2) 
			z0=atom( nchain(i,j) ,  3) 
			if(x0 > llx ) x0=x0-lxxx
			if(y0 > lly ) y0=y0-lyyy
			if(z0 > llz ) z0=z0-lzzz
			RcmABx_di = RcmABx_di + x0 ! 将原点（0,0,0）设为不动点
			RcmABY_di = RcmABY_di + y0
			RcmABZ_di = RcmABZ_di + z0
		enddo
			RcmAx_di = RcmAx_di / nnd_A
			RcmAY_di = RcmAY_di / nnd_A
			RcmAZ_di = RcmAZ_di / nnd_A
			RcmBx_di = RcmBx_di / nnd_B
			RcmBY_di = RcmBY_di / nnd_B
			RcmBZ_di = RcmBZ_di / nnd_B
			RcmABx_di = RcmABx_di / nnd_di
			RcmABY_di = RcmABY_di / nnd_di
			RcmABZ_di = RcmABZ_di / nnd_di
		do j=1,nnd_A 
			x0=atom( nchain(i,j) ,  1) 
			y0=atom( nchain(i,j) ,  2) 
			z0=atom( nchain(i,j) ,  3) 
			if(x0 > llx ) x0=x0-lxxx
			if(y0 > lly ) y0=y0-lyyy
			if(z0 > llz ) z0=z0-lzzz
			Rg2A_di = Rg2A_di + ( x0 - RcmAx_di )**2 + &
								 ( y0 - RcmAy_di )**2 + &
								  ( z0 - RcmAz_di )**2 
		enddo
		do j=nnd_A+1,nnd_di
			x0=atom( nchain(i,j) ,  1) 
			y0=atom( nchain(i,j) ,  2) 
			z0=atom( nchain(i,j) ,  3) 
			if(x0 > llx ) x0=x0-lxxx
			if(y0 > lly ) y0=y0-lyyy
			if(z0 > llz ) z0=z0-lzzz
			Rg2B_di = Rg2B_di + ( x0 - RcmBx_di )**2 + &
								 ( y0 - RcmBy_di )**2 + &
								  ( z0- RcmBz_di )**2 
		enddo
		DO	j=1,nnd_di
			x0=atom( nchain(i,j) ,  1) 
			y0=atom( nchain(i,j) ,  2) 
			z0=atom( nchain(i,j) ,  3) 
			if(x0 > llx ) x0=x0-lxxx
			if(y0 > lly ) y0=y0-lyyy
			if(z0 > llz ) z0=z0-lzzz
			Rg2AB_di = Rg2AB_di + ( x0- RcmABx_di )**2 + &
									( y0 - RcmABy_di )**2 + &
									  ( z0 - RcmABz_di )**2 
		enddo							

    enddo
		Rg2A_di = Rg2A_di / (ntotc1*nnd_A*(nnd_A-1))
		Rg2B_di = Rg2B_di / (ntotc1*nnd_B*(nnd_B-1))
		Rg2AB_di = Rg2AB_di /(ntotc1*nnd_di*(nnd_di-1))

		Rg2A_brush = Rg2A_di 
		Rg2B_brush = Rg2B_di 
		Rg2AB_brush = Rg2AB_di 


	write(*,*) ntotc, Rg2A_brush ,  Rg2B_brush,   Rg2AB_brush
 

!////////////////////// END 接枝链的Rg////////////////////////////////////////////


   !////////////// begin diblock  Rg   自由链的Rg  ///////////

			RcmAx_di = 0
			RcmAY_di = 0
			RcmAZ_di = 0
			RcmBx_di = 0
			RcmBY_di = 0
			RcmBZ_di = 0
			RcmABx_di = 0
			RcmABY_di = 0
			RcmABZ_di = 0
		Rg2A_di = 0
		Rg2B_di = 0
		Rg2AB_di = 0

	do i=ntotc1+1,ntotc
 			llx=lxxx
			lly=lyyy
			llz=lzzz

		kkk=1
		do while(kkk>0)
		kx=0
		ky=0
		kz=0
		do j=2,nnd_di
 		x1=atom( nchain(i,j) ,  1)
		y1=atom( nchain(i,j) ,  2)
		z1=atom( nchain(i,j) ,  3)
		x0=atom( nchain(i,j-1) ,  1) 
		y0=atom( nchain(i,j-1) ,  2) 
		z0=atom( nchain(i,j-1) ,  3) 
		if(x1 > llx ) x1=x1-lxxx
		if(y1 > lly ) y1=y1-lyyy
		if(z1 > llz ) z1=z1-lzzz
		if(x0 > llx ) x0=x0-lxxx
		if(y0 > lly ) y0=y0-lyyy
		if(z0 > llz ) z0=z0-lzzz
			rrx=x1-x0
			rry=y1-y0
			rrz=z1-z0
			
			if(abs(rrx) > XLXxxH) kx=kx+1 !数有多少个断点（即周期性过去的断点）
			if(abs(rry) > XLXxxH) ky=ky+1 !数有多少个断点（即周期性过去的断点）
			if(abs(rrz) > XLXxxH) kz=kz+1 !数有多少个断点（即周期性过去的断点）

		enddo

		if(kx /=0 ) llx=llx-1
		if(ky /=0 ) lly=lly-1
		if(kz /=0 ) llz=llz-1

		
		if(kx /=0 .or. ky /=0 .or. kz /=0 ) kkk=1
		if(kx ==0 .and. ky ==0 .and. kz ==0 ) kkk=0
		enddo
		write(*,*)i,llx,lly,llz
		if(llz<lzzz) then
			write(*,*)'llz<lzzz'
			pause
		endif



			RcmAx_di=0
			RcmAY_di=0
			RcmAZ_di=0	
			RcmBx_di=0
			RcmBY_di=0
			RcmBZ_di=0
			RcmABx_di=0
			RcmABY_di=0
			RcmABZ_di=0
		do j=1,nnd_A  
			x0=atom( nchain(i,j) ,  1) 
			y0=atom( nchain(i,j) ,  2) 
			z0=atom( nchain(i,j) ,  3) 
			if(x0 > llx ) x0=x0-lxxx
			if(y0 > lly ) y0=y0-lyyy
			if(z0 > llz ) z0=z0-lzzz
			RcmAx_di = RcmAx_di + x0 ! 将原点（0,0,0）设为不动点
			RcmAY_di = RcmAY_di + y0
			RcmAZ_di = RcmAZ_di + z0
		enddo
		do j=nnd_A+1,nnd_di
			x0=atom( nchain(i,j) ,  1) 
			y0=atom( nchain(i,j) ,  2) 
			z0=atom( nchain(i,j) ,  3) 
			if(x0 > llx ) x0=x0-lxxx
			if(y0 > lly ) y0=y0-lyyy
			if(z0 > llz ) z0=z0-lzzz
			RcmBx_di = RcmBx_di + x0 ! 将原点（0,0,0）设为不动点
			RcmBY_di = RcmBY_di + y0
			RcmBZ_di = RcmBZ_di + z0
		enddo
		do j=1,nnd_di
			x0=atom( nchain(i,j) ,  1) 
			y0=atom( nchain(i,j) ,  2) 
			z0=atom( nchain(i,j) ,  3) 
			if(x0 > llx ) x0=x0-lxxx
			if(y0 > lly ) y0=y0-lyyy
			if(z0 > llz ) z0=z0-lzzz
			RcmABx_di = RcmABx_di + x0 ! 将原点（0,0,0）设为不动点
			RcmABY_di = RcmABY_di + y0
			RcmABZ_di = RcmABZ_di + z0
		enddo
			RcmAx_di = RcmAx_di / nnd_A
			RcmAY_di = RcmAY_di / nnd_A
			RcmAZ_di = RcmAZ_di / nnd_A
			RcmBx_di = RcmBx_di / nnd_B
			RcmBY_di = RcmBY_di / nnd_B
			RcmBZ_di = RcmBZ_di / nnd_B
			RcmABx_di = RcmABx_di / nnd_di
			RcmABY_di = RcmABY_di / nnd_di
			RcmABZ_di = RcmABZ_di / nnd_di
		do j=1,nnd_A 
			x0=atom( nchain(i,j) ,  1) 
			y0=atom( nchain(i,j) ,  2) 
			z0=atom( nchain(i,j) ,  3) 
			if(x0 > llx ) x0=x0-lxxx
			if(y0 > lly ) y0=y0-lyyy
			if(z0 > llz ) z0=z0-lzzz
			Rg2A_di = Rg2A_di + ( x0 - RcmAx_di )**2 + &
								 ( y0 - RcmAy_di )**2 + &
								  ( z0 - RcmAz_di )**2 
		enddo
		do j=nnd_A+1,nnd_di
			x0=atom( nchain(i,j) ,  1) 
			y0=atom( nchain(i,j) ,  2) 
			z0=atom( nchain(i,j) ,  3) 
			if(x0 > llx ) x0=x0-lxxx
			if(y0 > lly ) y0=y0-lyyy
			if(z0 > llz ) z0=z0-lzzz
			Rg2B_di = Rg2B_di + ( x0 - RcmBx_di )**2 + &
								 ( y0 - RcmBy_di )**2 + &
								  ( z0- RcmBz_di )**2 
		enddo
		DO	j=1,nnd_di
			x0=atom( nchain(i,j) ,  1) 
			y0=atom( nchain(i,j) ,  2) 
			z0=atom( nchain(i,j) ,  3) 
			if(x0 > llx ) x0=x0-lxxx
			if(y0 > lly ) y0=y0-lyyy
			if(z0 > llz ) z0=z0-lzzz
			Rg2AB_di = Rg2AB_di + ( x0- RcmABx_di )**2 + &
									( y0 - RcmABy_di )**2 + &
									  ( z0 - RcmABz_di )**2 
		enddo							

    enddo
	if(ntotc2 /= 0)then
		Rg2A_di = Rg2A_di / (ntotc2*nnd_A*(nnd_A-1))
		Rg2B_di = Rg2B_di / (ntotc2*nnd_B*(nnd_b-1))
		Rg2AB_di = Rg2AB_di /(ntotc2*nnd_di*(nnd_di))
	


		Rg2A_free = Rg2A_di 
		Rg2B_free = Rg2B_di 
		Rg2AB_free = Rg2AB_di 
	endif

	write(*,*) ntotc,   Rg2A_free ,  Rg2B_free,   Rg2AB_free
 

!////////////////////// END 自由链的Rg////////////////////////////////////////////

 

	open(1,file='brush_free_Rg.txt',access='append')
	 
       ! write(1,'(100(1x,a20))') 'ntotc','Rg2A_brush ','Rg2B_brush',' Rg2AB_brush','Rg2A_free ','Rg2B_free',' Rg2AB_free'
		write(1,'(i20,100(1x,f20.6))') ntotc, Rg2A_brush ,  Rg2B_brush,   Rg2AB_brush,  Rg2A_free ,  Rg2B_free,   Rg2AB_free
 


	write(*,*)

 END