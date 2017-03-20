    
    program main
	Dimension natom(216000,3),icha(216000),nchain(1300,12),nna(216000,26),nn(216000)
	COMMON /CSEED/ ISEED
	  open(unit=5,file='d.dat')
	  open(unit=2,file='ao.dat')
	  ndx=2
	  nnd1=6
	  nnd=12 !
      nba=0
	  nbb=0
	  nbs=0
	  nbw=0
	  nndx=nnd-nnd1
	  lx=60
	  ly=60
	  lz=60
 ntot=lx*Ly*Lz
 ntotc1=lx*Ly/(ndx*ndx)
 ntotc2=400
 ntotc=ntotc1+ntotc2
 write(*,*)'ntotc=',ntotc
 ntotA=ntotc*nndx
      xxlx=60.
	  xxly=60.
	  xxlz=60.
xxlxb2=xxlx/2.
xxlyb2=xxly/2.
	  natom(1,1)=0
	  natom(1,2)=0
	  natom(1,3)=0
	  kkk=0
	  do 1 i1=1,lx
		 do 2 i2=1,ly
			do 3 i3=1,lz
			   kkk=kkk+1
			   natom(kkk,1)=natom(1,1)+(i1-1)
			   natom(kkk,2)=natom(1,2)+(i2-1)
			   natom(kkk,3)=natom(1,3)+(i3-1)
3			 continue
2         continue
1      continue  
  do 5 i=1,ntot-1
	    j=i
	    do j=i+1,ntot
	      
	      rrx=float(natom(i,1)-natom(j,1))
	      rry=float(natom(i,2)-natom(j,2))
	      rrz=float(natom(i,3)-natom(j,3))
	   	   if(rrx.gt.xxlxb2) then
	          rrx=-xxlx+rrx
	        else if(rrx.lt.-xxlxb2) then
	          rrx=xxlx+rrx
	        else
	        end if
	        if(rry.gt.xxlyb2) then
	          rry=-xxly+rry
	        else if(rry.lt.-xxlyb2) then
	          rry=xxly+rry
	        else
	        end if
		  rr=rrx*rrx+rry*rry+rrz*rrz
	      if(rr.lt.3.5)then
	   	    nn(i)=nn(i)+1
	        nn(j)=nn(j)+1
	        nna(i,nn(i))=j
	        nna(j,nn(j))=i
	      else
	      endif
	    enddo	
5	 continue
 write(*,*)'after 5'
 do 6 i1=1,ntot
	   	 if(natom(i1,3).eq.0)then
		   icha(i1)=4
		 else if(natom(i1,3).eq.(Lz-1))then
		   icha(i1)=5
		 else
		   icha(i1)=3
		 end if
6	  continue
   	  write(*,*)'after 6'     
	do 7 kkk=1,ntotc
	 read(5,*)i,(nchain(kkk,j),j=1,nnd)
  do 119 kki=1,nnd1
   icha(nchain(kkk,kki))=1
 119 continue
 do 129 kki=1+nnd1,nnd
   icha(nchain(kkk,kki))=2
 129 continue
7 continue
write(*,*)'after 7'
  do 9 kkk=1,ntotc
   do 521 kki=nnd1+1,nnd
    q=nchain(kkk,kki)
	do 120 kk=1,nn(q)
	nx=nna(q,kk)
    if(icha(nx).eq.1)then
     nba=nba+1
	 elseif(icha(nx).eq.2)then
     nbb=nbb+1
     elseif(icha(nx).eq.3)then
     nbs=nbs+1
     elseif(icha(nx).eq.4)then
	 nbw=nbw+1
	 else  
	endif
	120 continue
   521 continue
  9 continue
write(2,*) 'nba=',float(nba)/ntotA,'nbb=',float(nbb)/ntotA,'nbs=',float(nbs)/ntotA,'nbw=',float(nbw)/ntotA
  end

	  