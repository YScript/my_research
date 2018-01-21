module change_part
	use globalPara
! 	use caculate_e
	implicit none
contains
subroutine change1()

		e0 = e0+de
		nchain(nnx,nxx)=nx2
		icha(nx)=ichap(nx)
		icha(nx2)=ichap(nx2)
		accepted_MCS=accepted_MCS+ncc

end subroutine change1

subroutine change2()
		
		e0 = e0+de
	do 3301 i=1,ncc
			nccp=ichapp(i)
			icha(nccp)=ichap(nccp)
			mmp = nxx+(i-1)*mm
		    nchain(nnx,mmp) = nccp

3301		continue
			nccp=ichapp(ncc+1)
			icha(nccp)=ichap(nccp)
			accepted_MCS =accepted_MCS+ncc
			
end subroutine change2
			

end module change_part
