module replica
	use globalPara
	use random
	use attchange_snake 
	use change_part
	use caculate_e
	implicit none
contains

	
	subroutine monteC(mv_sol_cha1,mv_sol_cha2)

		integer, intent(in) ::mv_sol_cha1,mv_sol_cha2
		integer nstar, nend, iistend
		integer iik,njj,njj1
		double precision simAnnFactor,simAnnTemp 
		integer,parameter:: maxAttMo = 80 !*(numOffCh*lenOfft)+(numOfgCh*lenOfgch)
		integer simAnnSteps ,std_mcsu
		double precision,parameter::ann_coef = 0.94
		
		e0 = et
		simAnnTemp = 60.0
		write(*,*)'e0',e0,'et',et
		write(*,*)'entry the monomer motion process'
		write(*,*)


		if ( mv_sol_cha1 .eq. 3 ) then
			simAnnSteps  = 1
			std_mcsu = 8 !000
		else if ( mv_sol_cha1 .eq. 4 ) then
			simAnnSteps  = 80
			std_mcsu = 8 !000
		end if
		
		do 3001 n = 1,simAnnSteps

			simAnnFactor =1.0/ simAnnTemp
			do 3002 kkk = 1,std_mcsu
		        attMo = 0
				accepted_MCS = 0	
				do while(attMo .lt.maxAttMo)
					ncc = 1
					nnx = int(1.0+numOftCh*ranf())  !	find a chain randomly
					if (nnx.le.numOfgCh) then       ! judge the chain is grafted chain or free chain
						nxx=int(2.0+(lenOfgch-1)*ranf())  ! choose  a monomer randomly
						lenOfMoveC= lenOfgch    ! the moving chain's length is equal to grafted chain's；
					else
						nxx=int(1.0+lenOfft*ranf())
						lenOfMoveC = lenOfft   ! the moving chain's length is equal to freechain's；
					endif

					nx=nchain(nnx,nxx)            
 					nx1=int(1.0+nOfNeiP*ranf())       
 					nx2=idOfNearP(nx,nx1)
									
! 					do 3004 i =1, lenOfMoveC 
! 						if ( nchain(nnx,i).ne.nx2 ) then
! 							kk = 1
							
! 						else
! 							nx1 = int(1.0+nOfNeiP*ranf())
! 							nx2 = identifier(nx,nx1)
! 							kk = 0
! 						end if
! 3004 				continue
! 					write(*,*)'kk = ',kk 				


					if((icha(nx2).eq. mv_sol_cha1) .or.(icha(nx2).eq.mv_sol_cha2)) then		
						ichg = 1 			!ichg stand for the Number of chain disconnection, ichg = 1 ,chain is connected ,ichg = 2 , 1 broken ;
						if(nxx.eq.lenOfMoveC)then  ! footer
						    nstar=nxx-1
						    nend=nstar
						elseif(nxx.eq.1) then	!header
							nstar=nxx+1
							nend=nstar
						else    				!middle part
					        nstar=nxx-1
					        nend=nxx+1
						endif
						iik=0
						do 3003 ii = nstar,nend,2
							iistend = 0
							nni=nchain(nnx,ii) !the lattice point identifier of nstar/nend position in chain nnx 
							kk = 1
							do while((iistend.eq.0).and.(kk .le. nOfNeiP))
								if ( idOfNearP(nni,kk) .eq.nx2 )then
								!the chain is consecutive
								iik = iik +1
								iistend = 1
								else
									kk  = kk +1
								endif
							enddo
							if ( iistend .eq. 0 )then 
								! the chain has broken
				            	ichg=ichg+1
				            	nii = ii !nii is the order number of nstar and nend in a chain
				            	njj1 = nni  ! njj1
				            else
				            endif
3003 					continue
						njj = nx ! make the value of njj  beening equal to  identifier for  nx point
						de = 0.0
						if ( ichg .eq. 1 )then
							call attChange()

!!!!						end  call attChange()
						!	write(*,*) 'ichg = 1'
						else if ( ichg .eq.2)then      !!!! refer to 'if ( ichg .eq. 1 )then'
							
							if (nii .lt. nxx) then
								mm = -1
							else
								mm = 1
							endif

							ncc = 1
							
							do while(( ichg .eq.2 ) .and.(nii .ge.2).and.(nii .le. lenOfMoveC-1))
								nii=nii+mm
			 			       	nni=nchain(nnx,nii)
			                   	iii=0
			                   	kk=1
			                   	do while((ichg .ne.1).and.(kk .le.nOfNeiP))
			                   		if ( idOfNearP(nni,kk) .eq. njj )then
			                   			ichg = 1
			                   		else
			                   			kk = kk+1
			                   		endif
			                   	enddo

			                   	njj = njj1
			                   	njj1 = nni
			                   	ncc = ncc +1
			                enddo

			                if ( nii .eq. 1 )then
			                	if ( ichg .eq.2 ) then
			                		if ( nnx .le.numOfgCh )then
			                		!	write(*,*)'break' !,nni,nii,mm,ncc
			                			goto 9999
			                		else
			                			ncc = ncc+1
			                			ichg = 1
			                		endif
			                	else
			                	endif
			                else if ( nii .eq. lenOfMoveC )then
			                	if ( ichg .eq.2 ) then
			                		ncc = ncc +1
			                		ichg = 1
			                	else
			                	endif
			                else
			                endif

			            ! 						try tp snake motion  if the ichg == 2
				            
							call snake()
							!!!  end call snake()
				!			write(*,*) 'a snake'
! 							end snake motion

			            else 	!!!!	refer to 'if ( ichg .eq. 1 )then'  line 82
						endif	!!!!	refer to 'if ( ichg .eq. 1 )then'	line 82	
					else   		!!!!!  refer to 'if (icha(nx2).eq.3) then'
					endif		!!!!!  refer to 'if (icha(nx2).eq.3) then'		

9999        attMo = attMo +ncc
				
			enddo
3002		continue
		simAnnTemp =simAnnTemp*ann_coef
			write(*,*)'attMo-',attMo,'accepted_MCS-',accepted_MCS
			write(7,*)'attMo-',attMo,'accepted_MCS-',accepted_MCS
			call cacul_e()
			write(*,*)'n=',n,'et,e0,de',et,e0,de
			write(7,*)'n=',n,'et,e0,de',et,e0,de
			
3001 	continue
		
!		write(*,*)'enter check_snake'
		call check_snake()
!		write(*,*)'quit check_snake'

	end subroutine monteC




!	line 103 a goto!!!! be not allowed 

end module replica
	
