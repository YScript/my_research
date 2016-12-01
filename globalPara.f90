module globalPara

	implicit none
	integer,parameter:: lx=60, ly=60, lz=60
	integer,parameter:: lxyz = lx*ly*lz
	integer,parameter:: nOfNeiP=18
	real,parameter:: lenOfNeiP = 2 ! , numNP=18
	real rr
	integer i,j,k,ii,iii,jj,kk,kkk,l,m,n
	integer ncountNP
	double precision  e0, et, de 
	
	integer,dimension(216000,3)::box
	integer,dimension(216000)::icha,ichap
	double precision,dimension(5,5)::eab
	integer,allocatable,dimension(:,:):: nchain,nearP_coor,idOfNearP
	integer,allocatable,dimension(:)::nc,nearP,a,empty3
	integer nOftype3,nOftype4
 
!*************************************************************************************
!	the following parameters is used in initialization module
!************************************************************************************
	integer nct,nLattp,nAtom
	integer x,y,z
	integer testa,err
	integer nx,nx1,nx2  !*
	integer,parameter::ndx = 3
	integer,parameter:: numOffCh = 600
	integer,parameter:: lenOffas = 11,lenOfft = 17,lenOfgch= 11
	integer,parameter:: numOfgCh = (lx*ly)/(ndx*ndx),numOftCh=numOffCh+numOfgCh
!**************************************************************************************
!	the following parameter , replication module;
	integer nnx,nxx,nii,nni,ncc,ncc1,ncc2,mm,nccp  !***
	integer accepted_MCS,attMo
	integer ichg,mmp,lenOfMoveC,mv_cha
!	double precision pp
!	integer nearPofa ,nearPofb, ichaNPa,ichaNPb,tmpichaNPa,tmpichaNPb
!	integer ichanx,ichanx1,ichanx2,ichanccp,tmpichanccp
	double precision simAnnFactor,simAnnTemp 
	integer,dimension(:)::ichapp(lenOfft+1) 
!attemp_snake:  nx,nx2,ncc2,nnx,nxx,mm,ncc,icha(),ichap(),ichapp(),
!change_part:  nx,nx2,accepted_MCS,ncc,nchain(),nnx,nxx,nccp,mm,ncc1,ncc2 ~~must be global

!*************************************************************************************
!	the following parameters is used in random number generator function ,random module
!************************************************************************************
	INTEGER, PARAMETER :: K4B=selected_int_kind(9)
    INTEGER(K4B) iseed
    INTEGER(K4B), PARAMETER :: IA=16807, IM=2147483647, IQ=127773, IR=2836
    REAL, SAVE :: am
    INTEGER(K4B), SAVE :: ix=-1, iy=-1, kx
		

end module globalPara



!README
!
!
! lx,ly,lz					the simulation system box length of each orentation;
!i,j,k,ii,jj,kk,l,m,n,ll,mm,nneap
!							the counter number often used in loop;


!ranf_new					the random number return by module function;
!e0 						the initial system energy ;
!et							the current system energy ;
!de 						energy difference between old and new system config;
!lenOfNeiP 					length of neighbour points;
				
!nOfNeiP					number of neighbour points;

!array .box 				box lattice points 
!array .atom				the particles'coordinates in system;
!array .icha ,ichap 		atom character array ,and temporary array;
!array .eab 				interaction array of each pairs;
!array .nchain(m,n) 		chains array in simulation, m chains with n monomers per chain;
!ichapp 					temp array of monomer in per chain ,when chains moving;
!nc
!mmstar
!nearP(i,j)					neighbour(near) points array ,i atoms ,with j near points per point coordinate;
!nearP_coor(i,j)			coordinate of atom i 's jth neighbour points';
!							and the j is equal to numNP

!lenOffas					length of free chain, A segments
!lenOfft 					length of free chain, total length
!lenOfgch					length of grafted chain 


!nct 						used to counter ,initialize the value as 1;
!nLattp 					number of lattice points; conclude the top and bottum substrates;
!nAtom 						number of atoms between top abd bottum substrates;


!mv_sol_cha 					the character of choose  solvent