	program correlation

c inputs two ascii files, and computes their correlations

	implicit none

	integer i,j,k,N
	real A,B,SA,SB,SAB,SAA,SBB
	real num,den,cor

c	open(22,file=
c     1'/Users/monique/Desktop/ToTirion_1E0W_cache/1E0W_sbnma_Ca_m01.B')
	open(22,file=
     1'/Users/monique/Desktop/ML_1E0W/1e0w_H_GAP_tir.B')
	open(24,file=
     1'/Users/monique/Desktop/ML_1E0W/1e0w_nH_GAP_tir.B')
c	open(24,file=
c     1'/Users/monique/Desktop/WORK/NMStuff/ATMAN/DATA/Bfactor')


	N = 0
	SA = 0.0
	SB = 0.0
	SAB = 0.0
	SAA = 0.0
	SBB = 0.0
	

	do i = 1,10000


	  read(22,*,end=99)j,A
	  read(24,*,end=99)k,B


	  N = N + 1
	  SA = SA + A
	  SB = SB + B
	  SAB = SAB + A*B
	  SAA = SAA + A*A
	  SBB = SBB + B*B

	enddo

99	continue

	num = N*SAB - SA*SB

	den = sqrt(N*SAA-SA*SA)*sqrt(N*SBB-SB*SB)

	cor = num/den

	write(*,*)'Correlation:',cor
	

	stop
	end
