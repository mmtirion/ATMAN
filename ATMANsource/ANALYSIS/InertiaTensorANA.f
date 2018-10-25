	subroutine InertiaTensor(mass,xyz,Nprtn,Nres,Nlig,Rcm,evectors)


	implicit none

	integer Nprtn,Nres,Nlig
	integer i,j,N,INFO,LDZ,Mfound,icount,IL,IU,ITYPE
	character*1 JOBZ,RANGE,UPLO,IFAIL
	integer IWORK(5*3)
	real*8 mass(Nprtn),xyz(3,Nprtn)
	real*8 massT,Rcm(3)
	real*8 Evalues(3),Evectors(3,3)
	real*8 IT(3,3),Unit(3,3),ABSTOL
	real*8 A(6),B(6),WORK(8*3),VL,VU
	real*8 x,y,z
	integer rlst
	character*4 rsrt
	character*5 asrt
	character*6 atom


	do i = 1,3
	Rcm(1) = 0.d0
	Rcm(2) = 0.d0
	Rcm(3) = 0.d0
	do j = 1,3
	 IT(i,j) = 0.d0
	enddo
	enddo

c do everything wrt center of mass:
	massT = 0.d0
	do i = 1,Nprtn
	 massT = massT + mass(i)
	 Rcm(1) = Rcm(1) + mass(i)*xyz(1,i)
	 Rcm(2) = Rcm(2) + mass(i)*xyz(2,i)
	 Rcm(3) = Rcm(3) + mass(i)*xyz(3,i)
	enddo
	Rcm(1) = Rcm(1)/massT
	Rcm(2) = Rcm(2)/massT
	Rcm(3) = Rcm(3)/massT

c translate xyz to CoM coordinates:
	do i = 1,Nprtn
	 xyz(1,i) = xyz(1,i) - Rcm(1)
	 xyz(2,i) = xyz(2,i) - Rcm(2)
	 xyz(3,i) = xyz(3,i) - Rcm(3)
	enddo


	do i = 1,Nprtn
	 IT(1,1)= IT(1,1) + 
     1		(xyz(2,i)*xyz(2,i)+xyz(3,i)*xyz(3,i))*mass(i)
	 IT(2,2)= IT(2,2) +
     1		(xyz(1,i)*xyz(1,i)+xyz(3,i)*xyz(3,i))*mass(i)
	 IT(3,3)= IT(3,3) +
     1		(xyz(1,i)*xyz(1,i)+xyz(2,i)*xyz(2,i))*mass(i)
	 IT(1,2) = IT(1,2) - xyz(1,i)*xyz(2,i)*mass(i)
	 IT(1,3) = IT(1,3) - xyz(1,i)*xyz(3,i)*mass(i)
	 IT(2,3) = IT(2,3) - xyz(2,i)*xyz(3,i)*mass(i)
	enddo

	IT(2,1) = IT(1,2)
	IT(3,1) = IT(1,3)
	IT(3,2) = IT(2,3)


c since i only have subroutine DSPGVX from LAPACK, i will introduce unit matrix for B:
	do i = 1,3
	 Unit(i,i) = 1.d0
	do j = i+1,3
	 Unit(i,j) = 0.d0
	enddo
	enddo

c ITYPE is 1 for A*x = (lambda)*B*x
	ITYPE = 1
	JOBZ = 'V'
	RANGE = 'I'
	UPLO = 'L'
	N = 3
	
	icount = 0
	do i = 1,3
	do j = i,3
	 icount = icount+1
	 A(icount) = IT(i,j)
	 B(icount) = Unit(i,j)
	enddo
	enddo

	VL = 0.d0
	VU = 0.d0

	IL = 1
	IU = 3

	ABSTOL = 1.d-06

	LDZ = 3

        CALL DSPGVX(ITYPE,JOBZ,RANGE,UPLO,N,A,B,VL,VU,
     $   IL,IU,ABSTOL,Mfound,Evalues,Evectors,LDZ,WORK,IWORK,
     $               IFAIL,INFO )


	return
	end
