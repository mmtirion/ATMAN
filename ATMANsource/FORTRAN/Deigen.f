	subroutine Deigen(A,B,Nentry,N,Evalues,Evectors,Nmodes)

	implicit none

	character*1 jobz,range,uplo

	integer*4 il,info,itype,iu,ldz,Mfound,N,Nentry
	real*8 abstol, vl, vu
	real*8 eps,pi

	integer*4 ifail(N),iwork(5*N)
	real*8 A(Nentry),B(Nentry),WORK(8*N)
	real*8 Evalues(N),Evectors(N,N)

	integer i,j,Nmodes


	write(*,*)' '
	write(*,*)' try to open data file...'
c        open(55,status='unknown',file=
c     1'/Users/monique/Desktop/WORK/NMStuff/ATMAN/DATA/eigenvectors')
        open(55,status='unknown',file='../PDB/eigenvectors')


	pi = dacos(-1.d0)

	itype = 1  ! specifies problem type as: F x = lambda H x

	jobz = 'V'   ! compute eigenvalues AND eigenvectors
c       	jobz = 'N'   ! compute eigenvalues 


  	range ='I'   ! the ILth thru IUth eigenvalues will be found
c      	range ='A'   ! all eigenvalues will be found

	uplo = 'L' ! lower triangle of symmetric matrices are stored

c	n ! given as input parameter in subroutine call
c	AP and BP are given as input matrices A, B in subroutine call

	vl = 0.0d0
	vu = 0.0d0 ! eigen range to be computed: not used!

	il = 1
  	iu = Nmodes     ! eigenvalue numbers to be computed
c      	iu = N

	abstol = 0.0d0
	ldz =  N


	CALL DSPGVX(ITYPE,JOBZ,RANGE,UPLO,N,A,B,VL,VU,
     $   IL,IU,ABSTOL,Mfound,Evalues,Evectors,LDZ,WORK,IWORK,
     $               IFAIL,INFO )

	write(*,*)' For diagonalization:'
	write(*,*)' '
	write(*,*)' INFO = ',INFO

	write(*,*)' # of eigenvalues found: ', Mfound
	write(*,*)' '
	IF(info.eq.0)THEN

	 write(*,*)' Eigenvalues are:'
54	 format(d12.6)
	 do i = il, iu
	if(Evalues(i).le.0.d0)then
	  write(55,*)i,Evalues(i)
	  write(*,*)i,Evalues(i),'raw'
	else
	  Evalues(i) = dsqrt(Evalues(i))*2.0454828d13   ! in s^-1
	  write(55,*)i,Evalues(i)  !in s^-1
	  Evalues(i) = Evalues(i)/(2.d0*pi)  !frequency, not angular frequency
	  write(*,*)i,Evalues(i)/2.99792458d10   !'in cm^-1' 
	endif
	if(jobz.eq.'V')then
	 do j = 1, N
	  write(55,*)2.454d23*Evectors(j,i)    !1/sqrt(kg)/m
	 enddo
	endif
	 enddo
	ENDIF
	write(*,*)'... '
	write(*,*)' EIGENVALUES & EIGENVECTORS recorded'
	write(*,*)'... '

	close(55)

	return
	end
