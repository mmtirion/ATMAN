	subroutine hessian_check(A,N)

	implicit none
	integer N,i,j,iterms(N),ntot(N)
	real*8 A(N,N)
	real*8 max,avg



	do i = 1, N
	 max = 0.0d0
	 avg = 0.0d0
	 iterms(i) = 0

	do j = i, N
	 if(dabs(A(i,j)).gt.max)max=dabs(A(i,j))
	 avg = avg + dabs(A(i,j))
	 if(A(i,j).ne.0.d0)iterms(i) = iterms(i) + 1
	enddo

c         write(*,*)'Column ',i,' average: ',avg/real(n)
	 if(iterms(i).eq.0)write(*,*)i,iterms(i)
	enddo




	return
	end
