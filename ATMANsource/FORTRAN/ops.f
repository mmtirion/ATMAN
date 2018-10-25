	subroutine AmpA(A,Ahat)

	implicit none
	real*8 A(3),Ahat(3)
	real*8 mag


	mag = A(1)*A(1) + A(2)*A(2) + A(3)*A(3)
	mag = dsqrt(mag)

	Ahat(1) = A(1)/mag
	Ahat(2) = A(2)/mag
	Ahat(3) = A(3)/mag

	return
	end





	subroutine AcrossB(A,B,C)

	implicit none

	real*8 A(3),B(3),C(3)

c  A_vec  X  B_vec = C_vec

	C(1) =  A(2)*B(3) - A(3)*B(2)

	C(2) = -A(1)*B(3) + A(3)*B(1)
	
	C(3) =  A(1)*B(2) - A(2)*B(1)


	return
	end





	subroutine AdotB(A,B,C)

	implicit none

c A_vec dot B_vec = C

	real*8 A(3),B(3),C

	C = A(1)*B(1)+A(2)*B(2)+A(3)*B(3)

	return
	end


