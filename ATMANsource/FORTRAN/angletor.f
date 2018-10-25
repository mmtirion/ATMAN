	subroutine angletor(i1,i2,i3,i4,N,xyz,angle)

c for 4 consecutive atoms defining a dihdral stored in index,
c and at cartesian location xyz(3,index(i)), compute associated
c dihedral angle, positive if atom pair 'behind' rotated clockwise wrt front pair


	implicit none
	integer N,i1,i2,i3,i4
	real*8 xyz(3,N),angle

	real*8 A(3),B(3),C(3),Chat(3),D(3)
	real*8 V_A(3),V_D(3),Vhat_A(3),Vhat_D(3),V_AD(3)
	real*8 sina,cosb,pi


	pi = dacos(-1.d0)

	A(1) = xyz(1,i1)-xyz(1,i2)
	A(2) = xyz(2,i1)-xyz(2,i2)
	A(3) = xyz(3,i1)-xyz(3,i2)

	B(1) = xyz(1,i2)-xyz(1,i3)
	B(2) = xyz(2,i2)-xyz(2,i3)
	B(3) = xyz(3,i2)-xyz(3,i3)

	C(1) = xyz(1,i3)-xyz(1,i2)
	C(2) = xyz(2,i3)-xyz(2,i2)
	C(3) = xyz(3,i3)-xyz(3,i2)
	call AmpA(C,Chat)

	D(1) = xyz(1,i4)-xyz(1,i3)
	D(2) = xyz(2,i4)-xyz(2,i3)
	D(3) = xyz(3,i4)-xyz(3,i3)

	call AcrossB(A,B,V_A)
	call AmpA(V_A,Vhat_A)

	call AcrossB(C,D,V_D)
	call AmpA(V_D,Vhat_D)

	call AcrossB(Vhat_A,Vhat_D,V_AD)

	call AdotB(V_AD,Chat,sina)
	call AdotB(Vhat_A,Vhat_D,cosb)

	if(cosb.gt.0.d0)then

	 angle = dasin(sina)

	else
	 
	  if(sina.gt.0.d0)angle =  pi - dasin(sina)
	  if(sina.lt.0.d0)angle = -pi - dasin(sina)

	endif

	return
	end




