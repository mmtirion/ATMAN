	subroutine Enrgy2(xyz0,xyz,Natom,i,j,factor,energy)
c
c dated 9/12/16
c

	implicit none
c input:
	integer i,j,Natom
	real*8 xyz(3,Natom),xyz0(3,Natom)
c computation:
	real*8 dist0,dist,diff,factor
c output:
	real*8 energy


	dist0 = dsqrt ((xyz0(1,i) - xyz0(1,j))**2 +
     1	 	       (xyz0(2,i) - xyz0(2,j))**2 +
     2		       (xyz0(3,i) - xyz0(3,j))**2 )

	dist = dsqrt ( (xyz(1,i) - xyz(1,j))**2 +
     1	 	       (xyz(2,i) - xyz(2,j))**2 +
     2		       (xyz(3,i) - xyz(3,j))**2 )

	diff = dist - dist0

	energy = 0.5d0*factor*diff*diff


	return
	end



