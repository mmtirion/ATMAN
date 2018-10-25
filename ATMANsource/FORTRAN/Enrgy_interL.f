	subroutine Enrgy_interL(rLig0,rLig,j,k,factor,energy)

c  added 12/2016 to handle RB NBI energies in ATMAN


	implicit none
	integer j,k
	real*8 rLig0(3,2500),rLig(3,2500)
	real*8 dist0,dist,diff,energy,factor


	dist0 = dsqrt ((rLig0(1,j) - rLig0(1,k))**2 +
     1	 	       (rLig0(2,j) - rLig0(2,k))**2 +
     2		       (rLig0(3,j) - rLig0(3,k))**2 )


	dist = dsqrt ( (rLig(1,j) - rLig(1,k))**2 +
     1	 	       (rLig(2,j) - rLig(2,k))**2 +
     2		       (rLig(3,j) - rLig(3,k))**2 )

	diff = dist - dist0

	energy = 0.5d0*factor*diff*diff

	return
	end



