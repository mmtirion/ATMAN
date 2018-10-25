	subroutine Enrgy_inter(xyzP,xyzm,Nprtn,j,rLig0,rLig,
     1	NligNA,i,C,energy,type,vdWHet,vdWcut,Rhc,epsi)

c On 10/2/2016 this is amended to be consistent with FormII above
c however, L79 (parent potential) does not seem optimized to model
c inter-molecular forces, but better intra-molecular forces.
c This observation needs careful study. For ex. Wolynes claims
c that intra-molecular forces mediated not by hydrophobic forces,
c but by hydrophilic forces mediated by intervening H2O molecules.
c So at the very least, the cutoff distances here need to be larger.
c But more importantly, I may need to consider including hydrogens everywhere.


	implicit none
	integer i,j,Nprtn,NligNA
	character*1 type
	real*8 xyzP(3,Nprtn),xyzm(3,Nprtn)
	real*8 rLig0(3,2500),rLig(3,2500)
	real*8 dist0,dist,rij0(3),R(3),diff,energy,factor
	real*8 C,vdWcut,vdWHet,Rhc,epsi,epsj
	real*8 Rmin,Rrange,Rstar,num,ratio,x


	dist0 = dsqrt ( (rLig0(1,i) - xyzP(1,j))**2 +
     1	 	       (rLig0(2,i) - xyzP(2,j))**2 +
     2		       (rLig0(3,i) - xyzP(3,j))**2 )

	factor = C
     	call epsilons(type,epsj)
        factor = factor*12.d0*dsqrt(epsi*epsj)
    	Rmin = dsqrt(vdWHet*Rhc)
        if(dist0.le.Rmin)then
         x = dist0
        else
         Rstar = 1.10868341797d0*Rmin
         Rrange = Rstar + vdWcut
         num = (dist0-Rmin)*(Rstar-Rmin)
         x = Rmin + num/(Rrange-Rmin)
        endif
        factor = factor / (x*x)
        ratio = (Rmin/x)
        ratio = ratio*ratio
        ratio = ratio*ratio*ratio
        factor = factor*(13.d0*ratio*ratio - 7.d0*ratio)

	dist = dsqrt ( (rLig(1,i) - xyzm(1,j))**2 +
     1	 	       (rLig(2,i) - xyzm(2,j))**2 +
     2		       (rLig(3,i) - xyzm(3,j))**2 )

	diff = dist - dist0

	energy = 0.5d0*factor*diff*diff

	return
	end



