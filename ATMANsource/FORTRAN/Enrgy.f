	subroutine Enrgy(xyz0,xyz,Natom,i,j,C,energy,type,vdWcut,Rhc)
c
c see notes dated 10/20/2015
c

	implicit none
	character*1 type(Natom),type1
	integer i,j,Natom
	real*8 xyz(3,Natom),xyz0(3,Natom)
	real*8 dist0,dist,diff,energy
	real*8 epsi,epsj,SS,factor,C
	real*8 Rmin,Rrange,vdWcut,Rhc(Natom),ratio
	real*8 num,Rstar,x

	dist0 = dsqrt ( (xyz0(1,i) - xyz0(1,j))**2 +
     1	 	       (xyz0(2,i) - xyz0(2,j))**2 +
     2		       (xyz0(3,i) - xyz0(3,j))**2 )

	factor = C
	type1 = type(i)
     	call epsilons(type1,epsi)
	type1 = type(j)
    	call epsilons(type1,epsj)

     	factor = factor*12.d0*dsqrt(epsi*epsj)
    	Rmin = dsqrt(Rhc(i)*Rhc(j))
	x = dist0
c if stretching  x, use the following: ! note added 3/8/16:
c (note that somewhere along the way, vdWcut = 0.0 became required...)
c (also note that somewhere along way, C = 1.0 became the "norm", and
c only "enhance" in animations have an adjustable amplitude option)
c	if(dist0.le.Rmin)then
c	 x = dist0
c	else
c	 Rstar = 1.10868341797d0*Rmin
c	 Rrange = Rstar + vdWcut
c	 num = (dist0-Rmin)*(Rstar-Rmin)
c	 x = Rmin + num/(Rrange-Rmin)
c	endif
	factor = factor / (x*x)
	ratio = (Rmin/x)
	ratio = ratio*ratio
	ratio = ratio*ratio*ratio
	factor = factor*(13.d0*ratio*ratio - 7.d0*ratio)


c choose a relative strengh of disulfide bond, SS, compared to typical NCB
c to go much faster, block next four lines when there are NO DISULFIDE bonds!!
c 	SS = 10.0d0
c 	if(type(i).eq.'S'.and.type(j).eq.'S')then
c 	 if(dist0.le.2.2d0)factor=SS*factor
c 	endif

	dist = dsqrt ( (xyz(1,i) - xyz(1,j))**2 +
     1	 	       (xyz(2,i) - xyz(2,j))**2 +
     2		       (xyz(3,i) - xyz(3,j))**2 )

	diff = dist - dist0

	energy = 0.5d0*factor*diff*diff

	return
	end



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
	real*8 dist,rij0(3),R(3),dot,energy,factor
	real*8 C,vdWcut,vdWHet,Rhc,epsi,epsj
	real*8 Rmin,Rrange,Rstar,num,ratio,x


	dist = dsqrt ( (rLig0(1,i) - xyzP(1,j))**2 +
     1	 	       (rLig0(2,i) - xyzP(2,j))**2 +
     2		       (rLig0(3,i) - xyzP(3,j))**2 )

	factor = C
     	call epsilons(type,epsj)
        factor = factor*12.d0*dsqrt(epsi*epsj)
    	Rmin = dsqrt(vdWHet*Rhc)
        x = dist
c        if(dist.le.Rmin)then
c         x = dist
c        else
c         Rstar = 1.10868341797d0*Rmin
c         Rrange = Rstar + vdWcut
c         num = (dist-Rmin)*(Rstar-Rmin)
c         x = Rmin + num/(Rrange-Rmin)
c        endif
        factor = factor / (x*x)
        ratio = (Rmin/x)
        ratio = ratio*ratio
        ratio = ratio*ratio*ratio
        factor = factor*(13.d0*ratio*ratio - 7.d0*ratio)

	rij0(1) = ( rLig0(1,i)-xyzP(1,j) )/dist
	rij0(2) = ( rLig0(2,i)-xyzP(2,j) )/dist
	rij0(3) = ( rLig0(3,i)-xyzP(3,j) )/dist

	R(1) = (rLig(1,i)-xyzm(1,j)) - (rLig0(1,i)-xyzP(1,j))
	R(2) = (rLig(2,i)-xyzm(2,j)) - (rLig0(2,i)-xyzP(2,j))
	R(3) = (rLig(3,i)-xyzm(3,j)) - (rLig0(3,i)-xyzP(3,j))

	dot = rij0(1)*R(1) + rij0(2)*R(2) + rij0(3)*R(3)
	energy = factor*dot*dot

	return
	end



