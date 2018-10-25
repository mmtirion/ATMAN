	subroutine Enrgy0(xyz0,xyz,Natom,i,j,C,energy,type,vdWcut,Rhc)

	implicit none
	character*1 type(Natom),type1
	integer i,j,Natom
	real*8 xyz(3,Natom),xyz0(3,Natom)
	real*8 dist0,dist,diff,energy
	real*8 epsi,epsj,SS,factor,C
	real*8 vdWcut,Rhc(Natom)
	real*8 Rmin,Rrange,ratio

	dist0 = dsqrt ( (xyz0(1,i) - xyz0(1,j))**2 +
     1	 	       (xyz0(2,i) - xyz0(2,j))**2 +
     2		       (xyz0(3,i) - xyz0(3,j))**2 )

	factor = C
	type1 = type(i)
     	call epsilons(type1,epsi)
	type1 = type(j)
    	call epsilons(type1,epsj)

     	factor = factor*dsqrt(epsi*epsj)/(dist0*dist0)
    	Rmin = dsqrt(Rhc(i)*Rhc(j))
c      	Rrange = (Rhc(i) + Rhc(j))*0.6222d0 + vdWcut
	ratio = Rmin/dist0
	ratio = ratio*ratio
	ratio = ratio * ratio *ratio
	factor = 12.d0*factor*(13.*ratio*ratio -7.d0*ratio)
c    	factor = factor * (dist0-Rrange)*(dist0-Rrange)
c    	factor = factor / ((Rmin - Rrange)*(Rmin-Rrange))

c choose a relative strengh of disulfide bond, SS, compared to typical NCB
c	SS = 10.0d0
c	if(type(i).eq.'S'.and.type(j).eq.'S')then
c	 if(dist0.le.2.5d0)factor=SS*factor
c	endif

	dist = dsqrt ( (xyz(1,i) - xyz(1,j))**2 +
     1	 	       (xyz(2,i) - xyz(2,j))**2 +
     2		       (xyz(3,i) - xyz(3,j))**2 )


	diff = dist - dist0
 	energy = 0.5d0*factor*diff*diff

	return
	end



	subroutine Enrgy_inter(xyzP,xyzm,Nprtn,j,rLig0,rLig,
     1	NligNA,i,C,energy,type,vdWHet,vdWcut,Rhc,epsi)

	implicit none
	integer i,j,Nprtn,NligNA
	character*1 type
	real*8 xyzP(3,Nprtn),xyzm(3,Nprtn)
	real*8 rLig0(3,NligNA),rLig(3,NligNA)
	real*8 dist,rij0(3),R(3),dot,energy,factor
	real*8 C,vdWcut,vdWHet,Rhc,epsi,epsj
	real*8 Rmin,Rrange


	dist = dsqrt ( (rLig0(1,i) - xyzP(1,j))**2 +
     1	 	       (rLig0(2,i) - xyzP(2,j))**2 +
     2		       (rLig0(3,i) - xyzP(3,j))**2 )

	factor = C
     	call epsilons(type,epsj)
     	factor = factor*dsqrt(epsi*epsj)/(dist*dist)
    	Rmin = (vdWHet + Rhc)*0.5612d0
      	Rrange = (vdWHet + Rhc)*0.6222d0 + vdWcut
    	factor = factor * (dist-Rrange)*(dist-Rrange)
    	factor = factor / ((Rmin - Rrange)*(Rmin-Rrange))

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



