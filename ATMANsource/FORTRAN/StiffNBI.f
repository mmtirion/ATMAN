	subroutine StiffNBI(xyz,Natom,i,j,C,factor,type,vdWcut,Rhc)
c
c see notes dated 10/20/2015
c

	implicit none
	character*1 type(Natom),type1
	integer i,j,Natom
	real*8 xyz(3,Natom)
	real*8 dist,epsi,epsj,SS,factor,C
	real*8 Rmin,Rrange,vdWcut,Rhc(Natom),ratio
	real*8 num,Rstar,x

	dist = dsqrt ( (xyz(1,i) - xyz(1,j))**2 +
     1	 	       (xyz(2,i) - xyz(2,j))**2 +
     2		       (xyz(3,i) - xyz(3,j))**2 )

	factor = C
	type1 = type(i)
     	call epsilons(type1,epsi)
	type1 = type(j)
    	call epsilons(type1,epsj)

     	factor = factor*12.d0*dsqrt(epsi*epsj)
    	Rmin = dsqrt(Rhc(i)*Rhc(j))
	x = dist
c if stretching  x, use the following: ! note added 3/8/16:
c (note that somewhere along the way, vdWcut = 0.0 became required...)
c (also note that somewhere along way, C = 1.0 became the "norm", and
c only "enhance" in animations have an adjustable amplitude option)
c	if(dist.le.Rmin)then
c	 x = dist
c	else
c	 Rstar = 1.10868341797d0*Rmin
c	 Rrange = Rstar + vdWcut
c	 num = (dist-Rmin)*(Rstar-Rmin)
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
c 	 if(dist.le.2.2d0)factor=SS*factor
c 	endif


	return
	end


