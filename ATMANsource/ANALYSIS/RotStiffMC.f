	subroutine RotStiffMC(phi,psi,factorphiphi,
     1factorpsipsi,energy)

c Written 5/15 MMT. Use L79 potential (Levitt JMBiol 83) to derive
c stiffness constants, V''(NS), for dihedral angles for PDB-NMA

	implicit none

        real*8 phi,psi

	real*8 pi,conv
        real*8 Ei(4),PP0(4,2),W(4,2)
	integer i
	real*8 ai,bi,difphi,difpsi,fac,preA,preB

	real*8 factorphiphi,factorpsipsi,energy


    	Ei(1) = -4.d0
        Ei(2) = -2.d0
        Ei(3) = -4.d0
        Ei(4) = -2.d0

c Convert angles to radians!
        pi = dacos(-1.0d0)
        conv = pi/180.d0

        PP0(1,1) = -75.d0*conv
        PP0(1,2) = 0.d0

        PP0(2,1) = -60.d0*conv
        PP0(2,2) = 150.d0*conv

        PP0(3,1) = 75.d0*conv
        PP0(3,2) = 0.d0

        PP0(4,1) = 60.d0*conv
        PP0(4,2) =-150.d0*conv

        W(1,1) = 40.d0*conv
        W(1,2) = 40.d0*conv

        W(2,1) = 20.d0*conv
        W(2,2) = 30.d0*conv

        W(3,1) = 40.d0*conv
        W(3,2) = 40.d0*conv

        W(4,1) = 20.d0*conv
        W(4,2) = 30.d0*conv

	

	factorphiphi = 0.d0
	factorpsipsi = 0.d0
	energy = 0.d0

	do i = 1,4

	difphi = phi - PP0(i,1)
	difpsi = psi - PP0(i,2)

	ai = 0.693d0/dsin(W(i,1)/2.d0)
	bi = 0.693d0/dsin(W(i,2)/2.d0)

	fac = Ei(i)*dexp(-ai*(1.d0-dcos(difphi)))
	fac =   fac*dexp(-bi*(1.d0-dcos(difpsi)))

	energy = energy + fac

	preA = ai*ai*dsin(difphi)*dsin(difphi)
	preB = ai*dcos(difphi)
	factorphiphi = factorphiphi + (preA-preB)*fac

	preA = bi*bi*dsin(difpsi)*dsin(difpsi)
	preB = bi*dcos(difpsi)
	factorpsipsi = factorpsipsi + (preA-preB)*fac

	enddo

c "there is no "stretching" here,s o you should use these spring constants only for
c those phi,psi that yield a positive value":
	if(factorphiphi.lt.0.0d0)factorphiphi=0.d0
	if(factorpsipsi.lt.0.0d0)factorpsipsi=0.d0
	if(factorphiphi.lt.0.0d0)energy=0.d0
	if(factorpsipsi.lt.0.0d0)energy=0.d0

	return
	end
