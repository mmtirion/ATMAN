	subroutine RotStiffSC(angle,type1,type2,stiffness,energy)

c Coded 3/2015 by mmt.  Exact reproduction of Levitt JMB83 paper introducing L79.
c This computes the energy associated with each SC dihedral in kcal/mol.
c input xi angle  in radians. Parameters obtained from Levitt's JMB83

	implicit none
	character*1 type1,type2,hold1
	real*8 angle
	real*8 pi,K,n,del
	integer check
	real*8 stiffness,energy


	stiffness=0.d0
	K=0.d0
	n=0.d0
	del=0.d0
	pi = dacos(-1.d0)
	check = 0

	if(type1.eq.'C')then
		check=1
		goto25
	endif
	if(type2.eq.'C')then
		hold1=type1
		type1='C'
		type2=hold1
		check=1
		goto25
	endif

	if(type1.eq.'A')then
		check=2
		goto25
	endif
	if(type2.eq.'A')then
		hold1=type1
		type1='A'
		type2=hold1
		check=2
		goto25
	endif

25	continue
	IF(check.eq.1)THEN

	if(type2.eq.'C')then
		K=1.4d0
		n=3.d0
	endif
	if(type2.eq.'A'.or.type2.eq.'M')return
	if(type2.eq.'B')then
		K=0.1d0
		n=6.d0
	endif
	if(type2.eq.'S')then
		K=0.1d0
		n=3.d0
	endif
	
	ENDIF

	IF(check.eq.2)THEN

	if(type2.eq.'M')return
	if(type2.eq.'A')then
		K=20.d0
		n=2.0d0
		del=pi
	endif
	if(type2.eq.'N')then
		K=10.d0
		n=2.d0
		del=pi
	endif
	
	ENDIF

	if(check.eq.0)then
	write(*,*)'UNIDENTIFIED SIDE CHAIN DIHEDRAL'
	write(*,*)'between atom sorts:',type1,type2
	stop
	endif

	stiffness = -n*n*K*(dcos(n*angle+del))
	
	energy = K*(1.d0+dcos(n*angle+del))

	return
	end

