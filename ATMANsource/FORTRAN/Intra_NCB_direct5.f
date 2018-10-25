	subroutine Intra_NCB_direct5(Nprtn,xyz,Rhc,IntAct,Iint,vdWcut,C)

c Input xyz coordinates and identify interactions, InterActA(Nint,2) within 
c a cutoff distance Rint dependant on van der Waal radii and constant determined
c by dba in "improved_potentials.pdf".
c Subroutine Four_Bonds needed to filter out NBI closer
c than four bonds.  created from Intra_NCB_direct.f on 8/30/2016 mmt

	implicit none

	integer Nprtn
	real*8 xyz(3,Nprtn)
	real*8 C,vdWcut,Rhc(Nprtn)

	integer i,j,Iint,nearesti,nearestj
	real*8 dist,Rrange
	real*8 nearest,farthest

	integer IntAct(Nprtn,Nprtn)



	Iint=0
	nearest = 100.d0
	farthest = 0.d0
	DO i = 1, Nprtn
	 IntAct(i,i) = 0
	 Do j = i+1, Nprtn
	  IntAct(i,j)=0
	  IntAct(j,i)=0
	  dist= dsqrt( (xyz(1,i)-xyz(1,j))**2 +
     1	 	       (xyz(2,i)-xyz(2,j))**2 +
     1		       (xyz(3,i)-xyz(3,j))**2 )
	  Rrange = 1.10868341797d0*C*dsqrt(Rhc(i)*Rhc(j)) + vdWcut !included on 3/7/16
	  IF(dist.le.Rrange)THEN
	   if(dist.lt.nearest)then
		nearest=dist
		nearesti=i
		nearestj=j
	   endif
	   if(dist.gt.farthest)farthest = dist
	   Iint = Iint+1
	   IntAct(i,j) = 1
	   IntAct(j,i) = 1
	  ENDIF
	 EndDo  
	ENDDO  

	write(*,*)' '
	write(*,*)'In Intra_NCB_direct NBI=',Iint
	write(*,*)' '
	write(*,*)'Nearest any two NB atoms:',sngl(nearest),'btw',nearesti,nearestj
	write(*,*)'Farthest apart 2 NB atoms are is:',sngl(farthest)
	write(*,*)' '
 

	return
	end
