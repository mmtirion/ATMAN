	subroutine dRdTorsion(massP,xyzP,IndexMC,IndexXi,EndXi,Nres,
     1  Nxi,Npro,Nprtn,asort,delta,NDoF,Pyn,dRdTor)


	implicit none
c Input:
	integer Nprtn,NPro,Nres,Nxi,NDoF
	integer IndexMC(2*Nres,4),IndexXi(Nxi,4),EndXi(Nxi)
        real*8 xyzP(3,Nprtn),massP(Nprtn)
	character*5 asort(Nprtn)
c For computations:
	integer icount,itor,jtor,i,iatom,m,m2,m3
	real*8 delta
	integer ier,Pyn(2*Nres+Nxi)
        real*8 U(3,3),T(3),rms,mass(Nprtn)
	real*8 Tmass,Rcm(3),sum,dist
	real*8 nhatP(3,NDoF)
	real*8 x0,y0,z0
	real*8 xyzdP(3,Nprtn)
c output:
	real*8 dRdTor(Nprtn,NDoF,3)


c Calculate the center of mass: 
        Rcm(1) = 0.d00
        Rcm(2) = 0.d00
        Rcm(3) = 0.d00
        Tmass = 0.d00

        do iatom = 1,Nprtn
         Rcm(1) = Rcm(1) + massP(iatom)*xyzP(1,iatom) 
         Rcm(2) = Rcm(2) + massP(iatom)*xyzP(2,iatom) 
         Rcm(3) = Rcm(3) + massP(iatom)*xyzP(3,iatom) 
         Tmass = Tmass + massP(iatom)
         mass(iatom) = dsqrt(massP(iatom)) !!!!!!
        enddo

        Rcm(1) = Rcm(1)/Tmass
        Rcm(2) = Rcm(2)/Tmass
        Rcm(3) = Rcm(3)/Tmass

c Move coordinates to CoM of system:
        do iatom = 1,Nprtn    
         xyzP(1,iatom) = xyzP(1,iatom) - Rcm(1)
         xyzP(2,iatom) = xyzP(2,iatom) - Rcm(2)
         xyzP(3,iatom) = xyzP(3,iatom) - Rcm(3)
        enddo   

c  get nhat vector:
	do m = 1,2*Nres+Nxi
	 nhatP(1,m) = 0.0d00
	 nhatP(2,m) = 0.0d00
	 nhatP(3,m) = 0.0d00
	enddo

	do m = 1,2*Nres     
	 m3 = IndexMC(m,3)
	 m2 = IndexMC(m,2)
	 if(Pyn(m).eq.0)goto10
	 T(1) = xyzP(1,m3) - xyzP(1,m2)
	 T(2) = xyzP(2,m3) - xyzP(2,m2)
	 T(3) = xyzP(3,m3) - xyzP(3,m2)
	 dist = dsqrt( T(1)**2 + T(2)**2 + T(3)**2 )
	 nhatP(1,m) = T(1)/dist
	 nhatP(2,m) = T(2)/dist
	 nhatP(3,m) = T(3)/dist
10	 continue
	enddo

	do m = 1,Nxi
	 m3 = IndexXi(m,3)
	 m2 = IndexXi(m,2)
	 T(1) = xyzP(1,m3) - xyzP(1,m2)
	 T(2) = xyzP(2,m3) - xyzP(2,m2)
	 T(3) = xyzP(3,m3) - xyzP(3,m2)
	 dist = dsqrt( T(1)**2 + T(2)**2 + T(3)**2 )
	 nhatP(1,2*Nres+m) = T(1)/dist
	 nhatP(2,2*Nres+m) = T(2)/dist
	 nhatP(3,2*Nres+m) = T(3)/dist
	enddo


C================================================================
C ASSEMBLE THE dRdTor ARRAY FOR ALL DEGREES OF FREEDOM/ ALL ATOMS:
C================================================================

c Initialize:
	do itor = 1,NDoF
        do iatom = 1,Nprtn
         dRdTor(iatom,itor,1) = 0.d00
         dRdTor(iatom,itor,2) = 0.d00
         dRdTor(iatom,itor,3) = 0.d00
	enddo
	enddo


c DoF: MAINCHAIN TORSIONS OF PROTEIN:
c -----------------------------------
        DO itor = 1,2*Nres  

c ignore mainchain proline phi torsions ass well as phi(1):
	 if(Pyn(itor).eq.0)goto20

c update protein component with MC torsion # itor:
         call updateMC(xyzP,xyzdP,IndexMC,nhatP,itor,Nres,Nxi,
     1	 Nprtn,delta)

c compute overall translation/rotation about CoM due to update:
         call u3best(mass,xyzdP,xyzP,Nprtn,1,rms,U,T,ier)  

c eliminate overall translation/rotation about CoM:
	call MoveCoords(xyzdP,Nprtn,0,0,U,T)

c Fill in itor column in array dRdTor:
        do i = 1,Nprtn
          dRdTor(i,itor,1) = ( xyzdP(1,i) - xyzP(1,i) )/delta
          dRdTor(i,itor,2) = ( xyzdP(2,i) - xyzP(2,i) )/delta
          dRdTor(i,itor,3) = ( xyzdP(3,i) - xyzP(3,i) )/delta
        enddo   

20	continue
        ENDDO   ! end itor loop over MainChain DoF of protein



c DoF: SIDECHAIN TORSIONS OF PROTEIN:
c -----------------------------------
        DO itor = 1,Nxi  

c update protein component with SC torsion # itor:
         call updateSC(xyzP,xyzdP,IndexXi,EndXi,nhatP,itor,Nres,
     1    Nxi,Nprtn,delta)

c compute overall translation/rotation about CoM due to update:
         call u3best(mass,xyzdP,xyzP,Nprtn,1,rms,U,T,ier)  

c eliminate overall translation/rotation about CoM:
	call MoveCoords(xyzdP,Nprtn,0,0,U,T)

c Fill in itor column in array dRdTor:
         do i = 1,Nprtn
           dRdTor(i,2*Nres+itor,1) = 
     1		( xyzdP(1,i) - xyzP(1,i) )/delta
           dRdTor(i,2*Nres+itor,2) = 
     1		( xyzdP(2,i) - xyzP(2,i) )/delta
           dRdTor(i,2*Nres+itor,3) = 
     1		( xyzdP(3,i) - xyzP(3,i) )/delta
         enddo  

        ENDDO   ! end itor loop over SC of protein DoF


	return
	end



	subroutine MoveCoords(xyz,Nprtn,Nadp,NM,U,T)

	implicit none
	integer Nprtn,Nadp,NM
	real*8 xyz(3,Nprtn+Nadp+NM),T(3),U(3,3)

	integer iatom
	real*8 x,y,z


        do iatom = 1,Nprtn
         x = xyz(1,iatom)
         y = xyz(2,iatom)
         z = xyz(3,iatom)
         xyz(1,iatom) = T(1) + U(1,1)*x + U(1,2)*y + U(1,3)*z
         xyz(2,iatom) = T(2) + U(2,1)*x + U(2,2)*y + U(2,3)*z
         xyz(3,iatom) = T(3) + U(3,1)*x + U(3,2)*y + U(3,3)*z
        enddo

        do iatom = 1,Nadp
         x = xyz(1,Nprtn+iatom)
         y = xyz(2,Nprtn+iatom)
         z = xyz(3,Nprtn+iatom)
         xyz(1,Nprtn+iatom) = T(1) + U(1,1)*x + U(1,2)*y + U(1,3)*z
         xyz(2,Nprtn+iatom) = T(2) + U(2,1)*x + U(2,2)*y + U(2,3)*z
         xyz(3,Nprtn+iatom) = T(3) + U(3,1)*x + U(3,2)*y + U(3,3)*z
        enddo

        do iatom = 1,NM
         x = xyz(1,Nprtn+Nadp+iatom)
         y = xyz(2,Nprtn+Nadp+iatom)
         z = xyz(3,Nprtn+Nadp+iatom)
         xyz(1,Nprtn+Nadp+iatom)=T(1) + U(1,1)*x + U(1,2)*y + U(1,3)*z
         xyz(2,Nprtn+Nadp+iatom)=T(2) + U(2,1)*x + U(2,2)*y + U(2,3)*z
         xyz(3,Nprtn+Nadp+iatom)=T(3) + U(3,1)*x + U(3,2)*y + U(3,3)*z
        enddo

	return
	end



