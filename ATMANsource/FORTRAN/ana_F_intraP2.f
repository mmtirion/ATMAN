        subroutine ana_F_intraP2(xyz,IndexMC,IndexXi,EndXi,
     1vdWcut,NBI,Nres,Nprtn,Nxi,F,NF,type,Pyn,Rhc,P0,NBIk,IntF,C)

c upgraded 9/16 by mmt, to replace IntActA/B and labelA/B with NBI(N,N)
c and simplify do-loop constructs.

c Computes the analytic 2nd derivative, d2E/dq_idq_j, of a Hookian 
c potential:
c
c	  C 
c   E  =  -  Sum_lk { Rhat^0_lk dot (R_lk - R^0_lk) }**2
c         2
c
c The sum runs over all interacting atom pairs (l,k): these pairs are
c of 2 sorts: intra-protein pairs and inter-protein:ADP pairs. 
c Hence the program is divided into two parts: the first computes
c the contribution to the 2nd derivative matrix of intra-protein
c interactions, and the 2nd part computes the contribution of the 
c intra-protein:ADP interactions to the 2nd derivative matrix.
c The degrees of freedom, q_i (i=1,...,Ntor) include all torsion
c degrees of freedom within the protein component: Nphi + Npsi + Nxi.
c This software requires N (not N**2) execution cycles and is
c therefore faster than the numeric version.
c THe software, developed and written 6/95 by MMT, was tested against
c the numeric version (num_F_prtn), on July 17, 1995, using actin:ADP.



	implicit none
	integer Nprtn,Nres,Nxi,NF

c input:
c   main chain and side chain torsion indices:
	integer IndexMC(2*Nres,4),Pyn(2*Nres+Nxi)
	integer IndexXi(Nxi,4),EndXi(Nxi)
	character*1 type(Nprtn),type1,type2
	real*8 Rhc(Nprtn),vdWcut
c   intra NCB indices:
	integer NBI(Nprtn,Nprtn)
c   vector of DoF angles and coordinates of protein:
	real*8 P0(2*Nres+Nxi),xyz(3,Nprtn)

c computation:
	real*8 NBIk(Nprtn,Nprtn)
	real*8 rij0(3),rim(3),rjn(3),cross(3)
	real*8 energy,factor,dot,term,dist,C
	real*8 psi,fpsi,phi,fphi,xi,fxi
	real*8 RotStiff(2*Nres+Nxi)
	real*8 bmXi(3,Nxi),bmMC(3,2*Nres),sum(2*Nres)
	integer l,m,n,m2,m3,n3,i,j,ider
	integer i1,i2,sign
c   for computing intra-potential terms:
	real*8 dR_ldq_m(3),dr_ldq_n(3),Rlm(3),Rln(3)

c output:
	integer IntF(NF)
	real*8 F(NF)


c Obtain diagonal entries in Hessian due to dihedral stiffnesses:
        psi = P0(2)
        call RotStiffNter(psi,fpsi,energy)
        RotStiff(1) = 0.0d0
        RotStiff(2) = fpsi
        do i = 2, Nres
         phi = P0(2*i-1)
         psi = P0(2*i)
         if(Pyn(i).eq.0)then
          fphi = 0.d0
          call RotStiffPro(phi,psi,fpsi,energy)
         else
          call RotStiffMC(phi,psi,fphi,fpsi,energy)
         endif
         RotStiff(2*i-1) = fphi
         RotStiff(2*i) = fpsi
        enddo
c and then for the SC Xi angles:
        do i = 1, Nxi
         xi = P0(2*Nres+i)
         i1 = indexXi(i,2)
         i2 = indexXi(i,3)
         type1 = type(i1)
         type2 = type(i2)
         call RotStiffSC(xi,type1,type2,fxi,energy)
         if(fxi.le.1.0d0)fxi = 1.0d0  ! to stabilize protruding Args
         RotStiff(2*Nres+i)=fxi
        enddo
         write(*,*)' Side chain stiffness inclusion:'
         write(*,*)'(fxi.le.1.0d0)fxi = 1.0d0'


c set up bm vector:
	do i=1,2*Nres
	 if(Pyn(i).eq.0)goto44
	 m3 = IndexMC(i,3)
	 m2 = IndexMC(i,2)
	 dist = dsqrt( (xyz(1,m3)-xyz(1,m2))**2 +
     1                 (xyz(2,m3)-xyz(2,m2))**2 +
     2                 (xyz(3,m3)-xyz(3,m2))**2 )
	 bmMC(1,i) = ( xyz(1,m3)-xyz(1,m2) )/dist
	 bmMC(2,i) = ( xyz(2,m3)-xyz(2,m2) )/dist
	 bmMC(3,i) = ( xyz(3,m3)-xyz(3,m2) )/dist
44	 continue
	enddo

	do i=1,Nxi
	 m3 = IndexXi(i,3)
	 m2 = IndexXi(i,2)
	 dist = dsqrt( (xyz(1,m3)-xyz(1,m2))**2 +
     1                 (xyz(2,m3)-xyz(2,m2))**2 +
     2                 (xyz(3,m3)-xyz(3,m2))**2 )
	 bmXi(1,i) = ( xyz(1,m3)-xyz(1,m2) )/dist
	 bmXi(2,i) = ( xyz(2,m3)-xyz(2,m2) )/dist
	 bmXi(3,i) = ( xyz(3,m3)-xyz(3,m2) )/dist
	enddo

	do i = 1, NF
	 IntF(i) = 0
	enddo


c===================================================================
c Evaluate intra-protein potential contribution:
c===================================================================

c Main Loop:
	ider = 0
	DO m = 1,2*Nres ! OUTER loop over MAINCHAIN Degrees of Freedom (DOF)
	 m3 = IndexMC(m,3)

	DO n = m,2*Nres
	 n3 = IndexMC(n,3)

	 ider = ider + 1
	 F(ider) = 0.d0
	 if(m.eq.n)F(ider) = RotStiff(m)

	 if(Pyn(m).eq.0.OR.Pyn(n).eq.0)goto96


	sign = 1.d0
 	do i = 1,m3-1  ! locate all contributing i,j atom pairs
	do j = n3+1, Nprtn
	if(NBI(i,j).eq.1)then
	 IntF(ider) = IntF(ider) + 1

	 dist = dsqrt( (xyz(1,i)-xyz(1,j))**2 + 
     1                 (xyz(2,i)-xyz(2,j))**2 +
     2	               (xyz(3,i)-xyz(3,j))**2 )
 
	 rij0(1) = (xyz(1,i)-xyz(1,j))/dist
	 rij0(2) = (xyz(2,i)-xyz(2,j))/dist
	 rij0(3) = (xyz(3,i)-xyz(3,j))/dist

	 rim(1) = xyz(1,i)-xyz(1,m3) 
	 rim(2) = xyz(2,i)-xyz(2,m3) 
	 rim(3) = xyz(3,i)-xyz(3,m3) 

	 cross(1) = bmMC(2,m)*rim(3)-bmMC(3,m)*rim(2)
	 cross(2) = bmMC(3,m)*rim(1)-bmMC(1,m)*rim(3)
 	 cross(3) = bmMC(1,m)*rim(2)-bmMC(2,m)*rim(1)

	dot = rij0(1)*cross(1)+rij0(2)*cross(2)+rij0(3)*cross(3)

	 rjn(1) = xyz(1,j)-xyz(1,n3)
	 rjn(2) = xyz(2,j)-xyz(2,n3)
	 rjn(3) = xyz(3,j)-xyz(3,n3)

	 cross(1) = bmMC(2,n)*rjn(3)-bmMC(3,n)*rjn(2)
	 cross(2) = bmMC(3,n)*rjn(1)-bmMC(1,n)*rjn(3)
 	 cross(3) = bmMC(1,n)*rjn(2)-bmMC(2,n)*rjn(1)
	
	 factor = NBIk(i,j)

	 F(ider) = F(ider) + factor * dot * sign *
     1(rij0(1)*cross(1) + rij0(2)*cross(2) + rij0(3)*cross(3))


	endif
	enddo ! j loop
	enddo ! i loop

96	continue
	ENDDO ! n loop


	DO n = 1,Nxi ! CONTINUE INNER LOOP OVER ALL SIDECHAIN DOF
	 ider = ider+1
	 F(ider) = 0.d0
  	 term = 0.d0

	 if(Pyn(m).eq.0) goto200  ! n3 should never be zero...
	 n3 = indexXi(n,3)

	IF(m3.le.indexXi(n,2))THEN

	 sign = 1.d0
	 do i = 1, m3-1
	 do j = n3+1,EndXi(n)
	  if(NBI(i,j).eq.1)then
	 IntF(ider) = IntF(ider) + 1

	   dist = dsqrt( (xyz(1,i)-xyz(1,j))**2 + 
     1                   (xyz(2,i)-xyz(2,j))**2 +
     2	                 (xyz(3,i)-xyz(3,j))**2 )
 
	   rij0(1) = (xyz(1,i)-xyz(1,j))/dist  
	   rij0(2) = (xyz(2,i)-xyz(2,j))/dist  
	   rij0(3) = (xyz(3,i)-xyz(3,j))/dist  

	   rim(1) = xyz(1,i)-xyz(1,m3) 
	   rim(2) = xyz(2,i)-xyz(2,m3) 
	   rim(3) = xyz(3,i)-xyz(3,m3) 
 
	   cross(1) = bmMC(2,m)*rim(3)-bmMC(3,m)*rim(2)
	   cross(2) = bmMC(3,m)*rim(1)-bmMC(1,m)*rim(3)
 	   cross(3) = bmMC(1,m)*rim(2)-bmMC(2,m)*rim(1)

	 dot = rij0(1)*cross(1)+rij0(2)*cross(2)+rij0(3)*cross(3)

	   rjn(1) = xyz(1,j)-xyz(1,n3)
	   rjn(2) = xyz(2,j)-xyz(2,n3)
	   rjn(3) = xyz(3,j)-xyz(3,n3)

	   cross(1) = bmXi(2,n)*rjn(3)-bmXi(3,n)*rjn(2)
	   cross(2) = bmXi(3,n)*rjn(1)-bmXi(1,n)*rjn(3)
 	   cross(3) = bmXi(1,n)*rjn(2)-bmXi(2,n)*rjn(1)

	   factor = NBIk(i,j)

	   term = term + factor * dot * sign *
     1(rij0(1)*cross(1) + rij0(2)*cross(2) + rij0(3)*cross(3))


	 endif
	enddo ! j loop
	enddo ! i loop

	ELSE

	 sign = -1.d0
	 do i = n3+1, EndXi(n)
	 do j = m3+1, Nprtn
	  if(NBI(i,j).eq.1)then
	 IntF(ider) = IntF(ider) + 1

	   dist = dsqrt( (xyz(1,i)-xyz(1,j))**2 + 
     1                   (xyz(2,i)-xyz(2,j))**2 +
     2	                 (xyz(3,i)-xyz(3,j))**2 )
 
	   rij0(1) = (xyz(1,i)-xyz(1,j))/dist  
	   rij0(2) = (xyz(2,i)-xyz(2,j))/dist  
	   rij0(3) = (xyz(3,i)-xyz(3,j))/dist  
 
 	   rim(1) = xyz(1,i)-xyz(1,m3) 
	   rim(2) = xyz(2,i)-xyz(2,m3) 
	   rim(3) = xyz(3,i)-xyz(3,m3) 

 	   cross(1) = bmMC(2,m)*rim(3)-bmMC(3,m)*rim(2)
	   cross(2) = bmMC(3,m)*rim(1)-bmMC(1,m)*rim(3)
 	   cross(3) = bmMC(1,m)*rim(2)-bmMC(2,m)*rim(1)

	 dot = rij0(1)*cross(1)+rij0(2)*cross(2)+rij0(3)*cross(3)

	   rjn(1) = xyz(1,j)-xyz(1,n3)
	   rjn(2) = xyz(2,j)-xyz(2,n3)
	   rjn(3) = xyz(3,j)-xyz(3,n3)
 
	   cross(1) = bmXi(2,n)*rjn(3)-bmXi(3,n)*rjn(2)
	   cross(2) = bmXi(3,n)*rjn(1)-bmXi(1,n)*rjn(3)
 	   cross(3) = bmXi(1,n)*rjn(2)-bmXi(2,n)*rjn(1)

	   factor = NBIk(i,j)

	   term = term + factor * dot * sign *
     1(rij0(1)*cross(1) + rij0(2)*cross(2) + rij0(3)*cross(3))


	 endif
	enddo  ! j loop
	enddo  ! i loop

	ENDIF ! if m3 < indexXi(n,2)

200	continue
	F(ider) = F(ider) + term 

	ENDDO  ! n loop
	ENDDO  ! m loop

c this ends MC-MC block (no rigid body dofs)



	DO m = 1,Nxi  ! OUTER loop over all SC DOF
	 m3 = indexXi(m,3)

	DO n = m,Nxi ! INNER loop over SC DOF
	 n3 = IndexXi(n,3)

	 ider = ider+1
	 F(ider) = 0.d0
	 if(m.eq.n)F(ider)=RotStiff(2*Nres+m)
	 term = 0.0d00


	IF(EndXi(m).ne.EndXi(n))THEN

	sign = -1.d0
	do i = m3+1,EndXi(m)
	do j = n3+1,EndXi(n)
	 if(NBI(i,j).eq.1)then
	 IntF(ider) = IntF(ider) + 1

	  dist = dsqrt( (xyz(1,i)-xyz(1,j))**2 + 
     1                  (xyz(2,i)-xyz(2,j))**2 +
     2	 	        (xyz(3,i)-xyz(3,j))**2 )

	  rij0(1) = (xyz(1,i)-xyz(1,j))/dist
	  rij0(2) = (xyz(2,i)-xyz(2,j))/dist
	  rij0(3) = (xyz(3,i)-xyz(3,j))/dist

	  rim(1) = xyz(1,i)-xyz(1,m3) 
	  rim(2) = xyz(2,i)-xyz(2,m3) 
	  rim(3) = xyz(3,i)-xyz(3,m3) 

	  cross(1) = bmXi(2,m)*rim(3) - bmXi(3,m)*rim(2)
	  cross(2) = bmXi(3,m)*rim(1) - bmXi(1,m)*rim(3)
 	  cross(3) = bmXi(1,m)*rim(2) - bmXi(2,m)*rim(1)

	 dot = rij0(1)*cross(1)+rij0(2)*cross(2)+rij0(3)*cross(3)

	  rjn(1) = xyz(1,j)-xyz(1,n3)
	  rjn(2) = xyz(2,j)-xyz(2,n3)
	  rjn(3) = xyz(3,j)-xyz(3,n3)

	  cross(1) = bmXi(2,n)*rjn(3) - bmXi(3,n)*rjn(2)
	  cross(2) = bmXi(3,n)*rjn(1) - bmXi(1,n)*rjn(3)
 	  cross(3) = bmXi(1,n)*rjn(2) - bmXi(2,n)*rjn(1)

	  factor = NBIk(i,j)

	  term = term + factor * dot * sign *
     1(rij0(1)*cross(1) + rij0(2)*cross(2) + rij0(3)*cross(3))

	 endif
	enddo ! j loop 
	enddo ! i loop

	ELSE

	sign = 1.d0
	do i = 1, indexXi(m,2)
	do j = n3+1, EndXi(n)
	 if(NBI(i,j).eq.1)then
	 IntF(ider) = IntF(ider) + 1

	  dist = dsqrt( (xyz(1,i)-xyz(1,j))**2 + 
     1                  (xyz(2,i)-xyz(2,j))**2 +
     2	 	        (xyz(3,i)-xyz(3,j))**2 )

	  rij0(1) = (xyz(1,i)-xyz(1,j))/dist
	  rij0(2) = (xyz(2,i)-xyz(2,j))/dist
	  rij0(3) = (xyz(3,i)-xyz(3,j))/dist

	  rim(1) = xyz(1,i)-xyz(1,m3) 
	  rim(2) = xyz(2,i)-xyz(2,m3) 
	  rim(3) = xyz(3,i)-xyz(3,m3) 

	  cross(1) = bmXi(2,m)*rim(3) - bmXi(3,m)*rim(2)
	  cross(2) = bmXi(3,m)*rim(1) - bmXi(1,m)*rim(3)
 	  cross(3) = bmXi(1,m)*rim(2) - bmXi(2,m)*rim(1)

	 dot = rij0(1)*cross(1)+rij0(2)*cross(2)+rij0(3)*cross(3)

	  rjn(1) = xyz(1,j)-xyz(1,n3)
	  rjn(2) = xyz(2,j)-xyz(2,n3)
	  rjn(3) = xyz(3,j)-xyz(3,n3)

	  cross(1) = bmXi(2,n)*rjn(3)-bmXi(3,n)*rjn(2)
	  cross(2) = bmXi(3,n)*rjn(1)-bmXi(1,n)*rjn(3)
 	  cross(3) = bmXi(1,n)*rjn(2)-bmXi(2,n)*rjn(1)

	  factor = NBIk(i,j)

	  term = term + factor * dot * sign *
     1(rij0(1)*cross(1) + rij0(2)*cross(2) + rij0(3)*cross(3))


	 endif
	enddo  
	enddo  


	sign = 1.d0
	do i = n3+1,EndXi(n)
	do j = EndXi(n)+1, Nprtn
	 if(NBI(i,j).eq.1)then
	IntF(ider) = IntF(ider) + 1

	  dist = dsqrt( (xyz(1,i)-xyz(1,j))**2 + 
     1                  (xyz(2,i)-xyz(2,j))**2 +
     2	 	        (xyz(3,i)-xyz(3,j))**2 )

	  rij0(1) = (xyz(1,i)-xyz(1,j))/dist
	  rij0(2) = (xyz(2,i)-xyz(2,j))/dist
	  rij0(3) = (xyz(3,i)-xyz(3,j))/dist

	  rim(1) = xyz(1,i)-xyz(1,m3) 
	  rim(2) = xyz(2,i)-xyz(2,m3) 
	  rim(3) = xyz(3,i)-xyz(3,m3) 

	  cross(1) = bmXi(2,m)*rim(3) - bmXi(3,m)*rim(2)
	  cross(2) = bmXi(3,m)*rim(1) - bmXi(1,m)*rim(3)
 	  cross(3) = bmXi(1,m)*rim(2) - bmXi(2,m)*rim(1)

	 dot = rij0(1)*cross(1)+rij0(2)*cross(2)+rij0(3)*cross(3)

	  rjn(1) = xyz(1,j)-xyz(1,n3)
	  rjn(2) = xyz(2,j)-xyz(2,n3)
	  rjn(3) = xyz(3,j)-xyz(3,n3)

	  cross(1) = bmXi(2,n)*rjn(3)-bmXi(3,n)*rjn(2)
	  cross(2) = bmXi(3,n)*rjn(1)-bmXi(1,n)*rjn(3)
 	  cross(3) = bmXi(1,n)*rjn(2)-bmXi(2,n)*rjn(1)

	  factor = NBIk(i,j)

	  term = term + factor * dot * sign *
     1(rij0(1)*cross(1) + rij0(2)*cross(2) + rij0(3)*cross(3))


	 Endif
	enddo
	enddo  

	ENDIF

	 F(ider) = F(ider) +  term  
	ENDDO  ! END INNER LOOP OVER ALL SIDECHAIN DOF
	ENDDO  ! END OUTER LOOP OVER ALL SIDECHAIN DOF



	return
	end

