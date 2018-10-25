	subroutine ana_F_ATMAN(xyz,IndexMC,IndexXi,EndXi,NBIPP,Nres,
     1Nprtn,Nxi,Ninter,Nlig,Lig,mHET,xyzL,NBIPL,F,NF,type,Pyn,P0,NligNA,
     2CPP,CPL,NinterL,NBILL,CLL,IntF)

c upgraded 11/16 to include RB dofs. Block structure:
c
c        _  MC        SC         RB  _
c       |                             |
c   MC  |   I                         |
c       |                             |
c       |                             |
c   SC  |   II        III             |
c       |                             |
c   RBi |   IV         V        VI    |
c       |_                           _|
c
c Where MC are mainchain degrees of freedom,
c SC are sidechain degrees of freedom,
c and RB are rigid body degrees of freedom, running over all i=1,Nlig ligands, i
c adding at most 6*Nlig dofs, and fewer if some ligands consist of 1 atom 
c requiring only 3 rigid body degrees of freedom.



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
	integer Nprtn,Nres,Nxi,Ninter,NinterL,NligNA,Nlig,NF

c input:
c   main chain and side chain torsion indices:
	integer IndexMC(2*Nres,4),Pyn(2*Nres+Nxi)
	integer IndexXi(Nxi,4),EndXi(Nxi)
	character*1 type(Nprtn),type1,type2
c   intra NCB indices:
	integer NBIPP(Nprtn,Nprtn)
c   vector of DoF angles and coordinates of protein:
	real*8 P0(2*Nres+Nxi),xyz(3,Nprtn)
c Ligand stuff:
	integer Lig(1000,2),NBIPL(10000,2)
	integer NBILL(2500,2500)
	real*8 xyzL(3,2500),mHET(2500)
	real*8 Rcm_lig(3,1000),TmassL(1000)
	real*8 xyzLcom(3),unit(3),unitm(3),unitn(3)

c computation:
        real*8 CPP(Nprtn,Nprtn),CPL(10000),CLL(2500,2500)
	real*8 rij0(3),rim(3),rjn(3),cross(3),sign,vector(3)
	real*8 energy,factor,dot,value,term,dist
	real*8 psi,fpsi,phi,fphi,xi,fxi
	real*8 RotStiff(2*Nres+Nxi)
	real*8 bmXi(3,Nxi),bmMC(3,2*Nres),sum(2*Nres)
	real*8 avephi,avepsi,avexi
	integer l,ls,m,n,m2,m3,n3,i,j,k,k1,k2,icount,ider
	integer is,ks,k_m,k_n,nlo,i1,i2
	integer iLig,OnOff(2500),OnOff2(2500)
c   for computing intra-potential terms:
	real*8 dR_ldq_m(3),dr_ldq_n(3),Rlm(3),Rln(3)

c output:
	integer IntF(NF)
	real*8 F(NF)


	avexi =0.d0
	avepsi = 0.d0
	avephi = 0.d0
c Obtain diagonal entries in Hessian due to dihedral stiffnesses:
        psi = P0(2)
        call RotStiffNter(psi,fpsi,energy)
        RotStiff(1) = 0.0d0
        RotStiff(2) = fpsi
	avephi = avephi + 0.0d0
	avepsi = avepsi + fpsi
        do i = 2, Nres
         phi = P0(2*i-1)
         psi = P0(2*i)
         if(Pyn(i).eq.0)then
          fphi = 0.d0
          call RotStiffPro(phi,psi,fpsi,energy)
c	  if(fpsi.le.1.d0)fpsi=1.d0  !fine-tune MC torsion stiffness
         else
          call RotStiffMC(phi,psi,fphi,fpsi,energy)
         endif
	if(fphi.le.1.d0)fphi=1.d0 !fine-tune MC dihedral stiffnesses
	if(fpsi.le.1.d0)fpsi=1.d0 !fine-tune MC dihedral stiffnesses
         RotStiff(2*i-1) = fphi
         RotStiff(2*i) = fpsi
	avephi = avephi + fphi
	avepsi = avepsi + fpsi
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
c	 write(53,*)i,sngl(fxi)
	 avexi = avexi + fxi
        enddo
         write(*,*)' Side chain stiffness inclusion:'
         write(*,*)'(fxi.le.1.0d0)fxi = 1.0d0'
	write(*,*)' '
	write(*,*)'Ave Phi: ', sngl(avephi)/real(2*Nres-7-1)
	write(*,*)'Ave Psi: ', sngl(avepsi)/real(2*Nres)
	write(*,*)'Ave Xi: ', sngl(avexi)/real(Nxi)
	write(*,*)' '

c TEMPORARY OVERWRITE: TIRION TYPE POTENTIAL 3/20/18
c 	write(*,*)'WARNING: OVERWRITING TIRION TYPE POTENTIAL W/ fxi=0.0'
c 	do i = 1, 2*Nres+Nxi
c 	 RotStiff(i)=0.5d0
c  	enddo
c ENDTEMPORARY

c TEMPORARY CHECK ON EFFECT OF 1CEM's TAIL MOTILITY:
c	do i = 2*Nres-14,2*Nres  !Immobilize last 7 main chain residues' angles
c	do i = 1,14 ! N-terminal seven residues
c	do i = 15,28 !MC residues 40-46 of PDB restraint
c	do i = 83,96  ! MC residues 74-80 of PDB file
c	do i = 215,228 ! MC residues 140-146
c	do i = 159,172 ! MC residues 112-118
c	do i = 237,250 ! MC 151-157
c	do i = 469,482 ! MC 267-273
c	do i = 33,46 !MC 49-55
c	do i = 319,332  ! MC 192-198
c	do i = 537,550  !MC 301-307
c	do i = 725,726  !MC 363
c	do i = 1,2*Nres
c	 RotStiff(i) = 0.1d0*RotStiff(i) !Downsize L79 contribution to dihedral energy
c	 if(RotStiff(i).le.0.1d0)RotStiff(i)=0.1d0
c	 RotStiff(i)=0.1d0
c	enddo
c	do i = 2*Nres+1, 2*Nres+Nxi
c	 RotStiff(i) = 0.1d0*RotStiff(i) !Downsize L79 contribution to dihedral energy
c	 if(RotStiff(i).le.1.0d0)RotStiff(i)=1.0d0
c 	 RotStiff(i) = 0.1d0
c	enddo
c 	do i = 1, 2*Nres + Nxi
c 	 RotStiff(i) = 0.1d0
c 	enddo
c	do i = 1, Nprtn
c	do j = 1, Nprtn
c	 CPP(i,j) = 1.d0
c	enddo
c	enddo
c END TEMPORARY CHECK
	close(32)

c  Get Center of Mass of each Ligand (for RB rotations)
        i = 0
        DO k = 1, Nlig
        IF(Lig(k,1).gt.1)THEN
         TmassL(k) = 0.d0
         Rcm_lig(1,k) = 0.d0
         Rcm_lig(2,k) = 0.d0
         Rcm_lig(3,k) = 0.d0
        do ls = 1, Lig(k,1)
          i = i + 1
          TmassL(k) = TmassL(k) + mHET(i)
          Rcm_lig(1,k) = Rcm_lig(1,k) + mHET(i)*xyzL(1,i)
          Rcm_lig(2,k) = Rcm_lig(2,k) + mHET(i)*xyzL(2,i)
          Rcm_lig(3,k) = Rcm_lig(3,k) + mHET(i)*xyzL(3,i)
        enddo
         Rcm_lig(1,k) = Rcm_lig(1,k)/TmassL(k)
         Rcm_lig(2,k) = Rcm_lig(2,k)/TmassL(k)
         Rcm_lig(3,k) = Rcm_lig(3,k)/TmassL(k)
        ELSE
         i = i + 1
         Rcm_lig(1,k) = xyzL(1,i)
         Rcm_lig(2,k) = xyzL(2,i)
         Rcm_lig(3,k) = xyzL(3,i)
        ENDIF
        ENDDO  


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

c Main Loop 
	ider = 0
	sign = 1.d0
	DO m = 1,2*Nres ! BLOCK I: MC-MC
c TEMPORARY
	write(*,*)m
c END TEMPORARY
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
	if(NBIPP(i,j).eq.1)then
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
	
	 factor = CPP(i,j)

	 F(ider) = F(ider) + factor * dot * sign *
     1(rij0(1)*cross(1) + rij0(2)*cross(2) + rij0(3)*cross(3))


	endif
	enddo ! j loop
	enddo ! i loop

c include ligand interactions:
c
c	sign = 1.d0 !unknown
        do iLig = 1, Ninter
         j = NBIPL(iLig,2) ! belongs to protein subset
         if(j.gt.n3)then
          i = NBIPL(iLig,1) ! belongs to ligand subset
	  factor = CPL(iLig)
         
	  dist = dsqrt((xyzL(1,i)-xyz(1,j))**2 + 
     1                 (xyzL(2,i)-xyz(2,j))**2 +
     2	               (xyzL(3,i)-xyz(3,j))**2 )
 
	  rij0(1) = (xyzL(1,i)-xyz(1,j))/dist
	  rij0(2) = (xyzL(2,i)-xyz(2,j))/dist
	  rij0(3) = (xyzL(3,i)-xyz(3,j))/dist

	  rim(1) = xyzL(1,i)-xyz(1,m3) 
	  rim(2) = xyzL(2,i)-xyz(2,m3) 
	  rim(3) = xyzL(3,i)-xyz(3,m3) 

	  cross(1) = bmMC(2,m)*rim(3)-bmMC(3,m)*rim(2)
	  cross(2) = bmMC(3,m)*rim(1)-bmMC(1,m)*rim(3)
 	  cross(3) = bmMC(1,m)*rim(2)-bmMC(2,m)*rim(1)

	term = rij0(1)*cross(1)+rij0(2)*cross(2)+rij0(3)*cross(3)

	  rjn(1) = xyz(1,j)-xyz(1,n3)
	  rjn(2) = xyz(2,j)-xyz(2,n3)
	  rjn(3) = xyz(3,j)-xyz(3,n3)
 
	  cross(1) = bmMC(2,n)*rjn(3)-bmMC(3,n)*rjn(2)
	  cross(2) = bmMC(3,n)*rjn(1)-bmMC(1,n)*rjn(3)
 	  cross(3) = bmMC(1,n)*rjn(2)-bmMC(2,n)*rjn(1)
	
	value = rij0(1)*cross(1)+rij0(2)*cross(2)+rij0(3)*cross(3)

	  F(ider) = F(ider) + factor * term * sign * value

         endif
        enddo ! iLig loop
c
c end ligand contribution to block I

96	continue
	ENDDO ! n loop


	DO n = 1,Nxi  ! BLOCK II: MC-SC
	 ider = ider+1
	 F(ider) = 0.d0
  	 term = 0.d0

	 if(Pyn(m).eq.0) goto200  ! n3 should never be zero...
	 n3 = indexXi(n,3)

	IF(m3.lt.EndXi(n))THEN

	 sign = 1.d0
	 do i = 1, m3-1
	 do j = n3+1,EndXi(n)
	  if(NBIPP(i,j).eq.1)then
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

	   factor = CPP(i,j)

	   term = term + factor * dot * sign *
     1(rij0(1)*cross(1) + rij0(2)*cross(2) + rij0(3)*cross(3))


	 endif
	enddo ! j loop
	enddo ! i loop

c include ligand interactions:
c
c	sign = 1.d0 !unknown 
        do iLig = 1, Ninter
         j = NBIPL(iLig,2) ! belongs to protein subset
         if(j.gt.n3.and.j.le.EndXi(n))then
          i = NBIPL(iLig,1) ! belongs to ligand subset
	   factor = CPL(iLig)
          
	   dist = dsqrt( (xyzL(1,i)-xyz(1,j))**2 + 
     1                   (xyzL(2,i)-xyz(2,j))**2 +
     2	                 (xyzL(3,i)-xyz(3,j))**2 )
 
	   rij0(1) = (xyzL(1,i)-xyz(1,j))/dist  
	   rij0(2) = (xyzL(2,i)-xyz(2,j))/dist  
	   rij0(3) = (xyzL(3,i)-xyz(3,j))/dist  

	   rim(1) = xyzL(1,i)-xyz(1,m3) 
	   rim(2) = xyzL(2,i)-xyz(2,m3) 
	   rim(3) = xyzL(3,i)-xyz(3,m3) 
 
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

	 value = rij0(1)*cross(1)+rij0(2)*cross(2)+rij0(3)*cross(3)

	   term = term + factor * dot * sign * value

         endif
        enddo ! iLig loop
c
c ligand interaction

	ELSE

	 sign = -1.d0
	 do i = n3+1, EndXi(n)
	 do j = m3+1, Nprtn
	  if(NBIPP(i,j).eq.1)then
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

	   factor = CPP(i,j)

	   term = term + factor * dot * sign *
     1(rij0(1)*cross(1) + rij0(2)*cross(2) + rij0(3)*cross(3))


	 endif
	enddo  ! j loop
	enddo  ! i loop

c no ligand contributions here


	ENDIF ! if m3 < indexXi(n,2)


	F(ider) = F(ider) + term 
200	continue

	ENDDO  ! n loop

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                      added 11/16                             c
c                    Blocks in FIRST column                   c
c            for rigid-body D.o.F.s of ligands                c

        do k = 1, Nlig  ! BLOCK IV: MC-RB, first column for EACH ligand...

c select subset of Ligand array to activate by first initializing OnOff to zero
        do ks = 1, NligNA
         OnOff(ks) = 0
        enddo
c and then activating each ligand subset by setting appropriate OnOff entries to 1
        icount = 0
        do ks = 1, Nlig
        do ls = 1, Lig(ks,1)
         icount = icount + 1
         if(ks.eq.k)OnOff(icount) = 1
        enddo
        enddo

c now apply each rigid-body degree of freedom update to ligand-subset:
        ks = 3  ! for 1-atom RBs, only 3 RB dofs
        if(Lig(k,1).gt.1)ks=6
        do n = 1,ks
         ider = ider + 1
         F(ider) = 0.0d0
	 if(Pyn(m).eq.0) goto201  ! n3 should never be zero...
	 unit(1) = 0.d0
	 unit(2) = 0.d0
	 unit(3) = 0.d0
	 if(n.gt.3)unit(n-3)=1.d0 !this selects RB rotation axis
	 if(n-3.eq.2)unit(2)=-1.d0 !who knows, based on fit to numeric

        do is = 1, Ninter
         i = NBIPL(is,1)  ! ligand atom
         if(OnOff(i).eq.0)goto445
         j = NBIPL(is,2)    ! protein atom

	 if(j.gt.m3)then

	  dist = dsqrt( (xyzL(1,i)-xyz(1,j))**2 + 
     1                  (xyzL(2,i)-xyz(2,j))**2 +
     2	                (xyzL(3,i)-xyz(3,j))**2 )
 
	  rij0(1) = (xyzL(1,i)-xyz(1,j))/dist  
	  rij0(2) = (xyzL(2,i)-xyz(2,j))/dist  
	  rij0(3) = (xyzL(3,i)-xyz(3,j))/dist  

	  rim(1) = xyzL(1,i)-xyz(1,m3) 
	  rim(2) = xyzL(2,i)-xyz(2,m3) 
	  rim(3) = xyzL(3,i)-xyz(3,m3) 

	  vector(1) = bmMC(1,m)
	  vector(2) = bmMC(2,m)
	  vector(3) = bmMC(3,m)
 
	  call AcrossB(vector,rim,cross)

	  call AdotB(rij0,cross,term)

	  factor = CPL(is)

	  IF(n.le.3)THEN

	   value = -1.d0 * rij0(n) !atom b is ahead of dihedral angle...

	  ELSE

	   xyzLcom(1) = xyzL(1,i) - Rcm_lig(1,k)
           xyzLcom(2) = xyzL(2,i) - Rcm_lig(2,k)
           xyzLcom(3) = xyzL(3,i) - Rcm_lig(3,k)

	   call AcrossB(unit,xyzLcom,cross)

	   call AdotB(rij0,cross,value)

	  ENDIF

	  F(ider) = F(ider) + factor*term*sign*value

	 endif
445      continue
        enddo ! loop is = 1, Ninter
201	continue
        enddo  !loop n = 1,...6
        enddo  !loop k ligands

c                                                             c
c                                                             c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	ENDDO  ! m loop



	DO m = 1,Nxi  ! BLOCK III:  SC-SC ; second column
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
	 if(NBIPP(i,j).eq.1)then
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

	  factor = CPP(i,j)

	  term = term + factor * dot * sign *
     1(rij0(1)*cross(1) + rij0(2)*cross(2) + rij0(3)*cross(3))

	 endif
	enddo ! j loop 
	enddo ! i loop

c no ligand contributions here 

	ELSE

	sign = 1.d0
	do i = 1, indexXi(m,2)
	do j = n3+1, EndXi(n)
	 if(NBIPP(i,j).eq.1)then
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

	  factor = CPP(i,j)

	  term = term + factor * dot * sign *
     1(rij0(1)*cross(1) + rij0(2)*cross(2) + rij0(3)*cross(3))


	 endif
	enddo  
	enddo  


	sign = 1.d0
	do i = n3+1,EndXi(n)
	do j = EndXi(n)+1, Nprtn
	 if(NBIPP(i,j).eq.1)then
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

	  factor = CPP(i,j)

	  term = term + factor * dot * sign *
     1(rij0(1)*cross(1) + rij0(2)*cross(2) + rij0(3)*cross(3))


	 Endif
	enddo
	enddo  

c include ligand interactions:
c
c	sign = 1.d0 !unknown
        do iLig = 1, Ninter
         j = NBIPL(iLig,2) ! belongs to protein subset
         if(j.gt.n3.and.j.le.EndXi(n))then
          i = NBIPL(iLig,1) ! belongs to ligand subset
	  factor = CPL(iLig)
          
	  dist = dsqrt( (xyzL(1,i)-xyz(1,j))**2 + 
     1                  (xyzL(2,i)-xyz(2,j))**2 +
     2	                (xyzL(3,i)-xyz(3,j))**2 )
 
	  rij0(1) = (xyzL(1,i)-xyz(1,j))/dist  
	  rij0(2) = (xyzL(2,i)-xyz(2,j))/dist  
	  rij0(3) = (xyzL(3,i)-xyz(3,j))/dist  

	  rim(1) = xyzL(1,i)-xyz(1,m3) 
	  rim(2) = xyzL(2,i)-xyz(2,m3) 
	  rim(3) = xyzL(3,i)-xyz(3,m3) 
 
	  cross(1) = bmXi(2,m)*rim(3)-bmXi(3,m)*rim(2)
	  cross(2) = bmXi(3,m)*rim(1)-bmXi(1,m)*rim(3)
 	  cross(3) = bmXi(1,m)*rim(2)-bmXi(2,m)*rim(1)

	 dot = rij0(1)*cross(1)+rij0(2)*cross(2)+rij0(3)*cross(3)

	  rjn(1) = xyz(1,j)-xyz(1,n3)
	  rjn(2) = xyz(2,j)-xyz(2,n3)
	  rjn(3) = xyz(3,j)-xyz(3,n3)

	  cross(1) = bmXi(2,n)*rjn(3)-bmXi(3,n)*rjn(2)
	  cross(2) = bmXi(3,n)*rjn(1)-bmXi(1,n)*rjn(3)
 	  cross(3) = bmXi(1,n)*rjn(2)-bmXi(2,n)*rjn(1)

	 value = rij0(1)*cross(1)+rij0(2)*cross(2)+rij0(3)*cross(3)

	  term = term + factor * dot * sign * value

         endif
        enddo !iLig loop
c
c ligand loop


	ENDIF

	 F(ider) = F(ider) +  term  
	ENDDO  ! END INNER LOOP OVER ALL SIDECHAIN DOF
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                      added 11/16                             c
c                   Blocks in SECOND column                   c
c            for rigid-body D.o.F.s of ligands                c

        do k = 1, Nlig  ! BLOCK V: SC-RB; second column for EACH ligand...

c select subset of Ligand array to activate by first initializing OnOff to zero
        do ks = 1, NligNA
         OnOff(ks) = 0
        enddo
c and then activating each ligand subset by setting appropriate OnOff entries to 1
        icount = 0
        do ks = 1,Nlig
        do ls = 1,Lig(ks,1)
         icount = icount + 1
         if(ks.eq.k)OnOff(icount) = 1
        enddo
        enddo

c now apply each rigid-body degree of freedom update to ligand-subset:
        ks = 3
        if(Lig(k,1).gt.1)ks = 6
        do n = 1, ks
         ider = ider + 1
         F(ider) = 0.0d0
	 unit(1) = 0.d0
	 unit(2) = 0.d0
	 unit(3) = 0.d0
	 if(n.gt.3)unit(n-3)=1.d0 !this selects RB rotation axis
	 if(n-3.eq.2)unit(2)=-1.d0 !matching numeric

        do is = 1, Ninter
         i = NBIPL(is,1)  ! labels ligand atom
         if(OnOff(i).eq.0)goto447   !include only 'active' ligands...
         j = NBIPL(is,2)    ! j labels protein atoms

	 if(j.gt.m3.and.j.le.EndXi(m))then
          factor = CPL(is)
	  sign = 1.d0 !your guess as good as mine

	   dist = dsqrt( (xyzL(1,i)-xyz(1,j))**2 + 
     1                   (xyzL(2,i)-xyz(2,j))**2 +
     2	                 (xyzL(3,i)-xyz(3,j))**2 )
 
	   rij0(1) = (xyzL(1,i)-xyz(1,j))/dist  
	   rij0(2) = (xyzL(2,i)-xyz(2,j))/dist  
	   rij0(3) = (xyzL(3,i)-xyz(3,j))/dist  

	   rim(1) = xyzL(1,i)-xyz(1,m3) 
	   rim(2) = xyzL(2,i)-xyz(2,m3) 
	   rim(3) = xyzL(3,i)-xyz(3,m3) 

           vector(1) = bmXi(1,m)
           vector(2) = bmXi(2,m)
           vector(3) = bmXi(3,m)

	   call AcrossB(vector,rim,cross)

	   call AdotB(rij0,cross,term)

	   IF(n.le.3)THEN

	    value = -1.d0 * rij0(n) !atom b is ahead of dihedral angle,could be negative

	   ELSE
	
	    xyzLcom(1) = xyzL(1,i) - Rcm_lig(1,k)
            xyzLcom(2) = xyzL(2,i) - Rcm_lig(2,k)
            xyzLcom(3) = xyzL(3,i) - Rcm_lig(3,k)

	    call AcrossB(unit,xyzLcom,cross)

	    call AdotB(rij0,cross,value)

	   ENDIF

            F(ider) =  F(ider) + sign*factor*term*value 

	 endif
447      continue
        enddo ! loop is = 1, Ninter


        enddo  !loop n = 1,ks
        enddo  !loop k ligands
c                                                             c
c                                                             c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	ENDDO  ! END OUTER LOOP OVER ALL SIDECHAIN DOF

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                      added 11/16                            c
c                Blocks AFTER  MC & SC D.o.F.                 c
c            for rigid-body D.o.F.s of ligands                c

        DO k_m = 1, Nlig  !BLOCK VI: RB-RB; third column

c select subset of Ligand array to activate by first initializing OnOff to zero
        do ks = 1, NligNA
         OnOff(ks) = 0
        enddo
c and then activating each ligand subset by setting appropriate OnOff entries to 1
        icount = 0
        do ks = 1,Nlig
        do ls = 1,Lig(ks,1)
         icount = icount + 1
         if(ks.eq.k_m)OnOff(icount) = 1
        enddo
        enddo

        k1 = 3
        if(Lig(k_m,1).gt.1)k1=6
        Do m = 1,k1
	 unitm(1)=0.d0
	 unitm(2)=0.d0
	 unitm(3)=0.d0
	 if(m.gt.3)unitm(m-3)=-1.d0 !selects rotation axis
	 if(m-3.eq.2)unitm(2)=1.d0 ! due to check analytic vs numeric version

        DO k_n = k_m, Nlig
         nlo=1
         if(k_n.eq.k_m)nlo=m ! diagonal blocks...

c select subset of Ligand array to activate by first initializing OnOff to zero
        do ks = 1, NligNA
         OnOff2(ks) = 0
        enddo
c and then activating each ligand subset by setting appropriate OnOff entries to 1
        icount = 0
        do ks = 1,Nlig
        do ls = 1,Lig(ks,1)
         icount = icount + 1
         if(ks.eq.k_n)OnOff2(icount) = 1
        enddo
        enddo

        k2 = 3
        if(Lig(k_n,1).gt.1)k2=6
        Do n = nlo,k2
         ider = ider + 1
         F(ider) = 0.0d0
	 sign = 1.0d0
	 unitn(1) = 0.d0
	 unitn(2) = 0.d0
	 unitn(3) = 0.d0
	 if(n.gt.3)unitn(n-3)=-1.d0 !selects rotation axis
	 if(n-3.eq.2)unitn(2)=1.d0 

	if(k_m.eq.k_n)then  ! both dof belong to one ligand:

c include appropriate inter protein-ligand pairs...
        do i = 1, Ninter  ! over Ninter Interpairs 
         k = NBIPL(i,1) ! ligand atom
	 if(OnOff(k).eq.0)goto451  !include only ligand_m(n) atoms
c OnOff(k) = OnOff2(k) since k_m = k_n
         j = NBIPL(i,2) ! protein atom
c No restriction: use all inter-ligand protein NBI
         factor = CPL(i)

	 dist = dsqrt( (xyzL(1,k)-xyz(1,j))**2 + 
     1                 (xyzL(2,k)-xyz(2,j))**2 +
     2	               (xyzL(3,k)-xyz(3,j))**2 )
 
	 rij0(1) = (xyzL(1,k)-xyz(1,j))/dist  
	 rij0(2) = (xyzL(2,k)-xyz(2,j))/dist  
	 rij0(3) = (xyzL(3,k)-xyz(3,j))/dist  

	 IF(m.le.3)THEN

	  value =  rij0(m) !could be negative

	 ELSE
	
	  xyzLcom(1) = (xyzL(1,k) - Rcm_lig(1,k_m))
          xyzLcom(2) = (xyzL(2,k) - Rcm_lig(2,k_m))
          xyzLcom(3) = (xyzL(3,k) - Rcm_lig(3,k_m))

	  call AcrossB(unitm,xyzLcom,cross)

	  call AdotB(rij0,cross,value)

	 ENDIF

	 IF(n.le.3)THEN

	  term =  rij0(n) !could be negative

	 ELSE
	
	  xyzLcom(1) = xyzL(1,k) - Rcm_lig(1,k_n)
          xyzLcom(2) = xyzL(2,k) - Rcm_lig(2,k_n)
          xyzLcom(3) = xyzL(3,k) - Rcm_lig(3,k_n)

	  call AcrossB(unitn,xyzLcom,cross)

	  call AdotB(rij0,cross,term)

	 ENDIF

          F(ider) =  F(ider) + sign*factor*term*value 

451	continue
        enddo ! loop is = 1, Ninter

c and then also the appropriate inter ligand-ligand pairs:
         do k = 1, NligNA ! and include NBI of inter ligand-ligand NBI...
         do j = k, NligNA  !due to selection rule immediately below
         if(NBILL(k,j).ne.0)then
         
c select subset that contributes
          IF( (OnOff(k).eq.1.AND.OnOff(j).ne.1) .OR.
     1        (OnOff(j).eq.1.AND.OnOff(k).ne.1) ) THEN

	   factor = CLL(k,j)

           dist = dsqrt( (xyzL(1,k)-xyzL(1,j))**2 +
     1                   (xyzL(2,k)-xyzL(2,j))**2 +
     2                   (xyzL(3,k)-xyzL(3,j))**2 )

           rij0(1) = (xyzL(1,k)-xyzL(1,j))/dist
           rij0(2) = (xyzL(2,k)-xyzL(2,j))/dist
           rij0(3) = (xyzL(3,k)-xyzL(3,j))/dist

           IF(m.le.3)THEN

            value =  rij0(m) !could be negative

           ELSE

            xyzLcom(1) = (xyzL(1,k) - Rcm_lig(1,k_m))
            xyzLcom(2) = (xyzL(2,k) - Rcm_lig(2,k_m))
            xyzLcom(3) = (xyzL(3,k) - Rcm_lig(3,k_m))

            call AcrossB(unitm,xyzLcom,cross)

            call AdotB(rij0,cross,value)

           ENDIF

           IF(n.le.3)THEN

            term =   rij0(n) !could be negative

           ELSE

            xyzLcom(1) = xyzL(1,j) - Rcm_lig(1,k_n)
            xyzLcom(2) = xyzL(2,j) - Rcm_lig(2,k_n)
            xyzLcom(3) = xyzL(3,j) - Rcm_lig(3,k_n)

            call AcrossB(unitn,xyzLcom,cross)

            call AdotB(rij0,cross,term)

          ENDIF ! n.le.3

            F(ider) =  F(ider) + factor*term*value

	  ENDIF ! selection rule with OnOff

	 endif ! for contributing NBILL pairs
	 enddo
	 enddo



	else  ! two dof belong to two ligands
	
         do k = 1, NligNA ! and include NBI of inter ligand-ligand NBI...
         do j = k, NligNA  !due to selection rule immediately below
         if(NBILL(k,j).ne.0)then

	
c select subset that contributes
          IF((OnOff(k).eq.1.and.OnOff2(j).eq.1).OR. !possible problem with sign here: (k,j) vs ((j,k)
     1       (OnOff(j).eq.1.and.OnOff2(k).eq.1))THEN

	  
           factor = CLL(k,j)

	   dist = dsqrt( (xyzL(1,k)-xyzL(1,j))**2 + 
     1                   (xyzL(2,k)-xyzL(2,j))**2 +
     2	                 (xyzL(3,k)-xyzL(3,j))**2 )
 
	   rij0(1) = (xyzL(1,k)-xyzL(1,j))/dist  
	   rij0(2) = (xyzL(2,k)-xyzL(2,j))/dist  
	   rij0(3) = (xyzL(3,k)-xyzL(3,j))/dist  

	   IF(m.le.3)THEN

	    value =  rij0(m) !could be negative

	   ELSE
	
	    xyzLcom(1) = (xyzL(1,k) - Rcm_lig(1,k_m))
            xyzLcom(2) = (xyzL(2,k) - Rcm_lig(2,k_m))
            xyzLcom(3) = (xyzL(3,k) - Rcm_lig(3,k_m))

	    call AcrossB(unitm,xyzLcom,cross)

	    call AdotB(rij0,cross,value)

	   ENDIF
	   IF(n.le.3)THEN

	    term =  + rij0(n) !opposite of value above?

	   ELSE
	
	    xyzLcom(1) = xyzL(1,j) - Rcm_lig(1,k_n)
            xyzLcom(2) = xyzL(2,j) - Rcm_lig(2,k_n)
            xyzLcom(3) = xyzL(3,j) - Rcm_lig(3,k_n)

	    call AcrossB(unitn,xyzLcom,cross)

	    call AdotB(rij0,cross,term)

	  ENDIF ! end selection on n.le.3

            F(ider) =  F(ider) - factor*term*value 

	  ENDIF ! end selection on OnOff and OnOff2

	 endif ! end selection on active NBILL
	 enddo
	 enddo

	endif  ! end k_m.ne.k_n

        Enddo  !end of n dof
        ENDDO ! end of k_n ligand column

        Enddo ! end of m dof
        ENDDO ! of endo k_m ligand row
c                                                             c
c                                                             c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


	return
	end

