	subroutine num_F_prtn5(xyzP,Nprtn,IndexMC,Nres,IndexXi,EndXi,
     1Nxi,XiRes,NBIPP,P0,Nlig,NligNA,Lig,rLig0,NBIPL,Ninter,vdWHET,
     2mHET,delta,type,F,NF,Pyn,CPP,CPL,NBILL,CLL,NinterL)

	implicit none
c INPUT:
c set largest possible # atoms in all ligands(NligNA) to 2500,
c and set total possible # ligands(Nlig) to 1000
	integer Nprtn,Nres,Nxi,NF,Nlig,NligNA
	real*8 xyzP(3,Nprtn),rLig0(3,2500),vdWHET(2500)
	real*8 CPP(Nprtn,Nprtn),CPL(10000),CLL(2500,2500)
	integer IndexMC(2*Nres,4),Lig(1000,2)
	integer IndexXi(Nxi,4),EndXi(Nxi),XiRes(Nxi)
	integer NBIPP(Nprtn,Nprtn),NBIPL(10000,2),NBILL(2500,2500)
	integer Pyn(2*Nres+Nxi)
	character*1 type(Nprtn),type1,type2
	real*8 delta
	integer Ninter,NinterL
        real*8 P0(2*Nres+Nxi),RotStiff(2*Nres+Nxi)
        real*8 phi,psi,xi,fphi,fpsi,fxi,ProPhi
c   used in computation:
	integer m,n,m2,m3,n3,i,j,i1,i2,ij,ider,ier
	integer iLig,k,is,ks,ls,k_m,k_n,nlo,icount
	integer klimit,k0_m,k0_n,kl_m,kl_n,inc_m,inc_n,inc_e
	integer istart,jend,jlimit,iatom,atomafter
	integer OnOff(2500),OnOff2(2500)
	real*8 dist,energy,dEdq_m,dEdq_n,sum_mn,num(3)
	real*8 factor
        real*8 xyzm(3,Nprtn),xyzn(3,Nprtn),xyzmn(3,Nprtn)
	real*8 nhat(3,2*Nres+Nxi)
	real*8 rLign(3,2500),rLigm(3,2500),rLigmn(3,2500)
	real*8 mHET(2500),Rcm_lig(3,1000),TmassL(1000)
c OUTPUT:
	real*8 F(NF)



c 8/24/2016 P0 from index_prtn, includes proline//call Pvec(Nprtn,Nres,Nxi,IndexMC,IndexXi,xyzP,P0)
c get vector of dihedral stiffness constants:
c(for explanation see notes dated 10/20/2015)
c first for the MC phi/psi pairs:

c N terminal pair of psi/phi energies are done separately:
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
	  Rcm_lig(1,k) = Rcm_lig(1,k) + mHET(i)*rLig0(1,i)
	  Rcm_lig(2,k) = Rcm_lig(2,k) + mHET(i)*rLig0(2,i)
	  Rcm_lig(3,k) = Rcm_lig(3,k) + mHET(i)*rLig0(3,i)
	enddo
	 Rcm_lig(1,k) = Rcm_lig(1,k)/TmassL(k)
	 Rcm_lig(2,k) = Rcm_lig(2,k)/TmassL(k)
	 Rcm_lig(3,k) = Rcm_lig(3,k)/TmassL(k)
	ELSE
	 i = i + 1
	 Rcm_lig(1,k) = rLig0(1,i)
	 Rcm_lig(2,k) = rLig0(2,i)
	 Rcm_lig(3,k) = rLig0(3,i)
	ENDIF
	ENDDO




c  get nhat vector:
	do m = 1,2*Nres     
	 nhat(1,m) = 0.0d00
	 nhat(2,m) = 0.0d00
	 nhat(3,m) = 0.0d00
	enddo

	do m = 1,2*Nres     
	 if(Pyn(m).eq.0)goto10
	 m3 = IndexMC(m,3)
	 m2 = IndexMC(m,2)
	 num(1) = xyzP(1,m3) - xyzP(1,m2)
	 num(2) = xyzP(2,m3) - xyzP(2,m2)
	 num(3) = xyzP(3,m3) - xyzP(3,m2)
	 dist = dsqrt( num(1)**2 + num(2)**2 + num(3)**2 )
	 nhat(1,m) = num(1)/dist
	 nhat(2,m) = num(2)/dist
	 nhat(3,m) = num(3)/dist
10	 continue
	enddo
	
	do m = 1,Nxi
	 m3 = IndexXi(m,3)
	 m2 = IndexXi(m,2)
	 num(1) = xyzP(1,m3) - xyzP(1,m2)
	 num(2) = xyzP(2,m3) - xyzP(2,m2)
	 num(3) = xyzP(3,m3) - xyzP(3,m2)
	 dist = dsqrt( num(1)**2 + num(2)**2 + num(3)**2 )
	 nhat(1,2*Nres+m) = num(1)/dist
	 nhat(2,2*Nres+m) = num(2)/dist
	 nhat(3,2*Nres+m) = num(3)/dist
	enddo

c==================================================================
c Compute INTRA-prtn potential (varying protein torsions):
c==================================================================

c   MAIN LOOP:
	 ider = 0
	DO m = 1,2*Nres    ! BLOCK I DIAGONAL ELEMENTS: MC-MC, 1st column
	 m3 = IndexMC(m,3)

	 ider = ider+1
	 F(ider) = 0.0d00
	 dEdq_m = 0.0d00
	 sum_mn = 0.0d00  !edited in 9 Oct 07

	 if(Pyn(m).eq.0)goto108 ! for proline-Phi"s and Nterminal-Phi


         call updateMC(xyzP,xyzm,IndexMC,nhat,m,Nres,Nxi,
     1   Nprtn,delta)

c do case of F(m,m) seperately:
         call updateMC(xyzm,xyzmn,IndexMC,nhat,m,Nres,Nxi,
     1   Nprtn,delta)

	do i = 1, m3-1
	do j = m3+1, Nprtn
	 if(NBIPP(i,j).eq.1)then

	  factor = CPP(i,j)
	  call Enrgy2(xyzP,xyzmn,Nprtn,i,j,factor,energy)
	  sum_mn = sum_mn + energy  
	  call Enrgy2(xyzP,xyzm,Nprtn,i,j,factor,energy)
	  dEdq_m = dEdq_m + energy  

	 endif
	enddo  
	enddo  

c include all ligand interactions, if they exist:
	do iLig = 1, Ninter
	 i = NBIPL(iLig,1)  ! belongs to ligand subset
	 j = NBIPL(iLig,2)  ! belongs to protein subset
	 factor = CPL(iLig)
	 call Enrgy_inter2(xyzP,xyzmn,Nprtn,j,rLig0,rLig0,
     1	  i,factor,energy)
	 sum_mn = sum_mn + energy
	 call Enrgy_inter2(xyzP,xyzm,Nprtn,j,rLig0,rLig0,
     1	  i,factor,energy)
	 dEdq_m = dEdq_m + energy
	enddo

         F(ider) = (sum_mn - 2.0d00*dEdq_m)/(delta*delta)

 	 F(ider) = F(ider) + RotStiff(m)  !include stiffness of MC phi/psi dihedral dof

108	 continue


	DO n = m+1,2*Nres  ! REMAINDER OF BLOCK I: MC-MC, first column
	 n3 = IndexMC(n,3)

	 ider = ider+1
	 F(ider) = 0.0d00
	 dEdq_m = 0.0d00
	 dEdq_n = 0.0d00
	 sum_mn = 0.0d00  !edited in 9 Oct 07
	 if(Pyn(m).eq.0.or.Pyn(n).eq.0)goto110 ! for proline-Phi"s and Nterminal-Phi

         call updateMC(xyzP,xyzn,IndexMC,nhat,n,Nres,Nxi,
     1   Nprtn,delta)
         call updateMC(xyzm,xyzmn,IndexMC,nhat,n,Nres,Nxi,
     1   Nprtn,delta)

	do i = 1, m3-1
	do j = n3+1,Nprtn	
	 if(NBIPP(i,j).eq.1)then

	 factor = CPP(i,j)
	 call Enrgy2(xyzP,xyzmn,Nprtn,i,j,factor,energy)
	 sum_mn = sum_mn + energy  
	 call Enrgy2(xyzP,xyzm,Nprtn,i,j,factor,energy)
	 dEdq_m = dEdq_m + energy  
	 call Enrgy2(xyzP,xyzn,Nprtn,i,j,factor,energy)
	 dEdq_n = dEdq_n + energy  

	 endif
	enddo 
	enddo 

c include all ligand interactions, if Ligands exist:
	do iLig = 1, Ninter
	 i = NBIPL(iLig,1)
	 j = NBIPL(iLig,2)
	 factor = CPL(iLig)
	 call Enrgy_inter2(xyzP,xyzmn,Nprtn,j,rLig0,rLig0,
     1	  i,factor,energy)
	 sum_mn = sum_mn + energy
	 call Enrgy_inter2(xyzP,xyzm,Nprtn,j,rLig0,rLig0,
     1	  i,factor,energy)
	 dEdq_m = dEdq_m + energy
	call Enrgy_inter2(xyzP,xyzn,Nprtn,j,rLig0,rLig0,
     1	i,factor,energy)
	xi = energy
	 dEdq_n = dEdq_n + energy
	enddo


         F(ider) = (sum_mn - dEdq_m - dEdq_n)/(delta*delta)

110	 continue
	ENDDO  		! END INNER SUMMATION OVER ALL MAINCHAIN DOF


	DO n = 1,Nxi	! BLOCK II: MC-SC, 1st column
	 ider = ider+1
	 F(ider) = 0.0d00
	 dEdq_m = 0.0d00
	 dEdq_n = 0.0d00
	 sum_mn = 0.0d00  !edited in 9 Oct 07

	 if(Pyn(m).eq.0)goto210  ! n2,n3 should never be zero...
	 n3 = indexXi(n,3)

         call updateSC(xyzP,xyzn,IndexXi,EndXi,nhat,n,Nres,
     1    Nxi,Nprtn,delta)
         call updateSC(xyzm,xyzmn,IndexXi,EndXi,nhat,n,Nres,
     1    Nxi,Nprtn,delta)

	 IF(m3.le.IndexXi(n,2))THEN
	
	 do i = 1, m3-1
	 do j = n3+1,EndXi(n)
	  if(NBIPP(i,j).eq.1)then

	 factor = CPP(i,j)
	 call Enrgy2(xyzP,xyzmn,Nprtn,i,j,factor,energy)
	 sum_mn = sum_mn + energy  
	 call Enrgy2(xyzP,xyzm,Nprtn,i,j,factor,energy)
	 dEdq_m = dEdq_m + energy  
	 call Enrgy2(xyzP,xyzn,Nprtn,i,j,factor,energy)
	 dEdq_n = dEdq_n + energy  

	  endif 
	 enddo
	 enddo

	ELSE

	 do i = n3+1,EndXi(n)
	 do j = m3+1,Nprtn
	  if(NBIPP(i,j).eq.1)then

	 factor = CPP(i,j)
	 call Enrgy2(xyzP,xyzmn,Nprtn,i,j,factor,energy)
	 sum_mn = sum_mn + energy  
	 call Enrgy2(xyzP,xyzm,Nprtn,i,j,factor,energy)
	 dEdq_m = dEdq_m + energy  
	 call Enrgy2(xyzP,xyzn,Nprtn,i,j,factor,energy)
	 dEdq_n = dEdq_n + energy  

	  endif
	 enddo
	 enddo

	ENDIF


c include all ligand interactions, if Ligands exist:
	do iLig = 1, Ninter
	 i = NBIPL(iLig,1)
	 j = NBIPL(iLig,2)
	 factor = CPL(iLig)
	 call Enrgy_inter2(xyzP,xyzmn,Nprtn,j,rLig0,rLig0,
     1	  i,factor,energy)
	 sum_mn = sum_mn + energy
	 call Enrgy_inter2(xyzP,xyzm,Nprtn,j,rLig0,rLig0,
     1	  i,factor,energy)
	 dEdq_m = dEdq_m + energy
	call Enrgy_inter2(xyzP,xyzn,Nprtn,j,rLig0,rLig0,
     1	  i,factor,energy)
	 dEdq_n = dEdq_n + energy
	enddo

         F(ider) = (sum_mn - dEdq_m - dEdq_n)/(delta*delta)

210	 continue
	ENDDO  ! END INNER LOOP OVER ALL SIDECHAIN DOF
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                      added 9/14                             c 
c                    Blocks in FIRST column                   c
c            for rigid-body D.o.F.s of ligands                c 

	do k = 1, Nlig  ! BLOCKS IV: MC-RB, first column for EACH ligand...

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
	 dEdq_m = 0.0d0
	 dEdq_n = 0.d0
	 sum_mn = 0.d0
	 if(Pyn(m).eq.0)goto215

	 call updateRB(rLig0,rLign,n,OnOff,Nlig,NligNA,
     1   Rcm_lig,k,delta)
	
	do is = 1, Ninter
	 i = NBIPL(is,1)
	 if(OnOff(i).eq.0)goto445
	 j = NBIPL(is,2)    ! j labels protein atoms
	 factor = CPL(is)

	 call Enrgy_inter2(xyzP,xyzm,Nprtn,j,rLig0,rLign,i,
     1   factor,energy)
	 sum_mn = sum_mn + energy
	 call Enrgy_inter2(xyzP,xyzm,Nprtn,j,rLig0,rLig0,i,
     1	 factor,energy)
	 dEdq_m = dEdq_m + energy  
	 call Enrgy_inter2(xyzP,xyzP,Nprtn,j,rLig0,rLign,i,
     1	 factor,energy)
	 dEdq_n = dEdq_n + energy  

445	continue
	enddo ! loop is = 1, Ninter 

	F(ider) = (sum_mn - dEdq_m - dEdq_n)/(delta*delta)

215	continue
	enddo  ! loop over n = 1,...6
	enddo  !loop over k ligands
   
c                                                             c 
c                                                             c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 	ENDDO  ! END OUTER LOOP OVER ALL MAINCHAIN DOF


	DO m = 1, Nxi  		! BLOCK III: SC-SC, SECOND Column
	 m3 = indexXi(m,3)

         call updateSC(xyzP,xyzm,IndexXi,EndXi,nhat,m,Nres,
     1    Nxi,Nprtn,delta)

	DO n = m,Nxi 		! INNER LOOP OVER SC DOF
	 n3 = IndexXi(n,3)

         call updateSC(xyzP,xyzn,IndexXi,EndXi,nhat,n,Nres,
     1    Nxi,Nprtn,delta)
         call updateSC(xyzm,xyzmn,IndexXi,EndXi,nhat,n,Nres,
     1    Nxi,Nprtn,delta)
	 ider = ider+1
	 F(ider) = 0.0d00
	 dEdq_m = 0.0d00
	 dEdq_n = 0.0d00
	 sum_mn = 0.0d00  !edited in 9 Oct 07

	IF(EndXi(m).ne.EndXi(n))THEN

	 do i = m3+1,EndXi(m)
	 do j = n3+1,EndXi(n)
	  if(NBIPP(i,j).eq.1)then

	   factor = CPP(i,j)
	   call Enrgy2(xyzP,xyzmn,Nprtn,i,j,factor,energy)
	   sum_mn = sum_mn + energy  
	   call Enrgy2(xyzP,xyzm,Nprtn,i,j,factor,energy)
	   dEdq_m = dEdq_m + energy  
	   call Enrgy2(xyzP,xyzn,Nprtn,i,j,factor,energy)
	   dEdq_n = dEdq_n + energy  

	  endif
	 enddo 
	 enddo 

	ELSE

	 do i = 1, indexXi(m,2)
	 do j = n3+1,EndXi(n)
	  if(NBIPP(i,j).eq.1)then

	   factor = CPP(i,j)
	   call Enrgy2(xyzP,xyzmn,Nprtn,i,j,factor,energy)
	   sum_mn = sum_mn + energy  
	   call Enrgy2(xyzP,xyzm,Nprtn,i,j,factor,energy)
	   dEdq_m = dEdq_m + energy  
	   call Enrgy2(xyzP,xyzn,Nprtn,i,j,factor,energy)
	   dEdq_n = dEdq_n + energy  

	  endif
	 enddo 
	 enddo 

	 do i = n3+1,EndXi(n)
	 do j = EndXi(n)+1,Nprtn
	  if(NBIPP(i,j).eq.1)then

	   factor = CPP(i,j)
	   call Enrgy2(xyzP,xyzmn,Nprtn,i,j,factor,energy)
	   sum_mn = sum_mn + energy  
	   call Enrgy2(xyzP,xyzm,Nprtn,i,j,factor,energy)
	   dEdq_m = dEdq_m + energy  
	   call Enrgy2(xyzP,xyzn,Nprtn,i,j,factor,energy)
	   dEdq_n = dEdq_n + energy  

	  endif
	 enddo 
	 enddo 

	ENDIF

c include ligand interactions
	do i = 1, Ninter
	 k = NBIPL(i,1)
	 j = NBIPL(i,2)    ! j labels protein atoms
	 factor = CPL(i)

	 call Enrgy_inter2(xyzP,xyzmn,Nprtn,j,rLig0,rLig0,k,
     1   factor,energy)
	 sum_mn = sum_mn + energy
	 call Enrgy_inter2(xyzP,xyzm,Nprtn,j,rLig0,rLig0,k,
     1	 factor,energy)
	 dEdq_m = dEdq_m + energy  
	 call Enrgy_inter2(xyzP,xyzn,Nprtn,j,rLig0,rLig0,k,
     1	 factor,energy)
	 dEdq_n = dEdq_n + energy  
	enddo ! loop over Ninter prtn-Ligand interactions

         F(ider) = (sum_mn - dEdq_m - dEdq_n)/(delta*delta)

         if(m.eq.n)F(ider) = F(ider) + RotStiff(2*Nres+m) !MC dihedral stiffnesses

	ENDDO  ! END INNER LOOP OVER ALL SIDECHAIN DOF
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                      added 9/14                             c 
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
	if(Lig(k,1).gt.1)ks=6
	do n = 1,ks
	 ider = ider + 1
	 F(ider) = 0.0d0
	 dEdq_m = 0.0d0
	 dEdq_n = 0.d0
	 sum_mn = 0.d0

	 call updateRB(rLig0,rLign,n,OnOff,Nlig,NligNA,
     1Rcm_lig,k,delta)
	
	do is = 1, Ninter
	i = NBIPL(is,1)  ! labels ligand atom
	if(OnOff(i).eq.0)goto447   !include only 'active' ligands...
	j = NBIPL(is,2)    ! j labels protein atoms
	factor = CPL(is)

	call Enrgy_inter2(xyzP,xyzm,Nprtn,j,rLig0,rLign,i,
     1	factor,energy)
	sum_mn = sum_mn + energy
	call Enrgy_inter2(xyzP,xyzm,Nprtn,j,rLig0,rLig0,i,
     1	factor,energy)
	dEdq_m = dEdq_m + energy  
	call Enrgy_inter2(xyzP,xyzP,Nprtn,j,rLig0,rLign,i,
     1	factor,energy)
	dEdq_n = dEdq_n + energy  

447	continue
	enddo ! end of loop over all inter-protein-ligand pairs

	F(ider) = (sum_mn - dEdq_m - dEdq_n)/(delta*delta)

	enddo  ! loop over n = 1,...6
	enddo  !loop over k ligands
c                                                             c 
c                                                             c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	ENDDO  ! END OUTER LOOP OVER ALL SIDECHAIN DOF


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                      added 9/14                             c 
c                Blocks AFTER  MC & SC D.o.F.                 c
c            for rigid-body D.o.F.s of ligands                c 
c     12/2016 added inter-ligand contributions to Hessian

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
	
	ks = 3
	if(Lig(k_m,1).gt.1)ks=6
	Do m = 1,ks

	call updateRB(rLig0,rLigm,m,OnOff,Nlig,NligNA,
     1Rcm_lig,k_m,delta)


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

	ks = 3
	if(Lig(k_n,1).gt.1)ks=6
	Do n = nlo,ks
	ider = ider + 1

	call updateRB(rLig0,rLign,n,OnOff2,Nlig,NligNA,
     1Rcm_lig,k_n,delta)
	call updateRB(rLigm,rLigmn,n,OnOff2,Nlig,NligNA,
     1Rcm_lig,k_n,delta)

	F(ider) = 0.0d0
	dEdq_m = 0.0d0
	dEdq_n = 0.d0
	sum_mn = 0.d0

	if(k_m.eq.k_n)then !within single ligand, L:

	 do i = 1, Ninter ! include NBI of inter protein-ligands
	  k = NBIPL(i,1) ! ligand atom
	  IF(OnOff(k).eq.1)THEN !all active ligand atoms
	  j = NBIPL(i,2) ! protein atom
	  factor = CPL(i)

	  call Enrgy_inter2(xyzP,xyzP,Nprtn,j,rLig0,rLigmn,k,
     1	  factor,energy)
	  sum_mn = sum_mn + energy
	  call Enrgy_inter2(xyzP,xyzP,Nprtn,j,rLig0,rLigm,k,
     1	  factor,energy)
	  dEdq_m = dEdq_m + energy  
	  call Enrgy_inter2(xyzP,xyzP,Nprtn,j,rLig0,rLign,k,
     1	  factor,energy)
	  dEdq_n = dEdq_n + energy  
	  ENDIF
	 enddo ! end of loop over all inter-protein-ligand pairs


	 do k = 1, NligNA ! and include NBI of inter ligand-ligand NBI...
	 do j = k, NligNA  !due to selection rule immediately below
	 if(NBILL(k,j).ne.0)then
	 
c select subset that contributes
	  IF( (OnOff(k).eq.1.AND.OnOff(j).ne.1) .OR.
     1        (OnOff(j).eq.1.AND.OnOff(k).ne.1) ) THEN
	
	  factor = CLL(k,j)

          call Enrgy_interL(rLig0,rLigmn,j,k,factor,energy)
	 phi = energy
	  sum_mn = sum_mn + energy
	  call Enrgy_interL(rLig0,rligm,j,k,factor,energy)
	 psi = energy
	  dEdq_m = dEdq_m + energy  
	  call Enrgy_interL(rLig0,rlign,j,k,factor,energy)
	 xi = energy
	  dEdq_n = dEdq_n + energy  
	  ENDIF


	 endif
	 enddo
	 enddo

	 F(ider) = (sum_mn - dEdq_m - dEdq_n)/(delta*delta)
	
	else  ! k_m .NE. k_n ; or within separate ligands

	 do k = 1, NligNA ! and include NBI of inter ligand-ligand NBI...
	 do j = k, NligNA  !due to selection rule immediately below
	 if(NBILL(k,j).ne.0)then

	  factor = CLL(k,j)

c select subset that contributes
	  IF((OnOff(k).eq.1.and.OnOff2(j).eq.1).OR.
     1       (OnOff(j).eq.1.and.OnOff2(k).eq.1))THEN

          call Enrgy_interL(rLig0,rLigmn,j,k,factor,energy)
	phi = energy
	  sum_mn = sum_mn + energy
	  call Enrgy_interL(rLig0,rligm,j,k,factor,energy)
	psi = energy
	  dEdq_m = dEdq_m + energy  
	  call Enrgy_interL(rLig0,rlign,j,k,factor,energy)
	xi = energy
	  dEdq_n = dEdq_n + energy  
	  ENDIF

	 endif
	 enddo
	 enddo

	 F(ider) = (sum_mn - dEdq_m - dEdq_n)/(delta*delta)

	endif


	Enddo !end of n dof
	ENDDO ! end of k_n ligand column

	Enddo ! end of m dof
	ENDDO ! of endo k_m ligand row

c                                                             c 
c                                                             c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



	return
	end


