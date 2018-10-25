	PROGRAM ATMAN5

	implicit none
	integer Nder,itor,jtor,j0,i0,ij,Nzero,i,j,k
	integer check,icount,Nslow,i1,i2,i3,i4
	integer NDOF,NF,NFr,Nprtn,Nres,Nxi,NintraP
	integer Npro,Nmodes,power,mutate,Ninter,ijss
	integer Nlig,NligNA,NT,NI,NH,N_ss,NinterL,NChain,HOH
	integer rscale,NMAlig
	character*1 AnaNum

   	parameter ( Nprtn = 2598 ) !# of atoms in protein molecule
   	parameter ( Nres  =  302 ) !# of residues
   	parameter ( NPro  =  14 ) !# of nonaminoterminal-prolines 
   	parameter ( Nxi   =  478 ) !# of sidechain, Xi, torsions
 	parameter ( rscale = 0 ) ! ID's protein entry's first residue number
 	parameter ( NChain = 1 )  !currently 2/15, not developed yet, but will need to
 	parameter ( AnaNum = 'A' ) ! whether it's a numeric or analytic F computation

c LIGAND INFORMATION
	parameter ( Nlig  =  0 )  !# ligands, should equal # residues for HETATMs
	parameter ( HOH  =  0 )  !If present, labels first HOH entry's ligand #
	parameter ( NligNA  = 0 )  ! total # atoms in all ligands <= 2500, else up
	parameter ( NMAlig = 0 ) ! Number of Multi-Atom ligands
	parameter ( Ninter = 0 )  ! # of inter-protein-ligand NCB <=10000, else up
	parameter ( NinterL  = 0  )  ! #of inter ligand-ligand NBIs, <=5000, else up
c other params:
	parameter ( mutate = 0 ) ! manually block select D.o.F.s
	parameter ( Nmodes = 100  ) ! # of modes to compute
	parameter ( NDoF = 2*Nres+Nxi+NMAlig*6+(Nlig-NMAlig)*3 )  

	parameter ( NF  = NDoF*(NDoF+1)/2  )   
	parameter ( NFr = (NDoF-NPro-1-mutate)*(NDoF-Npro-mutate)/2 ) 

c variables to index  mainchain torsions:
	integer IndexMC(2*Nres,4),Block(100)
c variables to index sidechain torsions:
	integer IndexXi(Nxi,4),EndXi(Nxi),XiRes(Nxi)
c variables to identify atoms, residues, and coordinates:
	real*8 xyzP(3,Nprtn),massP(Nprtn),Rhc(Nprtn)
	real*8 Bhet(2500),B(Nprtn),P(2*Nres+Nxi)
	integer rlist(Nprtn),rnum(Nres),Pyn(2*Nres+Nxi)
	character*1 type(Nprtn),typei,typej
	character*4 rsort(Nprtn),rsrt
	character*5 asrt,asort(Nprtn)
c variables to index interatomic interactions:
	integer IntAct(Nprtn,Nprtn),NBIPP(Nprtn,Nprtn)
	real*8 CPP(Nprtn,Nprtn),epsi,epsj,eps
	real*8 sumF,sumH,Tmass
	real*8 C,C2,C3,SS
	real*8 vdWcut,delta
c HETATM parameters
	integer rlstHET(2500),Lig(1000,2)   !WARNING # LIGANDS LIMITED TO 1000 BY THIS
	character*5 aHET(2500)
	character*4 typeHET(2500)
	real*8 xyzL(3,2500),mHET(2500),vdWHET(2500)
	real*8 epsHET(2500)
	integer NBIPL(10000,2),NBILL(2500,2500)
	real*8 CPL(10000),CLL(2500,2500)
	integer LigNBI(2500)

 	parameter ( delta = 1.0d-05  )  ! numeric update parameter
 	parameter ( vdWcut = 0.00d00 )  !for subs
	parameter ( C = 1.d0 ) 		! intra Protein NBI adjuster
	parameter ( C2 = 1.d0 )	! inter Protein-Ligand adjuster
	parameter ( C3 = 1.d0 ) 	! inter Ligand-Ligand adjuster
	parameter ( SS = 0.d0 ) 	! SSBOND strength multiplier
c OUTPUT: 
	integer IntF(NF),IntBin(10)
	real*8 H(NF),F(NF)
	real*8 dRdTor(Nprtn+NligNA,NDoF,3)  ! this should be upped as i add dofs
	integer imax,N,Ncheck,irow,jmax,NligNAr,Nligr
	integer Ivalue,BinSize,Bbin(20)
	character*6 atom
	real*8 minH,minF,maxH,maxF
	real*8 Evalues(NDoF-Npro-1-mutate)
	real*8 Evectors(NDoF-Npro-1-mutate,NDoF-Npro-1-mutate)
        real*8 Ac(NDoF-Npro-1-mutate,NDoF-Npro-1-mutate),dave
	real*8 Bc(NDoF-Npro-1-mutate,NDoF-Npro-1-mutate)
	real*8 temp(NDoF,NDoF),Fr(NFr),Hr(NFr)
        real*8 lhs,rhs,Rmax,lhsmax,rhsmax,diff,ratio
	real*8 Rmin,dist,factor,Ortho,Sum(NFr)
	real*8 x,num,Rrange,Rstar
	integer Bin(20),BinStart,iBin


c CURRENTLY MUST ENTER LIGAND INFORMATION MANUALLY:
c enter #atoms for each of these ligands:
c this should be in parameter list too, above.
	if(Nlig.ne.0)then
 	 Lig(1,1) =  1     ! # of atoms in ligand 1 
 	 Lig(1,2) = 373  ! Residue # for ligand 1
c 	 Lig(2,1) =  31     ! # of atoms in ligand 2 
c 	 Lig(2,2) = 375  ! Residue # for ligand 2
c 	 Lig(3,1) =  1     ! # of atoms in ligand 2 
c 	 Lig(3,2) = 305  ! Residue # for ligand 2
c 	 do i = 3,436 ! # of waters in PDB file (not all may be used)
c 	  Lig(i,1) = 1
c 	  Lig(i,2) = Nres+rscale+2+i
c 	 enddo
	endif

	if(mutate.gt.0)then
c note to be careful w/ Phy-angles of Prolines!
	 do j=1,10
	  block(j)=78+j
	 enddo
	 do j=1,9
	  block(10+j)=89+j
	 enddo
	endif

	write(*,*)' '
	write(*,*)' '
	write(*,*)' 	Beginning computation...'
	write(*,*)' '
c===================================================================
c  index torsions; identify atoms, masses, residues; get coordinates
c===================================================================
	call index_prtn(indexMC,indexXi,massP,xyzP,asort,B,rlist,
     1rsort,rnum,type,EndXi,XiRes,Nres,Nxi,Nprtn,Rhc,ahet,xyzL,
     2rlstHET,mHET,vdWHET,typeHET,epsHET,NligNA,rscale,Bhet,P,Pyn)



	write(6,*)'Returned from index_prtn successfully'
c this translates coords to CoM so comment out normally
c	call InertiaTensor(massP,xyzP,Nprtn,Nres,Nlig)
	write(*,*)
	Tmass=0.d0
	do i = 1, Nprtn
	 if(type(i).ne.'H')Tmass = Tmass+massP(i)
	enddo
	write(*,*)' '
	write(*,*)'Total mass protein is',sngl(Tmass)

	power = 5
c===================================================================
c index all noncovalent (intra-protein chain) intra-atomic distances:
c===================================================================
c  	call Intra_NCB_direct4(Nprtn,xyzP,Rhc,IntAct,NintraP,vdWcut)
  	call Intra_NCB_direct5(Nprtn,xyzP,Rhc,IntAct,NintraP,vdWcut,C)

C TEMPORARY OVERWRITE 3/20/18
c        call Three_Bonds(asort,rsort,rnum,rlist,IntAct,Nres,
c     1  Nprtn,NintraP,Nxi,NBIPP,ijss,IndexXi,EndXi,IndexMC)
c END TEMPORARY OVERWRITE
 
        call Four_Bonds2(asort,rsort,rnum,rlist,IntAct,Nres,
     1  Nprtn,NintraP,Nxi,NBIPP,ijss,IndexXi,EndXi,IndexMC)

        write(*,*)' '
        write(*,*)' CHECKING INTRA-BONDS:'
        call Intra_check2(massp,asort,xyzp,rnum,rlist,rsort,NBIPP,
     1  ijss,nres,Nprtn,Rhc,C)

c obtain matrix R*8 CPP(i,j) of spring constants for  I*4 NBI(i,j)
	N_ss = 0
	Nh = 0
	NI =0
	NT = 0
	do i = 1, Nprtn
	do j = i+1, Nprtn
	if(NBIPP(i,j).eq.1)then

         dist = dsqrt ( (xyzP(1,i) - xyzP(1,j))**2 +
     1                  (xyzP(2,i) - xyzP(2,j))**2 +
     2                  (xyzP(3,i) - xyzP(3,j))**2 )

	 typei = type(i)
	 call epsilons(typei,epsi)
	 typej = type(j)
	 call epsilons(typej,epsj)
	 eps = dsqrt(epsi*epsj)

c	 Rmin = dsqrt(Rhc(i)*Rhc(j))
	 Rmin = C*dsqrt(Rhc(i)*Rhc(j))

	if(dist.lt.0.6d0*Rmin)then
	  NT = NT + 1
	  write(*,*)
     1'WARNING NBI CLOSE:',asort(i),rlist(i),rsort(i),'-',
     2asort(j),rlist(j),rsort(j)
	 if(rlist(i).eq.rlist(j))NI = NI + 1
	 if(asort(i).eq.' H   '.or.asort(j).eq.' H   ')NH=NH+1
	endif

	 if(dist.le.Rmin)then
	  x = dist
	 else
	  Rstar = 1.10868341797d0*Rmin
	  Rrange = Rstar + vdWcut
          num = (dist-Rmin)*(Rstar-Rmin)
          x = Rmin + num/(Rrange-Rmin)
         endif

	 ratio = Rmin/x
	 ratio = ratio * ratio
	 ratio = ratio * ratio * ratio
	 
c	 factor = C*12.d0*eps*(13.d0*ratio*ratio-7.d0*ratio)
c     1/(x*x)
	 factor = 12.d0*eps*(13.d0*ratio*ratio-7.d0*ratio)
     1/(x*x)

	   if(type(i).eq.'S'.and.type(j).eq.'S'.and.
     1dist.le.2.2d0)then
		factor = factor * SS
		write(*,*)' '
          	write(*,*)'S=S bond',asort(i),rlist(i),rsort(i),
     1'-',asort(j),rlist(j),rsort(j)
		write(*,*)' '
		N_ss = N_ss + 1
	   endif


	 CPP(i,j) = factor
	 CPP(j,i) = factor
c TEMPORARY FOR 2JDI SUBUNIT B which misses residues 400-409 in PDB file
c	if((i.ge.3218.and.i.le.3276).OR.(j.ge.3218.and.j.le.3276))then
c	 CPP(i,j)=0.01d0*CPP(i,j)
c	 CPP(j,i)=0.01d0*CPP(i,j)
c	endif
c END TEMPORARY 
 
c TEMPORARY BLOCK 3/20/18
c TIRION POTENTIAL
c	 factor = dist/5.d0
c 	 CPP(i,j) = dexp(-factor*factor)
c 	 CPP(j,i) = dexp(-factor*factor)
c 	 if(type(i).eq.'H'.or.type(j).eq.'H')then
c 		CPP(i,j)=0.d0
c 		CPP(j,i)=0.d0
c 	 endif
c END BLOCK	

	endif
	enddo
	enddo
c TEMPORARY
c	write(*,*)'WARNING: USING TIRION POTENTIAL WITH k=1 EVERYWHERE'
c	write(*,*)' WARNING: RESIDUE 400-409 CONTRIBUTIONS FOR SUB B '
c END TEMPORARY

	write(*,*)'Total near NBI:',NT
	write(*,*)'Total near NBI involving H:',NH
	write(*,*)'Total near NBI intra-residue:',NI
        write(*,*)'Total # of S=S bonds is:',N_ss
        write(*,*)'SS bond stiffness is:',SS
        write(*,*)' '


c===================================================================
c Bin all noncovalent inter-atomic distances:
c===================================================================
	if(HOH.ne.0)then
	k = 0
 	do i = 1,HOH-1
 	 k = k + Lig(i,1)
 	enddo
 	BinSize = 5
 	do i = 1,20
 	 Bbin(i) = 0
 	enddo
 	do i = k+1,NligNA
 	 Ivalue = Bhet(i)/BinSize + 1
 	 Bbin(Ivalue) = Bbin(Ivalue)+1
 	enddo
 	write(*,*)'Bin/Temp   #'
 	do i = 1,20
 	 write(*,*)(i-1)*BinSize,'-',i*BinSize,':',Bbin(i)
 	enddo
 	endif

c if present, cull PDB HOHs 
	Nligr = Nlig
	NligNAr = NligNA
	IF(HOH.ne.0)THEN

	 call CullHOH(xyzP,xyzL,NligNA,Nprtn,Lig,
     1rlist,Rhc,vdWHET,NligNAr,Nligr,HOH,Bhet,LigNBI)

	 if(HOH.ne.0)then
	 k = 0
 	 do i = 1,HOH-1
 	  k = k + Lig(i,1)
 	 enddo
 	 BinSize = 5
 	 do i = 1,20
 	  Bbin(i) = 0
 	 enddo
 	 do i = k+1,NligNAr
 	  Ivalue = LigNBI(i)/BinSize + 1
 	  Bbin(Ivalue) = Bbin(Ivalue)+1
 	 enddo
 	 write(*,*)'Bin/Temp   #'
 	 do i = 1,20
 	  write(*,*)(i-1)*BinSize,'-',i*BinSize,':',Bbin(i)
 	 enddo
 	 endif
c report temporarily:
	 write(*,*)' '
	 write(*,*)'CullHOH reports Nligr=',Nligr
	 write(*,*)' '

	 open(73,file='./DATA/eigen.pdb')
10       format(a6,x,i4,x,a5,a4,1x,i4,4x,3f8.3,2x,f4.2,x,f5.2)
	 atom = 'ATOM  '
	 do i = 1, Nprtn
	  write(73,10)atom,i,asort(i),rsort(i),rscale+rlist(i),
     1(xyzP(k,i),k=1,3),1.0,B(i)
	 enddo
	 atom='HETATM'
	 do i = 1,NligNAr
          write(73,10)atom,Nprtn+i,aHET(i),typeHET(i),
     1Nres+rscale+i,(xyzL(k,i),k=1,3),1.0,Bhet(i)
	 enddo
	 close(73)

	ENDIF

        write(*,*)'TOTAL # of D.o.F. in system (NDoF):',
     12*Nres+Nxi+6*NMAlig+3*(Nlig-NMAlig)

	Tmass=0.d0
	do i = 1, NligNAr
	 Tmass = Tmass + mHET(i)
	enddo
	write(*,*)'Total mass nonprotein is',sngl(Tmass)
	write(*,*)' '

        call Inter_Prtn_Lig(xyzP,xyzL,Nligr,NligNAr,Nprtn,Ninter,
     1NinterL,NBIPL,NBILL,rlist,Rhc,vdWHET,vdWcut)

c section added 11/29/2016
c obtain vector R*8 CPL(Ninter) of spring constants for interligand-protein NBIPL 
        do k = 1, Ninter

	i = NBIPL(k,1) !ligand
	j = NBIPL(k,2) !protein

         dist = dsqrt ( (xyzP(1,j) - xyzL(1,i))**2 +
     1                  (xyzP(2,j) - xyzL(2,i))**2 +
     2                  (xyzP(3,j) - xyzL(3,i))**2 )

         typej = type(j)
         call epsilons(typej,epsj)
	 epsi = epsHET(i)
         eps = dsqrt(epsi*epsj)

         Rmin = dsqrt(vdWHET(i)*Rhc(j))

         if(dist.le.Rmin)then
          x = dist
         else
          Rstar = 1.10868341797d0*Rmin
          Rrange = Rstar + vdWcut
          num = (dist-Rmin)*(Rstar-Rmin)
          x = Rmin + num/(Rrange-Rmin)
         endif

         ratio = Rmin/x
         ratio = ratio * ratio
         ratio = ratio * ratio * ratio

         factor = C2*12.d0*eps*(13.d0*ratio*ratio-7.d0*ratio)
     1/(x*x)

         CPL(k) = factor

        enddo
c section added 11/29/2016 ended

c section added 12/9/2016
c obtain vector R*8 CLL(NinterL) of spring constants for inter-ligand 
        do i = 1, NligNA
	do j = i, NligNA
	IF(NBILL(i,j).ne.0)THEN

         dist = dsqrt ( (xyzL(1,j) - xyzL(1,i))**2 +
     1                  (xyzL(2,j) - xyzL(2,i))**2 +
     2                  (xyzL(3,j) - xyzL(3,i))**2 )

         epsj = epsHET(j)
	 epsi = epsHET(i)
         eps = dsqrt(epsi*epsj)

         Rmin = dsqrt(vdWHET(i)*vdWHET(j))

         if(dist.le.Rmin)then
          x = dist
         else
          Rstar = 1.10868341797d0*Rmin
          Rrange = Rstar + vdWcut
          num = (dist-Rmin)*(Rstar-Rmin)
          x = Rmin + num/(Rrange-Rmin)
         endif

         ratio = Rmin/x
         ratio = ratio * ratio
         ratio = ratio * ratio * ratio

         factor = C3*12.d0*eps*(13.d0*ratio*ratio-7.d0*ratio)
     1/(x*x)

         CLL(i,j) = factor
         CLL(j,i) = factor

	ENDIF
        enddo
        enddo
c section added 12/9/2016 ended

c=========================================================================
c Compute 2nd Derivative Matrix of energy w.r.t. all dihedrals
c=========================================================================

	i0 = 0
	do itor = 1,NDoF
	do jtor = itor,NDoF
	 i0 = i0+1
	 F(i0) = 0.0d00
	 H(i0) = 0.0d00
	enddo
	enddo

	IF(AnaNum.eq.'N'.or.AnaNum.eq.'n')THEN

    	write(*,*)'Compute numeric F:'
    	write(*,*)' delta=',delta

      	call num_F_prtn5(xyzP,Nprtn,IndexMC,Nres,IndexXi,EndXi,
     1Nxi,XiRes,NBIPP,P,Nligr,NligNAr,Lig,xyzL,NBIPL,Ninter,vdWHET,
     2mHET,delta,type,F,NF,Pyn,CPP,CPL,NBILL,CLL,NinterL)

	ENDIF


	IF(AnaNum.eq.'A'.or.AnaNum.eq.'a')THEN

      	write(*,*)' Compute analytic F:'
	write(*,*)' '

        call ana_F_ATMAN(xyzP,IndexMC,IndexXi,EndXi,NBIPP,Nres,Nprtn,
     1Nxi,Ninter,Nlig,Lig,mHET,xyzL,NBIPL,F,NF,type,Pyn,P,NligNA,
     2CPP,CPL,NinterL,NBILL,CLL,IntF)

c bin the # of NBI contributing to each F(i,j):

	do i = 1,10
	 IntBin(i) = 0
	enddo
	do i = 1, NF
	 if(IntF(i).eq.0)IntBin(1)=IntBin(1) + 1
	enddo
	do i = 1, NF
	 if(IntF(i).lt.10)IntBin(2) = IntBin(2) + 1
	enddo
	do i = 1, NF
	 if(IntF(i).lt.100.and.IntF(i).ge.10)
     1IntBin(3) = IntBin(3) + 1
	enddo
	do i = 1, NF
	if(IntF(i).lt.500.and.IntF(i).ge.100)
     1IntBin(4) = IntBin(4) + 1
	enddo
	do i = 1, NF
	if(IntF(i).lt.1000.and.IntF(i).ge.500)
     1IntBin(5) = IntBin(5) + 1
	enddo
	do i = 1, NF
	 if(IntF(i).ge.1000)IntBin(6) = IntBin(6) + 1
	enddo
	write(*,*)' # NBI/F(i,j) histogram:'
	write(*,*)'# 0; 1<=N<10 ; 10<=N<100 ; 100<=N<500 ;
     1 500<=N<1000 ; 1000<=N'
	write(*,*)(IntBin(i),i=1,6)

	ENDIF

c=========================================================================
c Compute Mass Matrix, H
c=========================================================================

c in Oct '07, MMT pulled dRdTor out of subroutines H, since there
c was a segmentation fault when this large array was not announced
c in the main program...
	do i = 1, Nprtn
        do itor = 1,2*Nres+Nxi
         dRdTor(i,itor,1) = 0.0d0
         dRdTor(i,itor,2) = 0.0d0
         dRdTor(i,itor,3) = 0.0d0
        enddo
        enddo

 	write(6,*)'compute H...'
 	call num_H_prtn(massP,xyzP,IndexMC,IndexXi,EndXi,Nres,Nxi,
     1  Nprtn,xyzL,mHET,Nligr,NligNAr,lig,NF,H,delta,dRdTor,NDoF,
     2  Pyn)
 	write(6,*)' numeric H computed...'
	if(power.eq.1)stop

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


	minH = 1000.d0
	maxH = 0.d0
	minF = 1000.d0
	maxF = 0.d0
	do i0 = 1, NDoF*(NDoF+1)/2
	 if(F(i0).ne.0.d0.and.dabs(F(i0)).lt.minF)minF = dabs(F(i0))
	 if(H(i0).ne.0.d0.and.dabs(H(i0)).lt.minH)minH = dabs(H(i0))
	 if(F(i0).ne.0.d0.and.dabs(F(i0)).gt.maxF)maxF = dabs(F(i0))
	 if(H(i0).ne.0.d0.and.dabs(H(i0)).gt.maxH)maxH = dabs(H(i0))
	enddo
	write(6,*)' '
	write(6,*)' Min element of Fij and Hij: ',minF,minH
	write(6,*)' Max element of Fij and Hij: ',maxF,maxH
	write(6,*)' '

c Collapse prolines and phi(1) out of the temporary 2-d array, temp
        icount = 0
	sumF = 0.d0
        do i = 1, NDoF
        do j = i, NDoF
         icount = icount + 1
         temp(i,j) = F(icount)
	 sumF = sumF + temp(i,j)
        enddo
        enddo
	write(*,*)'Average Hessian entry:',sumF/dble(icount)

c Bin these values:
	BinStart = +9
	Bin(1) = 0  ! for F entries = 0.0
	do ij = 1, NF
	 if(F(ij).eq.0.d0)Bin(1) = Bin(1) + 1
	enddo
	DO iBin = 2,15  ! each bin spans 3 decades
	 Bin(iBin) = 0
	 BinStart = BinStart - 1
	do ij = 1, NF
	 if(dabs(F(ij)).ge.10.d0**(-3*BinStart).AND.
     1      dabs(F(ij)).lt.10.d0**(-3*(BinStart-1)))
     2      Bin(iBin) = Bin(iBin) + 1
	enddo
	ENDDO
	write(*,*)' iBin  ;  #  (bins of 3 decades, starts at 1.^-21)'
	icount = Bin(1)
	BinStart = +9
	 write(*,*)Bin(1),'# of zero entries'
	do iBin = 2,15
	 BinStart = BinStart-1
	 write(*,*)iBin,10.**(-3*BinStart),
     110.**(-3*(Binstart-1)),Bin(iBin)
	 icount = icount + Bin(iBin)
	enddo
	write(*,*)'Sum Bin:',icount,'(should equal:',NF,')'





	icount = 0
	DO i = 1, NDoF
	 If(i.le.2*Nres)Then
		if(Pyn(i).eq.0)goto20
		do k = 1, mutate
		 if(i.eq.Block(k))write(*,*)i
		 if(i.eq.Block(k))goto20
		enddo
	 EndIf
	Do j = i, NDoF
	 If(j.le.2*Nres)Then
		if(Pyn(j).eq.0)goto22
		do k = 1, mutate
		 if(j.eq.Block(k))goto22
		enddo
	 EndIf
	 icount = icount + 1
	 Fr(icount) = temp(i,j)
22	 continue
	EndDo
20	 continue
	ENDDO

        icount = 0
        do i = 1, NDoF
        do j = i, NDoF
         icount = icount + 1
         temp(i,j) = H(icount)
        enddo
        enddo

	icount = 0
	DO i = 1, NDoF
	 If(i.le.2*Nres)Then
		if(Pyn(i).eq.0)goto30
		do k = 1, mutate
		 if(i.eq.Block(k))goto30
		enddo
	 EndIf
	Do j = i, NDoF
	 If(j.le.2*Nres)Then
		if(Pyn(j).eq.0)goto32
		do k = 1, mutate
		 if(j.eq.Block(k))goto32
		enddo
	 EndIf
	 icount = icount + 1
	 Hr(icount) = temp(i,j)
32	 continue
	EndDo
30	 continue
	ENDDO

	write(6,*)' Original and reduced # ofarray elements: '
	write(6,*)NF, icount
	Ncheck = icount

c for diaganalization check:
        icount = 0
	N = NDoF - NPro - 1 - mutate
        do i = 1, N
        do j = i, N
         icount = icount + 1
         Ac(i,j) = Fr(icount)
         Ac(j,i) = Fr(icount)
         Bc(i,j) = Hr(icount)
         Bc(j,i) = Hr(icount)
        enddo
        enddo
	
c check out hessian properties:
	write(*,*)' '
	write(*,*)' Hessian row(s) with nul entries:'
	call hessian_check(Ac,N)
	write(*,*)' '
	write(*,*)' Mass matrix row(s) with nul entries:'
	call hessian_check(Bc,N)
	write(*,*)' '


c diagonalize the eq'n:  F Z = Lambda H Z using LAPACK:
	write(*,*)' Calling Deigen...'
 	call Deigen(Fr,Hr,Ncheck,N,Evalues,Evectors,Nmodes)
	if(power.eq.5)stop
	
c inserted 1/10/15: criterion that Enrgy.f's C-value is such that
c 14% of modes hve frequency value < 20 cm^-1

c	Nslow = 0
c	do i = 1, N
c	 if(Evalues(i)/2.9979d10.le.20.0d0)Nslow=Nslow+1
c	enddo
c	write(*,*)'# modes < 20: ',Nslow
c	write(*,*)' .14 (NDOFR) = ',14*N/100
c




c check diagonalization that F Z = Lambda H Z for each element:
c (i need to work in reduced space here, and use Ac and Bc )
         DO k = 1, Nmodes        !k loop (over N modes)
          Rmax = 0.0d0
 	  imax = 0
 	  jmax = 0
 	  rhsmax = 0.d0
 	  lhsmax =0.d0
 
         Do i = 1, N             !i loop (over rows)
 
          lhs = 0.0d0
          rhs = 0.0d0
          do j = 1, N             !j loop (over column)
 
           lhs = lhs + Ac(i,j) * Evectors(j,k) 
           rhs = rhs + Bc(i,j) * Evectors(j,k) * Evalues(k)
 
          enddo                   !j loop
 
 
          if(dabs(lhs - rhs).gt.Rmax)then
 	   Rmax = lhs - rhs
 	   imax = i
 	   jmax = j
 	   rhsmax = rhs
 	   lhsmax = lhs
 	  endif
 
         EndDo                   !i loop
         ENDDO                   !k loop


c check orthonormality: that Z**T H Z = I

	write(*,*)' '
 	write(6,*)' Checking Evectors orthogonality:'
 	rhs = 0.d0
 	do k = 1, Nmodes		!over k loop
 
 	do i = 1, N			!over i loop
 	 sum(i) = 0.d0
 	do j = 1, N			!over j loop
  	 sum(i) = sum(i) + Bc(i,j) * Evectors(j,k) 
 	enddo				!over j loop
 	enddo				!over i loop
 
 	do j = 1, Nmodes
 	Ortho = 0.0d0
 	do i = 1, N
 	 Ortho = Ortho + sum(i) * Evectors(i,j)
 	enddo
	if(i.eq.j.and.dabs(Ortho-1.d0).gt.0.001d0)
     1	write(6,*)' Evectors ',k,j,'th orthonormality: ',Ortho-1.d0
	if(i.ne.j.and.dabs(Ortho).gt.0.001d0)
     1	write(6,*)' Evectors ',k,j,'th orthonormality: ',Ortho
 	enddo
 
 	rhs = rhs + Ortho*Ortho
 	enddo				!over k loop (Nmodes)
 
c  	write(6,*)' Fit of animation-mode to Hookean-mode, using '
c  	write(6,*)Nmodes,' modes is: ',anint(100.*sqrt(rhs)),'%'
 	write(*,*)' '
 	write(*,*)' All done. Enjoy and have a nice day.'


	stop
	end
