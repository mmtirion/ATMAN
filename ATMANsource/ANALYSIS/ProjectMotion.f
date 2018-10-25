	program ProjectMotion

	implicit none
	character*80 file
	integer*4 Nprtn,Nres,Nxi,NPro,NDoF,mutate,Nlig
	integer NMAlig,NligNA,NDoFr,rscale
	real*8 Kelvin

	parameter ( Nprtn  = 2264 ) !# of atoms in protein molecule
	parameter ( Nres   = 263 ) !# of residues
	parameter ( rscale = 25  ) ! Residue count-start
	parameter ( Nxi    = 465 ) !# of sidechain, Xi, torsions
	parameter ( NPro   = 12 ) !# of non-Nterminal-prolines (no phi torsions)

	parameter ( Nlig   = 0 ) !# of ligands 
	parameter ( NligNA = 0  ) !# of ligand-atoms total
	parameter ( NMAlig = 0  )!# of Multi-Atom ligands

	parameter ( mutate = 0  )    !# of MC dihedral blocks

	parameter ( NDoF   = 2*Nres+Nxi+6*NMAlig+(Nlig-NMAlig)*3 ) 
	parameter ( NDoFr  = NDoF-Npro-1-mutate) ! reduced # of dofs
	parameter ( Kelvin =  293.d0 ) !for setting amplitude of oscillation


	integer*4 IndexMC(2*Nres,4),block(100)
	integer*4 IndexXi(Nxi,4),EndXi(Nxi),XiRes(Nxi)
	real*8 xyzP(3,Nprtn),massP(Nprtn),Bhet(Nprtn),B(Nprtn)
	integer*4 rlist(Nprtn),rnum(Nres),rlst,Pyn(2*Nres+Nxi)
	character*4 rsort(Nprtn),rsrt
	character*1 type(Nprtn)
	character*5 asort(Nprtn),asrt
	character*3 tor
	real*8 lambda,omega,alpha,enhance,amp0,amp,kb,Na
	real*8 A(NDoFr),occ(Nprtn),P(2*Nres+Nxi)
	real*8 xyzP0(3,Nprtn),angle,pi,max,RMS_Ca
c ligands:
	character*5 aHET(2500)
        character*4 typeHET(2500)
	integer is,ks,ls,OnOff(2500)
	integer rlstHET(2500),Lig(2500,2)  !WARNING: MAX # of LIGANDS NOW AT 2500
        real*8 xyzL0(3,2500),mHET(2500),vdWHET(2500)
	real*8 xyzL(3,2500),epsHET(2500)

c for computations:
	character*60 title
	integer itor,iLig,iCa,Loc,icount
	integer k,m2,m3,i,m,n,jj,j,ier,iatom,Nmode,iterm
        real*8 U(3,3),T(3),rms,mass(Nprtn+NLigNA)
	real*8 Tmass,sum,dist,AveDif,MaxDif,Rcm(3)
	real*8 nhatP(3,2*Nres+Nxi),evectors(3,3)
	real*8 TmassL(1000),Rcm_lig(3,1000),x0,y0,z0,delta
	real*8 xyzdP(3,Nprtn),Rhc(Nprtn),x,y,z
	real*8 xyz0(3,Nprtn+NLigNA),xyz(3,Nprtn+NLigNA)
c for output:
	character*6 atom
	character*4 ResSort
	integer ResNum,iSequence,Nsequences

c INPUT LIGAND DATA
c enter #atoms for each of these ligands:
 	if(Nlig.ne.0)then
           Lig(1,1) =  11     ! # of atoms in ligand 1
           Lig(1,2) = 262  ! Residue # for ligand 1
c           Lig(2,1) =  13    ! # of atoms in ligand 2
c           Lig(2,2) = 281   ! Residue # for ligand 2
c 	  do i = 1,Nlig
c 	   Lig(i,1) = 1
c 	   Lig(i,2) = Nres+rscale+i
c 	  enddo
  	endif



c===================================================================
c  index torsions; identify atoms, masses, residues; get coordinates
c===================================================================
	call index_prtn(indexMC,indexXi,massP,xyzP,asort,B,rlist,
     1	rsort,rnum,type,EndXi,XiRes,Nres,Nxi,Nprtn,Rhc,aHET,xyzL0,
     2  rlstHET,mHET,vdWHet,typeHET,epsHET,NligNA,rscale,Bhet,P,Pyn)

	do i = 1, Nprtn
	 xyzP0(1,i) = xyzP(1,i)
	 xyzP0(2,i) = xyzP(2,i)
	 xyzP0(3,i) = xyzP(3,i)
	enddo


c compute Center of Mass coords of ligands for RB rotations:
        i = 0
        DO k = 1, Nlig
	IF(lig(k,1).gt.1)THEN
         TmassL(k) = 0.d0
         Rcm_lig(1,k) = 0.d0
         Rcm_lig(2,k) = 0.d0
         Rcm_lig(3,k) = 0.d0
        do ls = 1, lig(k,1)
          i = i + 1
          TmassL(k) = TmassL(k) + mHET(i)
          Rcm_lig(1,k) = Rcm_lig(1,k) + mHET(i)*xyzL0(1,i)
          Rcm_lig(2,k) = Rcm_lig(2,k) + mHET(i)*xyzL0(2,i)
          Rcm_lig(3,k) = Rcm_lig(3,k) + mHET(i)*xyzL0(3,i)
        enddo
         Rcm_lig(1,k) = Rcm_lig(1,k)/TmassL(k)
         Rcm_lig(2,k) = Rcm_lig(2,k)/TmassL(k)
         Rcm_lig(3,k) = Rcm_lig(3,k)/TmassL(k)
	ELSE
	 i = i + 1
         Rcm_lig(1,k) = xyzL0(1,k)
         Rcm_lig(2,k) = xyzL0(2,k)
         Rcm_lig(3,k) = xyzL0(3,k)
	ENDIF
        ENDDO

c===================================================================
c  get nhat vector:
c===================================================================
        call Nhat(Nres,Nxi,Nprtn,IndexMC,IndexXi,xyzP,nhatP)

c===================================================================
c create single super-array of protein/ligand coordinates
c===================================================================
        do i = 1,Nprtn
         xyz0(1,i) = xyzP(1,i)
         xyz0(2,i) = xyzP(2,i)
         xyz0(3,i) = xyzP(3,i)
        enddo
	do i = 1, NligNA
	 xyz0(1,i+Nprtn) = xyzL0(1,i)
	 xyz0(2,i+Nprtn) = xyzL0(2,i)
	 xyz0(3,i+Nprtn) = xyzL0(3,i)
	enddo


c===================================================================
c set mass array, pointer vectors...
c===================================================================

        Tmass = 0.d00
	do iatom = 1, Nprtn
	 mass(iatom) = massP(iatom)
	 Tmass = Tmass + mass(iatom)
	enddo
        do iatom = 1,NligNA
	 mass(iatom+Nprtn) = mHET(iatom)
         Tmass = Tmass + mass(iatom+Nprtn)
        enddo

C================================================================
C  Read in eigenfrequencies and eigenvectors
C================================================================
	write(*,*)'Which eigenvector?'
	read(*,*)Nmode
        Na = 6.02205d23  ! avogadros #
	kb = 1.3806d-23  ! boltzmann's constant

c  read in eigenvectors A(NDoFr):
        open(81,file=
     1'/Users/monique/Desktop/WORK/NMStuff/ATMAN/DATA/eigenvectors')

        DO i = 1,Nmode
	 read(81,*)m, omega
	 alpha = dsqrt(2.d0*kb*Kelvin)/dabs(omega)
	
         do j = 1,NDoFr
          read(81,*)A(j)
         enddo
	ENDDO

C================================================================
C  Update according to chosen eigenvector:
C================================================================

	open(30,file='./project.pdb')

	enhance = 1.d0
	amp0 = alpha*enhance
	amp = amp0 

	iterm = 0
c DoF: MAINCHAIN TORSIONS OF PROTEIN:
        DO itor = 1,2*Nres  
	 IF(Pyn(itor).eq.1)THEN

	 iterm = iterm + 1
	 delta = amp*A(iterm)

         call updateMC(xyzP,xyzdP,IndexMC,nhatP,itor,Nres,Nxi,
     1                 Nprtn,delta)

	 do i = indexMC(itor,3)+1,Nprtn
	  xyzP(1,i) = xyzdP(1,i)
	  xyzP(2,i) = xyzdP(2,i)
	  xyzP(3,i) = xyzdP(3,i)
	 enddo
         call Nhat(Nres,Nxi,Nprtn,IndexMC,IndexXi,xyzP,nhatP)

	 ENDIF
	ENDDO


c DoF: SIDECHAIN TORSIONS OF PROTEIN:
        DO itor = 1,Nxi  
	  iterm = iterm + 1
	  delta = amp*A(iterm)
          call updateSC(xyzP,xyzdP,IndexXi,EndXi,nhatP,itor,Nres,Nxi,
     1    	      Nprtn,delta)

	  do i = indexXi(itor,3)+1,EndXi(itor)
	   xyzP(1,i) = xyzdP(1,i)
	   xyzP(2,i) = xyzdP(2,i)
	   xyzP(3,i) = xyzdP(3,i)
	  enddo
	ENDDO
	call Nhat(Nres,Nxi,Nprtn,IndexMC,IndexXi,xyzP,nhatP)


c DoF: LIGAND TRANSLATIONAL COORDINATES
	DO iLig = 1, Nlig

c select subset of Ligand array to activate by first initializing OnOff to zero
        do ks = 1, NligNA
         OnOff(ks) = 0
        enddo
c and then activating each ligand subset by setting appropriate OnOff entries to 1
        icount = 0
        do ks = 1,Nlig
        do ls = 1,lig(ks,1)
         icount = icount + 1
         if(ks.eq.iLig)OnOff(icount) = 1
        enddo
        enddo

	ks = 3
	if(lig(iLig,1).gt.1)ks=6
	do k = 1,ks
	 iterm = iterm + 1
	 delta = amp*A(iterm)
	 call updateRB(xyzL0,xyzL,k,OnOff,NLigNA,
     1Rcm_lig,iLig,delta)
	enddo
	ENDDO

c combine into single array:
	do i = 1, Nprtn
	 xyz(1,i) = xyzP(1,i)
	 xyz(2,i) = xyzP(2,i)
	 xyz(3,i) = xyzP(3,i)
	enddo
	do i = 1, NligNA
	 xyz(1,i+Nprtn)=xyzL(1,i)
	 xyz(2,i+Nprtn)=xyzL(2,i)
	 xyz(3,i+Nprtn)=xyzL(3,i)
	enddo

C================================================================
c Eliminate overall rotation and translations...
C================================================================
  	do i = 1, Nprtn+NligNA
  	 mass(i) = dsqrt(mass(i))
  	enddo
         call u3best(mass,xyz,xyz0,Nprtn+NligNA,1,rms,U,T,ier)  
 	call MoveCoords(xyz,Nprtn+NligNA,U,T)

C================================================================
c Accumulate statistics about update...
C================================================================
	 dist = sqrt ((xyz(1,iatom)-xyz0(1,iatom))**2
     1	  	 +    (xyz(2,iatom)-xyz0(2,iatom))**2
     2	  	 +    (xyz(3,iatom)-xyz0(3,iatom))**2)

	  dist = sqrt ((xyz(1,iatom)-xyz0(1,iatom))**2
     1	  	  +    (xyz(2,iatom)-xyz0(2,iatom))**2
     2	  	  +    (xyz(3,iatom)-xyz0(3,iatom))**2)









        stop
        end




	subroutine MoveCoords(xyz,Nprtn,U,T)

	implicit none
	integer Nprtn
	real*8 xyz(3,Nprtn),T(3),U(3,3)

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


	return
	end




