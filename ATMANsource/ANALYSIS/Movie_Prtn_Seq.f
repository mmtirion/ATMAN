	program Movie_Seq

	implicit none
	character*80 file
	integer*4 Nprtn,Nres,Nxi,NPro,NDoF,mutate,Nlig
	integer NMAlig,NligNA,NDoFr,rscale
	real*8 Kelvin
	character*1 amidehydrogen

	parameter ( Nprtn  = 2598 ) !# of atoms in protein molecule
	parameter ( Nres   = 302 ) !# of residues
	parameter ( rscale = 0  ) ! Residue count-start
	parameter ( Nxi    = 478 ) !# of sidechain, Xi, torsions
	parameter ( NPro   = 14 ) !# of non-Nterminal-prolines (no phi torsions)

	parameter ( Nlig   = 0 ) !# of ligands 
	parameter ( NligNA = 0  ) !# of ligand-atoms total
	parameter ( NMAlig = 0  )!# of Multi-Atom ligands

	parameter ( mutate = 0  )    !# of MC dihedral blocks

	parameter ( NDoF   = 2*Nres+Nxi+6*NMAlig+(Nlig-NMAlig)*3 ) 
	parameter ( NDoFr  = NDoF-Npro-1-mutate) ! reduced # of dofs
	parameter ( Kelvin =  293.d0 ) !for setting amplitude of oscillation
	parameter (amidehydrogen='n')


	integer*4 IndexMC(2*Nres,4),block(100)
	integer*4 IndexXi(Nxi,4),EndXi(Nxi),XiRes(Nxi)
	real*8 xyzP(3,Nprtn),massP(Nprtn),Bhet(Nprtn),B(Nprtn)
	integer*4 rlist(Nprtn),rnum(Nres),rlst,Pyn(2*Nres+Nxi)
	character*4 rsort(Nprtn),rsrt
	character*1 type(Nprtn)
	character*5 asort(Nprtn),asrt
	character*3 tor
	real*8 lambda,omega,alpha,enhance,amp0,amp,kb,Na
	real*8 cutofflow,cutoffmed
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
	real array(Nprtn)
	integer natoms,num(Nprtn),ires
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
           Lig(1,1) =  1     ! # of atoms in ligand 1
           Lig(1,2) = 373  ! Residue # for ligand 1
           Lig(2,1) =  31    ! # of atoms in ligand 2
           Lig(2,2) = 375   ! Residue # for ligand 2
c 	  do i = 1,Nlig
c 	   Lig(i,1) = 1
c 	   Lig(i,2) = Nres+rscale+i
c 	  enddo
  	endif



c===================================================================
c  index torsions; identify atoms, masses, residues; get coordinates
c===================================================================
c	if(amidehydrogen.eq.'y')then
c	call index_prtnH(indexMC,indexXi,massP,xyzP,asort,B,rlist,
c     1	rsort,rnum,type,EndXi,XiRes,Nres,Nxi,Nprtn,Rhc,aHET,xyzL0,
c     2  rlstHET,mHET,vdWHet,typeHET,epsHET,NligNA,rscale,Bhet,P,Pyn)
c	else
	call index_prtn(indexMC,indexXi,massP,xyzP,asort,B,rlist,
     1	rsort,rnum,type,EndXi,XiRes,Nres,Nxi,Nprtn,Rhc,aHET,xyzL0,
     2  rlstHET,mHET,vdWHet,typeHET,epsHET,NligNA,rscale,Bhet,P,Pyn)
c	endif

c TEMPORARY for creating plot with atom index plus residues index on plot2
c	icount = 0
c	do i = 1, Nprtn
c	 ires = rlist(i)
c	 if(type(i).ne.'H')then
c	  icount = icount + 1
c	  array(icount) = -2.0
c	  num(ires)=icount
c	 endif
c	enddo
c	natoms = icount
c	do i = 10, Nres, 10
c	 array(num(i)) = 1.30
c	enddo
cc	do i = 1, natoms
c	 write(18,*)i,array(i)
c	enddo
c	do i = 1, natoms
c	 array(i) = -2.0
c	enddo
c	do i = 50, Nres,50
c	 array(num(i)) = 1.30
c	enddo
c	do i = 1, natoms
c	 write(19,*)i,array(i)
c	enddo
c END TEMPORARY


c insert overwrite of theoretical Bfactors over experimental values:
c 	open(23,file=
c     1'/Users/monique/Desktop/WORK/NMStuff/ATMAN/DATA/Bfactor')
c impose cutoff
c 	cutofflow = 4.0d0
c 	cutoffmed = 7.5d0
c 	do i = 1, Nprtn
c 	 read(23,*)k,B(i)
c 	 if(B(i).gt.cutoffmed)then
c 		B(i)=40.0d0
c		goto876
c 	 endif
c	 if(B(i).le.cutofflow)then
c		B(i)=0.0d0
c		goto876
c	 endif
c	 B(i) = 15.d0
c876	 continue
c 	enddo
c 	close(23)
c end insert overwrite B factors



c output PDB coordinate file for CartEV listing ahead:
	atom='ATOM  '
c	open(31,file='./RESULTS/CartEV.pdb')
c        do i = 1, Nprtn ! note NO LIGAND
c           write(31,20)atom,i,asort(i),rsort(i),rlist(i)+rscale,
c     1          xyzP(1,i),xyzP(2,i),xyzP(3,i),occ(i),B(i)
c	enddo
c	close(31)

	do i = 1,Nprtn
	 xyzP0(1,i)=xyzP(1,i)
	 xyzP0(2,i)=xyzP(2,i)
	 xyzP0(3,i)=xyzP(3,i)
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

c in case of mutating main chain degrees of freedom list:
	if(mutate.ne.0)then
         do j=1,10
          block(j)=78+j
         enddo
         do j=1,9
          block(10+j)=89+j
         enddo
         do j=1,19
          block(19+j)=399+j
         enddo
         do j = 1, mutate
	  indexMC(Block(j),3)=0	  
          write(*,*)' Block(',j,') = ',block(j)
         enddo

	endif

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
	 mass(iatom) = dsqrt(massP(iatom))
	enddo
        do iatom = 1,NligNA
	 mass(iatom+Nprtn) = mHET(iatom)
         Tmass = Tmass + mass(iatom+Nprtn)
	 mass(iatom+Nprtn) = dsqrt(mHET(iatom))
        enddo

C================================================================
C  Read in eigenfrequencies and eigenvectors
C================================================================
	write(*,*)'Which eigenvector?'
	read(*,*)Nmode
        Na = 6.02205d23  ! avogadros #
	kb = 1.3806d-23  ! boltzmann's constant

c  read in eigenvectors A(NDoFr):
c        open(81,file=
c     1'/Users/monique/Desktop/WORK/NMStuff/ATMAN/DATA/eigenvectors')
        open(81,file='../DATA/eigenvectors')

        DO i = 1,Nmode
	 read(81,*)m, omega
	 alpha = dsqrt(2.d0*kb*Kelvin)/dabs(omega)
	 write(*,*)'Mode: ',m,' temperature: ',Kelvin,
     1' Amplitude: ',alpha
	
	 max = 0.0d0
         do j = 1,NDoFr
          read(81,*)A(j)
	  if(j.le.(2*Nres-Npro-1-mutate).and.
     1		dabs(A(j)).gt.max)max=dabs(A(j))
         enddo
	ENDDO


c Use "Occupancy value" to create eigenvector info:
	do i = 1, Nprtn
	 occ(i) = 0.0d0  ! set each atom to 0.0
	enddo

c initially, see only if the mainchain dihedrals can be conveyed like this...
	iterm = 0
	do itor = 1,2*Nres
	 IF(Pyn(itor).ne.0)THEN
	 iterm = iterm + 1
	 if(mod(itor,2).eq.0)occ(indexMC(itor,3))=
     1		dabs(A(iterm))/max
	 if(mod(itor,2).eq.1)occ(indexMC(itor,2))=
     1		dabs(A(iterm))/max
	 ENDIF
	enddo
	do i = 1, Nprtn
	  if(occ(i).ge.0.5)then
	    tor = 'PHI'
	    if(asort(i).eq.' C   ')tor = 'PSI' 
	    write(*,*)tor,' ',rsort(i),rlist(i),occ(i)
	  endif
	enddo
c bin the values for ease of visualization:
	do i = 1, Nprtn
	 if(occ(i).ge.0.5)occ(i)=1.0
	 if(occ(i).ge.0.35.and.occ(i).ne.1.0)occ(i)=0.5
	 if(occ(i).lt.0.35)occ(i)=0.0
	enddo


C================================================================
C  Update according to chosen eigenvector:
C================================================================

	
	write(*,*)' '
	write(*,*)' enhance = ?'
	read(*,*)enhance

	amp0 = alpha*enhance
	write(*,*)'amp used:',amp0
	RMS_Ca = 0.0d0

	Nsequences = 9
	pi = dacos(-1.0d0)
	angle = 2.*pi / (dble(Nsequences-1))
	open(30,file='./RESULTS/animation.pdb')
c	open(77,file='./CartEV.txt') ! note this requires reformatted filetemp2.pdb listing!
c	write(77,*)'1E0W Cartesian eigenvector for mode ',Nmode


 	DO iSequence = 1,Nsequences
	amp = -amp0 * dsin(angle*dble(iSequence-1))


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
        call u3best(mass,xyz,xyz0,Nprtn+NligNA,1,rms,U,T,ier)  
 	call MoveCoords(xyz,Nprtn+NligNA,U,T)

C================================================================
c Accumulate statistics about update...
C================================================================
	write(*,*)'Statistics on sequence:',iSequence
	write(*,*)'------------------------------------------------'

	MaxDif = 0.0
	AveDif = 0.0

        do iatom = 1,Nprtn+NligNA
	 dist = sqrt ((xyz(1,iatom)-xyz0(1,iatom))**2
     1	  	 +    (xyz(2,iatom)-xyz0(2,iatom))**2
     2	  	 +    (xyz(3,iatom)-xyz0(3,iatom))**2)
	 AveDif = AveDif + dist*dist
	 if(dist.gt.MaxDif)MaxDif=dist
        enddo
	write(*,26)sqrt(AveDif/real(iatom)),MaxDif
26	format(' Average and Maximum RMSD for any atom:',2F8.3)

	MaxDif = 0.0
	RMS_Ca = 0.0
	iCa = 0

        do iatom = 1,Nprtn

c	 if(iSequence.eq.1)write(25,*)(xyz(i,iatom),i=1,3)
c	 if(iSequence.eq.3)write(26,*)(xyz(i,iatom),i=1,3)

	 if(asort(iatom).eq.' CA  ')then
	  dist = sqrt ((xyz(1,iatom)-xyz0(1,iatom))**2
     1	  	  +    (xyz(2,iatom)-xyz0(2,iatom))**2
     2	  	  +    (xyz(3,iatom)-xyz0(3,iatom))**2)
	  iCa = iCa + 1
	  if(dist.gt.MaxDif)Loc = iCa
	  if(dist.gt.MaxDif)MaxDif=dist
	  RMS_Ca = RMS_Ca + dist*dist
	 endif

        enddo
	RMS_Ca = sqrt(RMS_Ca/real(2*Nres))
	write(*,29)RMS_Ca,MaxDif,Loc
29	format(' Average and Maximum RMSD for any Calpha:',2F8.3,
     1' at residue number:',I4)
	write(*,*)' '

C================================================================
c output sequence:
C================================================================

	if(iSequence.eq.1)file='Mode1_Seq1.pdb'
	if(iSequence.eq.2)file='Mode1_Seq2.pdb'
	if(iSequence.eq.3)file='Mode1_Seq3.pdb'
	if(iSequence.eq.4)file='Mode1_Seq4.pdb'
	if(iSequence.eq.5)file='Mode1_Seq5.pdb'
	if(iSequence.eq.6)file='Mode1_Seq6.pdb'
	if(iSequence.eq.7)file='Mode1_Seq7.pdb'
	if(iSequence.eq.8)file='Mode1_Seq8.pdb'
	if(iSequence.eq.9)file='Mode1_Seq9.pdb'
	write(30,30)file
30	format('MODEL  Output iset=     0  "',A14,'"')

20	format(a6,x,i4,x,a5,a4,1x,i4,4x,3f8.3,2f6.2)

	do i = 1,Nprtn+NligNA
	 if(i.le.Nprtn)then
	   write(30,20)atom,i,asort(i),rsort(i),rlist(i)+rscale,
     1		xyz(1,i),xyz(2,i),xyz(3,i),occ(i),B(i)
	 else
	   iLig = i - Nprtn
	   atom='HETATM'
	   write(30,20)atom,i,aHET(iLig),typeHET(iLig),
     1		Nres+rscale+iLig,xyz(1,i),xyz(2,i),xyz(3,i)
	 endif
	enddo

	file = "ENDMDL"
	write(30,31)file
31	format(a6)


c OUTPUT approximate EigenVector in Cartesian Coordinates
c	if(iSequence.eq.3)then
c	 do i = 1, Nprtn
c	  write(77,77)(xyz(j,i)-xyz0(j,i),j=1,3)
c	 enddo
c	endif
c output eigenvector angle file for some xyz set:
c	if(iSequence.eq.3)then
c	open(33,file='RESULTS/filename.pdb')
c	 do i = 1, Nprtn
c	  write(77,77)(xyz(j,i)-xyz0(j,i),j=1,3)
c	 enddo
c	endif



C================================================================
c Place original coordinates back in array:
C================================================================
	do i = 1, Nprtn
	 xyzP(1,i) = xyz0(1,i)
	 xyzP(2,i) = xyz0(2,i)
	 xyzP(3,i) = xyz0(3,i)
	enddo
	 

	ENDDO

77	format('{ ',F13.10,', ',F13.10,', ',F13.10,'},')


	close(30)

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




