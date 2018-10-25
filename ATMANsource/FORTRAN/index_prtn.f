	subroutine index_prtn(indexMC,indexXi,mass,xyz,asort,B,rlist,
     1rsort,rlist1,type,EndXi,XiRes,Nres,Nxi,Nprtn,Rhc,aHET,rHET,
     2rlstHET,mHET,vdWHET,typeHET,epsHET,NligNA,rscale,Bhet,P,Pyn)

c Written by MMT 4/95 to read in a pdb file, order atoms, and index 
c mainchain and sidechain dihedral angles. Also assigns masses for atoms.
c Residues must be numbered consecutively: N, N+1, N+2,..,N+Nres. 

	implicit none
c INPUT:
	character*80 file,ftemp
	integer*4 ires,Nprtn,n,Npro
	real*8 x,y,z
c OUTPUT:
	integer Nxi,i,j,Nres,ixi,ss,Natm
	integer rlst,rscale,lig
	integer NA,NAi,NB,NBi,NDE,NDEi,ND
	integer NligNA,Pyn(2*Nres+Nxi)
	character*6 atom
	character*5 asrt,asort(Nprtn)
	character*4 rsrt,rsort(Nprtn),rlist2(Nres)
	character*1 type(Nprtn)
	integer IndexMC(2*Nres,4),rlist1(Nres),rlist(Nprtn)
	integer IndexPsi(Nres,4),IndexPhi(Nres,4),ResN(21)
	integer IndexXi(Nxi,4),EndXi(Nxi),XiRes(nxi)
	real*8 xyz(3,Nprtn),mass(Nprtn),Rhc(Nprtn)
	real*8 Bhet(2500),B(Nprtn),P(2*Nres+Nxi),Occ,Temp
c HETATM !I may need to make array sizes exact here...
	character*5 aHET(2500)
	character*4 typeHET(2500)
	integer rlstHET(2500)
	real*8 rHET(3,2500),mHET(2500),vdWHET(2500)
	real*8 epsHET(2500),massi,vdWi,epsi


c     	file=
c     1'/Users/monique/Desktop/WORK/NMStuff/PDB/2JDI/2jdi_Bmodh.pdb'
c TEMPORARY OVERWRITE
c	file=
c     1'./1e0w_sb01frame3reH.pdb'
c END TEMPORARY
 	file='../PDB/1GOK_Ah.pdb'

	ftemp='/Users/monique/Desktop/WORK/NMStuff/ATMAN/filetemp2.pdb'
	open(30,file=file)
	open(40,file=ftemp)
	write(*,*)' '
	write(*,*)'Reading in data from: ',file
	write(*,*)'... '
10	format(a6,x,i4,x,a5,a4,1x,i4,3f22.17)  !for MINIMIZATION!!
12	format(a6,x,i4,x,a5,a4,1x,i4,4x,3f8.3,2x,f4.2,x,f5.2)

c==================================================================
c Read in coordinates, and reorder atom arrangement within residues
c==================================================================
c read in all except mc amide hydrogen atoms and ASN,GLN,HIS TRP amide hydrogens

	write(*,*)'Residue offset set to ',rscale
c rscale permits residue numbering to start at value other than 1.

	Natm = 0	
	lig = 0
	DO i = 1, 20000    !maximum atomic entries limited here
	 Temp =0.d0
	 Occ = 1.d0
	 read(30,12,end=900)atom,n,asrt,rsrt,rlst,x,y,z,Occ,Temp
c TEMPORARY
c	 write(*,12)atom,n,asrt,rsrt,rlst,x,y,z,Occ,Temp
c END TEMPORARY
c	 read(30,10,end=900)atom,n,asrt,rsrt,rlst,x,y,z  ! FOR MINIMIZATION
	 If(atom.eq.'ATOM  '.and.asrt(1:1).EQ.' ')Then
 	  if(asrt(2:2).eq.'H'.and.asrt(3:3).ne.' ')goto899
		Natm = Natm+1	
	 	asort(Natm) = asrt
	 	rsort(Natm) = rsrt
	 	rlist(Natm) = rlst-rscale
	 	xyz(1,Natm) = x
	 	xyz(2,Natm) = y
	 	xyz(3,Natm) = z
		B(Natm) = Temp  
	 	rlist1(rlist(Natm)) = Natm  !id's last entry w/in 1-d formula of res i
		rlist2(rlist(Natm)) = rsort(Natm)
899       continue	 
	 EndIf
	 If(atom.eq.'HETATM')Then !set up for general 3RB ligands on 8/14 and 6RB 2/15
	  lig = lig + 1
c	write(*,*)lig
	  if(lig.gt.2500)then
	   write(*,*)'Error'
	   write(*,*)'# of Ligand atoms > 2500 '
	   stop
	  endif
	  aHET(lig) = asrt
	  rlstHET(lig) = rlst - rscale 
	  rHET(1,lig) = x
	  rHET(2,lig) = y
	  rHET(3,lig) = z
	  typeHET(lig)=rsrt
	  call HETmass(asrt,rsrt,massi,vdWi,epsi)
	  mHET(lig) = massi
	  vdWHET(lig) = vdWi
	  epsHET(lig) = epsi	
	  Bhet(lig) = Temp  
	 EndIf
891	Continue
	ENDDO
900	close(30)
	write(*,*)' '
	write(*,*)'# of (unculled) ligand atoms:',lig


c count no. of PRO
	Npro = 0
	do i = 1, Nres
	 if(rlist2(i).eq.'PRO ')Npro=Npro + 1
	enddo
	write(*,*)'# of prolines in polymer: ',Npro
	write(*,*)'# of RES in protein:',rlist(Natm)  ! problem w/ HETATM!
	if(rlist2(1).eq.'PRO ')Npro = Npro - 1

	write(*,*)'# of atoms in SEQUENCE',Natm
	write(*,*)'	vs. Nprtn ',Nprtn
	if(Natm.ne.Nprtn)stop
	write(*,*)'# torsional DoF in protein:',2*Nres+Nxi
	write(*,*)' '

c count degrees of freedom
	NA = 0
	NB = 0
	NDE = 0
	do i = 1, 21
	 ResN(i) = 0
	enddo
	DO i = 1, Nres
	 call CountDoF(rlist2(i),NBi,NAi,NDEi,ResN)
	 NA = NA + NAi
	 NB = NB + NBi
	 NDE = NDE + NDEi
	ENDDO
 	NB = NB + (Nres-1)  ! Peptide Bonds
	NA = NA + 3*(Nres-1)
	NDE = NDE + Nres + 1  ! for MC carbonyl oxygens
	ND = NB - NDE 
	Natm = Nprtn - (Nres-Npro-1)

	write(*,*)' '
	write(*,*)'Degrees of Freedom:'
	write(*,*)'# of Cartesian of',Natm,'heavy atoms:'
	write(*,*)'         ',3*Natm
	write(*,*)' '
	write(*,*)'# of Bonds:', NB
	write(*,*)'# of Angles:',NA
	write(*,*)'# of Dihedrals:',ND
	write(*,*)'Total internal DoF:',NB+NA+ND+6
	write(*,*)' '
	 write(*,*)'ALA',ResN(1)
	 write(*,*)'GLY',ResN(2)
	 write(*,*)'ILE',ResN(3)
	 write(*,*)'LEU',ResN(4)
	 write(*,*)'PRO',ResN(5)
	 write(*,*)'PHE',ResN(6)
	 write(*,*)'VAL',ResN(7)
	 write(*,*)'ARG',ResN(8)
	 write(*,*)'ASP',ResN(9)
	 write(*,*)'GLU',ResN(10)
	 write(*,*)'SER',ResN(11)
	 write(*,*)'THR',ResN(12)
	 write(*,*)'CYS',ResN(13)
	 write(*,*)'ASN',ResN(14)
	 write(*,*)'GLN',ResN(15)
	 write(*,*)'HIS',ResN(16)
	 write(*,*)'LYS',ResN(17)
	 write(*,*)'TYR',ResN(18)
	 write(*,*)'MET',ResN(19)
	 write(*,*)'TRP',ResN(20)
	 write(*,*)'PCA',ResN(21)
c check if there's a terminal oxygen, OXT, and if so, place in last position:
	do i = 1, Nprtn
	 if(asort(i).eq.' OXT ')then
	  if(i.eq.Nprtn)goto898
	   asrt=asort(Nprtn)
	   rsrt=rsort(Nprtn)
	   rlst=rlist(Nprtn)
	   x=xyz(1,Nprtn)
	   y=xyz(2,Nprtn)
	   z=xyz(3,Nprtn)
	   Temp = B(Nprtn)
	   asort(Nprtn)=asort(i)
	   rsort(Nprtn)=rsort(i)
	   rlist(Nprtn)=rlist(i)
	   xyz(1,Nprtn)=xyz(1,i)
	   xyz(2,Nprtn)=xyz(2,i)
	   xyz(3,Nprtn)=xyz(3,i)
	   B(Nprtn) = B(i)
	   asort(i)=asrt
	   rsort(i)=rsrt
	   rlist(i)=rlst
	   xyz(1,i)=x
	   xyz(2,i)=y
	   xyz(3,i)=z
	   B(i) = Temp
	   rlist1(Nres)=Nprtn-1   ! this MAY have unforeseen consequences...8/'14
	   rlist2(Nres)=rsort(Nprtn)
	   goto898
	 endif
	enddo
898	continue

c now rearrange so that atomic-residue arrangement of proteins atoms is N-Ca-Cb-...-C-O  :
	write(*,*)'...now rearrange the atomic sequence of protein...'
       call rearrange(asort,rsort,rlist,xyz,rlist1,rlist2,
     1  B,Nres,Nprtn)
	atom='ATOM  '  !stop gap till i get HETATM corrected!! 8/14
	do i = 1, Nprtn
 	 write(40,12)atom,i,asort(i),rsort(i),rlist(i)+rscale,
     1	(xyz(j,i),j=1,3)
c 	 write(40,12)atom,i,asort(i),rsort(i),rlist(i)+rscale,
c     1	(xyz(j,i),j=1,3),Occ,B(i)
	enddo
	if(lig.gt.0)then
	atom='HETATM'
	do i = 1, lig
 	 write(40,12)atom,Nprtn+i,aHET(i),typeHET(i),Nres+i+rscale,
     1		(rHET(j,i),j=1,3)
	enddo
	endif
	close(40)

	write(*,*)'	...rearrangement complete.'
c==================================================================
c index MAIN-CHAIN torsions (Phi/Psi-dihedrals)
c==================================================================
c initialize:
	do i=1,nres
	 IndexPhi(i,1)=0
	 IndexPhi(i,2)=0
	 IndexPhi(i,3)=0
	 IndexPhi(i,4)=0
	 IndexPsi(i,1)=0
	 IndexPsi(i,2)=0
	 IndexPsi(i,3)=0
	 IndexPsi(i,4)=0
	enddo

	do i=1,Nprtn


c	write(*,*)i,asort(i),rsort(i)
	 call sort(asort(i),rsort(i),type(i))
	 call amu(type(i),mass(i),Rhc(i))


	 if(asort(i).eq.' N   ')then
		IndexPsi(rlist(i),1) = i
		if(rlist(i).ge.2)IndexPsi(rlist(i)-1,4) = i
		IndexPhi(rlist(i),2) = i
	 endif

	 if(asort(i).eq.' CA  ')then
		IndexPsi(rlist(i),2) = i
		IndexPhi(rlist(i),3) = i
	 endif
 
	 if(asort(i).eq.' C   ')then
		IndexPsi(rlist(i),3) = i
		IndexPhi(rlist(i),4) = i
		if(rlist(i).le.nres-1)IndexPhi(rlist(i)+1,1) = i
	 endif

c block (mainchain) Phi-torsions of prolines:
c right now (6/15) there's a problem if PCA is not a N terminal residue,
c since i've not blocked PCA's MC phi angles: it'll never have an H, for example!

c	if(rsort(i).eq.'PRO ')then
c		IndexPhi(rlist(i),1) = 0
c		IndexPhi(rlist(i),2) = 0
c		IndexPhi(rlist(i),3) = 0
c		IndexPhi(rlist(i),4) = 0
c	 endif

	enddo

c fix N- and C-terminal boundaries:
c block 1st Phi torsion:
	IndexPhi(1,1) = 0
	IndexPhi(1,2) = 0
	IndexPhi(1,3) = 0
	IndexPhi(1,4) = 0
c substitute OT for 4th-atom-index of final Psi torsion:
	IndexPsi(nres,4) = Nprtn

c rearrange into single array:
	do ires = 1,Nres
	 IndexMC(2*ires-1,1) = IndexPhi(ires,1)
	 IndexMC(2*ires-1,2) = IndexPhi(ires,2)
	 IndexMC(2*ires-1,3) = IndexPhi(ires,3)
	 IndexMC(2*ires-1,4) = IndexPhi(ires,4)
	 IndexMC(2*ires,1) = IndexPsi(ires,1)
	 IndexMC(2*ires,2) = IndexPsi(ires,2)
	 IndexMC(2*ires,3) = IndexPsi(ires,3)
	 IndexMC(2*ires,4) = IndexPsi(ires,4)
	enddo


C==================================================================
c Sections to index SIDE-CHAIN torsions (Xi-dihedrals)
c==================================================================
c initialize:
	do i=1,nxi
	 indexXi(i,1)=0
	 indexXi(i,2)=0
	 indexXi(i,3)=0
	 indexXi(i,4)=0
	enddo

	ixi=0
	ss=1  !ss stands for 'Start-Search'
	DO ires=1,Nres
	 if(rlist2(ires).eq.'ALA '.or.rlist2(ires).eq.'GLY '.or.
     1      rlist2(ires).eq.'PRO '.or.rlist2(ires).eq.'PCA ')
     2	    goto800
	 ixi=ixi+1
	 XiRes(ixi)=ires
	 if(ires.ne.1)ss=rlist1(ires-1)+1

	 if(rlist2(ires).eq.'VAL ')then
101	  if(asort(ss).eq.' N   ')IndexXi(ixi,1)=ss
	  if(asort(ss).eq.' CA  ')IndexXi(ixi,2)=ss
 	  if(asort(ss).eq.' CB  ')IndexXi(ixi,3)=ss
	  if(asort(ss).eq.' CG1 ')IndexXi(ixi,4)=ss
	  if(asort(ss).eq.' CG2 ')EndXi(ixi)=ss
	  ss=ss+1
	  if(ss.le.rlist1(ires))goto101
	 endif

	 if(rlist2(ires).eq.'CYS ')then
102	  if(asort(ss).eq.' N   ')IndexXi(ixi,1)=ss
	  if(asort(ss).eq.' CA  ')IndexXi(ixi,2)=ss
 	  if(asort(ss).eq.' CB  ')IndexXi(ixi,3)=ss
	  if(asort(ss).eq.' SG  ')IndexXi(ixi,4)=ss
	  if(asort(ss).eq.' SG  ')EndXi(ixi)=ss
	  ss=ss+1
	  if(ss.le.rlist1(ires))goto102
	 endif

	 if(rlist2(ires).eq.'SER ')then
103	  if(asort(ss).eq.' N   ')IndexXi(ixi,1)=ss
	  if(asort(ss).eq.' CA  ')IndexXi(ixi,2)=ss
 	  if(asort(ss).eq.' CB  ')IndexXi(ixi,3)=ss
	  if(asort(ss).eq.' OG  ')IndexXi(ixi,4)=ss
	  if(asort(ss).eq.' OG  ')EndXi(ixi)=ss
	  ss=ss+1
	  if(ss.le.rlist1(ires))goto103
	 endif

	 if(rlist2(ires).eq.'THR ')then
104	  if(asort(ss).eq.' N   ')IndexXi(ixi,1)=ss
	  if(asort(ss).eq.' CA  ')IndexXi(ixi,2)=ss
 	  if(asort(ss).eq.' CB  ')IndexXi(ixi,3)=ss
	  if(asort(ss).eq.' OG1 ')IndexXi(ixi,4)=ss
	  if(asort(ss).eq.' CG2 ')EndXi(ixi)=ss
	  ss=ss+1
	  if(ss.le.rlist1(ires))goto104
	 endif

	 if(rlist2(ires).eq.'LEU ')then
105	  if(asort(ss).eq.' N   ')IndexXi(ixi,1)=ss
	  if(asort(ss).eq.' CA  ')IndexXi(ixi,2)=ss
 	  if(asort(ss).eq.' CB  ')IndexXi(ixi,3)=ss
	  if(asort(ss).eq.' CG  ')IndexXi(ixi,4)=ss
	  if(asort(ss).eq.' CD2 ')EndXi(ixi)=ss

	  if(asort(ss).eq.' CA  ')IndexXi(ixi+1,1)=ss
	  if(asort(ss).eq.' CB  ')IndexXi(ixi+1,2)=ss
 	  if(asort(ss).eq.' CG  ')IndexXi(ixi+1,3)=ss
	  if(asort(ss).eq.' CD1 ')IndexXi(ixi+1,4)=ss
	  if(asort(ss).eq.' CD2 ')EndXi(ixi+1)=ss
	  ss=ss+1
	  if(ss.le.rlist1(ires))goto105
	  ixi=ixi+1
	  XiRes(ixi)=ires
	 endif

	 if(rlist2(ires).eq.'ILE ')then
106	  if(asort(ss).eq.' N   ')IndexXi(ixi,1)=ss
	  if(asort(ss).eq.' CA  ')IndexXi(ixi,2)=ss
 	  if(asort(ss).eq.' CB  ')IndexXi(ixi,3)=ss
	  if(asort(ss).eq.' CG1 ')IndexXi(ixi,4)=ss
	  if(asort(ss).eq.' CD1 ')EndXi(ixi)=ss
	  if(asort(ss).eq.' CD1 ')EndXi(ixi)=ss

	  if(asort(ss).eq.' CA  ')IndexXi(ixi+1,1)=ss
	  if(asort(ss).eq.' CB  ')IndexXi(ixi+1,2)=ss
 	  if(asort(ss).eq.' CG1 ')IndexXi(ixi+1,3)=ss
	  if(asort(ss).eq.' CD1 ')IndexXi(ixi+1,4)=ss
	  if(asort(ss).eq.' CD1 ')EndXi(ixi+1)=ss
	  ss=ss+1
	  if(ss.le.rlist1(ires))goto106
	  ixi=ixi+1
	  XiRes(ixi)=ires
	 endif

	 if(rlist2(ires).eq.'PHE ')then
107	  if(asort(ss).eq.' N   ')IndexXi(ixi,1)=ss
	  if(asort(ss).eq.' CA  ')IndexXi(ixi,2)=ss
 	  if(asort(ss).eq.' CB  ')IndexXi(ixi,3)=ss
	  if(asort(ss).eq.' CG  ')IndexXi(ixi,4)=ss
	  if(asort(ss).eq.' CZ  ')EndXi(ixi)=ss

	  if(asort(ss).eq.' CA  ')IndexXi(ixi+1,1)=ss
	  if(asort(ss).eq.' CB  ')IndexXi(ixi+1,2)=ss
 	  if(asort(ss).eq.' CG  ')IndexXi(ixi+1,3)=ss
	  if(asort(ss).eq.' CD1 ')IndexXi(ixi+1,4)=ss
	  if(asort(ss).eq.' CZ  ')EndXi(ixi+1)=ss
	  ss=ss+1
	  if(ss.le.rlist1(ires))goto107
	  ixi=ixi+1
	  XiRes(ixi)=ires
	 endif

	 if(rlist2(ires).eq.'TYR ')then
108	  if(asort(ss).eq.' N   ')IndexXi(ixi,1)=ss
	  if(asort(ss).eq.' CA  ')IndexXi(ixi,2)=ss
 	  if(asort(ss).eq.' CB  ')IndexXi(ixi,3)=ss
	  if(asort(ss).eq.' CG  ')IndexXi(ixi,4)=ss
	  if(asort(ss).eq.' OH  ')EndXi(ixi)=ss

	  if(asort(ss).eq.' CA  ')IndexXi(ixi+1,1)=ss
	  if(asort(ss).eq.' CB  ')IndexXi(ixi+1,2)=ss
 	  if(asort(ss).eq.' CG  ')IndexXi(ixi+1,3)=ss
	  if(asort(ss).eq.' CD1 ')IndexXi(ixi+1,4)=ss
	  if(asort(ss).eq.' OH  ')EndXi(ixi+1)=ss
	  ss=ss+1
	  if(ss.le.rlist1(ires))goto108
	  ixi=ixi+1
	  XiRes(ixi)=ires
	 endif

	 if(rlist2(ires).eq.'HIS ')then
109	  if(asort(ss).eq.' N   ')IndexXi(ixi,1)=ss
	  if(asort(ss).eq.' CA  ')IndexXi(ixi,2)=ss
 	  if(asort(ss).eq.' CB  ')IndexXi(ixi,3)=ss
	  if(asort(ss).eq.' CG  ')IndexXi(ixi,4)=ss
	  if(asort(ss).eq.' NE2 ')EndXi(ixi)=ss

	  if(asort(ss).eq.' CA  ')IndexXi(ixi+1,1)=ss
	  if(asort(ss).eq.' CB  ')IndexXi(ixi+1,2)=ss
 	  if(asort(ss).eq.' CG  ')IndexXi(ixi+1,3)=ss
	  if(asort(ss).eq.' ND1 ')IndexXi(ixi+1,4)=ss
	  if(asort(ss).eq.' NE2 ')EndXi(ixi+1)=ss
	  ss=ss+1
	  if(ss.le.rlist1(ires))goto109
	  ixi=ixi+1
	  XiRes(ixi)=ires
	 endif

	 if(rlist2(ires).eq.'HMS ')then
110	  if(asort(ss).eq.' N   ')IndexXi(ixi,1)=ss
	  if(asort(ss).eq.' CA  ')IndexXi(ixi,2)=ss
 	  if(asort(ss).eq.' CB  ')IndexXi(ixi,3)=ss
	  if(asort(ss).eq.' CG  ')IndexXi(ixi,4)=ss
	  if(asort(ss).eq.' CM  ')EndXi(ixi)=ss

	  if(asort(ss).eq.' CA  ')IndexXi(ixi+1,1)=ss
	  if(asort(ss).eq.' CB  ')IndexXi(ixi+1,2)=ss
 	  if(asort(ss).eq.' CG  ')IndexXi(ixi+1,3)=ss
	  if(asort(ss).eq.' ND1 ')IndexXi(ixi+1,4)=ss
	  if(asort(ss).eq.' CM  ')EndXi(ixi+1)=ss
	  ss=ss+1
	  if(ss.le.rlist1(ires))goto110
	  ixi=ixi+1
	  XiRes(ixi)=ires
	 endif

	 if(rlist2(ires).eq.'TRP ')then
111	  if(asort(ss).eq.' N   ')IndexXi(ixi,1)=ss
	  if(asort(ss).eq.' CA  ')IndexXi(ixi,2)=ss
 	  if(asort(ss).eq.' CB  ')IndexXi(ixi,3)=ss
	  if(asort(ss).eq.' CG  ')IndexXi(ixi,4)=ss
	  if(asort(ss).eq.' CH2 ')EndXi(ixi)=ss

	  if(asort(ss).eq.' CA  ')IndexXi(ixi+1,1)=ss
	  if(asort(ss).eq.' CB  ')IndexXi(ixi+1,2)=ss
 	  if(asort(ss).eq.' CG  ')IndexXi(ixi+1,3)=ss
	  if(asort(ss).eq.' CD1 ')IndexXi(ixi+1,4)=ss
	  if(asort(ss).eq.' CH2 ')EndXi(ixi+1)=ss
	  ss=ss+1
	  if(ss.le.rlist1(ires))goto111
	  ixi=ixi+1
	  XiRes(ixi)=ires
	 endif

	 if(rlist2(ires).eq.'ASP ')then
112	  if(asort(ss).eq.' N   ')IndexXi(ixi,1)=ss
	  if(asort(ss).eq.' CA  ')IndexXi(ixi,2)=ss
 	  if(asort(ss).eq.' CB  ')IndexXi(ixi,3)=ss
	  if(asort(ss).eq.' CG  ')IndexXi(ixi,4)=ss
	  if(asort(ss).eq.' OD2 ')EndXi(ixi)=ss

	  if(asort(ss).eq.' CA  ')IndexXi(ixi+1,1)=ss
	  if(asort(ss).eq.' CB  ')IndexXi(ixi+1,2)=ss
 	  if(asort(ss).eq.' CG  ')IndexXi(ixi+1,3)=ss
	  if(asort(ss).eq.' OD1 ')IndexXi(ixi+1,4)=ss
	  if(asort(ss).eq.' OD2 ')EndXi(ixi+1)=ss
	  ss=ss+1
	  if(ss.le.rlist1(ires))goto112
	  ixi=ixi+1
	  XiRes(ixi)=ires
	 endif

	 if(rlist2(ires).eq.'ASN ')then
118	  if(asort(ss).eq.' N   ')IndexXi(ixi,1)=ss
	  if(asort(ss).eq.' CA  ')IndexXi(ixi,2)=ss
 	  if(asort(ss).eq.' CB  ')IndexXi(ixi,3)=ss
	  if(asort(ss).eq.' CG  ')IndexXi(ixi,4)=ss
	  if(asort(ss).eq.' ND2 '.or.asort(ss).eq.'AD2 ')
     1		EndXi(ixi)=ss

	  if(asort(ss).eq.' CA  ')IndexXi(ixi+1,1)=ss
	  if(asort(ss).eq.' CB  ')IndexXi(ixi+1,2)=ss
 	  if(asort(ss).eq.' CG  ')IndexXi(ixi+1,3)=ss
	  if(asort(ss).eq.' OD1 '.or.asort(ss).eq.'AD1 ')
     1		IndexXi(ixi+1,4)=ss
	  if(asort(ss).eq.' ND2 '.or.asort(ss).eq.' AD2 ')
     1		EndXi(ixi+1)=ss
	  ss=ss+1
	  if(ss.le.rlist1(ires))goto118
	  ixi=ixi+1
	  XiRes(ixi)=ires
	 endif

	 if(rlist2(ires).eq.'GLU ')then
113	  if(asort(ss).eq.' N   ')IndexXi(ixi,1)=ss
	  if(asort(ss).eq.' CA  ')IndexXi(ixi,2)=ss
 	  if(asort(ss).eq.' CB  ')IndexXi(ixi,3)=ss
	  if(asort(ss).eq.' CG  ')IndexXi(ixi,4)=ss
	  if(asort(ss).eq.' OE2 ')EndXi(ixi)=ss

	  if(asort(ss).eq.' CA  ')IndexXi(ixi+1,1)=ss
	  if(asort(ss).eq.' CB  ')IndexXi(ixi+1,2)=ss
 	  if(asort(ss).eq.' CG  ')IndexXi(ixi+1,3)=ss
	  if(asort(ss).eq.' CD  ')IndexXi(ixi+1,4)=ss
	  if(asort(ss).eq.' OE2 ')EndXi(ixi+1)=ss

	  if(asort(ss).eq.' CB  ')IndexXi(ixi+2,1)=ss
	  if(asort(ss).eq.' CG  ')IndexXi(ixi+2,2)=ss
 	  if(asort(ss).eq.' CD  ')IndexXi(ixi+2,3)=ss
	  if(asort(ss).eq.' OE1 ')IndexXi(ixi+2,4)=ss
	  if(asort(ss).eq.' OE2 ')EndXi(ixi+2)=ss
	  ss=ss+1
	  if(ss.le.rlist1(ires))goto113
	  XiRes(ixi+1)=ires
	  XiRes(ixi+2)=ires
	  ixi=ixi+2
	 endif

	 if(rlist2(ires).eq.'GLN ')then
114	  if(asort(ss).eq.' N   ')IndexXi(ixi,1)=ss
	  if(asort(ss).eq.' CA  ')IndexXi(ixi,2)=ss
 	  if(asort(ss).eq.' CB  ')IndexXi(ixi,3)=ss
	  if(asort(ss).eq.' CG  ')IndexXi(ixi,4)=ss
	  if(asort(ss).eq.' NE2 '.or.asort(ss).eq.'AE2 ')
     1		EndXi(ixi)=ss

	  if(asort(ss).eq.' CA  ')IndexXi(ixi+1,1)=ss
	  if(asort(ss).eq.' CB  ')IndexXi(ixi+1,2)=ss
 	  if(asort(ss).eq.' CG  ')IndexXi(ixi+1,3)=ss
	  if(asort(ss).eq.' CD  ')IndexXi(ixi+1,4)=ss
	  if(asort(ss).eq.' NE2 '.or.asort(ss).eq.'AE2 ')
     1		EndXi(ixi+1)=ss

	  if(asort(ss).eq.' CB  ')IndexXi(ixi+2,1)=ss
	  if(asort(ss).eq.' CG  ')IndexXi(ixi+2,2)=ss
 	  if(asort(ss).eq.' CD  ')IndexXi(ixi+2,3)=ss
	  if(asort(ss).eq.' OE1 '.or.asort(ss).eq.'AE1 ')
     1		IndexXi(ixi+2,4)=ss
	  if(asort(ss).eq.' NE2 '.or.asort(ss).eq.' AE2 ')
     1		EndXi(ixi+2)=ss
	  ss=ss+1
	  if(ss.le.rlist1(ires))goto114
	  XiRes(ixi+1)=ires
	  XiRes(ixi+2)=ires
	  ixi=ixi+2
	 endif


	 if(rlist2(ires).eq.'MET ')then
115	  if(asort(ss).eq.' N   ')IndexXi(ixi,1)=ss
	  if(asort(ss).eq.' CA  ')IndexXi(ixi,2)=ss
 	  if(asort(ss).eq.' CB  ')IndexXi(ixi,3)=ss
	  if(asort(ss).eq.' CG  ')IndexXi(ixi,4)=ss
	  if(asort(ss).eq.' CE  ')EndXi(ixi)=ss

	  if(asort(ss).eq.' CA  ')IndexXi(ixi+1,1)=ss
	  if(asort(ss).eq.' CB  ')IndexXi(ixi+1,2)=ss
 	  if(asort(ss).eq.' CG  ')IndexXi(ixi+1,3)=ss
	  if(asort(ss).eq.' SD  ')IndexXi(ixi+1,4)=ss
	  if(asort(ss).eq.' CE  ')EndXi(ixi+1)=ss

	  if(asort(ss).eq.' CB  ')IndexXi(ixi+2,1)=ss
	  if(asort(ss).eq.' CG  ')IndexXi(ixi+2,2)=ss
 	  if(asort(ss).eq.' SD  ')IndexXi(ixi+2,3)=ss
	  if(asort(ss).eq.' CE  ')IndexXi(ixi+2,4)=ss
	  if(asort(ss).eq.' CE  ')EndXi(ixi+2)=ss
	  ss=ss+1
	  if(ss.le.rlist1(ires))goto115
	  XiRes(ixi+1)=ires
	  XiRes(ixi+2)=ires
	  ixi=ixi+2
	 endif

	 if(rlist2(ires).eq.'ARG ')then
116	  if(asort(ss).eq.' N   ')IndexXi(ixi,1)=ss
	  if(asort(ss).eq.' CA  ')IndexXi(ixi,2)=ss
 	  if(asort(ss).eq.' CB  ')IndexXi(ixi,3)=ss
	  if(asort(ss).eq.' CG  ')IndexXi(ixi,4)=ss
	  if(asort(ss).eq.' NH2 ')EndXi(ixi)=ss

	  if(asort(ss).eq.' CA  ')IndexXi(ixi+1,1)=ss
	  if(asort(ss).eq.' CB  ')IndexXi(ixi+1,2)=ss
 	  if(asort(ss).eq.' CG  ')IndexXi(ixi+1,3)=ss
	  if(asort(ss).eq.' CD  ')IndexXi(ixi+1,4)=ss
	  if(asort(ss).eq.' NH2 ')EndXi(ixi+1)=ss

	  if(asort(ss).eq.' CB  ')IndexXi(ixi+2,1)=ss
	  if(asort(ss).eq.' CG  ')IndexXi(ixi+2,2)=ss
 	  if(asort(ss).eq.' CD  ')IndexXi(ixi+2,3)=ss
	  if(asort(ss).eq.' NE  ')IndexXi(ixi+2,4)=ss
	  if(asort(ss).eq.' NH2 ')EndXi(ixi+2)=ss

	  if(asort(ss).eq.' CG  ')IndexXi(ixi+3,1)=ss
	  if(asort(ss).eq.' CD  ')IndexXi(ixi+3,2)=ss
 	  if(asort(ss).eq.' NE  ')IndexXi(ixi+3,3)=ss
	  if(asort(ss).eq.' CZ  ')IndexXi(ixi+3,4)=ss
	  if(asort(ss).eq.' NH2 ')EndXi(ixi+3)=ss

	  if(asort(ss).eq.' CD  ')IndexXi(ixi+4,1)=ss
	  if(asort(ss).eq.' NE  ')IndexXi(ixi+4,2)=ss
 	  if(asort(ss).eq.' CZ  ')IndexXi(ixi+4,3)=ss
	  if(asort(ss).eq.' NH1 ')IndexXi(ixi+4,4)=ss
	  if(asort(ss).eq.' NH2 ')EndXi(ixi+4)=ss
	  ss=ss+1
	  if(ss.le.rlist1(ires))goto116
	  XiRes(ixi+1)=ires
	  XiRes(ixi+2)=ires
	  XiRes(ixi+3)=ires
	  ixi=ixi+4
	 endif

	 if(rlist2(ires).eq.'LYS ')then
117	  if(asort(ss).eq.' N   ')IndexXi(ixi,1)=ss
	  if(asort(ss).eq.' CA  ')IndexXi(ixi,2)=ss
 	  if(asort(ss).eq.' CB  ')IndexXi(ixi,3)=ss
	  if(asort(ss).eq.' CG  ')IndexXi(ixi,4)=ss
	  if(asort(ss).eq.' NZ  ')EndXi(ixi)=ss

	  if(asort(ss).eq.' CA  ')IndexXi(ixi+1,1)=ss
	  if(asort(ss).eq.' CB  ')IndexXi(ixi+1,2)=ss
 	  if(asort(ss).eq.' CG  ')IndexXi(ixi+1,3)=ss
	  if(asort(ss).eq.' CD  ')IndexXi(ixi+1,4)=ss
	  if(asort(ss).eq.' NZ  ')EndXi(ixi+1)=ss

	  if(asort(ss).eq.' CB  ')IndexXi(ixi+2,1)=ss
	  if(asort(ss).eq.' CG  ')IndexXi(ixi+2,2)=ss
 	  if(asort(ss).eq.' CD  ')IndexXi(ixi+2,3)=ss
	  if(asort(ss).eq.' CE  ')IndexXi(ixi+2,4)=ss
	  if(asort(ss).eq.' NZ  ')EndXi(ixi+2)=ss

	  if(asort(ss).eq.' CG  ')IndexXi(ixi+3,1)=ss
	  if(asort(ss).eq.' CD  ')IndexXi(ixi+3,2)=ss
 	  if(asort(ss).eq.' CE  ')IndexXi(ixi+3,3)=ss
	  if(asort(ss).eq.' NZ  ')IndexXi(ixi+3,4)=ss
	  if(asort(ss).eq.' NZ  ')EndXi(ixi+3)=ss
	  ss=ss+1
	  if(ss.le.rlist1(ires))goto117
	  XiRes(ixi+1)=ires
	  XiRes(ixi+2)=ires
	  XiRes(ixi+3)=ires
	  ixi=ixi+3
	 endif

800	continue
	

	ENDDO

	write(*,*)' '
	write(*,*)' # of Xi torsions:',ixi
	write(*,*)' '

	do i = 1, ixi
	 IF(IndexXi(i,1).eq.0.or.indexXi(i,2).eq.0.or.
     1	 indexXi(i,3).eq.0.or.indexXi(i,4).eq.0)then
		write(*,*)' WARNING'
		write(*,*)
		write(*,*)' Residue: ',XiRes(i),rlist2(XiRes(i))
		write(*,*)' Stopping program...maybe more un-IDed atoms'
		stop
	 ENDIF
	enddo



	call Pvec(Nprtn,Nres,Nxi,IndexMC,IndexXi,xyz,P)

C TEMPORARY OUTPUT DIHEDRAL ANGLES
c13	format(i5,x,f6.2)
c        do i = 1, 2*Nres
c         write(61,13)i,P(i)
c        enddo
c        do i = 1, Nxi
c        write(63,13)i,P(2*Nres+i)
c        enddo
c        close(61)
c        close(63)
c END TEMPORARY




	Pyn(1) = 0  ! N terminus has no degree of freedom
	Pyn(2) = 1  
	do i = 2,Nres
	 Pyn(2*i) = 1 !all psi included
	 Pyn(2*i-1) = 1  ! exclude prolines phi
	 if(rlist2(i).eq.'PRO ')Pyn(2*i-1)=0
	enddo
	do i = 1, Nxi  !all xi included
	 Pyn(2*Nres+i) = 1
	enddo 







	return
	end





	subroutine rearrange(asort,rsort,rlist,xyz,rlist1,rlist2,
     1	B,Nres,Nprtn)

	integer Nres,Nprtn
	character*4 rsort(Nprtn),rlist2(Nres)
	character*5 asort(Nprtn)
	integer rlist1(Nres),rlist(Nprtn)
	real*8 xyz(3,Nprtn),B(Nprtn)

	integer i,j,imax,ires,ranksize,rank(15)
	character*4 resnam,temprsort(15)
	character*5 tempasort(15)
	real*8 tempxyz(3,15),tempB(15)
	

	DO ires=1,Nres

c initialize
	do j = 1, 15
	 rank(j)=-1
	enddo
	resnam=rlist2(ires)

	if(resnam.eq.'GLY ')ranksize=5
	if(resnam.eq.'ALA ')ranksize=6
	if(resnam.eq.'SER '.or.resnam.eq.'CYS '.OR.
     1	   resnam.eq.'PRO ')ranksize=7 
	if(resnam.eq.'THR '.or.resnam.eq.'VAL '.or.
     1resnam.eq.'PCA ')ranksize=8
	if(resnam.eq.'ASP '.or.resnam.eq.'LEU '.or.
     1     resnam.eq.'ILE '.or.resnam.eq.'MET '.or.
     2     resnam.eq.'ASN ')ranksize=9
	if(resnam.eq.'GLU '.or.resnam.eq.'LYS '.or.
     1     resnam.eq.'GLN ')ranksize=10
	if(resnam.eq.'HIS ')ranksize=11
	if(resnam.eq.'PHE '.or.resnam.eq.'ARG '.or.
     1     resnam.eq.'HMS ')ranksize=12
	if(resnam.eq.'TYR ')ranksize=13
	if(resnam.eq.'TRP ')ranksize=15

c determine ranking order of all atoms per residue
	i0=1
	if(ires.ne.1)i0=rlist1(ires-1)+1
	imax=rlist1(ires)
	do iatom=i0,imax
c if there's problems, i may want to check that ranksize=rlist1(ires)-i0

c exception PRO
	 If(resnam.eq.'PRO '.or.resnam.eq.'PCA ')Then   ! no MC amide hydrogen
	    if(asort(iatom).eq.' N   ')rank(1)=iatom
	    if(asort(iatom).eq.' CA  ')rank(2)=iatom
	    if(asort(iatom).eq.' CB  ')rank(3)=iatom
	    if(asort(iatom).eq.' CG  ')rank(4)=iatom
	    if(asort(iatom).eq.' CD  ')rank(5)=iatom
	    if(resnam.eq.'PRO '.and.asort(iatom).eq.' C   ')rank(6)=iatom
	    if(resnam.eq.'PRO '.and.asort(iatom).eq.' O   ')rank(7)=iatom
	    if(resnam.eq.'PCA '.and.asort(iatom).eq.' OE  ')rank(6)=iatom
	    if(resnam.eq.'PCA '.and.asort(iatom).eq.' C   ')rank(7)=iatom
	    if(resnam.eq.'PCA '.and.asort(iatom).eq.' O   ')rank(8)=iatom



	 Else
		if(asort(iatom).eq.' N   ') rank(1)=iatom
		if(asort(iatom).eq.' H   ') rank(2)=iatom
		if(asort(iatom).eq.' CA  ') rank(3)=iatom
		if(asort(iatom).eq.' CB  ') rank(4)=iatom
		if(asort(iatom).eq.' C   ') rank(ranksize-1)=iatom
		if(asort(iatom).eq.' O   ') rank(ranksize)=iatom
		if(asort(iatom)(3:4).eq.'G ') rank(5)=iatom
	if(resnam.ne.'ILE '.and.asort(iatom)(3:4).eq.'G1') rank(5)=iatom
	if(resnam.ne.'ILE '.and.asort(iatom)(3:4).eq.'G2') rank(6)=iatom

	 	if( (asort(iatom)(3:4).eq.'D '.or.asort(iatom)(3:4)
     1             .eq.'D1').AND.(resnam.ne.'ILE ') ) rank(6)=iatom

	 	if( (asort(iatom)(3:4).eq.'D1').AND.
     1              (resnam.eq.'ILE ') ) rank(7)=iatom

c 1/2017: i have to reverse CG1 and CG2 ordering for ILE!
                if(resnam.eq.'ILE ')then
                        if(asort(iatom)(3:4).eq.'G1')rank(6)=iatom
                        if(asort(iatom)(3:4).eq.'G2')rank(5)=iatom
                endif

	    	if(asort(iatom)(3:4).eq.'D2') rank(7)=iatom
		if(asort(iatom)(3:4).eq.'E ') rank(7)=iatom

		if( (asort(iatom)(3:4).eq.'E1').AND.
     1              (resnam(1:2).eq.'GL') ) rank(7)=iatom

		if( (asort(iatom)(3:4).eq.'E2').AND.
     1              (resnam(1:2).eq.'GL') ) rank(8)=iatom

		if( (asort(iatom)(3:4).eq.'E1').AND.
     1              (resnam.eq.'PHE '.or.resnam.eq.'HIS '.or.
     2               resnam.eq.'TRP '.or.resnam.eq.'TYR '.or.
     3  	     resnam.eq.'HMS ') )rank(8)=iatom

		if( (asort(iatom)(3:4).eq.'E2').AND.
     1              (resnam.eq.'PHE '.or.resnam.eq.'HIS '.or.
     2               resnam.eq.'TRP '.or.resnam.eq.'TYR '.or.
     3               resnam.eq.'HMS ') ) rank(9)=iatom

		if(asort(iatom)(3:4).eq.'E3') rank(10)=iatom  ! for TRP

		if(asort(iatom).eq.' CM  ') rank(10)=iatom  ! for HMS

		if( (asort(iatom)(3:4).eq.'Z ').AND.
     1 		(resnam.eq.'LYS '.or.resnam.eq.'ARG ') ) rank(8)=iatom

		if((asort(iatom)(3:4).eq.'Z ').AND.
     1 		(resnam.eq.'PHE '.or.resnam.eq.'TYR '))rank(10)=iatom

		if(asort(iatom)(3:4).eq.'Z2')rank(11)=iatom
		if(asort(iatom)(3:4).eq.'Z3')rank(12)=iatom
		if(asort(iatom)(2:3).eq.'OH')rank(11)=iatom
		if(asort(iatom)(3:4).eq.'H1')rank(9)=iatom

		if(asort(iatom)(3:4).eq.'H2'.AND.
     1 		resnam.eq.'ARG ')rank(10)=iatom

		if(asort(iatom)(3:4).eq.'H2'.AND.
     1 		resnam.eq.'TRP ')rank(13)=iatom
	 Endif
	enddo

c correct for absence of Residue 1 amide hydrogen:
	if(ires.eq.1.and.resnam.ne.'PRO '.and.
     1resnam.ne.'PCA ')then
	 ranksize=ranksize-1
	 do i = 2, ranksize
	  rank(i)=rank(i+1)
	 enddo
	endif

c now put arrays xyz, asort,rsort into correct order...
	do irank=1,ranksize
c	write(*,*)irank,ires,resnam,ranksize,i0,
c     1rlist1(ires),rank(irank)
	 tempxyz(1,irank) = xyz(1,rank(irank))
	 tempxyz(2,irank) = xyz(2,rank(irank))
	 tempxyz(3,irank) = xyz(3,rank(irank))
	 tempasort(irank) = asort(rank(irank))
	 temprsort(irank) = rsort(rank(irank))
	 tempB(irank)     = B(rank(irank))
	enddo
c and put back into premanent arrays:
	irank=0
	do iatom=i0,imax
	 irank = irank+1
c	write(*,*)iatom,irank,temprsort(irank),tempasort(irank)
	 xyz(1,iatom) = tempxyz(1,irank)
	 xyz(2,iatom) = tempxyz(2,irank)
	 xyz(3,iatom) = tempxyz(3,irank)
	 asort(iatom) = tempasort(irank)
	 rsort(iatom) = temprsort(irank)
	 B(iatom)     = tempB(irank)
	enddo


	ENDDO

c i may want to doublecheck that this leaves OXT in proper place and intact.

	return
	end

