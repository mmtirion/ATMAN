
	subroutine amu(TYPE,m,vdW)
c These params from M.Levitt's 1983 JMB paper 
c hydrogen atoms are included only on the oxygen and nitrogen atoms of peptide groups and
c amid groups.

	character*1 TYPE
	real*8 m,vdW


	if(TYPE.eq.'H')then  ! hydrogen
		m =1.008d0
		vdW = 2.8525d0
		return
	endif
	if(TYPE.eq.'O')then   ! oxygen
		m = 15.999d0
		vdW =3.1005d0
		return
	endif
	if(TYPE.eq.'Q')then  !oxygen in a carboxyl group
		m = 15.999d0
		vdW = 3.1005d0
		return
	endif
	if(TYPE.eq.'V')then  !An OH group treated as a single extended atom
		m = 15.999d0
		vdW = 3.1005d0
		return
	endif
	if(TYPE.eq.'N')then   ! nitrogen
		m = 14.007d0
		vdW = 3.8171d0
		return
	endif
	if(TYPE.eq.'M')then  ! nitrogen in a NH, NH2, NH3 group of Lys or Arg
		m = 14.007d0
		vdW = 3.8171d0
		return
	endif
	if(TYPE.eq.'C')then   !tetrahedral carbon
		m = 12.011d0
		vdW = 4.3150d0
		return
	endif
	if(TYPE.eq.'S')then  !sulfur
		m = 32.064d0
		vdW = 4.3150d0
		return
	endif
	if(TYPE.eq.'A')then  !trigonal carbon
		m = 12.011d0
		vdW = 4.2202d0
		return
	endif
	if(TYPE.eq.'B')then  !trigonal carbon
		m = 12.011d0
		vdW = 4.2202d0
		return
	endif
	if(TYPE.eq.'G')then  !trigonal carbon
		m = 12.011d0
		vdW = 4.2202d0
		return
	endif



		write(*,*)' Protein-atom mass, vdW radius unassigned'
		write(*,*)' TYPE: ',TYPE
		stop
 
	end


	subroutine sort(ATOM,RESI,TYPE)
	character*4 RESI
	character*5 ATOM
	character*1 TYPE

	if(ATOM.eq.' N   ')then
		TYPE='N' 
		return
	endif

	if(ATOM.eq.' H   ')then
		TYPE='H'
		return
	endif

	if(ATOM.eq.' C   ')then
		TYPE='A'
		return
	endif

	if(ATOM.eq.' O   ')then
		TYPE='O'
		return
	endif

	if(ATOM.eq.' OT  ')then
		TYPE='O'
		return
	endif

	if(ATOM.eq.' OXT ')then
		TYPE='O'
		return
	endif

	if(ATOM.eq.' OCT1')then
		TYPE='O'
		return
	endif

	if(ATOM.eq.' OCT2')then
		TYPE='O'
		return
	endif

	if(ATOM(1:1).eq.'H'.OR.ATOM(2:2).eq.'H')then
		TYPE='H'
		return
	endif



	if(RESI.eq.'ALA ')then
  		if(ATOM.eq.' CA  ')TYPE='C' 
  		if(ATOM.eq.' CB  ')TYPE='C' 
		return
	endif

	if(RESI.eq.'ARG ')then
  		if(ATOM.eq.' CA  ' )TYPE='C'
  		if(ATOM.eq.' CB  ' )TYPE='C'
  		if(ATOM.eq.' CG  ' )TYPE='C'
	  	if(ATOM.eq.' CD  ' )TYPE='C'
  		if(ATOM.eq.' NE  ' )TYPE='M'
  		if(ATOM.eq.' HE  ' )TYPE='H'
  		if(ATOM.eq.' CZ  ' )TYPE='A'
  		if(ATOM.eq.' NH1 ' )TYPE='M'
  		if(ATOM.eq.' HH11' )TYPE='H'
  		if(ATOM.eq.' HH12' )TYPE='H'
  		if(ATOM.eq.' NH2 ' )TYPE='M'
  		if(ATOM.eq.' HH21' )TYPE='H'
  		if(ATOM.eq.' HH22' )TYPE='H'
		return
	endif

	if(RESI.eq.'ASN ')then
  		if(ATOM.eq.' CA  ')TYPE='C'
  		if(ATOM.eq.' CB  ')TYPE='C'
  		if(ATOM.eq.' CG  ')TYPE='A'
  		if(ATOM.eq.' OD1 ')TYPE='O'
  		if(ATOM.eq.' ND2 ')TYPE='N'
  		if(ATOM.eq.' AD1 ')TYPE='O' !Oct 07
  		if(ATOM.eq.' AD2 ')TYPE='N' !Oct 07
  		if(ATOM.eq.' HD21')TYPE='H'
  		if(ATOM.eq.' HD22')TYPE='H'
		return
	endif

	if(RESI.eq.'ASP ')then
  		if(ATOM.eq.' CA  ')TYPE='C'
  		if(ATOM.eq.' CB  ')TYPE='C'
  		if(ATOM.eq.' CG  ')TYPE='A'
  		if(ATOM.eq.' OD1 ')TYPE='Q'
  		if(ATOM.eq.' OD2 ')TYPE='Q'
		return
	endif

	if(RESI.eq.'CYS ')then
  		if(ATOM.eq.' CA  ')TYPE='C'
  		if(ATOM.eq.' CB  ')TYPE='C'
  		if(ATOM.eq.' SG  ')TYPE='S'
		return
	endif

	if(RESI.eq.'GLN ')then
  		if(ATOM.eq.' CA  ')TYPE='C'
  		if(ATOM.eq.' CB  ')TYPE='C'
  		if(ATOM.eq.' CG  ')TYPE='C'
  		if(ATOM.eq.' CD  ')TYPE='A'
  		if(ATOM.eq.' OE1 ')TYPE='O'
  		if(ATOM.eq.' NE2 ')TYPE='N'
  		if(ATOM.eq.' AE1 ')TYPE='O' !Oct 07
  		if(ATOM.eq.' AE2 ')TYPE='N' !Oct 07
  		if(ATOM.eq.' HE21')TYPE='H'
  		if(ATOM.eq.' HE22')TYPE='H'
		return
	endif

	if(RESI.eq.'GLU ')then
  		if(ATOM.eq.' CA  ')TYPE='C'
  		if(ATOM.eq.' CB  ')TYPE='C'
  		if(ATOM.eq.' CG  ')TYPE='C'
  		if(ATOM.eq.' CD  ')TYPE='A'
  		if(ATOM.eq.' OE1 ')TYPE='Q'
  		if(ATOM.eq.' OE2 ')TYPE='Q'
		return
	endif

	if(RESI.eq.'GLY ')then
  		if(ATOM.eq.' CA  ')TYPE='C'
		return
	endif

	if(RESI.eq.'HIS ')then
  		if(ATOM.eq.' CA  ')TYPE='C'
  		if(ATOM.eq.' CB  ')TYPE='C'
  		if(ATOM.eq.' CG  ')TYPE='B'
  		if(ATOM.eq.' CD2 ')TYPE='B'
  		if(ATOM.eq.' ND1 ')TYPE='N'
  		if(ATOM.eq.' HD1 ')TYPE='H'
  		if(ATOM.eq.' CE1 ')TYPE='B'
  		if(ATOM.eq.' NE2 ')TYPE='N'
  		if(ATOM.eq.' HE2 ')TYPE='H'
		return
	endif

	if(RESI.eq.'HMS ')then
  		if(ATOM.eq.' CA  ')TYPE='C'
  		if(ATOM.eq.' CB  ')TYPE='C'
  		if(ATOM.eq.' CG  ')TYPE='B'
  		if(ATOM.eq.' CD2 ')TYPE='B'
  		if(ATOM.eq.' ND1 ')TYPE='N'
  		if(ATOM.eq.' HD1 ')TYPE='H'
  		if(ATOM.eq.' CE1 ')TYPE='B'
  		if(ATOM.eq.' NE2 ')TYPE='N'
  		if(ATOM.eq.' CM  ')TYPE='C'
		return
	endif

	if(RESI.eq.'ILE ')then
  		if(ATOM.eq.' CA  ')TYPE='C'
  		if(ATOM.eq.' CB  ')TYPE='C'
  		if(ATOM.eq.' CG2 ')TYPE='C'
  		if(ATOM.eq.' CG1 ')TYPE='C'
  		if(ATOM.eq.' CD  ')TYPE='C'
  		if(ATOM.eq.' CD1 ')TYPE='C' ! Oct 07
		return
	endif

	if(RESI.eq.'LEU ')then
  		if(ATOM.eq.' CA  ')TYPE='C'
  		if(ATOM.eq.' CB  ')TYPE='C'
  		if(ATOM.eq.' CG  ')TYPE='C'
  		if(ATOM.eq.' CD1 ')TYPE='C'
  		if(ATOM.eq.' CD2 ')TYPE='C'
		return
	endif

	if(RESI.eq.'LYS ')then
  		if(ATOM.eq.' CA  ')TYPE='C'
  		if(ATOM.eq.' CB  ')TYPE='C'
  		if(ATOM.eq.' CG  ')TYPE='C'
  		if(ATOM.eq.' CD  ')TYPE='C'
  		if(ATOM.eq.' CE  ')TYPE='C'
  		if(ATOM.eq.' NZ  ')TYPE='M'
  		if(ATOM.eq.' HZ1 ')TYPE='H'
  		if(ATOM.eq.' HZ2 ')TYPE='H'
  		if(ATOM.eq.' HZ3 ')TYPE='H'
		return
	endif

	if(RESI.eq.'PCA ')then
  		if(ATOM.eq.' CA  ')TYPE='C'
  		if(ATOM.eq.' CB  ')TYPE='C'
  		if(ATOM.eq.' CG  ')TYPE='C'
  		if(ATOM.eq.' CD  ')TYPE='C'
  		if(ATOM.eq.' OE  ')TYPE='V'
		return
	endif

	if(RESI.eq.'MET ')then
  		if(ATOM.eq.' CA  ')TYPE='C'
  		if(ATOM.eq.' CB  ')TYPE='C'
  		if(ATOM.eq.' CG  ')TYPE='C'
  		if(ATOM.eq.' SD  ')TYPE='S'
  		if(ATOM.eq.' CE  ')TYPE='C'
		return
	endif

	if(RESI.eq.'PHE ')then
  		if(ATOM.eq.' CA  ')TYPE='C'
  		if(ATOM.eq.' CB  ')TYPE='C'
  		if(ATOM.eq.' CG  ')TYPE='B'
  		if(ATOM.eq.' CD1 ')TYPE='A'
  		if(ATOM.eq.' CD2 ')TYPE='A'
  		if(ATOM.eq.' CE1 ')TYPE='A'
  		if(ATOM.eq.' CE2 ')TYPE='A'
  		if(ATOM.eq.' CZ  ')TYPE='A'
		return
	endif

	if(RESI.eq.'PRO ')then
  		if(ATOM.eq.' N   ')TYPE='N'
  		if(ATOM.eq.' CD  ')TYPE='C'
  		if(ATOM.eq.' CA  ')TYPE='C'
  		if(ATOM.eq.' CB  ')TYPE='C'
  		if(ATOM.eq.' CG  ')TYPE='C'
		return
	endif

	if(RESI.eq.'SER ')then
  		if(ATOM.eq.' CA  ')TYPE='C'
  		if(ATOM.eq.' CB  ')TYPE='C'
  		if(ATOM.eq.' OG  ')TYPE='V'
  		if(ATOM.eq.' HG  ')TYPE='H'
		return
	endif

	if(RESI.eq.'THR ')then
  		if(ATOM.eq.' CA  ')TYPE='C'
  		if(ATOM.eq.' CB  ')TYPE='C'
  		if(ATOM.eq.' OG1 ')TYPE='V'
  		if(ATOM.eq.' HG1 ')TYPE='H'
  		if(ATOM.eq.' CG2 ')TYPE='C'
		return
	endif

	if(RESI.eq.'TRP ')then
  		if(ATOM.eq.' CA  ')TYPE='C'
  		if(ATOM.eq.' CB  ')TYPE='C'
  		if(ATOM.eq.' CG  ')TYPE='B'
  		if(ATOM.eq.' CD2 ')TYPE='B'
  		if(ATOM.eq.' CE2 ')TYPE='G'
  		if(ATOM.eq.' CE3 ')TYPE='A'
  		if(ATOM.eq.' CD1 ')TYPE='A'
  		if(ATOM.eq.' NE1 ')TYPE='N'
  		if(ATOM.eq.' HE1 ')TYPE='H'
  		if(ATOM.eq.' CZ2 ')TYPE='A'
  		if(ATOM.eq.' CZ3 ')TYPE='A'
  		if(ATOM.eq.' CH2 ')TYPE='G'
		return
	endif

	if(RESI.eq.'TYR ')then
  		if(ATOM.eq.' CA  ')TYPE='C'
  		if(ATOM.eq.' CB  ')TYPE='C'
  		if(ATOM.eq.' CG  ')TYPE='B'
  		if(ATOM.eq.' CD1 ')TYPE='A'
  		if(ATOM.eq.' CE1 ')TYPE='A'
  		if(ATOM.eq.' CD2 ')TYPE='A'
  		if(ATOM.eq.' CE2 ')TYPE='A'
  		if(ATOM.eq.' CZ  ')TYPE='A'
  		if(ATOM.eq.' OH  ')TYPE='V'
  		if(ATOM.eq.' HH  ')TYPE='H'
		return
	endif

	if(RESI.eq.'VAL ')then
  		if(ATOM.eq.' CA  ')TYPE='C'
  		if(ATOM.eq.' CB  ')TYPE='C'
  		if(ATOM.eq.' CG1 ')TYPE='C' 
  		if(ATOM.eq.' CG2 ')TYPE='C'  
		return
	endif

	write(*,*)'ERROR: unassigned prtn-atom TYPE in subroutine sort:'
	write(*,*)'Residue:',RESI,' has atom type:',ATOM
	stop
	end



	subroutine CM(xyz,mass,N,Rcm)

	implicit none
	integer*4 N,i
	real*8 xyz(3,N),mass(N),Rcm(3),Tmass

        Tmass  = 0.0
        Rcm(1) = 0.0
        Rcm(2) = 0.0
        Rcm(3) = 0.0
        do i=1,N
        Rcm(1)=Rcm(1)+mass(i)*xyz(1,i) 
        Rcm(2)=Rcm(2)+mass(i)*xyz(2,i) 
        Rcm(3)=Rcm(3)+mass(i)*xyz(3,i) 
        Tmass=Tmass+mass(i)
        enddo
        Rcm(1)=Rcm(1)/Tmass
        Rcm(2)=Rcm(2)/Tmass
        Rcm(3)=Rcm(3)/Tmass

	return
	end



	subroutine HETmass(asort,rsort,mass,vdW,eps)

	character*5 asort
	character*4 rsort
	real*8 mass,vdW,eps

c these values i  made up, based on L79
c NOTE: I HAVE MADE NO EFFORT TO IDENTIFY VALUES FOR EPSILONS--THIS NEEDS
c SOME SORT OF DATA MINING TO RESOLVE: EACH CASE INDEPENDENTLY

	if(asort.eq.' CA  '.or.asort.eq.'CA   ')then  ! this is for Calcium!! not C-alpha
                mass = 40.08d0
                vdW = 4.3d0
		eps =0.2d0 
		return
	endif

	if(asort.eq.' NA  '.or.asort.eq.'NA   ')then
                mass = 22.99d0
                vdW = 4.7d0
		eps =0.2d0 
		return
	endif

	if(asort.eq.' FE  '.or.asort.eq.'FE   ')then
                mass = 55.85d0
                vdW = 3.0d0
		eps =0.2d0 
		return
	endif

	if(asort.eq.' ZN  '.or.asort.eq.'ZN   ')then
                mass = 65.38d0
                vdW = 3.0d0
		eps =0.2d0 
		return
	endif

	if(asort.eq.' MG  '.or.asort.eq.'MG   ')then
                mass = 24.31d0
                vdW = 4.3d0
		eps =0.2d0 
		return
	endif

	if(rsort.eq.'HOH ')then
		mass = 18.0d0
		vdW = 3.55322d0  !from ENCAD 95
		eps = 0.18479    !from ENCAD 95
		return
	endif

	if(asort.eq.' CL  '.or.asort.eq.'CL   ')then
		mass =35.45d0
		vdW = 4.3d0
		eps =0.2d0 
		return
	endif

	if(asort(2:2).eq.'O')then
		mass =15.999d0
		vdW = 3.7d0
		eps =0.185d0 
		return
	endif

	if(asort(1:1).eq.'C'.OR.asort(2:2).eq.'C')then
		mass =12.011d0
		vdW = 4.3d0
		eps =0.074d0 
		return
	endif

	if(asort(2:2).eq.'N')then
		mass =14.007d0
		vdW = 3.8d0
		eps =0.413d0 
		return
	endif

	if(asort(2:2).eq.'P')then
		mass =30.97d0
		vdW = 4.30d0
		eps =0.2d0 
		return
	endif
	if(asort(2:2).eq.'S')then
		mass = 32.064d0
		vdW = 4.3150d0
		eps =0.2d0 
		return
	endif

	write(*,*)'Unknown HETATM:',rsort
	stop
	end
