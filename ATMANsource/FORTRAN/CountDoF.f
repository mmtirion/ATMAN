	subroutine CountDoF(ResID,NB,NA,NDE,Nres)

c subroutine that counts internal degrees of freedom: # bonds, # bond angles plus # dihedrals and 
c compares them to # of Cartesian degrees of freedom. Includes only heavy atoms in analysis; ie no hydrogen atoms.

	implicit none
	character*4 ResID
	integer NB,NA,NDE,Nres(21)

	if(ResID.eq.'ALA ')then
		NB = 4
		NA = 4
		NDE = 1
		Nres(1) = Nres(1) + 1
	endif
	if(ResID.eq.'GLY ')then
		NB = 3
		NA = 2
		NDE = 0
		Nres(2) = Nres(2) + 1
	endif
	if(ResID.eq.'ILE ')then
		NB = 7
		NA = 8
		NDE = 2
		Nres(3) = Nres(3) + 1
	endif
	if(ResID.eq.'LEU ')then
		NB = 7
		NA = 8
		NDE = 2
		Nres(4) = Nres(4) + 1
	endif
	if(ResID.eq.'PRO ')then
		NB = 7
		NA = 3
		NDE = 0
		Nres(5) = Nres(5) + 1
	endif
	if(ResID.eq.'PHE ')then
		NB = 11
		NA = 7
		NDE = 0
		Nres(6) = Nres(6) + 1
	endif
	if(ResID.eq.'VAL ')then
		NB = 6
		NA = 7
		NDE = 2
		Nres(7) = Nres(7) + 1
	endif
	if(ResID.eq.'ARG ')then
		NB = 10
		NA = 11
		NDE = 2
		Nres(8) = Nres(8) + 1
	endif
	if(ResID.eq.'ASP ')then
		NB = 7
		NA = 8
		NDE = 2
		Nres(9) = Nres(9) + 1
	endif
	if(ResID.eq.'GLU ')then
		NB = 8
		NA = 9
		NDE = 2
		Nres(10) = Nres(10) + 1
	endif
	if(ResID.eq.'SER ')then
		NB = 5
		NA = 5
		NDE = 1
		Nres(11) = Nres(11) + 1
	endif
	if(ResID.eq.'THR ')then
		NB = 6
		NA = 7
		NDE = 2
		Nres(12) = Nres(12) + 1
	endif
	if(ResID.eq.'CYS ')then
		NB = 5
		NA = 5
		NDE = 1
		Nres(13) = Nres(13) + 1
	endif
	if(ResID.eq.'ASN ')then
		NB = 7
		NA = 8
		NDE = 2
		Nres(14) = Nres(14) + 1
	endif
	if(ResID.eq.'GLN ')then
		NB = 8
		NA = 9
		NDE = 2
		Nres(15) = Nres(15) + 1
	endif
	if(ResID.eq.'HIS ')then
		NB = 10
		NA = 6
		NDE = 0
		Nres(16) = Nres(16) + 1
	endif
	if(ResID.eq.'LYS ')then
		NB = 8
		NA = 8
		NDE = 1
		Nres(17) = Nres(17) + 1
	endif
	if(ResID.eq.'TYR ')then
		NB = 12 
		NA = 9
		NDE = 1
		Nres(18) = Nres(18) + 1
	endif
	if(ResID.eq.'MET ')then
		NB = 7
		NA = 7
		NDE = 1
		Nres(19) = Nres(19) + 1
	endif
	if(ResID.eq.'TRP ')then
		NB = 15
		NA = 8
		NDE = 0
		Nres(20) = Nres(20) + 1
	endif
	if(ResID.eq.'PCA ')then
		NB = 8 
		NA = 4
		NDE = 0  !  For terminal PCA only!
		Nres(21) = Nres(21) + 1
	endif


	return
	end
