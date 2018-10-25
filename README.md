
ATMAN — a first principle, i.e. analytic, computation of PDB-NMA in internal, dihedral, coordinates. 
Derivations by Dani ben-Avraham, coding by Monique Tirion.
Please direct software questions to Monique at mmtirion@clarkson.edu and theory questions to Dani benavraham@clarkson.edu.
Deposited to GitHub in October 2018.


This set of programs/subroutines, written by MMTirion in very basic Fortran77 (though a Fortran95 compiler works just as well. Just don’t compile with all the warnings), computes the NMA of whatever input structure is provided, minimized or otherwise, using dihedral dofs and using dihedral as well as non-bonded stiffness constants.  
The current parameters may be replaced by other values and/or functionals (see Tirion & ben-Avraham Phys Biol 2015).

The parent directory contains 4 folders: 

FORTRAN where the fortran files (endings .f) and compiler module called makeATMAN5 are located; 

PDB where the pdb structures to be analyzed are located;

DATA where the output data of the diagonalization, in file ‘eigenvectors’, is located;  

ANALYSIS where code to analyze the `eigenvector’ data belonging to PDB structure is found. 



FORTRAN
Use command module `makeATMAN5’ to compile/make executable (named a.out in working directory). Requires LAPACK.  
Use of internal, dihedral, coordinates, results in moving derivatives ( d R / d Torsion ), that are eliminated using u3best.f written by Wolfgang Kabsch.
The main program, ATMAN5.f, is for a single polypeptide chain (i.e. a different fortran program exists to compute the PDB-NMA of dimers).
Rigid ligands can be included: i.e. ligands with only translational & rotational DoFs. (Inclusion of ligand internal DoFs is currently #5 on my list of things to do) It may be best to avoid ligands until familiar with workings of ATMAN.

INPUT for a PDB-NMA requires the name of the PDB file (declared in subroutine index_prtn.f) as well as 5 parameter values (declared in main ATMAN5.f):

The PDB file name and location is declared in subroutine index_prtn.f, after the variable declaration list (around line 30)
The PDB-type file requires atleast the main-chain amide hydrogens built-in (no hydrogens on main-chain amides of residue 1 and prolines). 

The 5 parameters declared in the main program ATMAN5.f, within the variable-declaration list:

Nprtn=
Nres=
Npro=
Nxi=
rscale=


Nres is the number of (contiguous, i.e. no breaks) residues in chain.
Npro is the number of (non-residue 1) prolines in chain
rscale is the residue offset (equal to first residue number minus one; so if first residue starts at 1, rscale=0, if first residue=22, rscale=21)

Nprtn is the total number of atoms used in computation. This is equal to # heavy atoms in PDB coordinate file plus Nres - Npro - 1. However, just make a guess (overestimate your guess, otherwise you’ll run into dimension problems with Fortran complaining about array sizes being overwritten), and the program will terminate after declaring the correct value for Nprtn. Correct this value in ATMAN5.f, re-compile and re-run.

Nxi is the number of Xi angles for all the residue’s side chains.  If, for example, all residues are ALA, Nxi=0.  Again, overestimate this value, run the program, and it will terminate after declaring the correct value for Nxi.  Correct this value in parameter list of ATMAN5.f, re-compile and re-run.

This should do it. However, often there are A vs B conformations for isolated residues: dual declarations of atoms is not supported and will crash program. To ID any source of trouble, do a reverse-search (i.e. do a search starting from the end of the subroutine) for ‘write’ in index_prtn.f.  This will point you to two lines of code:

c       write(*,*)irank,ires,resnam,ranksize,i0,
c     1rlist1(ires),rank(irank)

and

c       write(*,*)iatom,irank,temprsort(irank),tempasort(irank)

These 3 lines are commented out: remove the letter c at the start of each of these lines to have them executed: recompile; now the program will report which residues are properly processed, and the final residue that is causing the program to crash.  Usually you’ll find an A and B conformation or some other idiosyncrasy that neats repair/remedy.
Btw, when you remove the letter c, make sure the line starting with the number 1 has that 1 starting at space #6…(a reminder of the days when coding was done using 72-column punch cards…and a reason Fortran went belly-up).

There are other parameter declarations that are not PDB-file dependent…except if you include ligands in analysis.
If rigid ligands are included, they need to be declared in ATMAN5 as well, under the protein parameter list PLUS (this is annoying!) the Lig(Nlig,2) array starting at line 90 of ATMAN5.f



PDB
The input structure needs to be in PDB format and reduced (Atman only uses the main-chain amide hydrogens however.) HETATM should only be for ligands. Any other HETATM in the residue list will not be supported.  Only the 20 standard RESidues are supported, but others may be (and have in past been) added. 



DATA
Deigen.f is the subroutine that takes F (Hessian matrix computed by ana_F_ATMAN.f ) and H (inertia matrix computed by num_H_prtn.f) and  solves the equation of motion  H A Lambda  =  F A .   We adopted notation of Levitt/Sander/Stern 1985, where A is eigenvector column matrix and Lambda is diagonal eigenvalue array. 
In Deigen.f, you can declare which eigenvalues to compute (variable range); whether to compute only eigenvalues or both eigenvalues PLUS eigenvectors (variable jobz); and eigenvalue numbers to be computed (variable iu, currently set to 100 via Nmodes).

The output of the diagonalization is stored in the file DATA/eigenvectors eigenvectors is an ASCII file and therefore potentially very large. If the calculation computed Nmode modes and used NDoF dofs, then the file eigenvectors consists of Nmode * [ NDoF + 1] entries (lines), with the first entry for the eigenFREQUENCY and the remaining NDoF entries for the dihedral angles (in weird units) per mode i, with i = 1, Nmode.
btw, NDoF = 2 x Nres + Nxi - Npro - 1 since each of N residues has two main-chain phi/psi angles as well as Nxi side-chain dofs, minus the proline and residue 1 phi angles.
Unfortunately, this file `eigenvectors’ is pretty useless: interpretation of these values depends entirely on internal software that re-connects each angular update back to Cartesian shifts. 


ANALYSIS
Several things may be done with the eigenvectors output file:

1) Generate sequences of frames for visual animation in PyMol or VMD. Use command module `makeANIM’. Set PDB filename in index_prtn.f and set parameter values in Movie_Prtn_Seq.f. Output file, animation.pdb, is stored in ANALYSIS/RESULTS. animation.db consists of 9 frames (this can be upped to make a smoother animation) suitable for animation sequences in PyMol.  

2) Compute averages (at kT with T currently set at 300K) using some or all of the `eigenvector’ data. Use command module `makeRMS’. Set PDB filename in index_prtn.f and set parameter values in Averages.f, including which modes to use in generating the output files. 
Output files stored in ANALYSIS/RESULTS and include: (a) Experimental B values for C_alpha atoms in Beep; (b) theoretical B values per C_alpha atom in file Btheo; (c) RMS_atom which are RMS per C_alpha; (d) RMS_mode which provide RMS values of all atoms PER mode, for all Nmodes.

3) Generate PDV file (similar to PDB but w/ OCC and B fields replaced by displacement vector due to select mode). This code does not yet exist. Please contact me if needed.

Units: Ahh. Don’t ask. They work and are self-consistent at this point. Everything proceeds from the fact that the energy expressions in Levitt/Sander/Stern are in units of kilcalories/mole and distance is expressed in Ångstroms. Final displacements are reported in Ångstroms, output frequencies reported (on-screen) in inverse centimeters.


EXAMPLE

ATMAN5.f is currently set up to compute, and analyze, the 100 slowest modes of PDB structure 1GOKA. Input file: PDB/1GOK_Ah.db. 
Hydrogens were added using Richardson’s REDUCE.
On my MacBook Pro w/ 2.8 GHz Inter Core i5, the executable requires a little under 12 cpu minutes for this 1082 DoF system. The reason it’s so slow, is that it’s currently set up with four levels of do-loops in the F matrix computation: NDoF, NDoF, Natoms, Natoms. (If I could correctly collapse the Natom x Natom array of non-bonded interactions into a vector plus pointer, it would go much faster. )  If this runs, the first 100 mode frequencies, starting at 5.1829 /cm , should appear on screen.

Please contact Monique at mmtirion@clarkson.edu to report problems. 
