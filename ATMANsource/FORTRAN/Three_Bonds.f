	subroutine Three_Bonds(asort,rsort,rnum,rlist,IntAct,Nres,
     1  Nprtn,Nintra,Nxi,NBij,ijss,index_xi,EndXi,index_mc)

c  Subroutine to filter and exclude (i,j) interactions between atom pairs either  one
c or two bond lengths (technically, w/in one degree of freedom) apart. 6/20/14
c on 9/2/16 altered to store data array in I*1 InterAct(Nprtn,Nprtn) and NBij(nprtn,nprtn)
c becoming Four_Bonds2.f

c Four_Bonds3 differs from Four_Bonds2 in that it ALLOWS NBI within 3 bond lengths.

	implicit none

	integer Nres,Nprtn,Nintra,Nxi,int,remark
        integer i,j,Nin,ijss,IDi,IDj,Idiff,Resk,Resl
	integer iatom,yn,index_xi(Nxi,4),EndXi(Nxi)
	integer index_mc(2*Nres,4)
	integer Acount,Bcount,k
	integer max,min,parse,hist(1000)
        integer*4 rlist(Nprtn),rnum(Nres)
	character*4 rsort(Nprtn)
	character*5 asort(Nprtn)
        integer IntAct(Nprtn,Nprtn),NBij(Nprtn,Nprtn)
	real avg

c Three cases (MC= Mainchain, SC=Sidechain):
c I:  {i} MC_k and {j} MC_l:
c   A. |k-l| > 2
c   B. |k-l| = 1
c   C. |k-l| = 0
c II:  {i} SC_k  and {j} SC_l
c   A. k .ne. l
c   B. k = l
c III: {i} MC_k  and {j} SC_l
c   A. k .ne. l
c   B. k = l




	do i=1, Nprtn
	do j=1, Nprtn
	 NBij(i,j)=0
	enddo
	enddo
	ijss = 0	!# {i,j} subset of Nintra

	do i = 1, Nprtn
	do j = i+1, Nprtn

	if(IntAct(i,j).ne.0)then

	 Resk = rlist(i)  ! IDs residue # of atom (i)
	 Resl = rlist(j)
	 Idiff = iabs(Resk - Resl)
	if(Resk.gt.Resl)then
		write(*,*)'Problem in subroutine Exclude_Neighbors:'
		write(*,*)'Res(k)=',Resk,' > Res(l)=',Resl
		stop
	endif

c  in order to parse interaction into categories:
	IDi = 0   ! atom i assumed to belong to SC (by definition)
	IDj = 0
	if(asort(i).eq.' N   '.or.asort(i).eq.' H   '.or.
     1	   asort(i).eq.' CA  '.or.asort(i).eq.' C   '.or.
     2	   asort(i).eq.' O   ')IDi=1
	if(asort(j).eq.' N   '.or.asort(j).eq.' H   '.or.
     1	   asort(j).eq.' CA  '.or.asort(j).eq.' C   '.or.
     2	   asort(j).eq.' O   ')IDj=1
	
c  Case I:
	IF(IDi.eq.1.and.IDj.eq.1)THEN
	  If(Idiff.ge.2)Then				!cases IA
		ijss = ijss + 1
		NBij(i,j)=1
		NBij(j,i)=1
	  Endif
	  If(Idiff.eq.1)then				!cases IB
cnew	 if(asort(i).eq.' N   '.and.asort(j).eq.' N   ')goto800
	 if(asort(i).eq.' CA  '.and.asort(j).eq.' N   ')goto800
cnew	 if(asort(i).eq.' CA  '.and.asort(j).eq.' H   ')goto800
cnew	 if(asort(i).eq.' CA  '.and.asort(j).eq.' CA  ')goto800
	 if(asort(i).eq.' C   '.and.asort(j).eq.' N   ')goto800
	 if(asort(i).eq.' C   '.and.asort(j).eq.' H   ')goto800
	 if(asort(i).eq.' C   '.and.asort(j).eq.' CA  ')goto800
cnew	 if(asort(i).eq.' C   '.and.asort(j).eq.' C   ')goto800
	 if(asort(i).eq.' O   '.and.asort(j).eq.' N   ')goto800
cnew	 if(asort(i).eq.' O   '.and.asort(j).eq.' H   ')goto800
cnew	 if(asort(i).eq.' O   '.and.asort(j).eq.' CA  ')goto800
c no special considerations here for prolines.
	 ijss = ijss + 1
	 NBij(i,j)=1
	 NBij(j,i)=1
	  Endif
	  If(Idiff.eq.0.and.asort(i).eq.' H   '.and.	!cases IC
     1			    asort(j).eq.' C   ')then
	 	ijss = ijss + 1
	 	NBij(i,j)=1
	 	NBij(j,i)=1
	  Endif
	  If(Idiff.eq.0.and.asort(i).eq.' H   '.and.	!new
     1			    asort(j).eq.' O   ')then    !new
	 	ijss = ijss + 1
	 	NBij(i,j)=1
	 	NBij(j,i)=1
	  Endif
	  If(Idiff.eq.0.and.asort(i).eq.' N   '.and.   !new
     1			    asort(j).eq.' O   ')then   !new
	 	ijss = ijss + 1
	 	NBij(i,j)=1
	 	NBij(j,i)=1
	  Endif
	ENDIF

c Case II:
	IF(IDi.eq.0.and.IDj.eq.0)THEN			
	 If(resk.ne.resl)then				!cases IIA
	  ijss = ijss + 1
	  NBij(i,j)=1
	  NBij(j,i)=1
	 Else						!cases IIB
 	  if(rsort(i).eq.'SER ')goto800
 	  if(rsort(i).eq.'ALA ')goto800
 	  if(rsort(i).eq.'VAL ')goto800
 	  if(rsort(i).eq.'PHE ')goto800
 	  if(rsort(i).eq.'PRO ')goto800
 	  if(rsort(i).eq.'PCA ')goto800
 	  if(rsort(i).eq.'LEU ')goto800
 	  if(rsort(i).eq.'ASP ')goto800
 	  if(rsort(i).eq.'THR ')goto800
 	  if(rsort(i).eq.'CYS ')goto800
 	  if(rsort(i).eq.'ASN ')goto800
 	  if(rsort(i).eq.'TRP ')goto800
 	  if(rsort(i).eq.'HIS ')goto800
 	  if(rsort(i).eq.'TYR ')goto800
 	  if(rsort(i).eq.'GLY ')goto800
cnew	  if(rsort(i).eq.'ILE ')goto800
cnew	  if(rsort(i).eq.'GLU ')goto800
cnew	  if(rsort(i).eq.'GLN ')goto800
cnew	  if(rsort(i).eq.'MET ')goto800

	  if(rsort(i).eq.'ARG ')then
	   if(
     6     ((asort(i).eq.' CB  ').and.(asort(j).eq.' NE  ')).OR.  !new
     1     ((asort(i).eq.' CB  ').and.(asort(j).eq.' CZ  ')).OR.
     2     ((asort(i).eq.' CB  ').and.(asort(j).eq.' NH1 ')).OR.
     3     ((asort(i).eq.' CB  ').and.(asort(j).eq.' NH2 ')).OR.
     7     ((asort(i).eq.' CG  ').and.(asort(j).eq.' CZ  ')).OR.  !new
     4     ((asort(i).eq.' CG  ').and.(asort(j).eq.' NH1 ')).OR.
     5     ((asort(i).eq.' CG  ').and.(asort(j).eq.' NH2 ')).OR.
     8     ((asort(i).eq.' CD  ').and.(asort(j).eq.' NH1 ')).OR.  !new
     9     ((asort(i).eq.' CD  ').and.(asort(j).eq.' NH2 ')))then  !new
     		ijss=ijss+1
	 	NBij(i,j)=1
	 	NBij(j,i)=1
	   endif	   
	  endif
c new
	  if(rsort(i).eq.'ILE '.and. 
     1     asort(i).eq.' CG2 '.and.asort(j).eq.' CD1 ')then
     		ijss=ijss+1
	 	NBij(i,j)=1
	 	NBij(j,i)=1
	  endif
c new
	  if(rsort(i).eq.'GLU '.and.
     1     (asort(i).eq.' CB  '.and.asort(j).eq.' OE1 ').OR. 
     1     (asort(i).eq.' CB  '.and.asort(j).eq.' OE2 '))then
     		ijss=ijss+1
	 	NBij(i,j)=1
	 	NBij(j,i)=1
	  endif
c new
	  if(rsort(i).eq.'GLN '.and. 
     1     (asort(i).eq.' CB  '.and.asort(j).eq.' OE1 ').OR. 
     1     (asort(i).eq.' CB  '.and.asort(j).eq.' OE2 '))then
     		ijss=ijss+1
	 	NBij(i,j)=1
	 	NBij(j,i)=1
	  endif
c new
	  if(rsort(i).eq.'MET '.and. 
     1     (asort(i).eq.' CB  '.and.asort(j).eq.' CE  '))then
     		ijss=ijss+1
	 	NBij(i,j)=1
	 	NBij(j,i)=1
	  endif
c new
	  if(rsort(i).eq.'LYS '.and.
     1     (asort(i).eq.' CB  '.and.asort(j).eq.' CE  ').OR.  
     1     (asort(i).eq.' CG  '.and.asort(j).eq.' NZ  ').OR. 
     1     (asort(i).eq.' CB  '.and.asort(j).eq.' NZ  '))then
     		ijss=ijss+1
	 	NBij(i,j)=1
	 	NBij(j,i)=1
	  endif

	 Endif
	ENDIF
	  
c Case III:
	IF((IDi.eq.1.and.IDj.eq.0).or.(IDi.eq.0.and.IDj.eq.1))THEN

	 If(Idiff.ge.2)Then				!cases IIIA
	  ijss = ijss + 1
	  NBij(i,j)=1
	  NBij(j,i)=1
	 Endif

	 If(Idiff.eq.1)Then                             !more cases IIIA
c	  if(asort(i).eq.' C   '.and.asort(j).eq.' CB  ')goto800         !new
	  if(rsort(i).eq.'PRO '.and.asort(i).eq.' CD  '.and.
     1asort(j).eq.' C   ')goto800
	  if(rsort(j).eq.'PRO '.and.asort(j).eq.' CD  '.and.
     1asort(i).eq.' C   ')goto800
	  ijss = ijss + 1
	  NBij(i,j)=1
	  NBij(j,i)=1
	 Endif

c cases IIIB: ! new
c proline exception:
	if(rsort(i).eq.'GLY '.or.rsort(i).eq.'PRO ')goto800
	if((asort(i).eq.' CB  '.and.asort(j).eq.' N   ').OR.
     1     (asort(j).eq.' CB  '.and.asort(i).eq.' N   '))goto800
	if((asort(i).eq.' CB  '.and.asort(j).eq.' C   ').OR.
     1     (asort(j).eq.' CB  '.and.asort(i).eq.' C   '))goto800
	if((asort(i).eq.' CB  '.and.asort(j).eq.' CA  ').OR.
     1     (asort(j).eq.' CB  '.and.asort(i).eq.' CA  '))goto800
	if((asort(i).eq.' CA  '.and.asort(j)(3:3).eq.'G').OR.
     1     (asort(j).eq.' CA  '.and.asort(i)(3:3).eq.'G'))goto800


	 ijss = ijss + 1
	 NBij(i,j)=1
	 NBij(j,i)=1
	ENDIF



800	continue 
	endif  ! ends if(IntAct(i,j).ne.0) statement
	enddo
	enddo

	write(*,*)' '
	write(*,*)' In Three-Bonds reduced to:' ,ijss
	write(*,*)' '



	return
	end
