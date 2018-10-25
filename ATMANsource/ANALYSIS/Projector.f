	program Projector

	implicit none
	integer*4 Nprtn,Nres,Npro,Nxi,NintraP
	integer*4 Nlig,NligNA,NDoF,Nmodes,rscale
	real*8 delta,pi,Kb,Kelvin,enhance,sum
c for protein:
	parameter ( Nprtn = 2264 ) 
	parameter ( Nres  = 263 ) 
	parameter ( Nxi   = 465 )
	parameter ( Nlig  = 0 )
	parameter ( NligNA  = 0 )
	parameter ( NDoF  = 2*Nres+Nxi ) 
	parameter ( rscale = 25 )
	parameter ( Nmodes = 50 )
	parameter ( delta = 1.d-05  ) 
	parameter ( Kelvin= 300.d0  )
	parameter ( enhance = 1.0 )
        parameter ( Kb    = 1.3806d-23 ) !boltzmann constant m^2kgs^-2K^-1
	parameter ( pi    = dacos(-1.d0) )

	integer*4 IndexMC(2*Nres,4)
	integer*4 IndexXi(Nxi,4),EndXi(Nxi),XiRes(Nxi)
	real*8 xyzP(3,Nprtn),massP(Nprtn),B(Nprtn)
	integer*4 rlist(Nprtn),rnum(Nres)
	character*4 rsort(Nprtn)
	character*5 asort(Nprtn)
	character*1 type(Nprtn)
	character*6 atom
c for ligands:
	character*5 aHET(2500)
	character*4 typeHET(2500)
	integer rlstHET(2500)
	real*8 rHet(3,2500),mHET(2500),vdWHET(2500)
	real*8 Rhc(Nprtn),epsHET(2500),Bhet(2500)
c for computations:
	real*8 dRdTor(Nprtn,NDoF,3),P(2*Nres+Nxi)
	real*8 vec(Nprtn,3),Occ,total,temp(Nprtn)
	real*8 Fluctuation,mii,mjj,mij,max,dist
	integer record,idist
	integer i,j,k,N,jj,istart,number,Pyn(2*Nres+Nxi)
	integer ARRAY(Nres,Nres)
	real*8 A(NDoF,Nmodes),lambda,omega,alpha(Nmodes)
	real*8 p1,p2,p3,m11,m22,m33

10      format(a6,x,i4,x,a5,a4,1x,i4,4x,3f8.3,2x,f4.2,x,f5.2)

c===================================================================
c  index torsions; identify atoms, masses, residues; get coordinates
c===================================================================
	call index_prtn(indexMC,indexXi,massP,xyzP,asort,B,rlist,
     1rsort,rnum,type,EndXi,XiRes,Nres,Nxi,Nprtn,Rhc,aHET,rHet,
     2rlstHET,mHET,vdWHET,typeHET,epsHET,NligNA,rscale,Bhet,P,Pyn)

c=========================================================================
c Compute dRdTor array for all Calpha atoms, and all DoF
c=========================================================================
	call dRdTorsion(massP,xyzP,IndexMC,IndexXi,EndXi,Nres,
     1  Nxi,Npro,Nprtn,asort,delta,NDoF,Pyn,dRdTor)

c=========================================================================
c  read in eigenvalues and eigenvectors A(NDoF,nmodes):
c=========================================================================
        open(81,status='unknown',file=
     1'/Users/monique/Desktop/WORK/NMStuff/ATMAN/DATA/eigenvectors')

	write(*,*)'ENHANCE=',real(enhance)
        DO i = 1, Nmodes
           read(81,*)k, omega  ! in cm^-1
           alpha(i) = sqrt(2.0*Kb*Kelvin) / omega  ! in [m*sqrt(kg)]

	   alpha(i) = enhance*alpha(i)  ! manual adjustment for now.

         Do j = 1, 2*Nres
	  if(Pyn(j).eq.0)then
	   A(j,i) = 0.0	  ! in 1/[m*sqrt(kg)]
	  else
           read(81,*)A(j,i)
	  endif
	 Enddo
         Do j = 2*Nres+1,2*Nres+Nxi
          read(81,*)A(j,i)
	 Enddo 
	if(Nlig.ne.0)then
         Do j = 1, 6  ! THIS MAY NEED TO BE 3!!
          read(81,*)sum
	 Enddo 
	endif

        ENDDO
        close(81)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Calculate the rms fluctuation of each atom i due to mode k:

	k = 1 ! to be set manually for each mode's projection!!
        DO i = 1, Nprtn
         vec(i,1) = 0.d0
         vec(i,2) = 0.d0
         vec(i,3) = 0.d0
         do j = 1, NDoF
          vec(i,1) =vec(i,1) + dRdTor(i,j,1)*A(j,k)*alpha(k)/dsqrt(2.d0)
          vec(i,2) =vec(i,2) + dRdTor(i,j,2)*A(j,k)*alpha(k)/dsqrt(2.d0)
          vec(i,3) =vec(i,3) + dRdTor(i,j,3)*A(j,k)*alpha(k)/dsqrt(2.d0)
         enddo
        ENDDO


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c compute momentum sum and A, B, C values

	m11 = 0.d0
	m22 = 0.d0
	m33 = 0.d0
	p1 = 0.d0
	p2 = 0.d0
	p3 = 0.d0

        DO i = 1, Nprtn

	 p1 = p1 + massP(i)*vec(i,1)
	 p2 = p2 + massP(i)*vec(i,2)
	 p3 = p3 + massP(i)*vec(i,3)

	 m11 = m11 + vec(i,1)*vec(i,1)
	 m22 = m22 + vec(i,2)*vec(i,2)
	 m33 = m33 + vec(i,3)*vec(i,3)

        ENDDO

	write(*,*)' '
	write(*,*)' '
	write(*,*)'p:',sngl(p1),sngl(p2),sngl(p3)
	write(*,*)' '
	write(*,*)' '
	write(*,*)'A,B,C:',sngl(m11),sngl(m22),sngl(m33)
	write(*,*)' '
	write(*,*)' '





	stop
	end
