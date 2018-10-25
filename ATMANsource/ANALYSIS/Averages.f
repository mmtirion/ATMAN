	program Averages

	implicit none
	integer*4 Nprtn,Nres,Npro,Nxi,NintraP
	integer*4 Nlig,NligNA,NDoF,Nmodes,rscale
	real*8 delta,pi,Kb,Kelvin,enhance,sum,Ave
c for protein:
	parameter ( Nprtn = 2598 ) 
	parameter ( Nres  = 302 ) 
	parameter ( Nxi   = 478 )
	parameter ( rscale = 0 )
	parameter ( Nlig  = 0 )
	parameter ( NligNA  = 0 )
	parameter ( NDoF  = 2*Nres+Nxi ) 
	parameter ( Nmodes = 100 )
	parameter ( delta = 1.d-05  ) 
	parameter ( Kelvin= 293.d0  )
        parameter ( Kb    = 1.3806d-23 ) !boltzmann constant m^2kgs^-2K^-1
	parameter ( pi    = dacos(-1.d0) )

	integer*4 IndexMC(2*Nres,4)
	integer*4 IndexXi(Nxi,4),EndXi(Nxi),XiRes(Nxi)
	real*8 xyzP(3,Nprtn),massP(Nprtn),B(Nprtn)
	integer*4 rlist(Nprtn),rnum(Nres),iatom
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
	real*8 vec(3),Occ,total,temp(Nprtn)
	integer i,j,k,N,jj,istart,number,Pyn(2*Nres+Nxi)
	real*8 A(NDoF,Nmodes),lambda,omega,alpha(Nmodes)
	real*8 RMSik(Nprtn,Nmodes),SigmaCa(Nprtn),Bfac
	real*8 num,den,Btheo,Bexp,correlation

10      format(a6,x,i4,x,a5,a4,1x,i4,4x,3f8.3,2x,f4.2,x,f5.2)

c===================================================================
c  index torsions; identify atoms, masses, residues; get coordinates
c===================================================================
	call index_prtn(indexMC,indexXi,massP,xyzP,asort,B,rlist,
     1rsort,rnum,type,EndXi,XiRes,Nres,Nxi,Nprtn,Rhc,aHET,rHet,
     2rlstHET,mHET,vdWHET,typeHET,epsHET,NligNA,rscale,Bhet,P,Pyn)

c record Exp. B factors
c	open(14,file='./Exp_Bfactors')
c	DO j = 1, Nres
c	 Temp = 0.0
c	 istart = 1
c	 if(j.ne.1)istart=rnum(j-1)+1
c
c	 do i = istart, rnum(j)
c	  Temp = Temp + B(i)
c	 enddo
c	 write(14,*)j,Temp/real(rnum(j)-istart+1)
c	ENDDO
c
c	sum = 0.0
c	number = 0
c	do i = 1, Nprtn
c
c   	 if(asort(i).eq.' CA  ')sum = sum + B(i)
c   	 if(asort(i).eq.' CA  ')number = number+1
c
c	enddo
c
c   	sum = sum /(rlist(Nprtn)-rlist(1)+1)
c	write(*,*)'Average experimental B value:',sum/real(number)
c	write(*,*)' '
c
	open(15,file='./RESULTS/Bexp')
	do i = 1, Nprtn
 	 if(asort(i).eq.' CA  ')then
 		write(15,*)rscale+rlist(i),',',B(i)
 	 endif
 	enddo
 	close(15)
c
c=========================================================================
c Compute dRdTor array for all Calpha atoms, and all DoF
c=========================================================================
	call dRdTorsion(massP,xyzP,IndexMC,IndexXi,EndXi,Nres,
     1  Nxi,Npro,Nprtn,asort,delta,NDoF,Pyn,dRdTor)


c=========================================================================
c  read in eigenvalues and eigenvectors A(NDoF,nmodes):
c=========================================================================
c        open(81,status='unknown',file=
c     1'/Users/monique/Desktop/WORK/NMStuff/ATMAN/DATA/eigenvectors')
        open(81,status='unknown',file='../DATA/eigenvectors')

c	enhance = 1.0d0
c	write(*,*)'ENHANCE=',real(enhance)
	write(*,*)' Enhance?'
	read(*,*)enhance
        DO i = 1,Nmodes
           read(81,*)k, omega  ! in cm^-1
           alpha(i) = sqrt(2.0*Kb*Kelvin) / omega  ! in [m*sqrt(kg)]

	   alpha(i) = enhance*alpha(i)  ! manual adjustment for now.

         Do j = 1,2*Nres
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
c  Calculate the rms fluctuation, RMSik, of atom i due to mode k:

        DO i = 1, Nprtn
        do k = 1, Nmodes
         vec(1) = 0.d0
         vec(2) = 0.d0
         vec(3) = 0.d0
        do j = 1, NDoF
         vec(1) = vec(1) + dRdTor(i,j,1)*A(j,k)*alpha(k)/dsqrt(2.d0)
         vec(2) = vec(2) + dRdTor(i,j,2)*A(j,k)*alpha(k)/dsqrt(2.d0)
         vec(3) = vec(3) + dRdTor(i,j,3)*A(j,k)*alpha(k)/dsqrt(2.d0)
        enddo
         RMSik(i,k) = sqrt( vec(1)*vec(1)+vec(2)*vec(2)+
     1			    vec(3)*vec(3) )
        enddo
        ENDDO


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Calculate the rms fluctuation of each Alpha carbon atom due to k nmodes:

	atom='ATOM  '
	Occ = 1.0
        open(81,status='unknown',file='./RESULTS/RMS_atom')


	Ave = 0.0
	iatom = 0
        DO i = 1, Nprtn
	 IF(asort(i).eq.' CA  ')THEN
	  iatom = iatom + 1
          SigmaCa(iatom) = 0.0

          do k = 1, Nmodes
           SigmaCa(iatom) = SigmaCa(iatom) + RMSik(i,k)**2
          enddo

          SigmaCa(iatom) = sqrt(SigmaCa(iatom))
	  Ave = Ave + SigmaCa(iatom)
	  write(81,*)rscale+rlist(i),sngl(SigmaCa(iatom))
	 ENDIF
        ENDDO

	write(*,*)'Average RMS value:',Ave/real(iatom)
	close(81)




        open(81,status='unknown',file='./RESULTS/Btheo')


        DO i = 1, Nprtn
         temp(i) = 0.0d0

         do k = 1, Nmodes
          temp(i) = temp(i) + (8.d0*pi*pi)*RMSik(i,k)**2
         enddo

           if(asort(i).eq.' CA  ')
     1write(81,*)rscale+rlist(i),sngl(temp(i))


	ENDDO

	close(81)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  Calculate the B value per Alpha carbon:
c        open(81,status='unknown',file=
c     1'/Users/monique/Desktop/WORK/NMStuff/ATMAN/DATA/Bfactors')
c
c        do i = 1,Nres
c         Bfac = (8.*pi*pi/3.)*SigmaCa(i)**2
c         write(81,*)i,Bfac
c        enddo
c        close(81)

        open(81,status='unknown',file='./RESULTS/RMS_mode')

        DO k = 1,Nmodes
         Bfac = 0.0
	 N = 0

         do i = 1, Nprtn
	  if(asort(i).eq.' CA  ')then
		Bfac = Bfac + RMSik(i,k)**2
		N = N + 1
	  endif
         enddo

         Bfac = sqrt(Bfac/real(N))
         write(81,*)k,',',Bfac
        ENDDO


	stop
	end
