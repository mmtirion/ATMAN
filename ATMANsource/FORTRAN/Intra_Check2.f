	subroutine Intra_check2(mass,asort,xyz,rnum,rlist,rsort,NBij,ijss,
     1  Nres,Nprtn,Rhc,C)


	implicit none

	character*5 space
	integer Nres,Nprtn,ijss,N0

	integer h,i,j,k,l,Irow,Icol
	integer mcmc,mcsc,scsc,num,nbi
	integer Iplus,Iminus,Rbin(5,5)
	real*8 mass(Nprtn),xyz(3,Nprtn)
	integer rlist(Nprtn),rnum(Nres)
	character*4 rsort(Nprtn)
	character*5 asort(Nprtn)
	real*8 Rmass(Nres),Rcm(3,Nres),C
	integer NBij(Nprtn,Nprtn)
	integer Imax,Imin,Iave,Nint(Nprtn),Hist(1000)
	real*8 Dist_ave,Dist_min,Dist_max,HD(5),dist,comp
	real*8 dmin(Nprtn),dave,Rhc(Nprtn),frac,AveDist
	integer Nhh,Nho,Nhn,Nhc,Nha,Noo,Non,Noc,Noa
	integer Nnn,Nnc,Nna,Ncc,Nca,Naa
	integer local,Qual(15,20000)
	real*4 distij(15,20000),Adist(15),AveNBI


c get a sense of the average, min and max # of interactions for any atom
	Imax = 0
	Imin = 10000
	Iave = 0
	do i = 1, Nprtn
	 Nint(i) = 0
	 do j = 1, Nprtn
	  if(NBij(i,j).eq.1)Nint(i) = Nint(i) + 1
	 enddo
c	write(51,*)i,Nint(i)
	enddo

	DO i = 1, Nprtn
	 if(Nint(i).gt.Imax)Imax = Nint(i)
	 if(Nint(i).lt.Imin)Imin = Nint(i)
	 Iave = Iave + Nint(i)
	ENDDO
	AveNBI = real(Iave)/real(Nprtn)

c	do i = 1, Nprtn
c	 write(88,*)i,real(Nint(i))/AveNBI
c	enddo
c	close(88)

c	write(*,*)' '
	write(*,*)' Unbound atoms at lines:'
	N0=0
	DO i = 1, Nprtn
	 if(Nint(i).eq.0)then
		write(*,*)asort(i),rlist(i),rsort(i)
		N0=N0+1
	 endif
	ENDDO
	write(*,*)'Total unbound atoms: ',N0

	DO h = 1, Imax+1
	 Hist(h) = 0
	Do i = 1, Nprtn
	 if(Nint(i).eq.h-1)Hist(h) = Hist(h) + 1
	Enddo
	ENDDO

	write(*,*)' '
	write(*,*)' Average # of interaction/atom is: ',
     1real(Iave)/real(Nprtn)
	write(*,*)' '
	write(*,*)' Maximum # of interactions for any atom is: ',Imax
	write(*,*)' '
	write(*,*)' Minimum # of interactions for any atom is: ',Imin
	write(*,*)' '
c	write(*,*)' # of inter-atomic interactions distribution'
c	write(*,*)'   # of interactions,       #'
c	space='     '
c 	do h = 1, Imax+1
c	 write(*,*)h-1,space,Hist(h)
c	enddo
c	write(*,*)' '

	do i = 1, Nprtn
 	 dmin(i) = 10.d0
 	enddo
 	AveDist = 0.0
 	num = 0

 	do i = 1, Nprtn
 	do j = i, Nprtn
	 if(NBij(i,j).ne.0)then
           dist = dsqrt( (xyz(1,i)-xyz(1,j))**2 +
     1                   (xyz(2,i)-xyz(2,j))**2 +
     1                   (xyz(3,i)-xyz(3,j))**2 )
 	   if(dmin(i).gt.dist)dmin(i) = dist
	 endif
 	enddo
 	enddo

 	AveDist = AveDist/dble(num)
 	dave = 0.d0
 	do i = 1, Nprtn
 	 dave = dave +dmin(i)
 	enddo
 	dave = dave / real(Nprtn)

c Make a distance histogram

	Dist_min = 15.0d0
	Dist_max = 0.0d0
	Dist_ave = 0.d0
	HD(1) = 0.d0
	HD(2) = 0.d0
	HD(3) = 0.d0
	HD(4) = 0.d0
	HD(5) = 0.d0
	mcmc = 0
	mcsc = 0
	scsc = 0

	do i = 1, Nprtn
	do j = i, Nprtn
	if(NBij(i,j).ne.0)then

	IF(asort(i).eq.' H   '.or.asort(i).eq.' N   '.or.
     1asort(i).eq.' CA  '.or.asort(i).eq.' C   '.or.
     2asort(i).eq.' O   ')THEN
	  if(asort(j).eq.' H   '.or.asort(j).eq.' N   '
     1.or.asort(j).eq.' CA  '.or.asort(j).eq.' C   '.or.
     2asort(j).eq.' O   ')then
		mcmc= mcmc + 1
	  else
		mcsc=mcsc+1
	  endif
	ELSE
	  if(asort(j).eq.' H   '.or.asort(j).eq.' N   '
     1.or.asort(j).eq.' CA  '.or.asort(j).eq.' C   '.or.
     2asort(j).eq.' O   ')then
		mcsc= mcsc + 1
	  else
	   	scsc = scsc + 1
	  endif
	ENDIF

         dist= dsqrt( (xyz(1,i)-xyz(1,j))**2 +
     1                (xyz(2,i)-xyz(2,j))**2 +
     1                (xyz(3,i)-xyz(3,j))**2 )
	if(dist.eq.0.d0)write(*,*)i,j,asort(i),rsort(i),rlist(i),
     1asort(j),rsort(j),rlist(j),(xyz(l,i),l=1,3),
     2(xyz(l,j),l=1,3)
	 if(dist.lt.Dist_min)Dist_min = dist
	 if(dist.gt.Dist_max)Dist_max = dist
	 Dist_Ave = Dist_Ave + dist
	 if(dist.le.1.d0)HD(1) = HD(1) + 1
	 if(dist.le.2.d0.and.dist.gt.1.d0)HD(2) = HD(2) + 1
	 if(dist.le.3.d0.and.dist.gt.2.d0)HD(3) = HD(3) + 1
	 if(dist.le.4.d0.and.dist.gt.3.d0)HD(4) = HD(4) + 1
	 if(dist.le.5.d0.and.dist.gt.4.d0)HD(5) = HD(5) + 1
	endif
	enddo
	enddo
	write(*,*)' '
	write(*,*)'Total # of interactions (ijss):',ijss
	write(*,*)'# MC-MC interactions:',mcmc
	write(*,*)'# MC-SC interactions:',mcsc
	write(*,*)'# SC-SC interactions:',scsc
	write(*,*)' '

	write(*,*)' Minimum distance of Separation: ',sngl(Dist_min)
	write(*,*)'Average minimum distance of separation 
     1over all atoms:',sngl(dave)
	write(*,*)' '
	write(*,*)' Maximum distance of Separation: ',sngl(Dist_max)
	
	write(*,*)' '
	write(*,*)'Separation(A), #'
	do k = 1,5
	 write(*,*)k-1,'-',k,'(A)',int(HD(k))
	enddo


c Check distribution characteristics of nonbonded interactions
	do i = 1,15
	 Adist(i)=0.0
	enddo
	Nhh = 0
	Nho = 0
	Nhn = 0
	Nhc = 0
	Nha = 0
	Noo = 0
	Non = 0
	Noc = 0
	Noa = 0
	Nnn = 0
	Nnc = 0
	Nna = 0
	Ncc = 0
	Nca = 0
	Naa = 0
	Iplus = 0
	Iminus = 0

	do i = 1, Nprtn
	do j = i, Nprtn
	if(NBij(i,j).ne.0)then

         dist= dsqrt( (xyz(1,i)-xyz(1,j))**2 +
     1                (xyz(2,i)-xyz(2,j))**2 +
     1                (xyz(3,i)-xyz(3,j))**2 )
	 comp = C*dsqrt(Rhc(i)*Rhc(j))
	 if(dist.gt.comp)Iplus=Iplus+1
	 if(dist.le.comp)Iminus=Iminus+1

	endif
	enddo
	enddo

	write(*,*)' '
	frac = Iplus/ijss
	write(*,*)'# NBI > C*sqrt(ri*rj) :',Iplus
	frac = Iminus/ijss
	write(*,*)'# NBI < C*sqrt(ri*rj) :',Iminus
	write(*,*)' '

	do i = 1,5
	do j = 1,5
	 Rbin(i,j)=0
	enddo
	enddo
	do i = 1, Nprtn
	do j = i, Nprtn
	if(NBij(i,j).ne.0)then


	 mcmc = 0
	 local = -1  !-1 means YES local

	 if(abs(rlist(i)-rlist(j)).gt.2)local=1 ! local=+1 means NOT local

	 IF(asort(i).eq.' H   '.or.asort(i).eq.' N   '.or.
     1asort(i).eq.' CA  '.or.asort(i).eq.' C   '.or.
     2asort(i).eq.' O   ')THEN
	   if(asort(j).eq.' H   '.or.asort(j).eq.' N   '
     1.or.asort(j).eq.' CA  '.or.asort(j).eq.' C   '.or.
     2asort(j).eq.' O   ')then
		mcmc = 1
	   else
		mcmc = 2
	   endif
	 ELSE
	   if(asort(j).eq.' H   '.or.asort(j).eq.' N   '
     1.or.asort(j).eq.' CA  '.or.asort(j).eq.' C   '.or.
     2asort(j).eq.' O   ')then
		mcmc = 2
	   else
	   	mcmc = 3
	   endif
	 ENDIF

	 if(Rhc(i).eq.2.8525d0)Icol=1
	 if(Rhc(i).eq.3.1005d0)Icol=2
	 if(Rhc(i).eq.3.8171d0)Icol=3
	 if(Rhc(i).eq.4.3150d0)Icol=4
	 if(Rhc(i).eq.4.2202d0)Icol=5

	 if(Rhc(j).eq.2.8525d0)Irow=1
	 if(Rhc(j).eq.3.1005d0)Irow=2
	 if(Rhc(j).eq.3.8171d0)Irow=3
	 if(Rhc(j).eq.4.3150d0)Irow=4
	 if(Rhc(j).eq.4.2202d0)Irow=5

	 Rbin(Icol,Irow)=Rbin(Icol,Irow)+1

         dist= dsqrt( (xyz(1,i)-xyz(1,j))**2 +
     1                (xyz(2,i)-xyz(2,j))**2 +
     1                (xyz(3,i)-xyz(3,j))**2 )

c       H      Q/V/O       N/M     C/S     A/B/G     
c    H 11(1)    12(2)    13(3)    14(4)    15(5)
c    O na       22(6)    23(7)    24(8)    25(9)
c    N na       na       33(10)   34(11)   35(12)
c    C na       na        na      44(13)   45(14)
c    A na       na        na       na      55(15) 

	 if(Icol.eq.1.and.Irow.eq.1)then
		Nhh=Nhh+1
		distij(1,Nhh)=dist
		Qual(1,Nhh) = local*mcmc
		Adist(1)=Adist(1)+dist
	 endif
	 if((Icol.eq.1.and.Irow.eq.2).or.
     1       (Icol.eq.2.and.Irow.eq.1))then
		Nho=Nho+1
		distij(2,Nho)=dist
		Qual(2,Nho) = local*mcmc
		Adist(2)=Adist(2)+dist
	 endif
	 if((Icol.eq.1.and.Irow.eq.3).or.
     1       (Icol.eq.3.and.Irow.eq.1))then
		Nhn=Nhn+1
		distij(3,Nhn)=dist
		Qual(3,Nhn) = local*mcmc
		Adist(3)=Adist(3)+dist
	 endif
	 if((Icol.eq.1.and.Irow.eq.4).or.
     1       (Icol.eq.4.and.Irow.eq.1))then
		Nhc=Nhc+1
		distij(4,Nhc)=dist
		Qual(4,Nhc) = local*mcmc
		Adist(4)=Adist(4)+dist
	 endif
	 if((Icol.eq.1.and.Irow.eq.5).or.
     1       (Icol.eq.5.and.Irow.eq.1))then
		Nha=Nha+1
		distij(5,Nha)=dist
		Qual(5,Nha) = local*mcmc
		Adist(5)=Adist(5)+dist
	 endif
	 if(Icol.eq.2.and.Irow.eq.2)then
		Noo=Noo+1
		distij(6,Noo)=dist
		Qual(6,Noo) = local*mcmc
		Adist(6)=Adist(6)+dist
	 endif
	 if((Icol.eq.2.and.Irow.eq.3).or.
     1       (Icol.eq.3.and.Irow.eq.2))then
		Non=Non+1
		distij(7,Non)=dist
		Qual(7,Non) = local*mcmc
		Adist(7)=Adist(7)+dist
	 endif
	 if((Icol.eq.2.and.Irow.eq.4).or.
     1       (Icol.eq.4.and.Irow.eq.2))then
		Noc=Noc+1
		distij(8,Noc)=dist
		Qual(8,Noc) = local*mcmc
		Adist(8)=Adist(8)+dist
	 endif
	 if((Icol.eq.2.and.Irow.eq.5).or.
     1       (Icol.eq.5.and.Irow.eq.2))then
		Noa=Noa+1
		distij(9,Noa)=dist
		Qual(9,Noa) = local*mcmc
		Adist(9)=Adist(9)+dist
	 endif
	 if(Icol.eq.3.and.Irow.eq.3)then
		Nnn=Nnn+1
		distij(10,Nnn)=dist
		Qual(10,Nnn) = local*mcmc
		Adist(10)=Adist(10)+dist
	 endif
	 if((Icol.eq.3.and.Irow.eq.4).or.
     1       (Icol.eq.4.and.Irow.eq.3))then
		Nnc=Nnc+1
		distij(11,Nnc)=dist
		Adist(11)=Adist(11)+dist
		Qual(11,Nnc) = local*mcmc
	 endif
	 if((Icol.eq.3.and.Irow.eq.5).or.
     1       (Icol.eq.5.and.Irow.eq.3))then
		Nna=Nna+1
		distij(12,Nna)=dist
		Adist(12)=Adist(12)+dist
		Qual(12,Nna) = local*mcmc
	 endif
	 if(Icol.eq.4.and.Irow.eq.4)then
		Ncc=Ncc+1
		distij(13,Ncc)=dist
		Adist(13)=Adist(13)+dist
		Qual(13,Ncc) = local*mcmc
	 endif
	 if((Icol.eq.4.and.Irow.eq.5).or.
     1       (Icol.eq.5.and.Irow.eq.4))then
		Nca=Nca+1
		distij(14,Nca)=dist
		Adist(14)=Adist(14)+dist
		Qual(14,Nca) = local*mcmc
	 endif
	 if(Icol.eq.5.and.Irow.eq.5)then
		Naa=Naa+1
		distij(15,Naa)=dist
		Adist(15)=Adist(15)+dist
		Qual(15,Naa) = local*mcmc
	 endif

	endif
	enddo
	enddo
	
	do i = 1,5
	do j = i+1,5
	 Rbin(i,j) = Rbin(i,j) + Rbin(j,i)
	enddo
	enddo
	do i =1,5
	do j =1,5
	if(j.lt.i)Rbin(i,j)=0
	enddo
	enddo


	write(*,*)' '
	write(*,*)' Distribution of TYPES of interaction:'
	write(*,*)'       H       O       N       C       A'
	Iplus = 0
	do i = 1,5
	 write(*,*)(Rbin(i,j),j=1,5)
	do j = 1,5
	 Iplus = Iplus + Rbin(i,j)
	enddo
	enddo

	write(*,*)'Total # NBI: ',Iplus

	write(*,*)
     1'Average distances of separation for 15 types of NBI:'
	k = 0
	do i=1,5
	do j=i,5
	 k = k+1
	 Adist(k) = Adist(k)/real(Rbin(i,j))
	 write(*,*)k,Adist(k)
	enddo
	enddo


c	open(21)
c	open(210)
cc	write(21,*)'H H interaction distances:',2.8525
c	do i = 1, Nhh
c	 write(21,*)i,distij(1,i)
c	 write(210,*)i,Qual(1,i)
c	enddo
c	close(21)
c	close(210)
c
c	open(22)
c	open(220)
c	comp = sqrt(2.8525*3.7005)
cc	write(22,*)'H O interaction distances:',comp
c	do i = 1, Nho
c	 write(22,*)i,distij(2,i)
c	 write(220,*)i,Qual(2,i)
c	enddo
c	close(22)
c	close(220)
c
c	open(23)
c	open(230)
c	comp = sqrt(2.8525*3.8171)
cc	write(23,*)'H N interaction distances:',comp
c	do i = 1, Nhn
c	 write(23,*)i,distij(3,i)
c	 write(230,*)i,Qual(3,i)
c	enddo
c	close(23)
c	close(230)
c
c	open(24)
c	open(240)
c	comp = sqrt(2.8525*4.315)
cc	write(24,*)'H C interaction distances:',comp
c	do i = 1, Nhc
c	 write(24,*)i,distij(4,i)
c	 write(240,*)i,Qual(4,i)
c	enddo
c	close(24)
c	close(240)
c
c	open(25)
cc	open(250)
c	comp = sqrt(2.8525*4.2202)
cc	write(25,*)'H A interaction distances:',comp
c	do i = 1, Nha
c	 write(25,*)i,distij(5,i)
c	 write(250,*)i,Qual(5,i)
c	enddo
c	close(25)
c	close(250)
c
c	open(26)
c	open(260)
cc	write(26,*)'O O interaction distances:',3.7005
cc	do i = 1, Noo
c	 write(26,*)i,distij(6,i)
c	 write(260,*)i,Qual(6,i)
c	enddo
c	close(26)
c	close(260)
c
c	open(27)
c	open(270)
c	comp = sqrt(3.7005*3.8171)
cc	write(27,*)'O N interaction distances:',comp
c	do i = 1, Non
c	 write(27,*)i,distij(7,i)
c	 write(270,*)i,Qual(7,i)
c	enddo
c	close(27)
c	close(270)
c
c	open(28)
c	open(280)
c	comp = sqrt(3.7005*4.315)
cc	write(28,*)'O C interaction distances:',comp
c	do i = 1, Noc
c	 write(28,*)i,distij(8,i)
c	 write(280,*)i,Qual(8,i)
c	enddo
c	close(28)
c	close(280)
c
c	open(29)
c	open(290)
c	comp = sqrt(3.7005*4.2202)
cc	write(29,*)'O A interaction distances:',comp
c	do i = 1, Noa
c	 write(29,*)i,distij(9,i)
c	 write(290,*)i,Qual(9,i)
c	enddo
c	close(29)
c	close(290)
c
c	open(30)
c	open(300)
cc	write(30,*)'N N interaction distances:',3.8171
c	do i = 1, Nnn
c	 write(30,*)i,distij(10,i)
c	 write(300,*)i,Qual(10,i)
c	enddo
c	close(30)
c	close(300)
c
cc	open(31)
c	open(310)
c	comp = sqrt(3.8171*4.315)
cc	write(31,*)'N C interaction distances:',comp
c	do i = 1, Nnc
c	 write(31,*)i,distij(11,i)
c	 write(310,*)i,Qual(11,i)
c	enddo
c	close(31)
c	close(310)
c
c	open(32)
c	open(320)
c	comp = sqrt(3.8171*4.2202)
cc	write(32,*)'N A interaction distances:',comp
c	do i = 1, Nna
c	 write(32,*)i,distij(12,i)
cc	 write(320,*)i,Qual(12,i)
c	enddo
c	close(32)
c	close(320)
c
c	open(33)
cc	open(330)
cc	write(33,*)'C C interaction distances:',4.315
c	do i = 1, Ncc
c	 write(33,*)i,distij(13,i)
c	 write(330,*)i,Qual(13,i)
c	enddo
c	close(33)
c	close(330)
c
c	open(34)
c	open(340)
c	comp = sqrt(4.315*4.2202)
ccc	write(34,*)i,'C A interaction distances:',comp
c	do i = 1, Nca
c	 write(34,*)distij(14,i)
c	 write(340,*)i,Qual(14,i)
c	enddo
c	close(34)
c	close(340)
c
c	open(35)
c	open(350)
cc	write(35,*)i,'A A interaction distances:',4.2202
c	do i = 1, Naa
c	 write(35,*)distij(15,i)
c	 write(350,*)i,Qual(15,i)
c	enddo
c	close(35)
c	close(350)
c

	return
	end
