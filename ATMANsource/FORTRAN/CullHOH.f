	subroutine CullHOH(xyzP,xyzL,NligNA,Nprtn,Lig,
     1rlist,vdWp,vdWlig,NligNAr,Nligr,HOH,Bhet,LigNBI)

	implicit none

	integer Nprtn,Nlig,NligNA
	integer Lig(1000,2)
	integer InterPairs(10000,2)
	integer rlist(Nprtn)
	integer i,j,k,l,iLig,Nint,LigNBIr(2500),LigNBI(2500)
	integer is,ks,ls,rmax

	real*8 xyzP(3,Nprtn),vdWP(Nprtn),Bhet(2500)
	real*8 xyzL(3,2500),vdWlig(2500),Bcull(2500)
	real*8 dist,CutOff,ave

	integer HOH,nonHOH,NligNAr,Nligr
	integer check,ki,ko,res1,res2
	real*8 xyzLr(3,2500),NBI(500,25)


	write(*,*)' '
	write(*,*)' Begin inter-protein-water search...'
	write(*,*)' '


c HOH ligands start after any other ligands:
	nonHOH = 0
	if(HOH.gt.1)then
	do j = 1,HOH-1
	  nonHOH = nonHOH + Lig(j,1)
	enddo
	write(*,*)'First HOH molecule at ligand atom #:',nonHOH+1
	endif

	Nint = 0
	DO i = 1,NligNA

	iLig = 0
	Do j = 1,Nprtn

	 CutOff = 1.10868341797d0*dsqrt(vdWlig(i)*vdWP(j)) !included on 3/7/16

	 dist = dsqrt((xyzL(1,i)-xyzP(1,j))**2 +
     1		      (xyzL(2,i)-xyzP(2,j))**2 +
     2		      (xyzL(3,i)-xyzP(3,j))**2 )

	 if(dist.le.CutOff)then
		iLig = iLig + 1
		Nint = Nint + 1
		NBI(i,iLig)=dist
		InterPairs(Nint,1)=i
		InterPairs(Nint,2)=j
	 endif

	EndDo ! j Nprtn loop

	LigNBI(i) = iLig
c	if(i.gt.nonHOH.and.iLig.eq.0)write(*,*)
c     1'Water ',i-HOH+1,'has zero interactions and is excluded'
	if(i.gt.nonHOH.and.iLig.le.18)then
	  iLig = 0
	  LigNBI(i) = iLig
c	  write(*,*)'Water ',i-nonHOH,'has =< 18 NBI and is excluded'
	 endif
	IF(i.gt.nonHOH.and.iLig.ne.0)THEN
	 write(*,*)'Water ',i-nonHOH,'has ',iLig,' NBI w/protein...'
	 check = 0
	 rmax = 0
	 do ki = Nint-iLig+1, Nint
	 do ko = ki, Nint
	  res1 = rlist(InterPairs(ki,2))
	  res2 = rlist(InterPairs(ko,2))
c	 write(*,*)'...interactions ',abs(res1-res2),'apart'
	  if(abs(res1-res2).gt.rmax)rmax=abs(res1-res2)
	  if(abs(res1-res2).gt.6)check=1 !should this be .gt.4, for the alpha helices?
	 enddo
	 enddo
	 ave = 0.d0
	 do k = 1, iLig
	  ave = ave + NBI(i,k)
	 enddo
	 write(*,*)'Average prtn-lig separation:',ave/dble(iLig)
	 write(*,*)'and has interactions with residues',rmax,'apart'
	 if(check.eq.0)LigNBI(i) = 0
	 if(ligNBI(i).eq.0)write(*,*)
     1' but is bound within 6 neighboring residues, excluded '
	 ENDIF
	ENDDO  ! i ligand loop


c collapse relevant HOHs out of coordinate, vdW arrays:
	NligNAr = 0
	do i = 1, nonHOH
	 NligNAr = NligNAr + 1
	 xyzLr(1,NligNAr) = xyzL(1,i)
	 xyzLr(2,NligNAr) = xyzL(2,i)
	 xyzLr(3,NligNAr) = xyzL(3,i)
	 Bcull(NligNAr) = Bhet(i)
	 LigNBIr(NligNAr) = LigNBI(i)
	enddo
	do i = nonHOH+1, NligNA
	 if(LigNBI(i).ne.0)then
	  NligNAr = NligNAr + 1
	  xyzLr(1,NligNAr) = xyzL(1,i)
	  xyzLr(2,NligNAr) = xyzL(2,i)
	  xyzLr(3,NligNAr) = xyzL(3,i)
	  Bcull(NligNAr) = Bhet(i)
	  LigNBIr(NligNAr) = LigNBI(i)
	 endif
	enddo

c note that since this is about HOH, neither aHET, typeHET, mHET, vdWHET,
c or epsHET need be updated: their values remain the same before/after collapse

c report and return updated coordinate files:
	write(*,*)' '
	check = 0
	do i = 1, HOH-1
	 check = check + Lig(i,1)
	enddo
	write(*,*)'PDB # of HOH: ',NligNA-check
	write(*,*)' Reduced included waters:',NligNAr-check
	Nligr = NligNAr - check + HOH - 1
	write(*,*)'Nligr =',Nligr

	do i = 1, NligNAr
	 xyzL(1,i) = xyzLr(1,i)
	 xyzL(2,i) = xyzLr(2,i)
	 xyzL(3,i) = xyzLr(3,i)
	 Bhet(i) = Bcull(i)
	 LigNBI(i) = LigNBIr(i)
	enddo
	do i = NligNAr+1,NligNA  ! this will create problems. Must never access these.
	 xyzL(1,i) = 0.d0	 
	 xyzL(2,i) = 0.d0	 
	 xyzL(3,i) = 0.d0	 
	 Bhet(i) = 0.d0
	 LigNBIr(i) = 0
	enddo



	return
	end
