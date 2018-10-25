	subroutine Inter_Prtn_Lig(xyzP,xyzL,Nlig,NligNA,Nprtn,Ninter,
     1NinterL,NBIPL,NBILL,rlist,vdWP,vdWlig,vdWcut)

	implicit none

	integer Nprtn,Nlig,NligNA,Ninter,NinterL
	integer NBIPL(10000,2)
	integer NBILL(2500,2500),rlist(Nprtn)
	integer i,j,k,l,iLig,Nint
	integer is,ks,ls

	real*8 xyzP(3,Nprtn),vdWP(Nprtn)
	real*8 xyzL(3,2500),vdWlig(2500)
	real*8 vdWcut,dist,CutOff
	real*8 xyzLr(3,2500),vdWligNAr(2500)


	write(*,*)' '
	write(*,*)' (Re)start inter-protein-ligand search...'
	write(*,*)' '

	Nint = 0
c note this algorithm ignores ligand-ligand interactions at this point...
	do i = 1,NligNA
	 iLig = 0
	do j = 1,Nprtn

c the value before next vdWcut may need increasing, due to inter-protein water modulation?
c the value needs to match the one 41 lines further down
	 CutOff = 1.10868341797d0*dsqrt(vdWlig(i)*vdWP(j)) + vdWcut !included on 3/7/16

	 dist = dsqrt((xyzL(1,i)-xyzP(1,j))**2 +
     1		      (xyzL(2,i)-xyzP(2,j))**2 +
     2		      (xyzL(3,i)-xyzP(3,j))**2 )

	 if(dist.le.CutOff)then
		iLig = iLig + 1
		Nint=Nint+1
		NBIPL(Nint,1)=i
		NBIPL(Nint,2)=j
	 endif

	enddo
c	 write(*,*)'Ligand atom #',i,'has ',iLig,' NBI w/protein'
	enddo

	write(*,*)' '
	write(*,*)'# of inter-Protein-Ligand interactions:',Nint
	write(*,*)' '

c include a 2nd search inter-ligand/ligand interactions, such as 
c btw Ca++ and ATP in g-Actin and add these into InterLigs 
c as of 1/17, this feature is not incorporated: no Hessian block
c includes inter-ligand NBI.

	if(Nlig.le.1)return
	iLig = 0  ! counter for inter ligand-ligand NBI


	 do k = 1, NligNA  
	  NBILL(k,k) = 0
	 do l = k+1, NligNA
	  NBILL(l,k) = 0
	  NBILL(k,l) = 0

	  CutOff = 1.10868341797d0*dsqrt(vdWlig(k)*vdWlig(l)) + vdWcut !included on 3/7/16
	
	  dist = dsqrt( (xyzL(1,k)-xyzL(1,l))**2 +
     1		        (xyzL(2,k)-xyzL(2,l))**2 +
     2		        (xyzL(3,k)-xyzL(3,l))**2 )

	  if(dist.le.CutOff)then
		iLig = iLig + 1
		NBILL(l,k) = 1
		NBILL(k,l) = 1
	  endif
	 enddo
	 enddo


	write(*,*)' '
	write(*,*)'# of inter ligand-ligand NBI:',iLig
	write(*,*)' '
	if(ilig.gt.NinterL)then
		write(*,*)'# ligand-ligand NBI>NinterL'
		stop
	endif

	return
	end
