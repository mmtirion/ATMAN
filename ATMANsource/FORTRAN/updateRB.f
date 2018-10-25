	subroutine updateRB(rLig0,rLig,n,OnOff,Nlig,NligNA,
     1Rcm_lig,k,delta)

	implicit none

	integer i,j,k,n
	integer Nlig,NligNA,OnOff(2500)
	real*8 rLig0(3,2500),rLig(3,2500)
	real*8 dx,dy,dz,delta
	real*8 x0,y0,z0
	real*8 Rcm_lig(3,1000)



c initiliaze rLig:
	do i = 1, NligNA
	 rLig(1,i) = rLig0(1,i)
	 rLig(2,i) = rLig0(2,i)
	 rLig(3,i) = rLig0(3,i)
	enddo


	dx = 0.d0
	dy = 0.d0
	dz = 0.d0

	if(n.eq.1)dx=delta
	if(n.eq.2)dy=delta
	if(n.eq.3)dz=delta

	do i = 1, NligNA
	 if(OnOff(i).eq.1)rLig(1,i) = rLig0(1,i) + dx
	 if(OnOff(i).eq.1)rLig(2,i) = rLig0(2,i) + dy
	 if(OnOff(i).eq.1)rLig(3,i) = rLig0(3,i) + dz
	enddo


        if(n.eq.4)then
        do i = 1, NligNA
	 if(OnOff(i).eq.1)then
          y0 = rLig0(2,i) - Rcm_lig(2,k)
          z0 = rLig0(3,i) - Rcm_lig(3,k)
          rLig(2,i) = y0*dcos(delta)+z0*dsin(delta) + Rcm_lig(2,k)
          rLig(3,i) =-y0*dsin(delta)+z0*dcos(delta) + Rcm_lig(3,k)
	 endif
        enddo
        endif

        if(n.eq.5)then
        do i = 1, NligNA
	 if(OnOff(i).eq.1)then
          x0 = rLig0(1,i) - Rcm_lig(1,k)
          z0 = rLig0(3,i) - Rcm_lig(3,k)
          rLig(1,i) = x0*dcos(delta)+z0*dsin(delta) + Rcm_lig(1,k)
          rLig(3,i) =-x0*dsin(delta)+z0*dcos(delta) + Rcm_lig(3,k)
	 endif
        enddo
        endif

        if(n.eq.6)then
        do i = 1, NligNA
	 if(OnOff(i).eq.1)then
          x0 = rLig0(1,i) - Rcm_lig(1,k)
          y0 = rLig0(2,i) - Rcm_lig(2,k)
          rLig(1,i) = x0*dcos(delta)+y0*dsin(delta) + Rcm_lig(1,k)
          rLig(2,i) =-x0*dsin(delta)+y0*dcos(delta) + Rcm_lig(2,k)
	 endif
        enddo
        endif


	return
	end

