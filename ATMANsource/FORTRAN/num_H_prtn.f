        subroutine num_H_prtn(massP,xyzP,IndexMC,IndexXi,EndXi,Nres,
     1	Nxi,Natom,xyzL,massL,Nlig,NligNA,Lig,Nder,H,delta,dRdTor,NDoF,
     2  Pyn)


c This program calculates the Hij=dr/dq(i)*dr/dq(j) matrix, computing
c moving derivatives to keep the CoM and rotation about CoM fixed.
c The derivatives are obtained numerically using W. Kabsch's program
c u3best to eliminate rotations/translations about the CoM.
c The numeric version is much faster than the analytic version, and
c requires large memory capacity (on order of Natom*(2*Nres+Nxi)*3 )
c The degrees of freedom (DoF) include all torsions within a protein.

c On 29 Aug. 1995, this routine was (re)checked against the analytic
c version, ana_H_prtn.f, using ken.pdb by MMT.  

c on 9 Sept. 2014 modified to include Nlig single atom ligands w/3 RB DoFs each
c on 11 Feb. 2015 modifief to include Nlig ligands w/ NligNA atoms total
c on 25 Feb, 2015 modified to include RB rotational DoF

	implicit none
c input:
	integer Natom,NDoF,Nres,Nxi,Nder,Nlig,NligNA,iLig
	integer IndexMC(2*Nres,4),IndexXi(Nxi,4),EndXi(Nxi)
        real*8 xyzP(3,Natom),massP(Natom),delta
	real*8 mass(Natom+NligNA),xyz(3,Natom+NligNA)
	real*8 xyz0(3,Natom+NligNA)
	real*8 xyzL(3,2500),massL(2500)
	integer Lig(1000,2),OnOff(2500),Pyn(2*Nres+Nxi)

c for computations:
	integer i,icount,itor,jtor,iatom,m,m2,m3,ier
	integer k,ks,ls,icheck
	real*8 Tmass,Rcm(3),rms,U(3,3),T(3),x,y,z,sum,dist
	real*8 TmassL(Nlig),Rcm_lig(3,Nlig)
	real*8 nhat(3,2*Nres+Nxi)
	real*8 dRdTor(Natom+NligNA,NDoF,3)
	real*8 xyzdP(3,Natom),xyzdL(3,2500)
c output:
	real*8 H(Nder)


c  Get Center of Mass of each Ligand (for RB rotations)
        i = 0
        DO k = 1, Nlig
	IF(Lig(k,1).gt.1)THEN
         TmassL(k) = 0.d0
         Rcm_lig(1,k) = 0.d0
         Rcm_lig(2,k) = 0.d0
         Rcm_lig(3,k) = 0.d0
        do ls = 1, Lig(k,1)
          i = i + 1
          TmassL(k) = TmassL(k) + massL(i)
          Rcm_lig(1,k) = Rcm_lig(1,k) + massL(i)*xyzL(1,i)
          Rcm_lig(2,k) = Rcm_lig(2,k) + massL(i)*xyzL(2,i)
          Rcm_lig(3,k) = Rcm_lig(3,k) + massL(i)*xyzL(3,i)
        enddo
         Rcm_lig(1,k) = Rcm_lig(1,k)/TmassL(k)
         Rcm_lig(2,k) = Rcm_lig(2,k)/TmassL(k)
         Rcm_lig(3,k) = Rcm_lig(3,k)/TmassL(k)
	ELSE
	 i = i + 1
         Rcm_lig(1,k) = xyzL(1,i)
         Rcm_lig(2,k) = xyzL(2,i)
         Rcm_lig(3,k) = xyzL(3,i)
	ENDIF
        ENDDO


c create super-arrays of both protein as well as ligand
	do iatom = 1,Natom
	 xyz0(1,iatom) = xyzP(1,iatom)
	 xyz0(2,iatom) = xyzP(2,iatom)
	 xyz0(3,iatom) = xyzP(3,iatom)
	 mass(iatom) = massP(iatom)
	enddo
	do iatom = Natom+1,Natom+NligNA
	 xyz0(1,iatom)=xyzL(1,iatom-Natom)
	 xyz0(2,iatom)=xyzL(2,iatom-Natom)
	 xyz0(3,iatom)=xyzL(3,iatom-Natom)
	 mass(iatom) = massL(iatom - Natom)
	enddo


c    calculate the center of mass:
        Rcm(1) = 0.d00
        Rcm(2) = 0.d00
        Rcm(3) = 0.d00
        Tmass = 0.d00

        do iatom = 1, Natom+NligNA
         Rcm(1) = Rcm(1) + mass(iatom)*xyz0(1,iatom) 
         Rcm(2) = Rcm(2) + mass(iatom)*xyz0(2,iatom) 
         Rcm(3) = Rcm(3) + mass(iatom)*xyz0(3,iatom) 
         Tmass = Tmass + mass(iatom)
         mass(iatom) = dsqrt(mass(iatom)) !!!!!!
        enddo

        Rcm(1) = Rcm(1)/Tmass
        Rcm(2) = Rcm(2)/Tmass
        Rcm(3) = Rcm(3)/Tmass

c Move coordinates to CM:
        do iatom = 1,Natom+NligNA
         xyz0(1,iatom) = xyz0(1,iatom) - Rcm(1)
         xyz0(2,iatom) = xyz0(2,iatom) - Rcm(2)
         xyz0(3,iatom) = xyz0(3,iatom) - Rcm(3)
        enddo   ! end iatom-loop
	do iatom = 1, Natom ! subset for updateMC/updateSC
	 xyzP(1,iatom) = xyz0(1,iatom)
	 xyzP(2,iatom) = xyz0(2,iatom)
	 xyzP(3,iatom) = xyz0(3,iatom)
	enddo
	do iatom = 1, NligNA  ! for updateRB
	 xyzL(1,iatom) = xyz0(1,iatom+Natom)
	 xyzL(2,iatom) = xyz0(2,iatom+Natom)
	 xyzL(3,iatom) = xyz0(3,iatom+Natom)
	enddo


c  get nhat vector for torsional updates:
	do m = 1,2*Nres + Nxi
	 nhat(1,m) = 0.0d00
	 nhat(2,m) = 0.0d00
	 nhat(3,m) = 0.0d00
	enddo

	do m = 1,2*Nres     
	 if(Pyn(m).eq.0)goto10
	 m3 = IndexMC(m,3)
	 m2 = IndexMC(m,2)
	 T(1) = xyzP(1,m3) - xyzP(1,m2)
	 T(2) = xyzP(2,m3) - xyzP(2,m2)
	 T(3) = xyzP(3,m3) - xyzP(3,m2)
	 dist = dsqrt( T(1)**2 + T(2)**2 + T(3)**2 )
	 nhat(1,m) = T(1)/dist
	 nhat(2,m) = T(2)/dist
	 nhat(3,m) = T(3)/dist
10	 continue
	enddo


	do m = 1,Nxi
	 m3 = IndexXi(m,3)
	 m2 = IndexXi(m,2)
	 T(1) = xyzP(1,m3) - xyzP(1,m2)
	 T(2) = xyzP(2,m3) - xyzP(2,m2)
	 T(3) = xyzP(3,m3) - xyzP(3,m2)
	 dist = dsqrt( T(1)**2 + T(2)**2 + T(3)**2 )
	 nhat(1,2*Nres+m) = T(1)/dist
	 nhat(2,2*Nres+m) = T(2)/dist
	 nhat(3,2*Nres+m) = T(3)/dist
	enddo


c Initialize:
	do itor = 1,NDoF  
        do iatom = 1,Natom+NligNA
         dRdTor(iatom,itor,1) = 0.d00
         dRdTor(iatom,itor,2) = 0.d00
         dRdTor(iatom,itor,3) = 0.d00
	enddo
	enddo

c--------------------------------------------------------------------
c Set up dRdTor(# of atoms,# of D.o.F.s,3) vector:


	icheck = 0
        DO itor = 1,2*Nres
	icheck = icheck + 1
	 if(Pyn(itor).eq.0)goto20
         call updateMC(xyzP,xyzdP,IndexMC,nhat,itor,Nres,Nxi,
     1			Natom,delta)
c recombine the coordinates
	do iatom = 1, Natom
	 xyz(1,iatom) = xyzdP(1,iatom)
	 xyz(2,iatom) = xyzdP(2,iatom)
	 xyz(3,iatom) = xyzdP(3,iatom)
	enddo
	do iatom = 1, NligNA
	 xyz(1,iatom+Natom) = xyz0(1,iatom+Natom)
	 xyz(2,iatom+Natom) = xyz0(2,iatom+Natom)
	 xyz(3,iatom+Natom) = xyz0(3,iatom+Natom)
	enddo

         call u3best(mass,xyz,xyz0,Natom+NligNA,1,rms,U,T,ier)  

         do iatom = 1,Natom+NligNA
          x = xyz(1,iatom)
          y = xyz(2,iatom)
          z = xyz(3,iatom)
          xyz(1,iatom) = T(1) + U(1,1)*x + U(1,2)*y + U(1,3)*z
          xyz(2,iatom) = T(2) + U(2,1)*x + U(2,2)*y + U(2,3)*z
          xyz(3,iatom) = T(3) + U(3,1)*x + U(3,2)*y + U(3,3)*z
         enddo

        do iatom=1,Natom+NligNA
         dRdTor(iatom,itor,1) = ( xyz(1,iatom)-xyz0(1,iatom) )/delta
         dRdTor(iatom,itor,2) = ( xyz(2,iatom)-xyz0(2,iatom) )/delta
         dRdTor(iatom,itor,3) = ( xyz(3,iatom)-xyz0(3,iatom) )/delta
        enddo   ! end iatom-loop
20	continue
        ENDDO   ! end itor loop



        DO itor = 1,Nxi
	icheck = icheck + 1
         call updateSC(xyzP,xyzdP,IndexXi,EndXi,nhat,itor,Nres,
     1    Nxi,Natom,delta)
c 	recombine the coordinates
	do iatom = 1, Natom
	 xyz(1,iatom) = xyzdP(1,iatom)
	 xyz(2,iatom) = xyzdP(2,iatom)
	 xyz(3,iatom) = xyzdP(3,iatom)
	enddo
	do iatom = 1, NligNA
	 xyz(1,iatom+Natom) = xyz0(1,iatom+Natom)
	 xyz(2,iatom+Natom) = xyz0(2,iatom+Natom)
	 xyz(3,iatom+Natom) = xyz0(3,iatom+Natom)
	enddo

         call u3best(mass,xyz,xyz0,Natom+NligNA,1,rms,U,T,ier)  

         do iatom=1,Natom+NligNA
          x = xyz(1,iatom)
          y = xyz(2,iatom)
          z = xyz(3,iatom)
          xyz(1,iatom) = T(1) + U(1,1)*x + U(1,2)*y + U(1,3)*z
          xyz(2,iatom) = T(2) + U(2,1)*x + U(2,2)*y + U(2,3)*z
          xyz(3,iatom) = T(3) + U(3,1)*x + U(3,2)*y + U(3,3)*z
         enddo

         do iatom = 1,Natom+NligNA
          dRdTor(iatom,2*Nres+itor,1) = ( xyz(1,iatom) - 
     1					   xyz0(1,iatom) )/delta
          dRdTor(iatom,2*Nres+itor,2) = ( xyz(2,iatom) - 
     1	 				   xyz0(2,iatom) )/delta
          dRdTor(iatom,2*Nres+itor,3) = ( xyz(3,iatom) - 
     1					   xyz0(3,iatom) )/delta
         enddo   ! end iatom-loop
        ENDDO   ! end itor loop


	itor = 0
	DO iLig = 1, Nlig

c select subset of Ligand array to activate by first initializing OnOff to zero
        do ks = 1, NligNA
         OnOff(ks) = 0
        enddo
c and then activating each ligand subset by setting appropriate OnOff entries to 1
        icount = 0
        do ks = 1,Nlig
        do ls = 1,lig(ks,1)
         icount = icount + 1
         if(ks.eq.iLig)OnOff(icount) = 1
        enddo
        enddo


	ks = 3
	if(lig(iLig,1).gt.1)ks=6
	do i = 1,ks
	icheck = icheck + 1
	itor = itor + 1

	call updateRB(xyzL,xyzdL,i,OnOff,Nlig,NligNA,
     1Rcm_lig,iLig,delta)

c recombine coords
	do iatom = 1, Natom
	 xyz(1,iatom) = xyz0(1,iatom)
	 xyz(2,iatom) = xyz0(2,iatom)
	 xyz(3,iatom) = xyz0(3,iatom)
	enddo
	do iatom = 1,NligNA
	 xyz(1,iatom+Natom) = xyzdL(1,iatom)
	 xyz(2,iatom+Natom) = xyzdL(2,iatom)
	 xyz(3,iatom+Natom) = xyzdL(3,iatom)
	enddo

         call u3best(mass,xyz,xyz0,Natom+NligNA,1,rms,U,T,ier)  

         do iatom=1,Natom+NligNA
          x = xyz(1,iatom)
          y = xyz(2,iatom)
          z = xyz(3,iatom)
          xyz(1,iatom) = T(1) + U(1,1)*x + U(1,2)*y + U(1,3)*z
          xyz(2,iatom) = T(2) + U(2,1)*x + U(2,2)*y + U(2,3)*z
          xyz(3,iatom) = T(3) + U(3,1)*x + U(3,2)*y + U(3,3)*z
         enddo

         do iatom = 1,Natom+NligNA
          dRdTor(iatom,2*Nres+Nxi+itor,1) = ( xyz(1,iatom) - 
     1					   xyz0(1,iatom) )/delta
          dRdTor(iatom,2*Nres+Nxi+itor,2) = ( xyz(2,iatom) - 
     1	 				   xyz0(2,iatom) )/delta
          dRdTor(iatom,2*Nres+Nxi+itor,3) = ( xyz(3,iatom) - 
     1					   xyz0(3,iatom) )/delta
         enddo   ! end iatom-loop

	enddo  ! end RB DoF (1-6)
	ENDDO  ! loop over all ligands (iLig)

	write(*,*)'icheck equals NDoF?',icheck,NDoF
	if(icheck.ne.NDoF)stop
c-------------------------------------------------------------------------


	icount = 0

        DO itor = 1,NDoF
        DO jtor = itor,NDoF
         sum = 0.d00

         do iatom = 1,Natom+NligNA
          sum = sum + (dRdTor(iatom,itor,1) * dRdTor(iatom,jtor,1)
     1              +  dRdTor(iatom,itor,2) * dRdTor(iatom,jtor,2)
     2              +  dRdTor(iatom,itor,3) * dRdTor(iatom,jtor,3))
     3              *  mass(iatom)*mass(iatom)
         enddo   ! iatom loop

	 icount = icount + 1
         H(icount) = sum

        ENDDO   ! jtor loop
        ENDDO   ! itor loop


        return
        end
