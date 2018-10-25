	subroutine Pvec(Nprtn,Nres,Nxi,IndexMC,IndexXi,xyz,P)

	implicit none
	integer Nprtn,Nres,Nxi
	integer i1,i2,i3,i4,itor
	integer IndexMC(2*Nres,4),IndexXi(Nxi,4)
	real*8 angle
	real*8 xyz(3,Nprtn),P(2*Nres+Nxi)


	P(1) = 0.d0 !Nterminal phi not defined
	do itor = 2, 2*Nres
	  i1 = IndexMC(itor,1)
	  i2 = IndexMC(itor,2)
	  i3 = IndexMC(itor,3)
	  i4 = IndexMC(itor,4)
	  call angletor(i1,i2,i3,i4,Nprtn,xyz,angle)
	  P(itor) = angle
	enddo

	do itor = 1, Nxi
         i1 = IndexXi(itor,1)
         i2 = IndexXi(itor,2)
         i3 = IndexXi(itor,3)
         i4 = IndexXi(itor,4)
         call angletor(i1,i2,i3,i4,Nprtn,xyz,angle)
	 P(2*Nres+itor) = angle
	enddo


        return
        end


