        subroutine Nhat(Nres,Nxi,Nprtn,IndexMC,IndexXi,xyz,hatVec)

	integer Nres,Nxi,Nprtn
	integer IndexMC(2*Nres,4),IndexXI(Nxi,4)
	real*8 xyz(3,Nprtn),hatVec(3,2*Nres+Nxi)

	integer m,m2,m3
	real*8 T(3),dist


        do m = 1,2*Nres+Nxi
         hatVec(1,m) = 0.0d00
         hatVec(2,m) = 0.0d00
         hatVec(3,m) = 0.0d00
        enddo



        do m = 1,2*Nres
         m3 = IndexMC(m,3)
         m2 = IndexMC(m,2)
         if(m3.eq.0)goto10
         T(1) = xyz(1,m3) - xyz(1,m2)
         T(2) = xyz(2,m3) - xyz(2,m2)
         T(3) = xyz(3,m3) - xyz(3,m2)
         dist = dsqrt( T(1)**2 + T(2)**2 + T(3)**2 )
         hatVec(1,m) = T(1)/dist
         hatVec(2,m) = T(2)/dist
         hatVec(3,m) = T(3)/dist
10       continue
        enddo

        do m = 1,Nxi
         m3 = IndexXi(m,3)
         m2 = IndexXi(m,2)
         T(1) = xyz(1,m3) - xyz(1,m2)
         T(2) = xyz(2,m3) - xyz(2,m2)
         T(3) = xyz(3,m3) - xyz(3,m2)
         dist = dsqrt( T(1)**2 + T(2)**2 + T(3)**2 )
         hatVec(1,2*Nres+m) = T(1)/dist
         hatVec(2,2*Nres+m) = T(2)/dist
         hatVec(3,2*Nres+m) = T(3)/dist
        enddo

	return
	end


