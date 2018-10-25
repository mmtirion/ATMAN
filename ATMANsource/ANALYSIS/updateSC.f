        subroutine updateSC(xyz,xyznew,IndexXi,EndXi,nhat,Ixi,Nres,
     1   Nxi,Natom,deldeg)

	IMPLICIT NONE
	integer*4 Natom,Nxi,Nres

	integer*4 i,j3,Ixi,Ientry
	integer*4 IndexXi(Nxi,4),EndXi(Nxi)
	real*8 xyz(3,Natom),xyznew(3,Natom)
	real*8 sum,nhat(3,2*Nres+Nxi),deldeg
        real*8 V1(3),V2(3),V3(Natom,3),rnk(Natom,3)

	do i=1,Natom
	 xyznew(1,i) = xyz(1,i)
	 xyznew(2,i) = xyz(2,i)
	 xyznew(3,i) = xyz(3,i)
	enddo

	j3 = IndexXi(Ixi,3)
	Ientry = 2*Nres + Ixi

	do i = j3+1,EndXi(Ixi)
	 rnk(i,1) = xyznew(1,i) - xyznew(1,j3)
	 rnk(i,2) = xyznew(2,i) - xyznew(2,j3)
	 rnk(i,3) = xyznew(3,i) - xyznew(3,j3)
	 sum = nhat(1,Ientry)*rnk(i,1) + 
     1	       nhat(2,Ientry)*rnk(i,2) + nhat(3,Ientry)*rnk(i,3)
	 V3(i,1) = nhat(1,Ientry) * sum
	 V3(i,2) = nhat(2,Ientry) * sum
	 V3(i,3) = nhat(3,Ientry) * sum
        enddo

	do i = j3+1,EndXi(Ixi)
	 V1(1) = nhat(2,Ientry)*rnk(i,3) - nhat(3,Ientry)*rnk(i,2)
	 V1(2) = nhat(3,Ientry)*rnk(i,1) - nhat(1,Ientry)*rnk(i,3)
	 V1(3) = nhat(1,Ientry)*rnk(i,2) - nhat(2,Ientry)*rnk(i,1)
	 V2(1) = -nhat(2,Ientry)*V1(3) + nhat(3,Ientry)*V1(2)
	 V2(2) = -nhat(3,Ientry)*V1(1) + nhat(1,Ientry)*V1(3)
	 V2(3) = -nhat(1,Ientry)*V1(2) + nhat(2,Ientry)*V1(1)
	 xyznew(1,i) = V1(1)*dsin(deldeg) + V2(1)*dcos(deldeg)
     1    	     + V3(i,1) + xyznew(1,j3)
	 xyznew(2,i) = V1(2)*dsin(deldeg) + V2(2)*dcos(deldeg)
     1     	     + V3(i,2) + xyznew(2,j3)
	 xyznew(3,i) = V1(3)*dsin(deldeg) + V2(3)*dcos(deldeg)
     1     	     + V3(i,3) + xyznew(3,j3)
	enddo     


	return
	end
