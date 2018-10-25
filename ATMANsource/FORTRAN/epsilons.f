	subroutine epsilons(type,eps)

	implicit none

	character*1 type
	real*8 eps

c values taken as closely as possible from L79 (precursor of ENCAD)


	if(TYPE.eq.' ')then
		write(*,*)'Error in subroutine epsilons'
		stop
	endif
	if(TYPE.eq.'H')eps=0.03800d0
	if(TYPE.eq.'O')eps=0.18479d0
	if(TYPE.eq.'Q')eps=0.18479d0
	if(TYPE.eq.'V')eps=0.18479d0
	if(TYPE.eq.'N')eps=0.41315d0
	if(TYPE.eq.'M')eps=0.41315d0
	if(TYPE.eq.'C')eps=0.07382d0
	if(TYPE.eq.'S')eps=0.07382d0
	if(TYPE.eq.'A')eps=0.03763d0
	if(TYPE.eq.'B')eps=0.03763d0
	if(TYPE.eq.'G')eps=0.03763d0

	return
	end


