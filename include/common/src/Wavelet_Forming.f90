
!***********************************************************************
	subroutine Wavelet_Forming(wavelet, dt, lt, fmain, nwt)

		implicit none
		!Dummy Variables
		integer ::	lt
		integer ::	nwt
		real	::	dt
		real	::	fmain
		real	::	wavelet(lt)

		!Local Variables
		real,parameter	::	pi=3.14159265359
		integer it
		real	tmain
		real	tp1, tp2
		real	time

		tmain = 1.0/fmain
		nwt = tmain/dt + 0.5
			
		do it=1, lt
			time = (it-1)*dt-tmain

			tp1 = pi*fmain*time
			tp2 = tp1*tp1
			wavelet(it) = (1.0-2.0*tp2)*exp(-tp2)
		enddo

		return
    end subroutine


!***********************************************************************
