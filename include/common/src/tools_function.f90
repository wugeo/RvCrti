!***********************************************************************
    subroutine normalization_f1d_array(v, num)

		implicit none
		integer	::	num
		real	::	v(num)

		real,parameter	::	eps=1.0e-8
		integer	::	it, itr
		real	::	tmp

		tmp = abs(maxval(v))
		if(tmp .lt. eps)	tmp=1.0
		v = v/tmp;

		return
    end	subroutine

!***********************************************************************
    subroutine normalization_f2d_array(v, n1, n2)

		implicit none
		integer	::	n1, n2
		real	::	v(n1, n2)

		real,parameter	::	eps=1.0e-8
		integer	::	it, itr
		real	::	tmp

		tmp = abs(maxval(v))
		if(tmp .lt. eps)	tmp=1.0
		v = v/tmp;

		return
    end	subroutine

!***********************************************************************

	subroutine fmaximum_1d(data, len, ii, amp)
		implicit none
		integer, intent(in)  :: len
		real   , intent(in)  :: data(len)
		integer, intent(out) :: ii
		real   , intent(out) :: amp

		integer :: i
		real    :: maxval

		maxval = data(1)
		ii = 1
		do i = 2, len
			if (data(i) > maxval) then
				maxval = data(i)
				ii = i
			end if
		end do

		amp = maxval
	end subroutine fmaximum_1d



	subroutine	cal_AdjointSource_crosscorrelation(shotcal, shotobs, shotres, misfit, &
					shotdt, lt, ntr, dt, nwt)

		use global
		implicit none
		real	::	misfit
		integer	::	lt, ntr, nwt
		real	::	dt

		real	::	shotcal(lt, ntr)
		real	::	shotobs(lt, ntr)
		real	::	shotres(lt, ntr)
		real	::	shotdt(ntr)

		integer	::	itr, it
		integer	::	flag_mod
		real	::	tmp1, tmp2, tmp3, gama
		real	::	sigma=0.00001
		
		real	::	coe1, coe2
		real	::	taper, taper1
		real	::	recgx, recgz
		real	::	offset, offset_min, offset_max
		real	::	tfb, t1, t2
		real	::	win
		integer	::	it1, it2, nwin
		real	::	dis
		real	::	coeff
		real,allocatable	::	trace_cal(:)
		real,allocatable	::	trace_obs(:)
		real,allocatable	::	trace(:)
		real,allocatable	::	trace_tmp(:)
		integer,allocatable	::	itp(:)
		real	::	tmin, tmax
		real	::	deltat
		real	::	tmp, tmp_max
		real	::	amp_max
		integer	::	iit, ip, ii
		integer	::	flag_itp

		character(len=para_char_flen)  ::  currt1, currt2, currt3, currt4
		character(len=para_char_flen)  ::  currtfile


		!===================================================================
		allocate(trace_cal(lt))
		allocate(trace_obs(lt))
		allocate(trace(lt))
		allocate(trace_tmp(lt))

		misfit=0.0
		do itr=1, ntr
			trace_cal=0.0
			trace_obs=0.0
			trace=0.0
			trace_tmp=0.0

			do it=1, lt
				trace_cal(it)=shotcal(it, itr)
				trace_obs(it)=shotobs(it, itr)
			enddo
			deltat=shotdt(itr)

			call gather_derivative_1d(trace_cal, trace, lt, dt)
			do it=1, lt
				shotres(it, itr)=trace(it)*deltat
			enddo

			call gather_derivative_2d(trace_cal, trace, lt, dt)
			tmp=0.0
			do it=1, lt
				tmp=tmp+ trace(it)*trace_cal(it)
			enddo
			do it=1, lt
				shotres(it, itr)=-1.0*shotres(it, itr)/(tmp+sigma)
			enddo

			misfit=misfit + deltat*deltat
		enddo

		deallocate(trace_cal)
		deallocate(trace_obs)
		deallocate(trace)
		deallocate(trace_tmp)

	end subroutine


!***********************************************************************
	subroutine gather_derivative_2d(trace_cal, trace, lt, dt)

		use global
		implicit none
		real	::	trace_cal(lt)
		real	::	trace(lt)
		real	::	dt
		integer	::	lt

		real	::	b(5)
		real	::	tmp
		integer	::	it, it1, it2
		integer	::	i

		b(1)=1.666667
		b(2)=-0.238095
		b(3)=0.03968254
		b(4)=-0.00496
		b(5)=0.00031746

		trace=0.0
		do it=6, lt-5
			tmp=0.0
			do i=1,5
				tmp=tmp+b(i)*(trace_cal(it-i)+trace_cal(it+i)-2.0*trace_cal(it))/(dt*dt)
			enddo
			trace(it)=tmp
		enddo
		

	end subroutine

!***********************************************************************
	subroutine gather_derivative_1d(trace, trace_cal, lt, dt)

		use global
		implicit none
		real	::	trace_cal(lt)
		real	::	trace(lt)
		real	::	dt
		integer	::	lt

		real	::	tmp1, tmp2
		integer	::	it, it1, it2

		real	::	poy_coe(5)

		poy_coe(1)=0.8333333
		poy_coe(2)=-0.2380953
		poy_coe(3)=5.9523813E-02
		poy_coe(4)=-9.9206381E-03
		poy_coe(5)=7.9365104E-04


		trace_cal=0.0
		do it=6, lt-5
			trace_cal(it)= &
				(	(trace(it+1)-trace(it-1))*poy_coe(1)+ &
					(trace(it+2)-trace(it-2))*poy_coe(2)+ &
					(trace(it+3)-trace(it-3))*poy_coe(3)+ &
					(trace(it+4)-trace(it-4))*poy_coe(4)+ &
					(trace(it+5)-trace(it-5))*poy_coe(5)  &
				)/dt
		enddo


	end subroutine











