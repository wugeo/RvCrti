!*************************************************************************
	subroutine	rv_crti_2d_shotgathers_image(shotobs, shotcal, shotres, &
					vv, image_r, misfit_loc, &
					sx, sz, gx, gz, index_recgx, index_recgz, ns_x, ns_z, &
                    nnz, nnx, lt, ntr, nwt, dt, dx, dz)

		use global
		implicit none
		integer	::	nnz, nnx
		integer	::	lt, ntr, nwt
		real	::	dt, dx, dz

		real	::	shotobs(lt, ntr), shotcal(lt, ntr)
		real	::	shotres(lt, ntr)
		real	::	vv(nnz, nnx), image_r(nnz, nnx)
		real	::	misfit_loc

		integer	::	ns_x, ns_z
		real	::	sx, sz
		real	::	gx(ntr), gz(ntr)
		integer	::	index_recgx(ntr), index_recgz(ntr)

		integer	::	flag_method
		character(len=para_char_flen)  :: 	fn_name
		character(len=para_char_flen)  :: 	currt1, currt2, currtsnap
		real	::	misfit_sig_sum
		integer	::	flag_itr
		integer	::	ix, inx, inz
		integer	::	it, itr, iit
		integer	::	ierr
		integer	::	crnum_cur
		integer	::	isok

		real,allocatable	::	shot_dt(:)
		real	::	amp_obs, amp_cal
		integer	::	it0_obs, it0_cal
		real	::	misfit

		allocate(shot_dt(ntr))
		shot_dt = 0.0

		do itr = 1, ntr
			call fmaximum_1d(shotobs(:,itr), lt, it0_obs, amp_obs)
			call fmaximum_1d(shotcal(:,itr), lt, it0_cal, amp_cal)
			shot_dt(itr) = 	(it0_obs - it0_cal)*dt
		end do

		shotres = 0.0
		CALL cal_AdjointSource_crosscorrelation(shotcal, shotobs, shotres, misfit, &
					shot_dt, lt, ntr, dt, nwt)

		deallocate(shot_dt)

		return
    end	subroutine


!*************************************************************************
	subroutine	rv_crti_2d_shotgathers_gradient(shotobs, shotcal, shotres, shotres_rv, &
					vv, image_r, misfit_loc, &
					sx, sz, gx, gz, index_recgx, index_recgz, ns_x, ns_z, &
                    nnz, nnx, lt, ntr, nwt, dt, dx, dz)

		use global
		implicit none
		integer	::	nnz, nnx
		integer	::	lt, ntr, nwt
		real	::	dt, dx, dz

		real	::	shotobs(lt, ntr), shotcal(lt, ntr)
		real	::	shotres(lt, ntr), shotres_rv(lt, ntr)
		real	::	vv(nnz, nnx), image_r(nnz, nnx)
		real	::	misfit_loc

		integer	::	ns_x, ns_z
		real	::	sx, sz
		real	::	gx(ntr), gz(ntr)
		integer	::	index_recgx(ntr), index_recgz(ntr)

		character(len=para_char_flen)	::	fn_debug
		integer	::	iflag_debug, debug_ishot
		integer	::	num_grad, ishot
		integer	::	isok

		integer	::	flag_method
		character(len=para_char_flen)  :: 	fn_name
		character(len=para_char_flen)  :: 	currt1, currt2, currtsnap
		real	::	misfit_sig_sum1, misfit_sig_sum2
		integer	::	flag_itr
		integer	::	ix, inx, inz
		integer	::	it, itr, iit
		integer	::	ierr
		integer	::	crnum_cur


		real,allocatable	::	shot_dt(:)
		real	::	amp_obs, amp_cal
		integer	::	it0_obs, it0_cal
		real	::	offset_limit 
		real	::	gx_target
		real	::	misfit


		allocate(shot_dt(ntr))
		shot_dt = 0.0

		gx_target = 5000.0
		do itr = 1, ntr
			if ( gx(itr) .eq. gx_target ) then
				call fmaximum_1d(shotobs(:,itr), lt, it0_obs, amp_obs)
				call fmaximum_1d(shotcal(:,itr), lt, it0_cal, amp_cal)
				shot_dt(itr) = 	(it0_obs - it0_cal)*dt
			end if
		end do

		shotres = 0.0
		CALL cal_AdjointSource_crosscorrelation(shotcal, shotobs, shotres, misfit, &
					shot_dt, lt, ntr, dt, nwt)

		shotres_rv = 0.0
		offset_limit = 20.0
		do itr = 1, ntr
			if ( abs(sx - gx(itr)) .le. offset_limit ) then
                shotres_rv(:, itr) = -1.0*shotcal(:, itr)
			end if
		end do

		deallocate(shot_dt)

		return
    end	subroutine

!***********************************************************************
