!***********************************************************************

	program	RV_CRTI_2D_ConstantDensityAcousticEquation

		use global
		use RV_CRTI_2d_shot_parameter 
		implicit none

		character(len=para_char_flen)	::	fn_par
		real,allocatable	::	wavelet(:)
		real,allocatable	::	gx(:), gz(:)
		integer,allocatable	::	index_gx(:), index_gz(:)
		real,allocatable	::	vel(:,:), vv(:,:)
		real,allocatable	::	image_r(:, :), image_i(:,:)
		real,allocatable	::	image_r_sig(:, :), image_i_sig(:,:)
		real,allocatable	::	image_r_loc(:, :), image_i_loc(:,:)
		real,allocatable	::	shotobs(:,:), shotcal(:,:), shotres(:,:), shotres_rv(:,:)
		real,allocatable	::	grad_to_sig(:, :), grad_rv_sig(:,:)
		integer		::	itr


		CALl Getarg(1, fn_par)
		CALL RV_CRTI2D_ReadPar(fn_par)

		nx = 1001
		nz = 401
		nx_with_apert = nx_apert_l + nx + nx_apert_r
		nz_with_apert = nz_apert_u + nz + nz_apert_d
		nnx = nx_bound_l +  nx_with_apert + nx_bound_r
		nnz = nz_bound_u +  nz_with_apert + nz_bound_d
		nvxx_shift = (cx_min-cvx_initial)/dvx - nx_apert_l - 1
		nvzz_shift = (cz_min-cvz_initial)/dvz - nz_apert_u - 1
		nx_shift = nx_apert_l + nx_bound_l
		nz_shift = nz_apert_u + nz_bound_u

		it_step = 1
		ns = lt
		dt = dt*0.001


		allocate(gx(ntr), gz(ntr), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"(gx, gz), Can not allocate working memory, stop!!!"
			stop
		endif
		gx = 0.0
		gz = 0.0

		allocate(index_gx(ntr), index_gz(ntr), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"(index_gx, index_gz), Can not allocate working memory, stop!!!"
			stop
		endif
		index_gx = 0
		index_gz = 0

		allocate(wavelet(lt), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"(wavelet), Can not allocate working memory, stop!!!"
			stop
		endif
		wavelet = 0.0

		allocate(image_r(nvz, nvx), image_i(nvz, nvx), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"(image_r, image_i), Can not allocate working memory, stop!!!"
			stop
		endif
		image_r = 0.0
		image_i = 0.0

		allocate(image_r_sig(nz_with_apert, nx_with_apert), image_i_sig(nz_with_apert, nx_with_apert), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"(image_r_sig, image_i_sig), Can not allocate working memory, stop!!!"
			stop
		endif
		image_r_sig = 0.0
		image_i_sig = 0.0

		allocate(image_r_loc(nnz, nnx), image_i_loc(nnz, nnx), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"(image_r_loc, image_i_loc), Can not allocate working memory, stop!!!"
			stop
		endif
		image_r_loc = 0.0
		image_i_loc = 0.0

		allocate(shotobs(lt, ntr), shotcal(lt, ntr), shotres(lt, ntr), shotres_rv(lt, ntr), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"(shotobs, shotcal, shotres), Can not allocate working memory, stop!!!"
			stop
		endif
		shotobs = 0.0
		shotcal = 0.0
		shotres = 0.0
		shotres_rv = 0.0

		allocate(vel(nvz, nvx), vv(nnz, nnx), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"(vv), Can not allocate working memory, stop!!!"
			stop
		endif
		vv = 0.0

		allocate(grad_to_sig(nz_with_apert, nx_with_apert), grad_rv_sig(nz_with_apert, nx_with_apert), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"(grad_to_sig, grad_rv_sig), Can not allocate working memory, stop!!!"
			stop
		endif
		grad_to_sig = 0.0
		grad_rv_sig = 0.0

		
		CALl read_float2d(fn_vel, vel, nvz, nvx, isok)
		CALl read_float2d(fn_image, image_r, nvz, nvx, isok)
		CALl read_float2d(fn_cs, shotobs, lt, ntr, isok)
		CALL Wavelet_Forming(wavelet, dt, lt, fmain, nwt)


		do itr = 1, ntr
			gx(itr) = (sx-ntr/2*dx) + (itr-1)*dx
			gz(itr) = 0.0
		end do

    	CALL Get_Current_Shot_Velocity_2d(vel, vv, &
					nvz, nvx, nnz, nnx, &
					nz_bound_u, nx_bound_l, &
					nvzz_shift, nvxx_shift, &
					isok)


    	CALL Get_Current_Shot_Velocity_2d(image_r, image_r_loc, &
					nvz, nvx, nnz, nnx, &
					nz_bound_u, nx_bound_l, &
					nvzz_shift, nvxx_shift, &
					isok)


		CALL Source_and_Receiver_Position_2d(sx, sz, gx, gz, &
				ns_x, ns_z, index_gx, index_gz, &
				ntr, cx_min, cz_min, &
				dx, dz, &
				nx_shift, nz_shift, &
				isok)
					
		CALL RV_CRTI_2D_image_single_shot(wavelet, vv, &
					shotobs, shotcal, shotres, &
					image_r_loc, image_r_sig, &
					gx, gz, index_gx, index_gz)


		CALL local_to_all_image_result_2d(image_i, image_r_sig, image_r, &
				nz_with_apert, nx_with_apert, nvz, nvx, &
				nx, nz, nx_apert_l, nz_apert_u, &
				nvxx_shift, nvzz_shift, &
				isok)


    	CALL Get_Current_Shot_Velocity_2d(image_i, image_i_loc, &
					nvz, nvx, nnz, nnx, &
					nz_bound_u, nx_bound_l, &
					nvzz_shift, nvxx_shift, &
					isok)
			

		CALL RV_CRTI_2D_gradient_single_shot(wavelet, vv, &
					shotobs, shotcal, shotres, shotres_rv, &
					image_r_loc, image_i_loc, grad_to_sig, grad_rv_sig, &
					gx, gz, index_gx, index_gz)


		deallocate(gx, gz)
		deallocate(index_gx, index_gz)
		deallocate(wavelet)
		deallocate(image_r, image_i)
		deallocate(image_r_sig, image_i_sig)
		deallocate(image_r_loc, image_i_loc)
		deallocate(vv, vel)
		deallocate(shotobs, shotcal, shotres)
		deallocate(shotres_rv)
		deallocate(grad_to_sig, grad_rv_sig)

	end 

!***********************************************************************












