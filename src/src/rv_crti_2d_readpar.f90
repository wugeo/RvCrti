

!***********************************************************************
    subroutine RV_CRTI2D_ReadPar(fn_par)

		use global
		use RV_CRTI_2d_shot_parameter 
		implicit none
		character(len=para_char_flen)	::	fn_par
		character(len=para_char_flen)	::	par_name

		!********************************************************************
		open(33, file=fn_par, status='old', iostat=ierr)
		if(ierr /= 0)then
			write(*,*)	'Parameter file cannot open right, please check it !!!'
			write(*,*)	'Shutting down the program'
			stop
		endif

		read(33, '(a)')	par_name
		read(33, '(a)') fn_vel

		read(33, '(a)')	par_name
		read(33, '(a)') fn_image

		read(33, '(a)')	par_name
		read(33, '(a)') fn_cs

		read(33, '(a)')	par_name
		read(33, *)	lt, dt

		read(33, '(a)')	par_name
		read(33, *)	nvx, nvz

		read(33, '(a)')	par_name
		read(33, *)	dvx, dvz

		read(33, '(a)')	par_name
		read(33, *)	ntr

		read(33, '(a)')	par_name
		read(33, *)	sx, sz

		read(33, '(a)')	par_name
		read(33, *)	dx, dz

		read(33, '(a)')	par_name
		read(33, *)	nx_apert_l, nx_apert_r, nz_apert_u, nz_apert_d

		read(33, '(a)')	par_name
		read(33, *)	nx_bound_l, nx_bound_r, nz_bound_u, nz_bound_d

		read(33, '(a)')	par_name
		read(33, *)	cx_min, cz_min

		read(33, '(a)')	par_name
		read(33, *)	cvx_initial, cvz_initial 

		read(33, '(a)')	par_name
		read(33, *)	fmain

		close(33)

		return
    end	subroutine

!***********************************************************************
