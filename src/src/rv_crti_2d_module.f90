

!***********************************************************************
	module	RV_CRTI_2d_shot_parameter
		use	global	
		implicit none

	  	integer,parameter	::	npad	=	5
	  	integer,parameter	::	order	=	11
	  	integer,parameter   ::  iflag_src = 1
	  	integer,parameter   ::  it_snap = 500
	  	real,parameter		::	r = 1.0e-4

		character(len=para_char_flen)	::	fn_vel
		character(len=para_char_flen)	::	fn_image
		character(len=para_char_flen)	::	fn_cs

		character(len=para_char_flen)	::	currt
		character(len=para_char_flen)	::	currtfile


		real	::	fmain
	  	integer	::	lt
	  	integer	::	ns
	  	integer	::	ntr
	  	integer	::	nwt
	  	integer	::	nvx, nvz
	  	integer	::	nx, nz
	  	integer	::	nnx, nnz
	  	integer	::	nx_with_apert, nz_with_apert
		integer ::	nx_apert_l, nx_apert_r
		integer ::	nz_apert_u, nz_apert_d
		integer ::	nx_bound_l, nx_bound_r
		integer ::	nz_bound_u, nz_bound_d
		integer ::	nvxx_shift, nvzz_shift
		integer	::	nx_shift, nz_shift
		integer	::	ns_x, ns_z
		integer	::	nr_x, nr_z

	  	real	::	dt
	  	real	::	dvx, dvz
		real	::	dx, dz
		real	::	cvx_initial, cvz_initial
		real	::	cx_min, cz_min
		real	::	sx, sz

		integer	::	iflag_order
	  	real	::	coe(18)
	  	integer ::  it_step

		integer	::	flag_sr
	  	integer	::	wflag
		integer	::	ierr
		integer	::	isok

	end module	RV_CRTI_2d_shot_parameter











