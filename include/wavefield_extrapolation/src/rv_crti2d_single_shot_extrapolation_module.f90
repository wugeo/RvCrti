



!***********************************************************************
	  module	module_rv_born_extrapolation_2d

	  	implicit none

	  	!=========================================================
    	!wavefield 
	  	type wavefield
	  		real,allocatable	::	u(:, :)
	  	end	type

		!* us0, ur0
		type(wavefield), target	::	us01, us02, us03
		type(wavefield), target	::	ur01, ur02, ur03

	  	type(wavefield), target	::	ux0_a1,  ux0_a2,  ux0_a3
	  	type(wavefield), target	::	ux0_b1,  ux0_b2
	  	type(wavefield), target	::	ux0_bp1, ux0_bp2, ux0_bp3
	  	type(wavefield), target	::	ux0_c1,  ux0_c2,  ux0_c3

	  	type(wavefield), target	::	uz0_a1,  uz0_a2,  uz0_a3
	  	type(wavefield), target	::	uz0_b1,  uz0_b2
	  	type(wavefield), target	::	uz0_bp1, uz0_bp2, uz0_bp3
	  	type(wavefield), target	::	uz0_c1,  uz0_c2,  uz0_c3


	  	type(wavefield), pointer	::	pus01, pus02, pus03
	  	type(wavefield), pointer	::	pur01, pur02, pur03

	  	type(wavefield), pointer	::	pux0_a1,  pux0_a2,  pux0_a3
	  	type(wavefield), pointer	::	pux0_b1,  pux0_b2
	  	type(wavefield), pointer	::	pux0_bp1, pux0_bp2, pux0_bp3
	  	type(wavefield), pointer	::	pux0_c1,  pux0_c2,  pux0_c3

	  	type(wavefield), pointer	::	puz0_a1,  puz0_a2,  puz0_a3
	  	type(wavefield), pointer	::	puz0_b1,  puz0_b2
	  	type(wavefield), pointer	::	puz0_bp1, puz0_bp2, puz0_bp3
	  	type(wavefield), pointer	::	puz0_c1,  puz0_c2,  puz0_c3


		!* usr, urr
		type(wavefield), target	::	usr1, usr2, usr3
		type(wavefield), target	::	urr1, urr2, urr3

	  	type(wavefield), target	::	uxr_a1,  uxr_a2,  uxr_a3
	  	type(wavefield), target	::	uxr_b1,  uxr_b2
	  	type(wavefield), target	::	uxr_bp1, uxr_bp2, uxr_bp3
	  	type(wavefield), target	::	uxr_c1,  uxr_c2,  uxr_c3

	  	type(wavefield), target	::	uzr_a1,  uzr_a2,  uzr_a3
	  	type(wavefield), target	::	uzr_b1,  uzr_b2
	  	type(wavefield), target	::	uzr_bp1, uzr_bp2, uzr_bp3
	  	type(wavefield), target	::	uzr_c1,  uzr_c2,  uzr_c3


	  	type(wavefield), pointer	::	pusr1, pusr2, pusr3
	  	type(wavefield), pointer	::	purr1, purr2, purr3

	  	type(wavefield), pointer	::	puxr_a1,  puxr_a2,  puxr_a3
	  	type(wavefield), pointer	::	puxr_b1,  puxr_b2
	  	type(wavefield), pointer	::	puxr_bp1, puxr_bp2, puxr_bp3
	  	type(wavefield), pointer	::	puxr_c1,  puxr_c2,  puxr_c3

	  	type(wavefield), pointer	::	puzr_a1,  puzr_a2,  puzr_a3
	  	type(wavefield), pointer	::	puzr_b1,  puzr_b2
	  	type(wavefield), pointer	::	puzr_bp1, puzr_bp2, puzr_bp3
	  	type(wavefield), pointer	::	puzr_c1,  puzr_c2,  puzr_c3


		!************************************************************
		type(wavefield), target	::	uv_r01, uv_r02, uv_r03

	  	type(wavefield), target	::	uv_x0_a1,  uv_x0_a2,  uv_x0_a3
	  	type(wavefield), target	::	uv_x0_b1,  uv_x0_b2
	  	type(wavefield), target	::	uv_x0_bp1, uv_x0_bp2, uv_x0_bp3
	  	type(wavefield), target	::	uv_x0_c1,  uv_x0_c2,  uv_x0_c3

	  	type(wavefield), target	::	uv_z0_a1,  uv_z0_a2,  uv_z0_a3
	  	type(wavefield), target	::	uv_z0_b1,  uv_z0_b2
	  	type(wavefield), target	::	uv_z0_bp1, uv_z0_bp2, uv_z0_bp3
	  	type(wavefield), target	::	uv_z0_c1,  uv_z0_c2,  uv_z0_c3

	  	type(wavefield), pointer	::	puv_r01, puv_r02, puv_r03

	  	type(wavefield), pointer	::	puv_x0_a1,  puv_x0_a2,  puv_x0_a3
	  	type(wavefield), pointer	::	puv_x0_b1,  puv_x0_b2
	  	type(wavefield), pointer	::	puv_x0_bp1, puv_x0_bp2, puv_x0_bp3
	  	type(wavefield), pointer	::	puv_x0_c1,  puv_x0_c2,  puv_x0_c3

	  	type(wavefield), pointer	::	puv_z0_a1,  puv_z0_a2,  puv_z0_a3
	  	type(wavefield), pointer	::	puv_z0_b1,  puv_z0_b2
	  	type(wavefield), pointer	::	puv_z0_bp1, puv_z0_bp2, puv_z0_bp3
	  	type(wavefield), pointer	::	puv_z0_c1,  puv_z0_c2,  puv_z0_c3

		!************************************************************
		type(wavefield), target	::	uv_sr1, uv_sr2, uv_sr3
		type(wavefield), target	::	uv_rr1, uv_rr2, uv_rr3

	  	type(wavefield), target	::	uv_xr_a1,  uv_xr_a2,  uv_xr_a3
	  	type(wavefield), target	::	uv_xr_b1,  uv_xr_b2
	  	type(wavefield), target	::	uv_xr_bp1, uv_xr_bp2, uv_xr_bp3
	  	type(wavefield), target	::	uv_xr_c1,  uv_xr_c2,  uv_xr_c3

	  	type(wavefield), target	::	uv_zr_a1,  uv_zr_a2,  uv_zr_a3
	  	type(wavefield), target	::	uv_zr_b1,  uv_zr_b2
	  	type(wavefield), target	::	uv_zr_bp1, uv_zr_bp2, uv_zr_bp3
	  	type(wavefield), target	::	uv_zr_c1,  uv_zr_c2,  uv_zr_c3

	  	type(wavefield), pointer	::	puv_sr1, puv_sr2, puv_sr3
	  	type(wavefield), pointer	::	puv_rr1, puv_rr2, puv_rr3

	  	type(wavefield), pointer	::	puv_xr_a1,  puv_xr_a2,  puv_xr_a3
	  	type(wavefield), pointer	::	puv_xr_b1,  puv_xr_b2
	  	type(wavefield), pointer	::	puv_xr_bp1, puv_xr_bp2, puv_xr_bp3
	  	type(wavefield), pointer	::	puv_xr_c1,  puv_xr_c2,  puv_xr_c3

	  	type(wavefield), pointer	::	puv_zr_a1,  puv_zr_a2,  puv_zr_a3
	  	type(wavefield), pointer	::	puv_zr_b1,  puv_zr_b2
	  	type(wavefield), pointer	::	puv_zr_bp1, puv_zr_bp2, puv_zr_bp3
	  	type(wavefield), pointer	::	puv_zr_c1,  puv_zr_c2,  puv_zr_c3

		!* pt
	  	type(wavefield), pointer	::	pt

	  	!=========================================================
	  	!wavefield3
	  	type wavefield3
	  		real,allocatable	::	u(:,:,:)
	 	end type

	  	type(wavefield3), target	::	top0, bot0, lef0, rig0
	  	type(wavefield3), target	::	ust0
	  	type(wavefield3), target	::	topr, botr, lefr, rigr
	  	type(wavefield3), target	::	ustr

	  	type(wavefield3), target	::	topv_r, botv_r, lefv_r, rigv_r
	  	type(wavefield3), target	::	ustv_r


	  	real,allocatable	::	image_tmp(:,:)
		real,allocatable	::	image_ms(:,:), image_dm(:,:)
		real,allocatable	::	image_rv_ms(:,:), image_rv_dm(:,:)
	  	
	  	
		real,allocatable	::	coe_1st(:), coe_2nd(:)
		real,allocatable	::	coe_1st_dx(:), coe_1st_dz(:)
		real,allocatable	::	coe_2nd_dx2(:), coe_2nd_dz2(:)
		real,allocatable	::	funa(:,:), dfuna(:,:)
		integer,allocatable	::	line_xl(:,:), line_xr(:,:)
		integer,allocatable	::	line_zu(:,:), line_zd(:,:)

	  	

	  end module 	module_rv_born_extrapolation_2d
!***********************************************************************
