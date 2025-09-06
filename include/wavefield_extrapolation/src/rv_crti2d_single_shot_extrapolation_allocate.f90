
!***********************************************************************
	subroutine rv_born_2d_single_shot_extrapolation_allocate(nnx, nnz, &
					nx_with_apert, nz_with_apert, &
					nx_bound_l, nx_bound_r, &
					nz_bound_u, nz_bound_d, &
					npad, order, ns, iflag_src)

		use module_rv_born_extrapolation_2d
		implicit none
		integer		::	nnx, nnz
		integer		::	nx_with_apert, nz_with_apert
		integer		::	nx_bound_l, nx_bound_r
		integer		::	nz_bound_u, nz_bound_d
		integer		::	npad
		integer		::	order, ns
		integer		::	iflag_src
		integer		::	ierr
	

		!******************************************************************
	!allocate wavefiled  us0, ur0
		allocate(us01%u(-npad+1:nnz+npad, -npad+1:nnx+npad), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about us01, stop!!!"
			stop
		endif

		allocate(us02%u(-npad+1:nnz+npad, -npad+1:nnx+npad), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about us02, stop!!!"
			stop
		endif

		allocate(us03%u(-npad+1:nnz+npad, -npad+1:nnx+npad), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about us03, stop!!!"
			stop
		endif

		allocate(ur01%u(-npad+1:nnz+npad, -npad+1:nnx+npad), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about ur01, stop!!!"
			stop
		endif

		allocate(ur02%u(-npad+1:nnz+npad, -npad+1:nnx+npad), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about ur02, stop!!!"
			stop
		endif

		allocate(ur03%u(-npad+1:nnz+npad, -npad+1:nnx+npad), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about ur03, stop!!!"
			stop
		endif


		!******************************************************************
		!allocate pml
		!x direction
		allocate(ux0_a1%u(nnz,nx_bound_l+nx_bound_r), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about ux0_a1, stop!!!"
			stop
		endif

		allocate(ux0_a2%u(nnz,nx_bound_l+nx_bound_r), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about ux0_a2, stop!!!"
			stop
		endif

		allocate(ux0_a3%u(nnz,nx_bound_l+nx_bound_r), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about ux0_a3, stop!!!"
			stop
		endif

		allocate(ux0_b1%u(nnz,nx_bound_l+nx_bound_r), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about ux0_b1, stop!!!"
			stop
		endif

		allocate(ux0_b2%u(nnz,nx_bound_l+nx_bound_r), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about ux0_b2, stop!!!"
			stop
		endif

		allocate(ux0_bp1%u(nnz,nx_bound_l+nx_bound_r), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about ux0_bp1, stop!!!"
			stop
		endif

		allocate(ux0_bp2%u(nnz,nx_bound_l+nx_bound_r), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about ux0_bp2, stop!!!"
			stop
		endif

		allocate(ux0_bp3%u(nnz,nx_bound_l+nx_bound_r), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about ux0_bp3, stop!!!"
			stop
		endif

		allocate(ux0_c1%u(nnz,nx_bound_l+nx_bound_r), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about ux0_c1, stop!!!"
			stop
		endif

		allocate(ux0_c2%u(nnz,nx_bound_l+nx_bound_r), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about ux0_c2, stop!!!"
			stop
		endif

		allocate(ux0_c3%u(nnz,nx_bound_l+nx_bound_r), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about ux0_c3, stop!!!"
			stop
		endif


		!z direction
		allocate(uz0_a1%u(nz_bound_u+nz_bound_d,nnx), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about uz0_a1, stop!!!"
			stop
		endif

		allocate(uz0_a2%u(nz_bound_u+nz_bound_d,nnx), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about uz0_a2, stop!!!"
			stop
		endif

		allocate(uz0_a3%u(nz_bound_u+nz_bound_d,nnx), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about uz0_a3, stop!!!"
			stop
		endif

		allocate(uz0_b1%u(nz_bound_u+nz_bound_d,nnx), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about uz0_b1, stop!!!"
			stop
		endif

		allocate(uz0_b2%u(nz_bound_u+nz_bound_d,nnx), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about uz0_b2, stop!!!"
			stop
		endif

		allocate(uz0_bp1%u(nz_bound_u+nz_bound_d,nnx), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about uz0_bp1, stop!!!"
			stop
		endif

		allocate(uz0_bp2%u(nz_bound_u+nz_bound_d,nnx), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about uz0_bp2, stop!!!"
			stop
		endif

		allocate(uz0_bp3%u(nz_bound_u+nz_bound_d,nnx), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about uz0_bp3, stop!!!"
			stop
		endif

		allocate(uz0_c1%u(nz_bound_u+nz_bound_d,nnx), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about uz0_c1, stop!!!"
			stop
		endif

		allocate(uz0_c2%u(nz_bound_u+nz_bound_d,nnx), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about uz0_c2, stop!!!"
			stop
		endif

		allocate(uz0_c3%u(nz_bound_u+nz_bound_d,nnx), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about uz0_c3, stop!!!"
			stop
		endif

		!******************************************************************
		!wavefields
		us01%u=0.0
		us02%u=0.0
		us03%u=0.0
		ur01%u=0.0
		ur02%u=0.0
		ur03%u=0.0

		!x direction
		ux0_a1%u=0.0
		ux0_a2%u=0.0
		ux0_a3%u=0.0
		ux0_b1%u=0.0
		ux0_b2%u=0.0
		ux0_bp1%u=0.0
		ux0_bp2%u=0.0
		ux0_bp3%u=0.0
		ux0_c1%u=0.0
		ux0_c2%u=0.0
		ux0_c3%u=0.0

		!z direction
		uz0_a1%u=0.0
		uz0_a2%u=0.0
		uz0_a3%u=0.0
		uz0_b1%u=0.0
		uz0_b2%u=0.0
		uz0_bp1%u=0.0
		uz0_bp2%u=0.0
		uz0_bp3%u=0.0
		uz0_c1%u=0.0
		uz0_c2%u=0.0
		uz0_c3%u=0.0

		!******************************************************************
		!wavefields
		pus01 => us01
		pus02 => us02
		pus03 => us03
		pur01 => ur01
		pur02 => ur02
		pur03 => ur03

		!x direction
		pux0_a1 => ux0_a1
		pux0_a2 => ux0_a2
		pux0_a3 => ux0_a3
		pux0_b1 => ux0_b1
		pux0_b2 => ux0_b2
		pux0_bp1 => ux0_bp1
		pux0_bp2 => ux0_bp2
		pux0_bp3 => ux0_bp3
		pux0_c1 => ux0_c1
		pux0_c2 => ux0_c2
		pux0_c3 => ux0_c3

		!z direction
		puz0_a1 => uz0_a1
		puz0_a2 => uz0_a2
		puz0_a3 => uz0_a3
		puz0_b1 => uz0_b1
		puz0_b2 => uz0_b2
		puz0_bp1 => uz0_bp1
		puz0_bp2 => uz0_bp2
		puz0_bp3 => uz0_bp3
		puz0_c1 => uz0_c1
		puz0_c2 => uz0_c2
		puz0_c3 => uz0_c3


		!******************************************************************
	!allocate wavefiled  usr, urr
		allocate(usr1%u(-npad+1:nnz+npad, -npad+1:nnx+npad), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about usr1, stop!!!"
			stop
		endif

		allocate(usr2%u(-npad+1:nnz+npad, -npad+1:nnx+npad), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about usr2, stop!!!"
			stop
		endif

		allocate(usr3%u(-npad+1:nnz+npad, -npad+1:nnx+npad), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about usr3, stop!!!"
			stop
		endif

		allocate(urr1%u(-npad+1:nnz+npad, -npad+1:nnx+npad), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about urr1, stop!!!"
			stop
		endif

		allocate(urr2%u(-npad+1:nnz+npad, -npad+1:nnx+npad), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about urr2, stop!!!"
			stop
		endif

		allocate(urr3%u(-npad+1:nnz+npad, -npad+1:nnx+npad), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about urr3, stop!!!"
			stop
		endif


		!******************************************************************
		!allocate pml
		allocate(uxr_a1%u(nnz,nx_bound_l+nx_bound_r), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about uxr_a1, stop!!!"
			stop
		endif

		allocate(uxr_a2%u(nnz,nx_bound_l+nx_bound_r), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about uxr_a2, stop!!!"
			stop
		endif

		allocate(uxr_a3%u(nnz,nx_bound_l+nx_bound_r), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about uxr_a3, stop!!!"
			stop
		endif

		allocate(uxr_b1%u(nnz,nx_bound_l+nx_bound_r), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about uxr_b1, stop!!!"
			stop
		endif

		allocate(uxr_b2%u(nnz,nx_bound_l+nx_bound_r), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about uxr_b2, stop!!!"
			stop
		endif

		allocate(uxr_bp1%u(nnz,nx_bound_l+nx_bound_r), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about uxr_bp1, stop!!!"
			stop
		endif

		allocate(uxr_bp2%u(nnz,nx_bound_l+nx_bound_r), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about uxr_bp2, stop!!!"
			stop
		endif

		allocate(uxr_bp3%u(nnz,nx_bound_l+nx_bound_r), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about uxr_bp3, stop!!!"
			stop
		endif

		allocate(uxr_c1%u(nnz,nx_bound_l+nx_bound_r), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about uxr_c1, stop!!!"
			stop
		endif

		allocate(uxr_c2%u(nnz,nx_bound_l+nx_bound_r), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about uxr_c2, stop!!!"
			stop
		endif

		allocate(uxr_c3%u(nnz,nx_bound_l+nx_bound_r), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about uxr_c3, stop!!!"
			stop
		endif


		!z direction
		allocate(uzr_a1%u(nz_bound_u+nz_bound_d,nnx), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about uzr_a1, stop!!!"
			stop
		endif

		allocate(uzr_a2%u(nz_bound_u+nz_bound_d,nnx), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about uzr_a2, stop!!!"
			stop
		endif

		allocate(uzr_a3%u(nz_bound_u+nz_bound_d,nnx), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about uzr_a3, stop!!!"
			stop
		endif

		allocate(uzr_b1%u(nz_bound_u+nz_bound_d,nnx), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about uzr_b1, stop!!!"
			stop
		endif

		allocate(uzr_b2%u(nz_bound_u+nz_bound_d,nnx), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about uzr_b2, stop!!!"
			stop
		endif

		allocate(uzr_bp1%u(nz_bound_u+nz_bound_d,nnx), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about uzr_bp1, stop!!!"
			stop
		endif

		allocate(uzr_bp2%u(nz_bound_u+nz_bound_d,nnx), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about uzr_bp2, stop!!!"
			stop
		endif

		allocate(uzr_bp3%u(nz_bound_u+nz_bound_d,nnx), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about uzr_bp3, stop!!!"
			stop
		endif

		allocate(uzr_c1%u(nz_bound_u+nz_bound_d,nnx), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about uzr_c1, stop!!!"
			stop
		endif

		allocate(uzr_c2%u(nz_bound_u+nz_bound_d,nnx), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about uzr_c2, stop!!!"
			stop
		endif

		allocate(uzr_c3%u(nz_bound_u+nz_bound_d,nnx), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about uzr_c3, stop!!!"
			stop
		endif

		!******************************************************************
		!wavefields
		usr1%u=0.0
		usr2%u=0.0
		usr3%u=0.0
		urr1%u=0.0
		urr2%u=0.0
		urr3%u=0.0

		!x direction
		uxr_a1%u=0.0
		uxr_a2%u=0.0
		uxr_a3%u=0.0
		uxr_b1%u=0.0
		uxr_b2%u=0.0
		uxr_bp1%u=0.0
		uxr_bp2%u=0.0
		uxr_bp3%u=0.0
		uxr_c1%u=0.0
		uxr_c2%u=0.0
		uxr_c3%u=0.0

		!z direction
		uzr_a1%u=0.0
		uzr_a2%u=0.0
		uzr_a3%u=0.0
		uzr_b1%u=0.0
		uzr_b2%u=0.0
		uzr_bp1%u=0.0
		uzr_bp2%u=0.0
		uzr_bp3%u=0.0
		uzr_c1%u=0.0
		uzr_c2%u=0.0
		uzr_c3%u=0.0

		!******************************************************************
		!wavefields
		pusr1 => usr1
		pusr2 => usr2
		pusr3 => usr3
		purr1 => urr1
		purr2 => urr2
		purr3 => urr3

		!x direction
		puxr_a1 => uxr_a1
		puxr_a2 => uxr_a2
		puxr_a3 => uxr_a3
		puxr_b1 => uxr_b1
		puxr_b2 => uxr_b2
		puxr_bp1 => uxr_bp1
		puxr_bp2 => uxr_bp2
		puxr_bp3 => uxr_bp3
		puxr_c1 => uxr_c1
		puxr_c2 => uxr_c2
		puxr_c3 => uxr_c3

		!z direction
		puzr_a1 => uzr_a1
		puzr_a2 => uzr_a2
		puzr_a3 => uzr_a3
		puzr_b1 => uzr_b1
		puzr_b2 => uzr_b2
		puzr_bp1 => uzr_bp1
		puzr_bp2 => uzr_bp2
		puzr_bp3 => uzr_bp3
		puzr_c1 => uzr_c1
		puzr_c2 => uzr_c2
		puzr_c3 => uzr_c3


		!******************************************************************

		allocate(uv_r01%u(-npad+1:nnz+npad, -npad+1:nnx+npad), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about uv_r01, stop!!!"
			stop
		endif

		allocate(uv_r02%u(-npad+1:nnz+npad, -npad+1:nnx+npad), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about uv_r02, stop!!!"
			stop
		endif

		allocate(uv_r03%u(-npad+1:nnz+npad, -npad+1:nnx+npad), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about uv_r03, stop!!!"
			stop
		endif

		!******************************************************************
		!allocate pml
		!x direction
		allocate(uv_x0_a1%u(nnz,nx_bound_l+nx_bound_r), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about uv_x0_a1, stop!!!"
			stop
		endif

		allocate(uv_x0_a2%u(nnz,nx_bound_l+nx_bound_r), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about uv_x0_a2, stop!!!"
			stop
		endif

		allocate(uv_x0_a3%u(nnz,nx_bound_l+nx_bound_r), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about uv_x0_a3, stop!!!"
			stop
		endif

		allocate(uv_x0_b1%u(nnz,nx_bound_l+nx_bound_r), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about uv_x0_b1, stop!!!"
			stop
		endif

		allocate(uv_x0_b2%u(nnz,nx_bound_l+nx_bound_r), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about uv_x0_b2, stop!!!"
			stop
		endif

		allocate(uv_x0_bp1%u(nnz,nx_bound_l+nx_bound_r), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about uv_x0_bp1, stop!!!"
			stop
		endif

		allocate(uv_x0_bp2%u(nnz,nx_bound_l+nx_bound_r), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about uv_x0_bp2, stop!!!"
			stop
		endif

		allocate(uv_x0_bp3%u(nnz,nx_bound_l+nx_bound_r), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about uv_x0_bp3, stop!!!"
			stop
		endif

		allocate(uv_x0_c1%u(nnz,nx_bound_l+nx_bound_r), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about uv_x0_c1, stop!!!"
			stop
		endif

		allocate(uv_x0_c2%u(nnz,nx_bound_l+nx_bound_r), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about uv_x0_c2, stop!!!"
			stop
		endif

		allocate(uv_x0_c3%u(nnz,nx_bound_l+nx_bound_r), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about uv_x0_c3, stop!!!"
			stop
		endif


		!z direction
		allocate(uv_z0_a1%u(nz_bound_u+nz_bound_d,nnx), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about uv_z0_a1, stop!!!"
			stop
		endif

		allocate(uv_z0_a2%u(nz_bound_u+nz_bound_d,nnx), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about uv_z0_a2, stop!!!"
			stop
		endif

		allocate(uv_z0_a3%u(nz_bound_u+nz_bound_d,nnx), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about uv_z0_a3, stop!!!"
			stop
		endif

		allocate(uv_z0_b1%u(nz_bound_u+nz_bound_d,nnx), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about uv_z0_b1, stop!!!"
			stop
		endif

		allocate(uv_z0_b2%u(nz_bound_u+nz_bound_d,nnx), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about uv_z0_b2, stop!!!"
			stop
		endif

		allocate(uv_z0_bp1%u(nz_bound_u+nz_bound_d,nnx), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about uv_z0_bp1, stop!!!"
			stop
		endif

		allocate(uv_z0_bp2%u(nz_bound_u+nz_bound_d,nnx), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about uv_z0_bp2, stop!!!"
			stop
		endif

		allocate(uv_z0_bp3%u(nz_bound_u+nz_bound_d,nnx), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about uv_z0_bp3, stop!!!"
			stop
		endif

		allocate(uv_z0_c1%u(nz_bound_u+nz_bound_d,nnx), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about uv_z0_c1, stop!!!"
			stop
		endif

		allocate(uv_z0_c2%u(nz_bound_u+nz_bound_d,nnx), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about uv_z0_c2, stop!!!"
			stop
		endif

		allocate(uv_z0_c3%u(nz_bound_u+nz_bound_d,nnx), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about uv_z0_c3, stop!!!"
			stop
		endif

		!******************************************************************
		!wavefields

		uv_r01%u=0.0
		uv_r02%u=0.0
		uv_r03%u=0.0

		!x direction
		uv_x0_a1%u=0.0
		uv_x0_a2%u=0.0
		uv_x0_a3%u=0.0
		uv_x0_b1%u=0.0
		uv_x0_b2%u=0.0
		uv_x0_bp1%u=0.0
		uv_x0_bp2%u=0.0
		uv_x0_bp3%u=0.0
		uv_x0_c1%u=0.0
		uv_x0_c2%u=0.0
		uv_x0_c3%u=0.0

		!z direction
		uv_z0_a1%u=0.0
		uv_z0_a2%u=0.0
		uv_z0_a3%u=0.0
		uv_z0_b1%u=0.0
		uv_z0_b2%u=0.0
		uv_z0_bp1%u=0.0
		uv_z0_bp2%u=0.0
		uv_z0_bp3%u=0.0
		uv_z0_c1%u=0.0
		uv_z0_c2%u=0.0
		uv_z0_c3%u=0.0

		!******************************************************************
		!wavefields

		puv_r01 => uv_r01
		puv_r02 => uv_r02
		puv_r03 => uv_r03

		!x direction
		puv_x0_a1 => uv_x0_a1
		puv_x0_a2 => uv_x0_a2
		puv_x0_a3 => uv_x0_a3
		puv_x0_b1 => uv_x0_b1
		puv_x0_b2 => uv_x0_b2
		puv_x0_bp1 => uv_x0_bp1
		puv_x0_bp2 => uv_x0_bp2
		puv_x0_bp3 => uv_x0_bp3
		puv_x0_c1 => uv_x0_c1
		puv_x0_c2 => uv_x0_c2
		puv_x0_c3 => uv_x0_c3

		!z direction
		puv_z0_a1 => uv_z0_a1
		puv_z0_a2 => uv_z0_a2
		puv_z0_a3 => uv_z0_a3
		puv_z0_b1 => uv_z0_b1
		puv_z0_b2 => uv_z0_b2
		puv_z0_bp1 => uv_z0_bp1
		puv_z0_bp2 => uv_z0_bp2
		puv_z0_bp3 => uv_z0_bp3
		puv_z0_c1 => uv_z0_c1
		puv_z0_c2 => uv_z0_c2
		puv_z0_c3 => uv_z0_c3

		!******************************************************************
	!allocate wavefiled  uv_sr, uv_rr
		allocate(uv_sr1%u(-npad+1:nnz+npad, -npad+1:nnx+npad), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about uv_sr1, stop!!!"
			stop
		endif

		allocate(uv_sr2%u(-npad+1:nnz+npad, -npad+1:nnx+npad), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about uv_sr2, stop!!!"
			stop
		endif

		allocate(uv_sr3%u(-npad+1:nnz+npad, -npad+1:nnx+npad), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about uv_sr3, stop!!!"
			stop
		endif

		allocate(uv_rr1%u(-npad+1:nnz+npad, -npad+1:nnx+npad), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about uv_rr1, stop!!!"
			stop
		endif

		allocate(uv_rr2%u(-npad+1:nnz+npad, -npad+1:nnx+npad), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about uv_rr2, stop!!!"
			stop
		endif

		allocate(uv_rr3%u(-npad+1:nnz+npad, -npad+1:nnx+npad), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about uv_rr3, stop!!!"
			stop
		endif

		!******************************************************************
		!allocate pml
		allocate(uv_xr_a1%u(nnz,nx_bound_l+nx_bound_r), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about uv_xr_a1, stop!!!"
			stop
		endif

		allocate(uv_xr_a2%u(nnz,nx_bound_l+nx_bound_r), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about uv_xr_a2, stop!!!"
			stop
		endif

		allocate(uv_xr_a3%u(nnz,nx_bound_l+nx_bound_r), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about uv_xr_a3, stop!!!"
			stop
		endif

		allocate(uv_xr_b1%u(nnz,nx_bound_l+nx_bound_r), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about uv_xr_b1, stop!!!"
			stop
		endif

		allocate(uv_xr_b2%u(nnz,nx_bound_l+nx_bound_r), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about uv_xr_b2, stop!!!"
			stop
		endif

		allocate(uv_xr_bp1%u(nnz,nx_bound_l+nx_bound_r), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about uv_xr_bp1, stop!!!"
			stop
		endif

		allocate(uv_xr_bp2%u(nnz,nx_bound_l+nx_bound_r), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about uv_xr_bp2, stop!!!"
			stop
		endif

		allocate(uv_xr_bp3%u(nnz,nx_bound_l+nx_bound_r), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about uv_xr_bp3, stop!!!"
			stop
		endif

		allocate(uv_xr_c1%u(nnz,nx_bound_l+nx_bound_r), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about uv_xr_c1, stop!!!"
			stop
		endif

		allocate(uv_xr_c2%u(nnz,nx_bound_l+nx_bound_r), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about uv_xr_c2, stop!!!"
			stop
		endif

		allocate(uv_xr_c3%u(nnz,nx_bound_l+nx_bound_r), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about uv_xr_c3, stop!!!"
			stop
		endif


		!z direction
		allocate(uv_zr_a1%u(nz_bound_u+nz_bound_d,nnx), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about uv_zr_a1, stop!!!"
			stop
		endif

		allocate(uv_zr_a2%u(nz_bound_u+nz_bound_d,nnx), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about uv_zr_a2, stop!!!"
			stop
		endif

		allocate(uv_zr_a3%u(nz_bound_u+nz_bound_d,nnx), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about uv_zr_a3, stop!!!"
			stop
		endif

		allocate(uv_zr_b1%u(nz_bound_u+nz_bound_d,nnx), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about uv_zr_b1, stop!!!"
			stop
		endif

		allocate(uv_zr_b2%u(nz_bound_u+nz_bound_d,nnx), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about uv_zr_b2, stop!!!"
			stop
		endif

		allocate(uv_zr_bp1%u(nz_bound_u+nz_bound_d,nnx), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about uv_zr_bp1, stop!!!"
			stop
		endif

		allocate(uv_zr_bp2%u(nz_bound_u+nz_bound_d,nnx), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about uv_zr_bp2, stop!!!"
			stop
		endif

		allocate(uv_zr_bp3%u(nz_bound_u+nz_bound_d,nnx), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about uv_zr_bp3, stop!!!"
			stop
		endif

		allocate(uv_zr_c1%u(nz_bound_u+nz_bound_d,nnx), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about uv_zr_c1, stop!!!"
			stop
		endif

		allocate(uv_zr_c2%u(nz_bound_u+nz_bound_d,nnx), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about uv_zr_c2, stop!!!"
			stop
		endif

		allocate(uv_zr_c3%u(nz_bound_u+nz_bound_d,nnx), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about uv_zr_c3, stop!!!"
			stop
		endif

		!******************************************************************
		!wavefields
		uv_sr1%u=0.0
		uv_sr2%u=0.0
		uv_sr3%u=0.0

		uv_rr1%u=0.0
		uv_rr2%u=0.0
		uv_rr3%u=0.0

		!x direction
		uv_xr_a1%u=0.0
		uv_xr_a2%u=0.0
		uv_xr_a3%u=0.0
		uv_xr_b1%u=0.0
		uv_xr_b2%u=0.0
		uv_xr_bp1%u=0.0
		uv_xr_bp2%u=0.0
		uv_xr_bp3%u=0.0
		uv_xr_c1%u=0.0
		uv_xr_c2%u=0.0
		uv_xr_c3%u=0.0

		!z direction
		uv_zr_a1%u=0.0
		uv_zr_a2%u=0.0
		uv_zr_a3%u=0.0
		uv_zr_b1%u=0.0
		uv_zr_b2%u=0.0
		uv_zr_bp1%u=0.0
		uv_zr_bp2%u=0.0
		uv_zr_bp3%u=0.0
		uv_zr_c1%u=0.0
		uv_zr_c2%u=0.0
		uv_zr_c3%u=0.0

		!******************************************************************
		!wavefields
		puv_sr1 => uv_sr1
		puv_sr2 => uv_sr2
		puv_sr3 => uv_sr3

		puv_rr1 => uv_rr1
		puv_rr2 => uv_rr2
		puv_rr3 => uv_rr3

		!x direction
		puv_xr_a1 => uv_xr_a1
		puv_xr_a2 => uv_xr_a2
		puv_xr_a3 => uv_xr_a3
		puv_xr_b1 => uv_xr_b1
		puv_xr_b2 => uv_xr_b2
		puv_xr_bp1 => uv_xr_bp1
		puv_xr_bp2 => uv_xr_bp2
		puv_xr_bp3 => uv_xr_bp3
		puv_xr_c1 => uv_xr_c1
		puv_xr_c2 => uv_xr_c2
		puv_xr_c3 => uv_xr_c3

		!z direction
		puv_zr_a1 => uv_zr_a1
		puv_zr_a2 => uv_zr_a2
		puv_zr_a3 => uv_zr_a3
		puv_zr_b1 => uv_zr_b1
		puv_zr_b2 => uv_zr_b2
		puv_zr_bp1 => uv_zr_bp1
		puv_zr_bp2 => uv_zr_bp2
		puv_zr_bp3 => uv_zr_bp3
		puv_zr_c1 => uv_zr_c1
		puv_zr_c2 => uv_zr_c2
		puv_zr_c3 => uv_zr_c3

		!******************************************************************
		if(iflag_src .eq. 1)then
			allocate(top0%u(order, nx_with_apert, ns), stat=ierr)
			if(ierr.ne.0)then
				write(*,*)"Can not allocate working memory about top0, stop!!!"
				stop
			endif

			allocate(bot0%u(order, nx_with_apert, ns), stat=ierr)
			if(ierr.ne.0)then
				write(*,*)"Can not allocate working memory about bot0, stop!!!"
				stop
			endif

			allocate(lef0%u(order, nz_with_apert, ns), stat=ierr)
			if(ierr.ne.0)then
				write(*,*)"Can not allocate working memory about lef0, stop!!!"
				stop
			endif

			allocate(rig0%u(order, nz_with_apert, ns), stat=ierr)
			if(ierr.ne.0)then
				write(*,*)"Can not allocate working memory about rig0, stop!!!"
				stop
			endif

			allocate(topr%u(order, nx_with_apert, ns), stat=ierr)
			if(ierr.ne.0)then
				write(*,*)"Can not allocate working memory about topr, stop!!!"
				stop
			endif

			allocate(botr%u(order, nx_with_apert, ns), stat=ierr)
			if(ierr.ne.0)then
				write(*,*)"Can not allocate working memory about botr, stop!!!"
				stop
			endif

			allocate(lefr%u(order, nz_with_apert, ns), stat=ierr)
			if(ierr.ne.0)then
				write(*,*)"Can not allocate working memory about lefr, stop!!!"
				stop
			endif

			allocate(rigr%u(order, nz_with_apert, ns), stat=ierr)
			if(ierr.ne.0)then
				write(*,*)"Can not allocate working memory about rigr, stop!!!"
				stop
			endif

			allocate(topv_r%u(order, nx_with_apert, ns), stat=ierr)
			if(ierr.ne.0)then
				write(*,*)"Can not allocate working memory about topv_r, stop!!!"
				stop
			endif

			allocate(botv_r%u(order, nx_with_apert, ns), stat=ierr)
			if(ierr.ne.0)then
				write(*,*)"Can not allocate working memory about botv_r, stop!!!"
				stop
			endif

			allocate(lefv_r%u(order, nz_with_apert, ns), stat=ierr)
			if(ierr.ne.0)then
				write(*,*)"Can not allocate working memory about lefv_r, stop!!!"
				stop
			endif

			allocate(rigv_r%u(order, nz_with_apert, ns), stat=ierr)
			if(ierr.ne.0)then
				write(*,*)"Can not allocate working memory about rigv_r, stop!!!"
				stop
			endif

		else if(iflag_src .eq. 0)then
			allocate(ust0%u(-npad+1:nnz+npad, -npad+1:nnx+npad, ns), stat=ierr)
			if(ierr.ne.0)then
				write(*,*)"Can not allocate working memory about ust0, stop!!!"
				stop
			endif

			allocate(ustr%u(-npad+1:nnz+npad, -npad+1:nnx+npad, ns), stat=ierr)
			if(ierr.ne.0)then
				write(*,*)"Can not allocate working memory about ustr, stop!!!"
				stop
			endif


			allocate(ustv_r%u(-npad+1:nnz+npad, -npad+1:nnx+npad, ns), stat=ierr)
			if(ierr.ne.0)then
				write(*,*)"Can not allocate working memory about ustv_r, stop!!!"
				stop
			endif

		endif

		if(iflag_src .eq. 1)then
			top0%u=0.0
			bot0%u=0.0
			lef0%u=0.0
			rig0%u=0.0
			topr%u=0.0
			botr%u=0.0
			lefr%u=0.0
			rigr%u=0.0

			topv_r%u=0.0
			botv_r%u=0.0
			lefv_r%u=0.0
			rigv_r%u=0.0

		else if(iflag_src .eq. 0)then
			ust0%u=0.0
			ustr%u=0.0

			ustv_r%u=0.0
		endif

		!******************************************************************
		allocate(image_tmp(nz_with_apert, nx_with_apert), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about uz_c3, stop!!!"
			stop
		endif
		allocate(image_ms(nz_with_apert, nx_with_apert), &
				 image_dm(nz_with_apert, nx_with_apert), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"(image_ms, image_dm), Can not allocate working memory, stop!!!"
			stop
		endif
		allocate(image_rv_ms(nz_with_apert, nx_with_apert), &
				 image_rv_dm(nz_with_apert, nx_with_apert), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"(image_rv_ms, image_rv_dm), Can not allocate working memory, stop!!!"
			stop
		endif

		image_tmp=0.0
		image_ms = 0.0
		image_dm = 0.0
		image_rv_ms = 0.0
		image_rv_dm = 0.0

		!******************************************************************
		allocate(coe_1st(npad), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about coe_1st, stop!!!"
			stop
		endif
		allocate(coe_2nd(npad), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about coe_2nd, stop!!!"
			stop
		endif
		allocate(coe_1st_dx(npad), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about coe_1st_dx, stop!!!"
			stop
		endif
		allocate(coe_1st_dz(npad), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about coe_1st_dz, stop!!!"
			stop
		endif
		allocate(coe_2nd_dx2(npad), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about coe_2nd_dx2, stop!!!"
			stop
		endif
		allocate(coe_2nd_dz2(npad), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about coe_2nd_dz2, stop!!!"
			stop
		endif
		allocate(funa(nnz, nnx), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about funa, stop!!!"
			stop
		endif
		allocate(dfuna(nnz, nnx), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about dfuna, stop!!!"
			stop
		endif
		allocate(line_xl(nx_bound_l, 2), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about line_xl, stop!!!"
			stop
		endif
		allocate(line_xr(nx_bound_r, 2), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about line_xr, stop!!!"
			stop
		endif
		allocate(line_zu(nz_bound_u, 2), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about line_zu, stop!!!"
			stop
		endif
		allocate(line_zd(nz_bound_d, 2), stat=ierr)
		if(ierr.ne.0)then
			write(*,*)"Can not allocate working memory about line_zd, stop!!!"
			stop
		endif

		coe_1st=0.0
		coe_2nd=0.0
		coe_1st_dx=0.0
		coe_1st_dz=0.0
		coe_2nd_dx2=0.0
		coe_2nd_dz2=0.0
		funa=0.0
		dfuna=0.0
		line_xl=0
		line_xr=0
		line_zu=0
		line_zd=0

	end subroutine

!***********************************************************************
	subroutine rv_born_2d_single_shot_extrapolation_deallocate(iflag_src)

		use module_rv_born_extrapolation_2d
		implicit none
		integer		::	iflag_src

		!deallocate wavefiled  us0, ur0
		deallocate(us01%u)
		deallocate(us02%u)
		deallocate(us03%u)
		deallocate(ur01%u)
		deallocate(ur02%u)
		deallocate(ur03%u)

		!deallocate pml
		!x direction
		deallocate(ux0_a1%u)
		deallocate(ux0_a2%u)
		deallocate(ux0_a3%u)
		deallocate(ux0_b1%u)
		deallocate(ux0_b2%u)
		deallocate(ux0_bp1%u)
		deallocate(ux0_bp2%u)
		deallocate(ux0_bp3%u)
		deallocate(ux0_c1%u)
		deallocate(ux0_c2%u)
		deallocate(ux0_c3%u)

		!z direction
		deallocate(uz0_a1%u)
		deallocate(uz0_a2%u)
		deallocate(uz0_a3%u)
		deallocate(uz0_b1%u)
		deallocate(uz0_b2%u)
		deallocate(uz0_bp1%u)
		deallocate(uz0_bp2%u)
		deallocate(uz0_bp3%u)
		deallocate(uz0_c1%u)
		deallocate(uz0_c2%u)
		deallocate(uz0_c3%u)

		!deallocate wavefiled usr, urr
		deallocate(usr1%u)
		deallocate(usr2%u)
		deallocate(usr3%u)
		deallocate(urr1%u)
		deallocate(urr2%u)
		deallocate(urr3%u)

		!deallocate pml
		!x direction
		deallocate(uxr_a1%u)
		deallocate(uxr_a2%u)
		deallocate(uxr_a3%u)
		deallocate(uxr_b1%u)
		deallocate(uxr_b2%u)
		deallocate(uxr_bp1%u)
		deallocate(uxr_bp2%u)
		deallocate(uxr_bp3%u)
		deallocate(uxr_c1%u)
		deallocate(uxr_c2%u)
		deallocate(uxr_c3%u)

		!z direction
		deallocate(uzr_a1%u)
		deallocate(uzr_a2%u)
		deallocate(uzr_a3%u)
		deallocate(uzr_b1%u)
		deallocate(uzr_b2%u)
		deallocate(uzr_bp1%u)
		deallocate(uzr_bp2%u)
		deallocate(uzr_bp3%u)
		deallocate(uzr_c1%u)
		deallocate(uzr_c2%u)
		deallocate(uzr_c3%u)

		deallocate(uv_r01%u)
		deallocate(uv_r02%u)
		deallocate(uv_r03%u)

		!x direction
		deallocate(uv_x0_a1%u)
		deallocate(uv_x0_a2%u)
		deallocate(uv_x0_a3%u)
		deallocate(uv_x0_b1%u)
		deallocate(uv_x0_b2%u)
		deallocate(uv_x0_bp1%u)
		deallocate(uv_x0_bp2%u)
		deallocate(uv_x0_bp3%u)
		deallocate(uv_x0_c1%u)
		deallocate(uv_x0_c2%u)
		deallocate(uv_x0_c3%u)

		!z direction
		deallocate(uv_z0_a1%u)
		deallocate(uv_z0_a2%u)
		deallocate(uv_z0_a3%u)
		deallocate(uv_z0_b1%u)
		deallocate(uv_z0_b2%u)
		deallocate(uv_z0_bp1%u)
		deallocate(uv_z0_bp2%u)
		deallocate(uv_z0_bp3%u)
		deallocate(uv_z0_c1%u)
		deallocate(uv_z0_c2%u)
		deallocate(uv_z0_c3%u)

		!deallocate wavefiled uv_sr, uv_rr
		deallocate(uv_sr1%u)
		deallocate(uv_sr2%u)
		deallocate(uv_sr3%u)
		deallocate(uv_rr1%u)
		deallocate(uv_rr2%u)
		deallocate(uv_rr3%u)

		!deallocate pml
		!x direction
		deallocate(uv_xr_a1%u)
		deallocate(uv_xr_a2%u)
		deallocate(uv_xr_a3%u)
		deallocate(uv_xr_b1%u)
		deallocate(uv_xr_b2%u)
		deallocate(uv_xr_bp1%u)
		deallocate(uv_xr_bp2%u)
		deallocate(uv_xr_bp3%u)
		deallocate(uv_xr_c1%u)
		deallocate(uv_xr_c2%u)
		deallocate(uv_xr_c3%u)

		!z direction
		deallocate(uv_zr_a1%u)
		deallocate(uv_zr_a2%u)
		deallocate(uv_zr_a3%u)
		deallocate(uv_zr_b1%u)
		deallocate(uv_zr_b2%u)
		deallocate(uv_zr_bp1%u)
		deallocate(uv_zr_bp2%u)
		deallocate(uv_zr_bp3%u)
		deallocate(uv_zr_c1%u)
		deallocate(uv_zr_c2%u)
		deallocate(uv_zr_c3%u)


		if(iflag_src .eq. 1)then
			deallocate(top0%u)
			deallocate(bot0%u)
			deallocate(lef0%u)
			deallocate(rig0%u)
			deallocate(topr%u)
			deallocate(botr%u)
			deallocate(lefr%u)
			deallocate(rigr%u)

			deallocate(topv_r%u)
			deallocate(botv_r%u)
			deallocate(lefv_r%u)
			deallocate(rigv_r%u)
		else if(iflag_src .eq. 0)then
			deallocate(ust0%u)
			deallocate(ustr%u)

			deallocate(ustv_r%u)
		endif

		deallocate(image_tmp)
		deallocate(image_ms, image_dm)
		deallocate(image_rv_ms, image_rv_dm)

		deallocate(coe_1st)
		deallocate(coe_2nd)
		deallocate(coe_1st_dx)
		deallocate(coe_1st_dz)
		deallocate(coe_2nd_dx2)
		deallocate(coe_2nd_dz2)
		deallocate(funa)
		deallocate(dfuna)
		deallocate(line_xl)
		deallocate(line_xr)
		deallocate(line_zu)
		deallocate(line_zd)


	end subroutine

!***********************************************************************
	subroutine rv_born_extrapolation_2d_forward_variable_update()

		use module_rv_born_extrapolation_2d
		implicit none
		!Data dictionary: declare variable types, definitons, & units
		!Dummy Variables
	
		!Local variables
		integer		::	ierr
		integer		::	ishot

		!===========================================================
	!update wavefield for us0
		pt     => pus01
		pus01  => pus02
		pus02  => pus03
		pus03  => pt

		!x direction
		pt      => pux0_a1
		pux0_a1 => pux0_a2
		pux0_a2 => pux0_a3
		pux0_a3 => pt

		pt		 => pux0_bp1
		pux0_bp1 => pux0_bp2
		pux0_bp2 => pux0_bp3
		pux0_bp3 => pt

		pt	    => pux0_b1
		pux0_b1 => pux0_b2
		pux0_b2 => pt

		pt	    => pux0_c1
		pux0_c1 => pux0_c2
		pux0_c2 => pux0_c3
		pux0_c3 => pt


		!z direction
		pt      => puz0_a1
		puz0_a1 => puz0_a2
		puz0_a2 => puz0_a3
		puz0_a3 => pt

		pt		 => puz0_bp1
		puz0_bp1 => puz0_bp2
		puz0_bp2 => puz0_bp3
		puz0_bp3 => pt

		pt	    => puz0_b1
		puz0_b1 => puz0_b2
		puz0_b2 => pt

		pt	    => puz0_c1
		puz0_c1 => puz0_c2
		puz0_c2 => puz0_c3
		puz0_c3 => pt


		!===========================================================
	!update wavefield for usr
		pt     => pusr1
		pusr1  => pusr2
		pusr2  => pusr3
		pusr3  => pt

		!x direction
		pt      => puxr_a1
		puxr_a1 => puxr_a2
		puxr_a2 => puxr_a3
		puxr_a3 => pt

		pt		 => puxr_bp1
		puxr_bp1 => puxr_bp2
		puxr_bp2 => puxr_bp3
		puxr_bp3 => pt

		pt	    => puxr_b1
		puxr_b1 => puxr_b2
		puxr_b2 => pt

		pt	    => puxr_c1
		puxr_c1 => puxr_c2
		puxr_c2 => puxr_c3
		puxr_c3 => pt


		!z direction
		pt      => puzr_a1
		puzr_a1 => puzr_a2
		puzr_a2 => puzr_a3
		puzr_a3 => pt

		pt		 => puzr_bp1
		puzr_bp1 => puzr_bp2
		puzr_bp2 => puzr_bp3
		puzr_bp3 => pt

		pt	    => puzr_b1
		puzr_b1 => puzr_b2
		puzr_b2 => pt

		pt	    => puzr_c1
		puzr_c1 => puzr_c2
		puzr_c2 => puzr_c3
		puzr_c3 => pt

		pt    		=> puv_sr1
		puv_sr1		=> puv_sr2
		puv_sr2		=> puv_sr3
		puv_sr3		=> pt

		!x direction
		pt      	=> puv_xr_a1
		puv_xr_a1	=> puv_xr_a2
		puv_xr_a2	=> puv_xr_a3
		puv_xr_a3	=> pt

		pt		 	=> puv_xr_bp1
		puv_xr_bp1	=> puv_xr_bp2
		puv_xr_bp2	=> puv_xr_bp3
		puv_xr_bp3	=> pt

		pt	   		=> puv_xr_b1
		puv_xr_b1	=> puv_xr_b2
		puv_xr_b2	=> pt

		pt	    	=> puv_xr_c1
		puv_xr_c1	=> puv_xr_c2
		puv_xr_c2	=> puv_xr_c3
		puv_xr_c3	=> pt

		!z direction
		pt      	=> puv_zr_a1
		puv_zr_a1	=> puv_zr_a2
		puv_zr_a2	=> puv_zr_a3
		puv_zr_a3	=> pt

		pt		 	=> puv_zr_bp1
		puv_zr_bp1	=> puv_zr_bp2
		puv_zr_bp2	=> puv_zr_bp3
		puv_zr_bp3	=> pt

		pt	    	=> puv_zr_b1
		puv_zr_b1	=> puv_zr_b2
		puv_zr_b2	=> pt

		pt	    	=> puv_zr_c1
		puv_zr_c1	=> puv_zr_c2
		puv_zr_c2	=> puv_zr_c3
		puv_zr_c3	=> pt


	    return
    end	subroutine

!***********************************************************************
	subroutine rv_born_extrapolation_2d_backward_variable_update(iflag_src)

		use module_rv_born_extrapolation_2d
		implicit none
		!Data dictionary: declare variable types, definitons, & units
		!Dummy Variables
		integer		::	iflag_src
	
		!Local variables

		!===========================================================
	!update wavefield  ur0
		pt    => pur01
		pur01 => pur02
		pur02 => pur03
		pur03 => pt

		!x direction
		pt      => pux0_a1
		pux0_a1 => pux0_a2
		pux0_a2 => pux0_a3
		pux0_a3 => pt

		pt		 => pux0_bp1
		pux0_bp1 => pux0_bp2
		pux0_bp2 => pux0_bp3
		pux0_bp3 => pt

		pt	    => pux0_b1
		pux0_b1 => pux0_b2
		pux0_b2 => pt

		pt	    => pux0_c1
		pux0_c1 => pux0_c2
		pux0_c2 => pux0_c3
		pux0_c3 => pt


		!z direction
		pt      => puz0_a1
		puz0_a1 => puz0_a2
		puz0_a2 => puz0_a3
		puz0_a3 => pt

		pt		 => puz0_bp1
		puz0_bp1 => puz0_bp2
		puz0_bp2 => puz0_bp3
		puz0_bp3 => pt

		pt	    => puz0_b1
		puz0_b1 => puz0_b2
		puz0_b2 => pt

		pt	    => puz0_c1
		puz0_c1 => puz0_c2
		puz0_c2 => puz0_c3
		puz0_c3 => pt

		!*===============================================
	!update wavefield urr
		pt    => purr1
		purr1 => purr2
		purr2 => purr3
		purr3 => pt

		!x direction
		pt      => puxr_a1
		puxr_a1 => puxr_a2
		puxr_a2 => puxr_a3
		puxr_a3 => pt

		pt		 => puxr_bp1
		puxr_bp1 => puxr_bp2
		puxr_bp2 => puxr_bp3
		puxr_bp3 => pt

		pt	    => puxr_b1
		puxr_b1 => puxr_b2
		puxr_b2 => pt

		pt	    => puxr_c1
		puxr_c1 => puxr_c2
		puxr_c2 => puxr_c3
		puxr_c3 => pt


		!z direction
		pt      => puzr_a1
		puzr_a1 => puzr_a2
		puzr_a2 => puzr_a3
		puzr_a3 => pt

		pt		 => puzr_bp1
		puzr_bp1 => puzr_bp2
		puzr_bp2 => puzr_bp3
		puzr_bp3 => pt

		pt	    => puzr_b1
		puzr_b1 => puzr_b2
		puzr_b2 => pt

		pt	    => puzr_c1
		puzr_c1 => puzr_c2
		puzr_c2 => puzr_c3
		puzr_c3 => pt


		!*===============================================
	!update wavefield  uv_r0
		pt    	=> puv_r01
		puv_r01 => puv_r02
		puv_r02 => puv_r03
		puv_r03 => pt

		!x direction
		pt     		=> puv_x0_a1
		puv_x0_a1	=> puv_x0_a2
		puv_x0_a2	=> puv_x0_a3
		puv_x0_a3 	=> pt

		pt		 	=> puv_x0_bp1
		puv_x0_bp1	=> puv_x0_bp2
		puv_x0_bp2	=> puv_x0_bp3
		puv_x0_bp3	=> pt

		pt	    	=> puv_x0_b1
		puv_x0_b1	=> puv_x0_b2
		puv_x0_b2	=> pt

		pt	    	=> puv_x0_c1
		puv_x0_c1	=> puv_x0_c2
		puv_x0_c2	=> puv_x0_c3
		puv_x0_c3	=> pt

		!z direction
		pt      	=> puv_z0_a1
		puv_z0_a1	=> puv_z0_a2
		puv_z0_a2	=> puv_z0_a3
		puv_z0_a3	=> pt

		pt		 	=> puv_z0_bp1
		puv_z0_bp1	=> puv_z0_bp2
		puv_z0_bp2	=> puv_z0_bp3
		puv_z0_bp3	=> pt

		pt	    	=> puv_z0_b1
		puv_z0_b1	=> puv_z0_b2
		puv_z0_b2	=> pt

		pt	    	=> puv_z0_c1
		puv_z0_c1	=> puv_z0_c2
		puv_z0_c2	=> puv_z0_c3
		puv_z0_c3	=> pt

		!*===============================================
	!update wavefield uv_rr
		pt    		=> puv_rr1
		puv_rr1		=> puv_rr2
		puv_rr2		=> puv_rr3
		puv_rr3		=> pt

		!x direction
		pt      	=> puv_xr_a1
		puv_xr_a1	=> puv_xr_a2
		puv_xr_a2	=> puv_xr_a3
		puv_xr_a3	=> pt

		pt		 	=> puv_xr_bp1
		puv_xr_bp1	=> puv_xr_bp2
		puv_xr_bp2	=> puv_xr_bp3
		puv_xr_bp3	=> pt

		pt	   		=> puv_xr_b1
		puv_xr_b1	=> puv_xr_b2
		puv_xr_b2	=> pt

		pt	    	=> puv_xr_c1
		puv_xr_c1	=> puv_xr_c2
		puv_xr_c2	=> puv_xr_c3
		puv_xr_c3	=> pt


		!z direction
		pt      	=> puv_zr_a1
		puv_zr_a1	=> puv_zr_a2
		puv_zr_a2	=> puv_zr_a3
		puv_zr_a3	=> pt

		pt		 	=> puv_zr_bp1
		puv_zr_bp1	=> puv_zr_bp2
		puv_zr_bp2	=> puv_zr_bp3
		puv_zr_bp3	=> pt

		pt	    	=> puv_zr_b1
		puv_zr_b1	=> puv_zr_b2
		puv_zr_b2	=> pt

		pt	    	=> puv_zr_c1
		puv_zr_c1	=> puv_zr_c2
		puv_zr_c2	=> puv_zr_c3
		puv_zr_c3	=> pt


		!*===============================================
		if(iflag_src .eq. 1)then
			!update source wavefield
			pt    => pus01
			pus01 => pus02
			pus02 => pus03
			pus03 => pt

			pt    => pusr1
			pusr1 => pusr2
			pusr2 => pusr3
			pusr3 => pt

			pt    	=> puv_sr1
			puv_sr1 => puv_sr2
			puv_sr2 => puv_sr3
			puv_sr3 => pt

		endif

	    return
    end	subroutine

!***********************************************************************
	subroutine rv_born_extrapolation_2d_backward_Reinitiallization()

		use module_rv_born_extrapolation_2d
		implicit none
		!Data dictionary: declare variable types, definitons, & units
		!Dummy Variables
	
		!Local variables

		!===========================================================
		pt    => pus02
		pus02 => pus01
		pus01 => pt

		pt    => pusr2
		pusr2 => pusr1
		pusr1 => pt

		pt    	=> puv_sr2
		puv_sr2 => puv_sr1
		puv_sr1 => pt

		!===========================================================
		!wavefield ur0
		!x direction
	    pux0_a1%u=0.0
	    pux0_a2%u=0.0
	    pux0_a3%u=0.0
	    pux0_b1%u=0.0
	    pux0_b2%u=0.0
  		pux0_bp1%u=0.0
	  	pux0_bp2%u=0.0
	    pux0_bp3%u=0.0
	    pux0_c1%u=0.0
	    pux0_c2%u=0.0
	    pux0_c3%u=0.0

	    !z direction
	    puz0_a1%u=0.0
	    puz0_a2%u=0.0
	    puz0_a3%u=0.0
	    puz0_b1%u=0.0
	    puz0_b2%u=0.0
	    puz0_bp1%u=0.0
	    puz0_bp2%u=0.0
	    puz0_bp3%u=0.0
	    puz0_c1%u=0.0
	    puz0_c2%u=0.0
	    puz0_c3%u=0.0

		!wavefield urr
		!x direction
	    puxr_a1%u=0.0
	    puxr_a2%u=0.0
	    puxr_a3%u=0.0
	    puxr_b1%u=0.0
	    puxr_b2%u=0.0
	    puxr_bp1%u=0.0
	    puxr_bp2%u=0.0
	    puxr_bp3%u=0.0
	    puxr_c1%u=0.0
	    puxr_c2%u=0.0
	    puxr_c3%u=0.0

	    !z direction
	    puzr_a1%u=0.0
	    puzr_a2%u=0.0
	    puzr_a3%u=0.0
	    puzr_b1%u=0.0
	    puzr_b2%u=0.0
	    puzr_bp1%u=0.0
	    puzr_bp2%u=0.0
	    puzr_bp3%u=0.0
	    puzr_c1%u=0.0
	    puzr_c2%u=0.0
	    puzr_c3%u=0.0

		!wavefield uv_r0
		!x direction
	    puv_x0_a1%u=0.0
	    puv_x0_a2%u=0.0
	    puv_x0_a3%u=0.0
	    puv_x0_b1%u=0.0
	    puv_x0_b2%u=0.0
  		puv_x0_bp1%u=0.0
	  	puv_x0_bp2%u=0.0
	    puv_x0_bp3%u=0.0
	    puv_x0_c1%u=0.0
	    puv_x0_c2%u=0.0
	    puv_x0_c3%u=0.0

	    !z direction
	    puv_z0_a1%u=0.0
	    puv_z0_a2%u=0.0
	    puv_z0_a3%u=0.0
	    puv_z0_b1%u=0.0
	    puv_z0_b2%u=0.0
	    puv_z0_bp1%u=0.0
	    puv_z0_bp2%u=0.0
	    puv_z0_bp3%u=0.0
	    puv_z0_c1%u=0.0
	    puv_z0_c2%u=0.0
	    puv_z0_c3%u=0.0

		!wavefield uv_rr
		!x direction
	    puv_xr_a1%u=0.0
	    puv_xr_a2%u=0.0
	    puv_xr_a3%u=0.0
	    puv_xr_b1%u=0.0
	    puv_xr_b2%u=0.0
	    puv_xr_bp1%u=0.0
	    puv_xr_bp2%u=0.0
	    puv_xr_bp3%u=0.0
	    puv_xr_c1%u=0.0
	    puv_xr_c2%u=0.0
	    puv_xr_c3%u=0.0

	    !z direction
	    puv_zr_a1%u=0.0
	    puv_zr_a2%u=0.0
	    puv_zr_a3%u=0.0
	    puv_zr_b1%u=0.0
	    puv_zr_b2%u=0.0
	    puv_zr_bp1%u=0.0
	    puv_zr_bp2%u=0.0
	    puv_zr_bp3%u=0.0
	    puv_zr_c1%u=0.0
	    puv_zr_c2%u=0.0
	    puv_zr_c3%u=0.0

		!*wavefiled
		pur01%u=0.0
		pur02%u=0.0
		pur03%u=0.0
		purr1%u=0.0
		purr2%u=0.0
		purr3%u=0.0

		puv_r01%u=0.0
		puv_r02%u=0.0
		puv_r03%u=0.0
		puv_rr1%u=0.0
		puv_rr2%u=0.0
		puv_rr3%u=0.0


	    return
    end	subroutine

!***********************************************************************


