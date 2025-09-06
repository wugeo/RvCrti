!*=======================================================================*
	subroutine RV_CRTI_2D_image_single_shot(wavelet, vv, &
					shotobs, shotcal, shotres, &
					image_r_loc, image_i_sig, &
					gx, gz, index_recgx, index_recgz)

		use	global
		use RV_CRTI_2d_shot_parameter 
		use module_born_extrapolation_2d
		implicit none
		!Data dictionary: declare variable types, definitons, & units
		!Dummy Variables
		real	::	wavelet(lt)
		real	::	vv(nnz, nnx)
		real	::	shotobs(lt, ntr)
		real	::	shotcal(lt, ntr)
		real	::	shotres(lt, ntr)
		real	::	image_r_loc(nnz, nnx)
		real	::	image_i_sig(nz_with_apert, nx_with_apert)	
		real	::	gx(ntr), gz(ntr)
		integer	::	index_recgx(ntr), index_recgz(ntr)


		!Local Variables
		integer	::	ix, iz
		integer ::	it, iit, iii, itt
		integer ::	iix, inx
		integer ::	iiz, inz
		integer	::	itr
		real	::	vsurf
		real	::	v2
		real	::	misfit_sig_sum
		real	::	time1, time2
		real	::	tmp
		integer	::	icount
		integer	::	ii, iter
		real	::	offset_set
		integer	::	misfit_loc



		CALL Form_Coe_extrapolation_2d(coe, dx, dz, dt)

		CALL born_2d_single_shot_extrapolation_allocate(nnx, nnz, &
						nx_with_apert, nz_with_apert, &
						nx_bound_l, nx_bound_r, &
						nz_bound_u, nz_bound_d, &
						npad, order, ns, iflag_src)

		CALL coefficient_2nd(npad, coe_2nd)
		CALL coefficient_1st(npad, coe_1st)
		do ii=1, npad
			coe_2nd_dx2(ii)=coe_2nd(ii)/(dx*dx)
			coe_2nd_dz2(ii)=coe_2nd(ii)/(dz*dz)
			coe_1st_dx(ii)=coe_1st(ii)/dx
			coe_1st_dz(ii)=coe_1st(ii)/dz
		enddo

		CALL absorbing_function_2d(funa, dfuna, vv, nnz, nnx, dx, dz, &
					line_xl, line_xr, line_zu, line_zd, &
					nx_bound_l, nx_bound_r, &
					nz_bound_u, nz_bound_d, &
					r)



		do it=1, lt + nwt

			itt=it-nwt

			if(mod(itt, it_snap).eq. 1) then
				write(*,*)'1--forward: it=',itt
			endif

			CALL  Extrapolation_2D_Add_PML_One_Step(vv, pus01%u, pus02%u, pus03%u, &
						pux0_a1%u, pux0_a2%u, pux0_a3%u, &
						pux0_b1%u, pux0_b2%u, pux0_bp1%u, pux0_bp2%u, pux0_bp3%u, &
						pux0_c1%u, pux0_c2%u, pux0_c3%u, &
						puz0_a1%u, puz0_a2%u, puz0_a3%u, &
						puz0_b1%u, puz0_b2%u, puz0_bp1%u, puz0_bp2%u, puz0_bp3%u, &
						puz0_c1%u, puz0_c2%u, puz0_c3%u, &
						funa, dfuna, &
						line_xl, line_xr, line_zu, line_zd, &
						coe_2nd_dx2, coe_2nd_dz2, &
						coe_1st_dx, coe_1st_dz, &
						npad, nnx, nnz, nx, nz, &
						nx_bound_l, nx_bound_r, &
						nz_bound_u, nz_bound_d, &
						dx, dz, dt, r, itt)

			if(it .ge. 1 .and. it .le. lt)then
				v2=vv(ns_z, ns_x)*vv(ns_z, ns_x)
				pus03%u(ns_z, ns_x)=pus03%u(ns_z, ns_x)+ wavelet(it)*dt*dt*v2
			endif

			iflag_order = 1
			CALL source_wavefield_rtm2d(pus01%u, pus02%u, pus03%u, &
					top0%u, bot0%u, lef0%u, rig0%u, &
					nnx, nnz, nx_with_apert, nz_with_apert, &
					nx_bound_l, nx_bound_r, &
					nz_bound_u, nz_bound_d, &
					npad, order, lt, ns, &
					itt, it_step, &
					iflag_src, iflag_order, &
					isok)

			if(iflag_src .eq. 0)then
				if(itt .ge. 1 .and. itt .le. lt) then
					if(mod(itt, it_step) .eq. 0) then
						iit=itt/it_step
						ust0%u(:, :, iit)=pus03%u
					endif
				endif
			endif

			CALL  Extrapolation_2D_Add_PML_One_Step(vv, pusr1%u, pusr2%u, pusr3%u, &
						puxr_a1%u, puxr_a2%u, puxr_a3%u, &
						puxr_b1%u, puxr_b2%u, puxr_bp1%u, puxr_bp2%u, puxr_bp3%u, &
						puxr_c1%u, puxr_c2%u, puxr_c3%u, &
						puzr_a1%u, puzr_a2%u, puzr_a3%u, &
						puzr_b1%u, puzr_b2%u, puzr_bp1%u, puzr_bp2%u, puzr_bp3%u, &
						puzr_c1%u, puzr_c2%u, puzr_c3%u, &
						funa, dfuna, &
						line_xl, line_xr, line_zu, line_zd, &
						coe_2nd_dx2, coe_2nd_dz2, &
						coe_1st_dx, coe_1st_dz, &
						npad, nnx, nnz, nx, nz, &
						nx_bound_l, nx_bound_r, &
						nz_bound_u, nz_bound_d, &
						dx, dz, dt, r, itt)

			do inx=1, nnx
				do inz=1, nnz
					v2=vv(inz, inx)*vv(inz, inx)
					tmp=pus01%u(inz, inx) + pus03%u(inz, inx) -2.0*pus02%u(inz, inx)
					tmp=tmp/(v2*dt*dt)
					tmp=-1.0*image_r_loc(inz, inx)*tmp
					pusr3%u(inz, inx)= pusr3%u(inz, inx)+tmp
				enddo
			enddo

			iflag_order = 1
			CALL source_wavefield_rtm2d(pusr1%u, pusr2%u, pusr3%u, &
					topr%u, botr%u, lefr%u, rigr%u, &
					nnx, nnz, nx_with_apert, nz_with_apert, &
					nx_bound_l, nx_bound_r, &
					nz_bound_u, nz_bound_d, &
					npad, order, lt, ns, &
					itt, it_step, &
					iflag_src, iflag_order, &
					isok)

			if(iflag_src .eq. 0)then
				if(itt .ge. 1 .and. itt .le. lt) then
					if(mod(itt, it_step) .eq. 0) then
						iit=itt/it_step
						ustr%u(:, :, iit)=pusr3%u
					endif
				endif
			endif

			if(itt .gt. 0 .and. itt .le. lt)then
				do itr=1, ntr
					nr_x=index_recgx(itr)
					nr_z=index_recgz(itr)
					if(nr_x .ge. 1   .and. nr_z .ge. 1   .and. &
					   nr_x .le. nnx .and. nr_z .le. nnz )then
						shotcal(itt, itr)=pusr3%u(nr_z, nr_x)
					endif
				enddo
			endif

			CALL born_extrapolation_2d_forward_variable_update()

		enddo

		CALL rv_crti_2d_shotgathers_image(shotobs, shotcal, shotres, &
					vv, image_r_loc, misfit_loc, &
					sx, sz, gx, gz, index_recgx, index_recgz, ns_x, ns_z, &
                    nnz, nnx, lt, ntr, nwt, dt, dx, dz)


		CALL born_extrapolation_2d_backward_Reinitiallization()

	    do it=lt, 1, -1		

			if(mod(it, it_snap).eq. 1) then
				write(*,*)'1--backward: it=',it
			endif

			iflag_order = -1
			CALL source_wavefield_rtm2d(pus01%u, pus02%u, pus03%u, &
					top0%u, bot0%u, lef0%u, rig0%u, &
					nnx, nnz, nx_with_apert, nz_with_apert, &
					nx_bound_l, nx_bound_r, &
					nz_bound_u, nz_bound_d, &
					npad, order, lt, ns, &
					it, it_step, &
					iflag_src, iflag_order, &
					isok)

			iflag_order = -1
			CALL source_wavefield_rtm2d(pusr1%u, pusr2%u, pusr3%u, &
					topr%u, botr%u, lefr%u, rigr%u, &
					nnx, nnz, nx_with_apert, nz_with_apert, &
					nx_bound_l, nx_bound_r, &
					nz_bound_u, nz_bound_d, &
					npad, order, lt, ns, &
					it, it_step, &
					iflag_src, iflag_order, &
					isok)

			if(iflag_src  .eq. 1)then
				if(it + nwt  .le. lt)then
					v2=vv(ns_z,ns_x)*vv(ns_z, ns_x)
					pus01%u(ns_z, ns_x)=pus01%u(ns_z, ns_x) - wavelet(it+nwt)*dt*dt*v2
				endif

				CALL Extrapolation_2D_Without_PML_One_Step(vv, pus01%u, pus02%u, pus03%u, &
						coe_2nd_dx2, coe_2nd_dz2, &
						npad, nnx, nnz, nx, nz, &
						nx_bound_l, nx_bound_r, &
						nz_bound_u, nz_bound_d, &
						dx, dz, dt, it)

				do ix=1, nx_with_apert
					inx=ix+nx_bound_l
					do iz=1, nz_with_apert
						inz=iz+nz_bound_u
						v2=vv(inz, inx)*vv(inz, inx)
						tmp=pus01%u(inz, inx) + pus03%u(inz, inx) -2.0*pus02%u(inz, inx)
						tmp=tmp/(v2*dt*dt)
						tmp=-1.0*image_r_loc(inz, inx)*tmp
						pusr1%u(inz, inx)= pusr1%u(inz, inx) - tmp
					enddo
				enddo

				CALL Extrapolation_2D_Without_PML_One_Step(vv, pusr1%u, pusr2%u, pusr3%u, &
						coe_2nd_dx2, coe_2nd_dz2, &
						npad, nnx, nnz, nx, nz, &
						nx_bound_l, nx_bound_r, &
						nz_bound_u, nz_bound_d, &
						dx, dz, dt, it)
			endif

			if(iflag_src .eq. 0)then
				if(it .ge. 1 .and. it .le. lt-2) then
					if(mod(it, it_step) .eq. 0) then
						iit=it/it_step
						pus01%u = ust0%u(:, :, iit+2)
						pus02%u = ust0%u(:, :, iit+1)
						pus03%u = ust0%u(:, :, iit)
					endif
				endif
			endif

			if(iflag_src .eq. 0)then
				if(it .ge. 1 .and. it .le. lt-2) then
					if(mod(it, it_step) .eq. 0) then
						iit=it/it_step
						pusr1%u = ustr%u(:, :, iit+2)
						pusr2%u = ustr%u(:, :, iit+1)
						pusr3%u = ustr%u(:, :, iit)
					endif
				endif
			endif

			CALL  Extrapolation_2D_Add_PML_One_Step(vv, pur01%u, pur02%u, pur03%u, &
						pux0_a1%u, pux0_a2%u, pux0_a3%u, &
						pux0_b1%u, pux0_b2%u, pux0_bp1%u, pux0_bp2%u, pux0_bp3%u, &
						pux0_c1%u, pux0_c2%u, pux0_c3%u, &
						puz0_a1%u, puz0_a2%u, puz0_a3%u, &
						puz0_b1%u, puz0_b2%u, puz0_bp1%u, puz0_bp2%u, puz0_bp3%u, &
						puz0_c1%u, puz0_c2%u, puz0_c3%u, &
						funa, dfuna, &
						line_xl, line_xr, line_zu, line_zd, &
						coe_2nd_dx2, coe_2nd_dz2, &
						coe_1st_dx, coe_1st_dz, &
						npad, nnx, nnz, nx, nz, &
						nx_bound_l, nx_bound_r, &
						nz_bound_u, nz_bound_d, &
						dx, dz, dt, r, itt)

			do itr=1, ntr
				nr_x=index_recgx(itr)
				nr_z=index_recgz(itr)
				if( nr_x .ge. 1   .and. nr_z .ge. 1   .and. &
					nr_x .le. nnx .and. nr_z .le. nnz )then
					v2=vv(nr_z, nr_x)*vv(nr_z, nr_x)
					pur03%u(nr_z, nr_x) = pur03%u(nr_z, nr_x) + shotres(it, itr)*v2*dt*dt
				endif
			enddo

			if(it .gt. 1 .and. it .lt. 1.0*lt)then
				CALL Cor_Imaging_Condition_2d(image_i_sig, vv, &
								pus01%u, pus02%u, pus03%u, &
								pur01%u, pur02%u, pur03%u, &
								nnz, nnx, npad, nz_with_apert, nx_with_apert, &
								nz_bound_u, nx_bound_l, dt)

			endif

			CALL born_extrapolation_2d_backward_variable_update(iflag_src)

		enddo		


		CALL rvcrti2d_image_source_position_abnormal_remove(image_i_sig, &
							nz_with_apert, nx_with_apert, &
							nz_bound_u, nx_bound_l, ns_z, ns_x)


		CALL born_2d_single_shot_extrapolation_deallocate(iflag_src)

		isok=0

	    return
    end	subroutine


!***********************************************************************
	subroutine rvcrti2d_image_source_position_abnormal_remove(grad, &
							nz_with_apert, nx_with_apert, &
							nz_bound_u, nx_bound_l, ns_z, ns_x)

		use global
		implicit none
		!Data dictionary: declare variable types, definitons, & units
		!Dummy Variables
		integer	::	nz_with_apert, nx_with_apert
		integer	::	nz_bound_u, nx_bound_l
		integer	::	ns_z, ns_x
		real	::	grad(nz_with_apert, nx_with_apert)

		!Local Variables
		integer	::	nwin 
		integer	::	ix, ixx
		integer ::	iz, izz 
		integer	::	ii
		real	::	tmp

		nwin = 2
		ix = ns_x - nx_bound_l
		iz = ns_z - nz_bound_u
		tmp = 0.0
		do ii=1, nwin
			tmp = tmp + grad(iz, ix+ii) + grad(iz, ix-ii)
		enddo
		grad(iz, ix) = tmp/(nwin*2)

		return
    end	subroutine


!***********************************************************************

