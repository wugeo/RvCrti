!***********************************************************************
	subroutine source_wavefield_rtm2d(u1, u2, u3, &
				top, bot, lef, rig, &
				nnx, nnz, &
				nx_with_apert, nz_with_apert, &
				nx_bound_l, nx_bound_r, &
				nz_bound_u, nz_bound_d, &
				npad, order, lt, ns, &
				it, it_step, &
				iflag_src, iflag, &
				isok)

		implicit none
		!Data dictionary: declare variable types, definitons, & units
		!Dummy Variables
		integer	::	nnx, nnz
		integer	::	nx_with_apert, nz_with_apert
		integer ::	nx_bound_l, nx_bound_r
		integer ::	nz_bound_u, nz_bound_d
		integer	::	npad, order, lt, ns
		integer	::	it, it_step
		integer	::	iflag_src, iflag
		integer	::	isok

		real	::	u1(-npad+1:nnz+npad, -npad+1:nnx+npad)
		real	::	u2(-npad+1:nnz+npad, -npad+1:nnx+npad)
		real	::	u3(-npad+1:nnz+npad, -npad+1:nnx+npad)
		real	::	top(order, nx_with_apert, ns)
		real	::	bot(order, nx_with_apert, ns)
		real	::	lef(order, nz_with_apert, ns)
		real	::	rig(order, nz_with_apert, ns)

		!Local Variables
		integer	::	ix, iy, iz

		if(iflag_src  .eq. 1)then
			if(iflag .eq. 1)then
				CALL wavefield_reconstruction_forward_rtm2d(u3, &
						top, bot, lef, rig, &
						nnx, nnz, &
						nx_with_apert, nz_with_apert, &
						nx_bound_l, nx_bound_r, &
						nz_bound_u, nz_bound_d, &
						npad, order, lt, ns, &
						it, it_step, &
						isok)
			endif
			if(iflag .eq. -1)then
				CALL wavefield_reconstruction_backward_rtm2d(u1, u2, u3, &
						top, bot, lef, rig, &
						nnx, nnz, &
						nx_with_apert, nz_with_apert, &
						nx_bound_l, nx_bound_r, &
						nz_bound_u, nz_bound_d, &
						npad, order, lt, ns, &
						it, it_step, &
						isok)
			endif
		endif

	end subroutine

!***********************************************************************
	subroutine wavefield_reconstruction_forward_rtm2d(u, &
				top, bot, lef, rig, &
				nnx, nnz, &
				nx_with_apert, nz_with_apert, &
				nx_bound_l, nx_bound_r, &
				nz_bound_u, nz_bound_d, &
				npad, order, lt, ns, &
				it, it_step, &
				isok)

		implicit none
		!Data dictionary: declare variable types, definitons, & units
		!Dummy Variables
		integer	::	nnx, nnz
		integer	::	nx_with_apert, nz_with_apert
		integer ::	nx_bound_l, nx_bound_r
		integer ::	nz_bound_u, nz_bound_d
		integer	::	npad, order, lt, ns

		real	::	u(-npad+1:nnz+npad, -npad+1:nnx+npad)
		real	::	top(order, nx_with_apert, ns)
		real	::	bot(order, nx_with_apert, ns)
		real	::	lef(order, nz_with_apert, ns)
		real	::	rig(order, nz_with_apert, ns)

		integer	::	it, it_step
		integer	::	isok

		!Local Variables
		integer	::	ix, iz
		integer	::	inx, inz
		integer	::	iii, iit


		if(mod(it, it_step) .eq. 0 .and. it .ge. 1 .and. it .le. lt) then
			iit=it/it_step
			do ix=1, nx_with_apert
				inx=ix+nx_bound_l
				do iii=1, order
					top(iii, ix, iit) = u(nz_bound_u+iii-order/2, inx)
					bot(iii, ix, iit) = u(nnz-nz_bound_d-iii+1+order/2, inx)
				enddo
			enddo
			do iz=1, nz_with_apert
				inz=iz+nz_bound_u
				do iii=1, order
					lef(iii, iz, iit) = u(inz, nx_bound_l+iii-order/2)
					rig(iii, iz, iit) = u(inz, nnx-nx_bound_r-iii+1+order/2)
				enddo
			enddo
		endif

		isok=0

	end subroutine

!***********************************************************************
	subroutine wavefield_reconstruction_backward_rtm2d(u1, u2, u3, &
				top, bot, lef, rig, &
				nnx, nnz, &
				nx_with_apert, nz_with_apert, &
				nx_bound_l, nx_bound_r, &
				nz_bound_u, nz_bound_d, &
				npad, order, lt, ns, &
				it, it_step, &
				isok)

		implicit none

		!Data dictionary: declare variable types, definitons, & units
		!Dummy Variables
		integer	::	nnx, nnz
		integer	::	nx_with_apert, nz_with_apert
		integer ::	nx_bound_l, nx_bound_r
		integer ::	nz_bound_u, nz_bound_d
		integer	::	npad, order, lt, ns

		real	::	u1(-npad+1:nnz+npad, -npad+1:nnx+npad)
		real	::	u2(-npad+1:nnz+npad, -npad+1:nnx+npad)
		real	::	u3(-npad+1:nnz+npad, -npad+1:nnx+npad)
		real	::	top(order, nx_with_apert, ns)
		real	::	bot(order, nx_with_apert, ns)
		real	::	lef(order, nz_with_apert, ns)
		real	::	rig(order, nz_with_apert, ns)

		integer	::	it, it_step
		integer	::	isok

		!Local Variables
		integer	::	ix, iz
		integer	::	inx, inz
		integer	::	iii, iit


		if(mod(it, it_step) .eq. 0 .and. it .ge. 1 .and. it .le. lt-1) then
			iit=it/it_step

			do ix=1, nx_with_apert
				inx=ix+nx_bound_l
				do iii=1, order
					u2(nz_bound_u+iii-order/2, inx) = top(iii, ix, iit+1)
					u2(nnz-nz_bound_d-iii+1+order/2, inx) = bot(iii, ix, iit+1)

					u2(nz_bound_u+iii-order/2, inx) = top(iii, ix, iit+1)
					u2(nnz-nz_bound_d-iii+1+order/2, inx) = bot(iii, ix, iit+1)
				enddo
			enddo
			do iz=1, nz_with_apert
				inz=iz+nz_bound_u
				do iii=1, order
					u2(inz, nx_bound_l+iii-order/2) = lef(iii, iz, iit+1)
					u2(inz, nnx-nx_bound_r-iii+1+order/2) = rig(iii, iz, iit+1)

					u2(inz, nx_bound_l+iii-order/2) = lef(iii, iz, iit+1)
					u2(inz, nnx-nx_bound_r-iii+1+order/2) = rig(iii, iz, iit+1)
				enddo
			enddo
		endif

		isok=0

	end subroutine

!***********************************************************************