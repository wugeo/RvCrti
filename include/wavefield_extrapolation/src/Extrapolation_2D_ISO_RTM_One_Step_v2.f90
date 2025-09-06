!*2020.04.28
!***********************************************************************
    subroutine  Extrapolation_2D_Add_PML_One_Step(vv, u1, u2, u3, &
			ux_a1, ux_a2, ux_a3, &
			ux_b1, ux_b2, ux_bp1, ux_bp2, ux_bp3, &
			ux_c1, ux_c2, ux_c3, &
			uz_a1, uz_a2, uz_a3, &
			uz_b1, uz_b2, uz_bp1, uz_bp2, uz_bp3, &
			uz_c1, uz_c2, uz_c3, &
			funa, dfuna, &
			line_xl, line_xr, line_zu, line_zd, &
			coe_2nd_dx2, coe_2nd_dz2, &
			coe_1st_dx, coe_1st_dz, &
			npad, nnx, nnz, nx, nz, &
			nx_bound_l, nx_bound_r, &
			nz_bound_u, nz_bound_d, &
			dx, dz, dt, r, it)

    implicit none
    integer ::	nnx, nnz
    integer ::	nx, nz
    integer ::	nx_bound_l, nx_bound_r
	integer ::	nz_bound_u, nz_bound_d
	integer	::	npad
    real    ::	dx, dz, dt
    real    ::	r
	integer ::	it
  
    real	::	vv(nnz, nnx)
    real	::	u1(-npad+1:nnz+npad, -npad+1:nnx+npad)
    real	::	u2(-npad+1:nnz+npad, -npad+1:nnx+npad)
    real	::	u3(-npad+1:nnz+npad, -npad+1:nnx+npad)

    real	::	ux_a1(nnz,nx_bound_l+nx_bound_r)
    real	::	ux_a2(nnz,nx_bound_l+nx_bound_r)
    real	::	ux_a3(nnz,nx_bound_l+nx_bound_r)
    real	::	ux_b1(nnz,nx_bound_l+nx_bound_r)
    real	::	ux_b2(nnz,nx_bound_l+nx_bound_r)
    real	::	ux_bp1(nnz,nx_bound_l+nx_bound_r)
    real	::	ux_bp2(nnz,nx_bound_l+nx_bound_r)
    real	::	ux_bp3(nnz,nx_bound_l+nx_bound_r)
    real	::	ux_c1(nnz,nx_bound_l+nx_bound_r)
    real	::	ux_c2(nnz,nx_bound_l+nx_bound_r)
    real	::	ux_c3(nnz,nx_bound_l+nx_bound_r)

    real	::	uz_a1(nz_bound_u+nz_bound_d,nnx)
    real	::	uz_a2(nz_bound_u+nz_bound_d,nnx)
    real	::	uz_a3(nz_bound_u+nz_bound_d,nnx)
    real	::	uz_b1(nz_bound_u+nz_bound_d,nnx)
    real	::	uz_b2(nz_bound_u+nz_bound_d,nnx)
    real	::	uz_bp1(nz_bound_u+nz_bound_d,nnx)
    real	::	uz_bp2(nz_bound_u+nz_bound_d,nnx)
    real	::	uz_bp3(nz_bound_u+nz_bound_d,nnx)
    real	::	uz_c1(nz_bound_u+nz_bound_d,nnx)
    real	::	uz_c2(nz_bound_u+nz_bound_d,nnx)
    real	::	uz_c3(nz_bound_u+nz_bound_d,nnx)

    real	::	funa(nnz, nnx)
    real	::	dfuna(nnz, nnx)
	real	::	coe_2nd_dx2(npad), coe_2nd_dz2(npad)
	real	::	coe_1st_dx(npad), coe_1st_dz(npad)
	integer	::	line_xl(nx_bound_l, 2), line_xr(nx_bound_r, 2)
	integer	::	line_zu(nz_bound_u, 2), line_zd(nz_bound_d, 2)


	integer	::	ii
    integer ::	ix, iz
    integer ::	inx, inz
    integer ::	iix, iiz
    real    ::	v2, dt2
	real	::	x1_deri, z1_deri
	real	::	x2_deri, z2_deri
	real	::	a, da

      
	dt2=dt*dt


	do ix = nx_bound_l+1, nnx-nx_bound_r
		inx = ix
		do iz = nz_bound_u+1, nnz-nz_bound_d
			inz = iz
			v2=vv(inz, inx)*vv(inz, inx)
			z2_deri = 0.0
			x2_deri = 0.0
			do ii = 1, npad
				z2_deri = z2_deri + coe_2nd_dz2(ii)*(u2(inz+ii, inx) + u2(inz-ii, inx) - 2.0*u2(inz, inx))
				x2_deri = x2_deri + coe_2nd_dx2(ii)*(u2(inz, inx+ii) + u2(inz, inx-ii) - 2.0*u2(inz, inx))
			enddo
			u3(inz, inx) = 2.0*u2(inz, inx) - u1(inz, inx) + dt2*v2*(z2_deri+x2_deri)
		enddo
	enddo

	do ix = 1, nx_bound_l
		iix = ix
		inx = ix
		do iz = line_xl(ix, 1), line_xl(ix, 2)
			iiz = iz
			inz = iz
			v2=vv(inz, inx)*vv(inz, inx)

			z2_deri = 0.0
			x2_deri = 0.0
			x1_deri = 0.0
			do ii = 1, npad
				z2_deri = z2_deri + coe_2nd_dz2(ii)*(u2(inz+ii, inx) + u2(inz-ii, inx) - 2.0*u2(inz, inx))
				x2_deri = x2_deri + coe_2nd_dx2(ii)*(u2(inz, inx+ii) + u2(inz, inx-ii) - 2.0*u2(inz, inx))
				x1_deri = x1_deri + coe_1st_dx(ii) *(u2(inz, inx+ii) - u2(inz, inx-ii)) 
			enddo

			a = funa(inz, inx)
			da = dfuna(inz, inx)

			ux_a3(iiz, iix) = 2.0*ux_a2(iiz, iix) - ux_a1(iiz, iix) + &
				dt2*(v2*x2_deri - 2.0*a*(ux_a2(iiz, iix) - ux_a1(iiz, iix))/dt - a*a*ux_a2(iiz, iix)) 
			ux_bp3(iiz, iix) = 2.0*ux_bp2(iiz, iix) - ux_bp1(iiz, iix) + &
				dt2*(-1.0*v2*da*x1_deri - 2.0*a*(ux_bp2(iiz, iix) - ux_bp1(iiz, iix))/dt - a*a*ux_bp2(iiz, iix)) 
			ux_b2(iiz, iix) = ux_b1(iiz, iix) + dt*(ux_bp2(iiz, iix) - a*ux_b1(iiz, iix))
			ux_c3(iiz, iix) = 2.0*ux_c2(iiz, iix) - ux_c1(iiz, iix) + dt2*v2*z2_deri
			u3(inz, inx) = ux_a3(iiz, iix) + ux_b2(iiz, iix) + ux_c3(iiz, iix)
		enddo
	enddo

	do ix = 1, nx_bound_r
		iix = nx_bound_l+ix
		inx = nnx-nx_bound_r+ix

		do iz=line_xr(ix, 1), line_xr(ix, 2)
			iiz = iz
			inz = iz
			v2=vv(inz, inx)*vv(inz, inx)
			z2_deri = 0.0
			x2_deri = 0.0
			x1_deri = 0.0
			do ii = 1, npad
				z2_deri = z2_deri + coe_2nd_dz2(ii)*(u2(inz+ii, inx) + u2(inz-ii, inx) - 2.0*u2(inz, inx))
				x2_deri = x2_deri + coe_2nd_dx2(ii)*(u2(inz, inx+ii) + u2(inz, inx-ii) - 2.0*u2(inz, inx))
				x1_deri = x1_deri + coe_1st_dx(ii) *(u2(inz, inx+ii) - u2(inz, inx-ii)) 
			enddo

			a = funa(inz, inx)
			da = dfuna(inz, inx)

			ux_a3(iiz, iix) = 2.0*ux_a2(iiz, iix) - ux_a1(iiz, iix) + &
				dt2*(v2*x2_deri - 2.0*a*(ux_a2(iiz, iix) - ux_a1(iiz, iix))/dt - a*a*ux_a2(iiz, iix)) 
			ux_bp3(iiz, iix) = 2.0*ux_bp2(iiz, iix) - ux_bp1(iiz, iix) + &
				dt2*(-1.0*v2*da*x1_deri - 2.0*a*(ux_bp2(iiz, iix) - ux_bp1(iiz, iix))/dt - a*a*ux_bp2(iiz, iix)) 
			ux_b2(iiz, iix) = ux_b1(iiz, iix) + dt*(ux_bp2(iiz, iix) - a*ux_b1(iiz, iix))
			ux_c3(iiz, iix) = 2.0*ux_c2(iiz, iix) - ux_c1(iiz, iix) + dt2*v2*z2_deri
			u3(inz, inx) = ux_a3(iiz, iix) + ux_b2(iiz, iix) + ux_c3(iiz, iix)
		enddo
	enddo

	do iz = 1, nz_bound_u
		iiz = iz
		inz = iz

		do ix = line_zu(iz, 1), line_zu(iz, 2)
			iix = ix
			inx = ix
			v2=vv(inz, inx)*vv(inz, inx)
			z2_deri = 0.0
			x2_deri = 0.0
			z1_deri = 0.0
			do ii = 1, npad
				z2_deri = z2_deri + coe_2nd_dz2(ii)*(u2(inz+ii, inx) + u2(inz-ii, inx) - 2.0*u2(inz, inx))
				x2_deri = x2_deri + coe_2nd_dx2(ii)*(u2(inz, inx+ii) + u2(inz, inx-ii) - 2.0*u2(inz, inx))
				z1_deri = z1_deri + coe_1st_dz(ii) *(u2(inz+ii, inx) - u2(inz-ii, inx)) 
			enddo

			a = funa(inz, inx)
			da = dfuna(inz, inx)

			uz_a3(iiz, iix) = 2.0*uz_a2(iiz, iix) - uz_a1(iiz, iix) + &
				dt2*(v2*z2_deri - 2.0*a*(uz_a2(iiz, iix) - uz_a1(iiz, iix))/dt - a*a*uz_a2(iiz, iix)) 
			uz_bp3(iiz, iix) = 2.0*uz_bp2(iiz, iix) - uz_bp1(iiz, iix) + &
				dt2*(-1.0*v2*da*z1_deri - 2.0*a*(uz_bp2(iiz, iix) - uz_bp1(iiz, iix))/dt - a*a*uz_bp2(iiz, iix)) 
			uz_b2(iiz, iix) = uz_b1(iiz, iix) + dt*(uz_bp2(iiz, iix) - a*uz_b1(iiz, iix))
			uz_c3(iiz, iix) = 2.0*uz_c2(iiz, iix) - uz_c1(iiz, iix) + dt2*v2*x2_deri
			u3(inz, inx) = uz_a3(iiz, iix) + uz_b2(iiz, iix) + uz_c3(iiz, iix)
		enddo
	enddo

	do iz = 1, nz_bound_d
		iiz = nz_bound_u+iz
		inz = nnz-nz_bound_d+iz

		do ix = line_zd(iz, 1), line_zd(iz, 2)
			iix = ix
			inx = ix
			v2=vv(inz, inx)*vv(inz, inx)
			z2_deri = 0.0
			x2_deri = 0.0
			z1_deri = 0.0
			do ii = 1, npad
				z2_deri = z2_deri + coe_2nd_dz2(ii)*(u2(inz+ii, inx) + u2(inz-ii, inx) - 2.0*u2(inz, inx))
				x2_deri = x2_deri + coe_2nd_dx2(ii)*(u2(inz, inx+ii) + u2(inz, inx-ii) - 2.0*u2(inz, inx))
				z1_deri = z1_deri + coe_1st_dz(ii) *(u2(inz+ii, inx) - u2(inz-ii, inx)) 
			enddo

			a = funa(inz, inx)
			da = dfuna(inz, inx)

			uz_a3(iiz, iix) = 2.0*uz_a2(iiz, iix) - uz_a1(iiz, iix) + &
				dt2*(v2*z2_deri - 2.0*a*(uz_a2(iiz, iix) - uz_a1(iiz, iix))/dt - a*a*uz_a2(iiz, iix)) 
			uz_bp3(iiz, iix) = 2.0*uz_bp2(iiz, iix) - uz_bp1(iiz, iix) + &
				dt2*(-1.0*v2*da*z1_deri - 2.0*a*(uz_bp2(iiz, iix) - uz_bp1(iiz, iix))/dt - a*a*uz_bp2(iiz, iix)) 
			uz_b2(iiz, iix) = uz_b1(iiz, iix) + dt*(uz_bp2(iiz, iix) - a*uz_b1(iiz, iix))
			uz_c3(iiz, iix) = 2.0*uz_c2(iiz, iix) - uz_c1(iiz, iix) + dt2*v2*x2_deri
			u3(inz, inx) = uz_a3(iiz, iix) + uz_b2(iiz, iix) + uz_c3(iiz, iix)
		enddo
	enddo



    return
    end subroutine

!***********************************************************************







