!*2020.04.28

	subroutine Extrapolation_2D_Without_PML_One_Step(vv, u1, u2, u3, &
			coe_2nd_dx2, coe_2nd_dz2, &
			npad, nnx, nnz, nx, nz, &
			nx_bound_l, nx_bound_r, &
			nz_bound_u, nz_bound_d, &
			dx, dz, dt, it)

    implicit none
    integer ::	nnx, nnz
    integer ::	nx, nz
    integer ::	nx_bound_l, nx_bound_r
	integer ::	nz_bound_u, nz_bound_d
	integer	::	npad
    real    ::	dx, dz, dt
	integer ::	it
  
    real	::	vv(nnz, nnx)
    real	::	u1(-npad+1:nnz+npad, -npad+1:nnx+npad)
    real	::	u2(-npad+1:nnz+npad, -npad+1:nnx+npad)
    real	::	u3(-npad+1:nnz+npad, -npad+1:nnx+npad)
	real	::	coe_2nd_dx2(npad), coe_2nd_dz2(npad)

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


    return

    end subroutine

!*************************************************************************

