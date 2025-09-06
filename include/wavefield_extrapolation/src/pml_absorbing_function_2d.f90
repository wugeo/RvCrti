




!***********************************************************************
	subroutine absorbing_function_2d(funa, dfuna, vv, nnz, nnx, dx, dz, &
				line_xl, line_xr, line_zu, line_zd, &
				nx_bound_l, nx_bound_r, &
				nz_bound_u, nz_bound_d, &
				r)

    implicit none
	!Data dictionary: declare variable types, definitons, & units
    !dummy variables
    integer	::	nnx, nnz
	integer ::	nx_bound_l, nx_bound_r
    integer ::	nz_bound_u, nz_bound_d
    real	::	dx, dz
    real	::	vv(nnz, nnx)
    real	::	funa(nnz, nnx)
    real	::	dfuna(nnz, nnx)
	integer	::	line_xl(nx_bound_l, 2), line_xr(nx_bound_r, 2)
	integer	::	line_zu(nz_bound_u, 2), line_zd(nz_bound_d, 2)
	real	::	r


    !local variables
    integer	::	ix, iz
    integer	::	iix, iiz
    integer	::	inx, inz
	real	::	tmp
	real	::	dis, dis3
	real	::	x, z
	real	::	kzuxl, kzdxl, kzuxr, kzdxr
	real	::	kxlzu, kxrzu, kxlzd, kxrzd


	tmp=3.0/2.0*log(1.0/r)


	kzuxl=(1.0*(nz_bound_u-1))/(1.0*(nx_bound_l-1))
	kzdxl=(1.0*(nz_bound_d-1))/(1.0*(nx_bound_l-1))
	kzuxr=(1.0*(nz_bound_u-1))/(1.0*(nx_bound_r-1))
	kzdxr=(1.0*(nz_bound_d-1))/(1.0*(nx_bound_r-1))
	kxlzu=(1.0*(nx_bound_l-1))/(1.0*(nz_bound_u-1))
	kxrzu=(1.0*(nx_bound_r-1))/(1.0*(nz_bound_u-1))
	kxlzd=(1.0*(nx_bound_l-1))/(1.0*(nz_bound_d-1))
	kxrzd=(1.0*(nx_bound_r-1))/(1.0*(nz_bound_d-1))


	do ix=1, nx_bound_l
		line_xl(ix,1) = int((ix-1)*kzuxl)
		line_xl(ix,2) = int(nnz-(ix-1)*kzdxl+1)
		if(line_xl(ix, 1) .lt. 1)	line_xl(ix, 1) = 1
		if(line_xl(ix, 2) .gt. nnz)	line_xl(ix, 2) = nnz
	enddo

	do ix=1, nx_bound_r
		line_xr(ix, 1) = int(nz_bound_u-(ix-1)*kzuxr-1) 
		line_xr(ix, 2) = int(nnz-nz_bound_d+(ix-1)*kzdxr+1)
		if(line_xr(ix, 1) .lt. 1)	line_xr(ix, 1) = 1
		if(line_xr(ix, 2) .gt. nnz)	line_xr(ix, 2) = nnz
	enddo

	do iz=1, nz_bound_u
		line_zu(iz, 1) = int((iz-1)*kxlzu)
		line_zu(iz, 2) = int(nnx-(iz-1)*kxrzu+1)
		if(line_zu(iz, 1) .lt. 1)	line_zu(iz, 1) = 1
		if(line_zu(iz, 2) .gt. nnx)	line_zu(iz, 2) = nnx
	enddo

	do iz=1, nz_bound_d
		line_zd(iz, 1) = int(nx_bound_l-(iz-1)*kxlzd-1)
		line_zd(iz, 2) = int(nnx-nx_bound_r +(iz-1)*kxrzd+1)
		if(line_zd(iz, 1) .lt. 1)	line_zd(iz, 1) = 1
		if(line_zd(iz, 2) .gt. nnx)	line_zd(iz ,2) = nnx
	enddo
		



	dis=nx_bound_l*dx
	dis3=dis*dis*dis
	do ix=1, nx_bound_l
		inx=ix
		do inz=line_xl(ix, 1), line_xl(ix, 2)
			x=(nx_bound_l-ix+1)*dx
			funa(inz, inx)=tmp*vv(inz, inx)*x*x/dis3
			dfuna(inz, inx)=tmp*2.0*vv(inz, inx)*x/dis3
		enddo	
	enddo

	dis=nx_bound_r*dx
	dis3=dis*dis*dis
	do ix=1, nx_bound_r
		inx= nnx-nx_bound_r+ix
		do inz=line_xr(ix, 1), line_xr(ix, 2) 
			x=ix*dx
			funa(inz, inx)=tmp*vv(inz, inx)*x*x/dis3
			dfuna(inz, inx)=tmp*2.0*vv(inz, inx)*x/dis3
		enddo	
	enddo

	dis=nz_bound_u*dz
	dis3=dis*dis*dis
	do iz=1, nz_bound_u
		inz=iz
		do inx=line_zu(iz, 1), line_zu(iz, 2)
			z=(nz_bound_u-iz+1)*dz
			funa(inz, inx)=tmp*vv(inz, inx)*z*z/dis3
			dfuna(inz, inx)=tmp*2.0*vv(inz, inx)*z/dis3
		enddo	
	enddo

	dis=nz_bound_d*dz
	dis3=dis*dis*dis
	do iz=1, nz_bound_d
		inz= nnz-nz_bound_d+iz
		do inx=line_zd(iz, 1), line_zd(iz, 2)
			z=iz*dz
			funa(inz, inx)=tmp*vv(inz, inx)*z*z/dis3
			dfuna(inz, inx)=tmp*2.0*vv(inz, inx)*z/dis3
		enddo
	enddo

 	do inx=nx_bound_l+1, nnx-nx_bound_r
        do inz=nz_bound_u+1, nnz-nz_bound_d
			funa(inz, inx)=0.0
			dfuna(inz, inx)=0.0
		enddo
    enddo



    end subroutine

!***********************************************************************
