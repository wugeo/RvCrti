
	subroutine Cor_Imaging_Condition_2d(image, vv, us1, us2, us3, ur1, ur2, ur3,&
							nnz, nnx, npad, nz_with_apert, nx_with_apert, &
							nz_bound_u, nx_bound_l ,dt)

		implicit none
		integer	::	nnz, nnx
		integer	::	npad
		integer	::	nz_with_apert, nx_with_apert
		integer	::	nz_bound_u, nx_bound_l
		real	::	dt


		real	::	image(nz_with_apert, nx_with_apert)
		real	::	vv(nnz, nnx)
		real	::	us1(-npad+1:nnz+npad, -npad+1:nnx+npad)
		real	::	us2(-npad+1:nnz+npad, -npad+1:nnx+npad)
		real	::	us3(-npad+1:nnz+npad, -npad+1:nnx+npad)
		real	::	ur1(-npad+1:nnz+npad, -npad+1:nnx+npad)
		real	::	ur2(-npad+1:nnz+npad, -npad+1:nnx+npad)
		real	::	ur3(-npad+1:nnz+npad, -npad+1:nnx+npad)


		!Local Variables
		integer	::	ix, ixx
		integer ::	iz, izz 
		real	::	tmp, tmp1


		do ix=1, nx_with_apert
			ixx = ix + nx_bound_l
			do iz=1, nz_with_apert
				izz=iz + nz_bound_u
				image(iz, ix)=image(iz, ix) + us3(izz, ixx)*ur3(izz, ixx)
	   	     enddo
		enddo

		return
    end	subroutine

