!***********************************************************************
	subroutine RV_CRTI_2D_Imaging_Condition_m1(image_ms, image_dm, vv, &
							us01, us02, us03, ur01, ur02, ur03, &
							usr1, usr2, usr3, urr1, urr2, urr3, &
							nnz, nnx, npad, nz_with_apert, nx_with_apert, &
							nz_bound_u, nx_bound_l, dt)
		use global
		implicit none
		!Data dictionary: declare variable types, definitons, & units
		!Dummy Variables
		integer	::	nnz, nnx
		integer	::	npad
		integer	::	nz_with_apert, nx_with_apert
		integer	::	nz_bound_u, nx_bound_l
		real	::	dt

		real	::	image_ms(nz_with_apert, nx_with_apert)
		real	::	image_dm(nz_with_apert, nx_with_apert)
		real	::	vv(nnz, nnx)
		real	::	us01(-npad+1:nnz+npad, -npad+1:nnx+npad)
		real	::	us02(-npad+1:nnz+npad, -npad+1:nnx+npad)
		real	::	us03(-npad+1:nnz+npad, -npad+1:nnx+npad)
		real	::	ur01(-npad+1:nnz+npad, -npad+1:nnx+npad)
		real	::	ur02(-npad+1:nnz+npad, -npad+1:nnx+npad)
		real	::	ur03(-npad+1:nnz+npad, -npad+1:nnx+npad)

		real	::	usr1(-npad+1:nnz+npad, -npad+1:nnx+npad)
		real	::	usr2(-npad+1:nnz+npad, -npad+1:nnx+npad)
		real	::	usr3(-npad+1:nnz+npad, -npad+1:nnx+npad)
		real	::	urr1(-npad+1:nnz+npad, -npad+1:nnx+npad)
		real	::	urr2(-npad+1:nnz+npad, -npad+1:nnx+npad)
		real	::	urr3(-npad+1:nnz+npad, -npad+1:nnx+npad)

		!Local Variables
		integer	::	ix, ixx
		integer ::	iz, izz 
		real	::	tmp, tmp1


			do ix=1, nx_with_apert
				ixx = ix + nx_bound_l
				do iz=1, nz_with_apert
					izz=iz + nz_bound_u
			
					tmp1=(us01(izz, ixx)+us03(izz, ixx)-2*us02(izz, ixx))
					tmp1=tmp1/(dt*dt)
					tmp=-1.0*tmp1*urr3(izz, ixx)
					image_ms(iz, ix)=image_ms(iz, ix) + tmp

					tmp1=(usr1(izz, ixx)+usr3(izz, ixx)-2*usr2(izz, ixx))
					tmp1=tmp1/(dt*dt)
					tmp=-1.0*tmp1*ur03(izz, ixx)
					image_dm(iz, ix)=image_dm(iz, ix) + tmp

		   	     enddo
			enddo

	    return
    end	subroutine

	subroutine RV_CRTI_2D_Imaging_Condition_m2(image_rv_ms, image_rv_dm, vv, &
							us01, us02, us03, &
							usr1, usr2, usr3, &
							uv_r01, uv_r02, uv_r03, &
							uv_rr1, uv_rr2, uv_rr3, &
							nnz, nnx, npad, nz_with_apert, nx_with_apert, &
							nz_bound_u, nx_bound_l, dt)
		use global
		implicit none
		!Data dictionary: declare variable types, definitons, & units
		!Dummy Variables
		integer	::	nnz, nnx
		integer	::	npad
		integer	::	nz_with_apert, nx_with_apert
		integer	::	nz_bound_u, nx_bound_l
		real	::	dt

		real	::	image_rv_ms(nz_with_apert, nx_with_apert)
		real	::	image_rv_dm(nz_with_apert, nx_with_apert)
		real	::	vv(nnz, nnx)

		real	::	us01(-npad+1:nnz+npad, -npad+1:nnx+npad)
		real	::	us02(-npad+1:nnz+npad, -npad+1:nnx+npad)
		real	::	us03(-npad+1:nnz+npad, -npad+1:nnx+npad)
		real	::	usr1(-npad+1:nnz+npad, -npad+1:nnx+npad)
		real	::	usr2(-npad+1:nnz+npad, -npad+1:nnx+npad)
		real	::	usr3(-npad+1:nnz+npad, -npad+1:nnx+npad)

		real	::	uv_r01(-npad+1:nnz+npad, -npad+1:nnx+npad)
		real	::	uv_r02(-npad+1:nnz+npad, -npad+1:nnx+npad)
		real	::	uv_r03(-npad+1:nnz+npad, -npad+1:nnx+npad)
		real	::	uv_rr1(-npad+1:nnz+npad, -npad+1:nnx+npad)
		real	::	uv_rr2(-npad+1:nnz+npad, -npad+1:nnx+npad)
		real	::	uv_rr3(-npad+1:nnz+npad, -npad+1:nnx+npad)

		!Local Variables
		integer	::	ix, ixx
		integer ::	iz, izz 
		real	::	tmp, tmp1


		do ix=1, nx_with_apert
			ixx = ix + nx_bound_l
			do iz=1, nz_with_apert
				izz=iz + nz_bound_u
			
				tmp1=(us01(izz, ixx)+us03(izz, ixx)-2*us02(izz, ixx))
				tmp1=tmp1/(dt*dt)
				tmp=-1.0*tmp1*uv_rr3(izz, ixx)
				image_rv_ms(iz, ix)=image_rv_ms(iz, ix) + tmp

				tmp1=(usr1(izz, ixx)+usr3(izz, ixx)-2*usr2(izz, ixx))
				tmp1=tmp1/(dt*dt)
				tmp=-1.0*tmp1*uv_r03(izz, ixx)
				image_rv_dm(iz, ix)=image_rv_dm(iz, ix) + tmp
		   	enddo
		enddo

	    return
    end	subroutine

!***********************************************************************

