!***********************************************************************
    subroutine local_to_all_image_result_2d(image_vel, image, image_r, &
					nz_with_apert, nx_with_apert, nvz, nvx, &
					nx, nz, nx_apert_l, nz_apert_u, &
					nvxx_shift, nvzz_shift, &
					isok)

		use global
		implicit none
		!Data dictionary: declare variable types, definitons, & units
		!Dummy Variables
		integer	::	nz, nx
		integer	::	nz_with_apert, nx_with_apert
		integer	::	nvz, nvx
		integer	::	nx_apert_l, nz_apert_u
		integer	::	nvxx_shift, nvzz_shift
		integer	::	isok
		real	::	image_vel(nvz, nvx)
		real	::	image_r(nvz, nvx)
		real	::	image(nz_with_apert, nx_with_apert)

		!Local Variables
		integer ::	ix, iz
		integer ::	ivx, ivz
		integer	::	irec
		integer	::	shift_x, shift_z
		integer	::	ix0, ishot, ds, nshot, ix_shift, iix1, iix2
		real,allocatable	:: image_tmp(:,:)
		character(len=para_char_flen)	::	currt
		integer	::	ntr
		integer	::	nwin


		ntr = 301
		nshot = 601
		image_tmp = 0.0
		do ishot=1, nshot
			ix_shift = (ishot-1)
			do ivx=-ntr, ntr
				iix1 = nx_with_apert/2 + ivx
				iix2 = ivx + ix_shift
				if(iix1 .ge. 1 .and. iix1 .le. nx_with_apert)then
					if(iix2 .ge. 1 .and. iix2 .le. nvx)then
						do ivz=1, nvz
							image_vel(ivz, iix2)=image_vel(ivz, iix2)+ image(ivz+nz_apert_u, iix1)
						enddo
					endif
				endif
			enddo
		enddo

		nwin = 20
		do ivx=1, nvx
			do ivz = 1, nvz
				if ( ivz .le. 225-nwin  .or. ivz .ge. 225+nwin ) then
					image_vel(ivz, ivx) = 0.0
				end if
			enddo
		enddo



		isok=0
		return

    end	subroutine


!***********************************************************************
    subroutine local_tovel_result_2d(image_vel, image,  &
					nz_with_apert, nx_with_apert, nvz, nvx, &
					nx, nz, nx_apert_l, nz_apert_u, &
					nvxx_shift, nvzz_shift, &
					isok)

		use global
		implicit none
		!Data dictionary: declare variable types, definitons, & units
		!Dummy Variables
		integer	::	nz, nx
		integer	::	nz_with_apert, nx_with_apert
		integer	::	nvz, nvx
		integer	::	nx_apert_l, nz_apert_u
		integer	::	nvxx_shift, nvzz_shift
		integer	::	isok
		real	::	image_vel(nvz, nvx)
		real	::	image(nz_with_apert, nx_with_apert)

		!Local Variables
		integer ::	ix, iz
		integer ::	ivx, ivz
		integer	::	irec
		integer	::	shift_x, shift_z
		character(len=para_char_flen)	::	currt
		integer	::	ntr
		integer	::	nwin

		
		shift_x = nvxx_shift + nx_apert_l
		shift_z = nvzz_shift + nz_apert_u

		do ix=1, nx
			ivx= ix + shift_x
			if(ivx .ge. 1 .and. ivx .le. nvx) then
				do iz=1, nz
					ivz= iz + shift_z
					if(ivz .ge. 1 .and. ivz .le. nvz)then
						image_vel(ivz, ivx)=image_vel(ivz, ivx) + image(iz+nz_apert_u, ix+nx_apert_l)
					endif
				enddo
			endif
		enddo

		isok=0
		return

    end	subroutine


!***********************************************************************


