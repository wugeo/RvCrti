	subroutine read_float2d(fn, d2, n1, n2, isok)
		use global
		implicit none
		!Data dictionary: declare variable types, definitons, & units
		!Dummy Variables
		character(len=para_char_flen)	::	fn
		integer		::	n1, n2
		real		::	d2(n1, n2)
		integer		::	isok

		!Local Variables
		integer		::	ierr
		integer		::	i1, i2

		open(93,file=trim(fn), access='direct', action='read', recl=n2*n1, status='old', iostat=ierr)
		if(ierr /= 0)then
			write(*,*)  'Parameter file cannot open right, please check it !!!'
			write(*,*)  'Shutting down the program'
			write(*,*)  'fn', trim(fn)
			stop
		endif
		read(93, rec=1)((d2(i1, i2),i1=1,n1), i2=1, n2)
		close(93)
		if(fdebug_print_flag .eq. 1)then
			write(*,*)'Success read_float2d: n1', n1, 'n2', n2
		endif

		isok=0
		return
	end	subroutine
	!=============================================================
    subroutine write_float2d(fn, d2, n1, n2, isok)
		use global
		implicit none
		!Data dictionary: declare variable types, definitons, & units
		!Dummy Variables
		character(len=para_char_flen)  	::  fn
		integer	::	n1, n2
		integer	::	isok
		real	::	d2(n1, n2)

		!Local Variables
		integer ::	i1, i2
		integer	::	ierr

		open(121, file=trim(fn), access='direct', status='replace', action='write', recl=lbyte*n1*n2, iostat=ierr)
		if(ierr /= 0)then
			write(*,*)  'Parameter file cannot open right, please check it !!!'
			write(*,*)  'Shutting down the program'
			write(*,*)  'fn', trim(fn)
			stop
		endif
		write(121, rec=1)((d2(i1, i2), i1=1, n1), i2=1, n2)
		close(121)
		if(fdebug_print_flag .eq. 1)then
			write(*,*)'Success write_float2d: n1', n1, 'n2', n2
		endif

		isok=0
		return
    end	subroutine
	!=============================================================
	!=============================================================
	subroutine read_image_2d_FloatBinary(fn_mig, image2d, nvz, nvx, isok)
		use global
		implicit none
		!Data dictionary: declare variable types, definitons, & units
		!Dummy Variables
		character(len=para_char_flen)	::	fn_mig
		integer		::	nvz, nvx
		real		::	image2d(nvz, nvx)
		integer		::	isok

		!Local Variables
		integer		::	ierr
		integer		::	ivz, ivx

		CALL read_float2d(fn_mig, image2d, nvz, nvx, isok)
		if(fdebug_print_flag .eq. 1)then
			write(*,*)'Success Read image2d: nvz', nvz, 'nvx', nvx
		endif

		isok=0
		return
	end	subroutine

	!=============================================================
	!=============================================================
	subroutine write_image_2d_FloatBinary(fn_mig, image2d, nvz, nvx, isok)
		use global
		implicit none
		!Data dictionary: declare variable types, definitons, & units
		!Dummy Variables
		character(len=para_char_flen)  	::  fn_mig
		integer	::	nvz, nvx
		real	::	image2d(nvz, nvx)
		integer	::	isok
		!Local Variables
		integer ::	ivz, ivx
		integer	::	ierr

		CALL write_float2d(fn_mig, image2d, nvz, nvx, isok)
		if(fdebug_print_flag .eq. 1)then
			write(*,*)'Success Write image2d: nvz', nvz, 'nvx', nvx
		endif

		isok=0
		return
    end	subroutine


!***********************************************************************
!***********************************************************************
