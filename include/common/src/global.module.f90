!*******************************************************************!
!*        Begin The Define of module Need In This Program           ! 
!*==================================================================!
      module global                                                !
       real   ,parameter :: PI = 3.14159265359                     !
	   real   ,parameter :: PI2 = 6.28318530718
       integer,parameter :: lbyte = 1                              !
	   integer,parameter :: para_int_fnum	=	100
	   integer,parameter :: para_long_fnum	=	50
	   integer,parameter :: para_float_fnum	=	50
	   integer,parameter :: para_double_fnum=	50
	   integer,parameter :: para_char_fnum	=	50
	   integer,parameter :: para_char_flen	=	256

       integer,parameter :: fdebug_print_flag = 1
    !    integer,parameter :: fdebug_print_flag = 0

       integer,parameter :: fdebug_info_flag = 1
      !  integer,parameter :: fdebug_info_flag = 0

       integer,parameter :: fdebug_io_flag = 1
      !  integer,parameter :: fdebug_io_flag = 0

       integer,parameter :: fdebug_error_flag = 1
      !  integer,parameter :: fdebug_error_flag = 0

      character(len=para_char_flen)  :: 	gl_method

      end module global                                            !

!*==================================================================!
!*        END Of The Define of module Need In This Program          ! 
!*******************************************************************!

