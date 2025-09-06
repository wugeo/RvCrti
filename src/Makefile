PackageDir_src :=	..
include $(PackageDir_src)/config/Makefile_inc
Debug_symbols	:=	FALSE
include $(PackageDir_src)/config/Makefile_compile
target_name		:=	rvcrti2d

lib_sub_name	:=	lib$(target_name)sub.a
lib_main_name	:=	lib$(target_name)main.a

target_c		:= 	$(target_name)
target_f		:= 	$(target_name)
CUR_DIR		= .
SRC_DIR		= $(CUR_DIR)/src
OBJ_DIR		= $(CUR_DIR)/obj
LIB_DIR		= $(CUR_DIR)/lib

INC 			:=	-I$(MPICH_INC)

INC_COMMON		:=	-I$(PackageDir_src)/include/common/obj \
					-I$(PackageDir_src)/include/common/src

INC_TARGET 		:=	\
					-I$(PackageDir_src)/include/wavefield_extrapolation/obj \
					-I$(PackageDir_src)/include/wavefield_extrapolation/src \
					-I$(PackageDir_src)/obj -I./obj

INCLUDE	    	:= 	 $(INC_TARGET) $(INC_COMMON) $(INC)

lib_package_obj	:=	-L$(PackageDir_src)/lib
LD_COMMON 		:= 	-L./lib 

lib_include_obj :=	\
					$(lib_package_obj) -lwavefield_extrasub \
					$(lib_package_obj) -linclude_wclsub

lib_common_obj	:=	$(lib_package_obj) $(LD_COMMON) $(lib_include_obj)

lib_target_obj	:=	$(LD_COMMON) -l$(target_name)main \
					$(LD_COMMON) -l$(target_name)sub \
					$(lib_common_obj)

FLAG_INCLUDE	:= $(INCLUDE) $(lib_common_obj)
FLAG_TARGET		:= $(INCLUDE) $(lib_target_obj)
f90_src := $(wildcard $(SRC_DIR)/*.f90)
f_src := $(wildcard $(SRC_DIR)/f/*.f)
c_src := $(wildcard $(SRC_DIR)/*.c)
cpp_src := $(wildcard $(SRC_DIR)/c/*.cpp)
include $(PackageDir_src)/config/Makefile_Rule
sub_lib : $(sub_obj) 
	$(AR) $(ARFUNCT) $(LIB_DIR)/$(lib_sub_name) $(sub_obj)
	cp $(LIB_DIR)/$(lib_sub_name) $(PackageDir_src)/lib

main_lib : $(main_obj)
	$(AR) $(ARFUNCT) $(LIB_DIR)/$(lib_main_name) $(main_obj)
	cp $(LIB_DIR)/$(lib_main_name) $(PackageDir_src)/lib

target	:=	
ifneq ($(f90_obj_dir_main), )
target	+=	$(target_f)
$(target_f) : $(f90_obj_dir_main)
	$(FC) -o $(target_f) $(f90_obj_dir_main) $(F90FLAG) $(FLAG_TARGET) $(LD_F90)

endif

ifneq ($(c_obj_dir_main), )
target	+=	$(target_c)
$(target_c) : $(main_obj)
	$(CC) -o $(target_c) $(c_obj_dir_main) $(CFLAG) $(FLAG_TARGET) $(LD_C)

endif

bin : $(target)
	@echo "bin: target", $(target)

all :
	cd $(PackageDir_src)/include/common; $(MAKE) sub_lib
	cd $(PackageDir_src)/include/wavefield_extrapolation; $(MAKE) sub_lib
	cd ./; $(MAKE) sub_lib
	cd ./; $(MAKE) main_lib
	cd ./; $(MAKE) bin
	mv $(target) $(PackageDir_src)/bin

.PHONY: clean
clean :	
	-rm $(OBJ_DIR)/*.o  $(OBJ_DIR)/*.mod
	-rm $(LIB_DIR)/*.a	
	-rm $(PackageDir_src)/lib/$(lib_sub_name)
	-rm $(PackageDir_src)/lib/$(lib_main_name)
	-rm $(target)
	-rm $(PackageDir_src)/bin/$(target)

.PHONY: cleanall
cleanall :	
	cd $(PackageDir_src)/include/common; $(MAKE) clean
	cd $(PackageDir_src)/include/wavefield_extrapolation; $(MAKE) clean
	-rm $(OBJ_DIR)/*.o  $(OBJ_DIR)/*.mod
	-rm $(LIB_DIR)/*.a	
	-rm $(PackageDir_src)/lib/$(lib_sub_name)
	-rm $(PackageDir_src)/lib/$(lib_main_name)
	-rm $(target)
	-rm $(PackageDir_src)/bin/$(target)

