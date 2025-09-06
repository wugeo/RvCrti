PackageDir_src :=	./
include $(PackageDir_src)/config/Makefile_inc
Debug_symbols	:=	FALSE
include $(PackageDir_src)/config/Makefile_compile
CUR_DIR		= .
SRC_DIR		= $(CUR_DIR)/src
OBJ_DIR		= $(CUR_DIR)/obj
INC_DIR		= $(CUR_DIR)/include
LIB_DIR		= $(CUR_DIR)/lib
BIN_DIR		= $(CUR_DIR)/bin



all:
	cd $(SRC_DIR); $(MAKE) all


.PHONY:clean
clean:
	cd $(SRC_DIR)/rv_crti2d; $(MAKE) clean
	-rm $(LIB_DIR)/*.a	
	-rm $(BIN_DIR)/*





