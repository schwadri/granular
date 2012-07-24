DEPTH = ../../../../..

include $(DEPTH)/make/openclsdkdefs.mk 

####
#
#  Targets
#
####

OPENCL			= 1
SAMPLE_EXE		= 1
INDEPENDENT_SAMPLE = 1
EXE_TARGET 		= Reduction
EXE_TARGET_INSTALL   	= Reduction

####
#
#  C/CPP files
#
####

vpath %.cpp : ../../../SDKUtil/

FILES 		= Reduction \
			  SDKApplication \
              SDKBitMap \
	          SDKCommon \
	          SDKCommandArgs \
	          SDKFile \
	          SDKThread

CLFILES	= Reduction_Kernels.cl

include $(DEPTH)/make/openclsdkrules.mk 

