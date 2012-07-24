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
EXE_TARGET 		= PrefixSum
EXE_TARGET_INSTALL   	= PrefixSum

####
#
#  C/CPP files
#
####

vpath %.cpp : ../../../SDKUtil/

FILES 		= PrefixSum \
			  SDKApplication \
              SDKBitMap \
	          SDKCommon \
	          SDKCommandArgs \
	          SDKFile \
	          SDKThread

CLFILES	= PrefixSum_Kernels.cl

include $(DEPTH)/make/openclsdkrules.mk 

