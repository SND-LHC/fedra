#------------------------------------------------------------------------------
include $(FEDRA_ROOT)/src/config/ProjectDef.mk
include $(FEDRA_ROOT)/src/config/RootDef.mk
#------------------------------------------------------------------------------

TARGET	=	$(BIN_DIR)/viewdist$(ExeSuf) $(BIN_DIR)/viewdist_old$(ExeSuf) $(BIN_DIR)/viewsideal$(ExeSuf) $(BIN_DIR)/moslink$(ExeSuf) 
TARGET	+=	$(BIN_DIR)/mosalign$(ExeSuf) $(BIN_DIR)/mosalignbeam$(ExeSuf) $(BIN_DIR)/mostag$(ExeSuf) $(BIN_DIR)/tagalign$(ExeSuf) 
TARGET	+=	$(BIN_DIR)/tagtra$(ExeSuf) 
##TARGET	+=	$(BIN_DIR)/mosvol$(ExeSuf) 

ROOTLIBS      += -lMLP
PROJECT_LIBS += -L$(FEDRA_ROOT)/lib -lScan -lMosaic -lEdb -lEbase -lEmath -lEphys -lvt -lEdr -lScan -lEdr -lEIO -lEMC -lAlignment -lDataConversion -lTreePlayer -lScan -lEdb -lEbase -lEmath -lEphys -lvt -lEdr -lScan -lEdr -lEIO -lEMC -lAlignment -lTreePlayer $(ROOTLIBS)

#------------------------------------------------------------------------------
include ../../config/TargetsDef.mk
#--------------------------------------------------------------------------
