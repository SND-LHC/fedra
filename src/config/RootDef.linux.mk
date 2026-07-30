ARCH         := $(shell root-config --arch)
#PLATFORM     := $(shell root-config --platform)
ObjSuf        = o
SrcSuf        = cxx
ExeSuf        =
DllSuf        = so
OutPutOpt     = -o

ROOTVERSION   = $(shell root-config --version | cut -b 1,2,3,4 | awk '{ if ($$1>5.22) print "NEWER522" }')
ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs)

# ifeq ($(ARCH),linux)
# Linux with egcs
CXX           = g++ 
CXXFLAGS      = -g -O2 -w -fPIC -DUSE_ROOT
LD            = g++
LDFLAGS       = -g
SOFLAGS       = -shared
# endif

CXXFLAGS     += $(ROOTCFLAGS)
LIBS          = $(ROOTLIBS)
GLIBS         = $(ROOTGLIBS) -lGeom
