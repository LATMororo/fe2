.DEFAULT_GOAL := extension

include python.mk

program   = material_bridge

OBJDIR    = OBJ_EXTENSION
extension = _fe2_material_native$(PYTHON_EXTENSION)
VPATH     = ../src
objects   = $(filter-out main.cpp PRNNMaterial.cpp, \
			$(notdir $(wildcard ../src/*.cpp))) # exclude such source files 

# Flags for Jem/Jive makefiles.

MY_CXX_STD_FLAGS = -std=c++17 -fPIC -Wno-unused-local-typedefs -Werror \
				   -DFE2_PYTHON_EXTENSION $(PYTHON_INCLUDES)
MY_LIBDIRS       = $(PYTHON_LIBDIR)
MY_LIBS          = python$(PYTHON_VERSION)
MY_CXX_LNK_FLAGS = -shared -Wl,-rpath,$(PYTHON_LIBDIR)

include $(JEMDIR)/packages/mt/Makedefs.mk
include $(JIVEDIR)/makefiles/prog.mk
include $(JIVEDIR)/makefiles/packages/*.mk

.PHONY: extension clean-all clean-shared print-python-vars
extension: $(extension)
$(extension): $(objects)
	$(LinkExe)

clean-all: clean-shared  # add one more rule for `make clean-all'

clean-shared:
	rm -f $(extension)

print-python-vars: # print python and pybind11 related variables
	$(info PYTHON_EXECUTABLE is $(PYTHON_EXECUTABLE))
	$(info PYTHON_INCLUDE    is $(PYTHON_INCLUDE))
	$(info PYBIND_INCLUDE    is $(PYBIND_INCLUDE))
	$(info PYTHON_VERSION    is $(PYTHON_VERSION))
	$(info PYTHON_LIBDIR     is $(PYTHON_LIBDIR))
	$(info PYTHON_EXTENSION  is $(PYTHON_EXTENSION))
	$(info PYTHON_INCLUDES   is $(PYTHON_INCLUDES))
