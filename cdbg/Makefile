CC     := g++ # This is the main compiler
SRCDIR := src

BUILDDIR := build
TARGET   := bin/cdbg
 
SRCEXT  := cpp
SOURCES := $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))
OBJDIRS := $(dir $(OBJECTS))
CFLAGS  := -O9 -std=c++11 -DNDEBUG -ffast-math -funroll-loops -msse4.2 -Wall
MACROS  := -DWTBV -DBV1BV -DBV3IL
#MARCOS  := -DWTBV -DBV1SD -DBV3SD
#MARCOS  := -DWTRRR15 -DBV1SD -DBV3SD
LIB     := -lsdsl -ldivsufsort -ldivsufsort64 -lcdbg
INC     := 

$(TARGET): $(OBJECTS)
	@echo " Linking..."
	@echo " $(CC) $^ -o $(TARGET) $(LIB)"; $(CC) $^ -o $(TARGET) $(LIB)

$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)
	@mkdir -p $(OBJDIRS)
	@echo " $(CC) $(CFLAGS) $(MACROS) $(INC) -c -o $@ $<"; $(CC) $(CFLAGS) $(MACROS) $(INC) -c -o $@ $<

clean:
	@echo " Cleaning..."
	@echo " $(RM) -r $(BUILDDIR) $(TARGET)"; $(RM) -r $(BUILDDIR) $(TARGET)

.PHONY: clean
