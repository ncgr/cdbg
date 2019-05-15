CC     := g++
AR     := ar
SRCDIR := src

BUILDDIR := build
TARGET   := lib/libcdbg.a
 
SRCEXT  := cpp
SOURCES := $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))
OBJDIRS := $(dir $(OBJECTS))
CFLAGS  := -c -O9 -std=c++11 -DNDEBUG -ffast-math -funroll-loops -msse4.2 -Wall
MACROS  := -DWTBV -DBV1BV -DBV3IL
#MARCOS  := -DWTBV -DBV1SD -DBV3SD
#MARCOS  := -DWTRRR15 -DBV1SD -DBV3SD
INC     := -Iinclude
AFLAGS  := rcs

$(TARGET): $(OBJECTS)
	@echo " Linking..."
	@echo " $(AR) $(AFLAGS) $(TARGET) $^"; $(AR) $(AFLAGS) $(TARGET) $^

$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)
	@mkdir -p $(OBJDIRS)
	@echo " $(CC) $(CFLAGS) $(MACROS) $(INC) -c -o $@ $<"; $(CC) $(CFLAGS) $(MACROS) $(INC) -c -o $@ $<

clean:
	@echo " Cleaning..."
	@echo " $(RM) -r $(BUILDDIR) $(TARGET)"; $(RM) -r $(BUILDDIR) $(TARGET)

.PHONY: clean
