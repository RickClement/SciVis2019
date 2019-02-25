OBJECTS     = fluids.o Quad.cpp
CFILES      = $(OBJECTS:.o=.c)
CPPFILES    = $(OBJECTS:.o=.cpp)
EXECFILE    = smoke
INCLUDEDIRS = -I./fftw-2.1.5/include/
LIBDIRS     = -L./fftw-2.1.5/lib/
LIBS        = -lrfftw -lfftw -lglut -lGL -lGLU -lGLEW -lm
#Possible flags for release (ffast-math uses less precision for floating-point numbers, check that your application can handle this)
#CFLAGS      = -O3 -march=x86-64 -mtune=generic -DNDEBUG -mfpmath=sse -ffast-math -Wall -pipe
#Debug flags
CFLAGS      = -ggdb -Wall -pipe
LINKFLAGS   =


.SILENT:

all: $(EXECFILE)

$(EXECFILE): $(OBJECTS)
		g++ $(LINKFLAGS) $(OBJECTS) -o $(EXECFILE) $(LIBDIRS) $(LIBS)

.c.o: $$@.c $$@.h
		g++ $(CFLAGS) $(INCLUDEDIRS) -c $<

clean:
		-rm -rf $(OBJECTS) $(EXECFILE)

depend:
		g++ -MM $(CFILES) > make.dep
		g++ $(CPPFILES) > make.dep

-include make.dep
