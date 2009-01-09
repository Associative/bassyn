# Make sure to set the PLATFORM environment variable to either 'linux' or 'osx'
# e.g.: $ PLATFORM=osx make

CC = gcc
CFLAGS = -O2 -std=c99 -pedantic -Wall `sdl-config --cflags`
LDFLAGS = -lm
OBJ = main.o apf.o apf-graph.o
DEMO = bassyn
TEST = 
PLATFORM = osx

# We assume a fink + sdl + xcode utilities setup on OS X
# Also note that we do a static link to SDL
ifeq (${PLATFORM}, osx)
	CFLAGS += -I/System/Library/Frameworks/GLUT.framework/Versions/A/Headers/
	LDFLAGS += `sdl-config --static-libs` -framework OpenGL -framework GLUT 
else
	LDFLAGS += `sdl-config --libs` -lGL -lGLU -lglut 
endif

test: ${DEMO}
	rm *.o
	./$(DEMO) ${TEST}

gt: graph
	rm *.o
	./graph

graph: apf.o apf-graph.o graphtest.o
	$(CC) -o graph apf.o apf-graph.o graphtest.o $(LDFLAGS)

bassyn: $(OBJ)
	$(CC) -o ${DEMO} $(OBJ) $(LDFLAGS)

%.o: %.c apf.h
	$(CC) $(CFLAGS) -c $<

clean:
	rm *.o ${DEMO}
