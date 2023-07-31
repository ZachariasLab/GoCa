CXX = g++
FLAGS = -O3 -I src/include/Eigen -fopenmp# -g -fsanitize=address -Wall
LDFLAGS = -fopenmp

EXECUTABLE = GoCa
SOURCEDIR = src
BUILDDIR = build
SOURCES = $(wildcard $(SOURCEDIR)/*.cpp)
DEPENDS = $(patsubst $(SOURCEDIR)/%.cpp, $(BUILDDIR)/%.d, $(SOURCES))
OBJECTS = $(patsubst $(SOURCEDIR)/%.cpp, $(BUILDDIR)/%.o, $(SOURCES))


.PHONY: all clean

all: dir $(EXECUTABLE)

dir:
	$(shell mkdir -p $(BUILDDIR))

$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(LDFLAGS) $^ -o $@

-include $(DEPENDS)

$(OBJECTS): $(BUILDDIR)/%.o : $(SOURCEDIR)/%.cpp makefile
	$(CXX) $(FLAGS) -MMD -MP -c $< -o $@

clean:
	rm -r $(BUILDDIR) $(EXECUTABLE)