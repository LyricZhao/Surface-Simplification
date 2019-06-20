APP = simplify
ENTRANCE = main

OBJECTS = main.o
HPPS = mesh.hpp
MESHES = meshes/*
OUTPUTS = outputs/*

CXXC = g++-9
CXXFLAGS = -std=c++11 -O3
LINKFLAGS = -O3

COMMENTS = ''
TEST_OBJ = horse.fine.90k.obj
TEST_RATIO = 0.001
TEST_T = 0.01
DEL = rm -rf

default:
	make $(APP)

$(APP): $(OBJECTS) Makefile
	echo 'Linking: $(APP)' && \
	$(CXXC) $(LINKFLAGS) $(OBJECTS) -o $(APP)

$(ENTRANCE).o: $(ENTRANCE).cpp $(HPPS) Makefile
	echo 'Compiling: $(ENTRANCE).o:'	&& \
	$(CXXC) $(CXXFLAGS) -c $*.cpp -o $*.o

%.o: %.cpp %.h Makefile
	echo 'Compiling: $*.o' && \
	$(CXXC) $(CXXFLAGS) -c $*.cpp -o $*.o

run:
	make
	time ./$(APP) meshes/$(TEST_OBJ) outputs/$(TEST_OBJ) $(TEST_RATIO) $(TEST_T)

clean:
	make clean_objs

clean_objs:
	echo 'Cleaning objects ...'
	-$(DEL) *.o
	-$(DEL) $(APP)

clean_outputs:
	echo 'Cleaning outputs ...'
	-$(DEL) $(OUTPUTS)

push:
	echo 'Comments: $(COMMENTS)'
	git add .gitignore
	git add *
	git commit -m "$(COMMENTS)"
	git push origin master