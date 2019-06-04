APP = simplify

OBJECTS = main.o
MESHES = meshes/*
OUTPUTS = outputs/*

CXXC = g++-9
CXXFLAGS = -std=c++11 -O3
LINKFLAGS = -O3

COMMENTS = ''
DEL = rm -rf

default:
	make $(APP)

$(APP): $(OBJECTS) Makefile
	echo 'Linking: $(APP)' && \
	$(CXXC) $(LINKFLAGS) $(OBJECTS) -o $(APP)

$(APP).o: $(APP).cpp Makefile
	echo 'Compiling: $(APP).o:'	&& \
	$(CXXC) $(CXXFLAGS) -c $*.cpp -o $*.o

%.o: %.cpp %.h Makefile
	echo 'Compiling: $*.o' && \
	$(CXXC) $(CXXFLAGS) -c $*.cpp -o $*.o

run:
	make
	./$(APP)

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