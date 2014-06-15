all: sepSAM fragment Assembly Astroid

OPTFLAGS = -O3

#LDFLAGS = -static

CFLAGS += $(OPTFLAGS)

#CFLAGS += $(LDFLAGS)


sepSAM :
	g++ $(CFLAGS) -o bin/sepSAM src/sepSAM.cpp 2>log

fragment:
	g++ $(CFLAGS) -o bin/fragment src/parseFragment.cpp 2>log

Assembly:
	g++ $(CFLAGS) -o bin/Assembly src/general_functions.cpp src/Assembly.cpp 2>log
	

Astroid:
	g++ -std=c++0x -m64 -o Astroid src/Astroid.cpp 2>log

clean: