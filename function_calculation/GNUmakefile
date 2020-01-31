
objects = function_calculation.o special_functions.o interpolation.o

function_calculation : $(objects)
		g++ -O2 -fopenmp -o function_calculation $(objects)

function_calculation.o : function_calculation.cpp special_functions.h integration.h utils.h
		g++ -O2 -fopenmp -c function_calculation.cpp
		
interpolation.o : interpolation.cpp interpolation.h
		g++ -O2 -fopenmp -c interpolation.cpp

special_functions.o : special_functions.cpp special_functions.h
		g++ -O2 -fopenmp -c special_functions.cpp
		
.PHONY :clean
clean :
	rm function_calculation $(objects)
