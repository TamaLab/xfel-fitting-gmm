###################################################
# this program is a part of GMM_DIFRACT_TOOLS. 
# Written by T. Nagai 
###################################################
CC=gcc
all:
	${CC}  -fPIC -Wall  -Q -O3 -mtune=native -c get_abs_struct_factor.c  -fopenmp -march=native  -ffast-math # with ffastmath ca. 20 % faster
	${CC}  -fPIC -shared -o get_abs_struct_factor.so get_abs_struct_factor.o  -lgomp
	${CC}  -fPIC -Wall  -Q -O3 -mtune=native -c Rotate_lib.c  -fopenmp -march=native 
	${CC}  -fPIC -shared -o Rotate_lib.so Rotate_lib.o 

clean:
	rm *.o  *so
