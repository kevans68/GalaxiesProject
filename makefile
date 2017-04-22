all: IV IVB galaxy

IV : 
	gcc -g -o randex randex.c nrutil.c
	
IVB : 
	gcc -g -o randexbulge randexbulge.c nrutil.c
	
galaxy:
	gcc -o galaxy Nbody_galaxy_1.c nrutil.c leapfrog_sa.c rk4_sa.c -O0 -lm
	
clean: 
	rm -f  randex randexbulge galaxy

