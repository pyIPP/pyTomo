CFLAGS = -m64 -Wall -Wextra -Wformat=2 -std=gnu89 -O2 -march=native
SPQRLIBS = -L../../SuiteSparse/lib -lspqr -lsuitesparseconfig -lcholmod -lamd -lcolamd -lm -L$MKL_HOME/lib/intel64 
spqr_wrapper.so: spqr_wrapper.c
	gcc $^ -o $@ $(CFLAGS) -fPIC -shared $(SPQRLIBS)
clean:
	rm -f spqr_wrapper.so 
