SHELL    = /bin/sh

f90	 = gfortran
f90FLAGS = -O3 -march=native -ffast-math -funroll-loops -Wno-unused
LIB	 = -I/opt/netcdf/netcdf-3.6.3/include/ -L/opt/netcdf/netcdf-3.6.3/lib/ -lnetcdf -lstdc++ -I/usr/X11R6/include -Dcimg_use_xshm -Dcimg_use_xrandr -L/usr/X11R6/lib -lpthread -lX11 -lXext -lXrandr -I../NetCDF/include -lnetcdf_c++ -lnetcdf -L../NetCDF/lib/
LD	 = gfortran
LDFLAGS  = -O3 -march=native -ffast-math -funroll-loops -Wno-unused /opt/acml5.3.0/gfortran64/lib/libacml.a
#LDFLAGS  = -O3 -march=native -ffast-math -funroll-loops -Wno-unused /opt/acml/acml-4.4.0.gfortran/gfortran64/lib/libacml.a
SRC1	 = dns1_2d.f dnspr1.f rgg.f subdns_2d.f vfft1.f interp.f90 controle.cpp
OBJ1	 = dns1_2d.o dnspr1.o rgg.o subdns_2d.o vfft1.o interp.o controle.o

dns2d:	$(OBJ1)
	@echo
	@echo "--- Linking ( $@ ) ---"
	$(LD) -o $@ $(OBJ1) $(LIB) $(LDFLAGS)

dns2d1:	$(OBJ1)
	@echo
	@echo "--- Linking ( $@ ) ---"
	$(LD) $(LDFLAGS) -o $@ $(OBJ1) $(LIB)
	@mv dns2d1 TOTO


.SUFFIXES: .f90 .f .o 

.f.o:
	@echo
	@echo "--- Compiling ( $< ) ---"
	@echo	
	$(f90) $(f90FLAGS) -c $<

.f90.o:
	@echo
	@echo "--- Compiling ( $< ) ---"
	@echo	
	$(f90) $(f90FLAGS) -c $<

clean:
	@echo "--- Clean Garbage ---"
	-rm -f  *.o 

vclean:
	@echo "--- Clean Garbage ---"
	-rm -f  dns2d *.o fort.*

vfft1.o: vfft1.f
	@echo
	@echo "--- Compiling ( $< ) ---"
	@echo
	$(f90) $(f90FLAGS) -fdefault-real-8 -c vfft1.f

dns1_2d.o : dns1_2d.f par.f 
dnspr1.o : dnspr1.f 
rgg.o : rgg.f 
subdns_2d.o : subdns_2d.f par.f
#Write_Ncdf.o : Write_Ncdf.f
#Util_Ncdf.o : Util_Ncdf.f 
interp.o : interp.f90
toto.o : controle.cpp
