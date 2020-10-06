#paramètres
F902=gfortran
F90=mpif90
##FFLAGS=-o -wall -ffpe-trap=invalid -pedantic -g
#ou
DEBUG=-o -wall -ffpe-trap=invalid -fbounds-check -pedantic -fdefault-real-8 -g
COMPILE= mpif90 $(DEBUG)
COMPILE2=gfortran $(DEBUG)

#edition de liens creation executable nommé run
run : main.o ModMesh.o ModMatrice.o ModParra.o ModTri.o ModPGC.o
	$(COMPILE) -o run main.o ModMesh.o ModMatrice.o ModParra.o ModTri.o ModPGC.o
	@echo compilation terminée

#compilation
main.o : main.f90 ModParra.o  ModMesh.o ModMatrice.o ModTri.o ModPGC.o
	$(COMPILE) -c main.f90
	$(F90) -c main.f90
ModParra.o : ModParra.f90
	$(COMPILE) -c ModParra.f90
	$(F90) -c  ModParra.f90
ModMesh.o : ModMesh.f90
	$(COMPILE) -c ModMesh.f90
	$(F90) -c ModMesh.f90
ModMatrice.o : ModMatrice.f90 ModMesh.o
	$(COMPILE) -c ModMatrice.f90
	$(F90) -c ModMatrice.f90
ModTri.o : ModTri.f90 ModMesh.o
	$(COMPILE) -c ModTri.f90
	$(F90) -c ModTri.f90
ModPGC.o : ModPGC.f90 ModMesh.o
	$(COMPILE) -c ModPGC.f90
	$(F90) -c ModPGC.f90
#destrution des fichiers objets et modules
clean:
	rm -f -- -wall ./-wall.o *.o *.dat *.mod  *~
#destruction des .txt
cleantxt:
	rm -f *.txt
#destruction des .txt dans les repertoires
cleanreptxt:
	rm -r ./DEF/*.txt ./DPL/*.txt ./MODULE_EVO/*.txt ./RES/*.txt ./SIGMA/*.txt ./TPS_CPU/*.txt
