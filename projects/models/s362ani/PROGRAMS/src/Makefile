SRC = ../src
OBJ = ../obj
BIN = ../bin

#--------------------------------------------------------------------

.f.o:
	 gfortran -c $*.f
 
all: 	velani.exe \
	topo.exe \
        test_subshsv.exe

#--------------------------------------------------------------------

LIB_velani =     $(OBJ)/velani.o \
                 $(OBJ)/evradker.o \
                 $(OBJ)/vbspl.o \
                 $(OBJ)/splconsngl.o \
                 $(OBJ)/ylm.o \
                 $(OBJ)/gt3dmodl.o \
                 $(OBJ)/rd3dmodl.o \
                 $(OBJ)/putcurrmodels.o \
                 $(OBJ)/getcurrmodels.o \
                 $(OBJ)/deck_subs.o \
                 $(OBJ)/drspln.o \
                 $(OBJ)/drsple.o 

velani.exe: $(BIN)/velani.exe

$(BIN)/velani.exe:	$(LIB_velani)
	gfortran $(LIB_velani) -o $(BIN)/velani.exe

$(OBJ)/velani.o: $(SRC)/velani.f
	gfortran -c $(SRC)/velani.f -o $@

$(OBJ)/splconsngl.o: $(SRC)/splconsngl.f
	gfortran -c $(SRC)/splconsngl.f -o $@

$(OBJ)/splgradcon.o: $(SRC)/splgradcon.f
	gfortran -c $(SRC)/splgradcon.f -o $@

$(OBJ)/ylm.o: $(SRC)/ylm.f
	gfortran -c $(SRC)/ylm.f -o $@

$(OBJ)/ylmgrad.o: $(SRC)/ylmgrad.f
	gfortran -c $(SRC)/ylmgrad.f -o $@

$(OBJ)/evradker.o: $(SRC)/evradker.f
	gfortran -c $(SRC)/evradker.f -o $@

$(OBJ)/vbspl.o: $(SRC)/vbspl.f
	gfortran -c $(SRC)/vbspl.f -o $@

$(OBJ)/rd3dmodl.o: $(SRC)/rd3dmodl.f $(SRC)/3dmodl.h
	gfortran -c $(SRC)/rd3dmodl.f -o $@

$(OBJ)/gt3dmodl.o: $(SRC)/gt3dmodl.f $(SRC)/3dmodl.h
	gfortran -c $(SRC)/gt3dmodl.f -o $@

$(OBJ)/putcurrmodels.o: $(SRC)/putcurrmodels.f
	gfortran -c $(SRC)/putcurrmodels.f -o $@

$(OBJ)/getcurrmodels.o: $(SRC)/getcurrmodels.f
	gfortran -c $(SRC)/getcurrmodels.f -o $@

$(OBJ)/deck_subs.o: $(SRC)/deck_subs.f
	gfortran -c $(SRC)/deck_subs.f -o $@

$(OBJ)/drspln.o: $(SRC)/drspln.f
	gfortran -c $(SRC)/drspln.f -o $@

$(OBJ)/drsple.o: $(SRC)/drsple.f
	gfortran -c $(SRC)/drsple.f -o $@


#--------------------------------------------------------------------

LIB_topo =       $(OBJ)/topo.o \
                 $(OBJ)/evradker.o \
                 $(OBJ)/vbspl.o \
                 $(OBJ)/splconsngl.o \
                 $(OBJ)/ylm.o \
                 $(OBJ)/gt3dmodl.o \
                 $(OBJ)/rd3dmodl.o \
                 $(OBJ)/putcurrmodels.o \
                 $(OBJ)/getcurrmodels.o 

topo.exe: $(BIN)/topo.exe

$(BIN)/topo.exe:	$(LIB_topo)
	gfortran $(LIB_topo) -o $(BIN)/topo.exe

$(OBJ)/topo.o: $(SRC)/topo.f
	gfortran -c $(SRC)/topo.f -o $@

#--------------------------------------------------------------------

LIB_test_subshsv =   	$(OBJ)/test_subshsv.o \
                        $(OBJ)/subshsv.o \
                        $(OBJ)/subtopo.o \
                        $(OBJ)/evradker.o \
                        $(OBJ)/vbspl.o \
                        $(OBJ)/splconsngl.o \
                        $(OBJ)/ylm.o \
                        $(OBJ)/gt3dmodl.o \
                        $(OBJ)/rd3dmodl.o \
                        $(OBJ)/putcurrmodels.o \
                        $(OBJ)/getcurrmodels.o 

test_subshsv.exe: $(BIN)/test_subshsv.exe

$(BIN)/test_subshsv.exe: $(LIB_test_subshsv)
	gfortran $(LIB_test_subshsv) -o $(BIN)/test_subshsv.exe

$(OBJ)/test_subshsv.o: $(SRC)/test_subshsv.f
	gfortran -c $(SRC)/test_subshsv.f -o $@

$(OBJ)/subshsv.o: $(SRC)/subshsv.f $(SRC)/mod.h
	gfortran -c $(SRC)/subshsv.f -o $@

$(OBJ)/subtopo.o: $(SRC)/subtopo.f $(SRC)/mod.h
	gfortran -c $(SRC)/subtopo.f -o $@

