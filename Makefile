
IFORT_FLAGS=-c -fpp -threads -recursive -DCLI_ONLY -fpe1 -O3

# OPENMP not working properly, disabling these flags for now
# OPENMP_FLAGS=-parallel -fopenmp -mkl  
PREPROCESS_DEFS=preprocessor_definitions.fpp
MODULE_SOURCES=w2modules.f90

# original sources may need to be renamed using: rename 'y/ /_/'

SOURCES= \
        CEMA_Bubbles_Code_01.f90 \
        CEMA_FFT_Layer_01.f90 \
        CEMA_Input_01.f90 \
        CEMA_Input_Files_Read_01.f90 \
        CEMA_Output_01.f90 \
        CEMA_Sediment_Flux_Model_04.f90 \
        CEMA_Sediment_Model_03.f90 \
        CEMA_Turbidity_01.f90 \
        aerate.f90 \
        az.f90 \
        balances.f90 \
        date.f90 \
        density.f90 \
        endsimulation.f90 \
        envir_perf.f90 \
        fishhabitat.f90 \
        gas-transfer.f90 \
        gate-spill-pipe.f90 \
        heat-exchange.f90 \
        hydroinout.f90 \
        init-cond.f90 \
        init-geom.f90 \
        init-u-elws.f90 \
        init.f90 \
        input.f90 \
        layeraddsub.f90 \
        macrophyte-aux.f90 \
        output.f90 \
        outputa2w2tools.f90 \
        outputinitw2tools.f90 \
        particle.f90 \
        restart.f90 \
        screen_output_intel.f90 \
        shading.f90 \
        tdg.f90 \
        temperature.f90 \
        time-varying-data.f90 \
        transport.f90 \
        update.f90 \
        w2_4_win.f90 \
        water-quality.f90 \
        waterbody.f90 \
        withdrawal.f90 \
        wqconstituents.f90

preprocess_def_objects=$(PREPROCESS_DEFS:.fpp=.o)
objects=$(SOURCES:.f90=.o)


.PHONY: renames
renames:
	rename 's/\.F90$$/.f90/' *.F90
	rename 'y/ /_/' *.[fF]90

.PHONY: w2modules clean
w2modules: $(MODULE_SOURCES)
	ifort $(IFORT_FLAGS) -o w2modules.o $<

t1: w2modules
	echo "done with modules"

w2_exe_linux: w2modules $(preprocess_def_objects) $(objects)
	ifort -o w2_exe_linux $(objects) w2modules.o preprocessor_definitions.o


clean:
	rm -rf *.o *.mod

%.o : %.f90
	ifort $(IFORT_FLAGS) -o $@ $<

%.o : %.fpp
	ifort $(IFORT_FLAGS) -o $@ $<
