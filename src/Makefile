# Makefile for RAJZEL

#--------------------------------------------------------
# edit here:

# source code for model generation

MODEL = 1D_start_BP2004.c
MODEL_GRID_SEARCH = 1D_start_gridsearch.c
EXEC= ../bin

# Compiler (LAM: CC=hcc, CRAY T3E: CC=cc)

# ON Linux cluster running LAM
#CC=hcc
#LFLAGS=-lm -lmpi 
#CFLAGS=-Wall -O4 

# On CRAY T3E
# CC=cc

# On Linux NEC-cluster with Intel-MPI
CC=mpiicc 
LFLAGS=-lm -ftz -lstdc++ -lumfpack
CFLAGS=-O3 -xAVX -fno-fnalias -restrict
SFLAGS=-L$//sfs/fs4/work-sh2/sungw331/suitesparse/SuiteSparse/lib
IFLAGS=-I$//sfs/fs4/work-sh2/sungw331/suitesparse/SuiteSparse/include

# On Desktop computer with LinuxMint 17, OpenMPI and gcc 4.8.2
#CC=mpicc
#LFLAGS=-lm -lfftw3 -lstdc++ -lumfpack
#CFLAGS=-O3 -w -fno-stack-protector -D_FORTIFY_SOURCE=0

# On HLRN system
#CC=mpcc
#LFLAGS=-lm  

# ALTIX
#CC=icc
#CFLAGS=-mp -O3 -ipo
#LFLAGS=-lmpi -lm -i-static

# after this line, no further editing should be necessary
# --------------------------------------------------------

.c.o:
	$(CC) $(CFLAGS) -c $< $(IFLAGS)

FDVEPS_UTIL = \
	absorb.c \
	exchange_par.c \
	info.c \
	initproc.c \
	$(MODEL_EL) \
	note.c \
	outseis.c \
	readdsk.c \
	read_par.c \
	readmod.c \
    	receiver.c \
	sources.c \
	util.c \
	writemod.c \
	write_par.c \
	writedsk.c 

GERMAINE= \
	alloc_fwiAC.c \
        alloc_matAC.c \
	alloc_seis_AC.c \
        alloc_waveAC.c \
	apply_hess_AC.c \
	ass_grad_AC.c \
	calc_res_AC.c \
        calc_mat_change.c \
        calc_mat_change_wolfe.c \
        calc_opt_step.c \
	calc_opt_step_test.c \
	calc_seis_AC.c \
        check_descent.c \
        cp_grad_frame.c \
	cp_vec.c \
        descent.c \
	dotp.c \
        dotp_matrix.c \
        exchange_grad_MPI.c \
	median.c \
        forward_AC.c \
        forward_shot_AC.c \
	free_seis_AC.c \
	fwi_FD_AC.c \
	gauss_filt.c \
	germaine.c \
	grad_obj_AC.c \
        grid_search.c \
	hess_shin_AC.c \
	info.c \
        info_mem.c \
        init_A_AC_9p_pml.c \
        init_grad.c \
        init_mat_AC.c \
        init_MPIshot.c \
        LBFGS.c \
	smooth2.c \
	smooth_model.c \
	smooth_grad.c \
	$(MODEL) \
        $(MODEL_GRID_SEARCH) \
	median_src.c \
	model_out.c \
        calc_nonzero.c \
	note.c \
	norm.c \
        norm1.c \
        norm_matrix.c \
	obj_AC.c \
        parabolicls_AC.c \
	PCG.c \
	pml_pro.c \
        precond.c \
	readdsk.c \
	read_par.c \
	read_par_inv.c \
	readmod.c \
    	receiver.c \
        RHS_source_AC.c \
        RHS_source_AC_adj.c \
	sources.c \
	solvelin.c \
        sum_grad_MPI.c \
        store_mat.c \
	taper_grad.c \
	taper_grad_hor.c \
	taper_grad_shot.c \
        taper_grad_shot1.c \
	util.c \
	vec2mat.c \
        wolfels.c \
        write_gridsearch.c \
	writemod.c \
	writemod_vec.c \
	writemod_true.c \
	write_par.c \
	write_seis_AC.c \
	writedsk.c \
	zero_LBFGS.c \
        zero_LBFGS1.c

GERMAINE_OBJ = $(GERMAINE:%.c=%.o)

germaine:		$(GERMAINE_OBJ)
	$(CC) $(SFLAGS) $(GERMAINE_OBJ) -o $(EXEC)/germaine $(LFLAGS)

clean:
	find . -name "*.o" -exec rm {} \; 
	find . -name "*.c%" -exec rm {} \;
	find . -name "*.bck" -exec rm {} \;
