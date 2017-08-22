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
#CC=mpiicc 
#LFLAGS=-lm -ftz -lstdc++ -lumfpack
#CFLAGS=-O3 -xAVX -fno-fnalias -restrict
#SFLAGS=-L$//sfs/fs4/work-sh2/sungw331/suitesparse/SuiteSparse/lib
#IFLAGS=-I$//sfs/fs4/work-sh2/sungw331/suitesparse/SuiteSparse/include -I./../include

# On Desktop computer with LinuxMint 17, OpenMPI and gcc 4.8.2
CC=mpicc
LFLAGS=-lm -lfftw3 -lstdc++ -lumfpack
CFLAGS=-O3 -w -fno-stack-protector -D_FORTIFY_SOURCE=0
SFLAGS=-L$//home/daniel-koehn/libs/SuiteSparse/lib
IFLAGS=-I$//home/daniel-koehn/libs/SuiteSparse/include -I./../include

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

VPATH = AC:SH:TE:models

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

GERMAINE := $(wildcard *.c)

# add AC source files
GERMAINE+= \
	alloc_fwiAC.c \
        alloc_matAC.c \
	alloc_seis_AC.c \
	alloc_seis_fwi_AC.c \
        alloc_waveAC.c \
	apply_hess_AC.c \
	ass_grad_AC.c \
	calc_res_AC.c \
	calc_seis_AC.c \
	extract_PCG_AC.c \
        forward_AC.c \
        forward_shot_AC.c \
	free_seis_AC.c \
	fwi_FD_AC.c \
	grad_obj_AC.c \
	hessian_AC.c \
	hessian_shin_AC.c \
        init_A_AC_9p_pml.c \
        init_mat_AC.c \
	$(MODEL) \
        $(MODEL_GRID_SEARCH) \
	obj_AC.c \
        parabolicls_AC.c \
	read_seis_AC.c \
        RHS_source_AC.c \
        RHS_source_AC_adj.c \
        RHS_source_AC_hess.c \
	RTM_FD_AC.c \
	RTM_AC_out.c \
	store_PCG_AC.c \
	wolfels_AC.c \
	write_seis_AC.c

# add SH source files
GERMAINE+= \
	alloc_matSH.c \
	forward_SH.c \
	forward_shot_SH.c \
	init_A_SH_9p_pml.c \
	pml_pro_SH.c \
	readmod_SH.c

# add TE source files
GERMAINE+= \
	alloc_fwiTE.c \
	alloc_matTE.c \
	ass_grad_TE.c \
	extract_PCG_TE.c \
	forward_TE.c \
	forward_shot_TE.c \
	fwi_FD_TE.c \
	grad_obj_TE.c \
	hessian_shin_TE.c \
	hessian_TE.c \
	init_A_TE_9p_pml.c \
	model_out_TE.c \
	obj_TE.c \
	parabolicls_TE.c \
	readmod_TE.c \
	store_PCG_TE.c \
	store_PCG_TE_old.c

GERMAINE_OBJ = $(GERMAINE:%.c=%.o)

germaine:		$(GERMAINE_OBJ)
	$(CC) $(SFLAGS) $(GERMAINE_OBJ) -o $(EXEC)/germaine $(LFLAGS)

clean:
	find . -name "*.o" -exec rm {} \; 
	find . -name "*.c%" -exec rm {} \;
	find . -name "*.bck" -exec rm {} \;
