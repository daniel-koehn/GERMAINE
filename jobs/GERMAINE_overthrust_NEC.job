#!/bin/bash
#PBS -l elapstim_req=48:00:00  # Walltime
#PBS -l cputim_job=384:00:00   # akkumulierte CPU-Zeit pro Knoten
#PBS -l memsz_job=50gb         # Hauptspeicherbedarf
#PBS -b 1                      # Anzahl der Knoten
#PBS -T intmpi                 # gibt Jobtyp an; intmpi fuer Intel-MPI 
#PBS -l cpunum_job=16          # Anzahl benoetigter CPUs pro Knoten 
#PBS -N GERMAINE               # Name des Batch-Jobs
#PBS -o GERMAINE.out           # Datei fuer die Standardausgabe
#PBS -j o                      # Standard- und Fehlerausgabe in eine Datei 
#PBS -q clmedium               # Batch-Klasse

# Initialisierung der Intel-Umgebung
module load intel17.0.4 intelmpi17.0.4

cd $WORK/GERMAINE/par
mpirun $NQSII_MPIOPTS -np 16 ../bin/germaine GERMAINE_overthrust.inp GERMAINE_workflow_overthrust.inp > GERMAINE.out

# Ausgabe der verbrauchten Ressourcen (Rechenzeit, Hauptspeicher) nach Jobende
/opt/nec/nqsv/bin/qstat -f ${PBS_JOBID/0:/}
