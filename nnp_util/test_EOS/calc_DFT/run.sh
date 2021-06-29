#!/bin/sh
#PBS -N c-GeTe
#PBS -l nodes=1:g11:ppn=28

cd $PBS_O_WORKDIR
cat $PBS_NODEFILE > nodefile
NPROC=`wc -l < $PBS_NODEFILE`
VASP='/vasp'

for step in `seq 1 61`
do
  mkdir sample_$step
  cp {INCAR,KPOINTS,POTCAR} sample_$step/
  cd sample_$step
  mv ../POSCAR_T$step .
  cp POSCAR_T$step POSCAR
  mpirun -np $NPROC $VASP >& stdout.x
  rm vasprun.xml PROCAR DOSCAR POTCAR CHGCAR CHG EIGENVAL PCDAT XDATCAR WAVECAR REPORT
  cd ..
done
