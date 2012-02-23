#$ -S /bin/bash
#$ -cwd
#$ -pe 8way 256
#$ -l h_rt=03:00:00
#$ -j y
#$ -N gadget2CR
#$ -o /work/00863/minerva/out.gadget2CR
#$ -e /work/00863/minerva/err.gadget2CR
#$ -q normal
#$ -V
#$ -A A-astro
## -M minerva@astro.as.utexas.edu
## -m be
   set -x
   ibrun /share/home/00863/minerva/programs/gadget2CR/g2 par.txt 2
