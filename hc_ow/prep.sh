#!/bin/bash

rm -f *.tpr \#* traj* confout* ener* fe*

rm *.pdb

grompp -f params.mdp -c struct_single.gro -p top.top -o run.tpr -maxwarn 2
mdrun -s run.tpr -v


rm \#*
