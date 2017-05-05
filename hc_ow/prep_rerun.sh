#!/bin/bash

rm -f \#* fe*

grompp -f params_fe.mdp -c struct_single.gro -p top.top -maxwarn 2 -o fe.tpr
mdrun -deffnm fe -rerun traj.trr -nt 1 -reprod -v

rm \#*

