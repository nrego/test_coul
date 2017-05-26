#!/bin/bash
rm \#*
rm step*
grompp -f params_run.mdp -c methanol.gro -p top_simple.top -maxwarn 1 -o simple.tpr
mdrun -deffnm simple -rerun traj.trr -reprod -nt 1
echo 5 6 0 | g_energy -f simple.edr -o simple.xvg
