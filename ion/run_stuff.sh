#!/bin/bash
rm \#*
rm *.{edr,xvg,tpr}

grompp -f params_run.mdp -c ion.gro -p top.top -o run.tpr -maxwarn 1
mdrun -s run.tpr -v 

grompp -f params_fe.mdp -c ion.gro -p top_fe.top -o fe_run.tpr -maxwarn 1
mdrun -deffnm fe_run -rerun traj.trr -reprod -nt 1 -v

# Output PME SR and LR
echo 2 3 4 0 | g_energy -f fe_run.edr -o ener_fe.xvg

rm \#*


