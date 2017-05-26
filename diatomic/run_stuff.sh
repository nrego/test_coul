#!/bin/bash
rm \#*
rm step*
rm *.xvg
grompp -f params_run.mdp -c diatomic.gro -p top.top -o run.tpr -maxwarn 1
mdrun -s run.tpr -v 

grompp -f params_run.mdp -c diatomic.gro -p top_excl.top -o excl.tpr -maxwarn 1
mdrun -deffnm excl -rerun traj.trr -reprod -nt 1 -v

grompp -f params_run.mdp -c diatomic.gro -p top_noexcl.top -o no_excl.tpr -maxwarn 1
mdrun -deffnm no_excl -rerunt traj.trr -reprod -nt 1 -v 

# Output PME SR and LR
echo 3 4 0 | g_energy -f excl.edr -o ener_excl.xvg
echo 3 4 0 | g_energy -f no_excl.edr -o ener_no_excl.xvg

rm \#*


