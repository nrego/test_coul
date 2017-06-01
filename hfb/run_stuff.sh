#!/bin/bash
rm \#*
rm step*
rm *.{edr,xvg,tpr}

grompp -f params_run.mdp -c diatomic.gro -p top.top -o run.tpr 
mdrun -s run.tpr -v 

grompp -f params_fe.mdp -c diatomic.gro -p top_excl.top -o excl.tpr -maxwarn 1
mdrun -deffnm excl -rerun traj.trr -reprod -nt 1 -v

grompp -f params_fe.mdp -c diatomic.gro -p top_pairs.top -o pairs.tpr -maxwarn 1
mdrun -deffnm pairs -rerun traj.trr -reprod -nt 1 -v

grompp -f params_fe.mdp -c diatomic.gro -p top_noexcl.top -o no_excl.tpr -maxwarn 2
mdrun -deffnm no_excl -rerunt traj.trr -reprod -nt 1 -v 

# Output PME SR and LR
echo 3 4 0 | g_energy -f excl.edr -o ener_excl.xvg
# Also include Coul (14)
echo 3 5 6 0 | g_energy -f pairs.edr -o ener_pairs.xvg
echo 3 4 0 | g_energy -f no_excl.edr -o ener_no_excl.xvg

rm \#*


