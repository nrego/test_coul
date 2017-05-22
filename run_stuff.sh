#!/bin/bash

rm \#*
rm *.tpr *.xvg *.trr *.edr
grompp -f params_run.mdp -c struct_single.gro -p top_run.top -maxwarn 2 -o run.tpr
mdrun -s run.tpr -v

echo 0 | trjconv -s run.tpr -f traj.trr -o traj_pbc.trr -pbc mol

grompp -f params_no_pme.mdp -c struct_single.gro -p top.top -maxwarn 2 -o no_pme.tpr
mdrun -deffnm no_pme -rerun traj.trr -reprod -nt 1

grompp -f params_pme.mdp -c struct_single.gro -p top.top -maxwarn 2 -o pme.tpr
mdrun -deffnm pme -rerun traj.trr -reprod -nt 1

grompp -f params_ewald.mdp -c struct_single.gro -p top.top -maxwarn 2 -o ewald.tpr
mdrun -deffnm ewald -rerun traj.trr -reprod -nt 1

echo 4 0 | g_energy -f pme.edr -o pme_energy.xvg
echo 3 0 | g_energy -f no_pme.edr -o no_pme_energy.xvg
echo 4 0 | g_energy -f ewald.edr -o ewald_energy.xvg
echo 2 0 | g_energy -f ewald.edr -o ewald_sr.xvg
echo 3 0 | g_energy -f ewald.edr -o ewald_lr.xvg

python calc_coul.py

rm \#*
