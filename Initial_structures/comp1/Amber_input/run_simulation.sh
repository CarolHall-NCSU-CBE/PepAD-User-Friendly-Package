#!/bin/bash
#BSUB -n 1
#BSUB -W 48:00
#BSUB -J amber
#BSUB -R "select[a100 || a30 || a10 || l40 || rtx2080 || gtx1080]"
#BSUB -gpu "num=1:mode=shared:mps=yes"
#BSUB -q gpu
#BSUB -o stdout.%J
#BSUB -e stderr.%J
pmemd_path=/rs1/researchers/h/hall/amber24/amber24_src/build/src/pmemd/src/pmemd.cuda_SPFP
sander_path=/rs1/researchers/h/hall/amber24/amber24_src/build/AmberTools/src/sander/sander
cpptraj_path=/rs1/researchers/h/hall/amber24/amber24_src/build/AmberTools/src/cpptraj/src/cpptraj

# module load amber
# module load cuda/12.0
# mkdir -p pdbfiles

tleap_file="tleap.in"
run_01="min_01"
run_02="heat_02"
run_03="equil_03"
run_04="equil_04"
run_05="min_05"
run_06="equil_06"
run_07="equil_07"
run_08="equil_08"
run_09="equil_09"
prod="prod"
center_file="center1.in"
cpptraj_file="getpdb.in"

pep="KLVFFAEini_12x2"
# saveamberparm pep KLVFFAEini_12x2_H2O.prmtop KLVFFAEini_12x2_H2O.inpcrd
module load amber
cat > "$tleap_file" <<EOF
    source leaprc.protein.ff14SB
    pep = loadpdb ${pep}.pdb
    source leaprc.water.tip3p
    solvateoct pep TIP3PBOX 9.0
    saveamberparm pep ${pep}_H2O.prmtop ${pep}_H2O.inpcrd

EOF

tleap -f "$tleap_file"
rm -f "$tleap_file"
module purge
prmtop_file="${pep}_H2O.prmtop"
init_inpcrd_file="${pep}_H2O.inpcrd"

$sander_path -O -i "${run_01}.in" -o "${run_01}.out" -p "$prmtop_file" -c "$init_inpcrd_file" -r "${run_01}.rst7" -ref "$init_inpcrd_file">> mdinfo
$pmemd_path -O -i "${run_02}.in" -o "${run_02}.out" -p "$prmtop_file" -c "${run_01}.rst7" -r "${run_02}.rst7" -ref "${run_01}.rst7" >> mdinfo
$pmemd_path -O -i "${run_03}.in" -o "${run_03}.out" -p "$prmtop_file" -c "${run_02}.rst7" -r "${run_03}.rst7" -ref "${run_02}.rst7" >> mdinfo
$pmemd_path -O -i "${run_04}.in" -o "${run_04}.out" -p "$prmtop_file" -c "${run_03}.rst7" -r "${run_04}.rst7" -ref "${run_03}.rst7" >> mdinfo
$sander_path -O -i "${run_05}.in" -o "${run_05}.out" -p "$prmtop_file" -c "${run_04}.rst7" -r "${run_05}.rst7" -ref "${run_04}.rst7" >> mdinfo
$pmemd_path -O -i "${run_06}.in" -o "${run_06}.out" -p "$prmtop_file" -c "${run_05}.rst7" -r "${run_06}.rst7" -ref "${run_05}.rst7" >> mdinfo
$pmemd_path -O -i "${run_07}.in" -o "${run_07}.out" -p "$prmtop_file" -c "${run_06}.rst7" -r "${run_07}.rst7" -ref "${run_06}.rst7" >> mdinfo
$pmemd_path -O -i "${run_08}.in" -o "${run_08}.out" -p "$prmtop_file" -c "${run_07}.rst7" -r "${run_08}.rst7" -ref "${run_07}.rst7" >> mdinfo
$pmemd_path -O -i "${run_09}.in" -o "${run_09}.out" -p "$prmtop_file" -c "${run_08}.rst7" -r "${run_09}.rst7" >> mdinfo

$pmemd_path -O -i "${prod}.in" -o "${prod}_1.out" -p "$prmtop_file" -c "${run_09}.rst7" -r "${prod}_1.rst7" -x "${prod}_1.mdcrd" >> mdinfo
$cpptraj_path -p "$prmtop_file" < "$center_file"

for i in {1..4}; do
    j=$((i + 1))
    $pmemd_path -O -i "${prod}.in" -o "${prod}_$j.out" -p "$prmtop_file" -c "${prod}_$i.rst7" -r "${prod}_$j.rst7" -x "${prod}_$j.mdcrd" >> mdinfo 
done

module purge
module load amber
for j in {1..5}; do
cat > "getpdb.in" <<EOF
parm KLVFFAEini_12x2_H2O.prmtop
trajin "prod_${j}.rst7"
center :1-168 mass origin
image origin center familiar
strip :WAT
trajout "pdbfiles/pep_$((5 * j))ns.pdb" pdb
run

EOF
cpptraj -i "getpdb.in" > pdb.out

done