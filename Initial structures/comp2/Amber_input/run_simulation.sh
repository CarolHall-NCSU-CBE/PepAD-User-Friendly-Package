#!/bin/bash
#BSUB -n 1
#BSUB -W 72:00
#BSUB -J amber
#BSUB -R "select[a100 || a30 || a10 || l40 || rtx2080 || gtx1080]"
#BSUB -gpu "num=1:mode=shared:mps=yes"
#BSUB -q gpu
#BSUB -o out.%J
#BSUB -e err.%J
################## need to select gpu model to prevent cuda error ####

pmemd_path=/rs1/researchers/h/hall/amber24/amber24_src/build/src/pmemd/src/pmemd.cuda_SPFP
sander_path=/rs1/researchers/h/hall/amber24/amber24_src/build/AmberTools/src/sander/sander
cpptraj_path=/rs1/researchers/h/hall/amber24/amber24_src/build/AmberTools/src/cpptraj/src/cpptraj
#tleap_path=/rs1/researchers/h/hall/amber24/amber24_src/AmberTools/src/leap/tleap
max=4
module load amber
tleap_file="tleap.in"
prmtop_file_1="top_crd/sheets_im_1.prmtop"
prmtop_file_2="top_crd/sheets_im_2.prmtop"
prmtop_file_3="top_crd/sheets_h2o.prmtop"
init_inpcrd_file_1="top_crd/sheets_im.inpcrd"
init_inpcrd_file_2="top_crd/sheets_h2o.inpcrd"
run_01="min_01"
run_02="heat_02"
run_03="equil_03"
run_04="equil_04"
run_05="equil_05"
run_06="min_06"
run_07="heat_07"
run_08="equil_08"
run_09="equil_09"
run_10="min_10"
run_11="equil_11"
run_12="equil_12"
run_13="equil_13"
run_14="equil_14"
prod="prod"
center_file="center.in"
getpdb="getpdb.in"
pep=E1209_11A

get_largest_prod_number() {
    # Find all files that match the pattern 'prodXX', extract the number, and sort it
    largest=$(ls output/prod[0-9][0-9].out 2>/dev/null | sed -E 's/output\/prod([0-9]{2}).out/\1/' | sort -n | tail -n 1)
    
    # If no files exist, return 0
    if [ -z "$largest" ]; then
        echo 0
    else
        echo $largest
    fi
}

k=$(get_largest_prod_number)
k=$((10#$k))  # Convert the string to a number

########################### Judge the output file ######################
recent_prod=$(printf "output/prod%02d.out" "$k")
echo "Recent prod file: $recent_prod"

if grep -q "5.  TIMINGS" "$recent_prod"; then
    echo "Found '5. TIMINGS' in $recent_prod. Proceeding to the next step."
else
    echo "'5. TIMINGS' not found in $recent_prod. Running the simulation at step $k."
    k=$((k - 1))
fi


########################### run simulation ############################
if [ "$k" -lt 1 ]; then
    echo "i is less than 2"

############# Use implicit water ###################
cat > "$tleap_file" <<EOF
    source leaprc.protein.ff14SBonlysc
    pep = loadpdb "${pep}.pdb"
    set default pbradii mbondi3
    saveamberparm pep "${prmtop_file_1}" "${init_inpcrd_file_1}"

EOF
    tleap -f "$tleap_file"
    rm -f "$tleap_file"

cat > parmed.in <<EOF
    hmassrepartition
    outparm "${prmtop_file_2}"
    quit

EOF
    parmed -i parmed.in -p "${prmtop_file_1}"

############ run equilibration with implicit water #################
    module purge
    module load openmpi-gcc/openmpi5.0.3-gcc11.3.1 cuda/10.0 gcc/13
    $sander_path -O -i "input/${run_01}.in" -o "output/${run_01}.out" -p "$prmtop_file_2" -c "$init_inpcrd_file_1" -r "rst7/${run_01}.rst7" >> mdinfo                           # minimize with restraints
    $pmemd_path -O -i "input/${run_02}.in" -o "output/${run_02}.out" -p "$prmtop_file_2" -c "rst7/${run_01}.rst7" -r "rst7/${run_02}.rst7" -ref "rst7/${run_01}.rst7" >> mdinfo # heat with restraints
    $pmemd_path -O -i "input/${run_03}.in" -o "output/${run_03}.out" -p "$prmtop_file_2" -c "rst7/${run_02}.rst7" -r "rst7/${run_03}.rst7" -ref "rst7/${run_02}.rst7" >> mdinfo # equilibrate with restraints
    $pmemd_path -O -i "input/${run_04}.in" -o "output/${run_04}.out" -p "$prmtop_file_2" -c "rst7/${run_03}.rst7" -r "rst7/${run_04}.rst7" -ref "rst7/${run_03}.rst7" >> mdinfo # equilibrate with restraints
    $pmemd_path -O -i "input/${run_05}.in" -o "output/${run_05}.out" -p "$prmtop_file_2" -c "rst7/${run_04}.rst7" -r "rst7/${run_05}.rst7" >> mdinfo                            # equilibrate 

############# Use explicit water ###################
    module purge
    module load amber

cat > getpdb.in <<EOF
    parm top_crd/sheets_im_2.prmtop
    trajin rst7/equil_05.rst7 
    center :1-384 mass origin
    image origin center familiar
    trajout "pdbfiles/${pep}_equil_05.pdb" pdb
    run

EOF
    cpptraj -i getpdb.in


cat > "$tleap_file" <<EOF
    source leaprc.protein.ff14SB
    pep = loadpdb "pdbfiles/${pep}_equil_05.pdb"
    source leaprc.water.tip3p
    solvateoct pep TIP3PBOX 12.0
    saveamberparm pep "${prmtop_file_3}" "${init_inpcrd_file_2}"

EOF
    tleap -f "$tleap_file"
    rm -f "$tleap_file"

############ run equilibration with explicit water #################
    module purge
    module load openmpi-gcc/openmpi5.0.3-gcc11.3.1 cuda/10.0 gcc/13
    $sander_path -O -i "input/${run_06}.in" -o "output/${run_06}.out" -p "$prmtop_file_3" -c "$init_inpcrd_file_2" -r "rst7/${run_06}.rst7" -ref "$init_inpcrd_file_2" >> mdinfo
    $pmemd_path -O -i "input/${run_07}.in" -o "output/${run_07}.out" -p "$prmtop_file_3" -c "rst7/${run_06}.rst7" -r "rst7/${run_07}.rst7" -ref "rst7/${run_06}.rst7" >> mdinfo
    $pmemd_path -O -i "input/${run_08}.in" -o "output/${run_08}.out" -p "$prmtop_file_3" -c "rst7/${run_07}.rst7" -r "rst7/${run_08}.rst7" -ref "rst7/${run_07}.rst7" >> mdinfo
    $pmemd_path -O -i "input/${run_09}.in" -o "output/${run_09}.out" -p "$prmtop_file_3" -c "rst7/${run_08}.rst7" -r "rst7/${run_09}.rst7" -ref "rst7/${run_08}.rst7" >> mdinfo
    $sander_path -O -i "input/${run_10}.in" -o "output/${run_10}.out" -p "$prmtop_file_3" -c "rst7/${run_09}.rst7" -r "rst7/${run_10}.rst7" -ref "rst7/${run_09}.rst7" >> mdinfo
    $pmemd_path -O -i "input/${run_11}.in" -o "output/${run_11}.out" -p "$prmtop_file_3" -c "rst7/${run_10}.rst7" -r "rst7/${run_11}.rst7" -ref "rst7/${run_10}.rst7" >> mdinfo
    $pmemd_path -O -i "input/${run_12}.in" -o "output/${run_12}.out" -p "$prmtop_file_3" -c "rst7/${run_11}.rst7" -r "rst7/${run_12}.rst7" -ref "rst7/${run_11}.rst7" >> mdinfo
    $pmemd_path -O -i "input/${run_13}.in" -o "output/${run_13}.out" -p "$prmtop_file_3" -c "rst7/${run_12}.rst7" -r "rst7/${run_13}.rst7" -ref "rst7/${run_12}.rst7" >> mdinfo
    $pmemd_path -O -i "input/${run_14}.in" -o "output/${run_14}.out" -p "$prmtop_file_3" -c "rst7/${run_13}.rst7" -r "rst7/${run_14}.rst7" >> mdinfo
    $pmemd_path -O -i "input/${prod}.in" -o "output/${prod}01.out" -p "$prmtop_file_3" -c "rst7/${run_09}.rst7" -r "rst7/${prod}_01.rst7" -x "coordinate/${prod}_01.mdcrd" >> mdinfo

    for i in $(seq 1 $max); do
        j=$((i + 1))

        if (( j < 10 )); then
            $pmemd_path -O -i "input/${prod}.in" -o "output/${prod}0${j}.out" -p "$prmtop_file_3" -c "rst7/${prod}_0${i}.rst7" -r "rst7/${prod}_0${j}.rst7" -x "coordinate/${prod}_0${j}.mdcrd" >> mdinfo 
        elif (( j == 10 )); then
            $pmemd_path -O -i "input/${prod}.in" -o "output/${prod}${j}.out" -p "$prmtop_file_3" -c "rst7/${prod}_0${i}.rst7" -r "rst7/${prod}_${j}.rst7" -x "coordinate/${prod}_${j}.mdcrd" >> mdinfo 
        else
            $pmemd_path -O -i "input/${prod}.in" -o "output/${prod}${j}.out" -p "$prmtop_file_3" -c "rst7/${prod}_${i}.rst7" -r "rst7/${prod}_${j}.rst7" -x "coordinate/${prod}_${j}.mdcrd" >> mdinfo 
        fi
    done

elif [ "$k" -le $max ]; then
    module purge
    module load openmpi-gcc/openmpi5.0.3-gcc11.3.1 cuda/10.0 gcc/13
    for i in $(seq $k $max); do
        j=$((i + 1))

        if (( j < 10 )); then
            $pmemd_path -O -i "input/${prod}.in" -o "output/${prod}0${j}.out" -p "$prmtop_file_3" -c "rst7/${prod}_0${i}.rst7" -r "rst7/${prod}_0${j}.rst7" -x "coordinate/${prod}_0${j}.mdcrd" >> mdinfo 
        elif (( j == 10 )); then
            $pmemd_path -O -i "input/${prod}.in" -o "output/${prod}${j}.out" -p "$prmtop_file_3" -c "rst7/${prod}_0${i}.rst7" -r "rst7/${prod}_${j}.rst7" -x "coordinate/${prod}_${j}.mdcrd" >> mdinfo 
        else
            $pmemd_path -O -i "input/${prod}.in" -o "output/${prod}${j}.out" -p "$prmtop_file_3" -c "rst7/${prod}_${i}.rst7" -r "rst7/${prod}_${j}.rst7" -x "coordinate/${prod}_${j}.mdcrd" >> mdinfo 
        fi
    done
fi