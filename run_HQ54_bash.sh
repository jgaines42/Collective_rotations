#!/bin/bash 
# Takes 5 parameters
# 1: pdb_name: 4 letter code for the PDB
# 2: babl: integer indicating which ba/bl variant to use
# 3: group: integer of which group/cluster in the PDB to use
# 4: folder_data: folder where all data will be stored
# 
pdb_name=$1
babl=$2
group=$3
folder_data=$4
folder_pdb=$5
echo $pdb_name
initial_E=1
finished=0
start_E=100

#Load matlab module
module load Apps/Matlab/R2014a


while [  $finished -le 0 ]
	do
		finished=1
		matlab -nosplash -nodisplay -nojvm -r "find_minimum_energy_toC_bash('$pdb_name' ,$babl, $group,$initial_E,'$folder_data','$folder_pdb'); quit" 

		result=`./find_Min $pdb_name $babl $group $folder_data $folder_pdb`
		echo $result
		if [ $result -eq "0" ];
			 then
			 	echo $initial_E
			 	initial_E=$(expr $initial_E + 1)
			 	finished=0
			 	echo 6
			 else
			 	finished=1
			 	echo 7
		fi
		echo 5
		
	done

matlab -nosplash -nodisplay -nojvm -r "process_C_recursive('$pdb_name' ,$babl, '$database', $group); quit" 

