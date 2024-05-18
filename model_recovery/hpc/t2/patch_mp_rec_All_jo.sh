#!/bin/bash
# source fc.sh # load myMod and myCeil functions

### This is for selective iterations that did not run
mp_rec_outputs/Fm1_Wm3_All_N30_mp_rec_3031222_29.out
mp_rec_outputs/Fm3_Wm1_All_N30_mp_rec_2933850_23.out
mp_rec_outputs/Fm3_Wm1_All_N30_mp_rec_3031197_36.out
mp_rec_outputs/Fm3_Wm2_All_N30_mp_rec_3031199_28.out
mp_rec_outputs/Fm4_Wm1_All_N30_mp_rec_3031204_50.out
mp_rec_outputs/Fm5_Wm2_All_N30_mp_rec_3031216_29.out
array_def=(m1 m3 29);
array_def=(m3 m1 23);
array_def=(m3 m2 28);
array_def=(m4 m1 50);
array_def=(m5 m2 29);

#from, with, iter



Nsubj=30;
g=All;
stanSpec=mp_rec_spec;
array_def=(m5 m2 29);
ls ../results/model_rec/${array_def[0]}/${array_def[1]}/*s${array_def[2]}*.csv;
# rm -r ../results/model_rec/${array_def[0]}/${array_def[1]}/*s${array_def[2]}*.csv;
if [ ${array_def[0]} == ${array_def[1]} ];
then
	ls ../results/param_rec/${array_def[0]}/*s${array_def[2]}*.csv;
	ls ../results/param_rec_extra/${array_def[0]}_extra/s${array_def[2]};
	# rm -r ../results/param_rec/${array_def[0]}/*s${array_def[2]}*.csv;
	# rm -r ../results/param_rec_extra/${array_def[0]}_extra/s${array_def[2]};
fi;

sbatch --array=${array_def[2]} --job-name=F${array_def[0]}_W${array_def[1]}_${g}_N$Nsubj t2/slrm_mp_recovery.slurm ${array_def[0]} ${array_def[1]} $g $Nsubj $stanSpec;

Nsubj=30;
g=All;
stanSpec=mp_rec_spec;
array_def=(m1 m3 29);
sbatch --array=${array_def[2]} --job-name=F${array_def[0]}_W${array_def[1]}_${g}_N$Nsubj t2/slrm_mp_recovery.slurm ${array_def[0]} ${array_def[1]} $g $Nsubj $stanSpec;
array_def=(m3 m1 23);
sbatch --array=${array_def[2]} --job-name=F${array_def[0]}_W${array_def[1]}_${g}_N$Nsubj t2/slrm_mp_recovery.slurm ${array_def[0]} ${array_def[1]} $g $Nsubj $stanSpec;
array_def=(m3 m2 28);
sbatch --array=${array_def[2]} --job-name=F${array_def[0]}_W${array_def[1]}_${g}_N$Nsubj t2/slrm_mp_recovery.slurm ${array_def[0]} ${array_def[1]} $g $Nsubj $stanSpec;
array_def=(m4 m1 50);
sbatch --array=${array_def[2]} --job-name=F${array_def[0]}_W${array_def[1]}_${g}_N$Nsubj t2/slrm_mp_recovery.slurm ${array_def[0]} ${array_def[1]} $g $Nsubj $stanSpec;
array_def=(m5 m2 29);
sbatch --array=${array_def[2]} --job-name=F${array_def[0]}_W${array_def[1]}_${g}_N$Nsubj t2/slrm_mp_recovery.slurm ${array_def[0]} ${array_def[1]} $g $Nsubj $stanSpec;

Nsubj=30;
g=All;
stanSpec=mp_rec_spec;
sbatch --array=36 --job-name=Fm3_Wm1_${g}_N$Nsubj t2/slrm_mp_recovery.slurm m3 m1 $g $Nsubj $stanSpec;


# sbatch --array=${sIterFrom}-${sIterTo} --job-name=F${model_from}_W${model_with}_${g}_N$Nsubj t2/slrm_mp_recovery.slurm $model_from $model_with $g $Nsubj $stanSpec
#
#
# m_from=(m1 m2 m3 m4 m5) # specify model names that you generate data from - has to be of the m1, m2,... form
# m_with=(m1 m2 m3 m4 m5) # specify model names that you fit data the generated data with - has to be of the m1, m2,... form
# # m_from=(m1 m2) # specify model names that you generate data from - has to be of the m1, m2,... form
# # m_with=(m1 m2) # specify model names that you fit data the generated data with - has to be of the m1, m2,... form
# Nsubj=30 #specify how many subjects you want to simulate data for
# sIterFrom=26 #specify the first of the number of elements in the job array on hpc
# sIterTo=50 #specify the last of the number of elements in the job array on hpc
# modelNo=${#m_from[@]}
# elementsNo=$(($modelNo*$modelNo))
# iTo_pr=$((${modelNo}-1)) # specify how many diagonal elements there are and subtract 1  (Dn -1) where Dn is the number of diagonal elements in the confusion matrix
# # iTo_mr=5 # specify how many off-diagonal elemets are and subtract 1  (oDn -1) where oDn is the number of off-diagonal elements in the confusion matrix
# iTo_mr=$(( ${elementsNo} -$modelNo - 1)) # specify how many off-diagonal elemets are and subtract 1  (oDn -1) where oDn is the number of off-diagonal elements in the confusion matrix
#
#
# a=($(seq 1 ${elementsNo})); #all model combs - specify the sequence indexing elements in the confusion matrix (flattened) - number goes from top left to right from the top to bottom
# b=($(seq 1 $((${modelNo}+1)) ${elementsNo})) #model combination for p.recoveyr - specify the indices of the diagnoal elements
# # b=(1 5 9); #model combination for p.recoveyr - specify the indices of the diagnoal elements
# c=(`echo ${a[@]} ${b[@]} | tr ' ' '\n' | sort | uniq -u `) #model combs for p.rec only
#
# ### optionally exclude some models
# # exc_m2=(2 7 12 17 22) # optional exclusion of a particular model that you fit data with
# # c=(`echo ${c[@]} ${exc_m2[@]} | tr ' ' '\n' | sort | uniq -u `) #all but excluded model fit combinations
#
# stanSpec=mp_rec_spec; # load stan spec
#
# echo --------------------Parameter recovery--------------------
# for i in $(seq 0 $iTo_pr);
# do
# 	g=All
# 	echo -----$g + $Nsubj-----
# 	m=${b[i]}; d=${modelNo}; # d - specify the number of models
# 	my_ceil; # ceil function which gets which model to generate data from <=> ceil(b/d)
# 	model_from=${m_from[(($?-1))]}
# 	my_mod; # mod function which gets which model to fit data with
# 	model_with=${m_with[(($?-1))]}
# 	echo $model_from
# 	echo $model_with
# 	sbatch --array=${sIterFrom}-${sIterTo} --job-name=F${model_from}_W${model_with}_${g}_N$Nsubj t2/slrm_mp_recovery.slurm $model_from $model_with $g $Nsubj $stanSpec
#
# 	# Rscript ../mp_recovery.R $model_from $model_with $g $Nsubj $stanSpec 1
# 	# $model_from $model_with $g $Nsubj
# 	# echo --
# done;
#
# echo --------------------Model recovery--------------------
# for i in $(seq 0 $iTo_mr);
# do
# 	g=All
# 	echo -----$g + $Nsubj-----
# 	m=${c[i]}; d=${modelNo}; # d - specify the number of models
# 	my_ceil; # ceil function which gets which model to generate data from <=> ceil(b/d)
# 	model_from=${m_from[(($?-1))]}
# 	my_mod; # mod function which gets which model to fit data with
# 	model_with=${m_with[(($?-1))]}
# 	echo $model_from
# 	echo $model_with
# 	sbatch --array=${sIterFrom}-${sIterTo} --job-name=F${model_from}_W${model_with}_${g}_N$Nsubj t2/slrm_mp_recovery.slurm $model_from $model_with $g $Nsubj $stanSpec
#
# 	# Rscript ../mp_recovery.R $model_from $model_with $g $Nsubj $stanSpec 1
# 	# $model_from $model_with $g $Nsubj
# 	# echo --
# done;
# unset my_mod
# unset my_ceil
