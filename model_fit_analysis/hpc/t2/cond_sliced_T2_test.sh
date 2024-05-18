stanSpec=spec_test;
g=All;
for m in specs/models/m*_*sliced*norm.txt;
do
  n=${m##*/}; p=${n%.*};
  sbatch --array=1-4 --job-name=${p}_${g} t2/slrm_cond_sliced_T2_test.slurm ${p} ${g} ${stanSpec} $1;
  # Rscript ../prefit_models.R ${p} ${g} ${stanSpec};

done;
