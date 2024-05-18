stanSpec=spec_test;
g=All;
i=1;
for m in specs/models/m4*_*sliced*.txt;
do
  n=${m##*/}; p=${n%.*};
  echo $i $p
  for c in $(seq 1 1 4)
  do
    echo $c;
    Rscript ../fit_models_cs.R ${p} ${g} ${stanSpec} ${c}
  done


done;
