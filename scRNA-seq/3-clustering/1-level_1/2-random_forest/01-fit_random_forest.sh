declare -a resolutions=("0.01" "0.03" "0.05" "0.07" "0.09" "0.11" "0.13" "0.15" "0.17" "0.19" "0.21" "0.23" "0.25" "0.27" "0.29" "0.31" "0.33" "0.35" "0.37" "0.39" "0.41" "0.43" "0.45" "0.47" "0.49" "0.51" "0.53" "0.55" "0.57" "0.59" "0.61" "0.63" "0.65" "0.67" "0.69")

for res in ${resolutions[@]}; do
   sbatch --reservation=massoni-reserv --err "log/rf/rf_${res}.err" --out "log/rf/rf_${res}.log" -J "rf_${res}" -c 4 --time 15:00:00 --wrap "module purge; module load gcc/6.3.0 hdf5/1.10.1 R/3.6.0; echo [`date "+%Y-%m-%d %T"`] starting job on $HOSTNAME; Rscript ./01-fit_random_forest.R ${res}; echo [`date "+%Y-%m-%d %T"`] job finished" 
done
