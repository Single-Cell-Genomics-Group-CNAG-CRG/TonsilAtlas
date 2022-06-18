for i in "0"  "1"  "2"  "3"  "4"  "5"  "6"  "7"  "8"  "9"  "10" "11" "12";
do
	sbatch  \
	    --reservation=massoni-reserv \
	    --err "log/markers/markers_k${i}.err" \
	    --out "log/markers/markers_k${i}.log" \
	    -J "markers_k${i}" \
	    -c 8 \
	    --time 23:00:00 \
	    --wrap "module purge; module load gcc/6.3.0 hdf5/1.10.1 R/3.6.0; echo [`date "+%Y-%m-%d %T"`] starting job on $HOSTNAME; Rscript ./01-find_markers.R ${i}; echo [`date "+%Y-%m-%d %T"`] job finished"
done
