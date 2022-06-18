PATH_OBJECTS=../results/Visualization/

for chr in ${PATH_OBJECTS}*.bedGraph 
do
LC_COLLATE=C sort -k1,1 -k2,2n $chr > tmp.bedGraph
bedgraphtobigwig tmp.bedGraph chrom.sizes $PATH_OBJECTS$(basename ${chr/bedGraph/bw})
done 
rm tmp.bedGraph
