#! /bin/bash

methods=('grad' 'fpt1' 'fpt2')
Kvals=(2 4 8 16)
pvals=(2)
nvals=(1000)

nrep=10
niter=5

for ((i_m=0; i_m<${#methods[@]} ;i_m++))
do
for ((i_K=0; i_K<${#Kvals[@]} ;i_K++))
do
for ((i_p=0; i_p<${#pvals[@]} ;i_p++))
do
for ((i_n=0; i_n<${#nvals[@]} ;i_n++))
do
  method=${methods[$i_m]}
  K=${Kvals[i_K]}
  n=${nvals[$i_p]}
  p=${pvals[$i_n]}
  Rscript lm-forest-tree-bench.R $method $K $p $n $nrep $niter
done
done
done
done
