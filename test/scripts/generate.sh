#! /bin/bash

FAMSA="../../bin/famsa"
N_THREADS=64

###############################################################################
# adeno-fiber - full tree
#
INPUT="../adeno_fiber/adeno_fiber"
REF_DIR="../adeno_fiber"

trees=("sl" "upgma" "slink") 
for tree in "${trees[@]}"; do
  set -x
  ${FAMSA} -gt ${tree} -gt_export ${INPUT} "${REF_DIR}/${tree}.dnd" 
  ${FAMSA} -gt ${tree} ${INPUT} "${REF_DIR}/${tree}.fasta" -t ${N_THREADS}
  set +x
done

###############################################################################
# adeno-fiber - duplicates
#
INPUT="../adeno_fiber_duplicates/adeno_fiber_duplicates"
REF_DIR="../adeno_fiber_duplicates"
tree="sl"

set -x
${FAMSA} -gt ${tree} -gt_export ${INPUT} "${REF_DIR}/${tree}.dnd"
${FAMSA} -gt ${tree} ${INPUT} "${REF_DIR}/${tree}.fasta" -t ${N_THREADS}
set +x

###############################################################################
# hemopexin - medoid tree
#
REF_DIR="../hemopexin"
INPUT="../hemopexin/hemopexin"
trees=("sl" "upgma" "slink" "nj") 

for tree in "${trees[@]}"; do
  set -x
  ${FAMSA}  -medoidtree -gt ${tree} -gt_export ${INPUT} "${REF_DIR}/medoid-${tree}.dnd"
  ${FAMSA}  -medoidtree -gt ${tree} ${INPUT} "${REF_DIR}/medoid-${tree}.fasta" -t ${N_THREADS}
  ${FAMSA}  -medoidtree -gt ${tree} -gt_export -subtree_size 10 -sample_size 100 -cluster_fraction 0.2 -cluster_iters 1 ${INPUT} "${REF_DIR}/medoid-${tree}-params.dnd"   
  set +x
done

###############################################################################
# hemopexin - duplicates
#
REF_DIR="../hemopexin_duplicates"
INPUT="../hemopexin_duplicates/hemopexin_duplicates"
tree="sl" 

set -x
${FAMSA}  -medoidtree -gt ${tree} -gt_export ${INPUT} "${REF_DIR}/medoid-${tree}.dnd"
${FAMSA}  -medoidtree -gt ${tree} ${INPUT} "${REF_DIR}/medoid-${tree}.fasta" -t ${N_THREADS}

${FAMSA} -keep-duplicates -medoidtree -gt ${tree} -gt_export ${INPUT} "${REF_DIR}/medoid-${tree}-dups.dnd"
${FAMSA} -keep-duplicates -medoidtree -gt ${tree} ${INPUT} "${REF_DIR}/medoid-${tree}-dups.fasta" -t ${N_THREADS}

set +x


###############################################################################
# adeno-fiber - other tests
#
INPUT="../adeno_fiber/adeno_fiber"
REF_DIR="../adeno_fiber"
tree="sl"

set -x
${FAMSA} -go 10 -ge 2 -tgo 0.5 -tge 1.0 -gsd 3 -gsl 30 ${INPUT} "${REF_DIR}/gaps.fasta" -t ${N_THREADS}
${FAMSA} -gt import "${REF_DIR}/upgma.dnd" -refine_mode off ${INPUT} "${REF_DIR}/upgma.no_refine.fasta" -t ${N_THREADS}
${FAMSA} -refine_mode off "${REF_DIR}/upgma.no_refine.part1.fasta" "${REF_DIR}/upgma.no_refine.part2.fasta" "${REF_DIR}/upgma.pp.fasta" -t ${N_THREADS}
   
     
${FAMSA} -dist_export ${INPUT} "${REF_DIR}/dist.csv"
${FAMSA} -dist_export -pid ${INPUT} "${REF_DIR}/pid.csv"
${FAMSA} -dist_export -square_matrix ${INPUT} "${REF_DIR}/dist_sq.csv"
${FAMSA} -dist_export -square_matrix -pid ${INPUT} "${REF_DIR}/pid_sq.csv"
        
${FAMSA} ../adeno_fiber_extra/adeno_fiber_extra ../adeno_fiber_extra/ref.fasta -t ${N_THREADS}
set +x  
  
###############################################################################
# LRR - 100k trees
#
REF_DIR="../LRR"
INPUT="../LRR/LRR"
tree="sl"

set -x
${FAMSA} -gt ${tree} -gt_export ${INPUT} "${REF_DIR}/${tree}.dnd"
set +x

