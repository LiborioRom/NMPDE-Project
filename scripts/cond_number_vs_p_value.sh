#!/bin/bash


PRECONDITIONERS=("identity" "jacobi" "ssor" "amg") # add to the vector new preconditioners
MESH="mesh-cube-20.msh"

for PRECONDITIONER in "${PRECONDITIONERS[@]}"
do
  for P_VALUE in {0..8}
  do
    for ((c=1; c<=1; c++))
    do
      pwd
      ../build/Project_parallel -m "$MESH"  -P "$PRECONDITIONER" -p "$P_VALUE"
    done
  done
done

python3 cond_number_vs_p_value.py