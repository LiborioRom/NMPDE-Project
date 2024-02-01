#!/bin/bash


DEGREE=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10")
PRECONDITIONERS=("ssor" "amg" "identity" "jacobi") # add to the vector new preconditioners
MESH="mesh-cube-20.msh"
P_VALUE=4
s="0"
PROCESSORS=("1" "2" "4" "8" "16" "32")
a="2"

for preco in "${PRECONDITIONERS[@]}"
do
  for d in "${DEGREE[@]}"
  do
    for ((c=1; c<=1; c++))
    do
      mpiexec -n 4 ../build/Project_parallel -m "$MESH" -P "$preco" -p "$P_VALUE" -r "$d" -s "$s"
    done
  done


for preco in "${PRECONDITIONERS[@]}"
do
  for n in "${PROCESSORS[@]}"
  do
    for ((c=1; c<=1; c++))
    do
      mpiexec -n n ../build/Project_parallel -m "$MESH" -P "$preco" -p "$P_VALUE" -r "$a" -s "$s"
    done
  done
done


# python3 cond_number_vs_p_value.py
# python3 iterations_vs_p_value.py