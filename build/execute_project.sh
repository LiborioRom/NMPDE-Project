#!/bin/bash




P_VALUE="6" # modify here the value of p
PRECONDITIONERS=("identity" "jacobi" "ssor" "sor") # add to the vector new preconditioners
MESH="mesh-cube-20.msh"




for PRECONDITIONER in "${PRECONDITIONERS[@]}"
do
for ((c=1; c<=10; c++))
do
  ./Project "$MESH" "$PRECONDITIONER" "$P_VALUE"
done
done