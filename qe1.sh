#!/bin/bash
nconf=$(grep -o 'Lattice=' db_W.xyz | wc -l)
#nconf=16
cd local/
run_and_echo() {
timeout 2s pw.x < $1.in > $1.out ; echo "File $1.out generated"
}
export -f run_and_echo
seq 1 $nconf | parallel --will-cite --keep-order -j8 run_and_echo {}
cd ..
