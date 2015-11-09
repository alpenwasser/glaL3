#!/usr/bin/env sh


# ------------------------------------------------------------------ #
# This script executes all python scripts in this directory          #
# ------------------------------------------------------------------ #

# -------------------------------------------------------- #
# First,  check if  Python3 is  actually installed  (there #
# might be more dependencies, obviosuly).                  #
# -------------------------------------------------------- #
python3 -c 'print("Running everything!")' || exit 1

for file in *py;do
    echo "Running ${file}"
    python3 "${file}"
done
