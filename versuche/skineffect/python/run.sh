#!/usr/bin/env sh


# ------------------------------------------------------------------ #
# This script executes all python scripts in this directory          #
# ------------------------------------------------------------------ #

# -------------------------------------------------------- #
# First,  check if  Python3 is  actually installed  (there #
# might be more dependencies, obviosuly).                  #
# -------------------------------------------------------- #
python3 -c 'print("Running everything!")' || exit 1

for file in hohl*py;do
    echo "Running ${file}"
    python3 "${file}"
done
for file in voll*py;do
    echo "Running ${file}"
    python3 "${file}"
done
for file in stue*py;do
    echo "Running ${file}"
    python3 "${file}"
done

echo "Running Error Analysis"
python3 error_analysis.py
