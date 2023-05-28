#!/bin/bash

rp_newton="$1"
loop_number="$2"
echo "loop number $loop_number"
echo "Running geodesic simulation for rp_newton=$rp_newton, loop=$loop_number"

# Generate the orbit.params file using Python script
# Path to the generate_orbit_params.py script
generate_script="./utils/write_orbit_params.py"
output_geodesics='./geodesic_outputs'

python3 "${generate_script}" ${rp_newton} 0.0

# Set the positions and velocities filenames
export positions_file="${output_geodesics}/positions_${loop_number}.dat"
export velocities_file="${output_geodesics}/velocities_${loop_number}.dat"
# Run grtest with orbit.params file and pass the filenames as arguments
./grtest < orbit.params
./grtest output_00000.dat

echo "Geodesic simulation complete for rp_newton=$rp_newton, loop=$loop_number"

