#!/bin/bash

# Array of rp_newton values
# Start and end values (logarithmic scale)
start=20.
end=80.

# Number of grid points
grid_points=100

# Calculate the logarithmic base
log_base=$(awk -v start="$start" -v end="$end" -v grid_points="$grid_points" 'BEGIN { print exp(log(end/start) / (grid_points-1)) }')

# Initialize the rp_newton_values array
rp_newton_values=()

# Loop to generate logarithmically spaced values
for ((i=0; i<grid_points; i++)); do
    value=$(awk -v start="$start" -v i="$i" -v log_base="$log_base" 'BEGIN { print start * (log_base^i) }')
    rp_newton_values+=($value)
done

# Print the rp_newton_valuesarray
for value in "${rp_newton_values[@]}"; do
    echo "$value"
done


# Array of inclination values
inc_parabola_values=(0.0)

# Path to the generate_orbit_params.py script
generate_script="/Users/megha/geodesic_inputs/write_orbit_params.py"

# Path to the new folder where position and velocity files need to be saved
output_geodesics='/Users/megha/geodesic_outputs'

# Remove the old rp_position_mapping.txt file, if it exists
rm -f "${output_geodesics}/rp_position_mapping.txt"

# Write the file headers
echo "rp_newton, positions_file" >> "${output_geodesics}/rp_position_mapping.txt"

# Loop through the arrays
for ((i=0; i<${#rp_newton_values[@]}; i++)); do
    rp_newton=${rp_newton_values[i]}
    inc_parabola=${inc_parabola_values[0]}

    # Generate orbit.params file using Python script
    python3.10 "${generate_script}" ${rp_newton} ${inc_parabola}

    echo " *******************************************************************"
    echo "Newtonian pericentre distance: ${rp_newton} code units"

    # Run grtest with orbit.params file
    ./grtest < orbit.params
    ./grtest output_00000.dat 

    # Move positions and velocities.dat files to output_geodesics
    mv positions.dat "${output_geodesics}/positions_${i}.dat"
    mv velocities.dat "${output_geodesics}/velocities_${i}.dat"

    # Write rp_newton and filename to a file
    echo "${rp_newton}, positions_${i}.dat" >> "${output_geodesics}/rp_position_mapping.txt"
    echo " *******************************************************************"

done




