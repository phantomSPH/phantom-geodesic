import sys

def generate_orbit_params(rp_newton, inc_parabola):
	output_file = '../orbit.params'
	with open(output_file, "w") as file:
		file.write("# tde setup file\n")
		file.write("rp_newton =      {:.3f}    ! newtonian rp\n".format(rp_newton))
		file.write("inc_parabola =         {:.3f}    ! inc of orbit and then update the numbers\n".format(inc_parabola))
	print("Finished generating a new orbits.params file with rp and inc as ",rp_newton,inc_parabola)
if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python generate_orbit_params.py <rp_newton> <inc_parabola>")
        sys.exit(1)

    rp_newton = float(sys.argv[1])
    inc_parabola = float(sys.argv[2])
    print(rp_newton,"rp newtonian in python script")
    generate_orbit_params(rp_newton, inc_parabola)

