import copy
import os
import time
from math import ceil, radians, cos, sin, sqrt
from file_conversion import convert_cssr_to_cif
"""
                           find_interp_structures.py

 This program generates potential interpenetrated structures for all CSSR files
 in the working directory. Attempts at interpenetrated structures are made
 using the displacement set by SHIFT_DISTANCE, and then the structure is
 checked for overlapping atoms using the hard-sphere radius set by
 OVERLAP_RADIUS. Set SAVE_TRIAL_STRUCTURES to True to save trial interpenetrated
 structures (SAVE_TRIAL_STRUCTURES == True recommended for debugging). Set
 FRACTIONAL_COORDS to True if CSSR input is in fractional coordinates, set to
 False if in Cartesian. Set CSSR_TO_CIF to True to convert the interpenetrated
 structures to CIF file format from the default CSSR.

"""
############################## SET VARIABLES HERE ##############################
SHIFT_DISTANCE = 1.0
OVERLAP_RADIUS = 1.53 # chosen because C covalent radius set to 0.76
SAVE_TRIAL_STRUCTURES = False
FRACTIONAL_COORDS = True
CSSR_TO_CIF = True
################################ END VARIABLES #################################

def translate_coords(x, y, z, old_list):
    """
    Translate atoms by set distance, wrapping them around if they go past
    the unit cell boundaries.
    """
    new_list = copy.deepcopy(old_list)

    for line in new_list:
        if float(line[2]) + x > 1:
            line[2] = str(float(line[2]) + x - 1)
        else:
            line[2] = str(float(line[2]) + x)

        if float(line[3]) + y > 1:
            line[3] = str(float(line[3]) + y - 1)
        else:
            line[3] = str(float(line[3]) + y)

        if float(line[4]) + z > 1:
            line[4] = str(float(line[4]) + z - 1)
        else:
            line[4] = str(float(line[4]) + z)
    return new_list

def get_boundary_atoms(all_atoms, a, b, c):
    boundary_atoms = []

    for line in all_atoms:
        if float(line[2]) >= 1 - OVERLAP_RADIUS/a:
            boundary_atoms.append(copy.deepcopy(line))
            boundary_atoms[-1][2] = str(float(boundary_atoms[-1][2]) - 1)
        elif float(line[3]) >= 1 - OVERLAP_RADIUS/b:
            boundary_atoms.append(copy.deepcopy(line))
            boundary_atoms[-1][3] = str(float(boundary_atoms[-1][3]) - 1)
        elif float(line[4]) >= 1 - OVERLAP_RADIUS/c:
            boundary_atoms.append(copy.deepcopy(line))
            boundary_atoms[-1][4] = str(float(boundary_atoms[-1][4]) - 1)
        elif (float(line[2]) <= OVERLAP_RADIUS/a
              or float(line[3]) <= OVERLAP_RADIUS/b
              or float(line[4]) <= OVERLAP_RADIUS/c):
            boundary_atoms.append(copy.deepcopy(line))
    return boundary_atoms

def check_overlap(list1, list2, a, b, c, alpha, beta, gamma):
    """
    Check for overlap between atoms; return True if no overlap, False otherwise.
    """
    # newlist1 and newlist2 are the positions of all the atoms
    newlist1 = copy.deepcopy(list1)
    newlist2 = copy.deepcopy(list2)

    # newlist3 and newlist4 are the positions of only the boundary atoms
    newlist3 = get_boundary_atoms(newlist1, a, b, c)
    newlist4 = get_boundary_atoms(newlist2, a, b, c)

    # convert x, y, z from fractional to Cartesian for all lists of coordinates
    for line in newlist1 + newlist2 + newlist3 + newlist4:
        a_x = float(line[2]) * a
        b_x = float(line[3]) * cos(gamma) * b
        b_y = float(line[3]) * sin(gamma) * b
        c_x = float(line[4]) * cos(beta) * c
        c_y = (
            float(line[4]) *
            (cos(alpha) - cos(beta) * cos(gamma)) /
            sin(gamma) * c
        )
        c_z = (
            float(line[4]) *
            sqrt(1 - cos(alpha)**2 - cos(beta)**2 - cos(gamma)**2 +
                 2 * cos(alpha) * cos(beta) * cos(gamma))
            / sin(gamma) * c
        )
        line[2] = str(a_x + b_x + c_x)
        line[3] = str(b_y + c_y)
        line[4] = str(c_z)

    # check for boundary overlaps.
    for line3 in newlist3:
        for line4 in newlist4:
            if (
                (float(line4[2]) - float(line3[2]))**2 +
                (float(line4[3]) - float(line3[3]))**2 +
                (float(line4[4]) - float(line3[4]))**2
                < 4 * OVERLAP_RADIUS**2
            ):
                print "Boundary overlap:", (
                    (float(line4[2]) - float(line3[2]))**2 +
                    (float(line4[3]) - float(line3[3]))**2 +
                    (float(line4[4]) - float(line3[4]))**2
                )
                print (
                    line4[0], line4[2], line4[3], line4[4],
                    line3[0], line3[2], line3[3], line3[4]
                )
                return False

    # check for regular overlaps.
    for line in newlist1:
        for line2 in newlist2:
            if (
                (float(line2[2]) - float(line[2]))**2 +
                (float(line2[3]) - float(line[3]))**2 +
                (float(line2[4]) - float(line[4]))**2
                < 4 * OVERLAP_RADIUS**2
            ):
                print "Regular overlap:", (
                    (float(line2[2]) - float(line[2]))**2 +
                    (float(line2[3]) - float(line[3]))**2 +
                    (float(line2[4]) - float(line[4]))**2
                )
                print (
                    line2[0], line2[2], line2[3], line2[4],
                    line[0], line[2], line[3], line[4]
                )
                return False
    return True

def get_lattice_vectors(filename):
    with open(filename) as data:
        elements = data.readline().split()
        a = float(elements[0])
        b = float(elements[1])
        c = float(elements[2])
    return a, b, c

def generate_twofold_interp(current_structure,
                            trial_indices,
                            coordinates_array,
                            outputfile):
    """
    Generate two-fold interpenetrated structures, where SHIFT_DISTANCE
    represents the distance by which we want the interpenetrated lattice
    to be translated. If two points on the structure are closer than
    2 * OVERLAP_RADIUS, then structure is rejected.
    """
    n_a = trial_indices[0]
    n_b = trial_indices[1]
    n_c = trial_indices[2]

    trial_structure = (
        current_structure[:-5] + "_trial_" + str(n_a) + "_" + str(n_b) +
        "_" + str(n_c) + ".cssr"
    )

    new_coordinates = (
        translate_coords(
            SHIFT_DISTANCE / float(coordinates_array[0][0]) * n_a,
            SHIFT_DISTANCE / float(coordinates_array[0][1]) * n_b,
            SHIFT_DISTANCE / float(coordinates_array[0][2]) * n_c,
            coordinates_array[4:]
        )
    )

    if (
        check_overlap(
            new_coordinates,
            coordinates_array[4:],
            float(coordinates_array[0][0]),
            float(coordinates_array[0][1]),
            float(coordinates_array[0][2]),
            radians(float(coordinates_array[1][0])),
            radians(float(coordinates_array[1][1])),
            radians(float(coordinates_array[1][2]))
        )
    ):
        print (
            "No overlap detected for 2-fold interpenetrated trial structure " +
            "with indices " + str(trial_indices) + "."
        )

        with open(outputfile, 'a') as text:
            text.write(
                "No overlap detected for 2-fold interpenetrated trial " +
                "structure with indices " + str(trial_indices) + "\n"
            )

        if SAVE_TRIAL_STRUCTURES:
            with open(trial_structure, 'w') as text:
                coordinates_array[2][0] = str(int(coordinates_array[2][0]) * 2)
                for line in coordinates_array: # non-coordinate information
                    text.write(' '.join(line) + "\n")

                coordinates_array[2][0] = str(int(coordinates_array[2][0]) / 2)
                for index, line in enumerate(new_coordinates):
                    line[0] = str(index + 1 + int(coordinates_array[2][0]))
                    text.write(' '.join(line) + "\n")
        return True
    else:
        print (
            "Overlap detected for 2-fold interpenetrated trial structure " +
            "with indices " + str(trial_indices) + "."
        )

        with open(outputfile, 'a') as text:
            text.write(
                "Overlap detected for 2-fold interpenetrated trial " +
                "structure with indices " + str(trial_indices) + "\n"
            )
        return False

def generate_higher_level_interp(current_structure,
                                 min_indices,
                                 max_lattice_vector,
                                 coordinates_array,
                                 outputfile):
    """
    Generate n-fold interpenetrated structures from the two-fold interpenetrated
    structure with minimum displacement vectors
    """
    n_a = min_indices[0]
    n_b = min_indices[1]
    n_c = min_indices[2]
    new_coordinates = list()

    for n in range(
        1,
        int(
            ceil(max_lattice_vector / (SHIFT_DISTANCE *
                 (n_a**2 + n_b**2 + n_c**2)**0.5)
            )
        )
    ):
        higher_level_indices_a = n_a * n
        higher_level_indices_b = n_b * n
        higher_level_indices_c = n_c * n
        higher_level_indices = [
            higher_level_indices_a,
            higher_level_indices_b,
            higher_level_indices_c
        ]

        new_coordinates += (
            translate_coords(
                SHIFT_DISTANCE / float(coordinates_array[0][0]) *
                higher_level_indices_a,
                SHIFT_DISTANCE / float(coordinates_array[0][1]) *
                higher_level_indices_b,
                SHIFT_DISTANCE / float(coordinates_array[0][2]) *
                higher_level_indices_c,
                coordinates_array[4:]
            )
        )

        if check_overlap(
                new_coordinates,
                coordinates_array[4:],
                float(coordinates_array[0][0]),
                float(coordinates_array[0][1]),
                float(coordinates_array[0][2]),
                radians(float(coordinates_array[1][0])),
                radians(float(coordinates_array[1][1])),
                radians(float(coordinates_array[1][2]))
        ):
            print (
                "No overlap detected for " + str(n+1) +
                 "-fold interpenetrated trial structure with indices " +
                 str(higher_level_indices) + "!"
            )

            with open(outputfile, 'a') as text:
                text.write(
                    "No overlap detected for " + str(n + 1) +
                    "-fold interpenetrated trial structure with indices " +
                    str(higher_level_indices) + "!\n"
                )

            coordinates_array[2][0] = str(int(coordinates_array[2][0])*(n+1))

            higher_level_interp_structure = (
                current_structure[:-5] + "_interp_" + str(n+1) + "_indices_" +
                str(higher_level_indices_a) + "_" + str(higher_level_indices_b) +
                "_" + str(higher_level_indices_c) + ".cssr"
            )

            with open(higher_level_interp_structure, 'w') as text:
                for line in coordinates_array: # non-coordinate information
                    text.write(' '.join(line) + "\n")

                coordinates_array[2][0] = str(
                    int(coordinates_array[2][0]) /
                    (n + 1)
                )

                for index, line in enumerate(new_coordinates):
                    line[0] = str(index + 1 + int(coordinates_array[2][0]))
                    text.write(' '.join(line) + "\n")

            if CSSR_TO_CIF:
                convert_cssr_to_cif(higher_level_interp_structure)
        else:
            print (
                "Overlap detected for " + str(n + 1) +
                "-fold interpenetrated trial structure with indices " +
                str(higher_level_indices) + "."
            )

            with open(outputfile, 'a') as text:
                text.write(
                    "Overlap detected for " + str(n + 1) +
                    "-fold interpenetrated trial structure with indices " +
                    str(higher_level_indices) + "\n"
                )

def main():
    files = [f for f in os.listdir(".") if f.endswith(".cssr")]
    for f in files:
        output_file = str(f[:-5]) + "_output.txt"

        with open(output_file, 'w') as text:
            text.write(
                "Beginning the program " +
                 str(time.strftime("%A %d %B %Y, %I:%M%p")) + "\n"
            )
            text.write(
                "Parameters: SHIFT_DISTANCE, OVERLAP_RADIUS = " +
                str(SHIFT_DISTANCE) + ", " + str(OVERLAP_RADIUS) + "\n"
            )
            text.write("Structure: " + f + "\n")

        with open(f, 'r') as text:
            # make a nested array of text elements; textfile[line][word]
            coordinates_array = list(
                list(line.split()) for line in text.read().split('\n')
            )
            # remove newline at the end of the file
            coordinates_array.pop()

        lattice_vec_a, lattice_vec_b, lattice_vec_c = get_lattice_vectors(f)

        with open(output_file, 'a') as text:
            text.write("Lattice vectors:" + " \n")
            text.write("a = " + str(lattice_vec_a) + " A \n")
            text.write("b = " + str(lattice_vec_b) + " A \n")
            text.write("c = " + str(lattice_vec_c) + " A \n")

        t0 = time.time()

        with open(output_file, 'a') as text:
            text.write(
                "Generating interpenetrated structures and checking overlap\n"
            )

        min_indices = [
            int(ceil(lattice_vec_a / SHIFT_DISTANCE)),
            int(ceil(lattice_vec_b / SHIFT_DISTANCE)),
            int(ceil(lattice_vec_c / SHIFT_DISTANCE))
        ]

        min_indices_magnitude = sqrt(
            (min_indices[0] * SHIFT_DISTANCE)**2 +
            (min_indices[1] * SHIFT_DISTANCE)**2 +
            (min_indices[2] * SHIFT_DISTANCE)**2
        )

        for n_c in range(int(ceil(lattice_vec_c / SHIFT_DISTANCE))):
            for n_b in range(int(ceil(lattice_vec_b / SHIFT_DISTANCE))):
                for n_a in range(int(ceil(lattice_vec_a / SHIFT_DISTANCE))):
                    trial_indices = [n_a, n_b, n_c]

                    indices_magnitude = sqrt(
                        (n_a * SHIFT_DISTANCE)**2 +
                        (n_b * SHIFT_DISTANCE)**2 +
                        (n_c * SHIFT_DISTANCE)**2
                    )

                    if indices_magnitude > min_indices_magnitude:
                        break
                    elif (
                        (indices_magnitude == min_indices_magnitude) and (
                            generate_twofold_interp(
                                f,
                                trial_indices,
                                coordinates_array,
                                output_file
                            )
                        )
                    ):
                        with open(output_file, 'a') as text:
                            text.write(
                                "Same magnitude of displacement for trial " +
                                "structures " + str(trial_indices) + " and " +
                                str(min_indices) + " (magnitude: " +
                                str(indices_magnitude) + " A) \n"
                            )
                            min_indices = list(trial_indices)
                            min_indices_magnitude = sqrt(
                                (min_indices[0] * SHIFT_DISTANCE)**2 +
                                (min_indices[1] * SHIFT_DISTANCE)**2 +
                                (min_indices[2] * SHIFT_DISTANCE)**2
                            )
                    elif (
                        (indices_magnitude < min_indices_magnitude) and (
                            generate_twofold_interp(
                                f,
                                trial_indices,
                                coordinates_array,
                                output_file
                            )
                        )
                    ):
                        min_indices = list(trial_indices)
                        min_indices_magnitude = sqrt(
                            (min_indices[0] * SHIFT_DISTANCE)**2 +
                            (min_indices[1] * SHIFT_DISTANCE)**2 +
                            (min_indices[2] * SHIFT_DISTANCE)**2
                        )
                        with open(output_file, 'a') as text:
                            text.write(
                                "New indices for best trial structure " +
                                str(min_indices) + " \n"
                            )
                            text.write(
                                "min_indices_magnitude " +
                                str(min_indices_magnitude) + " \n"
                            )
                    else:
                        pass

        t1 = time.time()

        with open(output_file, 'a') as text:
            text.write(
                "Time elapsed for finding the best 2-fold interpenetrated " +
                "structure: " + str((t1 - t0) / 60) + " min\n"
            )
        if (
            min_indices != [
                int(ceil(lattice_vec_a)),
                int(ceil(lattice_vec_b)),
                int(ceil(lattice_vec_c))
            ]
        ):
            with open(output_file, 'a') as text:
                text.write(
                    "Generating higher-level interpenetrated structures " +
                    "and checking overlap\n"
                )
            max_lattice_vector = max(
                lattice_vec_a,
                lattice_vec_b,
                lattice_vec_c
            )
            generate_higher_level_interp(
                f,
                min_indices,
                max_lattice_vector,
                coordinates_array,
                output_file
            )
        else:
            print (
                "No successful interpenetrated structures found for " +
                str(f) + "\n"
            )
            with open(output_file, 'a') as text:
                text.write(
                    "No successful interpenetrated structures found for " +
                    str(f) + "\n"
                )
        with open(output_file, 'a') as text:
            text.write(
                "Program complete " +
                str(time.strftime("%A %d %B %Y, %I:%M%p")) + "\n"
            )

if __name__ == "__main__":
    main()
