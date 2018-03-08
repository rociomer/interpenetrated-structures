import copy, glob, os, subprocess, time
from math import ceil, radians, cos, sin, sqrt
from fileConversion import convertCssrToCif

############################## SET VARIABLES HERE ##############################

#                      findInterpenetratedStructures.py
# This program generates potential interpenetrated structures for all CSSR files
#   in the working directory. Attempts at interpenetrated structures are made 
#   using the displacement set by shiftDistance, and then the structure is
#   checked for overlapping atoms using the hard-sphere radius set by
#   overlapRadius. Set saveTrialStructures to True to save trial interpenetrated 
#   structures (saveTrialStructures == True recommended for debugging). Set 
#   fractionalCoords to True if CSSR input is in fractional coordinates, set to 
#   False if in Cartesian. Set cssrToCif to True to convert the interpenetrated 
#   structures to CIF file format from the default CSSR.

shiftDistance = 1.0
overlapRadius = 1.53 # chosen because C covalent radius set to 0.76
saveTrialStructures = False
fractionalCoords = True 
cssrToCif = True

################################ END VARIABLES #################################

def coordTranslation(x, y, z, oldList): 
# translate atoms by set distance, wrapping them around if they go past the 
#  unit cell boundaries.
    newList = copy.deepcopy(oldList)

    for line in newList:
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
    return newList

def getBoundaryAtoms(allAtoms, a, b, c):
    boundaryAtoms = [] 

    for line in allAtoms:
        if float(line[2]) >= 1-overlapRadius/a:
            boundaryAtoms.append(copy.deepcopy(line))
            boundaryAtoms[-1][2] = str(float(boundaryAtoms[-1][2]) - 1)
        elif float(line[3]) >= 1-overlapRadius/b:
            boundaryAtoms.append(copy.deepcopy(line))
            boundaryAtoms[-1][3] = str(float(boundaryAtoms[-1][3]) - 1)
        elif float(line[4]) >= 1-overlapRadius/c:
            boundaryAtoms.append(copy.deepcopy(line))
            boundaryAtoms[-1][4] = str(float(boundaryAtoms[-1][4]) - 1)
        elif (
            float(line[2]) <= overlapRadius/a 
            or float(line[3]) <= overlapRadius/b 
            or float(line[4]) <= overlapRadius/c
        ):
            boundaryAtoms.append(copy.deepcopy(line))
    return boundaryAtoms

def checkOverlap(overlapRadius, list1, list2, a, b, c, alpha, beta, gamma): 
# check for overlap between atoms; return True if no overlap, False otherwise.
    # newlist1 and newlist2 are the positions of all the atoms
    newlist1 = copy.deepcopy(list1)
    newlist2 = copy.deepcopy(list2)

    # newlist3 and newlist4 are the positions of only the boundary atoms
    newlist3 = getBoundaryAtoms(newlist1, a, b, c)
    newlist4 = getBoundaryAtoms(newlist2, a, b, c)

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
            sqrt(
                1 - cos(alpha)**2 - cos(beta)**2 - cos(gamma)**2 + 
                2 * cos(alpha) * cos(beta) * cos(gamma)
            ) / sin(gamma) * c
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
                < 4 * overlapRadius**2
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
                < 4 * overlapRadius**2
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

def getLatticeVectors(file):
     with open(file) as data:
         elements = data.readline().split()
         a = float(elements[0])
         b = float(elements[1])
         c = float(elements[2])
     return a, b, c 

def generateInterpenetrated(shiftDistance, trialIndices, overlapRadius, \
coordinatesArray): 
# generate two-fold interpenetrated structures, where shiftDistance 
#  represents the distance by which we want the interpenetrated lattice 
#  to be translated, then if two points on the structure are closer than 
#  2 * overlapRadius, the structure is rejected
    n_a = trialIndices[0]
    n_b = trialIndices[1]
    n_c = trialIndices[2]

    trialStructure = (
        f[:-5] + "_trial_" + str(n_a) + "_" + str(n_b) + 
        "_" + str(n_c) + ".cssr"
    )

    newCoordinates = (
        coordTranslation(
            shiftDistance / float(coordinatesArray[0][0]) * n_a, 
            shiftDistance / float(coordinatesArray[0][1]) * n_b, 
            shiftDistance / float(coordinatesArray[0][2]) * n_c, 
            coordinatesArray[4:]
        )
    )

    if (
        checkOverlap(
            overlapRadius, newCoordinates, coordinatesArray[4:], 
            float(coordinatesArray[0][0]), 
            float(coordinatesArray[0][1]), 
            float(coordinatesArray[0][2]), 
            radians(float(coordinatesArray[1][0])), 
            radians(float(coordinatesArray[1][1])), 
            radians(float(coordinatesArray[1][2]))
        )
    ):
        print (
            "No overlap detected for 2-fold interpenetrated trial structure " +
            "with indices " + str(trialIndices) + "."
        )

        with open(outputfile, 'a') as text: 
            text.write(
                "No overlap detected for 2-fold interpenetrated trial " +
                "structure with indices " + str(trialIndices) + "\n"
            )

        if saveTrialStructures == True:
            with open(trialStructure, 'w') as text: 
                coordinatesArray[2][0] = str(int(coordinatesArray[2][0]) * 2)
                for line in coordinatesArray: # non-coordinate information
                    text.write(' '.join(line) + "\n")

                coordinatesArray[2][0] = str(int(coordinatesArray[2][0]) / 2)
                for index, line in enumerate(newCoordinates):
                    line[0] = str(index + 1 + int(coordinatesArray[2][0]))
                    text.write(' '.join(line) + "\n")
        return True
    else:
        print (
            "Overlap detected for 2-fold interpenetrated trial structure " +
            "with indices " + str(trialIndices) + "."
        )

        with open(outputfile, 'a') as text: 
            text.write(
                "Overlap detected for 2-fold interpenetrated trial " +
                "structure with indices " + str(trialIndices) + "\n"
            )
        return False

def generateHigherLevelInterpenetrated(
    shiftDistance, minIndices, maxLatticeVector, 
    overlapRadius, coordinatesArray
): 
# generate n-fold interpenetrated structures from the two-fold interpenetrated 
#  structure with minimum displacement vectors
    n_a = minIndices[0]
    n_b = minIndices[1]
    n_c = minIndices[2]
    newCoordinates = list()

    for n in range(
        1, int(
            ceil(
                maxLatticeVector / 
                (shiftDistance * (n_a**2 + n_b**2 + n_c**2)**0.5)
            )
        )
    ):
        higherLevelIndices_a = n_a * n
        higherLevelIndices_b = n_b * n
        higherLevelIndices_c = n_c * n
        higherLevelIndices = [
            higherLevelIndices_a, 
            higherLevelIndices_b, 
            higherLevelIndices_c
        ]

        newCoordinates += (
            coordTranslation(
                shiftDistance / float(coordinatesArray[0][0]) * 
                higherLevelIndices_a, 
                shiftDistance / float(coordinatesArray[0][1]) * 
                higherLevelIndices_b, 
                shiftDistance / float(coordinatesArray[0][2]) * 
                higherLevelIndices_c, 
                coordinatesArray[4:]
            )
        ) 

        if (
            checkOverlap(
                overlapRadius, newCoordinates, coordinatesArray[4:], 
                float(coordinatesArray[0][0]), 
                float(coordinatesArray[0][1]), 
                float(coordinatesArray[0][2]), 
                radians(float(coordinatesArray[1][0])), 
                radians(float(coordinatesArray[1][1])), 
                radians(float(coordinatesArray[1][2]))
            )
        ):
            print (
                "No overlap detected for " + str(n+1) + 
                 "-fold interpenetrated trial structure with indices " + 
                 str(higherLevelIndices) + "!"
            )

            with open(outputfile, 'a') as text: 
                text.write(
                    "No overlap detected for " + str(n + 1) + 
                    "-fold interpenetrated trial structure with indices " + 
                    str(higherLevelIndices) + "!\n"
                )

            coordinatesArray[2][0] = str(int(coordinatesArray[2][0])*(n+1))

            higherLevelInterpStructure = (
                f[:-5] + "_interp_" + str(n+1) + "_indices_" + 
                str(higherLevelIndices_a) + "_" + str(higherLevelIndices_b) + 
                "_" + str(higherLevelIndices_c) + ".cssr"
            )

            with open(higherLevelInterpStructure, 'w') as text:
                for line in coordinatesArray: # non-coordinate information
                    text.write(' '.join(line) + "\n")

                coordinatesArray[2][0] = str(
                    int(coordinatesArray[2][0]) / 
                    (n + 1)
                )

                for index, line in enumerate(newCoordinates):
                    line[0] = str(index + 1 + int(coordinatesArray[2][0]))
                    text.write(' '.join(line) + "\n")

            if cssrToCif == True:
                convertCssrToCif(higherLevelInterpStructure)
        else:
            print (
                "Overlap detected for " + str(n + 1) + 
                "-fold interpenetrated trial structure with indices " + 
                str(higherLevelIndices) + "."
            )

            with open(outputfile, 'a') as text: 
                text.write(
                    "Overlap detected for " + str(n + 1) + 
                    "-fold interpenetrated trial structure with indices " + 
                    str(higherLevelIndices) + "\n"
                )

# carry out routine
files = [f for f in os.listdir(".") if f.endswith(".cssr")] 
for f in files:
    outputfile = str(f[:-5]) + "_output.txt"

    with open(outputfile, 'w') as text: 
        text.write(
            "Beginning the program " + 
             str(time.strftime("%A %d %B %Y, %I:%M%p")) + "\n"
        )
        text.write(
            "Parameters: shiftDistance, overlapRadius = " + 
            str(shiftDistance) + ", " + str(overlapRadius) + "\n"
        )
        text.write("Structure: " + f + "\n")

    with open(f, 'r') as text: 
        # make a nested array of text elements; textfile[line][word]
        coordinatesArray = list(
            list(line.split()) for line in text.read().split('\n')
        ) 
        # remove newline at the end of the file
        coordinatesArray.pop() 

    latticeVector_a, latticeVector_b, latticeVector_c = getLatticeVectors(f)

    with open(outputfile, 'a') as text:
        text.write("Lattice vectors:" + " \n")
        text.write("a = " + str(latticeVector_a) + " A \n")
        text.write("b = " + str(latticeVector_b) + " A \n")
        text.write("c = " + str(latticeVector_c) + " A \n")

    t0 = time.time()

    with open(outputfile, 'a') as text: 
        text.write(
            "Generating interpenetrated structures and checking overlap\n"
        )

    minIndices = [
        int(ceil(latticeVector_a / shiftDistance)), 
        int(ceil(latticeVector_b / shiftDistance)), 
        int(ceil(latticeVector_c / shiftDistance))
    ]

    minIndicesMagnitude = sqrt(
        (minIndices[0] * shiftDistance)**2 + 
        (minIndices[1] * shiftDistance)**2 + 
        (minIndices[2] * shiftDistance)**2
    )

    for n_c in range(int(ceil(latticeVector_c / shiftDistance))):
        for n_b in range(int(ceil(latticeVector_b / shiftDistance))):
            for n_a in range(int(ceil(latticeVector_a / shiftDistance))):
                trialIndices = [n_a, n_b, n_c]

                IndicesMagnitude = sqrt(
                    (n_a * shiftDistance)**2 + 
                    (n_b * shiftDistance)**2 + 
                    (n_c * shiftDistance)**2
                )

                if IndicesMagnitude > minIndicesMagnitude:
                    break
                elif (
                    (IndicesMagnitude == minIndicesMagnitude) and (
                        generateInterpenetrated(
                            shiftDistance, 
                            trialIndices, 
                            overlapRadius, 
                            coordinatesArray
                        ) == True
                    )
                ):
                    with open(outputfile, 'a') as text:
                        text.write(
                            "Same magnitude of displacement for trial " + 
                            "structures " + str(trialIndices) + " and " + 
                            str(minIndices) + " (magnitude: " + 
                            str(IndicesMagnitude) + " A) \n"
                        )
                        minIndices = list(trialIndices)
                        minIndicesMagnitude = sqrt(
                            (minIndices[0] * shiftDistance)**2 + 
                            (minIndices[1] * shiftDistance)**2 + 
                            (minIndices[2] * shiftDistance)**2
                        )
                elif (
                    (IndicesMagnitude < minIndicesMagnitude) and (
                        generateInterpenetrated(
                            shiftDistance, 
                            trialIndices, 
                            overlapRadius, 
                            coordinatesArray
                        ) == True
                    )
                ):
                    minIndices = list(trialIndices)
                    minIndicesMagnitude = sqrt(
                        (minIndices[0] * shiftDistance)**2 + 
                        (minIndices[1] * shiftDistance)**2 + 
                        (minIndices[2] * shiftDistance)**2
                    )
                    with open(outputfile, 'a') as text:
                        text.write(
                            "New indices for best trial structure " + 
                            str(minIndices) + " \n"
                        )
                        text.write(
                            "minIndicesMagnitude " + 
                            str(minIndicesMagnitude) + " \n"
                        )
                else:
                    pass

    t1 = time.time()

    with open(outputfile, 'a') as text: 
        text.write(
            "Time elapsed for finding the best 2-fold interpenetrated " +
            "structure: " + str((t1 - t0) / 60) + " min\n"
        )
    if (
        minIndices != [
            int(ceil(latticeVector_a)), 
            int(ceil(latticeVector_b)), 
            int(ceil(latticeVector_c))
        ]
    ):
        with open(outputfile, 'a') as text: 
            text.write(
                "Generating higher-level interpenetrated structures " + 
                "and checking overlap\n"
            )
        maxLatticeVector = max(
            latticeVector_a, 
            latticeVector_b, 
            latticeVector_c
        )
        generateHigherLevelInterpenetrated(
            shiftDistance, 
            minIndices, 
            maxLatticeVector, 
            overlapRadius, 
            coordinatesArray
        )
    else:
        print (
            "No successful interpenetrated structures found for " + 
            str(f) + "\n"
        )
        with open(outputfile, 'a') as text: 
            text.write(
                "No successful interpenetrated structures found for " + 
                str(f) + "\n"
            )
    with open(outputfile, 'a') as text:
        text.write(
            "Program complete " + 
            str(time.strftime("%A %d %B %Y, %I:%M%p")) + "\n"
        )
