from math import ceil, radians, cos, sin, sqrt
import copy, glob, os, subprocess, time

### TO DO: lose dependency on Zeo++
### TO DO: make code versatile, so that it could read fractional coordinates or Cartesian from a .cssr

# Program assumes that fractional coordinates are given.
# n goes up to half the diameter of the largest included sphere.
# Each axis is traversed independently.

### TO DO: Fix the description, as this program is not for using the surface area or pore volume

# SET VARIABLES HERE
# dist is the displacement for the interpenetrated structure
# radius is the hard-sphere radius used to determine if a trial structure has overlapping atoms
# proberad is the radius used in computing the surface area of each trial structure
# trials is the number of insertions used in computing the surface area of each structure

# 7/27: Tested for interpenetrating structures using both the free and included sphere diameters.
# Similar results were observed in both cases, except when the two were separated by the 20 A
# boundary, which happened often. This suggests that the 20 A boundary is inappropriate. The
# boundary has thus been changed to 25 A.
# Changing the radius from 0.5 A to 0.6 A led to one structure, compared to eight,
# for dist = 1 A, linker3_dia.

# 8/22: Modified findVol and findSA so that if a .sa or a .vol file is missing/empty, that 
# structure is assigned an arbitrarily large SA or PoreVol so that it is not chosen.
# This is because in cases where Zeo++ runs into an error with 
# the Voronoi cell ("Error: Voronoi cell of sampled point does not have any nodes") then
# Zeo++ does not create an output file (either .vol or .sa) then the functions just return a
# large number (1000000000000) so that indexMinVol and indexMinSA can still function (before,
# the array "values" was not being created if a single .vol or .sa file was missing from the
# trial structures

#distance, radius, proberad, trials, useVol = 0.25, 0.4, 0.5, 2500, True
distance, radius, proberad, trials, useVol = 0.25, 1.05, 5.0, 250, False

def coordTranslation(x, y, z, oldList): # Translates atoms by set distance, wrapping them around if they go past the unit cell boundaries.
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

def checkOverlap(radius, list1, list2, a, b, c, alpha, beta, gamma): # Checks for overlap between atoms; returns True if no overlap and False otherwise.
    newlist1, newlist2, newlist3, newlist4 = copy.deepcopy(list1), copy.deepcopy(list2), list(), list()
    # newlist3 and newlist4 - boundary atoms
    for line in newlist1:
        if float(line[2]) >= 1-radius/a:
            newlist3.append(copy.deepcopy(line))
            newlist3[-1][2] = str(float(newlist3[-1][2])-1)
        elif float(line[3]) >= 1-radius/b:
            newlist3.append(copy.deepcopy(line))
            newlist3[-1][3] = str(float(newlist3[-1][3])-1)
        elif float(line[4]) >= 1-radius/c:
            newlist3.append(copy.deepcopy(line))
            newlist3[-1][4] = str(float(newlist3[-1][4])-1)
        elif float(line[2]) <= radius/a or float(line[3]) <= radius/b or float(line[4]) <= radius/c:
            newlist3.append(copy.deepcopy(line))
    for line in newlist2:
        if float(line[2]) >= 1-radius/a:
            newlist4.append(copy.deepcopy(line))
            newlist4[-1][2] = str(float(newlist4[-1][2])-1)
        elif float(line[3]) >= 1-radius/b:
            newlist4.append(copy.deepcopy(line))
            newlist4[-1][3] = str(float(newlist4[-1][3])-1)
        elif float(line[4]) >= 1-radius/c:
            newlist4.append(copy.deepcopy(line))
            newlist4[-1][4] = str(float(newlist4[-1][4])-1)
        elif float(line[2]) <= radius/a or float(line[3]) <= radius/b or float(line[4]) <= radius/c:
            newlist4.append(copy.deepcopy(line))
    for line in newlist1 + newlist2 + newlist3 + newlist4:
        a_x, b_x, b_y = float(line[2])*a, float(line[3])*cos(gamma)*b, float(line[3])*sin(gamma)*b
        c_x, c_y = float(line[4])*cos(beta)*c, float(line[4])*(cos(alpha)-cos(beta)*cos(gamma))/sin(gamma)*c
        c_z = float(line[4])*sqrt(1-cos(alpha)**2-cos(beta)**2-cos(gamma)**2+2*cos(alpha)*cos(beta)*cos(gamma))/sin(gamma)*c
        line[2], line[3], line[4] = str(a_x+b_x+c_x), str(b_y+c_y), str(c_z)
    # Checks for boundary overlaps.
    for line3 in newlist3:
        for line4 in newlist4:
            if (float(line4[2])-float(line3[2]))**2 + (float(line4[3])-float(line3[3]))**2 + (float(line4[4])-float(line3[4]))**2 < 4*radius**2:
                print "Boundary overlap:", (float(line4[2])-float(line3[2]))**2 + (float(line4[3])-float(line3[3]))**2 + (float(line4[4])-float(line3[4]))**2
                print line4[0], line4[2], line4[3], line4[4], line3[0], line3[2], line3[3], line3[4]
                return False
    # Checks for regular overlaps.
    for line in newlist1:
        for line2 in newlist2:
            if (float(line2[2])-float(line[2]))**2 + (float(line2[3])-float(line[3]))**2 + (float(line2[4])-float(line[4]))**2 < 4*radius**2:
                print "Regular overlap:", (float(line2[2])-float(line[2]))**2 + (float(line2[3])-float(line[3]))**2 + (float(line2[4])-float(line[4]))**2
                print line2[0], line2[2], line2[3], line2[4], line[0], line[2], line[3], line[4]
                return False
    return True

def findMinDisplacement(filenames): # Returns the indices of the framework with the minimum surface area.
    values = [findSA(file) for file in filenames]
    with open("output.txt", 'a') as text: #Rocio
        text.write("List of surface areas for trial structures: " + str(values) + "\n")
        text.write("Smallest surface area: " + str(min(values)) + " A^2 \n")
    fileSA = filenames[values.index(min(values))]
    ultUnd = len(fileSA)-fileSA[::-1].find('_')-1
    penultUnd = len(fileSA[0:ultUnd])-fileSA[0:ultUnd][::-1].find('_')-1
    antepenultUnd = len(fileSA[0:penultUnd])-fileSA[0:penultUnd][::-1].find('_')-1
    return [fileSA[antepenultUnd+1:penultUnd], fileSA[penultUnd+1:ultUnd], fileSA[ultUnd+1:-5]]

def getLatticeVectors(file): # Returns the diameter of the largest included sphere in a framework.
     with open(file) as data:
         elements = data.readline().split()
         a = float(elements[0])
         b = float(elements[1])
         c = float(elements[2])
     return a, b, c

def generateInterpenetrated(dist, n_a, n_b, n_c, radius, coordinates): # Generates 2-interpenetrated structures.
    # dist represents the distance by which we want the interpenetrated lattice to be translated
    # if two points on the structure are closer than 2*radius, then they are overlapping
    newCoords = coordTranslation(dist/float(coordinates[0][0])*n_a, dist/float(coordinates[0][1])*n_b, dist/float(coordinates[0][2])*n_c, coordinates[4:]) # coordinate information
    if checkOverlap(radius, newCoords, coordinates[4:], float(coordinates[0][0]), float(coordinates[0][1]), float(coordinates[0][2]), radians(float(coordinates[1][0])), radians(float(coordinates[1][1])), radians(float(coordinates[1][2]))):
        print "No overlap detected for 2-fold interpenetrated trial structure " + str(n_a) + "_" + str(n_b) + "_" + str(n_c) + "."
        with open("output.txt", 'a') as text: 
            text.write("No overlap detected for 2-fold interpenetrated trial structure " + str(n_a) + "_" + str(n_b) + "_" + str(n_c) + "\n")
        with open(f[:-5] + "_trial_" + str(n_a) + "_" + str(n_b) + "_" + str(n_c) + ".cssr", 'w') as text: # writing into file
            coordinates[2][0] = str(int(coordinates[2][0])*2)
            for line in coordinates: # non-coordinate information
                text.write(' '.join(line) + "\n")
            coordinates[2][0] = str(int(coordinates[2][0])/2)
            for index, line in enumerate(newCoords):
                line[0] = str(index + 1 + int(coordinates[2][0]))
                text.write(' '.join(line) + "\n")
    else:
        print "Overlap detected for 2-fold interpenetrated trial structure " + str(n_a) + "_" + str(n_b) + "_" + str(n_c) + "."
        with open("output.txt", 'a') as text: 
            text.write("Overlap detected for 2-fold interpenetrated trial structure " + str(n_a) + "_" + str(n_b) + "_" + str(n_c) + "\n")

def generateHigherInterpenetrated(dist, n_a, n_b, n_c, diamInclSphere, radius, coordinates): # Generates n-interpenetrated structures for the structure with minimal surface area.
    newCoords = list()
    for n in range(1, int(ceil(diamInclSphere/(dist*(n_a**2+n_b**2+n_c**2)**0.5)))):
        newCoords += coordTranslation(dist/float(coordinates[0][0])*n_a*n, dist/float(coordinates[0][1])*n_b*n, dist/float(coordinates[0][2])*n_c*n, coordinates[4:]) # coordinate information
        if checkOverlap(radius, newCoords, coordinates[4:], float(coordinates[0][0]), float(coordinates[0][1]), float(coordinates[0][2]), radians(float(coordinates[1][0])), radians(float(coordinates[1][1])), radians(float(coordinates[1][2]))):
            print "No overlap detected for " + str(n+1) + "-fold interpenetrated trial structure " + str(n_a) + "_" + str(n_b) + "_" + str(n_c) + "."
            coordinates[2][0] = str(int(coordinates[2][0])*(n+1))
            with open("output.txt", 'a') as text:
                text.write("Generating higher-level interpenetrated structures, where the optimal displacement is " + str(n_a) + "_" + str(n_b) + "_" + str(n_c) + "\n")
            with open(f[:-5] + "_interp_" + str(n+1) + "_trial_" + str(n_a) + "_" + str(n_b) + "_" + str(n_c) + ".cssr", 'w') as text: # writing into file
                for line in coordinates: # non-coordinate information
                    text.write(' '.join(line) + "\n")
                coordinates[2][0] = str(int(coordinates[2][0])/(n+1))
                for index, line in enumerate(newCoords):
                    line[0] = str(index + 1 + int(coordinates[2][0]))
                    text.write(' '.join(line) + "\n")
            subprocess.call("~/Dropbox\ \(LSMO\)/Research/Zeo++/zeo/trunk/network -cif " + f[:-5] + "_interp_" + str(n+1) + "_trial_" + str(n_a) + "_" + str(n_b) + "_" + str(n_c) + ".cssr", shell=True) #Rocio
        else:
            print "Overlap detected for " + str(n+1) + "-fold interpenetrated trial structure " + str(n_a) + "_" + str(n_b) + "_" + str(n_c) + "."
            with open("output.txt", 'a') as text: 
                text.write("Overlap detected for " + str(n+1) + "-fold interpenetrated trial structure " + str(n_a) + "_" + str(n_b) + "_" + str(n_c) + "\n")

files = [f for f in os.listdir(".") if f.endswith(".cssr")] # cssr files
with open("output.txt", 'w') as text: # Prints time taken into a text file, 'output.txt'
    text.write("Beginning the program at " + str(time.localtime(None)) + "\n")
    text.write("distance, radius, proberad, trials, useVol = " + str(distance) + ", " + str(radius) + ", " + str(proberad) + ", " + str(trials) + ", " + str(useVol) + "\n")
for f in files:
    with open("output.txt", 'a') as text: # Prints time taken into a text file, 'output.txt'
        text.write(f + "\n")
    with open(f, 'r') as text: # reading file
        textfile = list(list(line.split()) for line in text.read().split('\n')) # makes a nested array of text elements; textfile[line][word]
        textfile.pop() # removes newline at the end of the file
#    diamInclSphere = findDiam(f)
    with open("output.txt", 'a') as text: 
        text.write("Computing the diameter of the largest included sphere: " + str(diamInclSphere) + " A \n")
    if diamInclSphere > 20:
        dist = 2*distance
    else:
        dist = distance
    t0 = time.time()
    with open("output.txt", 'a') as text: 
        text.write("Generating interpenetrated structures and checking overlap:\n")
### TO DO: change this loop
    for n_a in range(int(ceil(diamInclSphere/dist/2))):
        for n_b in range(int(ceil(diamInclSphere/dist/2))):
            for n_c in range(int(ceil(diamInclSphere/dist/2))):
                generateInterpenetrated(dist, n_a, n_b, n_c, radius, textfile)
    t1 = time.time()
    with open("output.txt", 'a') as text: # Prints time taken into a text file, 'output.txt'
        text.write(str((t1-t0)/60) + " min for checking overlap of 2-fold interpenetrated structures.\n")
    try:
        if useVol:
            indices = indexMinVol([g for g in os.listdir(".") if g.startswith(f[:-5]) and g.endswith(".cssr") and g != f])
            with open("output.txt", 'a') as text: #Rocio
                text.write("Indices of the best 2-fold interp trial structure: " + str(indices) + "\n")
        else:
            indices = indexMinSA([g for g in os.listdir(".") if g.startswith(f[:-5]) and g.endswith(".cssr") and g != f])
            with open("output.txt", 'a') as text: #Rocio
                text.write("Indices of the best 2-fold interp trial structure: " + str(indices) + "\n")
        t2 = time.time()
        with open("output.txt", 'a') as text:
            text.write(str((t2-t1)/60) + " min for computing pore volume/surface area for 2-interpenetrated structures.\n")
            text.write("Generating higher-level interpenetrated structures.\n")
        generateInterpenetrated2(dist, int(indices[0]), int(indices[1]), int(indices[2]), diamInclSphere, radius, textfile)
    except:
        print "No successful interpenetrated structures for " + str(f) + "\n"
        with open("output.txt", 'a') as text:
            text.write("NO SUCCESSFULLY GENERATED INTERPENETRATED STRUCTURES\n")

### FOR DEBUGGING: Converts all .cssr files for trial structures to cif format
#with open("output.txt", 'a') as text: # Prints time taken into a text file, 'output.txt'
#    text.write("Converting all cssr files into cif files\n")
#files = [f for f in os.listdir(".") if f.endswith(".cssr")]
#for f in files: # Converts all cssr files into cif files.
#    subprocess.call("~/Dropbox\ \(LSMO\)/Research/Zeo++/zeo/trunk/network -cif " + f, shell=True)

with open("output.txt", 'a') as text: # Prints time taken into a text file, 'output.txt'
    text.write("Program complete at " + str(time.localtime(None)) + "\n")
