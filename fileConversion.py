import subprocess

# update path to Zeo++ network binary below
pathToNetworkBinary = "~/zeo++-0.3/network"
 
def convertCssrToCif(filename):
    subprocess.call(pathToNetworkBinary +" -cif " + filename, shell=True)

