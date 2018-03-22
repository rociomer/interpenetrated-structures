import subprocess

# update path to Zeo++ network binary below
PATH_TO_NETWORK_BINARY = "~/zeo++-0.3/network"

def convert_cssr_to_cif(filename):
    subprocess.call(PATH_TO_NETWORK_BINARY +" -cif " + filename, shell=True)

