import subprocess

### Update path to Zeo++ network binary below
pathToNetworkBinary = '/hoem/rocio/Zeo++_source_code/zeo/trunk/network'
 
def convertCssrToCif(filename):
    subprocess.call(pathToNetworkBinary +" -cif " + filename, shell=True)

