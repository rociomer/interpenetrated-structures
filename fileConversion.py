import subprocess

### Update path to Zeo++ network binary below
pathToNetworkBinary = '~/Dropbox\ \(LSMO\)/Research/Zeo++/zeo/trunk/network'
 
def convertCssrToCif(filename):
    subprocess.call(pathToNetworkBinary +" -cif " + filename, shell=True)

