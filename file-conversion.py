import subprocess

### Update path to Zeo++ network binary below
pathToNetworkBinary = '~/Dropbox\ \(LSMO\)/Research/Zeo++/zeo/trunk/'
 
def convertCssrToCif(filename):
    subprocess.call(pathToNetworkBinary +"network -cif " + filename, shell=True)

