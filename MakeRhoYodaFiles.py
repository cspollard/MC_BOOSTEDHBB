import yoda
import itertools 
import os, optparse
import re
from sys import stdout


op = optparse.OptionParser()

##add args to name plot so we have meaningful legends later
##add axis labels by default


opts, args = op.parse_args()
name = os.path.splitext(args[0])

aodict = yoda.core.read(args[0], True)

yfiledict = {}

findrho = re.compile(r"_Rho([0-9]*)_")
for k, ao in aodict.iteritems():
    # TODO
    # should be unnecessary soon
    k = k.replace('-', '_')
    ao.path = ao.path.replace('-', '_')

    m = findrho.search(k)
    if m:
        rho = int(m.group(1))
    else:
        continue

    # replace "_RhoXXX_" with "_"
    ao.path = findrho.sub("_", ao.path)
    if rho in yfiledict:
        yfiledict[rho].append(ao)
    else:
        yfiledict[rho] = [ao]

    continue

# loop over rho, aos
for rho, aos in yfiledict.iteritems():
    yoda.write(aos, "Rho%03d.yoda" % rho)
