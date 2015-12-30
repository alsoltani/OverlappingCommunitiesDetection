import os
import pandas as pd


"""
A short Pandas script to reformat the list of the nodes and their membership.
"""

# Find all files.

path = os.getcwd() + "/Data/DERMemberships/"

# Read file.
with open(path + "N1000_K20_MAXK50_MU01.dat", "rb") as f:
    data = map(lambda r: r.split(" "), f.readlines())

data = map(lambda r: [i for i in r if len(i) > 1], data)
data = map(lambda r: " ".join(r), data)

with open(path + "N1000_K20_MAXK50_MU01_Reformatted.dat", "wb") as f:
    for item in data:
        f.write("%s\n" % item)
