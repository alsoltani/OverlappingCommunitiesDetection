import os
import pandas as pd


"""
A short Pandas script to reformat the list of the nodes and their membership.
"""

# Find all files.

path = os.getcwd() + "/Data/DERMemberships/"

for root, dirs, files in os.walk(path):

    print root, dirs, files

    for name in files:
        if "N100_" in name and "Formatted" not in name:

            # Read file.
            with open(os.path.join(root, name), "rb") as f:
                data = map(lambda r: r.split(" "), f.readlines())

            data = map(lambda r: [i for i in r if len(i) > 1], data)
            data = map(lambda r: " ".join(r).rstrip(), data)

            print data

            file_prefix = name.split(".")[0]
            print file_prefix

            with open("Data/DERMemberships/" + file_prefix + "_Formatted.dat", "wb") as f:
                for item in data:
                    f.write("%s\n" % item)
