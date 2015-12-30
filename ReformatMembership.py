import os
import pandas as pd


"""
A short Pandas script to reformat the list of the nodes and their membership.
"""

# Find all files.

path = os.getcwd() + "/Data/Networks/"

for root, dirs, files in os.walk(path):
    for name in files:
        if "community.dat" in name:

            print os.path.join(root, name)

            # Read file.
            A = pd.read_csv(os.path.join(root, name), sep="\t")

            A.columns = ["Node", "Membership"]

            # Remove space at the end of membership list.
            A["Membership"] = A["Membership"].apply(lambda r: r.rstrip())

            # Reformat as two-column dataframe.
            A = pd.concat([pd.Series(row["Node"], row["Membership"].split(" "), )
                           for _, row in A.iterrows()]).reset_index()
            A.columns = ["Membership", "Node"]
            A["Node"] = A["Node"].apply(str)

            # Group by membership values.
            grouped_A = A.groupby("Membership")["Node"].apply(lambda x: x.tolist())
            grouped_A = grouped_A.apply(lambda r: " ".join(r))

            # Save.
            file_prefix = root.split("/")[-1]
            grouped_A.to_csv("Data/Memberships/" + file_prefix + ".dat", index=False)
