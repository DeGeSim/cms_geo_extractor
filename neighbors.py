# %%
import yaml
import uproot
import os
import pickle
import numpy as np
import pandas as pd
import awkward as ak
import matplotlib.pyplot as plt
import seaborn as sns


fngeopic = "output/geometry.pickle"
if os.path.isfile(fngeopic):
    with open(fngeopic, "rb") as f:
        geoD = pickle.load(f)
else:
    with open("output/geometry.yaml", "r") as f:
        geoD = yaml.load(f)
    with open(fngeopic, "wb") as f:
        pickle.dump(geoD, f)

# %%
rf = uproot.open("output/DetIdLUT.root")
arr = rf["analyzer/tree"].arrays()
keydf = ak.to_pandas(arr[0])
keydf = keydf.set_index("globalid")

# %%
def getcellfromrow(row):
    return geoD[row["detectorid"]][row["layerid"]][
        (row["waferortileid.first"], row["waferortileid.second"])
    ][(row["cellid.first"], row["cellid.second"])]


# %%
# Cells per detetor
sns.catplot(x="detectorid", data=keydf, kind="count")

# Types of the detector cells
sns.catplot(x="type", data=keydf, kind="count", hue="detectorid")
# %%

# Check how frequent the neigbors are
ids, counts = np.unique(keydf["next"], return_counts=True)
ids, counts = list(ids), list(counts)

countsL = [(i, c) for i, c in zip(ids, counts) if i != 0 and c != 1]
plt.hist(list(zip(*countsL))[1])


# %%
# Number of cells per player
d8df = keydf[keydf["detectorid"] == 8]
d8df["layerid"].value_counts().plot(kind="bar")
# %%

def issiliconFromId(id):
    if id == 0:
        return None
    else:
        return keydf.loc[id]["issilicon"]


booltab = keydf[["n0", "n1", "n2", "n3", "n4", "n5", "n6", "n7"]].applymap(
    issiliconFromId
)
# %%
def detectorjump(row):
    issilicon = issiliconFromId(row.name)
    data = [e for e in row if e is not None]
    same = 0
    diff = 0
    for e in data:
        if e == issilicon:
            same += 1
        else:
            diff += 1
    return (same, diff)


booltab["nsame"], booltab["ndiff"] = list(zip(*booltab.apply(detectorjump, axis=1)))

# %%
def recurkeys(a):
    return [e for e in a.keys() if "keys" in dir(a[e])]


# %%
# Search for the elements witht the properties
def getproprec(prop: str, iterd: dict):
    if prop in iterd:
        return iterd[prop]
    else:
        return [getproprec(prop, iterd[key]) for key in recurkeys(iterd)]


def fullflatten(flist):
    flat_list = []
    recurflat(flist, flat_list)
    return flat_list


def recurflat(l, flat_list):
    for element in l:
        if type(element) == list:
            recurflat(element, flat_list)
        else:
            flat_list.append(element)


def getprop(prop: str, dx):
    return fullflatten(getproprec(prop, dx))


# %%
propsL = ["issilicon", "type", "z"]
detL = [8, 9, 10]
fig, axs = plt.subplots(len(propsL), 1, figsize=(5, 4 * len(propsL)))

for detectorid in detL:
    for prop, ax in zip(propsL, axs):
        ax.hist(getprop(prop, geoD[detectorid]))
        ax.title.set_text(prop)

# %%
def getnextcellsfreq(detector, layer, threshold=1):
    foo = [
        geoD[detector][layer][wafer][cell]["next"]
        for wafer in recurkeys(geoD[detector][layer])
        for cell in recurkeys(geoD[detector][layer][wafer])
    ]
    return [
        b for a, b in np.array(np.unique(foo, return_counts=True)).T if b > threshold
    ]


def countcells(detector, layer):
    foo = [
        1  # geoD[detector][layer][wafer][cell]
        for wafer in recurkeys(geoD[detector][layer])
        for cell in recurkeys(geoD[detector][layer][wafer])
    ]
    return len(foo)


# %%
for detector in recurkeys(geoD):
    for layer in recurkeys(geoD[detector]):
        print(f"{detector}:\t{layer}:\t{getnextcellsfreq(detector,layer)}")
    print("")
# %%
for detector in recurkeys(geoD):
    for layer in recurkeys(geoD[detector]):
        print(f"{detector}:\t{layer}:\t{countcells(detector,layer)}")
    print("")
# %%
