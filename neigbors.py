#%%
import yaml
import uproot
import os
import pickle
import numpy as np

rf = uproot.open("output/DetIdLUT.root")

fngeopic = "output/geometry.pickle"
if os.path.isfile(fngeopic):
    with open(fngeopic, "rb") as f:
        d = pickle.load(f)
else:
    with open("output/geometry.yaml", "r") as f:
        d = yaml.load(f)
    with open(fngeopic, "wb") as f:
        pickle.dump(d, f)

#%%


def recurkeys(a):
    return [e for e in a.keys() if hasattr(a[e], "__getitem__")]


# %%
def getnextcellsfreq(detector,layer, threshold=1):
    foo = [
        d[detector][0][layer][wafer][cell]["next"]
        for wafer in recurkeys(d[detector][0][layer])
        for cell in recurkeys(d[detector][0][layer][wafer])
    ]
    return [
        b
        for a, b in np.array(np.unique(foo, return_counts=True)).T
        if b > threshold
    ]

def countcells(detector,layer):
    foo = [
        1 #d[detector][0][layer][wafer][cell]
        for wafer in recurkeys(d[detector][0][layer])
        for cell in recurkeys(d[detector][0][layer][wafer])
    ]
    return len(foo)
# %%
for detector in recurkeys(d):
    for layer in recurkeys(d[detector][0]):
        print(f"{detector}:\t{layer}:\t{getnextcellsfreq(detector,layer)}")
    print("")
# %%
for detector in recurkeys(d):
    for layer in recurkeys(d[detector][0]):
        print(f"{detector}:\t{layer}:\t{countcells(detector,layer)}")
    print("")
# %%
