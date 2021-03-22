# %%
import yaml
import uproot
import os
import pickle
import numpy as np
import pandas as pd
import awkward as ak

# %%
import matplotlib as mpl

# mpl.use("pgf")
import matplotlib.pyplot as plt
import seaborn as sns


# plt.rcParams.update({
#     "text.usetex": True,
#     "pgf.texsystem": "lualatex",
#     "pgf.rcfonts": False,
#     "font.family": "serif",
#     "font.serif": [],
#     "font.sans-serif": [],
#     "font.monospace": [],
#     # "figure.figsize": [default_width, default_width * default_ratsio],
#     "pgf.preamble": "\\usepackage{mymplsetup}"
# })
# plt.rcParams = plt.rcParamsDefault


# %%
rf = uproot.open("output/DetIdLUT.root")
arr = rf["analyzer/tree"].arrays()
keydf = ak.to_pandas(arr[0])
keydf = keydf.set_index("globalid")
keydf.head()
# %%
# Debug code to see the if the arrays are filled correctly
index = [
    "globalid",
    "detectorid",
    "subdetid",
    "layerid",
    "waferortileid.first",
    "waferortileid.second",
    "cellid.first",
    "cellid.second",
    "x",
    "y",
    "celltype",
    "issilicon",
    "next",
    "previous",
    "nneighbors",
    "ngapneighbors",
    "n0",
    "n1",
    "n2",
    "n3",
    "n4",
    "n5",
    "n6",
    "n7",
]
for key in keydf.columns:
    foo = ak.to_pandas(arr[0][key])
    print(key, foo.shape)

# %%
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
# # Types of the detector cells
# sns.catplot(
#     x="celltype",
#     data=keydf,
#     kind="count",
#     hue="detectorid",
# )
# plt.tight_layout()
# plt.savefig("CellTypes.pdf")
# # %%
# # Number of neighbors
# ax = sns.catplot(
#     x="layerid",
#     data=keydf,
#     kind="count",
#     hue="detectorid",
#     legend_out=True,
#     height=4,
#     aspect=2.5,
# )
# locs, labels = plt.xticks()
# plt.xticks(locs[::3], labels[::3])
# plt.savefig("CellsPerLayer.pdf")
# %%
# sns.color_palette("tab10")
# sns.set()
# fig , axes = plt.subplots(4,7,figsize=(35, 20),sharex=True,squeeze=True)
# for i in range(len(axes)):
#     for j in range(len(axes[0])):
#         layerid=i*6+j+1
#         if layerid>=23:
#             continue
#         dfsel=keydf[(keydf["layerid"]==layerid) & (keydf["detectorid"]==8)]
#         print(f"pos {i} {j} => {layerid} ({len(dfsel)})")
#         sns.histplot(
#             x="nneighbors",
#             data=dfsel,
#             hue="detectorid",
#             ax=axes[i][j]
#         )
#         axes[i][j].set_title(f"layer {layerid}")
#         axes[i][j].set_yscale('log')
# plt.title("Number of Neighbors per Cell")
# plt.yscale('log')
# plt.savefig("neighborsHist8.pdf")
# # %%


# plt.cla()
# plt.clf()
# plt.close()

# sns.color_palette("tab10")
# sns.set()
# fig , axes = plt.subplots(4,7,figsize=(35, 20),sharex=True,squeeze=True)
# for i in range(len(axes)):
#     for j in range(len(axes[0])):
#         layerid=i*6+j+1
#         if layerid>=23:
#             continue
#         dfsel=keydf[(keydf["layerid"]==layerid) & (keydf["detectorid"]!=8)]
#         print(f"pos {i} {j} => {layerid} ({len(dfsel)})")
#         sns.histplot(
#             x="nneighbors",
#             data=dfsel,
#             hue="detectorid",
#             ax=axes[i][j]
#         )
#         axes[i][j].set_title(f"layer {layerid}")
#         axes[i][j].set_yscale('log')
# plt.title("Number of Neighbors per Cell")
# plt.yscale('log')
# plt.savefig("neighborsHist.pdf")
# %%

plt.cla()
plt.clf()
plt.close()

sns.color_palette("tab10")
sns.set()
fig, axes = plt.subplots(4, 7, figsize=(35, 20), sharex=True, squeeze=True)
for i in range(len(axes)):
    for j in range(len(axes[0])):
        layerid = i * 6 + j + 1
        if layerid >= 23:
            continue
        dfsel = keydf[
            (keydf["layerid"] == layerid)
            & (keydf["detectorid"] != 8)
            & (keydf["ngapneighbors"] != 0)
        ]
        print(f"pos {i} {j} => {layerid} ({len(dfsel)})")
        sns.histplot(x="ngapneighbors", data=dfsel, hue="detectorid", ax=axes[i][j])
        axes[i][j].set_title(f"layer {layerid}")
        # axes[i][j].set_yscale("log")
plt.title("Number of Neighbors per Cell")
# plt.yscale("log")
plt.savefig("gapneighbors.pdf")


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
