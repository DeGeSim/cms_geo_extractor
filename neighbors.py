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
plt.rcParams = plt.rcParamsDefault
mpl.rcParams.update({"font.size": 16})
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
sns.catplot(
    x="celltype",
    data=keydf,
    kind="count",
    hue="detectorid",
)
plt.tight_layout()
plt.savefig("plots/CellTypes.pdf")
# %%
# Number of neighbors
ax = sns.catplot(
    x="layerid",
    data=keydf,
    kind="count",
    hue="detectorid",
    legend_out=True,
    height=4,
    aspect=2.5,
)
locs, labels = plt.xticks()
plt.xticks(locs[::3], labels[::3])
plt.savefig("plots/CellsPerLayer.pdf")
# %%
plt.cla()
plt.clf()
plt.close()

fig, axes = plt.subplots(4, 7, figsize=(35, 20), sharex=True, sharey=False)
for i in range(len(axes)):
    for j in range(len(axes[0])):
        layerid = i * (len(axes[0]) - 1) + j + 1
        if layerid >= 23:
            continue
        dfsel = keydf[
            (keydf["layerid"] == layerid)
            & (keydf["detectorid"] == 8)
        ]
        print(f"pos {i} {j} => {layerid} ({len(dfsel)})")
        sns.histplot(
            x="nneighbors",
            data=dfsel,
            discrete=True,
            multiple="dodge",
            ax=axes[i][j],
        )
        axes[i][j].set_title(f"layer {layerid}")
        axes[i][j].set_yscale('log')
fig.suptitle("Number of Neighbors per Cell in HgCalEE")

plt.tight_layout()
fig.subplots_adjust(top=0.95)
plt.savefig("plots/neighborsHistEE.pdf")
# # %%

# Hadronic Neighbors

plt.cla()
plt.clf()
plt.close()

fig, axes = plt.subplots(4, 7, figsize=(35, 20), sharex=True, sharey=False)
for i in range(len(axes)):
    for j in range(len(axes[0])):
        layerid = i * (len(axes[0]) - 1) + j + 1
        if layerid >= 23:
            continue
        dfsel = keydf[
            (keydf["layerid"] == layerid)
            & (keydf["detectorid"] != 8)
        ]
        print(f"pos {i} {j} => {layerid} ({len(dfsel)})")
        sns.histplot(
            x="nneighbors",
            data=dfsel,
            discrete=True,
            multiple="dodge",
            hue="detectorid",
            ax=axes[i][j],
        )
        axes[i][j].set_title(f"layer {layerid}")
        axes[i][j].set_yscale('log')
fig.suptitle("Number of Neighbors per Cell in HSi (9) and HSc(10)")

plt.tight_layout()
fig.subplots_adjust(top=0.95)
plt.savefig("plots/neighborsHistH.pdf")


# %%
# Added neighbors
plt.cla()
plt.clf()
plt.close()

fig, axes = plt.subplots(4, 7, figsize=(35, 20), sharex=True, sharey=True)
for i in range(len(axes)):
    for j in range(len(axes[0])):
        layerid = i * (len(axes[0]) - 1) + j + 1
        if layerid >= 23:
            continue
        dfsel = keydf[
            (keydf["layerid"] == layerid)
            & (keydf["detectorid"] != 8)
            & (keydf["ngapneighbors"] != 0)
        ]
        print(f"pos {i} {j} => {layerid} ({len(dfsel)})")
        sns.histplot(
            x="ngapneighbors",
            data=dfsel,
            discrete=True,
            multiple="dodge",
            hue="detectorid",
            ax=axes[i][j],
        )
        axes[i][j].set_title(f"layer {layerid}")
fig.suptitle("Number of added Neighbors from the other subdetector in the HGCalH")

plt.tight_layout()
fig.subplots_adjust(top=0.95)
plt.savefig("plots/gapneighbors.pdf")

# %%
# Hadronic Neighbors with gapfixing

plt.cla()
plt.clf()
plt.close()

fig, axes = plt.subplots(4, 7, figsize=(35, 20), sharex=True, sharey=False)
for i in range(len(axes)):
    for j in range(len(axes[0])):
        layerid = i * (len(axes[0]) - 1) + j + 1
        if layerid >= 23:
            continue
        dfsel = keydf[
            (keydf["layerid"] == layerid)
            & (keydf["detectorid"] != 8)
        ]
        dfsel=dfsel.assign(totalneighbors=dfsel["nneighbors"]+dfsel["ngapneighbors"])
        print(f"pos {i} {j} => {layerid} ({len(dfsel)})")
        sns.histplot(
            x="totalneighbors",
            data=dfsel,
            discrete=True,
            multiple="dodge",
            hue="detectorid",
            ax=axes[i][j],
        )
        axes[i][j].set_title(f"layer {layerid}")
        axes[i][j].set_yscale('log')
fig.suptitle("Number of Neighbors per Cell in HSi (9) and HSc(10) after fixing the gaps.")

plt.tight_layout()
fig.subplots_adjust(top=0.95)
plt.savefig("plots/neighborsHistH_with_gapfixing.pdf")

# %%
# Cell Scatterplot
plt.cla()
plt.clf()
plt.close()

fig, axes = plt.subplots(4, 7, figsize=(35, 20), sharex=True, sharey=True, squeeze=True)
for i in range(len(axes)):
    for j in range(len(axes[0])):
        layerid = i * (len(axes[0]) - 1) + j + 1
        if layerid >= 23:
            continue
        dfsel = keydf[(keydf["layerid"] == layerid) & (keydf["detectorid"] != 8)]

        dfsel = pd.concat(
            [
                dfsel[(keydf["detectorid"] == 9)][::3],
                dfsel[(keydf["detectorid"] == 10)],
            ]
        )
        print(f"pos {i} {j} => {layerid} ({len(dfsel)})")
        sns.scatterplot(
            x="x",
            y="y",
            data=dfsel,
            hue="detectorid",
            ax=axes[i][j],
        )
        axes[i][j].set_title(f"layer {layerid}")
fig.suptitle("Cells Scatterplot")

plt.tight_layout()
fig.subplots_adjust(top=0.95)
plt.savefig("plots/scatter.pdf")

# %%
### Arrow plot
# Cell Scatterplot
for layerid in range(1, 23):
    plt.cla()
    plt.clf()
    plt.close()

    fig = plt.figure(figsize=(10, 12))

    dfsel = keydf[(keydf["layerid"] == layerid) & (keydf["detectorid"] != 8)]
    print(f"layer {layerid}")

    sns.scatterplot(
        x="x",
        y="y",
        data=dfsel,
        hue="detectorid",
        # markers=["o","o"]
        style="detectorid",
    )
    for originid, row in dfsel.iterrows():
        if row.ngapneighbors == 0:
            continue
        for i in range(row.nneighbors, row.nneighbors + row.ngapneighbors):
            targetid = row["n" + str(i)]
            dx = keydf.loc[targetid].x - row.x
            dy = keydf.loc[targetid].y - row.y
            plt.arrow(row.x, row.y, dx, dy, head_width=0.4,length_includes_head=True)

    plt.gca().set_title(f"layer {layerid}")
    plt.title("Cells Scatterplot")

    plt.tight_layout()
    fig.subplots_adjust(top=0.95)
    plt.savefig(f"plots/gaparrows-{layerid}.pdf")

# %%
exit(0)
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
