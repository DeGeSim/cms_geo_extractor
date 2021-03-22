def getcellfromrow(row):
    return geoD[row["detectorid"]][row["layerid"]][
        (row["waferortileid.first"], row["waferortileid.second"])
    ][(row["cellid.first"], row["cellid.second"])]

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


def issiliconFromId(id):
    if id == 0:
        return None
    else:
        return keydf.loc[id]["issilicon"]


# booltab = keydf[["n0", "n1", "n2", "n3", "n4", "n5", "n6", "n7"]].applymap(
#     issiliconFromId
# )
# # %%
# def detectorjump(row):
#     issilicon = issiliconFromId(row.name)
#     data = [e for e in row if e is not None]
#     same = 0
#     diff = 0
#     for e in data:
#         if e == issilicon:
#             same += 1
#         else:
#             diff += 1
#     return (same, diff)


# booltab["nsame"], booltab["ndiff"] = list(zip(*booltab.apply(detectorjump, axis=1)))