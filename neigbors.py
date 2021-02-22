import yaml 
import uproot
rf = uproot.open("output/DetIdLUT.root")

with open("output/geometry.yaml", "r") as f:
    d=yaml.load(f)

