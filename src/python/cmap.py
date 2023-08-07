import os

import numpy as np
from matplotlib.colors import ListedColormap

CUSTOM_CMAP: dict[str, ListedColormap] = {}

directory = "input" if os.path.isdir("input") else "../../input"
for filename in os.listdir(directory):
	if filename.startswith("cmap_") and filename.endswith(".csv"):
		cmap_data = np.loadtxt(os.path.join(directory, filename), delimiter=",")
		cmap_name = filename[5:-4]
		CUSTOM_CMAP[cmap_name] = ListedColormap(cmap_data, name=cmap_name)
