import numpy as np

RAGGED_SPECTRUM = "C:/Users/justi/Downloads/time_resolved_neutron_spectra_from Radha.xlsx - シート1.csv"

energy = np.loadtxt("../../input/energy.txt")
energy = (energy[1:] + energy[:-1])/2
spectrum = []
with open(RAGGED_SPECTRUM, "r") as f:
	for i, line in enumerate(f.readlines()):
		row = []
		for cell in line.strip().split(","):
			if len(cell) == 0:
				break
			else:
				row.append(float(cell))
		if i%2 == 0:
			input_energy = row
		else:
			input_values = row
			spectrum.append(np.interp(energy, input_energy, input_values))

print(np.shape(spectrum))
np.savetxt("../../input/omega spectrum.txt", np.array(spectrum).T, delimiter="\t")
