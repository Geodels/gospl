from gospl.model import Model as sim
from scripts import fuzzyCarb as carbfuzz

input = "input/normal.yml"

# Carbonate controller
carbctrl = carbfuzz.fuzzyCarb()

# Initialise model
normal = sim(input, True, False, carbctrl)

# Run model
normal.runProcesses()

# Destroy model
normal.destroy()
