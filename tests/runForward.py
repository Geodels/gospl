from gospl.model import Model as sim

input = "input/forward.yml"

# Initialise model
forward = sim(input, True, False)

# Run model
forward.runProcesses()

# Destroy model
forward.destroy()
