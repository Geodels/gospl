from gospl.model import Model as sim

input = "input/backward.yml"

# Initialise model
normal = sim(input, True, False)

# Run model
normal.runProcesses()

# Destroy model
normal.destroy()
