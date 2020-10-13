import argparse
from gospl.model import Model as sim

# Parsing command line arguments
parser = argparse.ArgumentParser(
    description="This is a simple entry to run goSPL model.", add_help=True
)
parser.add_argument("-i", "--input", help="Input file name (YAML file)", required=True)
parser.add_argument(
    "-v",
    "--verbose",
    help="True/false option for verbose",
    required=False,
    action="store_true",
    default=False,
)

parser.add_argument(
    "-l",
    "--log",
    help="True/false option for PETSC log",
    required=False,
    action="store_true",
    default=False,
)

args = parser.parse_args()

carbctrl = None

# Reading input file
model = sim(args.input, args.verbose, args.log, carbctrl)

# Running forward model
model.runProcesses()

# Cleaning model
model.destroy()
