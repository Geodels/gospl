import pytest
import warnings

warnings.filterwarnings("ignore", category=RuntimeWarning)


def test_gospl_init_destroy():
    import gospl
    from gospl.model import Model

    input = "input/normal.yml"

    # Initialise model
    normal = Model(input, True, False)

    # Destroy model
    normal.destroy()


def test_gospl_run():
    import gospl
    from gospl.model import Model

    input = "input/normal.yml"

    # Initialise model
    normal = Model(input, True, False)

    # Run model
    normal.runProcesses()

    # Destroy model
    normal.destroy()
