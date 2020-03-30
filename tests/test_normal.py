import pytest
import warnings

warnings.filterwarnings("ignore", category=RuntimeWarning)


def test_gospl_normal_run():
    import gospl
    from gospl.model import Model

    input = "/live/lib/gospl/tests/input/normal.yml"

    # Initialise model
    normal = Model(input, True, False)

    # Run model
    normal.runProcesses()

    # Destroy model
    normal.destroy()
