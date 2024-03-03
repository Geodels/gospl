import pytest
import warnings

warnings.filterwarnings("ignore", category=RuntimeWarning)


def test_gospl_normal_run():
    import gospl
    from gospl.model import Model
    import warnings

    warnings.filterwarnings("ignore", category=RuntimeWarning)

    input = "/live/lib/gospl/tests/input/normal.yml"

    # Initialise model
    normal = Model(input, True, False)

    # Run model
    normal.runProcesses()

    # Reinitialise model
    normal.reInitialiseZ()

    # Destroy model
    normal.destroy()


def test_gospl_backward_run():
    import gospl
    from gospl.model import Model
    import warnings

    warnings.filterwarnings("ignore", category=RuntimeWarning)

    input = "/live/lib/gospl/tests/input/backward.yml"

    # Initialise model
    normal = Model(input, True, False)

    # Run model
    normal.runProcesses()

    # Destroy model
    normal.destroy()


def test_gospl_strat_run():
    import gospl
    from gospl.model import Model
    import warnings

    warnings.filterwarnings("ignore", category=RuntimeWarning)

    input = "/live/lib/gospl/tests/input/strat.yml"

    # Initialise model
    normal = Model(input, True, False)

    # Run model
    normal.runProcesses()

    # Reinitialise model
    normal.reInitialiseZ()

    # Destroy model
    normal.destroy()


def test_gospl_forward_run():
    import gospl
    from gospl.model import Model
    import warnings

    warnings.filterwarnings("ignore", category=RuntimeWarning)

    input = "/live/lib/gospl/tests/input/forward.yml"

    # Initialise model
    normal = Model(input, True, False)

    # Run model
    normal.runProcesses()

    # Destroy model
    normal.destroy()


def test_gospl_restart_run():
    import gospl
    from gospl.model import Model
    import warnings

    warnings.filterwarnings("ignore", category=RuntimeWarning)

    input = "/live/lib/gospl/tests/input/restart.yml"

    # Initialise model
    normal = Model(input, True, False)

    # Run model
    normal.runProcesses()

    # Destroy model
    normal.destroy()


def test_gospl_normalnl_run():
    import gospl
    from gospl.model import Model
    import warnings

    warnings.filterwarnings("ignore", category=RuntimeWarning)

    input = "/live/lib/gospl/tests/input/normalNl.yml"

    # Initialise model
    normal = Model(input, True, False)

    # Run model
    normal.runProcesses()

    # Reinitialise model
    normal.reInitialiseZ()

    # Destroy model
    normal.destroy()

