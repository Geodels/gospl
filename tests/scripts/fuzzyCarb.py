import numpy as np
import skfuzzy as fuzz
from skfuzzy import control as ctrl


class fuzzyCarb:
    def __init__(self):
        """
        Initialization function.
        """

        self.controlKeys = []
        self.controlBounds = []

        # Accomodation space
        self.controlKeys.append("depth")
        self.controlBounds.append([1.0e-3, 39.0])
        self.waterdepth = np.linspace(-10.0, 40.0, num=101, endpoint=True)
        self.accomodation_control()

        # Carbonate growth ranges:
        self.carbgrowth = np.linspace(0.0, 2.0, num=101, endpoint=True)  # m/kyr
        self.carbonate_control()

        # Build carbonate fuzzy controller
        self.fuzzy_system()

        return

    def accomodation_control(self):

        shallow_depth_sigma = 0.1
        shallow2 = fuzz.gaussmf(self.waterdepth, 0.5, shallow_depth_sigma)
        shallow_depth_sigma = 3
        shallow = fuzz.gaussmf(self.waterdepth, 5, shallow_depth_sigma)
        id1 = np.where(self.waterdepth < 0.5)[0]
        id2 = np.where(self.waterdepth > 5)[0]
        shallow[id1[-1] : id2[0]] = 1.0
        shallow[: id1[-1]] = shallow2[: id1[-1]]
        id = np.where(self.waterdepth <= 0)[0]
        shallow[id] = 0
        shallow[shallow < 0.00001] = 0.0

        medium_depth_sigma = 2
        medium2 = fuzz.gaussmf(self.waterdepth, 9, medium_depth_sigma)
        medium_depth_sigma = 4
        medium = fuzz.gaussmf(self.waterdepth, 15, medium_depth_sigma)
        id1 = np.where(self.waterdepth < 9)[0]
        id2 = np.where(self.waterdepth > 15)[0]
        medium[id1[-1] : id2[0]] = 1.0
        medium[: id1[-1]] = medium2[: id1[-1]]
        id = np.where(self.waterdepth <= 0)[0]
        medium[id] = 0
        medium[medium < 0.00001] = 0.0

        deep_sigma = 5
        deep = 1.0 - fuzz.gaussmf(self.waterdepth, 15, deep_sigma)
        id = np.where(self.waterdepth < 15)[0]
        deep[id] = 0

        # Antecedent object
        self.depth = ctrl.Antecedent(self.waterdepth, "depth")

        # Membership functions population
        self.depth["shallow"] = shallow
        self.depth["medium"] = medium
        self.depth["deep"] = deep

        return

    def carbonate_control(self):

        rate_sigma = 0.15
        low = fuzz.gaussmf(self.carbgrowth, 0.2, rate_sigma)
        id = np.where(self.carbgrowth < 0.2)[0]
        low[id] = 1
        id = np.where(self.carbgrowth > 1)[0]
        low[id] = 0

        mid_sigma = 0.15
        mid2 = fuzz.gaussmf(self.carbgrowth, 0.75, mid_sigma)
        mid = fuzz.gaussmf(self.carbgrowth, 1.25, mid_sigma)
        id1 = np.where(self.carbgrowth < 0.75)[0]
        id2 = np.where(self.carbgrowth > 1.25)[0]
        mid[id1[-1] : id2[0]] = 1.0
        mid[: id1[-1]] = mid2[: id1[-1]]
        id = np.where(self.carbgrowth <= 0)[0]
        mid[id] = 0
        mid[mid < 0.00001] = 0.0

        high_sigma = 0.2
        high = 1.0 - fuzz.gaussmf(self.carbgrowth, 1.2, high_sigma)
        id = np.where(self.carbgrowth < 1.2)[0]
        high[id] = 0

        # Consequent object
        self.growth = ctrl.Consequent(self.carbgrowth, "growth")

        # Membership functions population
        self.growth["low"] = low
        self.growth["mid"] = mid
        self.growth["high"] = high

        return

    def fuzzy_system(self):

        rule1 = ctrl.Rule(self.depth["shallow"], self.growth["high"])
        rule2 = ctrl.Rule(self.depth["medium"], self.growth["mid"])
        rule3 = ctrl.Rule(self.depth["deep"], self.growth["low"])

        growth_ctrl = ctrl.ControlSystem([rule1, rule2, rule3])

        self.carbControlSystem = ctrl.ControlSystemSimulation(growth_ctrl)

        return

    def get_crisp(self, depth=None):

        self.carbControlSystem.input["depth"] = depth

        self.carbControlSystem.compute()

        return self.carbControlSystem.output["growth"]
