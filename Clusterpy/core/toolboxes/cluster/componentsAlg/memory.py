from copy import deepcopy
import numpy as np

class BasicMemory:
    """
    Keeps essential information about a solution:
        - Objective function value (self.objInfo)
        - Region assignment for each area (self.regions)
        - Seeds used (self.seeds)
    """

    def __init__(self, objInfo=99999999E10, regions={}, seeds=[]):
        self.objInfo = objInfo
        self.regions = regions
        self.seeds = seeds

    def updateBasicMemory(self, rm):
        """
        Updates BasicMemory from another BasicMemory or ExtendedMemory object.
        """
        self.objInfo = rm.objInfo
        self.regions = rm.regions  # Directly access the 'regions' attribute
        self.seeds = rm.seeds


class ExtendedMemory(BasicMemory):
    """
    Extends BasicMemory to store additional information for backtracking:
        - area2Region: Mapping of areas to their regions
        - region2Area: Mapping of regions to their areas
        - intraBorderingAreas: Areas on region borders
    """

    def __init__(self, objInfo=99999999E10, area2Region={}, region2Area={}, intraBorderingAreas={}):
        super().__init__(objInfo, {})  # Call parent's constructor
        self.area2Region = area2Region
        self.region2Area = region2Area
        self.intraBorderingAreas = intraBorderingAreas

    def updateExtendedMemory(self, rm):
        """
        Updates ExtendedMemory from another BasicMemory or ExtendedMemory object.
        """
        super().updateBasicMemory(rm)  # Update base attributes
        self.area2Region = deepcopy(rm.area2Region)
        self.region2Area = deepcopy(rm.region2Area)
        self.intraBorderingAreas = deepcopy(rm.intraBorderingAreas)

    def returnRegions(self):
        """
        Return regions created
        """
        areasId = self.area2Region.keys()
        areasId = np.sort(areasId).tolist()
        return [self.area2Region[area] for area in areasId]