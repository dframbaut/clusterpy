from Clusterpy.core.toolboxes.cluster.componentsAlg.dist2Regions import distanceStatDispatcher
from Clusterpy.core.toolboxes.cluster.componentsAlg.areacl import AreaCl

class AreaManager:
    """
    This class contains operations at areal level, including the generation of
    instances of areas, a wide range of area2area and area2region distance
    functions.
    """
    def __init__(self, w, y, distanceType="EuclideanSquared", variance="false"):
        """
        @type w: dictionary
        @param w: With B{key} = area Id, and B{value} = list with Ids of neighbours of
        each area.

        @type y: dictionary
        @param y: With B{key} = area Id, and B{value} = list with attribute
        values.

        @type distanceType: string
        @keyword distanceType: Function to calculate the distance between areas. Default value I{distanceType = 'EuclideanSquared'}.

        @type variance: boolean
        @keyword variance: Boolean indicating if the data have variance matrix. Default value I{variance = 'false'}.
        """
        self.y = y
        self.areas = {}
        self.noNeighs = set([])
        self.variance = variance
        self.distanceType = distanceType
        self.createAreas(w, y)
        self.distanceStatDispatcher = distanceStatDispatcher

    def createAreas(self, w, y):
        """
        Creates instances of areas based on a sparse weights matrix (w) and a
        data array (y).
        """
        n = len(self.y)
        self.distances = {}
        noNeighs = []
        for key in range(n):
            data = y[key]
            try:
                neighbours = w[key]
            except:
                neighbours = {}
                w[key] = {}
            if len(w[key]) == 0:
                self.noNeighs = self.noNeighs | set([key])
            a = AreaCl(key, neighbours, data, self.variance)
            self.areas[key] = a
        if len(self.noNeighs) > 0:
            print("Disconnected areas neighs: ", list(self.noNeighs))