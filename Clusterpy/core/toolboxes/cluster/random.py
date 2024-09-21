import time as tm
from componentsAlg import AreaManager

def execRandom(y, w, regions):
    """Generate random regions
    
    This algorithm aggregates, at random, a set of areas into a predefined
    number of spatially contiguous regions. ::

    :keyword vars: Area attribute(s) (e.g. ['SAR1','SAR2']) 
    :type vars: list
    :keyword regions: Number of regions 
    :type regions: integer
    :keyword wType: Type of first-order contiguity-based spatial matrix: 'rook' or 'queen'. Default value wType = 'rook'. 
    :type wType: string
    :keyword dissolve: If = 1, then you will get a "child" instance of the layer that contains the new regions. Default value = 0. Note:. Each child layer is saved in the attribute layer.results. The first algorithm that you run with dissolve=1 will have a child layer in layer.results[0]; the second algorithm that you run with dissolve=1 will be in layer.results[1], and so on. You can export a child as a shapefile with layer.result[<1,2,3..>].exportArcData('filename')
    :type dissolve: binary
    :keyword dataOperations: Dictionary which maps a variable to a list of operations to run on it. The dissolved layer will contains in it's data all the variables specified in this dictionary. Be sure to check the input layer's fieldNames before use this utility.
    :type dataOperations: dictionary

    The dictionary structure must be as showed bellow.

    >>> X = {}
    >>> X[variableName1] = [function1, function2,....]
    >>> X[variableName2] = [function1, function2,....]
    """

    if regions >= len(y):
        message = "\n WARNING: You are aggregating "+str(len(y))+" into"+\
        str(regions)+" regions!!. The number of regions must be an integer"+\
        " number lower than the number of areas being aggregated"
        raise Exception, message 

    distanceType = "EuclideanSquared" 
    distanceStat = "Centroid"
    objectiveFunctionType = "SS"
    selectionType = "FullRandom"
    am = AreaManager(w, y, distanceType)
    start = tm.time()