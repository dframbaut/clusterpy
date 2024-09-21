from componentsAlg import somManager
import inputs

import time as tm
import copy

__all__ = ['originalSOM']

def originalSOM(y,w,
        nRows=10,
        nCols=10,
        iters=1000,
        alphaType='linear',
        initialDistribution='Uniform',
        wType='rook',
        fileName=None):
    """
    :keyword vars: Area attribute(s)
    :type vars: list

    :keyword nRows: Number of rows in the lattice
    :type nRows: list

    :keyword nCols: Number of columns in the lattice
    :type nCols: list

    :keyword wType: Type of first-order contiguity-based spatial matrix: 'rook' or 'queen'. Default value wType = 'rook'.
    :type wType: string

    :keyword iters: Number of iterations for the SOM algorithm. Default value iters = 1000.
    :type iters: integer

    :keyword alphaType: Name of the scalar-valued decreasing function which maps iterations onto (0,1) float values. This function is used to define how much modify the BMU neighborhood areas. In clusterPy we have to possible functions: 'linear' (linear decreasing function), or 'quadratic' (quadratic decreasing function). Default value alphaType = 'linear'.
    :type alphaType: string

    :keyword initialDistribution: Data generator process to initialize the neural wights. Default value initialDistribution = 'uniform'.
    :type initialDistribution: string

    :keyword fileName: Parameter used to export neural output layer topology as a shapefile. Default value fileName = None.
    :type fileName: string
    """

    print ("Original Self Organizing Maps")
    start = tm.time()
    print ("---Generating SOM topology---")
    oLayer = inputs.createGrid(nRows, nCols)
    manager = somManager(y,
                 iters,
                 oLayer,
                 alphaType,
                 initialDistribution,
                 wType)
    print ("Done")
    for iter in range(iters):
        manager.clusters = copy.deepcopy(manager.emptyClusters)
        for areaId in manager.order:
            bmu = manager.findBMU(areaId)
            manager.clusters[bmu] += [areaId]
            manager.modifyUnits(bmu, areaId, iter)
        solution = manager.addSolution(iter)
    time = tm.time() - start
    Sol = manager.compressSolution(solution)
    Of = 0
    print ("FINAL SOLUTION: ", Sol)
    print ("FINAL O.F.: ", Of)
    output = { "objectiveFunction": Of,
    "runningTime": time,
    "algorithm": "originalSOM",
    "regions": len(Sol),
    "r2a": Sol,
    "distanceType": None,
    "distanceStat": None,
    "selectionType": None,
    "ObjectiveFuncionType": None,
    "SOMOutputLayer": manager.outputLayer}
    print ("Done")
    if fileName != None:
        manager.outputLayer.exportArcData(fileName)
    return output