import time as tm
import numpy as np
from multiprocessing import Pool, cpu_count
from concurrent.futures import ProcessPoolExecutor, as_completed


from Clusterpy.core.toolboxes.cluster.componentsAlg.areamanagernodes import AreaManagerNodes 
from Clusterpy.core.toolboxes.cluster.componentsAlg.memory import ExtendedMemory 
from Clusterpy.core.toolboxes.cluster.componentsAlg.regionmakernodes import RegionMakerNodes  

def constructPossible(am, pRegions, initialSolution, distanceType, distanceStat,
                      selectionType, objectiveFunctionType, Dij, Cio):
    """Create one instance of a region maker"""
    
    rm = RegionMakerNodes(am, pRegions,
                     initialSolution = initialSolution,
                     distanceType = distanceType,
                     distanceStat = distanceStat,
                     selectionType = selectionType,
                     objectiveFunctionType = objectiveFunctionType,
                     Dij = Dij,
                     Cio = Cio)
    return rm

# def execpExhaustiveFunctionalRegions(y, w, pRegions, inits = 3, initialSolution = [],
#                convTabu = 0, tabuLength = 10, Dij={}, Cio={}):
#     """
#     :keyword vars: Area attribute(s) (e.g. ['SAR1','SAR2'])
#     :type vars: list
#     :keyword regions: Number of regions
#     :type regions: integer
#     :keyword wType: Type of first-order contiguity-based spatial matrix: 'rook'
#     or 'queen'. Default value wType = 'rook'.
#     :type wType: string
#     :keyword std: If = 1, then the variables will be standardized.
#     :type std: binary
#     :keyword inits: number of initial feasible solutions to be constructed
#     before applying Tabu Search.
#     :type inits: integer. Default value inits = 5.
#     :keyword initialSolution: List with a initial solution vector. It is useful
#     when the user wants a solution that is not very different from a preexisting
#     solution (e.g. municipalities,districts, etc.). Note that the number of
#     regions will be the same as the number of regions in the initial feasible
#     solution (regardless the value you assign to parameter "regions").
#     IMPORTANT: make sure you are entering a feasible solution and according to
#     the W matrix you selected, otherwise the algorithm will not converge.
#     :type initialSolution: list
#     :keyword convTabu: Stop the search after convTabu nonimproving moves
#     (nonimproving moves are those moves that do not improve the current
#     solution.
#     Note that "improving moves" are different to "aspirational moves").
#     If convTabu=0 the algorithm will stop after Int(M/N) nonimproving moves.
#     Default value convTabu = 0.
#     :type convTabu: integer
#     :keyword tabuLength: Number of times a reverse move is prohibited. Default
#     value *tabuLength = 10*.
#     :type tabuLength: integer
#     :keyword dissolve: If = 1, then you will get a "child" instance of the layer
#     that contains the new regions. Default value *dissolve = 0*.  **Note:**.
#     Each child layer is saved in the attribute *layer.results*.  The first
#     algorithm that you run with *dissolve=1* will have a child layer in
#     *layer.results[0]*; the second algorithm that you run with *dissolve=1* will
#     be in *layer.results[1]*, and so on. You can export a child as a shapefile
#     with *layer.result[<1,2,3..>].exportArcData('filename')*
#     :type dissolve: binary
#     :keyword dataOperations: Dictionary which maps a variable to a list of
#     operations to run on it. The dissolved layer will contains in it's data all
#     the variables specified in this dictionary. Be sure to check the input
#     layer's fieldNames before use this utility.
#     :type dataOperations: dictionary

#     The dictionary structure must be as showed bellow.

#     >>> X = {}
#     >>> X[variableName1] = [function1, function2,....]
#     >>> X[variableName2] = [function1, function2,....]
#     """
#     lenY = len(y) # Quien es y,w pRegions, si execp... es llamada con los argumentos *args y **kargs, entonces
#                     # y,w y pRegions deben alguno de los siguientes: 'pExhaustiveFunctionalRegions', ['dn'], p    
#     lenCio = np.sum(y.values())   # ninguno de los elementos antes mendionados admite el uso de .values()
#     start = 0.0
#     time2 = 0.0

#     listResults=[]

#     print("Running original p-Exhaustive-Functional-Regions algorithm")
#     print("Number of areas: ", lenY)
#     #print 'Numer of development poles: ', int(lenCio)
#     if initialSolution:
#         print("Number of regions: ", len(np.unique(initialSolution)))
#         pRegions = len(set(initialSolution))
#     else:
#         print("Number of regions: ", pRegions)
#     if pRegions >= lenY:
#         message = "\n WARNING: You are aggregating "+str(lenY)+" into"+\
#         str(pRegions)+" regions!!. The number of regions must be an integer"+\
#         " number lower than the number of areas being aggregated"
#         raise Exception(message)

#     if convTabu <= 0:
#         convTabu = lenY/pRegions  #   convTabu = 230*numpy.sqrt(pRegions)
#     distanceType = "EuclideanSquared" #'NOTA: Distances in this case are exogenous, given by the user'
#     distanceStat = "Functional"
#     objectiveFunctionType = "Functional"
#     #selectionType = "Minimum"
#     selectionType = "FullRandom"
#     am = AreaManagerNodes(w, y, Dij, Cio, distanceType)
#     print("Termino de calcular las areas")
#     extendedMemory = ExtendedMemory()
#     print("Termino de extender las areas")

#     pool = Pool(processes = cpu_count())
#     procs = []

#     start = tm.time()

#     for dummy in range(inits):
#         ans = pool.apply_async(constructPossible, [am, pRegions,
#                                                    initialSolution,
#                                                    distanceType,
#                                                    distanceStat,
#                                                    selectionType,
#                                                    objectiveFunctionType,
#                                                    Dij,
#                                                    Cio])
#         procs.append(ans)
 

#     print("Termino for de dummy")
#     results = []
#     for p in procs:
#         try:
#             results.append(p.get())
#         except Exception as e:
#             print(f"Error retrieving result from process: {e}")

#     print("Termino for p")
    
#     tmp_ans = extendedMemory
#     print('TMP', tmp_ans.objInfo)
#     for rm in results:
#         print('Antes de entrar')
#         print (rm.objInfo)
#         print (tmp_ans.objInfo)
#         if rm.objInfo < tmp_ans.objInfo:  #objInfo pertenece a la clase BasicMemory
#             print( 'MEJOR')
#             print( type(tmp_ans))
#             tmp_ans = rm
#             print ('DESPUES')
#             print( type(tmp_ans))
#     rm = tmp_ans

#     extendedMemory.updateExtendedMemory(rm)
#     rm.recoverFromExtendedMemory(extendedMemory)
#     print('----------------')
#     print("INITIAL SOLUTION: ", rm.returnRegions(), "\nINITIAL OF: ", rm.objInfo)
    
#     time2 = tm.time() - start
#     #Sol = rm.regions
#     Of = rm.objInfo
#     seeds = rm.seeds

def execpExhaustiveFunctionalRegions(y, w, pRegions, inits=3, initialSolution=[],
                                     convTabu=0, tabuLength=10, Dij={}, Cio={}):
    lenY = len(y)
    lenCio = np.sum(y.values())
    start = 0.0
    time2 = 0.0

    print("Running original p-Exhaustive-Functional-Regions algorithm")
    print("Number of areas:", lenY)
    
    if initialSolution:
        print("Number of regions:", len(np.unique(initialSolution)))
        pRegions = len(set(initialSolution))
    else:
        print("Number of regions:", pRegions)
    
    if pRegions >= lenY:
        message = f"\n WARNING: You are aggregating {lenY} into {pRegions} regions!!. The number of regions must be an integer number lower than the number of areas being aggregated"
        raise Exception(message)

    if convTabu <= 0:
        convTabu = lenY / pRegions

    distanceType = "EuclideanSquared"
    distanceStat = "Functional"
    objectiveFunctionType = "Functional"
    selectionType = "FullRandom"

    am = AreaManagerNodes(w, y, Dij, Cio, distanceType)
    extendedMemory = ExtendedMemory()

    start = tm.time()

    with ProcessPoolExecutor(max_workers=cpu_count()) as executor:
        futures = [executor.submit(constructPossible, am, pRegions, initialSolution, distanceType,
                                   distanceStat, selectionType, objectiveFunctionType, Dij, Cio) for _ in range(inits)]

        results = []
        print(len(futures))
        for future in as_completed(futures):
            try:
                results.append(future.result(timeout=60))  # Agregar un tiempo de espera
            except Exception as e:
                print(f"Error retrieving result from process: {e}")

    print("Parte problematica")
    tmp_ans = extendedMemory

    for rm in results:
        if rm.objInfo < tmp_ans.objInfo:
            tmp_ans = rm

    rm = tmp_ans
    extendedMemory.updateExtendedMemory(rm)
    rm.recoverFromExtendedMemory(extendedMemory)

    print('----------------')
    print("INITIAL SOLUTION:", rm.returnRegions(), "\nINITIAL OF:", rm.objInfo)
    
    time2 = tm.time() - start
    Of = rm.objInfo
    seeds = rm.seeds

    return Of, seeds


