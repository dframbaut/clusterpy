# Custome libraries, located in dissolve.py
from Clusterpy.core.data.CreateVariable import fieldOperation

from Clusterpy.core.toolboxes.pExhaustiveFunctionalRegions import execpExhaustiveFunctionalRegions

from Clusterpy.core.outputs import dbfWriter, shpWriterDis

import numpy as np
import itertools
import time
import copy
import os

class Layer():
    """Main class in clusterPy

    It is an object that represents an original map and all the
    other maps derived from it after running any algorithm.

    The layer object can be also represented as an inverse tree
    with an upper root representing the original map and the
    different branches representing other layers related to the
    root.

    """
    def __init__(self):
        """
        **Attributes**

        * Y: dictionary (attribute values of each feature)
        * fieldNames: list (fieldNames List of attribute names)
        * areas: list (list containing the coordinates of each feature)
        * region2areas: list (list of lenght N (number of areas) with the ID of the region to which each area has been assigned during the last algorithm run)
        * Wqueen: dictionary (spatial contiguity based on queen criterion)
        * Wrook: dictionary (spatial contiguity based on rook criterion)
        * Wcustom: dictionary (custom spatial contiguity based on any other criterion)
        * type: string (layer's geometry type ('polygons','lines','points'))
        * results: list (repository of layer instances from running an algorithm)
        * outputCluster: dictionary (information about different characteristics of a solution (time, parameters, OS, among others))
        * name: string (layer's name; default is 'root')
        * outputDissolve: dictionary (keep information from which the current layer has been created)
        * father: Layer (layer from which the current layer has been generated)
        * bbox: tuple (bounding box)
        """
        self.Y = {}
        self.fieldNames = []
        self.areas = []
        self.region2areas = []
        self.Wqueen = {}
        self.Wrook = {}
        self.customW = {}
        self.shpType = ''
        self.results = []
        self.name = ""
        self.outputCluster = {}
        self.outputCluster['r2a'] = []
        self.outputCluster['r2aRoot'] = []
        self.outputDissolve = {}
        self.outputDissolve['r2a'] = []
        self.outputDissolve['r2aRoot'] = []
        self.father = []
        self.bbox = []
        self.tStats = []

    def addVariable(self, names, values):
        """Adding new variables
        
        :param names: field name
        :type names: list
        :param values: data
        :type values: dictionary

        **Description**

        On the example below the population of China in 1990 is multiplied by 10 and stored on the layer as "10Y1900". Note that using the power of Python and clusterPy together the number of possible new variables is unlimited.

        **Examples**

        **Example 1**::

            import clusterpy
            china = clusterpy.importArcData("clusterpy/data_examples/china")
            Y1990 = china.getVars(["Y1990"])
            MY1990 = {}
            for area_i,pop in enumerate(Y1990):
               MY1990[area_i] = pop*10
            china.addVariable(['10Y1990'], MY1990)

        **Example 2** ::

            import clusterpy
            china = clusterpy.importArcData("clusterpy/data_examples/china")
            chinaData = clusterpy.importCSV("clusterpy/data_examples/china.csv")
            china.addVariable(chinaData[1],chinaData[0])
        """
        #print "Adding variables"
        self.fieldNames += (names)
        for area in range(len(values)):
            if area in self.Y:
                if type(values[area]) is not list:
                    self.Y[area] += [values[area]]
                else:
                    self.Y[area] += values[area]
            else:
                self.Y[area] = [values[area]]

    def getVars(self, *args):
        """Getting subsets of data
        
        :param args: subset data fieldNames.
        :type args: tuple
        :rtype: Dictionary (Data subset)

        **Description**

        This function allows the user to extract a subset of variables from a
        layer object.
        
        **Examples**

        Getting Y1998 and Y1997 from China ::
        
            import clusterpy
            china = clusterpy.importArcData("clusterpy/data_examples/china")
            subset = china.getVars(["Y1998", "Y1997"])
        """
        #print "Getting variables"
        fields = []
        for argument in args:
            if isinstance(argument, list):
                for argumentIn in argument:
                    fields.append(argumentIn)
            else:
                fields.append(argument)
        labels = self.fieldNames
        count = 0
        subY = {}
        #print fields
        #print labels
        for i in self.Y.keys():
            subY[i] = []
        for j in fields:
            for i in range(len(labels)):
                if labels[i] == j:
                    for j in self.Y.keys():
                        #print 'i,j', i, j
                        subY[j] = subY[j] + [self.Y[j][i]]
        #print "Variables successfully extracted"
        return subY

    def cluster(self, *args, **kargs):
        """
        :param args: Basic parameters.
        :type args: tuple
        :param kargs: Optional parameter keywords.
        :type kargs: dictionary
        
        The dataOperations dictionary used by 'dissolveMap <dissolveMap>' could be 
        passed in order to specify which data should be calculated for the dissolved
        layer. The dataOperations dictionary must be:

        >>> X = {}
        >>> X[variableName1] = [function1, function2,....]
        >>> X[variableName2] = [function1, function2,....]
        """

        # *args = 'pExhaustiveFunctionalRegions', ['dn'], p
        # **kargs = convTabu=int(len(areas)/p)*200, tabuLength=10, 
        #   wType='rook', inits=1000, dissolve=0, Dij=Dij, Cio=Cio
        
        algorithm = args[0]   

        # Extracting W type from arguments
        wType = kargs['wType']
        kargs.pop('wType')
        
        # Extracting W according to requirement
        algorithmW = self.Wrook
        
        # Extracting standardize variables
        std = 0

        # Setting dissolve according to requirement
        dissolve = kargs.pop('dissolve')
        
        # Extracting dataOperations
        dataOperations = {}

        # Construction of parameters            
        # print(args[1])
        fieldNames = tuple(args[1]) ## Error: args[2] no se puede convertir en una lista
        algorithmY = self.getVars(*fieldNames)
        if std==1:
            for nn,name in enumerate(fieldNames):
                values = [i[0] for i in self.getVars(name).values()]
                mean_value = np.mean(values)
                std_value = np.std(values)
                newVar = fieldOperation("( " + name + " - " + str(mean_value) + ")/float(" + str(std_value) + ")", algorithmY, fieldNames)
                for nv,val in enumerate(newVar):
                    algorithmY[nv][nn] = val
            if algorithm == "maxpTabu":
                population = fieldNames[-1]
                populationY = self.getVars(population)
                for key in populationY:
                    algorithmY[key][-1] = populationY[key][0]
        args = (algorithmY,algorithmW) + args[2:]  # No existen más de 2 argumentos

        name = algorithm + "_" +  time.strftime("%Y%m%d%H%M%S")
        print("LLamando al algoritmo")
        self.outputCluster[name] = {
            "pExhaustiveFunctionalRegions": lambda *args, **kargs: execpExhaustiveFunctionalRegions(*args, **kargs),
            }[algorithm](*args, **kargs)
        print("Terminando argoritmo")
        self.outputCluster[name]["weightType"] = wType
        self.outputCluster[name]["aggregationVariables"] = fieldNames 
        self.outputCluster[name]["OS"] = os.name
        self.outputCluster[name]["proccesorArchitecture"] = os.getenv('PROCESSOR_ARCHITECTURE')
        self.outputCluster[name]["proccesorIdentifier"] = os.getenv('PROCESSOR_IDENTIFIER')
        self.outputCluster[name]["numberProccesor"] = os.getenv('NUMBER_OF_PROCESSORS')
        sol = self.outputCluster[name]["r2a"]
        self.region2areas = sol
        self.addVariable([name], sol)
        self.outputCluster[name]["fieldName"] = self.fieldNames[-1]

    def exportArcData(self, filename):
        """
        Creates an ESRI shapefile from a clusterPy's layer.
        
        :param filename: shape file name to create, without ".shp"
        :type filename: string 

        **Examples** ::
        
            import clusterpy
            china = clusterpy.importArcData("clusterpy/data_examples/china")
            china.exportArcData("china")
        """
        #print "Writing ESRI files"
        shpWriterDis(self.areas, filename, self.shpType)
        self.exportDBFY(filename)
        #print "ESRI files created"

    def exportDBFY(self, fileName, *args):    
        """Exports the database file

        :param fileName: dbf file name to create, without ".dbf"
        :type fileName: string 
        :param args: variables subset to be exported
        :type args: tuple 

        **Examples** ::

            import clusterpy
            clusterpy.importArcData("clusterpy/data_examples/china")
            china.exportDBFY("china")
        """
        #print "Writing DBF file"
        if args != ():
            Y = self.getVars(self, *args) 
            fieldNames = args
        else:
            Y = self.Y
            fieldNames = self.fieldNames
        fieldspecs = []
        types = Y[0] 
        for i in types:
            itype = str(type(i))
            if 'str' in itype:
                fieldspecs.append(('C', 10, 0))
            else:
                fieldspecs.append(('N', 10, 3))
        records = range(len(Y))
        for i in range(len(Y)):
            if len(fieldNames) == 2:
                records[i] = []
                records[i] = records[i] + Y.values()[i]
            else:
                records[i] = []
                records[i] = records[i] + Y.values()[i]
        dbfWriter(fieldNames, fieldspecs, records, fileName + '.dbf')
        #print "Done"

    def calculate_bbox(self):
            # Supongamos que este método calcula y actualiza self.bbox basado en self.areas
            self.bbox = [min(area) for area in self.areas]