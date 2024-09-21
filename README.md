# Clusterpy:_python2-to-python3

Los siguientes diagramas presentan de manera visual la organización de la librería clusterpy.


En el archivo model_script.py, se utilizan las siguientes 
funciones provenientes de la librería clusterpy

```
. model_script.py
│
├── importArcData
│
├── cluster 
│
└── exportArcData
```

La función importArcData se encuentra en el archivo inputs.py,
de clusterpy\core. Las funciones utilizadas para que la función
importArcData funcionen, son:

```
. inputs (PATH: core\inputs.py)
└── importArcData
    ├── Class Layer
    │   └── Modules used:
    │       ├── (name), (fieldNames), (Y), _defBbox (bbox)
    │       └── (areas), (Wqueen), (Wrook), (shpType)
    │   
    ├── importDBF
    └── importShape
        ├── readShape
        │   ├── readPoints
        │   ├── readPolylines
        │   └── readPolygons
        └── weightsFromAreas (PATH: core\contiguity\weightsFromAreas.py)
```

El atributo **_defBbox** no aparece en la clase Layer

===============================================================

```
. cluster (PATH: core\layer.py)
├── fieldOperation() **
│
├── execpExhaustiveFunctionalRegions() (core\toolboxes)
│   │    
│   ├── AreaManagerNodes
│   │    
│   ├── ExtendedMemory (PATH: toolboxes\componentsAlg)
│   │   │    
│   │   └── updateExtendedMemory
│   │
│   └── constructPossible
│
│       └──RegionMakerNodes
│           ├── objInfo
│           ├── recoverFromExtendedMemory
│           ├── returnRegions
│           └── seeds
│   
├── getVars()
└── addVariable()
```

============================================================
```
. exportArcData (PATH: core\layer.py)
├── shpWriterDis() (PATH: core\outputs.py)
│   ├── shpWriter2()
│   └── shpWriter() 
│
├── exportDBFY()
    │
    ├── getVars()
    └── dbfWriter() (PATH: core\outputs.py)

REMOVE
├── dissolveMap
    │   
    ├── dissolveLayer (PATH: core\geometry\dissolve.py)
    │    │
    │    ├── getScaleRatio (check)
    │    │
    │    ├── scale (check)
    │    │
    │    ├── readAreas 
    │    │     
    │    │
    │    ├── dissolveMap
    │    │   │ 
    │    │   ├── MultiRingCluster
    │    │   ├── isHole (check)
    │    │   │
    │    │   └── Modules: 
    │    │          ├── (Ring Class): ringNeighbors  
    │    │          │
    │    │          └── (MultiRingCluster Class):
    │    │          independentRings, holeRings,
    │    │          contiguosRings, unlinkedSegments
    │    │
    │    ├── reconnect (check)
    │    │   │
    │    │   ├── connector
    │    │   └── getArea
    │    │
    │    └── invertScale (check)
    │     
    └── dissolveData (PATH: core\data\dissolvedata.py)
```