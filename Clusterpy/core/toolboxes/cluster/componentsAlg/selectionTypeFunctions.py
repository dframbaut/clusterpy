import numpy as np
from random import randint

def minimumSelection(RegionMaker):
    """
    Select and assign the nearest area to a region
    """
    nInd = 0
    minIndex = 0
    idx = 0
    it = 0 #iterator
    rid = 0
    aid = 0
    val = 0.0
    minVal = 0.0
    values = []
    indicesMin = []
    keys = RegionMaker.candidateInfo.keys()

    if keys:
        #values = [ RegionMaker.candidateInfo[i] for i in keys ]
        #minVal = min(values)
        # random selection for ties
        #indicesMin = [ l[0] for l in enumerate(values) if l[1] == minVal ]
        indicesMin = []
        minVal =  float('Inf')
        for it, key in enumerate(keys):
            val = RegionMaker.candidateInfo[key]
            if val < minVal:
                minVal = val
                indicesMin = [it]
                nInd = 1
            elif val == minVal:
                indicesMin.append(it)
                nInd += 1

        idx = randint(0, nInd - 1)
        minIndex = indicesMin[idx]
        aid = keys[minIndex][0]
        rid = keys[minIndex][1]

        for key in keys:
            if key[0] == aid:
                RegionMaker.candidateInfo.pop(key)
        RegionMaker.assignArea(aid, rid)

def fullRandom(RegionMaker):
    """
    Select and assign randomly an area
    """
    keys = list(RegionMaker.candidateInfo.keys())
    #print('keys content: ', keys, '\n')

# #===========================================================
#     print('Arguments in RegionMaker')
#     for attr_name in dir(RegionMaker):
#         if not attr_name.startswith("__"):
#             try:
#                 attr_value = getattr(RegionMaker, attr_name)
#                 print(f"Attribute {attr_name}: {attr_value}", '\n')
#             except AttributeError:
#                 print(f"Attribute {attr_name} could not be accessed.")
# #===========================================================

    # regionmakernodes line 605, candidateInfo is full of 0's
    values = [ RegionMaker.candidateInfo[i] for i in keys ]
    #print('values content: ', values, '\n')
    if len(values) > 0:
        randomIndex = np.random.randint(0, len(values))
        aid,rid = keys[randomIndex]
        [RegionMaker.candidateInfo.pop(key) for key in keys if key[0] == aid]
        RegionMaker.assignArea(aid, rid)


def indexMultiple(x,value):
    """
    Return indexes in x with multiple values.
    """
    return [ i[0] for i in enumerate(x) if i[1] == value ]



selectionTypeDispatcher = {}
selectionTypeDispatcher["Minimum"] = minimumSelection
selectionTypeDispatcher["FullRandom"] = fullRandom