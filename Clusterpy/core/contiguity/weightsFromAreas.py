

def weightsFromAreas(AREAS):
    """Generates contiguity matrix

    This Function returns the Wrook and Wqueen dictionaries from a set of areas.
   
    :param AREAS: Set of polygons to calculate de weights 
    :type AREAS: list
    
    The example below shows how the user can print the Wrook and the
    Wqueen of a layer.

    >>> import clusterpy
    >>> lay = clusterpy.new()
    >>> lay.Wrook
    >>> lay.Wqueen

    :rtype: (dictionary, dictionary), tuple with queen and rook dictionaries
    """
    areas = AREAS
    ring2areas = {}
    segment2areas = {}
    point2areas = {}
    point2segment = {}
    Wqueen = {}
    Wrook = {}
    for idx in range(len(AREAS)):
        Wqueen[idx] = []
        Wrook[idx] = []
    for a, area in enumerate(areas):
        for r, ring in enumerate(area):
            for p, point in enumerate(ring[0: -1]):
                p1 = tuple(map(lambda x: round(x,3),point))
                p2 = tuple(map(lambda x: round(x,3),ring[p + 1]))
                segment = [p1, p2]
                segment.sort(key=lambda x: x[0]**2 + x[1]**2)
                sortSegment = tuple(segment)
                if sortSegment in segment2areas:
                    segment2areas[sortSegment] += [a]
                    areasRook = segment2areas[sortSegment]
                    for area1 in areasRook:
                        for area2 in areasRook:
                            if area2 not in Wrook[area1] and area2 != area1:
                                Wrook[area1] += [area2]
                                Wrook[area2] += [area1]
                else:
                    segment2areas[sortSegment] = [a]
                if p1 in point2areas:
                    point2areas[p1] += [a]
                    areasQueen = point2areas[p1]
                    for area1 in areasQueen:
                        for area2 in areasQueen:
                            if area2 not in Wqueen[area1] and area2 != area1:
                                Wqueen[area1] += [area2]
                                Wqueen[area2] += [area1]
                else:
                    point2areas[p1] = [a]
    return Wqueen, Wrook