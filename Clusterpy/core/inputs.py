import struct
import geopandas as gpd
from shapely.geometry import Polygon, Point

from Clusterpy.core.layer import Layer
from Clusterpy.core.contiguity.weightsFromAreas import weightsFromAreas
from tqdm import tqdm

def importArcData(filename):
    """Creates a new Layer from a shapefile (<file>.shp)
    
    :param filename: filename without extension 
    :type filename: string
    :rtype: Layer (CP project)

    **Description**

    `ESRI <http://www.esri.com/>`_ shapefile is a binary file used to
    save and transport maps. During the last times it has become
    the most used format for the spatial scientists around the world.

    On clusterPy's "data_examples" folder you can find some shapefiles. To
    load a shapefile in clusterPy just follow the example bellow.

    **Example** ::

        import clusterpy
        china = clusterpy.importArcData("clusterpy/data_examples/china")

    """
    layer = Layer()
    layer.name = filename.split("/")[-1]
    #print("Loading " + filename.split("/")[-1] + ".dbf")
    data, fields, specs = importDBF(filename + '.dbf')
    #print "Loading " + filename + ".shp"
    layer.areas, layer.Wqueen, layer.Wrook, layer.shpType = importShape(filename + '.shp')
    
    #print('layer.Wrook:', layer.Wrook, '\n'),
    #print('layer.Queen:', layer.Wqueen)
    #print('area content: ', layer.areas, '\n')

    if fields[0] != "ID":
        fields = ["ID"] + fields
        for y in data.keys():
            data[y] = [y] + data[y] # -- INFORMATION STORED ON EACH REGION
            data[y] = [value for value in data[y] if isinstance(value, (int, float))]
    # Check 
            
    layer.fieldNames = fields
    layer.Y = data
    return layer
    
def importDBF(filename):
    """Get variables from a dbf file.
    
    :param filename: name of the file (String) including ".dbf"
    :type filename: string
    :rtype: tuple (dbf file Data, fieldNames and fieldSpecs).

    **Example** ::

        import clusterpy
        chinaData = clusterpy.importDBF("clusterpy/data_examples/china.dbf")
    """
    #Y = {}
    fieldNames = []
    fieldSpecs = []
    fileBytes = open(filename, 'rb')
    fileBytes.seek(4, 1)
    numberOfRecords = struct.unpack('i', fileBytes.read(4))[0]
    #print('\n','numberOfRecords: ', numberOfRecords)
    firstDataRecord = struct.unpack('h', fileBytes.read(2))[0]
    #print('firstDataRecord: ', firstDataRecord)
    lenDataRecord = struct.unpack('h', fileBytes.read(2))[0]
    #print('lenDataRecord: ', lenDataRecord, '\n')
    fileBytes.seek(20, 1)
    while fileBytes.tell() < firstDataRecord - 1:
        name_bytes = fileBytes.read(11)
        name = name_bytes.decode('ascii').replace("\x00", "") # name = ''.join(struct.unpack(11 * 'c', fileBytes.read(11))).replace("\x00", "")
        #print("nombre ", name)
        typ = fileBytes.read(1).decode('ascii')
        fileBytes.seek(4, 1)
        siz = struct.unpack('B', fileBytes.read(1))[0]
        dec = struct.unpack('B', fileBytes.read(1))[0]
        spec = (typ, siz, dec)
        fieldNames += [name] # <======================================
        fieldSpecs += [spec] # <======================================
        fileBytes.seek(14, 1)
    fileBytes.seek(1, 1)
    Y = {} 
    for nrec in range(numberOfRecords):
        record = fileBytes.read(lenDataRecord)
        start = 0
        first = 0
        Y[nrec] = [] # <======================================
        #print('\n' , 'fieldSpecs: ', fieldSpecs)
        for nf, field in enumerate(fieldSpecs):
            l = field[1] + 1
            #print('\n' , 'l: ', l)
            dec = field[2]
            #print('\n' , 'dec: ', dec)
            end = start + l + first
            #print('\n' , 'end: ', end)
            value = record[start: end]
            # print('\n', 'record: ', record)
            while value.find(b"  ") != -1:
                value = value.replace(b"  ", b" ")
            if value.startswith(b" "):
                value = value[1:]
            if value.endswith(b" "):
                value = value[:-1]
            if field[0] in ["N", "F", "B", "I", "O"]:
                if dec == 0:
                    value = int(float(value))
                else:
                    value = float(value)
            start = end
            first = -1
            Y[nrec] += [value] # <======================================
    return (Y, fieldNames, fieldSpecs)

def importShape(shapefile):
    """Reads the geographic information stored in a shape file and returns
    them in python objects.
    
    :param shapefile: path to shapefile including the extension ".shp"
    :type shapefile: string
    :rtype: tuple (coordinates(List), Wqueen(Dict), Wrook(Dict)).

    **Example** ::

        import clusterpy
        chinaAreas = clusterpy.importShape("clusterpy/data_examples/china.shp")
    """

    INFO, areas = readShape(shapefile)
    count = 0
    for i in range(len(areas)):
        count += len(areas[i])
    print('length of AREAS: ', count)
    # INFO['type'] = 5
    if INFO['type'] == 5:
        Wqueen, Wrook = weightsFromAreas(areas)
        #Wqueen = {}
        #Wrook = {}
        shpType = 'polygon'
    elif INFO['type'] == 3:
        shpType = 'line'
        Wrook = {}
        Wqueen = {}
    elif INFO['type'] == 1:
        shpType = 'point'
        Wrook = {}
        Wqueen = {}
    return areas, Wqueen, Wrook, shpType

def readShape(filename):
    """ This function automatically detects the type of the shape and then reads an ESRI shapefile of polygons, polylines or points.
    :param filename: name of the file to be read
    :type filename: string
    :rtype: tuple (information about the layer and areas coordinates).
    """
    with open(filename, 'rb') as fileObj:
        # Leer la cabecera para obtener el tipo de forma
        fileObj.seek(32)
        shape_type = struct.unpack('<i', fileObj.read(4))[0]
        #print('tipo de la forma: ', shape_type)
        # Dependiendo del tipo de forma, leer los datos apropiados
        if shape_type == 1:  # Points
            INFO, areas = readPoints(fileObj)
        elif shape_type == 3:  # PolyLine
            INFO, areas = readPolylines(fileObj)
        elif shape_type == 5:  # Polygon
            INFO, areas = readPolygons(filename)
        else:
            raise ValueError("Unsupported shape type")
        
        return INFO, areas

def readPoints(bodyBytes):
    """This function reads an ESRI shapefile of points.

    :param bodyBytes: bytes to be processed
    :type bodyBytes: string
    :rtype: tuple (information about the layer and area coordinates).
    """
    INFO = {}
    INFO['type'] = 1
    AREAS = []
    id = 0
    bb0 = struct.unpack('>d', bodyBytes.read(8))[0]
    bb1 = struct.unpack('>d', bodyBytes.read(8))[0]
    bb2 = struct.unpack('>d', bodyBytes.read(8))[0]
    bb3 = struct.unpack('>d', bodyBytes.read(8))[0]
    bb4 = struct.unpack('>d', bodyBytes.read(8))[0]
    bb5 = struct.unpack('>d', bodyBytes.read(8))[0]
    bb6 = struct.unpack('>d', bodyBytes.read(8))[0]
    bb7 = struct.unpack('>d', bodyBytes.read(8))[0]
    while bodyBytes.read(1) != "":
        bodyBytes.seek(11, 1)
        x = struct.unpack('<d', bodyBytes.read(8))[0]
        y = struct.unpack('<d', bodyBytes.read(8))[0]
        area = [x, y] 
        AREAS = AREAS + [[[tuple(area)]]]
    return INFO, AREAS

def readPolylines(bodyBytes):
    """This function reads a ESRI shape file of lines.

    :param bodyBytes: bytes to be processed
    :type bodyBytes: string
    :rtype: tuple (information about the layer and areas coordinates). 
    """
    INFO = {}
    INFO['type'] = 3
    AREAS=[]
    id = 0
    pos = 100
    bb0 = struct.unpack('>d', bodyBytes.read(8))[0]
    bb1 = struct.unpack('>d', bodyBytes.read(8))[0]
    bb2 = struct.unpack('>d', bodyBytes.read(8))[0]
    bb3 = struct.unpack('>d', bodyBytes.read(8))[0]
    bb4 = struct.unpack('>d', bodyBytes.read(8))[0]
    bb5 = struct.unpack('>d', bodyBytes.read(8))[0]
    bb6 = struct.unpack('>d', bodyBytes.read(8))[0]
    bb7 = struct.unpack('>d', bodyBytes.read(8))[0]
    while bodyBytes.read(1) != "":
        bodyBytes.seek(7, 1)
        bodyBytes.seek(36, 1)
        nParts = struct.unpack('<i', bodyBytes.read(4))[0]
        nPoints = struct.unpack('<i', bodyBytes.read(4))[0]
        r = 1
        parts = []
        while r <= nParts:
            parts += [struct.unpack('<i', bodyBytes.read(4))[0]]
            r += 1
        ring = []
        area = []
        l = 0
        while l < nPoints:
            if l in parts[1:]:
                area += [ring]
                ring = []
            x = struct.unpack('<d', bodyBytes.read(8))[0]
            y = struct.unpack('<d', bodyBytes.read(8))[0]
            l += 1
            ring = ring + [(x, y)]
        area += [ring]
        AREAS = AREAS + [area]
        id += 1
    return INFO, AREAS

def extract_rings(AREAS, geometry):
    area = []
    outer_ring = []
    inner_ring = []

    if isinstance(geometry, Polygon):
        for i in range(len(geometry.exterior.coords[:])):
            outer_ring += [geometry.exterior.coords[i]]
        inner_ring = [interior.coords[:] for interior in geometry.interiors]

    area.append(outer_ring)
    AREAS.append(area)

def readPolygons(filename):
    INFO = {}
    INFO['type'] = 5
    AREAS = []
    
    gdf = gpd.read_file(filename)    
    
    total_iterations = gdf.shape[0]
    pbar = tqdm(total=total_iterations, desc="Constructing of AREAS")
    
    for i in range(gdf.shape[0]):
        geometry = gdf.loc[i, 'geometry']
        extract_rings(AREAS, geometry)
        pbar.update(1)
    pbar.close()
    return INFO, AREAS
# def readPolygons(shapefile):
#     """This function reads an ESRI shapefile of polygons starting from a given position."""
#     INFO = {'type': 5}
#     AREAS = []
#     centroids = {}

#     gdf = gpd.read_file(shapefile)
#     gdf['centroid'] = gdf.geometry.centroid

#     gdf['centroid_coords'] = gdf['centroid'].apply(lambda x: (x.x, x.y))
#     for key, coord in gdf['centroid_coords'].items():
#         centroids[key] = coord

#     for area in gdf['Area_Km']:
#         AREAS.append(area)
    
#     minx, miny, maxx, maxy = gdf.total_bounds
#     bbox = [minx, miny, maxx, maxy]
#     return INFO, AREAS, bbox, centroids

