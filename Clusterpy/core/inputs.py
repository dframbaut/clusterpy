import struct

from Clusterpy.core.layer import Layer
from Clusterpy.core.contiguity.weightsFromAreas import weightsFromAreas

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
    if fields[0] != "ID":
        fields = ["ID"] + fields
        for y in data.keys():
            data[y] = [y] + data[y]
    # Check 
            
    layer.fieldNames = fields
    layer.Y = data
    layer.areas, layer.Wqueen, layer.Wrook, layer.shpType = importShape(filename + '.shp')
    layer.calculate_bbox()
    # print(layer.bbox)  # para imprimir toda la lista
    # print(layer.bbox[0]) 
      # Eliminé .defBbox por bbox pq defBbox no aparece dentro  
    #print "Done"    de los atributos de la Clase Layer
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
        fileObj.seek(32)  # Saltar los primeros 32 bytes para llegar al tipo de forma
        shape_type = struct.unpack('<i', fileObj.read(4))[0]  # Leer el tipo de forma

        # Dependiendo del tipo de forma, leer los datos apropiados
        if shape_type == 1:  # Points
            INFO, areas = readPoints(fileObj)
        elif shape_type == 3:  # PolyLine
            INFO, areas = readPolylines(fileObj)
        elif shape_type == 5:  # Polygon
            INFO, areas = readPolygons(fileObj, fileObj.tell())
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

def readPolygons(bodyBytes, start_position):
    """This function reads an ESRI shape file of polygons starting from a given position."""
    bodyBytes.seek(start_position)  # Set the start position
    INFO = {'type': 5}
    AREAS = []
    
    try:
        while True:
            if bodyBytes.read(1) == b"":  
                break  # Check for the end of the file
            bodyBytes.seek(-1, 1)  # Go back one byte to start reading the polygon record
            
            # Skip header parts that are not needed for reading numParts and numPoints
            bodyBytes.seek(43, 1)  # Adjust according to your file's structure
            
            parts_bytes = bodyBytes.read(4)
            points_bytes = bodyBytes.read(4)
            if len(parts_bytes) < 4 or len(points_bytes) < 4:
                break  # Not enough data to read, exit loop
            
            numParts = struct.unpack('<i', parts_bytes)[0]
            numPoints = struct.unpack('<i', points_bytes)[0]
            
            # Read parts indices
            parts = []
            for _ in range(numParts):
                part_index_bytes = bodyBytes.read(4)
                if len(part_index_bytes) < 4:
                    break  # Not enough data to read, exit loop
                parts.append(struct.unpack('<i', part_index_bytes)[0])
                
            # Read coordinates for each point
            area = []
            ring = []
            for i in range(numPoints):
                point_data = bodyBytes.read(16)
                if len(point_data) < 16:
                    break  # Not enough data to read, exit loop
                x, y = struct.unpack('<2d', point_data)
                if i in parts and i != 0:
                    area.append(ring)
                    ring = []
                ring.append((x, y))
            area.append(ring)  # Add the last ring
            AREAS.append(area)

    except struct.error as e:
        print(f"Error during unpacking data: {e}")

    return INFO, AREAS
