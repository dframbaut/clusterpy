import struct

def read_shp_header(filename):
    """Leer la cabecera de un archivo SHP para obtener el tipo de forma."""
    with open(filename, 'rb') as file:
        file.seek(24)  # Saltar los primeros 24 bytes
        file_length = struct.unpack('>i', file.read(4))[0] * 2  # Longitud del archivo en bytes
        
        file.seek(32)  # Ir al byte 32 donde comienza la descripción de la forma
        shape_type = struct.unpack('<i', file.read(4))[0]  # Leer el tipo de forma
        
        # Leer los límites espaciales (mínimo y máximo de x, y, z, m)
        bounds = struct.unpack('<4d', file.read(32))  # Leer los 4 dobles para los límites Xmin, Ymin, Xmax, Ymax
        
        return {
            'file_length': file_length,
            'shape_type': shape_type,
            'bounds': bounds
        }

# Uso de la función
filename = 'Oficial_Antioquia\Municipios_RIMISP.shp'
header_info = read_shp_header(filename)
print(header_info)
