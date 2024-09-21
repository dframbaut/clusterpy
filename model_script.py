#"from sets import Set"   Sets are built-in functions in Python3
import numpy as np
import Clusterpy
from pandas import *

if __name__ == '__main__':
    # Total Municipios
    colombia = Clusterpy.importArcData(
        "Oficial_Antioquia/Oficial_Antioquia") #"Oficial_Antioquia/Municipios_RIMISP"
    data = len(colombia.Y.keys())
    print('Total de municipios en el mapa: ', data)

    # Leemos la matriz Dij
    xls = ExcelFile('Oficial_Antioquia/diccionario_numeros.xlsx')
    df = xls.parse(xls.sheet_names[0])
    ind = df.set_index(['m']).to_dict()['i']

    # csv.reader(open('Entradas/Oficial_Colombia/dij_def.csv'))
    reader = read_csv('Oficial_Antioquia/dij_def_ant.csv')
    # reader = reader.iloc[1:].reset_index(drop=True)

    Dij = dict()
    for _, row in reader.iterrows():
        if len(row) >= 3:
            try:
                i_int = int(row.iloc[0])
                j_int = int(row.iloc[1])
                d1 = float(row.iloc[2])

                if (
                    i_int not in {27086, 88001, 88564}
                    and j_int not in {27086, 88001, 88564}
                ):
                    i = ind[i_int] 
                    j = ind[j_int]
                    Dij[i, j] = d1 if i < j else 0
            except (ValueError, KeyError):
                continue



    # Leemos la matriz Cio
    #xls = ExcelFile('Entradas/Oficial_Colombia/Cio_Colombia.xlsx')
    #df = xls.parse(xls.sheet_names[0])
        #Cio = df.set_index(['i', 'o']).to_dict()['C']

    reader = read_csv('Oficial_Antioquia/cij_def_ant.csv')
    # reader = reader.iloc[1:].reset_index(drop=True)
    # print("Matriz cij")
    # print(reader.head())

    Cio = dict()
    for _, row in reader.iterrows():
        if len(row) >= 3:
            try:
                i_int = int(row.iloc[0])
                j_int = int(row.iloc[1])
                d1 = float(row.iloc[2])

                if (
                    i_int not in {27086, 88001, 88564}
                    and j_int not in {27086, 88001, 88564}
                ):
                    i = ind[i_int] 
                    j = ind[j_int]
                    Cio[i, j] = d1 if i != j else 0 
            except (ValueError, KeyError):
                continue  # Skip rows with invalid data or missing indices


    items = []
    mpios = []
    for item in Cio.keys():
        items.append(item[1])
        mpios.append(item[0])
    items = set(items)
    mpios = set(mpios)
    print('Municipios', len(mpios))
    if data != len(mpios):
        print('La informacion de las distancias Cio no conicide con la informacion del mapa')

    areas = []
    for area in Dij.keys():
        areas.append(area[1])
    areas = set(areas)
    if data != len(areas):
        print('La informacion de las distancias Dij no conicide con la informacion del mapa')

    nd = np.zeros(len(areas)+1)
    for potential in items:
        nd[potential] = 1

    dn = {}
    for _ in range(len(areas)):
        dn[_] = nd[_]
    colombia.addVariable(['dn'], dn)

    print("Ejecutando Modelos")
    p = 2
    colombia.cluster('pExhaustiveFunctionalRegions', ['dn'], p, convTabu=int(len(
        areas)/p)*200, tabuLength=10, wType='rook', inits=10, dissolve=0, Dij=Dij, Cio=Cio)
    # colombia.exportArcData('Salidas/Colombia_int_'+str(p))
