import time as tm
import numpy as np
from multiprocessing import Pool, cpu_count
import Clusterpy

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

def execpExhaustiveFunctionalRegions(y, w, pRegions, inits = 3, initialSolution = [],
               convTabu = 0, tabuLength = 10, Dij={}, Cio={}):
    """
    :keyword vars: Area attribute(s) (e.g. ['SAR1','SAR2'])
    :type vars: list
    :keyword regions: Number of regions
    :type regions: integer
    :keyword wType: Type of first-order contiguity-based spatial matrix: 'rook'
    or 'queen'. Default value wType = 'rook'.
    :type wType: string
    :keyword std: If = 1, then the variables will be standardized.
    :type std: binary
    :keyword inits: number of initial feasible solutions to be constructed
    before applying Tabu Search.
    :type inits: integer. Default value inits = 5.
    :keyword initialSolution: List with a initial solution vector. It is useful
    when the user wants a solution that is not very different from a preexisting
    solution (e.g. municipalities,districts, etc.). Note that the number of
    regions will be the same as the number of regions in the initial feasible
    solution (regardless the value you assign to parameter "regions").
    IMPORTANT: make sure you are entering a feasible solution and according to
    the W matrix you selected, otherwise the algorithm will not converge.
    :type initialSolution: list
    :keyword convTabu: Stop the search after convTabu nonimproving moves
    (nonimproving moves are those moves that do not improve the current
    solution.
    Note that "improving moves" are different to "aspirational moves").
    If convTabu=0 the algorithm will stop after Int(M/N) nonimproving moves.
    Default value convTabu = 0.
    :type convTabu: integer
    :keyword tabuLength: Number of times a reverse move is prohibited. Default
    value *tabuLength = 10*.
    :type tabuLength: integer
    :keyword dissolve: If = 1, then you will get a "child" instance of the layer
    that contains the new regions. Default value *dissolve = 0*.  **Note:**.
    Each child layer is saved in the attribute *layer.results*.  The first
    algorithm that you run with *dissolve=1* will have a child layer in
    *layer.results[0]*; the second algorithm that you run with *dissolve=1* will
    be in *layer.results[1]*, and so on. You can export a child as a shapefile
    with *layer.result[<1,2,3..>].exportArcData('filename')*
    :type dissolve: binary
    :keyword dataOperations: Dictionary which maps a variable to a list of
    operations to run on it. The dissolved layer will contains in it's data all
    the variables specified in this dictionary. Be sure to check the input
    layer's fieldNames before use this utility.
    :type dataOperations: dictionary

    The dictionary structure must be as showed bellow.

    >>> X = {}
    >>> X[variableName1] = [function1, function2,....]
    >>> X[variableName2] = [function1, function2,....]
    """

    #print("y: ", y)
    #print("w: ", w)
    lenY = len(y)     
    #lenCio = np.sum(y.values())  
    start = 0.0
    time2 = 0.0
    listResults=[]

    print("Running original p-Exhaustive-Functional-Regions algorithm")
    print("Number of areas: ", lenY)
    #print('Numer of development poles: ', int(lenCio))
    if initialSolution:
        print("Number of regions: ", len(np.unique(initialSolution)))
        pRegions = len(set(initialSolution))
    else:
        print("Number of regions: ", pRegions)
    if pRegions >= lenY:
        message = "\n WARNING: You are aggregating "+str(lenY)+" into"+\
        str(pRegions)+" regions!!. The number of regions must be an integer"+\
        " number lower than the number of areas being aggregated"
        raise Exception(message)

    # if convTabu <= 0:
    #     convTabu = lenY/pRegions  #   convTabu = 230*numpy.sqrt(pRegions)
    distanceType = "EuclideanSquared" #'NOTA: Distances in this case are exogenous, given by the user'
    distanceStat = "Functional"
    objectiveFunctionType = "Functional"
    #selectionType = "Minimum"
    selectionType = "FullRandom"
    am = AreaManagerNodes(w, y, Dij, Cio, distanceType)
    extendedMemory = ExtendedMemory()


    #pool = Pool(processes = cpu_count())
    #procs = []

    start = tm.time()

    ans = constructPossible(am, pRegions, initialSolution,
                            distanceType, distanceStat,
                            selectionType, objectiveFunctionType,
                            Dij, Cio)
    
    results = [ans]
    # for dummy in range(inits):
    #     ans = pool.apply_async(constructPossible, [am, pRegions,
    #                                                initialSolution,
    #                                                distanceType,
    #                                                distanceStat,
    #                                                selectionType,
    #                                                objectiveFunctionType,
    #                                                Dij,
    #                                                Cio])
    #     procs.append(ans)

    # results = []
    # for p in procs:
    #     try:
    #         results.append(p.get())
    #     except Exception as e:
    #         print(f"Error retrieving result from process: {e}")
    
    
    tmp_ans = extendedMemory
    print('TMP', extendedMemory.objInfo)
    for rm in results:
        print('Anteds de entrar')
        print ('Object Info from construct possible: ', rm.objInfo)
        print ('Object Info from Extended memory: ', tmp_ans.objInfo)

        # print('Arguments in obj Info from Extended Memory')
        # for attr_name in dir(tmp_ans.objInfo):
        #     if not attr_name.startswith("__"):
        #         try:
        #             attr_value = getattr(tmp_ans.objInfo, attr_name)
        #             print(f"Attribute {attr_name}: {attr_value}", '\n')
        #         except AttributeError:
        #             print(f"Attribute {attr_name} could not be accessed.")
    
        if rm.objInfo < tmp_ans.objInfo:  
            print( 'MEJOR')
            print( type(tmp_ans))
            tmp_ans = rm
            print ('DESPUES')
            print( type(tmp_ans))
    rm = tmp_ans

    extendedMemory.updateExtendedMemory(rm)
    rm.recoverFromExtendedMemory(extendedMemory)
    print('----------------')
    print("INITIAL SOLUTION: ", rm.returnRegions(), "\nINITIAL OF: ", rm.objInfo)
    
    time2 = tm.time() - start
    #Sol = rm.regions
    Of = rm.objInfo
    seeds = rm.seeds

    # print('Arguments in rm from Region Maker')
    # for attr_name in dir(rm):
    #     if not attr_name.startswith("__"):
    #         try:
    #             attr_value = getattr(rm, attr_name)
    #             print(f"Attribute {attr_name}: {attr_value}", '\n')
    #         except AttributeError:
    #             print(f"Attribute {attr_name} could not be accessed.")

    #print('Elements fro Tabu Search: ', '\n', 'tabu length: ', tabuLength, '\n', 'ConvTabu: ', convTabu)
    rm.tabuMove(tabuLength=tabuLength, convTabu=convTabu, typeTabu = 'exact', Dij = Dij, Cio = Cio, seeds = rm.seeds)
    time2 = tm.time() - start
    Sol = rm.regions
    Of = rm.objInfo
    seeds = rm.seeds
    
    #rm = impactSeeds(rm)
    time2 = tm.time() - start
    Sol = rm.regions
    Of = rm.objInfo
    seeds = rm.seeds

    for res in rm.resList:
        listResults.append(res)

    dicc = {0: u'MEDELL\xc3\x8dN', 1: u'ABEJORRAL', 2: u'ABRIAQU\xc3\x8d', 3: u'ALEJANDR\xc3\x8dA', 4: u'AMAG\xc3\x81', 5: u'AMALFI', 6: u'ANDES', 7: u'ANGEL\xc3\u201cPOLIS', 8: u'ANGOSTURA', 9: u'ANOR\xc3\x8d', 10: u'SANTA F\xc3\u2030 DE ANTIOQUIA', 11: u'ANZ\xc3\x81', 12: u'APARTAD\xc3\u201c', 13: u'ARBOLETES', 14: u'ARGELIA', 15: u'ARMENIA', 16: u'BARBOSA', 17: u'BELMIRA', 18: u'BELLO', 19: u'BETANIA', 20: u'BETULIA', 21: u'CIUDAD BOL\xc3\x8dVAR', 22: u'BRICE\xc3\u2018O', 23: u'BURITIC\xc3\x81', 24: u'C\xc3\x81CERES', 25: u'CAICEDO', 26: u'CALDAS', 27: u'CAMPAMENTO', 28: u'CA\xc3\u2018ASGORDAS', 29: u'CARACOL\xc3\x8d', 30: u'CARAMANTA', 31: u'CAREPA', 32: u'EL CARMEN DE VIBORAL', 33: u'CAROLINA', 34: u'CAUCASIA', 35: u'CHIGOROD\xc3\u201c', 36: u'CISNEROS', 37: u'COCORN\xc3\x81', 38: u'CONCEPCI\xc3\u201cN', 39: u'CONCORDIA', 40: u'COPACABANA', 41: u'DABEIBA', 42: u'DONMAT\xc3\x8dAS', 43: u'EB\xc3\u2030JICO', 44: u'EL BAGRE', 45: u'ENTRERR\xc3\x8dOS', 46: u'ENVIGADO', 47: u'FREDONIA', 48: u'FRONTINO', 49: u'GIRALDO', 50: u'GIRARDOTA', 51: u'G\xc3\u201cMEZ PLATA', 52: u'GRANADA', 53: u'GUADALUPE', 54: u'GUARNE', 55: u'GUATAP\xc3\u2030', 56: u'HELICONIA', 57: u'HISPANIA', 58: u'ITAG\xc3\u0153\xc3\x8d', 59: u'ITUANGO', 60: u'JARD\xc3\x8dN', 61: u'JERIC\xc3\u201c', 62: u'LA CEJA', 63: u'LA ESTRELLA', 64: u'LA PINTADA', 65: u'LA UNI\xc3\u201cN', 66: u'LIBORINA', 67: u'MACEO', 68: u'MARINILLA', 69: u'MONTEBELLO', 70: u'MURIND\xc3\u201c', 71: u'MUTAT\xc3\x81', 72: u'NARI\xc3\u2018O', 73: u'NECOCL\xc3\x8d', 74: u'NECH\xc3\x8d', 75: u'OLAYA', 76: u'PE\xc3\u2018OL', 77: u'PEQUE', 78: u'PUEBLORRICO', 79: u'PUERTO BERR\xc3\x8dO', 80: u'PUERTO NARE', 81: u'PUERTO TRIUNFO', 82: u'REMEDIOS', 83: u'RETIRO', 84: u'RIONEGRO', 85: u'SABANALARGA', 86: u'SABANETA', 87: u'SALGAR', 88: u'SAN ANDR\xc3\u2030S DE CUERQU\xc3\x8dA', 89: u'SAN CARLOS', 90: u'SAN FRANCISCO', 91: u'SAN JER\xc3\u201cNIMO', 92: u'SAN JOS\xc3\u2030 DE LA MONTA\xc3\u2018A', 93: u'SAN JUAN DE URAB\xc3\x81', 94: u'SAN LUIS', 95: u'SAN PEDRO DE LOS MILAGROS', 96: u'SAN PEDRO DE URAB\xc3\x81', 97: u'SAN RAFAEL', 98: u'SAN ROQUE', 99: u'SAN VICENTE FERRER', 100: u'SANTA B\xc3\x81RBARA', 101: u'SANTA ROSA DE OSOS', 102: u'SANTO DOMINGO', 103: u'EL SANTUARIO', 104: u'SEGOVIA', 105: u'SONS\xc3\u201cN', 106: u'SOPETR\xc3\x81N', 107: u'T\xc3\x81MESIS', 108: u'TARAZ\xc3\x81', 109: u'TARSO', 110: u'TITIRIB\xc3\x8d', 111: u'TOLEDO', 112: u'TURBO', 113: u'URAMITA', 114: u'URRAO', 115: u'VALDIVIA', 116: u'VALPARA\xc3\x8dSO', 117: u'VEGACH\xc3\x8d', 118: u'VENECIA', 119: u'VIG\xc3\x8dA DEL FUERTE', 120: u'YAL\xc3\x8d', 121: u'YARUMAL', 122: u'YOLOMB\xc3\u201c', 123: u'YOND\xc3\u201c', 124: u'ZARAGOZA', 125: u'BARRANQUILLA', 126: u'BARANOA', 127: u'CAMPO DE LA CRUZ', 128: u'CANDELARIA', 129: u'GALAPA', 130: u'JUAN DE ACOSTA', 131: u'LURUACO', 132: u'MALAMBO', 133: u'MANAT\xc3\x8d', 134: u'PALMAR DE VARELA', 135: u'PIOJ\xc3\u201c', 136: u'POLONUEVO', 137: u'PONEDERA', 138: u'PUERTO COLOMBIA', 139: u'REPEL\xc3\u201cN', 140: u'SABANAGRANDE', 141: u'SABANALARGA', 142: u'SANTA LUC\xc3\x8dA', 143: u'SANTO TOM\xc3\x81S', 144: u'SOLEDAD', 145: u'SUAN', 146: u'TUBAR\xc3\x81', 147: u'USIACUR\xc3\x8d', 148: u'BOGOT\xc3\x81, D.C.', 149: u'CARTAGENA DE INDIAS', 150: u'ACH\xc3\x8d', 151: u'ALTOS DEL ROSARIO', 152: u'ARENAL', 153: u'ARJONA', 154: u'ARROYOHONDO', 155: u'BARRANCO DE LOBA', 156: u'CALAMAR', 157: u'CANTAGALLO', 158: u'CICUCO', 159: u'C\xc3\u201cRDOBA', 160: u'CLEMENCIA', 161: u'EL CARMEN DE BOL\xc3\x8dVAR', 162: u'EL GUAMO', 163: u'EL PE\xc3\u2018\xc3\u201cN', 164: u'HATILLO DE LOBA', 165: u'MAGANGU\xc3\u2030', 166: u'MAHATES', 167: u'MARGARITA', 168: u'MAR\xc3\x8dA LA BAJA', 169: u'MONTECRISTO', 170: u'MOMP\xc3\u201cS', 171: u'MORALES', 172: u'NOROS\xc3\x8d', 173: u'PINILLOS', 174: u'REGIDOR', 175: u'R\xc3\x8dO VIEJO', 176: u'SAN CRIST\xc3\u201cBAL', 177: u'SAN ESTANISLAO', 178: u'SAN FERNANDO', 179: u'SAN JACINTO', 180: u'SAN JACINTO DEL CAUCA', 181: u'SAN JUAN NEPOMUCENO', 182: u'SAN MART\xc3\x8dN DE LOBA', 183: u'SAN PABLO', 184: u'SANTA CATALINA', 185: u'SANTA ROSA', 186: u'SANTA ROSA DEL SUR', 187: u'SIMIT\xc3\x8d', 188: u'SOPLAVIENTO', 189: u'TALAIGUA NUEVO', 190: u'TIQUISIO', 191: u'TURBACO', 192: u'TURBAN\xc3\x81', 193: u'VILLANUEVA', 194: u'ZAMBRANO', 195: u'TUNJA', 196: u'ALMEIDA', 197: u'AQUITANIA', 198: u'ARCABUCO', 199: u'BEL\xc3\u2030N', 200: u'BERBEO', 201: u'BET\xc3\u2030ITIVA', 202: u'BOAVITA', 203: u'BOYAC\xc3\x81', 204: u'BRICE\xc3\u2018O', 205: u'BUENAVISTA', 206: u'BUSBANZ\xc3\x81', 207: u'CALDAS', 208: u'CAMPOHERMOSO', 209: u'CERINZA', 210: u'CHINAVITA', 211: u'CHIQUINQUIR\xc3\x81', 212: u'CHISCAS', 213: u'CHITA', 214: u'CHITARAQUE', 215: u'CHIVAT\xc3\x81', 216: u'CI\xc3\u2030NEGA', 217: u'C\xc3\u201cMBITA', 218: u'COPER', 219: u'CORRALES', 220: u'COVARACH\xc3\x8dA', 221: u'CUBAR\xc3\x81', 222: u'CUCAITA', 223: u'CU\xc3\x8dTIVA', 224: u'CH\xc3\x8dQUIZA', 225: u'CHIVOR', 226: u'DUITAMA', 227: u'EL COCUY', 228: u'EL ESPINO', 229: u'FIRAVITOBA', 230: u'FLORESTA', 231: u'GACHANTIV\xc3\x81', 232: u'G\xc3\x81MEZA', 233: u'GARAGOA', 234: u'GUACAMAYAS', 235: u'GUATEQUE', 236: u'GUAYAT\xc3\x81', 237: u'G\xc3\u0153IC\xc3\x81N', 238: u'IZA', 239: u'JENESANO', 240: u'JERIC\xc3\u201c', 241: u'LABRANZAGRANDE', 242: u'LA CAPILLA', 243: u'LA VICTORIA', 244: u'LA UVITA', 245: u'VILLA DE LEYVA', 246: u'MACANAL', 247: u'MARIP\xc3\x8d', 248: u'MIRAFLORES', 249: u'MONGUA', 250: u'MONGU\xc3\x8d', 251: u'MONIQUIR\xc3\x81', 252: u'MOTAVITA', 253: u'MUZO', 254: u'NOBSA', 255: u'NUEVO COL\xc3\u201cN', 256: u'OICAT\xc3\x81', 257: u'OTANCHE', 258: u'PACHAVITA', 259: u'P\xc3\x81EZ', 260: u'PAIPA', 261: u'PAJARITO', 262: u'PANQUEBA', 263: u'PAUNA', 264: u'PAYA', 265: u'PAZ DE R\xc3\x8dO', 266: u'PESCA', 267: u'PISBA', 268: u'PUERTO BOYAC\xc3\x81', 269: u'QU\xc3\x8dPAMA', 270: u'RAMIRIQU\xc3\x8d', 271: u'R\xc3\x81QUIRA', 272: u'ROND\xc3\u201cN', 273: u'SABOY\xc3\x81', 274: u'S\xc3\x81CHICA', 275: u'SAMAC\xc3\x81', 276: u'SAN EDUARDO', 277: u'SAN JOS\xc3\u2030 DE PARE', 278: u'SAN LUIS DE GACENO', 279: u'SAN MATEO', 280: u'SAN MIGUEL DE SEMA', 281: u'SAN PABLO DE BORBUR', 282: u'SANTANA', 283: u'SANTA MAR\xc3\x8dA', 284: u'SANTA ROSA DE VITERBO', 285: u'SANTA SOF\xc3\x8dA', 286: u'SATIVANORTE', 287: u'SATIVASUR', 288: u'SIACHOQUE', 289: u'SOAT\xc3\x81', 290: u'SOCOT\xc3\x81', 291: u'SOCHA', 292: u'SOGAMOSO', 293: u'SOMONDOCO', 294: u'SORA', 295: u'SOTAQUIR\xc3\x81', 296: u'SORAC\xc3\x81', 297: u'SUSAC\xc3\u201cN', 298: u'SUTAMARCH\xc3\x81N', 299: u'SUTATENZA', 300: u'TASCO', 301: u'TENZA', 302: u'TIBAN\xc3\x81', 303: u'TIBASOSA', 304: u'TINJAC\xc3\x81', 305: u'TIPACOQUE', 306: u'TOCA', 307: u'TOG\xc3\u0153\xc3\x8d', 308: u'T\xc3\u201cPAGA', 309: u'TOTA', 310: u'TUNUNGU\xc3\x81', 311: u'TURMEQU\xc3\u2030', 312: u'TUTA', 313: u'TUTAZ\xc3\x81', 314: u'\xc3\u0161MBITA', 315: u'VENTAQUEMADA', 316: u'VIRACACH\xc3\x81', 317: u'ZETAQUIRA', 318: u'MANIZALES', 319: u'AGUADAS', 320: u'ANSERMA', 321: u'ARANZAZU', 322: u'BELALC\xc3\x81ZAR', 323: u'CHINCHIN\xc3\x81', 324: u'FILADELFIA', 325: u'LA DORADA', 326: u'LA MERCED', 327: u'MANZANARES', 328: u'MARMATO', 329: u'MARQUETALIA', 330: u'MARULANDA', 331: u'NEIRA', 332: u'NORCASIA', 333: u'P\xc3\x81CORA', 334: u'PALESTINA', 335: u'PENSILVANIA', 336: u'RIOSUCIO', 337: u'RISARALDA', 338: u'SALAMINA', 339: u'SAMAN\xc3\x81', 340: u'SAN JOS\xc3\u2030', 341: u'SUP\xc3\x8dA', 342: u'VICTORIA', 343: u'VILLAMAR\xc3\x8dA', 344: u'VITERBO', 345: u'FLORENCIA', 346: u'ALBANIA', 347: u'BEL\xc3\u2030N DE LOS ANDAQU\xc3\x8dES', 348: u'CARTAGENA DEL CHAIR\xc3\x81', 349: u'CURILLO', 350: u'EL DONCELLO', 351: u'EL PAUJ\xc3\x8dL', 352: u'LA MONTA\xc3\u2018ITA', 353: u'MIL\xc3\x81N', 354: u'MORELIA', 355: u'PUERTO RICO', 356: u'SAN JOS\xc3\u2030 DEL FRAGUA', 357: u'SAN VICENTE DEL CAGU\xc3\x81N', 358: u'SOLANO', 359: u'SOLITA', 360: u'VALPARA\xc3\x8dSO', 361: u'POPAY\xc3\x81N', 362: u'ALMAGUER', 363: u'ARGELIA', 364: u'BALBOA', 365: u'BOL\xc3\x8dVAR', 366: u'BUENOS AIRES', 367: u'CAJIB\xc3\x8dO', 368: u'CALDONO', 369: u'CALOTO', 370: u'CORINTO', 371: u'EL TAMBO', 372: u'FLORENCIA', 373: u'GUACHEN\xc3\u2030', 374: u'GUAP\xc3\x8d', 375: u'INZ\xc3\x81', 376: u'JAMBAL\xc3\u201c', 377: u'LA SIERRA', 378: u'LA VEGA', 379: u'L\xc3\u201cPEZ DE MICAY', 380: u'MERCADERES', 381: u'MIRANDA', 382: u'MORALES', 383: u'PADILLA', 384: u'P\xc3\x81EZ', 385: u'PAT\xc3\x8dA', 386: u'PIAMONTE', 387: u'PIENDAM\xc3\u201c', 388: u'PUERTO TEJADA', 389: u'PURAC\xc3\u2030', 390: u'ROSAS', 391: u'SAN SEBASTI\xc3\x81N', 392: u'SANTANDER DE QUILICHAO', 393: u'SANTA ROSA', 394: u'SILVIA', 395: u'SOTARA', 396: u'SU\xc3\x81REZ', 397: u'SUCRE', 398: u'TIMB\xc3\x8dO', 399: u'TIMBIQU\xc3\x8d', 400: u'TORIB\xc3\x8dO', 401: u'TOTOR\xc3\u201c', 402: u'VILLA RICA', 403: u'VALLEDUPAR', 404: u'AGUACHICA', 405: u'AGUST\xc3\x8dN CODAZZI', 406: u'ASTREA', 407: u'BECERRIL', 408: u'BOSCONIA', 409: u'CHIMICHAGUA', 410: u'CHIRIGUAN\xc3\x81', 411: u'CURUMAN\xc3\x8d', 412: u'EL COPEY', 413: u'EL PASO', 414: u'GAMARRA', 415: u'GONZ\xc3\x81LEZ', 416: u'LA GLORIA', 417: u'LA JAGUA DE IBIRICO', 418: u'MANAURE BALC\xc3\u201cN DEL CESAR', 419: u'PAILITAS', 420: u'PELAYA', 421: u'PUEBLO BELLO', 422: u'R\xc3\x8dO DE ORO', 423: u'LA PAZ', 424: u'SAN ALBERTO', 425: u'SAN DIEGO', 426: u'SAN MART\xc3\x8dN', 427: u'TAMALAMEQUE', 428: u'MONTER\xc3\x8dA', 429: u'AYAPEL', 430: u'BUENAVISTA', 431: u'CANALETE', 432: u'CERET\xc3\u2030', 433: u'CHIM\xc3\x81', 434: u'CHIN\xc3\u0161', 435: u'CI\xc3\u2030NAGA DE ORO', 436: u'COTORRA', 437: u'LA APARTADA', 438: u'LORICA', 439: u'LOS C\xc3\u201cRDOBAS', 440: u'MOMIL', 441: u'MONTEL\xc3\x8dBANO', 442: u'MO\xc3\u2018ITOS', 443: u'PLANETA RICA', 444: u'PUEBLO NUEVO', 445: u'PUERTO ESCONDIDO', 446: u'PUERTO LIBERTADOR', 447: u'PUR\xc3\x8dSIMA DE LA CONCEPCI\xc3\u201cN', 448: u'SAHAG\xc3\u0161N', 449: u'SAN ANDR\xc3\u2030S DE SOTAVENTO', 450: u'SAN ANTERO', 451: u'SAN BERNARDO DEL VIENTO', 452: u'SAN CARLOS', 453: u'SAN JOS\xc3\u2030 DE UR\xc3\u2030', 454: u'SAN PELAYO', 455: u'TIERRALTA', 456: u'TUCH\xc3\x8dN'
    , 457: u'VALENCIA', 458: u'AGUA DE DIOS', 459: u'ALB\xc3\x81N', 460: u'ANAPOIMA', 461: u'ANOLAIMA', 462: u'ARBEL\xc3\x81EZ', 463: u'BELTR\xc3\x81N', 464: u'BITUIMA', 465: u'BOJAC\xc3\x81', 466: u'CABRERA', 467: u'CACHIPAY', 468: u'CAJIC\xc3\x81', 469: u'CAPARRAP\xc3\x8d', 470: u'C\xc3\x81QUEZA', 471: u'CARMEN DE CARUPA', 472: u'CHAGUAN\xc3\x8d', 473: u'CH\xc3\x8dA', 474: u'CHIPAQUE', 475: u'CHOACH\xc3\x8d', 476: u'CHOCONT\xc3\x81', 477: u'COGUA', 478: u'COTA', 479: u'CUCUNUB\xc3\x81', 480: u'EL COLEGIO', 481: u'EL PE\xc3\u2018\xc3\u201cN', 482: u'EL ROSAL', 483: u'FACATATIV\xc3\x81', 484: u'F\xc3\u201cMEQUE', 485: u'FOSCA', 486: u'FUNZA', 487: u'F\xc3\u0161QUENE', 488: u'FUSAGASUG\xc3\x81', 489: u'GACHAL\xc3\x81', 490: u'GACHANCIP\xc3\x81', 491: u'GACHET\xc3\x81', 492: u'GAMA', 493: u'GIRARDOT', 494: u'GRANADA', 495: u'GUACHET\xc3\x81', 496: u'GUADUAS', 497: u'GUASCA', 498: u'GUATAQU\xc3\x8d', 499: u'GUATAVITA', 500: u'GUAYABAL DE S\xc3\x8dQUIMA', 501: u'GUAYABETAL', 502: u'GUTI\xc3\u2030RREZ', 503: u'JERUSAL\xc3\u2030N', 504: u'JUN\xc3\x8dN', 505: u'LA CALERA', 506: u'LA MESA', 507: u'LA PALMA', 508: u'LA PE\xc3\u2018A', 509: u'LA VEGA', 510: u'LENGUAZAQUE', 511: u'MACHET\xc3\x81', 512: u'MADRID', 513: u'MANTA', 514: u'MEDINA', 515: u'MOSQUERA', 516: u'NARI\xc3\u2018O', 517: u'NEMOC\xc3\u201cN', 518: u'NILO', 519: u'NIMAIMA', 520: u'NOCAIMA', 521: u'VENECIA', 522: u'PACHO', 523: u'PAIME', 524: u'PANDI', 525: u'PARATEBUENO', 526: u'PASCA', 527: u'PUERTO SALGAR', 528: u'PUL\xc3\x8d', 529: u'QUEBRADANEGRA', 530: u'QUETAME', 531: u'QUIPILE', 532: u'APULO', 533: u'RICAURTE', 534: u'SAN ANTONIO DEL TEQUENDAMA', 535: u'SAN BERNARDO', 536: u'SAN CAYETANO', 537: u'SAN FRANCISCO', 538: u'SAN JUAN DE RIOSECO', 539: u'SASAIMA', 540: u'SESQUIL\xc3\u2030', 541: u'SIBAT\xc3\u2030', 542: u'SILVANIA', 543: u'SIMIJACA', 544: u'SOACHA', 545: u'SOP\xc3\u201c', 546: u'SUBACHOQUE', 547: u'SUESCA', 548: u'SUPAT\xc3\x81', 549: u'SUSA', 550: u'SUTATAUSA', 551: u'TABIO', 552: u'TAUSA', 553: u'TENA', 554: u'TENJO', 555: u'TIBACUY', 556: u'TIBIRITA', 557: u'TOCAIMA', 558: u'TOCANCIP\xc3\x81', 559: u'TOPAIP\xc3\x8d', 560: u'UBAL\xc3\x81', 561: u'UBAQUE', 562: u'VILLA DE SAN DIEGO DE UBAT\xc3\u2030', 563: u'UNE', 564: u'\xc3\u0161TICA', 565: u'VERGARA', 566: u'VIAN\xc3\x8d', 567: u'VILLAG\xc3\u201cMEZ', 568: u'VILLAPINZ\xc3\u201cN', 569: u'VILLETA', 570: u'VIOT\xc3\x81', 571: u'YACOP\xc3\x8d', 572: u'ZIPAC\xc3\u201cN', 573: u'ZIPAQUIR\xc3\x81', 574: u'QUIBD\xc3\u201c', 575: u'ACAND\xc3\x8d', 576: u'ALTO BAUD\xc3\u201c', 577: u'ATRATO', 578: u'BAGAD\xc3\u201c', 579: u'BAH\xc3\x8dA SOLANO', 580: u'BAJO BAUD\xc3\u201c', 581: u'BOJAY\xc3\x81', 582: u'EL CANT\xc3\u201cN DEL SAN PABLO', 583: u'CARMEN DEL DARI\xc3\u2030N', 584: u'C\xc3\u2030RTEGUI', 585: u'CONDOTO', 586: u'EL CARMEN DE ATRATO', 587: u'EL LITORAL DEL SAN JUAN', 588: u'ISTMINA', 589: u'JURAD\xc3\u201c', 590: u'LLOR\xc3\u201c', 591: u'MEDIO ATRATO', 592: u'MEDIO BAUD\xc3\u201c', 593: u'MEDIO SAN JUAN', 594: u'N\xc3\u201cVITA', 595: u'NUQU\xc3\x8d', 596: u'R\xc3\x8dO IR\xc3\u201c', 597: u'R\xc3\x8dO QUITO', 598: u'RIOSUCIO', 599: u'SAN JOS\xc3\u2030 DEL PALMAR', 600: u'SIP\xc3\x8d', 601: u'TAD\xc3\u201c', 602: u'UNGU\xc3\x8dA', 603: u'UNI\xc3\u201cN PANAMERICANA', 604: u'NEIVA', 605: u'ACEVEDO', 606: u'AGRADO', 607: u'AIPE', 608: u'ALGECIRAS', 609: u'ALTAMIRA', 610: u'BARAYA', 611: u'CAMPOALEGRE', 612: u'COLOMBIA', 613: u'EL\xc3\x8dAS', 614: u'GARZ\xc3\u201cN', 615: u'GIGANTE', 616: u'GUADALUPE', 617: u'HOBO', 618: u'\xc3\x8dQUIRA', 619: u'ISNOS', 620: u'LA ARGENTINA', 621: u'LA PLATA', 622: u'N\xc3\x81TAGA', 623: u'OPORAPA', 624: u'PAICOL', 625: u'PALERMO', 626: u'PALESTINA', 627: u'PITAL', 628: u'PITALITO', 629: u'RIVERA', 630: u'SALADOBLANCO', 631: u'SAN AGUST\xc3\x8dN', 632: u'SANTA MAR\xc3\x8dA', 633: u'SUAZA', 634: u'TARQUI', 635: u'TESALIA', 636: u'TELLO', 637: u'TERUEL', 638: u'TIMAN\xc3\x81', 639: u'VILLAVIEJA', 640: u'YAGUAR\xc3\x81', 641: u'RIOHACHA', 642: u'ALBANIA', 643: u'BARRANCAS', 644: u'DIBULLA', 645: u'DISTRACCI\xc3\u201cN', 646: u'EL MOLINO', 647: u'FONSECA', 648: u'HATONUEVO', 649: u'LA JAGUA DEL PILAR', 650: u'MAICAO', 651: u'MANAURE', 652: u'SAN JUAN DEL CESAR', 653: u'URIBIA', 654: u'URUMITA', 655: u'VILLANUEVA', 656: u'SANTA MARTA', 657: u'ALGARROBO', 658: u'ARACATACA', 659: u'ARIGUAN\xc3\x8d', 660: u'CERRO DE SAN ANTONIO', 661: u'CHIVOLO', 662: u'CI\xc3\u2030NAGA', 663: u'CONCORDIA', 664: u'EL BANCO', 665: u'EL PI\xc3\u2018\xc3\u201cN', 666: u'EL RET\xc3\u2030N', 667: u'FUNDACI\xc3\u201cN', 668: u'GUAMAL', 669: u'NUEVA GRANADA', 670: u'PEDRAZA', 671: u'PIJI\xc3\u2018O DEL CARMEN', 672: u'PIVIJAY', 673: u'PLATO', 674: u'PUEBLOVIEJO', 675: u'REMOLINO', 676: u'SABANAS DE SAN \xc3\x81NGEL', 677: u'SALAMINA', 678: u'SAN SEBASTI\xc3\x81N DE BUENAVISTA', 679: u'SAN ZEN\xc3\u201cN', 680: u'SANTA ANA', 681: u'SANTA B\xc3\x81RBARA DE PINTO', 682: u'SITIONUEVO', 683: u'TENERIFE', 684: u'ZAPAY\xc3\x81N', 685: u'ZONA BANANERA', 686: u'VILLAVICENCIO', 687: u'ACAC\xc3\x8dAS', 688: u'BARRANCA DE UP\xc3\x8dA', 689: u'CABUYARO', 690: u'CASTILLA LA NUEVA', 691: u'SAN LUIS DE CUBARRAL', 692: u'CUMARAL', 693: u'EL CALVARIO', 694: u'EL CASTILLO', 695: u'EL DORADO', 696: u'FUENTE DE ORO', 697: u'GRANADA', 698: u'GUAMAL', 699: u'MAPIRIP\xc3\x81N', 700: u'MESETAS', 701: u'LA MACARENA', 702: u'URIBE', 703: u'LEJAN\xc3\x8dAS', 704: u'PUERTO CONCORDIA', 705: u'PUERTO GAIT\xc3\x81N', 706: u'PUERTO L\xc3\u201cPEZ', 707: u'PUERTO LLERAS', 708: u'PUERTO RICO', 709: u'RESTREPO', 710: u'SAN CARLOS DE GUAROA', 711: u'SAN JUAN DE ARAMA', 712: u'SAN JUANITO', 713: u'SAN MART\xc3\x8dN', 714: u'VISTAHERMOSA', 715: u'PASTO', 716: u'ALB\xc3\x81N', 717: u'ALDANA', 718: u'ANCUY\xc3\x81', 719: u'ARBOLEDA', 720: u'BARBACOAS', 721: u'BEL\xc3\u2030N', 722: u'BUESACO', 723: u'COL\xc3\u201cN', 724: u'CONSAC\xc3\x81', 725: u'CONTADERO', 726: u'C\xc3\u201cRDOBA', 727: u'CUASP\xc3\u0161D', 728: u'CUMBAL', 729: u'CUMBITARA', 730: u'CHACHAG\xc3\u0153\xc3\x8d', 731: u'EL CHARCO', 732: u'EL PE\xc3\u2018OL', 733: u'EL ROSARIO', 734: u'EL TABL\xc3\u201cN DE G\xc3\u201cMEZ', 735: u'EL TAMBO', 736: u'FUNES', 737: u'GUACHUCAL', 738: u'GUAITARILLA', 739: u'GUALMAT\xc3\x81N', 740: u'ILES', 741: u'IMU\xc3\u2030S', 742: u'IPIALES', 743: u'LA CRUZ', 744: u'LA FLORIDA', 745: u'LA LLANADA', 746: u'LA TOLA', 747: u'LA UNI\xc3\u201cN', 748: u'LEIVA', 749: u'LINARES', 750: u'LOS ANDES', 751: u'MAG\xc3\u0153\xc3\x8d', 752: u'MALLAMA', 753: u'MOSQUERA', 754: u'NARI\xc3\u2018O', 755: u'OLAYA HERRERA', 756: u'OSPINA', 757: u'FRANCISCO PIZARRO', 758: u'POLICARPA', 759: u'POTOS\xc3\x8d', 760: u'PROVIDENCIA', 761: u'PUERRES', 762: u'PUPIALES', 763: u'RICAURTE', 764: u'ROBERTO PAY\xc3\x81N', 765: u'SAMANIEGO', 766: u'SANDON\xc3\x81', 767: u'SAN BERNARDO', 768: u'SAN LORENZO', 769: u'SAN PABLO', 770: u'SAN PEDRO DE CARTAGO', 771: u'SANTA B\xc3\x81RBARA', 772: u'SANTACRUZ', 773: u'SAPUYES', 774: u'TAMINANGO', 775: u'TANGUA', 776: u'SAN ANDR\xc3\u2030S DE TUMACO', 777: u'T\xc3\u0161QUERRES', 778: u'YACUANQUER', 779: u'C\xc3\u0161CUTA', 780: u'\xc3\x81BREGO', 781: u'ARBOLEDAS', 782: u'BOCHALEMA', 783: u'BUCARASICA', 784: u'C\xc3\x81COTA', 785: u'C\xc3\x81CHIRA', 786: u'CHIN\xc3\x81COTA', 787: u'CHITAG\xc3\x81', 788: u'CONVENCI\xc3\u201cN', 789: u'CUCUTILLA', 790: u'DURANIA', 791: u'EL CARMEN', 792: u'EL TARRA', 793: u'EL ZULIA', 794: u'GRAMALOTE', 795: u'HACAR\xc3\x8d', 796: u'HERR\xc3\x81N', 797: u'LABATECA', 798: u'LA ESPERANZA', 799: u'LA PLAYA', 800: u'LOS PATIOS', 801: u'LOURDES', 802: u'MUTISCUA', 803: u'OCA\xc3\u2018A', 804: u'PAMPLONA', 805: u'PAMPLONITA', 806: u'PUERTO SANTANDER', 807: u'RAGONVALIA', 808: u'SALAZAR', 809: u'SAN CALIXTO', 810: u'SAN CAYETANO', 811: u'SANTIAGO', 812: u'SARDINATA', 813: u'SILOS', 814: u'TEORAMA', 815: u'TIB\xc3\u0161', 816: u'TOLEDO', 817: u'VILLA CARO', 818: u'VILLA DEL ROSARIO', 819: u'ARMENIA', 820: u'BUENAVISTA', 821: u'CALARC\xc3\x81', 822: u'CIRCASIA', 823: u'C\xc3\u201cRDOBA', 824: u'FILANDIA', 825: u'G\xc3\u2030NOVA', 826: u'LA TEBAIDA', 827: u'MONTENEGRO', 828: u'PIJAO', 829: u'QUIMBAYA', 830: u'SALENTO', 831: u'PEREIRA', 832: u'AP\xc3\x8dA', 833: u'BALBOA', 834: u'BEL\xc3\u2030N DE UMBR\xc3\x8dA', 835: u'DOSQUEBRADAS', 836: u'GU\xc3\x81TICA', 837: u'LA CELIA', 838: u'LA VIRGINIA', 839: u'MARSELLA', 840: u'MISTRAT\xc3\u201c', 841: u'PUEBLO RICO', 842: u'QUINCH\xc3\x8dA', 843: u'SANTA ROSA DE CABAL', 844: u'SANTUARIO', 845: u'BUCARAMANGA', 846: u'AGUADA', 847: u'ALBANIA', 848: u'ARATOCA', 849: u'BARBOSA', 850: u'BARICHARA', 851: u'BARRANCABERMEJA', 852: u'BETULIA', 853: u'BOL\xc3\x8dVAR', 854: u'CABRERA', 855: u'CALIFORNIA', 856: u'CAPITANEJO', 857: u'CARCAS\xc3\x8d', 858: u'CEPIT\xc3\x81', 859: u'CERRITO', 860: u'CHARAL\xc3\x81', 861: u'CHARTA', 862: u'CHIMA', 863: u'CHIPAT\xc3\x81', 864: u'CIMITARRA', 865: u'CONCEPCI\xc3\u201cN', 866: u'CONFINES', 867: u'CONTRATACI\xc3\u201cN', 868: u'COROMORO', 869: u'CURIT\xc3\x8d', 870: u'EL CARMEN DE CHUCUR\xc3\x8d', 871: u'EL GUACAMAYO', 872: u'EL PE\xc3\u2018\xc3\u201cN', 873: u'EL PLAY\xc3\u201cN', 874: u'ENCINO', 875: u'ENCISO', 876: u'FLORI\xc3\x81N', 877: u'FLORIDABLANCA', 878: u'GAL\xc3\x81N', 879: u'G\xc3\x81MBITA', 880: u'GIR\xc3\u201cN', 881: u'GUACA', 882: u'GUADALUPE', 883: u'GUAPOT\xc3\x81', 884: u'GUAVAT\xc3\x81', 885: u'G\xc3\u0153EPSA', 886: u'HATO', 887: u'JES\xc3\u0161S MAR\xc3\x8dA', 888: u'JORD\xc3\x81N', 889: u'LA BELLEZA', 890: u'LAND\xc3\x81ZURI', 891: u'LA PAZ', 892: u'LEBRIJA', 893: u'LOS SANTOS', 894: u'MACARAVITA', 895: u'M\xc3\x81LAGA', 896: u'MATANZA', 897: u'MOGOTES', 898: u'MOLAGAVITA', 899: u'OCAMONTE', 900: u'OIBA', 901: u'ONZAGA', 902: u'PALMAR', 903: u'PALMAS DEL SOCORRO', 904: u'P\xc3\x81RAMO', 905: u'PIEDECUESTA', 906: u'PINCHOTE', 907: u'PUENTE NACIONAL', 908: u'PUERTO PARRA', 909: u'PUERTO WILCHES', 910: u'RIONEGRO', 911: u'SABANA DE TORRES', 912: u'SAN ANDR\xc3\u2030S', 913: u'SAN BENITO', 914: u'SAN GIL', 915: u'SAN JOAQU\xc3\x8dN', 916: u'SAN JOS\xc3\u2030 DE MIRANDA', 917: u'SAN MIGUEL', 918: u'SAN VICENTE DE CHUCUR\xc3\x8d',
    919: u'SANTA B\xc3\x81RBARA', 920: u'SANTA HELENA DEL OP\xc3\u201cN', 921: u'SIMACOTA', 922: u'SOCORRO', 923: u'SUAITA', 924: u'SUCRE', 925: u'SURAT\xc3\x81', 926: u'TONA', 927: u'VALLE DE SAN JOS\xc3\u2030', 928: u'V\xc3\u2030LEZ', 929: u'VETAS', 930: u'VILLANUEVA', 931: u'ZAPATOCA', 932: u'SINCELEJO', 933: u'BUENAVISTA', 934: u'CAIMITO', 935: u'COLOSO', 936: u'COROZAL', 937: u'COVE\xc3\u2018AS', 938: u'CHAL\xc3\x81N', 939: u'EL ROBLE', 940: u'GALERAS', 941: u'GUARANDA', 942: u'LA UNI\xc3\u201cN', 943: u'LOS PALMITOS', 944: u'MAJAGUAL', 945: u'MORROA', 946: u'OVEJAS', 947: u'PALMITO', 948: u'SAMPU\xc3\u2030S', 949: u'SAN BENITO ABAD', 950: u'SAN JUAN DE BETULIA', 951: u'SAN MARCOS', 952: u'SAN ONOFRE', 953: u'SAN PEDRO', 954: u'SAN LUIS DE SINC\xc3\u2030', 955: u'SUCRE', 956: u'SANTIAGO DE TOL\xc3\u0161', 957: u'TOL\xc3\u0161 VIEJO', 958: u'IBAGU\xc3\u2030', 959: u'ALPUJARRA', 960: u'ALVARADO', 961: u'AMBALEMA', 962: u'ANZO\xc3\x81TEGUI', 963: u'ARMERO GUAYABAL', 964: u'ATACO', 965: u'CAJAMARCA', 966: u'CARMEN DE APICAL\xc3\x81', 967: u'CASABIANCA', 968: u'CHAPARRAL', 969: u'COELLO', 970: u'COYAIMA', 971: u'CUNDAY', 972: u'DOLORES', 973: u'ESPINAL', 974: u'FALAN', 975: u'FLANDES', 976: u'FRESNO', 977: u'GUAMO', 978: u'HERVEO', 979: u'HONDA', 980: u'ICONONZO', 981: u'L\xc3\u2030RIDA', 982: u'L\xc3\x8dBANO', 983: u'SAN SEBASTI\xc3\x81N DE MARIQUITA', 984: u'MELGAR', 985: u'MURILLO', 986: u'NATAGAIMA', 987: u'ORTEGA', 988: u'PALOCABILDO', 989: u'PIEDRAS', 990: u'PLANADAS', 991: u'PRADO', 992: u'PURIFICACI\xc3\u201cN', 993: u'RIOBLANCO', 994: u'RONCESVALLES', 995: u'ROVIRA', 996: u'SALDA\xc3\u2018A', 997: u'SAN ANTONIO', 998: u'SAN LUIS', 999: u'SANTA ISABEL', 1000: u'SU\xc3\x81REZ', 1001: u'VALLE DE SAN JUAN', 1002: u'VENADILLO', 1003: u'VILLAHERMOSA', 1004: u'VILLARRICA', 1005: u'CALI', 1006: u'ALCAL\xc3\x81', 1007: u'ANDALUC\xc3\x8dA', 1008: u'ANSERMANUEVO', 1009: u'ARGELIA', 1010: u'BOL\xc3\x8dVAR', 1011: u'BUENAVENTURA', 1012: u'GUADALAJARA DE BUGA', 1013: u'BUGALAGRANDE', 1014: u'CAICEDONIA', 1015: u'CALIMA', 1016: u'CANDELARIA', 1017: u'CARTAGO', 1018: u'DAGUA', 1019: u'EL \xc3\x81GUILA', 1020: u'EL CAIRO', 1021: u'EL CERRITO', 1022: u'EL DOVIO', 1023: u'FLORIDA', 1024: u'GINEBRA', 1025: u'GUACAR\xc3\x8d', 1026: u'JAMUND\xc3\x8d', 1027: u'LA CUMBRE', 1028: u'LA UNI\xc3\u201cN', 1029: u'LA VICTORIA', 1030: u'OBANDO', 1031: u'PALMIRA', 1032: u'PRADERA', 1033: u'RESTREPO', 1034: u'RIOFR\xc3\x8dO', 1035: u'ROLDANILLO', 1036: u'SAN PEDRO', 1037: u'SEVILLA', 1038: u'TORO', 1039: u'TRUJILLO', 1040: u'ULLOA', 1041: u'VERSALLES', 1042: u'VIJES', 1043: u'YOTOCO', 1044: u'YUMBO', 1045: u'ZARZAL', 1046: u'ARAUCA', 1047: u'ARAUQUITA', 1048: u'CRAVO NORTE', 1049: u'FORTUL', 1050: u'PUERTO ROND\xc3\u201cN', 1051: u'SARAVENA', 1052: u'TAME', 1053: u'YOPAL', 1054: u'AGUAZUL', 1055: u'CH\xc3\x81MEZA', 1056: u'HATO COROZAL', 1057: u'LA SALINA', 1058: u'MAN\xc3\x8d', 1059: u'MONTERREY', 1060: u'NUNCH\xc3\x8dA', 1061: u'OROCU\xc3\u2030', 1062: u'PAZ DE ARIPORO', 1063: u'PORE', 1064: u'RECETOR', 1065: u'SABANALARGA', 1066: u'S\xc3\x81CAMA', 1067: u'SAN LUIS DE PALENQUE', 1068: u'T\xc3\x81MARA', 1069: u'TAURAMENA', 1070: u'TRINIDAD', 1071: u'VILLANUEVA', 1072: u'MOCOA', 1073: u'COL\xc3\u201cN', 1074: u'ORITO', 1075: u'PUERTO AS\xc3\x8dS', 1076: u'PUERTO CAICEDO', 1077: u'PUERTO GUZM\xc3\x81N', 1078: u'PUERTO LEGU\xc3\x8dZAMO', 1079: u'SIBUNDOY', 1080: u'SAN FRANCISCO', 1081: u'SAN MIGUEL', 1082: u'SANTIAGO', 1083: u'VALLE DEL GUAMUEZ', 1084: u'VILLAGARZ\xc3\u201cN', 1085: u'LETICIA', 1086: u'EL ENCANTO', 1087: u'LA CHORRERA', 1088: u'LA PEDRERA', 1089: u'LA VICTORIA', 1090: u'MIRIT\xc3\x8d - PARAN\xc3\x81', 1091: u'PUERTO ALEGR\xc3\x8dA', 1092: u'PUERTO ARICA', 1093: u'PUERTO NARI\xc3\u2018O', 1094: u'PUERTO SANTANDER', 1095: u'TARAPAC\xc3\x81', 1096: u'IN\xc3\x8dRIDA', 1097: u'BARRANCO MINAS', 1098: u'MAPIRIPANA', 1099: u'SAN FELIPE', 1100: u'PUERTO COLOMBIA', 1101: u'LA GUADALUPE', 1102: u'CACAHUAL', 1103: u'PANA PANA', 1104: u'MORICHAL', 1105: u'SAN JOS\xc3\u2030 DEL GUAVIARE', 1106: u'CALAMAR', 1107: u'EL RETORNO', 1108: u'MIRAFLORES', 1109: u'MIT\xc3\u0161', 1110: u'CARUR\xc3\u0161', 1111: u'PACOA', 1112: u'TARAIRA', 1113: u'PAPUNAUA', 1114: u'YAVARAT\xc3\u2030', 1115: u'PUERTO CARRE\xc3\u2018O', 1116: u'LA PRIMAVERA', 1117: u'SANTA ROSAL\xc3\x8dA', 1118: u'CUMARIBO', 1119: u'TULUA'}

    print('----------------')
    print("FINAL SOLUTION: ", Sol, "\nFINAL OF: ", Of)
    output = { "objectiveFunction": Of,
               "runningTime": time2,
               "algorithm": "p-Exhaustive-Functional-Regions",
               "regions": len(Sol),
               "seeds" : seeds,
               "r2a": Sol,
               "distanceType": distanceType,
               "distanceStat": distanceStat,
               "selectionType": selectionType,
               "ObjectiveFuncionType": objectiveFunctionType}
    
    # instance = Clusterpy.importArcData("Entradas/Oficial_Colombia/Municipios_RIMISP")
    # '''
    # if len(rm.area2Region)> 50:
    #     instance = clusterpy.importArcData("Entradas/Antioquia_2")
    # elif len(rm.area2Region)> 130:
    #     instance = clusterpy.importArcData("Entradas/mpio")
    # else:
    #     instance = clusterpy.createGrid(rows, cols)
    # '''
    # NDs = rm.seeds#[5,8,17] #nodos dinamizadores
    # regions = rm.returnRegions() #[2,2,2,0,0,2,2,2,0,0,2,2,1,0,0,1,1,1,1,0,1,1,1,1,1] #solución
    # OF = rm.getObj(Dij, Cio, rm.seeds)#4300 #Función objetivo
    # it = r #número de la iteración que se grafica (este número se usa para nombrar la figura que se exporta al folder outputMaps)

    return output

