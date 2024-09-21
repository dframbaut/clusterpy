from contiguity import weightsFromAreas
from mrpolygons import mrpolygon, scalePolygon, polarPolygon2cartesian, transportPolygonGeometry, transportPolygon

from Polygon import Polygon

class rimap():
    def __init__(self,n=3600,N=30,alpha=[0.1,0.5],sigma=[1.1,1.4],dt=0.1,pg=0.01,pu=0.05,su=0.315,boundary=""):
        """Creates an irregular maps

        :param n: number of areas 
        :type n: integer
        :param N: number of points sampled from each irregular polygon (MR-Polygon) 
        :type N: integer
        :param alpha: min and max value to sampled alpha; default is (0.1,0.5)
        :type alpha: List
        :param sigma: min and max value to sampled sigma; default is (1.2,1.5)
        :type sigma: List
        :param dt: time delta to be used to create irregular polygons (MR-Polygons)
        :type dt: Float
        :param pg: parameter to define the scaling factor of each polygon before being introduced as part of the irregular map
        :type pg: Float
        :param pu: parameter to define the probability of increase the number of areas of each polygon before being introduced into the irregular map
        :type pu: Float
        :param su: parameter to define how much is increased the number of areas of each polygon before being introduced into the irregular map
        :type su: Float
        :param boundary: Initial irregular boundary to be used into the recursive irregular map algorithm
        :type boundary: Layer

        :rtype: Layer
        :return: RI-Map instance 
        """
        self.carteAreas = []
        self.carteExternal = []
        self.n = n
        self.pg = pg
        self.pu = pu
        self.su = su
        # Initializing area parameters
        self.alpha = alpha
        self.sigma = sigma
        self.N = N + 1
        self.mu = 10
        self.X_0 = 10
        self.dt = dt
        alp = 0.4
        sig = 1.2
        self.lAreas = 0
        if boundary == "":
            a,r,sa,sr,X1,times = mrpolygon(alp,sig,self.mu,self.X_0,self.dt,self.N)
            sa,sr = scalePolygon(sa,sr,1000)
            polygon = polarPolygon2cartesian(zip(sa,sr))
        else:
            layer = boundary
            polygon = layer.areas[0][0]
        polygon = Polygon(polygon)
        self.areasPerLevel = {}
        self.oPolygon = polygon
        areas, coveredArea = self.dividePolygon(polygon,Polygon(),self.n,0.97)
        areaUnion = Polygon()
        if len(areas) > n:
            areas = self.postCorrectionDissolve(areas,n)
        for a in areas:
            a = a[0]
            if a[-1] != a[0]:
                a.append(a[0])
            self.carteAreas.append([a])
        print ("closing: " + str(len(self.carteAreas)))

    def postCorrectionDissolve(self,areas,nAreas):
        def deleteAreaFromW(areaId,newId,W):
            neighs = W[areaId]
            W.pop(areaId)
            for n in neighs:
                W[n].remove(areaId)
                if n != newId:
                    W[n].append(newId)
                    W[n] = list(set(W[n]))
            W[newId].extend(neighs)
            W[newId] = list(set(W[newId]))
            W[newId].remove(newId)
            return W
        pos = 0
        Wrook, Wqueen = weightsFromAreas(areas)
        aIds = filter(lambda x: len(Wrook[x])>0,Wrook.keys())
        aIds0 = filter(lambda x: len(Wrook[x])==0,Wrook.keys())
        areas = [areas[x] for x in aIds]
        areas.sort(key = lambda x: len(x[0]))
        Wrook, Wqueen = weightsFromAreas(areas)
        availableAreas = Wrook.keys()
        end = False
        pos = availableAreas.pop(0)
        id2pos = Wrook.keys()
        while len(areas) > nAreas and not end:
            area = areas[id2pos.index(pos)]
            if len(Wrook[pos]) > 0:
                neighs = Wrook[pos]
                neighs.sort(key=lambda x: areas[id2pos.index(x)].area(),reverse=True)
                for k,nneigh in enumerate(neighs):
                    neigh = areas[id2pos.index(nneigh)]
                    narea = area | neigh
                    if len(narea) == 1:
                        areas[id2pos.index(nneigh)] = narea
                        areas.pop(id2pos.index(pos))
                        id2pos.remove(pos)
                        Wrook = deleteAreaFromW(pos,nneigh,Wrook) # El error puede estar aca
                        break
                    else:
                        pass
                if len(areas) == nAreas:
                    end = True
                    break
                else:
                    if len(availableAreas) > 0:
                        pos = availableAreas.pop(0)
                    else:
                        end = True
            else:
                Wrook.pop(pos)
                area = areas.pop(id2pos.index(pos))
                id2pos.remove(pos)
                if len(availableAreas) > 0:
                    pos = availableAreas.pop(0)
                else:
                    end = True
        return areas