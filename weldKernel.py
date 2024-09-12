# --*coding=UTF-8*--
from abaqus import *
from abaqusConstants import *
from weldModule import WeldSequence
from fileModule import FileModule
   
   
def weldKernel(weldParameters):
    """
    #    
    """  
    ws = WeldSequence(weldParameters)  
    fm = FileModule(ws)   
    
    
# Test	
if __name__ == '__main__':	
    session.journalOptions.setValues(replayGeometry=INDEX,recoverGeometry=INDEX)
    
    # --------------------------------------------------------------------------------     
    #p1 = mdb.models['tee'].parts['tee']
    #v1, n1, f1, e1 = p1.vertices, p1.nodes, p1.faces, p1.edges   
    #
    #p2 = mdb.models['sweep'].parts['sweep']
    #v2, n2, f2, e2 = p2.vertices, p2.nodes, p2.faces, p2.edges     
    #
    #weldParameters = [{'modelName':'tee',    'partName':'tee',    'crossFace':   (f1[32],),    'curveType':'ArbitraryCurve', 'startPoint':n1[2271],'alongPoint': n1[2251],'toePoint':n1[2229], 'segNumber': 4, 'weldVoltage': 10.2,  'weldCurrent': 300.0,'weldVelocty': 5.0,  'weldYita': 0.8, 'coolTime':2.0, 'deltTemp': 1500,'heatType':'DoubleEllipsoid', 'heatPara':(2.0, 4.0, 0.9, 3.0) },
    #                  {'modelName':'sweep',  'partName':'sweep',  'crossFace':(f2[22],f2[23]), 'curveType':'Circle'        , 'startPoint':n2[0],   'alongPoint': n2[96],  'toePoint':n2[6],    'segNumber': 3, 'weldVoltage': 10.2,  'weldCurrent': 300.0,'weldVelocty': 5.0,  'weldYita': 0.8, 'coolTime':2.0, 'deltTemp': 1500,'heatType':'DoubleEllipsoid', 'heatPara':(3.0, 4.0, 0.9, 3.0) },
    #                  {'modelName':'sweep',  'partName':'sweep',  'crossFace':(f2[9],f2[10]),  'curveType':'StraightLine'  , 'startPoint':n2[26],  'alongPoint': n2[407], 'toePoint':n2[61],   'segNumber': 1, 'weldVoltage': 10.2,  'weldCurrent': 300.0,'weldVelocty': 5.0,  'weldYita': 0.8, 'coolTime':2.0, 'deltTemp': 1500,'heatType':'DoubleEllipsoid', 'heatPara':(9.0, 4.0, 0.9, 3.0) }]
    
    p1 = mdb.models['tee'].parts['tee']
    v1, n1, f1, e1 = p1.vertices, p1.nodes, p1.faces, p1.edges   

    p2 = mdb.models['segn'].parts['Part-1']
    v2, n2, f2, e2 = p2.vertices, p2.nodes, p2.faces, p2.edges     
    
    weldParameters = [{'modelName':'tee',    'partName':'tee',    'crossFace':(f1[32],), 'curveType':'ArbitraryCurve', 'startPoint':n1[7],   'alongPoint': n1[46], 'toePoint':n1[6],   'segNumber': 16, 'weldVoltage': 9.5,  'weldCurrent': 300.0,'weldVelocty': 3.9,   'weldYita': 0.9, 'coolTime':2.0, 'deltTemp': 1500,'heatType':'DoubleEllipsoid', 'heatPara':(1.62, 3.24, 1.62, 0.77) },
                      {'modelName':'segn',  'partName':'Part-1',  'crossFace':(f2[28],), 'curveType':'Circle'        , 'startPoint':n2[0],   'alongPoint': n2[42], 'toePoint':n2[1],   'segNumber': 27, 'weldVoltage': 9.5,  'weldCurrent': 300.0,'weldVelocty': 10.0,  'weldYita': 0.9, 'coolTime':0.0, 'deltTemp': 1500,'heatType':'DoubleEllipsoid', 'heatPara':(3.2, 6.4, 3.2, 3.2) },
                      {'modelName':'segn',  'partName':'Part-1',  'crossFace':(f2[5],),  'curveType':'StraightLine'  , 'startPoint':n2[582], 'alongPoint': n2[392],'toePoint':n2[570], 'segNumber': 25, 'weldVoltage': 9.5,  'weldCurrent': 300.0,'weldVelocty': 10.0,  'weldYita': 0.9, 'coolTime':5.0, 'deltTemp': 1500,'heatType':'DoubleEllipsoid', 'heatPara':(3.2, 6.4, 3.2, 3.2) }]
   
 
    weldKernel(weldParameters)  

    
p = mdb.models['tlp20-3tobecontinue'].parts['Part-1']
v = p.queryGeometry(printResults=False)['volume']
l = v**(1.0/3)
e = p.edges
edge = e[5]
length = edge.getSize(printResults=False)
if length > 1:
    num = int(length)*10
else:
    num = 10    
oldFeature = list(p.features.keys())
oldDatums = list(p.datums.keys())
coords = []    
for i in range(num):
    p.DatumPointByEdgeParam(edge=edge, parameter=float(i)/(num-1))
    
    newFeature = list(p.features.keys())
    keys = [f for f in newFeature if f not in oldFeature]
    
    newDatums = list(p.datums.keys())
    indx = [d for d in newDatums if d not in oldDatums]
    
    coords.append(p.datums[indx[0]].pointOn)
    del p.features[keys[0]]
    
def getCurveType(edge):  
    try:
        ra = edge.getRadius() 
        curveType = 'Circle'
    except:
        try:
            ra = getCurvature(parameter=0.0)['radius'] 
            curveType = 'ArbitraryCurve'
        except:
            curveType = 'StraightLine'
            
    return curveType        