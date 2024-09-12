# --*coding=UTF-8*--
from abaqus import *
from abaqusConstants import *

def getCurveType(edge):  
    try:
        ra = edge.getRadius() 
        curveType = 'CIRCLE'
    except:
        try:
            ra = edge.getCurvature(parameter=0.0)['radius'] 
            curveType = 'ARBITRARY_CURVE'
        except:
            curveType = 'STRAIGHT_LINE'
            
    return curveType 
    
# Test	
if __name__ == '__main__':	
    session.journalOptions.setValues(replayGeometry=INDEX,recoverGeometry=INDEX)    
    p = mdb.models['tee'].parts['tee']    
    e = p.edges
    edge1 = e[9]  # StraightLine
    highlight(edge1)
    edge2 = e[34]  # Circle
    edge3 = e[36]  # ArbitraryCurve 
    t1 = getCurveType(edge1)
    t2 = getCurveType(edge2)
    t3 = getCurveType(edge3)
    
    a = mdb.models['tee'].rootAssembly
    e1 = a.instances['tee-1'].edges
    edge4 = e1[9]
    edge5 = e1[34]  
    edge6 = e1[36]
    t4 = getCurveType(edge4)
    t5 = getCurveType(edge5)
    t6 = getCurveType(edge6)  

    r1 = e1.getClosest(coordinates=((15.0,22.0,0.0),)) 
    r2 = e1.getClosest(coordinates=((14.369817,22.498789,-4.302135),)) 