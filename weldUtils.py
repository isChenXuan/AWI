# --*coding=UTF-8*--
from abaqus import *
from abaqusConstants import *
import numpy as np
import numpy.linalg as nla
import re

 
# -------------------------------------------------------------------------------------------------------------------------------------------
def getNodeFromPoint(point): 
    """
    # Get element node associated with input point. Note: coordinates of the same point may be different between abaqus part and assembly.
    # Inputs: 
    #        point : Abaqus part point object.
    # Outputs: 
    #        pNode : Part node object associated with input point. If no node found, return None.
    """   
    if hasattr(point, 'pointOn'):
        # Input point is of type geometry Vertex.
        pNode = point.getNodes()[0]
    elif hasattr(point, 'coordinates'): 
        # Input point is of type element node.
        pNode = point
    else:
        return
     
    return pNode    


# -------------------------------------------------------------------------------------------------------------------------------------------
def getElemFacesFromFace(face): 
    """
    # Get all element faces associated with input face and save them in elemFaces.
    # Inputs: 
    #        face: Tuple or list of abaqus face object input by user, can consist of multi element faces or geometry faces.
    # Outputs: 
    #        elemFaces: List of element faces found with input faceObject. If no element face found, return None.
    """
    elemFaces = []    
    for obj in face:
        if hasattr(obj, 'getElementsViaTopology'):
            # Input face is of type element face.
            elemFaces.append(obj)
        elif hasattr(obj, 'getAdjacentFaces'):
            # Input face is of type geometry face.
            elemFaces = elemFaces + list(obj.getElementFaces())
        else:
            return        
     
    return elemFaces   
    
    
# -------------------------------------------------------------------------------------------------------------------------------------------
def getAdjacentNodes(node):
    """
    # Get adjacent nodes of input node. As shown below.
    #                                 ×----×---1----×----×
    #                                /|       /|        /|  
    #                               ×-|--2---0-|--4----× |
    #                               | |      | |       | |  
    #                               | ×---×--|-×----×--|-×
    #                               |/       |/        |/  
    #                               ×---×----3----×----×
    # Inputs: 
    #        node: Node object. Node 0 in above figure.
    # Outputs: 
    #        adjaNodes: Adjacent nodes of input node. Nodes 1,2,3 and 4 in above figure. 
    """
    adjaNodes = []
    elemEdges = list(node.getElemEdges())
    if len(elemEdges) == 1:  # Input node is on element edge. 
        nodes = list(elemEdges[0].getNodes())
        nodes.remove(node)
        adjaNodes = nodes
    else:                    # Input node is on element vertex.
        for edge in elemEdges:
            nodes = list(edge.getNodes())
            nodes.remove(node)
            if len(nodes) == 1:
                adjaNodes.append(nodes[0])
            else:
                for item in nodes:
                    if len(item.getElemEdges()) == 1:
                        adjaNodes.append(item)    
    
    return adjaNodes
    

# -------------------------------------------------------------------------------------------------------------------------------------------
def getRowOfNodes(currElem, nextElem): 
    """
    # Get row of nodes from two adjacent elements. Input elments must be of 3-D HEX type(linear or quadratic). 
    # Four paths will be found, each path passes through one vertex node of common nodes between adjacent elements. As shown below:
    #                                 1----×---2----×----3
    #                                /|       /|        /|  
    #                               4-|--×---5-|--×----6 |
    #                               | |   I  | |   II  | |    ------> Sweep direction.   Common vertex nodes are 2,5,8 and 11. 
    #                               | 7---×--|-8----×--|-9
    #                               |/       |/        |/  
    #                              10---×---11----×---12
    # Inputs: 
    #        currElem: Current element. Element I in figure above.
    #        nextElem: Next element. Element II in figure above.    
    # Outputs: 
    #        nodeList: List of 4 element node lists associated with input elements.
    #                  In figure above, 4 paths we get are node list [1,2,3],[4,5,6],[7,8,9] and [10,11,12]. 
    #                  If elements are quadratic, additional nodes × will be inserted into node list, like [1, ×, 2, ×, 3].
    """  
    if ((re.search('3D8' ,str(currElem.type)) == None) and 
        (re.search('3D20',str(currElem.type)) == None) and
        (re.search('3D8' ,str(nextElem.type)) == None) and
        (re.search('3D20',str(nextElem.type)) == None)) :
        return 
    else: 
        # Get common VERTEX nodes of input elements.    
        nodes1 = list(currElem.getNodes())
        nodes2 = list(nextElem.getNodes())
        commonNodes = [node for node in nodes1 if node in nodes2]
        edgeNumbers = [len(node.getElemEdges()) for node in commonNodes]
        zipNodes = list(zip(commonNodes,edgeNumbers))
        sortNodes = sorted(zipNodes, key=lambda x:x[1], reverse=True)[:4]  # Nodes on vertex always share more element edges than nodes on edge.
        vertNodes = list(zip(*sortNodes))[0]
        
        # Get center of common face of input elements.
        avgX = 0.0
        avgY = 0.0
        avgZ = 0.0
        for node in vertNodes:
            avgX = avgX + node.coordinates[0]
            avgY = avgY + node.coordinates[1]
            avgZ = avgZ + node.coordinates[2]
        faceCent = [avgX/4.0, avgY/4.0, avgZ/4.0]
        
        # Get row of nodes.
        vNorm = lambda x : np.array([0.0,0.0,0.0]) if(nla.norm(x)==0) else x/nla.norm(x)
        dirBack = np.array(currElem._getCentroid()) - np.array(faceCent) 
        dirBack = vNorm(dirBack)  # Direction from common face center to current element center.
        dirFrnt = np.array(nextElem._getCentroid()) - np.array(faceCent) 
        dirFrnt = vNorm(dirFrnt)  # Direction from common face center to next element center.
        nodeList = []
        for vNode in vertNodes:
            adjaNodes = getAdjacentNodes(vNode)
            direction = [vNorm(np.array(aNode.coordinates) - np.array(vNode.coordinates)) for aNode in adjaNodes]             
            dotBack = [np.dot(dir, dirBack) for dir in direction]
            dotFrnt = [np.dot(dir, dirFrnt) for dir in direction]
            zipNodes = list(zip(adjaNodes, dotBack, dotFrnt))
            sortBack = sorted(zipNodes, key=lambda x:x[1], reverse=True)
            sortFrnt = sorted(zipNodes, key=lambda x:x[2], reverse=True)
            lastNode = list(zip(*sortBack))[0][0]
            nextNode = list(zip(*sortFrnt))[0][0]
            if (re.search('3D8' ,str(currElem.type)) != None):
                nodes = [lastNode, vNode, nextNode]
            else:
                backEdges = list(lastNode.getElemEdges())
                backNodes = list(backEdges[0].getNodes())
                backNodes.remove(lastNode)
                backNodes.remove(vNode)
                frntEdges = list(nextNode.getElemEdges())
                frntNodes = list(frntEdges[0].getNodes())
                frntNodes.remove(nextNode)
                frntNodes.remove(vNode)  
                nodes = [backNodes[0], lastNode, vNode, nextNode, frntNodes[0]]
            nodeList.append(nodes)  

    return nodeList   

 
# -------------------------------------------------------------------------------------------------------------------------------------------
def mergeList(list1, list2): 
    """
    # Merge two lists which have common items in inner list. Inner list should be quite different from each other.
    # For example: 
    #       list1 = [[1,2,3],  list2 = [[3,10,11],  then mergList = [[1,2,3,10,11],
    #                [4,5,6],           [9,12,13],                   [4,5,6,14,15],
    #                [7,8,9]]           [6,14,15]]                   [7,8,9,12,13]].
    # Inputs: 
    #        Two lists.    
    # Outputs: 
    #        One list.
    """
    if not len(list1):  # List1 is empty.
        return list2  
    elif not len(list2): # List2 is empty.
        return list1
    else:
        mergList = []
        for lst1 in list1:
            # Search inner lists which have common items.
            for lst2 in list2:
                if any([i in lst1 for i in lst2]):
                    lst1 = lst1 + lst2
            # Remove duplicated items.
            tmp = [] 
            for item in lst1:
                if item not in tmp:
                    tmp.append(item)                
            mergList.append(tmp)
            
        return mergList    
   
   
# -------------------------------------------------------------------------------------------------------------------------------------------   
def getTopologyAttributes(crossSection): 
    """
    # Get all element nodes on input face object and save them in crossFaceNodes.
    # Get all element faces on input face object and save them in crossElemFaces.
    # Get all elements from input face object via topology and save them in topoElements.
    # Get other topology information.
    # Inputs: 
    #        crossSection    : Weld line cross section.
    # Outputs:               
    #        crossFaceNodes  : List of element nodes found with input crossSection. 
    #        crossElemFaces  : List of element faces found with input crossSection.
    #        topoElements    : List of elements found with input crossSection via topology.
    #        isClosedCurve   : Is weld line a closed curve?
    #        topoDirection   : Direction vector of topology, from topology start to topology end.
    #        iniTopoFaceNodes: Nodes on start face of topology.
    #        endTopoFaceNodes: Nodes on end face of topology.
    """ 
    topoInfo = {}
    
    # Get element faces on cross face.
    crossElemFaces = getElemFacesFromFace(crossSection)
    topoInfo['crossElemFaces'] = crossElemFaces
    
    # Get vertex nodes on cross face.  
    crossFaceNodes = []  
    for elemFace in crossElemFaces:
        for node in list(elemFace.getNodes()):
            if (node not in crossFaceNodes) and (len(node.getElemEdges())!=1):
                crossFaceNodes.append(node)
    topoInfo['crossFaceNodes'] = crossFaceNodes
    
    # Get elements via topology.
    topoElements = []
    for elemFace in crossElemFaces:
        elements = list(elemFace.getElementsViaTopology()) 
        if elements not in topoElements:
            topoElements.append(elements) 
    # Let sweep direction of topoElements[0] be topology direction.        
    topoDirection = np.array(topoElements[0][-1]._getCentroid()) - np.array(topoElements[0][0]._getCentroid())
    topoDirection = topoDirection/nla.norm(topoDirection)
    # If start topology element shares common nodes with end topology element, the welding line is a closed curve. 
    isClosedCurve = any([item in topoElements[0][-1].getNodes() for item in topoElements[0][0].getNodes()]) 
    # If welding line cross section includes more than one face, topology sweep direction may be different between each topology obtained by using getElementsViaTopology() method.
    # Rearrange all topologies by topoDirection.     
    for i in range(1,len(topoElements)):
        direction = np.array(topoElements[i][-1]._getCentroid()) - np.array(topoElements[i][0]._getCentroid())
        direction = direction/nla.norm(direction)
        if np.dot(topoDirection, direction) < 0.0:
            topoElements[i].reverse()
    topoInfo['isClosedCurve'] = isClosedCurve
    topoInfo['topoDirection'] = topoDirection    
    topoInfo['topoElements'] = topoElements
    
    # Get row of nodes by topology sweep path.
    topoNodes = []
    for elements in topoElements:
        # Get row of node in single elment topology.
        fullRow = []
        for i in range(len(elements)-1):
            partRow = getRowOfNodes(elements[i], elements[i+1])
            fullRow = mergeList(fullRow, partRow)
 
        for lst in fullRow:
            if lst not in topoNodes:         
                topoNodes.append(lst)        
    topoInfo['topoNodes'] = topoNodes   

    # Get nodes on start and end topology face. 
    topoInfo['iniTopoFaceNodes'] = [nodes[0] for nodes in topoNodes] 
    if isClosedCurve:    
        topoInfo['endTopoFaceNodes'] = topoInfo['iniTopoFaceNodes']
    else:
        topoInfo['endTopoFaceNodes'] = [nodes[-1] for nodes in topoNodes]        
   
    return topoInfo      
    

# -------------------------------------------------------------------------------------------------------------------------------------------       
def displayView(weld, arrowName, arrowLength, arrowColor, headStyle, textOffset, textName, textContent, textColor, bgColor, fontSize):
    """
    # Plot welding direction arrow and text for welding line in abaqus viewport.
    # Inputs: 
    #        arrowName  : Optional argument, name of arrow.
    #        arrowLength: Optional argument, length of arrow, default value is 15.
    #        arrowColor : Optional argument, color of arrow, default value is #0000FF, i.e. blue.
    #        headStyle  : Optional argument, style of arrow start, default value is HOLLOW_CIRCLE.
    #
    #        textName   : Optional argument, name of text object.
    #        textContent: Optional argument, contents of text.
    #        textColor  : Optional argument, color of text, default value is #FF0000, i.e. red.
    #        bgColor    : Optional argument, color of text background, default value is #CCCCCC.
    #        fontSize   : Optional argument, font size of text, default value is 14.
    """ 
    p = weld.p
    p.setValues(geometryRefinement=EXTRA_FINE)
    
    # Create arrow and text.        
    arrowStart = weld.weldNodeList[0].coordinates
    arrowDirec = weld.iniVector
    arrowEnd = tuple(np.array(arrowStart) + arrowLength*arrowDirec)
    
    pArrow = mdb.Arrow(name=arrowName, startAnchor=arrowStart, endAnchor=arrowEnd, color=arrowColor,
                       startHeadStyle=headStyle, lineStyle=DASHED, lineThickness=MEDIUM) 
                       
    pText = mdb.Text(name=textName, text=textContent, offset=textOffset, anchor=arrowEnd, 
                     font='-*-verdana-bold-r-normal-*-*-{}-*-*-p-*-*-*'.format(fontSize*10), color=textColor, 
                     backgroundStyle=OTHER, backgroundColor=bgColor, box=ON, justification=CENTER) 
    
    # Create viewport.      
    viewportName = 'Viewport: {}'.format(weld.modelName)
    if viewportName not in session.viewports.keys():
        session.Viewport(name=viewportName, origin=(0.0, 0.0), width=100, height=100)

    session.viewports[viewportName].makeCurrent()
    session.viewports[viewportName].maximize() 
    
    # Display arrow and text.        
    session.viewports[viewportName].setValues(displayedObject=p)
    session.viewports[viewportName].view.setProjection(projection=PARALLEL)
    session.viewports[viewportName].plotAnnotation(annotation=pArrow)
    session.viewports[viewportName].plotAnnotation(annotation=pText)        