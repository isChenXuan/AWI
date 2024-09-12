# --*coding=UTF-8*--
from abaqus import *
from abaqusConstants import *
from csv import *
from dataModule import DataProcess, Circle
from weldUtils import getTopologyAttributes, getNodeFromPoint, displayView
import numpy as np
import numpy.linalg as nla


class WeldLine():

    def __init__(self, modelName, partName, crossFace, curveType):
        """
        # WeldLine class defines a single welding line object, which has attributes of geometry parameters and welding process parameters.
        #    crossFace : Ccross section face of welding line, selected in abaqus GUI.
        #    curveType : Type of weld line, coud be 'StraightLine', 'Circle' and 'ArbitraryCurve'.
        """
        self.modelName = modelName
        self.partName  = partName
        self.crossFace = crossFace
        self.curveType = curveType
        
        # Fetch abaqus part.
        p = mdb.models[modelName].parts[partName]
        # Check mesh control in welding zone.
        try:
            cells = p.sets['weld'].cells
        except:
            raise ValueError, "Error: no weld set found, please create a set named 'weld' which includes all cells of welding zone under part {}!".format(partName)
        else:
            if any(p.getMeshControl(region=cell, attribute=TECHNIQUE) != SWEEP for cell in cells):
                raise ValueError, "Error: meshing technique in welding zone must be SWEEP!"
            else:
                self.p = p
                
        # Fetch abaqus assembly.     
        a = mdb.models[modelName].rootAssembly
        instName = list(a.instances.keys())
        nameList = [key for key in instName if a.instances[key].partName == partName]        
        if len(nameList) >= 2:
            raise ValueError, "Error: more than one instance associated with part {} found, plese delete redundant instances!".format(partName)
        elif len(nameList) == 0:
            raise ValueError, "Error: no instance associated with part {} found, plese make Assembly first!".format(partName)
        else:
            instance = a.instances[nameList[0]]
            if not instance.dependent:
                raise ValueError, "Error: please make instance Dependent!"  
            else:
                self.instance = instance             
        
        # Get all elements in weld region. Unpack self.topoInfo['topoElements'], which is list of list of elements.
        self.topoInfo  = getTopologyAttributes(crossFace)
        self.__elementInRegion = [elem for lst in self.topoInfo['topoElements'] for elem in lst]  
        self.__layer = len(self.topoInfo['topoElements'][0])
        
    
    def setRegionParameters(self, startPoint, alongPoint, toePoint, segNumber):  
        """ 
        # Geometry parameters:
        #    startPoint : Point at which welding process starts, selected in abaqus GUI. 
        #    alongPoint : Combined with start point to establish welding forward direction, selected in abaqus GUI.
        #    toePoint   : Start point of welding toe line, selected in abaqus GUI.
        #    segNumber  : Segement number that a welding line splited into, also abaqus steps needed to finish welding that line, integer. 
        """      
        # Geometry parameters.       
        startNode = getNodeFromPoint(startPoint)
        toeNode   = getNodeFromPoint(toePoint)
        # Check input point.         
        if (len(startNode.getElemEdges()) == 1) or (len(toeNode.getElemEdges()) == 1) :
            raise ValueError, "Error: start point cannot be middle node of a quadratic element edge!"
        #elif (self.__layer//segNumber != 0) or (segNumber <= 0) or (segNumber > self.__layer) :
        #    raise ValueError, "Error: there are {} layers of element along sweep direction, segment number should be divisor of {}, minimum segment number is 1 and maximum segment number is {}!".format(self.__layer,self.__layer,self.__layer)             
        else:
            self.startNode = startNode
            self.toeNode   = toeNode  
            self.alongNode = getNodeFromPoint(alongPoint) 
            alongVector = np.array(self.alongNode.coordinates)-np.array(self.startNode.coordinates)
            self.alongVector = alongVector/nla.norm(alongVector)
            self.segNumber = segNumber
            
            highlight(self.startNode)
            highlight(self.toeNode)
            highlight(self.alongNode)
            
            # Update attributes.
            self.__sortItems__()
            self.__groupItems__()
            self.__calWeldLength__()
            
        
    def setWeldingParameters(self, weldVoltage, weldCurrent, weldVelocty, weldYita, coolTime, deltTemp, heatType, heatPara):
        """
        # Welding parameters:
        #    weldVoltage : Voltage used during welding process, in V, float.   
        #    weldCurrent : Current used during welding process, in A, float.
        #    weldVelocty : Velocity of welding process, in mm/s, float.  
        #    weldYita    : Heat efficiency during welding process, float.   
        #    coolTime    : Cooling down time after finishing a welding line, in second, float.       
        #    deltTemp    : Maximum allowable temperature change during welding process, in C degree, float.
        #    heatType    : Heat source type. 
        #    heatPara    : Heat source parameters, tuple of floats.       
        """
        # Welding parameters.
        self.weldVoltage = weldVoltage   
        self.weldCurrent = weldCurrent
        self.weldVelocty = weldVelocty  
        self.weldYita    = weldYita 
        if coolTime <= 0.0:
            self.coolTime = 1.0e-7
        else:            
            self.coolTime = coolTime     
        self.deltTemp = deltTemp
        self.__time = [item/self.weldVelocty for item in self.__length]  # Time spent from start node to current node.  
        self.weldTime = self.__time[-1]
        self.__totaTime = self.weldTime + self.coolTime # Total time spent for current welding line, i.e., weld time plus cool time. 
        
        # Heat source settings.
        self.heatType = heatType
        if self.heatType == 'DoubleEllipsoid':
            self.a1 = heatPara[0]                            # a1  ：热源前进方向椭球半轴长,毫米,建议a1=b=0.45*焊缝宽度
            self.a2 = heatPara[1]                            # a2  ：热源后方方向椭球半轴长,毫米,建议a2=2*a1 
            self.b  = heatPara[2]                            # b   ：焊缝宽度方向椭球半轴长,毫米,建议a1=b=0.45*焊缝宽度
            self.c  = heatPara[3]                            # c   ：焊缝深度方向椭球半轴长,毫米,建议c=0.90*焊缝深度
            POWR = 1000.0*weldYita*weldVoltage*weldCurrent   # POWR：输入功率,毫瓦       
            f1 = 2.0/(1 + self.a2/self.a1)                   # f1  ：前半椭球能量系数,
            f2 = 2.0 - f1                                    # f2  ：后半椭球能量系数 
            self.qm1 = 1.8663*f1*POWR/self.a1/self.b/self.c  # qm1 ：前半椭球能量峰值   
            self.qm2 = 1.8663*f2*POWR/self.a2/self.b/self.c  # qm2 ：后半椭球能量峰值        
        elif self.heatType == 'PlanarGauss': 
            self.a = heatPara[0]
            self.b = heatPara[1]
            self.c = heatPara[2]
        else:
            pass  


    def __sortItems__(self): 
        """
        # Rearrange nodes and elements in topology.
        """        
        # Find weld line node list and toe line node list in topoNodes.
        for nodeList in self.topoInfo['topoNodes']:
            if self.startNode in nodeList:
                topoWeldNodeList = nodeList
            if self.toeNode in nodeList:
                topoToeNodeList = nodeList  
        
        # Rearrange nodes and elements.
        if (self.startNode in self.topoInfo['iniTopoFaceNodes']):
            weldNodeList = topoWeldNodeList
            toeNodeList  = topoToeNodeList
            elementList  = self.topoInfo['topoElements'] 
        elif (self.startNode in self.topoInfo['endTopoFaceNodes']):
            weldNodeList = list(reversed(topoWeldNodeList))
            toeNodeList  = list(reversed(topoToeNodeList))  
            for elems in self.topoInfo['topoElements']:
                elems.reverse()
            elementList  = self.topoInfo['topoElements'] 
            if self.topoInfo['isClosedCurve']:
                weldNodeList = weldNodeList[-1:] + weldNodeList[:-1]
                toeNodeList  = toeNodeList[-1:] + toeNodeList[:-1]
        else:  
            # Find start node index in topology node list.
            for node in topoWeldNodeList:
                if self.startNode == node:
                    startNodeIndex = topoWeldNodeList.index(node)
            
            frntDir = np.array(topoWeldNodeList[startNodeIndex+1].coordinates) - np.array(topoWeldNodeList[startNodeIndex].coordinates)
            frntDir = frntDir/nla.norm(frntDir)
            backDir = np.array(topoWeldNodeList[startNodeIndex-1].coordinates) - np.array(topoWeldNodeList[startNodeIndex].coordinates)
            backDir = backDir/nla.norm(backDir)
        
            # Rearrange nodes. 
            if self.topoInfo['isClosedCurve']:
                weldNodeList = topoWeldNodeList[startNodeIndex:] + topoWeldNodeList[:startNodeIndex]
                toeNodeList  = topoToeNodeList[startNodeIndex:] + topoToeNodeList[:startNodeIndex]
            else:
                if np.dot(self.alongVector, frntDir) > np.dot(self.alongVector, backDir):
                    weldNodeList = topoWeldNodeList[startNodeIndex:]
                    toeNodeList  = topoToeNodeList[startNodeIndex:]
                else:
                    weldNodeList = list(reversed(topoWeldNodeList[:startNodeIndex+1]))
                    toeNodeList  = list(reversed(topoToeNodeList[:startNodeIndex+1]))  
        
            # Find start element index in topology element list.        
            iniElements = [elem for elem in weldNodeList[0].getElements() if elem in weldNodeList[1].getElements()]
            for elemList in self.topoInfo['topoElements']:
                for elem in elemList:
                    if iniElements[0] == elem:
                        startElemIndex = elemList.index(elem)  
        
            # Rearrange elements. 
            if self.topoInfo['isClosedCurve']:
                elementList  = [item[startElemIndex:]+item[:startElemIndex] for item in self.topoInfo['topoElements']]  
            else:
                if np.dot(self.alongVector, frntDir) > np.dot(self.alongVector, backDir):
                    elementList  = [item[startElemIndex:] for item in self.topoInfo['topoElements']] 
                else:
                    elementList  = [item[:startElemIndex+1].reverse() for item in self.topoInfo['topoElements']] 
        
        # Whether to flip arranging direction for closed curve.
        tmpDirection = np.array(weldNodeList[1].coordinates) - np.array(weldNodeList[0].coordinates)
        iniDirection = tmpDirection/nla.norm(tmpDirection)
        if np.dot(self.alongVector, iniDirection) < 0.0 and self.topoInfo['isClosedCurve']:
            weldNodeList.reverse()
            weldNodeList = weldNodeList[-1:] + weldNodeList[:-1]
            toeNodeList.reverse()
            toeNodeList = toeNodeList[-1:] + toeNodeList[:-1]
            for item in elementList:
                item.reverse()
 
        # Get element list and labels.
        self.elementList = elementList
        self.__elementLabels  = [[item.label for item in lst] for lst in self.elementList]
        
        # Get node list and labels.
        self.weldNodeList = weldNodeList
        self.toeNodeList = toeNodeList
        
        weldNodeLabels = [item.label for item in self.weldNodeList]
        toeNodeLabels = [item.label for item in self.toeNodeList] 
        
        # Get node coordinates under assembly.
        if any(n.instanceName==None for n in self.weldNodeList): # There is part node in weldNodeList.
            weldNodeCoords = []        
            for node in self.weldNodeList:
                indx = node.label
                aNode = self.instance.nodes[indx-1]
                weldNodeCoords.append(aNode.coordinates) 
        else:
            weldNodeCoords = [item.coordinates for item in self.weldNodeList]
                
        if any(n.instanceName==None for n in self.toeNodeList): # There is part node in toeNodeList.  
            toeNodeCoords = []        
            for node in self.toeNodeList:
                indx = node.label
                aNode = self.instance.nodes[indx-1]
                toeNodeCoords.append(aNode.coordinates) 
        else:
            toeNodeCoords  = [item.coordinates for item in self.toeNodeList]

        # If weld line is closed curve, duplicate firt node to node list end.
        if self.topoInfo['isClosedCurve']:
            weldNodeLabels.append(weldNodeLabels[0])
            toeNodeLabels.append(toeNodeLabels[0])
            weldNodeCoords.append(weldNodeCoords[0]) 
            toeNodeCoords.append(toeNodeCoords[0])        
        
        # Add attributes.         
        self.__weldNodeLabels = weldNodeLabels
        self.__toeNodeLabels  = toeNodeLabels
        
        self.__weldNodeCoords = np.array(weldNodeCoords)
        self.__toeNodeCoords  = np.array(toeNodeCoords )
        
    
    def __groupItems__(self):
        """
        # Regroup elements by layer and segment.
        # Get elements in each layer and segment, used in abaqus model change setting when activating welding elements segment by segment.
        # Each segment may include one or more element layers.
        """  
        # Get elements by layer.
        self.__elementByLayer = [[self.elementList[i][j] for i in range(len(self.elementList))] for j in range(self.__layer)] 
        # Get elements by segment.
        elementBySegment = []
        # Divide layer into segments by segNumber.
        # Note: for integer m > n, m can be divided into m = n*q + r = r*(q+1) + (n-r)*q.
        q = self.__layer//self.segNumber
        r = self.__layer%self.segNumber
        for i in range(r):                                 # i = r
            elems = []                                    
            for j in range(i*(q+1), (i+1)*(q+1)):          # j = i*(q+1)
                elems = elems + self.__elementByLayer[j]
            elementBySegment.append(elems)
            
        for i in range(self.segNumber - r):                # i = n-r
            elems = []
            for j in range(r*(q+1)+i*q, r*(q+1)+(i+1)*q):  # j = r*(q+1)+i*q
                elems = elems + self.__elementByLayer[j]
            elementBySegment.append(elems)             
        
        self.elementBySegment = elementBySegment
        
        
    def __calWeldLength__(self): 
        """
        # Calculate weld length.
        """
        vNorm = lambda x : np.array([0.0,0.0,0.0]) if(nla.norm(x)==0) else x/nla.norm(x)
        
        if self.curveType == 'StraightLine':
            dir = [self.__weldNodeCoords[i] - self.__weldNodeCoords[0] for i in range(len(self.__weldNodeCoords))]
            self.__length = list(map(nla.norm, dir))
            self.weldLength = self.__length[-1]
            self.iniVector = vNorm(dir[1])
            
        elif self.curveType == 'Circle':
            # Get circle points.
            n = len(self.__weldNodeCoords)
            n2 = n//3
            n3 = 2*n2
            
            # Get weld line circle.
            PW1 = self.__weldNodeCoords[0]
            PW2 = self.__weldNodeCoords[n2]
            PW3 = self.__weldNodeCoords[n3] 
            self.weldCircle = Circle(PW1, PW2, PW3)
            
            # Get weld toe line circle.
            PT1 = self.__toeNodeCoords[0]
            PT2 = self.__toeNodeCoords[n2]
            PT3 = self.__toeNodeCoords[n3] 
            self.toeCircle = Circle(PT1, PT2, PT3)            
            
            # Get tangent vector at start node.
            self.iniVector = self.weldCircle.vy
            
            # Get angle between point[i] and point[0].
            VN = [vNorm(self.__weldNodeCoords[i] - self.weldCircle.PC) for i in range(len(self.__weldNodeCoords))] # Vectors of circle center to weld nodes.
            VC = [vNorm(np.cross(self.weldCircle.vx, V)) for V in VN]
            VD = [np.dot(self.weldCircle.vx, V) for V in VN]
            cosTheta = np.clip(np.array(VD), -1, 1)
            theta = []
            for i in range(len(VC)):
                if nla.norm(VC[i]+VC[1]) < 1.0e-6:
                    tmp = 2.0*np.pi - np.arccos(cosTheta[i]) 
                else:
                    tmp = np.arccos(cosTheta[i])  
                theta.append(tmp)                 

            self.__length = [self.weldCircle.RA*t for t in theta]
            self.weldLength = self.__length[-1]
                        
        elif self.curveType == 'ArbitraryCurve':    
            self.weldDp = DataProcess(self.__weldNodeCoords)
            self.__length = self.weldDp.arcLength   # Length from start node to current node.
            self.weldLength = self.__length[-1]
            self.iniVector = self.weldDp.iniVector
            self.toeDp = DataProcess(self.__toeNodeCoords)
        
 
class WeldSequence():

    def __init__(self, weldParameters):
        """
        # WeldSequence class defines a sequence of welding line object.
        # Welding parameters inputs example:
        # Note : POINT is selected by user in GUI, maybe geometry vertex or element node, wehre NODE is node associated with input point.     
        # Model input:
         ====================================================================
        | Line NO. |  Model Name |  Part  Name |   Cross Face  | Curve Type  |
        |------------------------------------------------------|-------------|
        |     1    | pick in GUI | pick in GUI |  pick in GUI  |  selected   |
        |------------------------------------------------------|-------------|
        |     2    | pick in GUI | pick in GUI |  pick in GUI  |  selected   |
        |------------------------------------------------------|-------------|
        |    ...   |             |             |               |             |
         ====================================================================   
        # Region parameters:
         =================================================================
        | Line NO. | Start Point | Along Point |  Toe Point  |  Segments  |
        |-----------------------------------------------------------------|
        |     1    | pick in GUI | pick in GUI | pick in GUI | user enter |
        |-----------------------------------------------------------------|
        |     2    | pick in GUI | pick in GUI | pick in GUI | user enter |
        |-----------------------------------------------------------------|
        |    ...   |             |             |             |            |
         ================================================================= 
        # Welding parameters:
         =============================================================================================================================
        | Line NO. | Weld Voltage | Weld Current | Weld Velocty |   Weld Yita   |  Cool Time   | delta Temp | Heat Type  | Heat Para  |
        |---------------------------------------------------------------------------------------------------|--------------------------
        |     1    |  user enter  |  user enter  |  user enter  |   user enter  |  user enter  | user enter | user enter | user enter |
        |---------------------------------------------------------------------------------------------------|--------------------------
        |     2    |  user enter  |  user enter  |  user enter  |   user enter  |  user enter  | user enter | user enter | user enter |
        |---------------------------------------------------------------------------------------------------|--------------------------
        |    ...   |              |              |              |               |              |            |            |            |
         =============================================================================================================================   
        """ 
        self.weldParameters = weldParameters
     
        # Create welds.
        weldList = []
        for wp in self.weldParameters:
            try:
                weld = WeldLine(modelName = wp['modelName'], 
                                partName  = wp['partName'], 
                                crossFace = wp['crossFace'],
                                curveType = wp['curveType'])  
                weld.setRegionParameters(startPoint = wp['startPoint'], 
                                         alongPoint = wp['alongPoint'], 
                                         toePoint   = wp['toePoint'], 
                                         segNumber  = wp['segNumber'])        
                weld.setWeldingParameters(weldVoltage = wp['weldVoltage'], 
                                          weldCurrent = wp['weldCurrent'], 
                                          weldVelocty = wp['weldVelocty'], 
                                          weldYita    = wp['weldYita'], 
                                          coolTime    = wp['coolTime'], 
                                          deltTemp    = wp['deltTemp'], 
                                          heatType    = wp['heatType'],
                                          heatPara    = wp['heatPara'])
            except ValueError, error:
                print error
                return
            else: 
                weldList.append(weld)
        self.weldList = weldList 
        
        # Call functions in order.
        global firstStep
        firstStep = 'Step-1-removeAll'
        
        self.__getWeldSequence__()
        
        self.createSets()
        self.createStep()
        self.createInts()
        
        self.__getTimeSequence__()
        
        self.abaqConfig() 
        self.display()         

        
    def __getWeldSequence__(self):
        """
        # Get weld and instance sequence.
        """  
        # Get model list.
        modelList = []
        for w in self.weldList:  
            model = mdb.models[w.modelName]
            if model not in modelList:
                modelList.append(model)
        self.modelList = modelList 
        
        # Get sequence.
        weldSeq = {} 
        partSeq = {}
        instSeq = {}       
        for model in self.modelList:
            name = model.name
            welds = [w for w in self.weldList if w.modelName==name] 
            weldSeq[name] = welds
            # Get parts and instances sequence.
            partList = []
            instList = []
            for w in welds:
                part = w.p
                inst = w.instance 
                if part not in partList:
                    partList.append(part)  
                if inst not in instList:
                    instList.append(inst)                
            partSeq[name] = partList
            instSeq[name] = instList
        
        self.weldSeq = weldSeq
        self.partSeq = partSeq
        self.instSeq = instSeq        
        
        
    def createSets(self):
        """
        # Create sets for cells, faces, nodes and elements under assembly in abaqus.
        """
        for name in list(self.weldSeq.keys()):
            a = mdb.models[name].rootAssembly
            instanceList = self.instSeq[name]
            
            # Create set for all surfaces, cells and weld elements. Surfaces need to manually update by user.
            faces = instanceList[0].faces
            cells = instanceList[0].cells
            elems = instanceList[0].sets['weld'].elements
            if len(instanceList) > 1:
                for inst in instanceList[1:]:
                    faces = faces + inst.faces
                    cells = cells + inst.cells
                    elems = elems + inst.sets['weld'].elements
            a.Surface(side1Faces=faces, name='allSurfaces')
            a.Set(cells=cells, name='allCells')
            a.Set(elements=elems, name='allWeldElems')
            
            # Create sets for nodes and elements.
            weldList = self.weldSeq[name]
            for indx,weld in enumerate(weldList, 1):
                # Create set for all nodes in welding and toe line.
                a.SetFromNodeLabels(name='weldLineNodes-{}'.format(indx), nodeLabels=((weld.instance.name, weld._WeldLine__weldNodeLabels),))   
                a.SetFromNodeLabels(name= 'weldToeNodes-{}'.format(indx), nodeLabels=((weld.instance.name, weld._WeldLine__toeNodeLabels ),)) 
                
                # Create sets for elements in welding zone by segment.
                segElemLabels = [[elem.label for elem in lst] for lst in weld.elementBySegment] 
                for j in range(weld.segNumber):
                    segSetName = "w{}-{}".format(indx,j+1)
                    a.SetFromElementLabels(name=segSetName, elementLabels=((weld.instance.name, segElemLabels[j]),))       
                
                # Create set for elements in each welding line.
                regElemLabels = [elem.label for elem in weld._WeldLine__elementInRegion]
                a.SetFromElementLabels(name="weldRegion-{}".format(indx), elementLabels=((weld.instance.name, regElemLabels),)) 
        
        print "Surface for film condition setting have been collected in 'allSurfaces', please manually update it before run analysis!"
        print "Set 'allCells' for all cells in each model has been created!" 
        print "Set 'allWeldElems' for all weld elements in each model has been created!"          
        print "Sets 'weldLineNodes' for all welding nodes in each welding line have been created!" 
        print "Sets 'weldToeNodes' for all welding toe nodes in each welding line have been created!"
        print "Sets for all welding segments in each welding line have been created!"
        print "Sets for all elements in each welding line has been created!"         
        
    
    def createStep(self):    
        """
        # Create abaqus analysis steps.
        """ 
        for name in list(self.weldSeq.keys()):
            model = mdb.models[name]
            weldList = self.weldSeq[name]
            
            for indx,weld in enumerate(weldList, 1):
                # Get welding time for each welding segment.
                timePerSeg = weld.weldTime/weld.segNumber
                
                # Set welding step name.
                weldStep = ['Step-w{}-{}'.format(indx, i+1) for i in range(weld.segNumber)]
        
                # Define first Abaqus step, in which all elements of welding lines were removed.
                if firstStep not in model.steps.keys():
                    model.HeatTransferStep(name=firstStep, previous='Initial', timePeriod=1.0e-7, maxNumInc=1000, 
                        initialInc=1.0e-7, minInc=1.0e-7, maxInc=1.0e-7, deltmx=weld.deltTemp) 
        
                # Create welding step for each welding line.
                for i in range(weld.segNumber):
                    if i == 0:
                        if indx == 1:
                            preStep = firstStep
                        else:
                            preStep = 'Step-coolDown-{}'.format(indx-1)
                    else:
                        preStep = weldStep[i-1]  
                
                    model.HeatTransferStep(name=weldStep[i], previous=preStep, timePeriod=timePerSeg, 
                        maxNumInc=1000, initialInc=timePerSeg/10**4, minInc=timePerSeg/10**5, maxInc=timePerSeg, 
                        deltmx= weld.deltTemp)     
                        
                # Create cooling down step for each welding line. 
                model.HeatTransferStep(name='Step-coolDown-{}'.format(indx), previous=weldStep[-1],
                    timePeriod=weld.coolTime,maxNumInc=1000, initialInc=weld.coolTime/10**4, minInc=weld.coolTime/10**5, 
                    maxInc=weld.coolTime, deltmx= weld.deltTemp)
                             
        print "Steps in each model have been created!"     
 
 
    def createInts(self):
        """
        # Create interactions to reactivate elements set of welding segments.
        """ 
        for name in list(self.weldSeq.keys()):
            model = mdb.models[name]
            a = mdb.models[name].rootAssembly
            weldList = self.weldSeq[name]
            
            # Remove all welding elements.
            model.ModelChange(name='Int-1-removeAll', createStepName=firstStep, region=a.sets['allWeldElems'], 
                              regionType=ELEMENTS, activeInStep=False, includeStrain=False)
            
            for indx,weld in enumerate(weldList, 1):                  
                # Reactive element segments.                      
                for i in range(weld.segNumber):
                    intName = "int-w{}-{}".format(indx,i+1)
                    stpName = "Step-w{}-{}".format(indx, i+1)
                    segName = "w{}-{}".format(indx,i+1)
                    segRegion = a.sets[segName]
                    model.ModelChange(name=intName, createStepName=stpName, region=segRegion, regionType=ELEMENTS, activeInStep=True, 
                        includeStrain=False)
                    
        print "Interactions in each model have been created!"     
    

    def __getTimeSequence__(self):
        """
        # Get step and time sequence.
        """     
        stepSeq = {}
        timeSeq = {}        
        for name in list(self.weldSeq.keys()):
            model = mdb.models[name]
            weldList = self.weldSeq[name]

            # Get abaqus KSTEP for each model. 
            #  For example, there are 2 models in abaqus, each model has 2 weld lines:
            #      model-1 line1 is diveded into 4 segs, line2 is diveded into 4 segs.
            #      model-2 line1 is diveded into 3 segs, line2 is diveded into 5 segs.
            #    For model-1, there will be 12 abaqus Steps(including 'Initial'):
            #      line1 start KSTEP is 2, end with 5. line2 start KSTEP is 7, end with 10.
            #    For model-2, there will be 12 abaqus Steps(including 'Initial'):
            #      line1 start KSTEP is 2, end with 4. line2 start KSTEP is 6, end with 10. 
            kstep = []
            a, b = 2, 2+weldList[0].segNumber-1
            kstep.append((a, b))
            if len(weldList) > 1:
                for weld in weldList[1:]:
                    a = b + 2
                    b = a + weld.segNumber-1
                    kstep.append((a, b))
            stepSeq[name] = kstep
            
            # Get total time spent before current welding line begins.
            time = []
            keys = list(model.steps.keys())
            for i in range(len(kstep)):
                indx = kstep[i][0]
                currStep = keys[indx]
                t = 0.0
                while True:
                    preStep = model.steps[currStep].previous
                    if preStep == 'Initial':
                        break
                    t = t + model.steps[preStep].timePeriod    
                    currStep = preStep
                time.append(t)               
            timeSeq[name] = time          
        
        self.stepSeq = stepSeq         
        self.timeSeq = timeSeq 
    
    
    def abaqConfig(self, tZero=-273.15, stefanBoltzmann=5.669E-011, filmCoeff=0.025, Tambient=21.1):
        """
        # Set parameters for abaqus model and create steps and interactions.
        # Inputs:   
        #        tZero          : Absolute zero temperature, default value is -273.15C.
        #        stefanBoltzmann: Stefan-Boltzmann constant, default value is 5.669E-011mJ/s/mm^2/K^4.
        #        filmCoeff      : Film coefficient, default value is 0.025mJ/s/mm^2/K.
        #        Tambient       : Ambient temperature, default value is 21.1C.    
        # Outputs: 
        #        None.
        """ 
        for name in list(self.weldSeq.keys()):
            model = mdb.models[name]
            a = mdb.models[name].rootAssembly
            
            # Set physical constants for abaqus model.
            model.setValues(absoluteZero=tZero, stefanBoltzmann=stefanBoltzmann)
            
            # Create abaqus field output requests.
            model.FieldOutputRequest(name='F-Output-1',createStepName=firstStep, variables=PRESELECT)   

            # Create constant amplitude curve.
            model.TabularAmplitude(name='Amp-1', timeSpan=STEP,smooth=SOLVER_DEFAULT, data=((0.0, 1.0), (1.0, 1.0)))  
            
            # Set surface film coefficient for all surfaces.
            model.FilmCondition(name='Int-film',createStepName=firstStep, surface=a.surfaces['allSurfaces'], 
                                    definition=EMBEDDED_COEFF, filmCoeff=filmCoeff, filmCoeffAmplitude='Amp-1', 
                                    sinkTemperature=Tambient, sinkAmplitude='Amp-1', sinkDistributionType=UNIFORM, 
                                    sinkFieldName='')             
            
            # Set predefined temperature field and heat loading for all cells.
            model.Temperature(name='Predefined Field-1',createStepName='Initial',region=a.sets['allCells'],
                distributionType=UNIFORM,crossSectionDistribution=CONSTANT_THROUGH_THICKNESS, magnitudes=(Tambient, ))
            model.BodyHeatFlux(name='Load-1-bodyFluxHeat',createStepName=firstStep, region=a.sets['allCells'], 
                magnitude=1.0, distributionType=USER_DEFINED) 
                
             
        for indx,model in enumerate(self.modelList, 1):     
            # Backup models.
            bakName = model.name + '-bak'
            if bakName not in list(mdb.models.keys()):
                mdb.Model(name=bakName, objectToCopy=model)
            
            # Create jobs.
            mdb.Job(name='Job-{}'.format(indx), model=model, description='', type=ANALYSIS, 
                atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=70, memoryUnits=PERCENTAGE, 
                getMemoryFromAnalysis=True, explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, 
                echoPrint=OFF, modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, 
                userSubroutine="DFLUX-{}.for".format(model.name), scratch='', resultsFormat=ODB, 
                multiprocessingMode=DEFAULT, numCpus=2, numDomains=2, numGPUs=0) 
                
    
    def display(self):
        for name in list(self.weldSeq.keys()):
            weldList = self.weldSeq[name] 
            for indx,weld in enumerate(weldList, 1):
                displayView(weld, arrowName='{}-Arrow-{}'.format(weld.modelName,indx), arrowLength=15.0, arrowColor='#0000FF', 
                            headStyle=HOLLOW_CIRCLE, textOffset=(0,0), textName='{}-Text-{}'.format(weld.modelName, indx), 
                            textContent='WeldLine-{}'.format(indx), textColor='#FF0000', bgColor='#CCCCCC', fontSize=14)  
    

# Test	
if __name__ == '__main__':	
    session.journalOptions.setValues(replayGeometry=INDEX,recoverGeometry=INDEX)
    
    # --------------------------------------------------------------------------------     
    p1 = mdb.models['tee'].parts['tee']
    v1, n1, f1, e1 = p1.vertices, p1.nodes, p1.faces, p1.edges   

    p2 = mdb.models['sweep'].parts['sweep']
    v2, n2, f2, e2 = p2.vertices, p2.nodes, p2.faces, p2.edges     
    
    #               modelName  | partName  |  crossSect     |  curveType      | startPoint |  alongPoint |  toePoint | segNumber | weldVoltage | weldCurrent | weldVelocty | weldYita | coolTime  |  deltTemp  |    heatType         |    heatPara
    weldParameters = [{'modelName':'tee',    'partName':'tee',    'crossFace':   (f1[32],),    'curveType':'ArbitraryCurve', 'startPoint':n1[2271],'alongPoint': n1[2251],'toePoint':n1[2229], 'segNumber': 4, 'weldVoltage': 10.2,  'weldCurrent': 300.0,'weldVelocty': 5.0,  'weldYita': 0.8, 'coolTime':2.0, 'deltTemp': 1500,'heatType':'DoubleEllipsoid', 'heatPara':(2.0, 4.0, 0.9, 3.0) },
                      {'modelName':'sweep',  'partName':'sweep',  'crossFace':(f2[22],f2[23]), 'curveType':'Circle'        , 'startPoint':n2[0],   'alongPoint': n2[96],  'toePoint':n2[6],    'segNumber': 3, 'weldVoltage': 10.2,  'weldCurrent': 300.0,'weldVelocty': 5.0,  'weldYita': 0.8, 'coolTime':2.0, 'deltTemp': 1500,'heatType':'DoubleEllipsoid', 'heatPara':(3.0, 4.0, 0.9, 3.0) },
                      {'modelName':'sweep',  'partName':'sweep',  'crossFace':(f2[9],f2[10]),  'curveType':'StraightLine'  , 'startPoint':n2[26],  'alongPoint': n2[407], 'toePoint':n2[61],   'segNumber': 1, 'weldVoltage': 10.2,  'weldCurrent': 300.0,'weldVelocty': 5.0,  'weldYita': 0.8, 'coolTime':2.0, 'deltTemp': 1500,'heatType':'DoubleEllipsoid', 'heatPara':(9.0, 4.0, 0.9, 3.0) }]
   
    ws = WeldSequence(weldParameters)  
