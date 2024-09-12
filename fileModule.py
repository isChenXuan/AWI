# --*coding=UTF-8*--
from abaqus import *
from csv import *
from datetime import *
import os

     
class FileModule():
    
    def __init__ (self, weldSequence):  
        self.weldSequence = weldSequence  
        
        # Call functions.
        self.writeCsvs()
        self.writeMODULES()  
        self.writeDFLUX() 

    def fileStatus(self, fileName):
        """
        # Checking whether file exists. If yes, let user choose to overwrite file or not.
        """ 
        file = "./" + fileName
        if os.path.exists(file) == True:
            reply = getWarningReply("The file {} exists!\nOverwrite?".format(fileName), (YES, NO))
            if reply == YES:
                return True
            elif reply == NO:
                return False
        else:
            return True
            
        
    def writeCsvs(self):
        """
        # Write weld nodes and elements info to csv file.
        """ 
        dir = os.getcwd()
        
        for name in list(self.weldSequence.weldSeq.keys()):
            weldList = self.weldSequence.weldSeq[name]
            for indx,weld in enumerate(weldList, 1):
                # write node list csv.
                fileName = "{}-weldLineNodes-{}.csv".format(name, indx)
                with open(fileName, 'w') as f:
                    f.write("No.,weldNodeLabel,x,y,z,distance,time\n")
                    for i in range(len(weld._WeldLine__weldNodeLabels)):
                        csv_line = "{},{},{},{},{},{},{}\n".format(i+1, weld._WeldLine__weldNodeLabels[i], 
                                                                        weld._WeldLine__weldNodeCoords[i][0], 
                                                                        weld._WeldLine__weldNodeCoords[i][1], 
                                                                        weld._WeldLine__weldNodeCoords[i][2], 
                                                                        weld._WeldLine__length[i], 
                                                                        weld._WeldLine__time[i])
                        f.write(csv_line)
                #print "Coordinates of weld line nodes have been written into file {} successfully!".format(dir+"\\"+fileName)  
                #os.remove(fileName)
                
                fileName = "{}-toeLineNodes-{}.csv".format(name, indx)
                with open(fileName, 'w') as f:
                    f.write("No.,toeNodeLabel,x,y,z\n")
                    for i in range(len(weld._WeldLine__toeNodeLabels)):
                        csv_line = "{},{},{},{},{}\n".format(i+1, weld._WeldLine__toeNodeLabels[i], 
                                                                  weld._WeldLine__toeNodeCoords[i][0], 
                                                                  weld._WeldLine__toeNodeCoords[i][1], 
                                                                  weld._WeldLine__toeNodeCoords[i][2])
                        f.write(csv_line)
                #print "Coordinates of toe line nodes have been written into file {} successfully!".format(dir+"\\"+fileName)         
                #os.remove(fileName)
                
                # write element list csv.
                fileName = "{}-rowOfElements-{}.csv".format(name, indx)
                with open(fileName, 'w') as f:
                    for i in range(len(weld._WeldLine__elementLabels)):
                        f.write("the {}th sweep path".format(i+1))
                        f.write("No.,elemLabel\n")
                        for no, label in enumerate(weld._WeldLine__elementLabels[i], 1):
                            csv_line = "{},{}\n".format(no, label)
                            f.write(csv_line)
                #print "Row of elements has been written in file {} successfully!".format(dir+"\\"+fileName) 
                #os.remove(fileName)
                
        
    def writeMODULES(self):
        """
        # Write FORTRAN math modules file used in abaqus subroutine. Shouldn't be modified!!!
        """ 
        with open("./MODULES.INC", 'w') as f:
            f.writelines("      MODULE WELDMATH                                                       \n")
            f.writelines("      CONTAINS                                                              \n")
            f.writelines("C   本模块包含了热源子程序计算使用到的函数                                  \n")
            f.writelines("C   !!! 以下向量默认都是行向量 !!!                                          \n")
            f.writelines("                                                                            \n")
            f.writelines("C   向量点乘V1·V2                                                           \n")
            f.writelines("      FUNCTION VDOT(V1, V2) RESULT(DOT)                                     \n")
            f.writelines("          REAL:: V1(3), V2(3), DOT                                          \n")
            f.writelines("          DOT = V1(1)*V2(1)+V1(2)*V2(2)+V1(3)*V2(3)                         \n")
            f.writelines("      END FUNCTION VDOT                                                     \n")
            f.writelines("                                                                            \n")
            f.writelines("C   向量叉乘V1×V2                                                           \n")
            f.writelines("      FUNCTION VCROSS(V1, V2) RESULT(V3)                                    \n")
            f.writelines("          REAL:: V1(3), V2(3), V3(3)                                        \n")
            f.writelines("          V3(1) = V1(2)*V2(3)-V1(3)*V2(2)                                   \n")
            f.writelines("          V3(2) = V1(3)*V2(1)-V1(1)*V2(3)                                   \n")
            f.writelines("          V3(3) = V1(1)*V2(2)-V1(2)*V2(1)                                   \n")
            f.writelines("      END FUNCTION VCROSS                                                   \n")
            f.writelines("                                                                            \n")
            f.writelines("C   向量的模                                                                \n")
            f.writelines("      FUNCTION VNORM(V) RESULT(VM)                                          \n")
            f.writelines("          REAL:: V(3), VM                                                   \n")
            f.writelines("          VM = SQRT(VDOT(V, V))                                             \n")
            f.writelines("      END FUNCTION VNORM                                                    \n")
            f.writelines("                                                                            \n")
            f.writelines("C   坐标变换                                                                \n")
            f.writelines("    ! r = Xi*Ei = xi*ei                                                     \n")
            f.writelines("    ! 旧坐标Xi,新坐标xi,均为行向量                                          \n")
            f.writelines("    ! 旧坐标基Ei,新坐标基ei                                                 \n")
            f.writelines("    ! =>  Xi*Ei·ej = xi*ei·ej = xj                                          \n")
            f.writelines("    ! 或  xi*ei·Ej = Xi*Ei·Ej = Xj                                          \n")
            f.writelines("    !                  | E1*e1  E1*e2  E1*e3 |                              \n")
            f.writelines("    ! 令RMTX = Ei·ej = | E2*e1  E2*e2  E2*e3 |                              \n")
            f.writelines("    !                  | E3*e1  E3*e2  E3*e3 |                              \n")
            f.writelines("    ! 则x = X*RMTX,或X = x*RMTX逆 = x*RMTX转置                              \n")
            f.writelines("    !                                        | E1 |                         \n")
            f.writelines("    ! 特别地,若新坐标系为全局坐标架,则RMTX = | E2 |                         \n")
            f.writelines("    !                                        | E3 |                         \n") 
            f.writelines("      FUNCTION PROT(RMTX, POLD) RESULT(PNEW)                                \n")
            f.writelines("          REAL:: RMTX(3, 3), POLD(3), PNEW(3)                               \n")
            f.writelines("          PNEW(1) = VDOT(POLD, RMTX(:,1))                                   \n")
            f.writelines("          PNEW(2) = VDOT(POLD, RMTX(:,2))                                   \n")
            f.writelines("          PNEW(3) = VDOT(POLD, RMTX(:,3))                                   \n")
            f.writelines("      END FUNCTION PROT                                                     \n")
            f.writelines("                                                                            \n")
            f.writelines("C   设置局部坐标系单位向量                                                  \n")
            f.writelines("    ! P1为局部坐标原点,P1P2为局部1轴,局部3轴为P1P2×P1P3                     \n")
            f.writelines("      FUNCTION SETLOC(P1, P2, P3) RESULT(ELOC)                              \n")
            f.writelines("          REAL:: P1(3), P2(3), P3(3), ELOC(3, 3)                            \n")
            f.writelines("          ELOC(1,:) = (P2-P1)/VNORM(P2-P1)                                  \n")
            f.writelines("          ELOC(3,:) = VCROSS(P2-P1,P3-P1)/VNORM(VCROSS(P2-P1,P3-P1))        \n")
            f.writelines("          ELOC(2,:) = VCROSS(ELOC(3,:), ELOC(1,:))                          \n")
            f.writelines("      END FUNCTION SETLOC                                                   \n")
            f.writelines("                                                                            \n")
            f.writelines("C   计算任意点PA与当前时刻热源平面的距离                                    \n")
            f.writelines("    ! P1为当前时刻焊接点,P2为焊接前进方向点,P3为焊趾点,PA为空间中任意一点   \n")
            f.writelines("      FUNCTION DISTANCE (P1, P2, P3, PA) RESULT(DIS)                        \n")
            f.writelines("          REAL:: P1(3), P2(3), P3(3), PA(3)                                 \n")
            f.writelines("          REAL:: DIS(3)                                                     \n")
            f.writelines("          REAL:: VX(3), VY(3), VZ(3)                                        \n")
            f.writelines("          REAL:: VTMP(3, 3)                  ! 临时变量                     \n")
            f.writelines("          VTMP = SETLOC(P1, P2, P3)          ! 以P1为原点,搭建局部直角坐标系\n")
            f.writelines("          VX = VTMP(1,:)                     ! 焊接前进方向单位向量         \n")
            f.writelines("          VY = VTMP(2,:)                     ! 焊缝宽度方向单位向量         \n")
            f.writelines("          VZ = VTMP(3,:)                     ! 焊缝深度方向单位向量         \n")
            f.writelines("    ! 计算点PA到热源平面的距离(带符号)：d = |P1PA|*COS<V,P1PA> = V·P1PA/|V| \n")
            f.writelines("          DIS(1) = VDOT(VX, PA-P1)                                          \n")
            f.writelines("          DIS(2) = VDOT(VY, PA-P1)                                          \n")
            f.writelines("          DIS(3) = VDOT(VZ, PA-P1)                                          \n")
            f.writelines("      END FUNCTION DISTANCE                                                 \n")
            f.writelines("                                                                            \n")
            f.writelines("      END MODULE WELDMATH                                                   \n")  

            
    def writeDFLUX(self):
        """
        # Write FORTRAN dflux subroutine file used in abaqus.
        """ 
        now = datetime.now()
        tim = now.strftime("%Y-%m-%d %H:%M:%S")
        dir = os.getcwd()
        
        for name in list(self.weldSequence.weldSeq.keys()):
            fileName = "DFLUX-{}.for".format(name)
            weldList = self.weldSequence.weldSeq[name]
            stepList = self.weldSequence.stepSeq[name]
            timeList = self.weldSequence.timeSeq[name]
            if self.fileStatus(fileName):    
                with open(fileName, 'w') as f:
                    # Write file head and parameters definition.
                    f.writelines("C=====================================================================\n")
                    f.writelines("C          This subroutine is generated by ABAQUS weld plugin        C\n")
                    f.writelines("C                 on " + tim +  "                             C       \n")
                    f.writelines("C                    CopyRight © 2024 ChenXuan                       C\n")
                    f.writelines("C                   Email: isChenXuan@outlook.com                    C\n")
                    f.writelines("C            本子程序用于生成沿任意空间曲线运动的热源函数.           C\n")
                    f.writelines("C=====================================================================\n")
                    f.writelines("      INCLUDE 'MODULES.INC'                                           \n")
                    f.writelines("      SUBROUTINE DFLUX(FLUX,SOL,KSTEP,KINC,TIME,NOEL,NPT,COORDS,JLTYP,\n")
                    f.writelines("     1TEMP,PRESS,SNAME)                                               \n")
                    f.writelines("      USE WELDMATH                                                    \n")
                    f.writelines("      INCLUDE 'ABA_PARAM.INC'                                         \n")
                    f.writelines("      DIMENSION COORDS(3),FLUX(2),TIME(2)                             \n")
                    f.writelines("      CHARACTER*80 SNAME                                              \n")
                    f.writelines("C   以上声明除INCLUDE 'MODULES.INC'及USE WELDMATH外,均为ABAQUS设置的接\n")
                    f.writelines("C   口参数,修改就会报错,勿动!!!                                       \n")
                    f.writelines("                                                                      \n")        
                    f.writelines("C   以下为用户变量声明                                                \n") 
                    f.writelines("C   焊接工艺参数                                                      \n") 
                    f.writelines("      ! WI：焊接电流,安培                                             \n")
                    f.writelines("      ! WU：焊接电压,伏特                                             \n")
                    f.writelines("      ! YI：焊接热效率,<1.0                                           \n") 
                    f.writelines("      ! VE：焊接线速度,毫米/秒                                        \n")                      
                    f.writelines("      REAL :: WI, WU, YI, VE                                          \n")
                    f.writelines("                                                                      \n")            
                    f.writelines("C   双椭球热源参数                                                    \n")
                    f.writelines("      ! a1  ：热源前进方向椭球半轴长,毫米,建议a1=b=0.45*焊缝宽度      \n")
                    f.writelines("      ! a2  ：热源后方方向椭球半轴长,毫米,建议a2=2a1                  \n")
                    f.writelines("      ! b   ：焊缝宽度方向椭球半轴长,毫米,建议a1=b=0.45*焊缝宽度      \n")
                    f.writelines("      ! c   ：焊缝深度方向椭球半轴长,毫米,建议c=0.90*焊缝深度         \n")
                    f.writelines("      ! POWR：输入功率,毫瓦                                           \n")
                    f.writelines("      ! f1  ：前半椭球能量系数,f1+f2=2,f1/f2=a1/a2                    \n")
                    f.writelines("      ! f2  ：后半椭球能量系数                                        \n")
                    f.writelines("      ! qm1 ：前半椭球能量峰值                                        \n")
                    f.writelines("      ! qm2 ：后半椭球能量峰值                                        \n")
                    f.writelines("      REAL :: a1, a2, b, c, POWR, f1, f2, qm1, qm2                    \n")
                    f.writelines("                                                                      \n")            
                    f.writelines("C   坐标计算                                                          \n")
                    f.writelines("      ! P1  : 焊接点                                                  \n")
                    f.writelines("      ! P2  : 焊接前进点                                              \n")
                    f.writelines("      ! P3  : 焊趾点                                                  \n")
                    f.writelines("      ! DIS ：点到热源面的距离                                        \n")
                    f.writelines("      REAL :: P1(3), P2(3), P3(3)                                     \n")
                    f.writelines("      REAL :: DIS(3)                                                  \n")
                    f.writelines("                                                                      \n")
                    f.writelines("C   接口参数                                                          \n")
                    f.writelines("      ! TCLD：本条焊缝开始之前的累积消耗时间,需要覆盖本条焊缝之前的所 \n")
                    f.writelines("      !       有焊接及非焊接时间.                                     \n")
                    f.writelines("      ! TCRT：本条焊缝已焊接时间,TCRT=TIME(2)-TCLD,若本条焊缝在一个   \n")
                    f.writelines("      !       STEP内完成,则TCRT=TIME(1).                              \n")
                    f.writelines("      ! PA  ：任意空间点坐标,即ABAQUS接口中的COORDS,由于COORDS参数无  \n")
                    f.writelines("      !       REAL修饰,直接调用会报错.                                \n")
                    f.writelines("      REAL :: TCLD, TCRT, PA(3)                                       \n")
                    f.writelines("                                                                      \n") 

                    # Write KSTEP.
                    for i in range(len(weldList)):
                        weld = weldList[i]
                        kstep = stepList[i]
                        f.writelines("C  ========================第" + str(i+1) + "条焊缝开始============================                                      \n")                        
                        if i == 0:
                            f.writelines("      IF (KSTEP >= " + str(kstep[0]) + " .AND. KSTEP <= " + str(kstep[1]) + ") THEN                                  \n")
                        else:                                                                                                                                  
                            f.writelines("      ELSEIF (KSTEP >= " + str(kstep[0]) + " .AND. KSTEP <= " + str(kstep[1]) + ") THEN                              \n")
                        f.writelines("C  焊接参数设置                                                                                                          \n")
                        f.writelines("          WI = " + str(weld.weldCurrent) + "                                                                             \n")
                        f.writelines("          WU = " + str(weld.weldVoltage) + "                                                                             \n")
                        f.writelines("          YI = " + str(weld.weldYita) + "                                                                                \n")
                        f.writelines("          VE = " + str(weld.weldVelocty) + "                                                                             \n")
                        f.writelines("C  双椭球热源参数设置                                                                                                    \n")
                        f.writelines("          ! 输入值                                                                                                       \n")
                        f.writelines("          a1 = " + str(weld.a1) + "                                                                                      \n")
                        f.writelines("          a2 = " + str(weld.a2) + "                                                                                      \n")
                        f.writelines("          b  = " + str(weld.b) + "                                                                                       \n")
                        f.writelines("          c  = " + str(weld.c) + "                                                                                       \n")
                        f.writelines("          qm1  = " + str(weld.qm1) + "                                                                                   \n")
                        f.writelines("          qm2  = " + str(weld.qm2) + "                                                                                   \n")
                        f.writelines("C  计算坐标                                                                                                              \n") 
                        
                        # Time spent before this weld line.
                        Ts = timeList[i] 
                        # TCRT = TIME(2) - Ts, i.e., current time.
                        f.writelines("          TCRT = TIME(2) - " + str(Ts) + "                                                                               \n") 
                        
                        # Write coordinates.
                        if weld.curveType == 'StraightLine':
                            self.writeStraightLine(weld, f)
                        elif weld.curveType == 'Circle':
                            self.writeCircle(weld, f)
                        elif weld.curveType == 'ArbitraryCurve':
                            self.writeArbitraryCurve(weld, f)
                            
                        f.writelines("                                                                                                                            \n")
                        f.writelines("          P1 = [XW, YW, ZW]                                                                                                 \n") 
                        f.writelines("          P2 = [XF, YF, ZF]                                                                                                 \n") 
                        f.writelines("          P3 = [XT, YT, ZT]                                                                                                 \n") 
                        f.writelines("                                                                                                                            \n")
                        
                        # Calculate heat.
                        f.writelines("C  计算热源                                                                                                                \n")
                        f.writelines("          PA = COORDS                                                                                                      \n")
                        f.writelines("          DIS = DISTANCE(P1, P2, P3, PA)                                                                                   \n")
                        f.writelines("          xx = DIS(1)**2                                                                                                   \n")
                        f.writelines("          yy = DIS(2)**2                                                                                                   \n")
                        f.writelines("          zz = DIS(3)**2                                                                                                   \n")
                        f.writelines("          IF (DIS(1) .GE. 0.0) THEN                                                                                        \n")
                        f.writelines("              FLUX(1) = qm1*EXP(-3*(xx/a1**2 + yy/b**2 + zz/c**2))                                                         \n")
                        f.writelines("          ELSE                                                                                                             \n")
                        f.writelines("              FLUX(1) = qm2*EXP(-3*(xx/a2**2 + yy/b**2 + zz/c**2))                                                         \n")
                        f.writelines("          END IF                                                                                                           \n")
                        f.writelines("C  ========================第" + str(i+1) + "条焊缝结束============================                                        \n")
                     
                    # Write file tail. 
                    f.writelines("      ELSE                                                                                                               \n")
                    f.writelines("          FLUX(1) = 0.0                                                                                                  \n")
                    f.writelines("      END IF                                                                                                             \n")
                    f.writelines("                                                                                                                         \n")
                    f.writelines("      RETURN                                                                                                             \n")
                    f.writelines("      END                                                                                                                \n")
                
                print "DFLUX subroutine for model {} has been written in file {} successfully!".format(name, dir+"\\"+fileName) 
                
    
    def writeStraightLine(self, weld, f):
        to4f = lambda x : float(format(x, '.4f'))
        XW0 = to4f(weld._WeldLine__weldNodeCoords[0][0])
        YW0 = to4f(weld._WeldLine__weldNodeCoords[0][1])
        ZW0 = to4f(weld._WeldLine__weldNodeCoords[0][2])
             
        XT0 = to4f(weld._WeldLine__toeNodeCoords[0][0])
        YT0 = to4f(weld._WeldLine__toeNodeCoords[0][1])
        ZT0 = to4f(weld._WeldLine__toeNodeCoords[0][2])
        
        DX = to4f(weld.iniVector[0])
        DY = to4f(weld.iniVector[1])
        DZ = to4f(weld.iniVector[2])

        VX = weld.weldVelocty*DX
        VY = weld.weldVelocty*DY
        VZ = weld.weldVelocty*DZ
        
        f.writelines("          XW0 = " + str(XW0) + "                                                                     \n")
        f.writelines("          YW0 = " + str(YW0) + "                                                                     \n")
        f.writelines("          ZW0 = " + str(ZW0) + "                                                                     \n")
                                                                                                                           
        f.writelines("          XT0 = " + str(XT0) + "                                                                     \n")
        f.writelines("          YT0 = " + str(YT0) + "                                                                     \n")
        f.writelines("          ZT0 = " + str(ZT0) + "                                                                     \n")
                                                                                                                           
        f.writelines("          DX = " + str(DX) + "                                                                       \n")
        f.writelines("          DY = " + str(DY) + "                                                                       \n")
        f.writelines("          DZ = " + str(DZ) + "                                                                       \n")
                                                                                                                           
        f.writelines("          VX = " + str(VX) + "                                                                       \n")
        f.writelines("          VY = " + str(VY) + "                                                                       \n")
        f.writelines("          VZ = " + str(VZ) + "                                                                       \n")
        
        # Weld point.
        f.writelines("          XW = XW0 + VX*TCRT                                                                         \n")
        f.writelines("          YW = YW0 + VY*TCRT                                                                         \n")
        f.writelines("          ZW = ZW0 + VZ*TCRT                                                                         \n")
                                                                                                                                                   
        # Front point.                                                                                                                             
        f.writelines("          XF = XW + DX                                                                               \n")
        f.writelines("          YF = YW + DY                                                                               \n")
        f.writelines("          ZF = ZW + DZ                                                                               \n")
                                                                                                                                                   
        # Toe point.                                                                                                                               
        f.writelines("          XT = XT0 + VX*TCRT                                                                         \n")
        f.writelines("          YT = YT0 + VY*TCRT                                                                         \n")
        f.writelines("          ZT = ZT0 + VZ*TCRT                                                                         \n")       
                                                                                                                                                   

    def writeCircle(self, weld, f):
        to4f = lambda x : float(format(x, '.4f'))
        RW = to4f(weld.weldCircle.RA)
        RT = to4f(weld.toeCircle.RA )
        omiga = to4f(weld.weldVelocty/RW)
        f.writelines("          RW = " + str(RW) + "                                                                                         \n")
        f.writelines("          RT = " + str(RT) + "                                                                                         \n")
        f.writelines("          OMEG = " + str(omiga) + "                                                                                    \n")
        # Write local coordinates.
        f.writelines("          ! 局部坐标                                                                                                   \n")
        # Weld point.
        f.writelines("          XWL = RW*COS(OMEG*TCRT)                                                                                      \n")
        f.writelines("          YWL = RW*SIN(OMEG*TCRT)                                                                                      \n")
        f.writelines("          ZWL = 0.0                                                                                                    \n")
        
        # Front point.
        f.writelines("          XFL = XWL - RW*OMEG*SIN(OMEG*TCRT)                                                                           \n")
        f.writelines("          YFL = YWL + RW*OMEG*COS(OMEG*TCRT)                                                                           \n")
        f.writelines("          ZFL = 0.0                                                                                                    \n") 
        
        # Toe point.
        f.writelines("          XTL = RT*COS(OMEG*TCRT)                                                                                      \n")
        f.writelines("          YTL = RT*SIN(OMEG*TCRT)                                                                                      \n")
        f.writelines("          ZTL = 0.0                                                                                                    \n") 
        
        # Transfer to global coordinates.
        f.writelines("          ! 全局坐标                                                                                                   \n")
        # Local to global tansformation matrix.
        MTXW = weld.weldCircle.local2Global
        MTXT = weld.toeCircle.local2Global
        # Local coordinates origin point.
        XOW = to4f(weld.weldCircle.PC[0])
        YOW = to4f(weld.weldCircle.PC[1])
        ZOW = to4f(weld.weldCircle.PC[2])
        XOT = to4f(weld.toeCircle.PC[0])
        YOT = to4f(weld.toeCircle.PC[1])
        ZOT = to4f(weld.toeCircle.PC[2])
        # Weld point.
        f.writelines("          XW = " + str(XOW) + "+" + str(MTXW[0][0]) + "*XWL+" + str(MTXW[0][1]) + "*YWL+" + str(MTXW[0][2]) +"*ZWL        \n")
        f.writelines("          YW = " + str(YOW) + "+" + str(MTXW[1][0]) + "*XWL+" + str(MTXW[1][1]) + "*YWL+" + str(MTXW[1][2]) +"*ZWL        \n")
        f.writelines("          ZW = " + str(ZOW) + "+" + str(MTXW[2][0]) + "*XWL+" + str(MTXW[2][1]) + "*YWL+" + str(MTXW[2][2]) +"*ZWL        \n")  
                                                                                                                                                  
        # Front point.                                                                                                                            
        f.writelines("          XF = " + str(XOW) + "+" + str(MTXW[0][0]) + "*XFL+" + str(MTXW[0][1]) + "*YFL+" + str(MTXW[0][2]) +"*ZFL        \n")
        f.writelines("          YF = " + str(YOW) + "+" + str(MTXW[1][0]) + "*XFL+" + str(MTXW[1][1]) + "*YFL+" + str(MTXW[1][2]) +"*ZFL        \n")
        f.writelines("          ZF = " + str(ZOW) + "+" + str(MTXW[2][0]) + "*XFL+" + str(MTXW[2][1]) + "*YFL+" + str(MTXW[2][2]) +"*ZFL        \n") 
                                                       
        # Toe point.                                   
        f.writelines("          XT = " + str(XOT) + "+" + str(MTXT[0][0]) + "*XTL+" + str(MTXT[0][1]) + "*YTL+" + str(MTXT[0][2]) +"*ZTL        \n")
        f.writelines("          YT = " + str(YOT) + "+" + str(MTXT[1][0]) + "*XTL+" + str(MTXT[1][1]) + "*YTL+" + str(MTXT[1][2]) +"*ZTL        \n")
        f.writelines("          ZT = " + str(ZOT) + "+" + str(MTXT[2][0]) + "*XTL+" + str(MTXT[2][1]) + "*YTL+" + str(MTXT[2][2]) +"*ZTL        \n")        
        

    def writeArbitraryCurve(self, weld, f):
        time = [float(l)/weld.weldVelocty for l in weld._WeldLine__length]
        
        # Weld nodes.
        fw = weld.weldDp.symbolCalculate('TCRT', time, flag=0)
        xw = fw[0]
        yw = fw[1]
        zw = fw[2]
        df = weld.weldDp.symbolCalculate('TCRT', time, flag=1)
        dx = df[0]
        dy = df[1]
        dz = df[2]
        # Toe nodes.
        ft = weld.toeDp.symbolCalculate('TCRT', time, flag=0)
        xt = ft[0]
        yt = ft[1]
        zt = ft[2]
        
        # Write points coordinate.
        for j in range(1, len(time)):
            if j == 1:
                f.writelines("          IF (TCRT <= " + str(float(format(time[j], '.4f'))) + ") THEN                                                 \n")   
            else:                                                                                                                                    
                f.writelines("          ELSEIF (TCRT <= " + str(float(format(time[j], '.4f'))) + ") THEN                                             \n") 
            
            # Weld point.    
            f.writelines("              XW = " + str(xw[j-1]) + "                                                                                    \n")
            f.writelines("              YW = " + str(yw[j-1]) + "                                                                                    \n") 
            f.writelines("              ZW = " + str(zw[j-1]) + "                                                                                    \n") 
                                                                                                                                                     
            # Front point.                                                                                                                           
            f.writelines("              XF = XW+" + str(dx[j-1]) + "                                                                                 \n") 
            f.writelines("              YF = YW+" + str(dy[j-1]) + "                                                                                 \n") 
            f.writelines("              ZF = ZW+" + str(dz[j-1]) + "                                                                                 \n") 
                                                                                                                                                     
            # Toe points                                                                                                                             
            f.writelines("              XT = " + str(xt[j-1]) + "                                                                                    \n") 
            f.writelines("              YT = " + str(yt[j-1]) + "                                                                                    \n") 
            f.writelines("              ZT = " + str(zt[j-1]) + "                                                                                    \n")
            
        f.writelines("          END IF                                                                                                            \n") 

        
# Test	
if __name__ == '__main__':	
    session.journalOptions.setValues(replayGeometry=INDEX,recoverGeometry=INDEX)
    
    # --------------------------------------------------------------------------------     
    p1 = mdb.models['tee'].parts['tee']
    v1, n1, f1, e1 = p1.vertices, p1.nodes, p1.faces, p1.edges   

    p2 = mdb.models['sweep'].parts['sweep']
    v2, n2, f2, e2 = p2.vertices, p2.nodes, p2.faces, p2.edges     
    
    ##               modelName  | partName  |  crossSect |  startPoint |  alongPoint |  toePoint | segNumber | weldVoltage | weldCurrent | weldVelocty | weldYita | coolTime  |  deltTemp      |         heatSource
    #weldParameters = [('tee',    'tee',       (f1[32],),    n1[2271],      n1[2251],    n1[2229],     4,         10.2,          300.0,         5.0,          0.8,      2.0,         1500,   {'DoubleEllipsoid':(2.0, 4.0, 0.9, 3.0)} ),
    #                  ('sweep',  'sweep',  (f2[22],f2[23]),  n2[0],         n2[96],     n2[6],        3,         10.2,          300.0,         5.0,          0.8,      2.0,         1500,   {'DoubleEllipsoid':(3.0, 4.0, 0.9, 3.0)} ),
    #                  ('sweep',  'sweep',  (f2[9],f2[10]),   n2[26],        n2[407],    n2[61],       1,         10.2,          300.0,         5.0,          0.8,      2.0,         1500,   {'DoubleEllipsoid':(9.0, 4.0, 0.9, 3.0)} )]
   
    weldParameters = [{'modelName':'tee',    'partName':'tee',    'crossFace':   (f1[32],),    'curveType':'ArbitraryCurve', 'startPoint':n1[2271],'alongPoint': n1[2251],'toePoint':n1[2229], 'segNumber': 4, 'weldVoltage': 10.2,  'weldCurrent': 300.0,'weldVelocty': 5.0,  'weldYita': 0.8, 'coolTime':2.0, 'deltTemp': 1500,'heatType':'DoubleEllipsoid', 'heatPara':(2.0, 4.0, 0.9, 3.0) },
                      {'modelName':'sweep',  'partName':'sweep',  'crossFace':(f2[22],f2[23]), 'curveType':'Circle'        , 'startPoint':n2[0],   'alongPoint': n2[96],  'toePoint':n2[6],    'segNumber': 3, 'weldVoltage': 10.2,  'weldCurrent': 300.0,'weldVelocty': 5.0,  'weldYita': 0.8, 'coolTime':2.0, 'deltTemp': 1500,'heatType':'DoubleEllipsoid', 'heatPara':(3.0, 4.0, 0.9, 3.0) },
                      {'modelName':'sweep',  'partName':'sweep',  'crossFace':(f2[9],f2[10]),  'curveType':'StraightLine'  , 'startPoint':n2[26],  'alongPoint': n2[407], 'toePoint':n2[61],   'segNumber': 1, 'weldVoltage': 10.2,  'weldCurrent': 300.0,'weldVelocty': 5.0,  'weldYita': 0.8, 'coolTime':2.0, 'deltTemp': 1500,'heatType':'DoubleEllipsoid', 'heatPara':(9.0, 4.0, 0.9, 3.0) }]
   
    ws = WeldSequence(weldParameters)  
    fm = FileModule(ws)   