C=====================================================================
C          This subroutine is generated by ABAQUS weld plugin        C
C                 on 2024-06-05 16:04:20                             C       
C                    CopyRight © 2024 ChenXuan                       C
C                   Email: isChenXuan@outlook.com                    C
C            本子程序用于生成沿任意空间曲线运动的热源函数.           C
C=====================================================================
      INCLUDE 'MODULES.INC'                                           
      SUBROUTINE DFLUX(FLUX,SOL,KSTEP,KINC,TIME,NOEL,NPT,COORDS,JLTYP,
     1TEMP,PRESS,SNAME)                                               
      USE WELDMATH                                                    
      INCLUDE 'ABA_PARAM.INC'                                         
      DIMENSION COORDS(3),FLUX(2),TIME(2)                             
      CHARACTER*80 SNAME                                              
C   以上声明除INCLUDE 'MODULES.INC'及USE WELDMATH外,均为ABAQUS设置的接
C   口参数,修改就会报错,勿动!!!                                       
                                                                      
C   以下为用户变量声明                                                
C   焊接工艺参数                                                      
      ! WI：焊接电流,安培                                             
      ! WU：焊接电压,伏特                                             
      ! YI：焊接热效率,<1.0                                           
      ! VE：焊接线速度,毫米/秒                                        
      REAL :: WI, WU, YI, VE                                          
                                                                      
C   双椭球热源参数                                                    
      ! a1  ：热源前进方向椭球半轴长,毫米,建议a1=b=0.45*焊缝宽度      
      ! a2  ：热源后方方向椭球半轴长,毫米,建议a2=2a1                  
      ! b   ：焊缝宽度方向椭球半轴长,毫米,建议a1=b=0.45*焊缝宽度      
      ! c   ：焊缝深度方向椭球半轴长,毫米,建议c=0.90*焊缝深度         
      ! POWR：输入功率,毫瓦                                           
      ! f1  ：前半椭球能量系数,f1+f2=2,f1/f2=a1/a2                    
      ! f2  ：后半椭球能量系数                                        
      ! qm1 ：前半椭球能量峰值                                        
      ! qm2 ：后半椭球能量峰值                                        
      REAL :: a1, a2, b, c, POWR, f1, f2, qm1, qm2                    
                                                                      
C   坐标计算                                                          
      ! P1  : 焊接点                                                  
      ! P2  : 焊接前进点                                              
      ! P3  : 焊趾点                                                  
      ! DIS ：点到热源面的距离                                        
      REAL :: P1(3), P2(3), P3(3)                                     
      REAL :: DIS(3)                                                  
                                                                      
C   接口参数                                                          
      ! TCLD：本条焊缝开始之前的累积消耗时间,需要覆盖本条焊缝之前的所 
      !       有焊接及非焊接时间.                                     
      ! TCRT：本条焊缝已焊接时间,TCRT=TIME(2)-TCLD,若本条焊缝在一个   
      !       STEP内完成,则TCRT=TIME(1).                              
      ! PA  ：任意空间点坐标,即ABAQUS接口中的COORDS,由于COORDS参数无  
      !       REAL修饰,直接调用会报错.                                
      REAL :: TCLD, TCRT, PA(3)                                       
                                                                      
C  ========================第1条焊缝开始============================                                      
      IF (KSTEP >= 2 .AND. KSTEP <= 28) THEN                                  
C  焊接参数设置                                                                                                          
          WI = 300.0                                                                             
          WU = 9.5                                                                             
          YI = 0.9                                                                                
          VE = 10.0                                                                             
C  双椭球热源参数设置                                                                                                    
          ! 输入值                                                                                                       
          a1 = 3.2                                                                                      
          a2 = 6.4                                                                                      
          b  = 3.2                                                                                       
          c  = 3.2                                                                                       
          qm1  = 97392.9748535                                                                                   
          qm2  = 97392.9748535                                                                                   
C  计算坐标                                                                                                              
          TCRT = TIME(2) - 1e-07                                                                               
          RW = 12.4651                                                                                         
          RT = 10.0                                                                                         
          OMEG = 0.8022                                                                                    
          ! 局部坐标                                                                                                   
          XWL = RW*COS(OMEG*TCRT)                                                                                      
          YWL = RW*SIN(OMEG*TCRT)                                                                                      
          ZWL = 0.0                                                                                                    
          XFL = XWL - RW*OMEG*SIN(OMEG*TCRT)                                                                           
          YFL = YWL + RW*OMEG*COS(OMEG*TCRT)                                                                           
          ZFL = 0.0                                                                                                    
          XTL = RT*COS(OMEG*TCRT)                                                                                      
          YTL = RT*SIN(OMEG*TCRT)                                                                                      
          ZTL = 0.0                                                                                                    
          ! 全局坐标                                                                                                   
          XW = -19.976+0.2588*XWL+-0.9659*YWL+0.0*ZWL        
          YW = -12.4651+0.0*XWL+0.0*YWL+1.0*ZWL        
          ZW = 40.4006+-0.9659*XWL+-0.2588*YWL+0.0*ZWL        
          XF = -19.976+0.2588*XFL+-0.9659*YFL+0.0*ZFL        
          YF = -12.4651+0.0*XFL+0.0*YFL+1.0*ZFL        
          ZF = 40.4006+-0.9659*XFL+-0.2588*YFL+0.0*ZFL        
          XT = -19.976+0.2588*XTL+-0.9659*YTL+0.0*ZTL        
          YT = -10.0+0.0*XTL+0.0*YTL+1.0*ZTL        
          ZT = 40.4006+-0.9659*XTL+-0.2588*YTL+0.0*ZTL        
                                                                                                                            
          P1 = [XW, YW, ZW]                                                                                                 
          P2 = [XF, YF, ZF]                                                                                                 
          P3 = [XT, YT, ZT]                                                                                                 
                                                                                                                            
C  计算热源                                                                                                                
          PA = COORDS                                                                                                      
          DIS = DISTANCE(P1, P2, P3, PA)                                                                                   
          xx = DIS(1)**2                                                                                                   
          yy = DIS(2)**2                                                                                                   
          zz = DIS(3)**2                                                                                                   
          IF (DIS(1) .GE. 0.0) THEN                                                                                        
              FLUX(1) = qm1*EXP(-3*(xx/a1**2 + yy/b**2 + zz/c**2))                                                         
          ELSE                                                                                                             
              FLUX(1) = qm2*EXP(-3*(xx/a2**2 + yy/b**2 + zz/c**2))                                                         
          END IF                                                                                                           
C  ========================第1条焊缝结束============================                                        
C  ========================第2条焊缝开始============================                                      
      ELSEIF (KSTEP >= 30 .AND. KSTEP <= 54) THEN                              
C  焊接参数设置                                                                                                          
          WI = 300.0                                                                             
          WU = 9.5                                                                             
          YI = 0.9                                                                                
          VE = 10.0                                                                             
C  双椭球热源参数设置                                                                                                    
          ! 输入值                                                                                                       
          a1 = 3.2                                                                                      
          a2 = 6.4                                                                                      
          b  = 3.2                                                                                       
          c  = 3.2                                                                                       
          qm1  = 97392.9748535                                                                                   
          qm2  = 97392.9748535                                                                                   
C  计算坐标                                                                                                              
          TCRT = TIME(2) - 5.87404254933                                                                               
          XW0 = 27.1348                                                                     
          YW0 = -12.4651                                                                     
          ZW0 = 42.0687                                                                     
          XT0 = 25.0                                                                     
          YT0 = -10.0                                                                     
          ZT0 = 43.3013                                                                     
          DX = -0.5                                                                       
          DY = 0.0                                                                       
          DZ = -0.866                                                                       
          VX = -5.0                                                                       
          VY = 0.0                                                                       
          VZ = -8.66                                                                       
          XW = XW0 + VX*TCRT                                                                         
          YW = YW0 + VY*TCRT                                                                         
          ZW = ZW0 + VZ*TCRT                                                                         
          XF = XW + DX                                                                               
          YF = YW + DY                                                                               
          ZF = ZW + DZ                                                                               
          XT = XT0 + VX*TCRT                                                                         
          YT = YT0 + VY*TCRT                                                                         
          ZT = ZT0 + VZ*TCRT                                                                         
                                                                                                                            
          P1 = [XW, YW, ZW]                                                                                                 
          P2 = [XF, YF, ZF]                                                                                                 
          P3 = [XT, YT, ZT]                                                                                                 
                                                                                                                            
C  计算热源                                                                                                                
          PA = COORDS                                                                                                      
          DIS = DISTANCE(P1, P2, P3, PA)                                                                                   
          xx = DIS(1)**2                                                                                                   
          yy = DIS(2)**2                                                                                                   
          zz = DIS(3)**2                                                                                                   
          IF (DIS(1) .GE. 0.0) THEN                                                                                        
              FLUX(1) = qm1*EXP(-3*(xx/a1**2 + yy/b**2 + zz/c**2))                                                         
          ELSE                                                                                                             
              FLUX(1) = qm2*EXP(-3*(xx/a2**2 + yy/b**2 + zz/c**2))                                                         
          END IF                                                                                                           
C  ========================第2条焊缝结束============================                                        
      ELSE                                                                                                               
          FLUX(1) = 0.0                                                                                                  
      END IF                                                                                                             
                                                                                                                         
      RETURN                                                                                                             
      END                                                                                                                
