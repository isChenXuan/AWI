      MODULE WELDMATH                                                       
      CONTAINS                                                              
C   本模块包含了热源子程序计算使用到的函数                                  
C   !!! 以下向量默认都是行向量 !!!                                          
                                                                            
C   向量点乘V1·V2                                                           
      FUNCTION VDOT(V1, V2) RESULT(DOT)                                     
          REAL:: V1(3), V2(3), DOT                                          
          DOT = V1(1)*V2(1)+V1(2)*V2(2)+V1(3)*V2(3)                         
      END FUNCTION VDOT                                                     
                                                                            
C   向量叉乘V1×V2                                                           
      FUNCTION VCROSS(V1, V2) RESULT(V3)                                    
          REAL:: V1(3), V2(3), V3(3)                                        
          V3(1) = V1(2)*V2(3)-V1(3)*V2(2)                                   
          V3(2) = V1(3)*V2(1)-V1(1)*V2(3)                                   
          V3(3) = V1(1)*V2(2)-V1(2)*V2(1)                                   
      END FUNCTION VCROSS                                                   
                                                                            
C   向量的模                                                                
      FUNCTION VNORM(V) RESULT(VM)                                          
          REAL:: V(3), VM                                                   
          VM = SQRT(VDOT(V, V))                                             
      END FUNCTION VNORM                                                    
                                                                            
C   坐标变换                                                                
    ! r = Xi*Ei = xi*ei                                                     
    ! 旧坐标Xi,新坐标xi,均为行向量                                          
    ! 旧坐标基Ei,新坐标基ei                                                 
    ! =>  Xi*Ei·ej = xi*ei·ej = xj                                          
    ! 或  xi*ei·Ej = Xi*Ei·Ej = Xj                                          
    !                  | E1*e1  E1*e2  E1*e3 |                              
    ! 令RMTX = Ei·ej = | E2*e1  E2*e2  E2*e3 |                              
    !                  | E3*e1  E3*e2  E3*e3 |                              
    ! 则x = X*RMTX,或X = x*RMTX逆 = x*RMTX转置                              
    !                                        | E1 |                         
    ! 特别地,若新坐标系为全局坐标架,则RMTX = | E2 |                         
    !                                        | E3 |                         
      FUNCTION PROT(RMTX, POLD) RESULT(PNEW)                                
          REAL:: RMTX(3, 3), POLD(3), PNEW(3)                               
          PNEW(1) = VDOT(POLD, RMTX(:,1))                                   
          PNEW(2) = VDOT(POLD, RMTX(:,2))                                   
          PNEW(3) = VDOT(POLD, RMTX(:,3))                                   
      END FUNCTION PROT                                                     
                                                                            
C   设置局部坐标系单位向量                                                  
    ! P1为局部坐标原点,P1P2为局部1轴,局部3轴为P1P2×P1P3                     
      FUNCTION SETLOC(P1, P2, P3) RESULT(ELOC)                              
          REAL:: P1(3), P2(3), P3(3), ELOC(3, 3)                            
          ELOC(1,:) = (P2-P1)/VNORM(P2-P1)                                  
          ELOC(3,:) = VCROSS(P2-P1,P3-P1)/VNORM(VCROSS(P2-P1,P3-P1))        
          ELOC(2,:) = VCROSS(ELOC(3,:), ELOC(1,:))                          
      END FUNCTION SETLOC                                                   
                                                                            
C   计算任意点PA与当前时刻热源平面的距离                                    
    ! P1为当前时刻焊接点,P2为焊接前进方向点,P3为焊趾点,PA为空间中任意一点   
      FUNCTION DISTANCE (P1, P2, P3, PA) RESULT(DIS)                        
          REAL:: P1(3), P2(3), P3(3), PA(3)                                 
          REAL:: DIS(3)                                                     
          REAL:: VX(3), VY(3), VZ(3)                                        
          REAL:: VTMP(3, 3)                  ! 临时变量                     
          VTMP = SETLOC(P1, P2, P3)          ! 以P1为原点,搭建局部直角坐标系
          VX = VTMP(1,:)                     ! 焊接前进方向单位向量         
          VY = VTMP(2,:)                     ! 焊缝宽度方向单位向量         
          VZ = VTMP(3,:)                     ! 焊缝深度方向单位向量         
    ! 计算点PA到热源平面的距离(带符号)：d = |P1PA|*COS<V,P1PA> = V·P1PA/|V| 
          DIS(1) = VDOT(VX, PA-P1)                                          
          DIS(2) = VDOT(VY, PA-P1)                                          
          DIS(3) = VDOT(VZ, PA-P1)                                          
      END FUNCTION DISTANCE                                                 
                                                                            
      END MODULE WELDMATH                                                   
