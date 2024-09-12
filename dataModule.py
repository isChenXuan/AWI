# --*coding=UTF-8*--
import numpy as np
import numpy.linalg as nla
#from scipy.interpolate import CubicSpline
#from scipy.integrate import quad


def SimpsonIntegrate(f, a, b, n=100):
    dx = (b-a)/n
    x = np.linspace(a, b, n)
    y = np.array([f(i) for i in x])
    xm = (x[:-1] + x [1:])/2.0
    ym = np.array([f(i) for i in xm])
    s = np.sum(y*dx)/3 - (f(a)+f(b))*dx/6 + 2/3*np.sum(ym*dx) 
    return s  


class RombergIntegrate():
    def __init__(self, f, a, b, eps=10**-7):
        self.f = f
        self.a = a
        self.b = b
        self.eps = eps
    
    
    def Trap(self, Iold, k):
        if k == 1:
            Inew = (self.f(self.a)+self.f(self.b))*(self.b-self.a)/2.0
        else:
            n = np.power(2,k-2)
            h = (self.b-self.a)/n 
            x = self.a + (h/2) 
            sum = 0.0
            for i in range(n):
                sum = sum + self.f(x) 
                x = x + h 
            Inew = (Iold + h*sum)/2.0 
        return Inew


    def Richardson(self, R, k):
        for i in range(k-1,0,-1):
            #c = np.power(2,2*(k-i))  # c can be a big number, change c to ic.
            ic = np.power(2,2*(i-k))
            R[i] = (R[i+1]-ic*R[i])/(1-ic) 
        return R


    def Romberg(self):
        T = {} 
        k = 1
        T[1] = self.Trap(0.0, 1)
        former_R = T[1]
        while True:
            k += 1
            T[k] = self.Trap(T[k-1], k)
            T = self.Richardson(T,k)
            if abs(T[1] - former_R) < self.eps:
                return T[1]
            former_R = T[1] 
            

class NAKCubicSpline():
    """
    # Same as 'not-a-knot' CubicSpline of scipy. 
    """
    def __init__(self, x, y):
        self.x = np.array(x)  # discrete x
        self.y = np.array(y)  # discrete y
        self.n = len(x)
        self.h = 1.0*self.x[1:] - self.x[:-1]  # n-1 dimension, dx=x[i+1]-x[i]
        self.dy = 1.0*self.y[1:] - self.y[:-1] # n-1 dimension, dy=y[i+1]-y[i]
        
        # Generating equation coefficient matrix.
        if self.n > 3:
            v1 = self.h.copy()     # v1: left lower part of H.     
            v1[-1] = v1[-1] + self.h[-2]

            v2 = 2.0*np.ones(self.n) # v2: diagonal of H. 
            v2[1:-1] = 2*(self.h[:-1]+self.h[1:])
            v2[0] =- self.h[0]
            v2[-1] =-self.h[-2]

            v3 = self.h.copy()     # v3: right upper part of H. 
            v3[0] = v3[0] + self.h[1]

            H = np.diag(v1,-1) + np.diag(v2,0) + np.diag(v3,1)
            H[0,2] = -self.h[0]
            H[-1,-3] = -self.h[-1]

            d = np.zeros(self.n)
            d[1:-1] = 6.0*(self.dy[1:]/self.h[1:] - self.dy[:-1]/self.h[:-1])
        elif self.n==3:
            H = np.array([           0,         0,   1,            0,        0,   0,
                          self.h[0]**2, self.h[0],   1,            0,        0,   0,
                                     0,         0,   0,            0,        0,   1,
                                     0,         0,   0, self.h[1]**2,self.h[1],   1,
                           2*self.h[0],         1,   0,            0,       -1,   0,
                                     1,         0,   0,           -1,        0,   0]).reshape((6,6))
            d = np.array([self.y[0], self.y[1], self.y[1], self.y[2], 0, 0])
        elif self.n==2:
            H = np.array([        0,  1,
                          self.h[0],  1]).reshape((2,2))
            d = np.array([self.y[0], self.y[1]])
        else:
            return
        
        # Fetch coefficient ai，bi，ci，di.
        try:
            m = np.linalg.solve(H,d)
        except:
            return        
            
        if self.n==2:
            ai = [0.0 ]
            bi = [0.0 ]
            ci = [m[0]]
            di = [m[1]] 
        elif self.n==3:
            ai = [0.0 , 0.0 ] 
            bi = [m[0], m[3]]
            ci = [m[1], m[4]]
            di = [m[2], m[5]] 
        else:    
            ai = 1.0/6/self.h*(m[1:] - m[:-1])
            bi = m[:-1]/2.0
            ci = self.dy/self.h - 1/2.0*self.h*m[:-1] - 1/6.0*self.h*(m[1:] - m[:-1])
            di = self.y[:-1]*1.0
            
        coef = np.array([ai,bi,ci,di])   # In accordance with scipy CubicSpline.
        self.c = coef
        

class DataProcess():
    """
    # Handling with discreted point coordinates data, like
    #      np.array([(x1, y1, z1),
    #                (x2, y2, z2),
    #                (   ...    )]).
    """
    def __init__(self, data, method='user'):
        self.data = data.T
        self.x = self.data[0]
        self.y = self.data[1]
        self.z = self.data[2]
        self.n = len(self.x)
        self.t = np.arange(self.n)
        self.method = method
        
        if self.method == 'user':
            # Using user defined NAK CubicSpline.
            self.csx = NAKCubicSpline(self.t, self.x)
            self.csy = NAKCubicSpline(self.t, self.y)
            self.csz = NAKCubicSpline(self.t, self.z) 
        elif self.method == 'scipy':
            # Using scipy NAK CubicSpline.
            self.csx = CubicSpline(self.t, self.x, bc_type='not-a-knot')
            self.csy = CubicSpline(self.t, self.y, bc_type='not-a-knot')
            self.csz = CubicSpline(self.t, self.z, bc_type='not-a-knot')                    
        
        self.arcLength = self.calArcLength(method='Approximate')
        self.iniVector = self.fetchTangent(r=0.0)
        
    
    def polyFunc(self, c, x, flag):
        """
        # Calculate fuction value.
        # For CubicSpline cs=CubicSpline(t, f), the i-th curve segement was expressed as poly-function of t:
        #     f[i](t) =   a[i]*(t-t[i])**3 +   b[i]*(t-t[i])**2 + c[i]*(t-t[i]) + d[i], i=0,1,...,len(t)-2.
        #    df[i](t) = 3*a[i]*(t-t[i])**2 + 2*b[i]*(t-t[i])    + c[i]
        #   d2f[i](t) = 6*a[i]*(t-t[i])    + 2*b[i]
        # where:
        #   a[i]=cs.c[0,i]
        #   b[i]=cs.c[1,i]
        #   c[i]=cs.c[2,i]
        #   d[i]=cs.c[3,i]
        """
        m = len(c[0])  # Segements of function curve.
        t = np.arange(m+1)
        
        try:
            if x < t[0] or x > t[-1]:
                return
        except:
            if x.all() < t[0] or x.all() > t[-1]:
                return 
        
        f   = [lambda p, i=i:   c[0,i]*(p-t[i])**3 +   c[1,i]*(p-t[i])**2 + c[2,i]*(p-t[i]) + c[3,i] for i in range(m)]
        df  = [lambda p, i=i: 3*c[0,i]*(p-t[i])**2 + 2*c[1,i]*(p-t[i])    + c[2,i]                   for i in range(m)]
        d2f = [lambda p, i=i: 6*c[0,i]*(p-t[i])    + 2*c[1,i]                                        for i in range(m)]
        
        try:
            indx = np.where(x >= t)[0][-1]
        except:
            indx = np.where(x.all() >= t.all())[0][-1]
            
        if indx < m:
            y   =   f[indx](x)
            dy  =  df[indx](x)
            d2y = d2f[indx](x)
        else:
            y   =   f[-1](x)
            dy  =  df[-1](x)
            d2y = d2f[-1](x)
             
        if flag == 0:
            return y
        elif flag == 1:
            return dy
        elif flag == 2:
            return d2y  

    
    def evaluate(self, p, flag): 
        """
        # Fetch function value.
        """ 
        if flag == 0:   # Function value. 
            return np.array([self.polyFunc(self.csx.c, p, flag=0),
                             self.polyFunc(self.csy.c, p, flag=0),
                             self.polyFunc(self.csz.c, p, flag=0)])
        elif flag == 1: # First-order derivative. 
            return np.array([self.polyFunc(self.csx.c, p, flag=1),
                             self.polyFunc(self.csy.c, p, flag=1),
                             self.polyFunc(self.csz.c, p, flag=1)])   
        elif flag == 2: # Second-order derivative.
            return np.array([self.polyFunc(self.csx.c, p, flag=2),
                             self.polyFunc(self.csy.c, p, flag=2),
                             self.polyFunc(self.csz.c, p, flag=2)]) 
                

    def calArcLength(self, method='Approximate'):
        """
        # Calculate arc length.
        """    
        ds = lambda x : nla.norm(self.evaluate(x,flag=1))
        
        arcLength = [0.0]
        for i in range(self.n-1):
            if method=='Quad':
                # Using scipy quad
                length = quad(ds, self.t[i], self.t[i+1])
                arcLength.append(arcLength[i] + length[0]) 
            else:
                if method=='Romberg':    # Best.
                    # Using RombergIntegrate
                    r = RombergIntegrate(ds, self.t[i], self.t[i+1])
                    length = r.Romberg()
                elif method=='Simpson':   # Worst.
                    # Using SimpsonIntegrate 
                    length = SimpsonIntegrate(ds, self.t[i], self.t[i+1])
                elif method=='Approximate':
                    # Approximate value.  # Middle.
                    length = nla.norm(self.data.T[i+1,:] - self.data.T[i,:])                                          
                arcLength.append(arcLength[i] + length) 	
    
        return np.array(arcLength)
    

    def fetchTangent(self, r):
        # r : Relative positon of current point at weld line, should be within [0,1].
        #  r=0: start point; r=1: end point.
        p = r*self.arcLength[-1]
        df = self.evaluate(p, flag=1)
        return df/nla.norm(df)


    def symbolCalculate(self, t, time, flag=0):
        # time : Time list from start point to current position.
        time = np.array(time)
        csx = NAKCubicSpline(time, self.x)
        csy = NAKCubicSpline(time, self.y)
        csz = NAKCubicSpline(time, self.z)
        
        time = np.around(time, decimals=4)
        cx = np.around(csx.c, decimals=4)
        cy = np.around(csy.c, decimals=4)
        cz = np.around(csz.c, decimals=4)
        
        fx = ["{}*({}-{})**3+{}*({}-{})**2+{}*({}-{})+{}".format(cx[0,i], t, time[i], cx[1,i], t, time[i], cx[2,i], t, time[i], cx[3,i]) for i in range(self.n-1)]
        fy = ["{}*({}-{})**3+{}*({}-{})**2+{}*({}-{})+{}".format(cy[0,i], t, time[i], cy[1,i], t, time[i], cy[2,i], t, time[i], cy[3,i]) for i in range(self.n-1)]
        fz = ["{}*({}-{})**3+{}*({}-{})**2+{}*({}-{})+{}".format(cz[0,i], t, time[i], cz[1,i], t, time[i], cz[2,i], t, time[i], cz[3,i]) for i in range(self.n-1)]
        
        dfx = ["{}*3*({}-{})**2+{}*2*({}-{})+{}".format(cx[0,i], t, time[i], cx[1,i], t, time[i], cx[2,i]) for i in range(self.n-1)]
        dfy = ["{}*3*({}-{})**2+{}*2*({}-{})+{}".format(cy[0,i], t, time[i], cy[1,i], t, time[i], cy[2,i]) for i in range(self.n-1)]
        dfz = ["{}*3*({}-{})**2+{}*2*({}-{})+{}".format(cz[0,i], t, time[i], cz[1,i], t, time[i], cz[2,i]) for i in range(self.n-1)]    
        
        if flag==0:
            return fx,fy,fz 
        elif flag==1:
            return dfx,dfy,dfz


class Circle():
    """
    # Create circle from three points.        
    """
    def __init__(self, P1, P2, P3):
        self.P1 = np.array(P1)
        self.P2 = np.array(P2)
        self.P3 = np.array(P3)

        self.getCircle()
        self.getLocalAxis()
        
        
    def getCircle(self):
        V12 = self.P2 - self.P1
        V13 = self.P3 - self.P1
        VN = np.cross(V12, V13)    # Normal vector of circle plane.
        
        M1 = (self.P1 + self.P2)/2.0
        M2 = (self.P1 + self.P3)/2.0 
        
        b1 = np.dot(V12, M1)
        b2 = np.dot(V13, M2)
        b3 = np.dot(VN, self.P1)
        
        A = np.array([V12, V13, VN])
        b = np.array([b1, b2, b3])
        
        PC = np.linalg.solve(A, b)  
        RA = nla.norm(self.P1 - PC) 

        self.PC = PC         # Center of circle.
        self.RA = RA         # Radius of circle.


    def getLocalAxis(self):
        vNorm = lambda x : np.array([0.0,0.0,0.0]) if(nla.norm(x)==0) else x/nla.norm(x)
        self.vz = vNorm(np.cross(self.P2 - self.P1, self.P3 - self.P1)) # Normal vector of circle plane, also local z-axis.
        self.vx = vNorm(self.P1 - self.PC)                              # Local x-axis.
        self.vy = np.cross(self.vz, self.vx)                            # Local y-axis.   
        self.local2Global = np.round(np.array([self.vx, self.vy, self.vz]).T, decimals=4)  # Transformation matrix from loacal coordinates to global coordinates.      

        
# Test	
if __name__ == '__main__':					   
    # read node list csv.
    df = pd.read_csv('tee-weldLineNodes-1.csv', encoding='utf-8')
    darray = np.array(df)
    data = np.array(darray[:,2:5])
    
    dp = DataProcess(data, method='user')
    print('quad:',dp.calArcLength(method='Quad')[-1])
    print('romberg:',dp.calArcLength(method='Romberg')[-1])
    #print('simpson:',dp.calArcLength(method='Simpson')[-1])
    print('approx:',dp.calArcLength(method='Approximate')[-1])
    
    df = pd.DataFrame(dp.csx.c.T)
    df.to_csv('coeffx.csv')
    
    time = np.arange(dp.n)
    f = dp.symbolCalculate('r', time, flag=0)
    df = dp.symbolCalculate('r', time, flag=1)
    
    print(f[0])
    print(df[0])
    
    P1 = (2, 2, 1)
    P2 = (1.414, 3.414, 1)
    P3 = (0, 0, 1)
    c = Circle(P1, P2, P3)    