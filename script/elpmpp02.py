#!/usr/bin/env python

import os
import math
import struct
import argparse
import urllib.request
import numpy as np

PATH_SELF = os.path.realpath(os.path.dirname(__file__))

np_dtype = np.float64
rad = 648000.0 / math.pi
deg = math.pi / 180.0


def DMS(ideg:float, imin:float, sec:float)->float:
    deg = math.pi / 180.0
    return (ideg + imin/60.0 + sec/3600.0) * deg


class ELPMPP02:
    def __init__(self, icor: int, path_base: str='.'):
        """
        icor     : Index of the corrections to the constants (integer).

        The nominal values of some constants have to be corrected.
        There are 2 sets of corrections. The corresponding choice is 
        marked by the parameter 'icor' (used in INITIAL).
        icor=0, the constants are fitted to LLR observations provided 
                from 1970 to 2001; it is the default value;
        icor=1, the constants are fitted to DE405 ephemeris over one 
                century (1950-1960); the lunar angles W1, W2, W3 receive 
                also additive corrections to the secular coefficients.
 
        --- NOTATIONS --------------------------------------------------------
       
        Moon elements (polynomials coefficients until order 4):
        w[0,0:4]  : mean longitude of the Moon                        (W1)
        w[1,0:4]  : mean longitude of the lunar perigee               (W2)
        w[2,0:4]  : mean longitude of the lunar ascending node        (W3)
        zeta(0:4) : mean longitude of the Moon + precession        (W1+pt)
                    p is the precession rate and t is the time
       
        Earth-Moon (EMB) elements (polynomials coefficients until order 4):
        eart(0:4) : mean longitude of EMB                             (Te)
        peri(0:4) : mean longitude of the perihelion of EMB          (Pip)
       
        Delaunay arguments (polynomials coefficients until order 4):
        delu[0,0:4]: D  =  W1 - Te + 180 degrees                       (D)
        delu[1,0:4]: F  =  W1 - W3                                     (F)
        delu[2,0:4]: l  =  W1 - W2   mean anomaly of the Moon          (l)
        delu[3,0:4]: l' =  Te - Pip  mean anomaly of EMB              (l')
       
        Planetary arguments (mean longitudes at J2000 and mean motions):
        p(0,0:1)  : mean longitude of Mercury
        p(1,0:1)  : mean longitude of Venus
        p(2,0:1)  : mean longitude of EMB (eart(0:1))
        p(3,0:1)  : mean longitude of Mars 
        p(4,0:1)  : mean longitude of Jupiter
        p(5,0:1)  : mean longitude of Saturn
        p(6,0:1)  : mean longitude of Uranus
        p(7,0:1)  : mean longitude of Neptune
       
        Moon constants:
        nu        : mean motion of the Moon (W1(1,1))                 (Nu)
        g         : half coefficient of sin(F) in latitude         (Gamma)
        e         : half coefficient of sin(l) in longitude            (E)
        np        : mean motion of EMB (eart(1))                      (n')
        ep        : eccentricity of EMB                               (e')
        alpha     : Ratio of the semi-major axis (Moon / EMB)
        am        : Ratio of the mean motions (EMB / Moon)
        dtasm     : (2*alpha) / (3*am)     
       
        Corrections to the constants Nu, Gamma, E, n', e':
        delnu     : to the mean motion of the Moon 
        delg      : to the half coefficient of sin(F) in latitude
        dele      : to the half coefficient of sin(l) in longitude
        delnp     : to the mean motion of EMB 
        delep     : to the eccentricity of EMB
       
        Precession of the longitude of the ascending node of the mean 
        ecliptic of date on fixed ecliptic J2000:
        pi[i=1,5] : sine coefficients
        qi[i=1,5] : cosine coefficients
       
        """        

        self.sc = 36525.0
        self.a405=384747.9613701725
        self.aelp=384747.980674318

        self.filename = [
            'ELP_MAIN.S1',
            'ELP_MAIN.S2',
            'ELP_MAIN.S3',
            'ELP_PERT.S1',
            'ELP_PERT.S2',
            'ELP_PERT.S3'
        ]

        self.icor = icor

        self.bp = np.array([
            +0.311079095e+0, -0.448239800e-2, -0.110248500e-2, +0.105606200e-2, +0.509280000e-4, 
            -0.103837907e+0, +0.668287000e-3, -0.129807200e-2, -0.178028000e-3, -0.373420000e-4,
        ]).reshape(2, 5).T

        
# --- Value of the correction to the constant of precession ------------

        self.Dprec=-0.29965e0                 # ***** IAU 2000A

# --- Constants for the evaluation of the partial derivatives ----------

        self.am     =  0.074801329e0
        self.alpha  =  0.002571881e0
        self.dtasm  =  (2.0 * self.alpha)/(3.0 * self.am)
        self.xa     =  (2.0 * self.alpha)/3.0

# --- Corrections to constants -----------------------------------------

#     Values of the corrections to the constants fitted to LLR.
#     Fit 13-05-02 (2 iterations) except Phi and eps w2_1 et w3_1
        if self.icor == 0:
            self.Dw1_0   = -0.10525e0
            self.Dw2_0   =  0.16826e0
            self.Dw3_0   = -0.10760e0
            self.Deart_0 = -0.04012e0
            self.Dperi   = -0.04854e0
            self.Dw1_1   = -0.32311e0
            self.Dgam    =  0.00069e0
            self.De      = +0.00005e0
            self.Deart_1 =  0.01442e0
            self.Dep     =  0.00226e0
            self.Dw2_1   =  0.08017e0
            self.Dw3_1   = -0.04317e0
            self.Dw1_2   = -0.03794e0
        else:
#     Values of the corrections to the constants fitted to DE405
#     over the time interval (1950-2060)
            self.Dw1_0   = -0.07008e0
            self.Dw2_0   =  0.20794e0
            self.Dw3_0   = -0.07215e0
            self.Deart_0 = -0.00033e0
            self.Dperi   = -0.00749e0
            self.Dw1_1   = -0.35106e0
            self.Dgam    =  0.00085e0
            self.De      = -0.00006e0
            self.Deart_1 =  0.00732e0
            self.Dep     =  0.00224e0
            self.Dw2_1   =  0.08017e0
            self.Dw3_1   = -0.04317e0
            self.Dw1_2   = -0.03743e0
            
# --- Fundamental arguments (Moon and EMB) -----------------------------

        self.w = np.zeros([3, 5], dtype=np_dtype)
        self.eart = np.zeros([5], dtype=np_dtype)
        self.peri = np.zeros([5], dtype=np_dtype)
        
        self.w[0,0] = DMS(218,18,59.95571+self.Dw1_0)         #***** ELP
        self.w[0,1] = (1732559343.73604+self.Dw1_1)/rad       #***** ELP
        self.w[0,2] = (        -6.8084 +self.Dw1_2)/rad       #***** DE405
        self.w[0,3] =          0.66040e-2/rad                 #***** ELP
        self.w[0,4] =         -0.31690e-4/rad                 #***** ELP

        self.w[1,0] = DMS( 83,21,11.67475+self.Dw2_0)         #***** ELP
        self.w[1,1] = (  14643420.3171 +self.Dw2_1)/rad       #***** DE405
        self.w[1,2] = (       -38.2631)/rad                   #***** DE405
        self.w[1,3] =         -0.45047e-1/rad                 #***** ELP
        self.w[1,4] =          0.21301e-3/rad                 #***** ELP

        self.w[2,0] = DMS(125, 2,40.39816+self.Dw3_0)         #***** ELP
        self.w[2,1] = (  -6967919.5383 +self.Dw3_1)/rad       #***** DE405
        self.w[2,2] = (         6.3590)/rad                   #***** DE405
        self.w[2,3] =          0.76250e-2/rad                 #***** ELP
        self.w[2,4] =         -0.35860e-4/rad                 #***** ELP

        self.eart[0]= DMS(100,27,59.13885+self.Deart_0)       #***** VSOP2000
        self.eart[1]= (129597742.29300 +self.Deart_1)/rad     #***** VSOP2000
        self.eart[2]=         -0.020200/rad                   #***** ELP
        self.eart[3]=          0.90000e-5/rad                 #***** ELP
        self.eart[4]=          0.15000e-6/rad                 #***** ELP

        self.peri[0]= DMS(102,56,14.45766+self.Dperi)         #***** VSOP2000
        self.peri[1]=       1161.24342/rad                    #***** VSOP2000
        self.peri[2]=          0.529265/rad                   #***** VSOP2000
        self.peri[3]=         -0.11814e-3/rad                 #***** VSOP2000
        self.peri[4]=          0.11379e-4/rad                 #***** VSOP2000

        if icor == 1:
            self.w[0,3] += -0.00018865e0/rad
            self.w[0,4] += -0.00001024e0/rad
            self.w[1,2] += +0.00470602e0/rad
            self.w[1,3] += -0.00025213e0/rad
            self.w[2,2] += -0.00261070e0/rad
            self.w[2,3] += -0.00010712e0/rad
#
# --- Corrections to the mean motions of the Moon angles W2 and W3 -----
#     infered from the modifications of the constants

        x2     =   self.w[1,1]/self.w[0,1]
        x3     =   self.w[2,1]/self.w[0,1]
        y2     =   self.am*self.bp[0,0]+self.xa*self.bp[4,0]
        y3     =   self.am*self.bp[0,1]+self.xa*self.bp[4,1]

        d21    =   x2-y2
        d22    =   self.w[0,1]*self.bp[1,0]
        d23    =   self.w[0,1]*self.bp[2,0]
        d24    =   self.w[0,1]*self.bp[3,0]
        d25    =   y2/self.am
       
        d31    =   x3-y3
        d32    =   self.w[0,1]*self.bp[1,1]
        d33    =   self.w[0,1]*self.bp[2,1]
        d34    =   self.w[0,1]*self.bp[3,1]
        d35    =   y3/self.am

        Cw2_1  =  d21*self.Dw1_1+d25*self.Deart_1+d22*self.Dgam+d23*self.De+d24*self.Dep
        Cw3_1  =  d31*self.Dw1_1+d35*self.Deart_1+d32*self.Dgam+d33*self.De+d34*self.Dep

        self.w[1,1] +=  Cw2_1/rad
        self.w[2,1] +=  Cw3_1/rad

#
# --- Arguments of Delaunay --------------------------------------------
#
        self.delu = np.zeros([4, 5], dtype=np_dtype)
        for i in range(5):
            self.delu[0,i] = self.w[0,i]-self.eart[i]                  #***** D 
            self.delu[1,i] = self.w[0,i]-self.w[2,i]                   #***** F 
            self.delu[2,i] = self.w[0,i]-self.w[1,i]                   #***** l 
            self.delu[3,i] = self.eart[i]-self.peri[i]                 #***** l'
        self.delu[0,0] += math.pi

#
# --- Planetary arguments (mean longitudes and mean motions) -----------
#
        self.p = np.zeros([8, 5], dtype=np_dtype)
        self.p[0,0] = DMS(252,15, 3.216919e0)                   #***** VSOP2000
        self.p[1,0] = DMS(181,58,44.758419e0)                   #***** VSOP2000
        self.p[2,0] = DMS(100,27,59.138850e0)                   #***** VSOP2000
        self.p[3,0] = DMS(355,26, 3.642778e0)                   #***** VSOP2000
        self.p[4,0] = DMS( 34,21, 5.379392e0)                   #***** VSOP2000
        self.p[5,0] = DMS( 50, 4,38.902495e0)                   #***** VSOP2000
        self.p[6,0] = DMS(314, 3, 4.354234e0)                   #***** VSOP2000
        self.p[7,0] = DMS(304,20,56.808371e0)                   #***** VSOP2000

        self.p[0,1] = 538101628.66888e0/rad                     #***** VSOP2000
        self.p[1,1] = 210664136.45777e0/rad                     #***** VSOP2000
        self.p[2,1] = 129597742.29300e0/rad                     #***** VSOP2000
        self.p[3,1] =  68905077.65936e0/rad                     #***** VSOP2000
        self.p[4,1] =  10925660.57335e0/rad                     #***** VSOP2000
        self.p[5,1] =   4399609.33632e0/rad                     #***** VSOP2000
        self.p[6,1] =   1542482.57845e0/rad                     #***** VSOP2000
        self.p[7,1] =    786547.89700e0/rad                     #***** VSOP2000

#
# --- Zeta : Mean longitude W1 + Rate of the precession ----------------
#
        self.zeta = np.zeros([5], dtype=np_dtype)
        self.zeta[0]= self.w[0,0]
        self.zeta[1]= self.w[0,1] + (5029.0966e0+self.Dprec)/rad
        self.zeta[2]= self.w[0,2]
        self.zeta[3]= self.w[0,3]
        self.zeta[4]= self.w[0,4]
#
# --- Corrections to the parameters: Nu, E, Gamma, n' et e' ------------
#
        self.delnu  = (+0.55604e0+self.Dw1_1)/rad/self.w[0,1]        #***** ELP
        self.dele   = (+0.01789e0+self.De)/rad                       #***** ELP
        self.delg   = (-0.08066e0+self.Dgam)/rad                     #***** ELP
        self.delnp  = (-0.06424e0+self.Deart_1)/rad/self.w[0,1]      #***** ELP
        self.delep  = (-0.12879e0+self.Dep)/rad                      #***** ELP
#
# --- Precession coefficients for P and Q (Laskar, 1986) ---------------
#
        self.p1 =  0.10180391e-04
        self.p2 =  0.47020439e-06
        self.p3 = -0.5417367e-09
        self.p4 = -0.2507948e-11
        self.p5 =  0.463486e-14

        self.q1 = -0.113469002e-03
        self.q2 =  0.12372674e-06
        self.q3 =  0.1265417e-08
        self.q4 = -0.1371808e-11
        self.q5 = -0.320334e-14

        self.max1 = 2645
        self.max2 = 33256

        pis2 = math.pi / 2.0
        dpi = 2.0 * math.pi

        b    = np.zeros([5], dtype=np_dtype)
        ilu  = np.zeros([4], dtype=np.int32)
        ipla = np.zeros([11], dtype=np_dtype)
        ifi  = np.zeros([16], dtype=np_dtype)

        self.cmpb = np.zeros([self.max1], dtype=np_dtype)
        self.fmpb = np.zeros([5, self.max1], dtype=np_dtype)
        self.nmpb = np.zeros([3, 3], dtype=np.int32)

        self.cper = np.zeros([self.max2], dtype=np_dtype)
        self.fper = np.zeros([5, self.max2], dtype=np_dtype)
        self.nper = np.zeros([3, 4, 3], dtype=np.int32)

#
# --- Reading Main Problem series --------------------------------------
#

        ir = 0
        for iv in range(3):
            f = open(os.path.join(path_base, self.filename[iv]))
            line = f.readline()
            if line == '':
                return
            
            # format 1001
            self.nmpb[iv, 0] = int(line[25:35].strip())
            self.nmpb[iv, 1] = ir
            self.nmpb[iv, 2] = self.nmpb[iv, 0] + self.nmpb[iv, 1] 

            for n in range(self.nmpb[iv, 0].item()):
                line = f.readline()
                # format 1002
                values = []
                values.extend(int(line[i:i+3]) for i in range(0, 12, 3))
                values.append(float(line[14:27]))
                values.extend(float(line[i:i+12]) for i in range(27, 87, 12))

                ilu[:] = values[:4]
                a = values[4]
                b[:] = values[5:]

                tgv = b[0] + self.dtasm * b[4]
                if iv == 2:
                    a -= 2.0 * a * self.delnu/3.0

                self.cmpb[ir] = a \
                    + tgv * (self.delnp - self.am * self.delnu) \
                    + b[1] * self.delg + b[2] * self.dele + b[3] * self.delep
                
                for k in range(5):
                    for i in range(4):
                        self.fmpb[k, ir] += ilu[i] * self.delu[i, k]
                
                if iv == 2:
                    self.fmpb[0, ir] += pis2
                
                ir += 1
            f.close()

#
# --- Reading Perturbations series -------------------------------------
#
        ir = 0
        for iv in range(3):
            f = open(os.path.join(path_base, self.filename[iv + 3]))
            for it in range(4):
                line = f.readline()
                if line == '': return

                self.nper[iv, it, 0], ipt = int(line[25:35]), int(line[35:45])
                self.nper[iv, it, 1] = ir 
                self.nper[iv, it, 2] = self.nper[iv, it, 0] + self.nper[iv, it, 1] 
                if self.nper[iv, it, 0] == 0:
                    continue
                
                for n in range(self.nper[iv, it, 0].item()):
                    line = f.readline().replace('D', 'e')

                    values = []
                    values.append(int(line[:5]))
                    values.extend(float(line[i:i+20]) for i in range(5, 45, 20))
                    values.extend(int(line[i:i+3]) for i in range(45, 93, 3))
                    icount = values[0]
                    s, c = values[1:3]
                    ifi = values[3:]


                    self.cper[ir] = (c*c + s*s) ** 0.5
                    pha = math.atan2(c, s)
                    if pha < 0.0:
                        pha += math.pi * 2.0
                    
                    for k in range(5):
                        self.fper[k, ir] = 0.0
                        if k == 0: self.fper[k, ir] = pha

                        for i in range(0, 4):
                            self.fper[k, ir] += ifi[i] * self.delu[i, k]

                        for i in range(4, 12):
                            self.fper[k, ir] += ifi[i] * self.p[i-4, k]

                        self.fper[k, ir] += ifi[12] * self.zeta[k]
                    ir += 1

    def evaluate(self, tj: float)->np.ndarray:
#
# --- Initialization of time powers ------------------------------------
#
        t = np.power(tj/self.sc, np.arange(0, 5, dtype=np_dtype))

# --- Evaluation of the series: substitution of time in the series -----
#     Variable iv=1 : Longitude
#     Variable iv=2 : Latitude
#     Variable iv=3 : Distance
#
        v = np.zeros([6], dtype=np_dtype)

        for iv in range(3):

#
# ------ Main Problem series -------------------------------------------
#
            for n in range(self.nmpb[iv, 1], self.nmpb[iv, 2]):
                x = self.cmpb[n]
                y = self.fmpb[0, n]
                yp = 0.0
                for k in range(1, 5):
                    y  += self.fmpb[k, n] * t[k]
                    yp += k * self.fmpb[k, n] * t[k-1]
                v[iv]   += x * math.sin(y)
                v[iv+3] += x * yp * math.cos(y)
#
# ------ Perturbations series  -----------------------------------------
#
            for it in range(4):
                for n in range(self.nper[iv, it, 1], self.nper[iv, it, 2]):
                    x = self.cper[n]
                    y = self.fper[0, n]
                    xp = 0.0
                    yp = 0.0
                    if it != 0:
                        xp = it * x * t[it-1]
                    for k in range(1, 5):
                        y  += self.fper[k, n] * t[k]
                        yp += k * self.fper[k, n] * t[k-1]

                    v[iv]   += x * t[it] * math.sin(y)
                    v[iv+3] += xp * math.sin(y) + x * t[it] * yp * math.cos(y)

#
# --- Computation of the rectangular coordinates (Epoch J2000) ---------           
#
        v[0]   = v[0]/rad + self.w[0,0] + self.w[0,1]*t[1] + self.w[0,2]*t[2] + self.w[0,3]*t[3] + self.w[0,4]*t[4]
        v[1]   = v[1]/rad
        v[2]   = v[2]*self.a405/self.aelp
        v[3]   = v[3]/rad + self.w[0,1] + 2.0*self.w[0,2]*t[1] + 3.0*self.w[0,3]*t[2] + 4.0*self.w[0,4]*t[3]
        v[4]   = v[4]/rad

        xyz = np.zeros([6], dtype=np_dtype)
     
        clamb  = math.cos(v[0])
        slamb  = math.sin(v[0])
        cbeta  = math.cos(v[1])
        sbeta  = math.sin(v[1])
        cw     = v[2]*cbeta
        sw     = v[2]*sbeta
     
        x1     = cw*clamb
        x2     = cw*slamb
        x3     = sw
        xp1    = (v[5]*cbeta-v[4]*sw)*clamb-v[3]*x2
        xp2    = (v[5]*cbeta-v[4]*sw)*slamb+v[3]*x1
        xp3    = v[5]*sbeta+v[4]*cw
     
        pw     = (self.p1+self.p2*t[1]+self.p3*t[2]+self.p4*t[3]+self.p5*t[4])*t[1]
        qw     = (self.q1+self.q2*t[1]+self.q3*t[2]+self.q4*t[3]+self.q5*t[4])*t[1]
        ra     = 2.0*((1-pw*pw-qw*qw)**0.5)
        pwqw   = 2.0*pw*qw
        pw2    = 1.0-2.0*pw*pw
        qw2    = 1.0-2.0*qw*qw
        pwra   = pw*ra
        qwra   = qw*ra
     
        xyz[0] = pw2*x1+pwqw*x2+pwra*x3
        xyz[1] = pwqw*x1+qw2*x2-qwra*x3
        xyz[2] = -pwra*x1+qwra*x2+(pw2+qw2-1)*x3
     
        ppw    = self.p1+(2.0*self.p2+3.0*self.p3*t[1]+4.0*self.p4*t[2]+5.0*self.p5*t[3])*t[1]
        qpw    = self.q1+(2.0*self.q2+3.0*self.q3*t[1]+4.0*self.q4*t[2]+5.0*self.q5*t[3])*t[1]
        ppw2   = -4.0*pw*ppw
        qpw2   = -4.0*qw*qpw
        ppwqpw = 2.0*(ppw*qw+pw*qpw)
        rap    = (ppw2+qpw2)/ra
        ppwra  = ppw*ra+pw*rap
        qpwra  = qpw*ra+qw*rap
     
        xyz[3] = (pw2*xp1+pwqw*xp2+pwra*xp3 + ppw2*x1+ppwqpw*x2+ppwra*x3)/self.sc
        xyz[4] = (pwqw*xp1+qw2*xp2-qwra*xp3 + ppwqpw*x1+qpw2*x2-qpwra*x3)/self.sc
        xyz[5] = (-pwra*xp1+qwra*xp2+(pw2+qw2-1.0)*xp3 - ppwra*x1+qpwra*x2+(ppw2+qpw2)*x3)/self.sc

        return xyz

    def export(self):
        # main problem
        coeff_main = {}
        for iv in range(3):
            coeffs = []
            n_start, n_end = self.nmpb[iv, 1], self.nmpb[iv, 2]
            for n in range(n_start, n_end):
                coeffs.append(tuple(
                    [self.cmpb[n]] + self.fmpb[:,n].tolist()
                ))
            coeffs.sort(key=lambda e: abs(e[0]))
            coeff_main[iv] = coeffs
        
        # perturbations
        coeff_perturb = {}

        for iv in range(3):
            for it in range(4):
                n_start, n_end = self.nper[iv, it, 1], self.nper[iv, it, 2]
                coeffs = []
                for n in range(n_start, n_end):
                    coeffs.append(tuple(
                        [self.cper[n]] + self.fper[:, n].tolist()
                    ))
                coeffs.sort(key=lambda e: abs(e[0]))
                coeff_perturb[(iv, it)] = coeffs
        
        return {
            'coeff_main': coeff_main,
            'coeff_perturb': coeff_perturb,
            'a405': self.a405,
            'aelp': self.aelp,
            'w': self.w.reshape(-1).tolist(),
            'p_poly': [ self.p1, self.p2, self.p3, self.p4, self.p5, ],
            'q_poly': [ self.q1, self.q2, self.q3, self.q4, self.q5, ],
            'sc': self.sc,
        }

def main_orig(args: argparse.Namespace):

    f = open('ELPMPP02.PY.TXT', 'w')

#-----------------------------------------------------------------------
#     COMPUTE LUNAR COORDINATES (POSITIONS ANS VELOCITIES):
#     5 DATES FROM JD 2444239.5 TO 2452239.5,
#     TIME INTERVAL = 2000 DAYS,
#     CONSTANTS FITTED TO LLR OBSERVATIONS.
#-----------------------------------------------------------------------

    dj2000=2451545.0

    e=ELPMPP02(0, path_base=args.data)
    tj0=2444239.5
    dtj=2000.e0
    nda=5

    for n in range(1, nda+1):
        t = tj0 + (n-1) * dtj 
        f.write(f'  JD {t:9.1f}\n')

        res = e.evaluate(t - dj2000)
        res = res.tolist()
        f.write(f'  X = {res[0]:17.7f}')
        f.write(f'  Y = {res[1]:17.7f}')
        f.write(f'  Z = {res[2]:17.7f}  km\n')
        f.write(f'  X\'= {res[3]:17.7f}')
        f.write(f'  Y\'= {res[4]:17.7f}')
        f.write(f'  Z\'= {res[5]:17.7f}  km/day\n')

#-----------------------------------------------------------------------
#     COMPUTE LUNAR COORDINATES (POSITIONS ANS VELOCITIES):
#     5 DATES FROM JD 2500000.5 TO 1700000.5,
#     TIME INTERVAL = -200000 DAYS,
#     CONSTANTS FITTED TO JPL EPHEMERIS DE405.
#-----------------------------------------------------------------------
    e=ELPMPP02(1, path_base=args.data)
    tj0=2500000.5
    dtj=-200000.
    nda=5

    for n in range(1, nda+1):
        t = tj0 + (n-1) * dtj 
        f.write(f'  JD {t:9.1f}\n')

        res = e.evaluate(t - dj2000)
        res = res.tolist()
        f.write(f'  X = {res[0]:17.7f}')
        f.write(f'  Y = {res[1]:17.7f}')
        f.write(f'  Z = {res[2]:17.7f}  km\n')
        f.write(f'  X\'= {res[3]:17.7f}')
        f.write(f'  Y\'= {res[4]:17.7f}')
        f.write(f'  Z\'= {res[5]:17.7f}  km/day\n')
    f.close()



def main_export(args: argparse.Namespace):
    path_export = args.output
    if path_export is None:
        raise ValueError('--output is required')
    
    os.makedirs(path_export, exist_ok=True)
    os.makedirs(os.path.join(path_export, 'bins'), exist_ok=True)
    e=ELPMPP02(int(args.correction), path_base=args.data)

    export = e.export()
    size_main = {}

    for iv in range(3):
        data = []
        for e in export['coeff_main'][iv]:
            data.extend(e)

        size_main[iv] = len(data)
        with open(os.path.join(path_export, 'bins', f'coeff_main_{iv}.bin'), 'wb') as f:
            for e in data:
                f.write(struct.pack('d', e))

    size_perturb = {}
    for iv in range(3):
        for it in range(4):
            data = []
            for e in export['coeff_perturb'][(iv, it)]:
                data.extend(e)
            size_perturb[(iv, it)] = len(data) 
            if len(data) == 0: continue
            with open(os.path.join(path_export, 'bins', f'coeff_perturb_{iv}_{it}.bin'), 'wb') as f:
                for e in data:
                    f.write(struct.pack('d', e))

    code = []
    code.append('use include_bytes_plus::include_bytes;')
    code.append('')

    for iv in range(3):
        sz = size_main[iv]
        if sz == 0: continue
        code.append(f'pub static DATA_COEFF_MAIN_{iv}: [u64; {sz}] = include_bytes!("src/bins/coeff_main_{iv}.bin" as u64le);')

    for iv in range(3):
        for it in range(4):
            sz = size_perturb[(iv ,it)]
            if sz == 0: continue
            code.append(f'pub static DATA_COEFF_PERTURB_{iv}_{it}: [u64; {sz}] = include_bytes!("src/bins/coeff_perturb_{iv}_{it}.bin" as u64le);')
    
    for k, v in export.items():
        if k in ['coeff_perturb', 'coeff_main']:
            continue
        if isinstance(v, float):
            code.append(f'pub const CONST_{k.upper()}: f64 = {v};');
        elif isinstance(v, list):
            code.append(f'pub const CONST_{k.upper()}: [f64;{len(v)}] = [');
            v = [f'{e:27.20e}' for e in v]
            for idx_start in range(0, len(v), 3):
                idx_end = min(idx_start + 3, len(v))
                line = ' ' * 2
                line += ', '.join(v[idx_start:idx_end])
                if idx_end < len(v):
                    line += ','
                code.append(line)
            code.append('];')
    
    with open(os.path.join(path_export, 'coeff.rs'), 'w') as f:
        f.write('\n'.join(code)) 

def main_download(args: argparse.Namespace):
    base_url = 'ftp://cyrano-se.obspm.fr/pub/2_lunar_solutions/2_elpmpp02/'
    fnames = [
        'README.TXT',
        'ELPMPP02.for',
        'elpmpp02.pdf',
        'ELP_MAIN.S1',
        'ELP_MAIN.S2',
        'ELP_MAIN.S3',
        'ELP_PERT.S1',
        'ELP_PERT.S2',
        'ELP_PERT.S3',
    ]

    for fname in fnames:
        url = base_url + fname
        path_dst = os.path.join(args.data, fname)
        os.makedirs(os.path.dirname(path_dst), exist_ok=True)

        if os.path.exists(path_dst):
            print(f'! skipping {url}: {path_dst} exists')
        else:
            print(f'* downloading {url} -> {path_dst}')
            urllib.request.urlretrieve(url, path_dst)
    print(f'* done')



def main_check(args: argparse.Namespace):
    import subprocess

    main_download(args)
    subprocess.call([
        'gfortran', 'ELPMPP02.for', '-o', 'ELPMPP02', '-O3'
    ], cwd=args.data)
    subprocess.call([
        './ELPMPP02'
    ], cwd=args.data)

    main_orig(args)
    
    with open('ELPMPP02.PY.TXT') as f:
        res_ported = [x.strip() for x in f.readlines() if x.strip() != '']
    with open(os.path.join(args.data, 'ELPMPP02.TXT')) as f:
        res_orig = [x.strip() for x in f.readlines() if x.strip() != '']
    
    assert len(res_ported) == len(res_orig)
    for line_ported, line_orig in zip(res_ported, res_orig):
        assert line_ported == line_orig

    print('* check passed!')



 
def parse_args()->argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument('mode', choices=['main', 'export', 'download', 'check'])
    parser.add_argument('-d,--data', help='path to the directory containing ELP_{MAIN, PERT}.S{1..3} files', 
                        dest='data', 
                        default=os.path.join(PATH_SELF, 'data'))
    parser.add_argument('-o,--output', help='path to output', 
                        dest='output', 
                        default=os.path.join(PATH_SELF, '..', 'src')
                        )
    parser.add_argument('-c,--correction', choices=['0', '1'], help='correction mode', dest='correction', default='0')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = parse_args()

    mode = args.mode
    if mode == 'main':
        main_orig(args)
    elif mode == 'export':
        main_export(args)
    elif mode == 'download':
        main_download(args)
    elif mode == 'check':
        main_check(args)
    else:
        raise ValueError('unsupported mode {args.mode}')