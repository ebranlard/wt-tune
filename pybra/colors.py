def rgb2hex(C,g=None,b=None):
    if len(C)==3 :
        r=C[0]
        g=C[1]
        b=C[2]
    if r<1.1 and g<1.1 and b<1.1:
        r=r*255
        g=g*255
        b=b*255
        
    return '#%02X%02X%02X' % (r,g,b)


def fColrs_hex(*args):
    return rgb2hex(fColrs(*args))

def fColrs(i=-1,n=-1,bBW=True):
    # Possible calls
# M=fColrs()  : returns a nx3 matrix of RBG colors
# C=fColrs(i) : cycle through colors, modulo the number of color
# G=fColrs(i,n) : return a grey color (out of n), where i=1 is black
# % Thrid argument add a switch possibility between black and white or colors:
# % G=fColrs(i,n,1) : return a grey color (out of n), where i=1 is black
# % G=fColrs(i,n,0) : cycle through colors
    import numpy as np
# 
    MathematicaBlue       = np.array([63,63,153])/255.     ; 
    MathematicaRed        = np.array([153,61,113])/255.    ; 
    MathematicaGreen      = np.array([61,153,86])/255.     ; 
    MathematicaYellow     = np.array([152,140,61])/255.    ; 
    MathematicaLightBlue  = np.array([159,159,204 ])/255.  ; 
    MathematicaLightRed   = np.array([204,158,184 ])/255.  ; 
    MathematicaLightGreen = np.array([158,204,170  ])/255. ; 
    # 
    ManuDarkBlue    = np.array([0,0,0.7 ])       ; 
    ManuDarkRed     = np.array([138,42,93])/255.  ; 
    ManuDarkOrange  = np.array([245,131,1])/255.  ; 
    ManuDarkOrange  = np.array([198,106,1])/255.  ; 
    ManuLightOrange = np.array([255.,212,96])/255. ; 
    # 
    Red    = np.array([1,0,0])     ; 
    Blue   = np.array([0,0,1])     ; 
    Green  = np.array([0,0.6,0])   ; 
    Yellow = np.array([0.8,0.8,0]); 

    MatlabGreen   = np.array([0,0.5,1                             ] );
    MatlabCyan    = np.array([0.0e+0  , 750.0e-03 ,  750.0e-03    ] );
    MatlabMagenta = np.array([ 750.0e-03 ,  0.0e+0 ,  750.0e-03   ] );
    MatlabYellow  = np.array([750.0e-03 ,  750.0e-03 ,  0.0e+0    ] );
    MatlabGrey    = np.array([250.0e-03  , 250.0e-03 ,  250.0e-03 ] );


    # Table of Color used
    mcolrs=np.array([
            MathematicaBlue,
            MathematicaGreen,
            ManuDarkRed,
            ManuDarkOrange,
            MathematicaLightBlue,
            MathematicaLightGreen,
            MathematicaLightRed,
            ManuLightOrange,
            Blue,
            Green,
            Red,
            Yellow,
            MatlabCyan,
            MatlabMagenta ]);
    # 
    if i==-1:
        return mcolrs
    elif (i!=-1 and n==-1):
        return mcolrs[np.mod(i-1,len(mcolrs)),:];
    elif (i!=-1 and n!=-1):
        if bBW:
            if n==1:
                return [0,0,0]
            else:
                return [0.55,0.55,0.55]*(v-1)/(n-1); #grayscale
        else:
            return mcolrs[mod(i-1,len(mcolrs,1)),:];
    else:
        return mcolrs
        
# 
def test_colrs():

    #from ebra.plot import *
    from numpy import linspace,sin,pi
    from matplotlib import pyplot as plt

    x=linspace(0,2*pi,100);
    plt.figure()
    plt.title('fig_python')
    plt.grid()
    for i in range(30):
        plt.plot(x,sin(x)+i,'-',color=fColrs(i))
    plt.xlabel('x coordinate [m]')
    plt.ylabel('Velocity  U_i [m/s]')
    plt.xlim([0,2*pi])

    plt.show()
if __name__ == "__main__":
    test_colrs()
