from screeninfo import get_monitors
import matplotlib.pyplot as plt

def getRightScreenArea(LeftPanel=0,BottomPanel=0,TopPanel=0,RightPanel=0):
    Monitors=get_monitors()
    if len(Monitors)>1:
        for m in get_monitors():
            if m.__repr__().find('+0+0')<=0:
                SA = m
    else:
        SA=Monitors[0]
    #try:
    #except:
    #    SA         = get_monitors()[0]
    SA.y      += LeftPanel
    SA.x      += TopPanel
    SA.width  = SA.width-LeftPanel-RightPanel
    SA.height = SA.height-BottomPanel-TopPanel
    return SA

def getLeftScreenArea(LeftPanel=105,BottomPanel=0,TopPanel=0,RightPanel=0):
    SA         = get_monitors()[0]
    SA.y      += LeftPanel
    SA.x      += TopPanel
    SA.width  = SA.width-LeftPanel-RightPanel
    SA.height = SA.height-BottomPanel-TopPanel
    return SA

def fig_grid(nX=1,nY=1,AreaName=None,ScreenName='rightscreen',Area=None):
    if Area is None:
        # Determining Area from a acombination of area name and screen name
        if ScreenName.lower()=='rightscreen':
            ScreenArea=getRightScreenArea()
        elif ScreenName.lower()=='leftscreen':
            ScreenArea=getLeftScreenArea()
        else :
            raise Exception('Unknown ScreenName `{}`'.format(ScreenName))

        #print('Screen Position  : X={} Y={}'.format(ScreenArea.x,ScreenArea.y))
        #print('Screen Dimensions: Width={} Height={}'.format(ScreenArea.width,ScreenArea.height))
        if AreaName is not None:
            if AreaName.lower()=='leftscreen':
                Area=getLeftScreenArea()
            elif AreaName.lower()=='rightscreen':
                Area=getRightScreenArea()
            elif AreaName.lower()=='top':
                Area=ScreenArea
                Area.height /=2
            elif AreaName.lower()=='bottom':
                Area=ScreenArea
                Area.height/=2
                Area.y +=Area.height
            elif AreaName.lower()=='left':
                Area=ScreenArea
                Area.width/=2
            elif AreaName.lower()=='right':
                Area=ScreenArea
                Area.width/=2
                Area.x +=Area.width
            elif AreaName.lower()=='topright':
                Area=ScreenArea
                Area.width/=2
                Area.height /=2
                Area.x +=Area.width
            elif AreaName.lower()=='topleft':
                Area=ScreenArea
                Area.width/=2
                Area.height /=2
            elif AreaName.lower()=='bottomright':
                Area=ScreenArea
                Area.width/=2
                Area.height /=2
                Area.x +=Area.width
                Area.y +=Area.height
            elif AreaName.lower()=='bottomleft':
                Area=ScreenArea
                Area.width/=2
                Area.height /=2
                Area.y +=Area.height
            else:
                raise Exception('Screen area name not supported `{}`'.format(Name))
        else: 
            Area = ScreenArea

    W=int(Area.width)
    H=int(Area.height)
    #print('Area  Position  : X={} Y={}'.format(Area.x,Area.y))
    #print('Area  Dimensions: Width={} Height={}'.format(W,H))

    # Creating figures
    figs=[]
    fig1=plt.figure()
    #bKeepFirstFig=fig1.number!=1
    bKeepFirstFig=False
    for i in range(nX):
        for j in range(nY):
            if i==0 and j==0 and bKeepFirstFig:
                figs.append(fig1)
            else:
                figs.append(plt.figure())
            
    fW=int(W/nY)
    fH=int(H/nX)
    #print('Window Dimensions: Width={} Height={}'.format(fW,fH))
    c=0
    for i in range(nX):
        for j in range(nY):
            fx=int(i*fH+Area.y)
            fy=int(j*fW+Area.x)
            s='{:d}x{:d}+{:d}+{:d}'.format(fW,fH,fy,fx)
            figs[c].canvas.manager.window.wm_geometry(s) # NOTE for TK, QT TODO
            c += 1
    if not bKeepFirstFig:
        plt.close(fig1)
    return figs



if __name__ == "__main__":
    from numpy.linalg import inv
    import numpy as np
    
    plt.ion()
    figs=fig_grid(2,2)
    # figs=[]
    # figs+=fig_grid(AreaName='topright')
    # figs+=fig_grid(AreaName='topleft')
    # figs+=fig_grid(AreaName='bottomleft')
    # figs+=fig_grid(AreaName='bottomright')
    #figs=create_fig_grid(2,2,Name='rightscreen')
    plt.show()

    for f in figs:
        f.add_subplot(111)

    n=1000
    for i in range(1):
        # Plotting 
        x=np.random.rand(1,10)
        y=np.random.rand(1,10)
        for f in figs:
            f.axes[0].clear()
            f.axes[0].scatter(x,y,s=10)
        plt.pause(0.001)
        # Computing something heavy
        a = np.random.rand(n,n); ainv = inv(a)
        # Plotting additional data
        x=np.random.rand(1,10)
        y=np.random.rand(1,10)
        for f in figs:
            f.axes[0].scatter(x,y,color='k',marker='+',s=80)
        plt.pause(0.001)
        # Computing something heavy
        a = np.random.rand(n,n); ainv = inv(a)

    plt.pause(3)    
