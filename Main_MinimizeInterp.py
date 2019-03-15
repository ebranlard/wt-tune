import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.optimize import minimize

from mpl_toolkits.mplot3d import Axes3D 
import matplotlib.pyplot as plt

def f(x, y):
    return  (x-2)**2 + (y-5)**2 

x = np.linspace(1, 4, 11)
y = np.linspace(4, 7, 22)
data = f(*np.meshgrid(x, y, indexing='ij', sparse=True))
print(data.shape)
print(type(data))

# data is now a 3D array with data[i,j,k] = f(x[i], y[j], z[k]). Next, define an interpolating function from this data:

fInterp = RegularGridInterpolator((x, y), data)

# Evaluate the interpolating function at the two points (x,y,z) = (2.1, 6.2, 8.3) and (3.3, 5.2, 7.1):
# >>>


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

X, Y = np.meshgrid(x, y)
# zs = np.array([fun(x,y) for x,y in zip(np.ravel(X), np.ravel(Y))])
zs = np.array([f(x,y) for x,y in zip(np.ravel(X), np.ravel(Y))])
print(zs.shape)
print(type(zs.shape))
Z = zs.reshape(X.shape)
print(Z.shape)
print(type(Z.shape))

ax.plot_surface(X, Y, Z)
plt.show()


pts = np.array([[2.1, 6.2], [3.3, 5.2]])
print(fInterp(pts))
bnds = ((1, 4), (4, 7))
print('Minimizing')
res = minimize(fInterp, (2, 5), bounds=bnds)
print(res)


