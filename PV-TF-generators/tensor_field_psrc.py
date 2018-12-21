import vtk
import numpy as np
from vtk.util import numpy_support
from vtk.numpy_interface import dataset_adapter as dsa
import numpy as np
import vtk
import vtk.util.numpy_support as ns

# set user parameters
Dimensions = [11,11,11]
T = 10.0
A = 0.1
eps = 0.1

def double_gyre_vectorized( x, y, t):
	a = eps * np.sin(2.0 * np.pi / T * t)
	b = 1.0 - 2.0 * eps * np.sin(2.0 * np.pi / T * t)
	f = a * x**2 + b * x
	df = 2.0 * a * x + b

	xx = -np.pi * A * np.sin(np.pi * f) * np.cos(np.pi * y)
	yx = np.pi * A * np.cos(np.pi * f) * np.sin(np.pi * y) * df
	zx = np.pi * A * np.cos(np.pi * f)
	xy = -np.pi * A * np.sin(np.pi * f) * np.cos(np.pi * y)
	yy = np.pi * A * np.cos(np.pi * f) * np.sin(np.pi * y) * df
	zy = np.pi * A * np.sin(np.pi * f)
	xz = -np.pi * A * np.sin(np.pi * f) * np.cos(np.pi * y)
	yz = np.pi * A * np.cos(np.pi * f) * np.sin(np.pi * y) * df
	zz = np.pi * A * np.sin(np.pi * f) * np.cos(np.pi * f)

	return np.dstack([xx, yx, zx, xy, yy, zy, xz, yz, zz])

# create meshgrid and spacing
xs = np.linspace(0.0, 2.0, Dimensions[0])
ys = np.linspace(0.0, 1.0, Dimensions[1])
ts = np.linspace(0.0, 10.0, Dimensions[2])
x, y, t = np.meshgrid(xs, ys, ts, indexing='ij')

# set data sink and image attributes
img = self.GetImageDataOutput()
img.SetDimensions(Dimensions)
dims = Dimensions
img.SetSpacing(abs(xs[1] - xs[0]), abs(ys[1] - ys[0]), abs(ts[1] - ts[0]))
img.SetOrigin(xs[0], ys[0], ts[0])
img = dsa.WrapDataObject(output)

v = double_gyre_vectorized(x, y, t)

img.PointData.append(v.reshape((-1, 9), order='F'), 'v')
img.VTKObject.GetPointData().SetActiveScalars('v')

# Request Information Script 1:
# from paraview import util
# Dimensions = [101,101,101]
# op = self.GetOutput()
# util.SetOutputWholeExtent(self, [0, Dimensions[0] - 1, 0, Dimensions[1] - 1, 0, Dimensions[2] - 1])

# //////

# Request Information Script 2
# Dimensions = [101,101,101]
# executive = self.GetExecutive()
# outInfo = executive.GetOutputInformation(0)
# outInfo.Set(executive.WHOLE_EXTENT(), [0, Dimensions[0] - 1, 0, Dimensions[1] - 1, 0, Dimensions[2] - 1], 6)

