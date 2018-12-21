import numpy as np
from vtk.util.vtkAlgorithm import * # use (standalone) vtk module for sake of completeness
from vtk.vtkCommonDataModel import vtkImageData
from vtk.numpy_interface import dataset_adapter as dsa
from timeit import default_timer as timer
from paraview.pythonalgorithm import smproxy, smhint, smdomain, smproperty # import custom pythonalgorithm to generate XMLs via decorators in deprecated versions (<5.6.0)


@smproxy.source(label="Double Gyre")
@smhint.xml('''<ShowInMenu category="My Sources" />''')
class DoubleGyre(VTKPythonAlgorithmBase):
    def __init__(self):
        self.Dimensions = [101, 101, 101]
        self.T = 10.0
        self.A = 0.1
        self.eps = 0.1
        self.Vectorize = True
        VTKPythonAlgorithmBase.__init__(self, nInputPorts=0, nOutputPorts=1, outputType='vtkImageData')

    @smproperty.intvector(name='Dimensions', default_values=[101, 101, 101])
    def SetDimensions(self, dx, dy, dz):
        self.Dimensions = [dx, dy, dz]
        self.Modified()

    @smproperty.doublevector(name='T', default_values=10.0)
    def SetT(self, T):
        self.T = T
        self.Modified()

    @smproperty.doublevector(name='A', default_values=0.1)
    def SetA(self, A):
        self.A = A
        self.Modified()

    @smproperty.doublevector(name='eps', default_values=0.1)
    def Seteps(self, eps):
        self.eps = eps
        self.Modified()

    @smproperty.intvector(name='Vectorize', default_values=1)
    @smdomain.xml('''<BooleanDomain name="bool" />''')
    def SetVectorize(self, vectorize):
        self.Vectorize = vectorize
        self.Modified()

    def RequestInformation(self, request, inInfo, outInfo):
        executive = self.GetExecutive()
        info = outInfo.GetInformationObject(0)
        info.Set(executive.WHOLE_EXTENT(),
                 [0, self.Dimensions[0] - 1, 0, self.Dimensions[1] - 1, 0, self.Dimensions[2] - 1], 6)
        return 1

    def double_gyre_vectorized(self, x, y, t):
        a = self.eps * np.sin(2.0 * np.pi / self.T * t)
        b = 1.0 - 2.0 * self.eps * np.sin(2.0 * np.pi / self.T * t)
        f = a * x**2 + b * x
        df = 2.0 * a * x + b

        dx = -np.pi * self.A * np.sin(np.pi * f) * np.cos(np.pi * y)
        dy = np.pi * self.A * np.cos(np.pi * f) * np.sin(np.pi * y) * df
        dt = np.full_like(t, 1.0)

        return np.stack([dx, dy, dt], axis=-1)

    def double_gyre_loop(self, x, y, t):
        dx = np.empty_like(x)
        dy = np.empty_like(y)
        dt = np.empty_like(t)

        for k in range(self.Dimensions[2]):
            for j in range(self.Dimensions[1]):
                for i in range(self.Dimensions[0]):
                    px = x[i, j, k]
                    py = y[i, j, k]
                    pt = t[i, j, k]

                    a = self.eps * np.sin(2.0 * np.pi / self.T * pt)
                    b = 1.0 - 2.0 * self.eps * np.sin(2.0 * np.pi / self.T * pt)
                    f = a * px ** 2 + b * px
                    df = 2.0 * a * px + b

                    dx[i, j, k] = -np.pi * self.A * np.sin(np.pi * f) * np.cos(np.pi * py)
                    dy[i, j, k] = np.pi * self.A * np.cos(np.pi * f) * np.sin(np.pi * py) * df
                    dt[i, j, k] = 1.0

        return np.stack([dx, dy, dt], axis=-1)

    def RequestData(self, request, inInfo, outInfo=None):
        xs = np.linspace(0.0, 2.0, self.Dimensions[0])
        ys = np.linspace(0.0, 1.0, self.Dimensions[1])
        ts = np.linspace(0.0, 10.0, self.Dimensions[2])
        x, y, t = np.meshgrid(xs, ys, ts, indexing='ij')

        output = vtkImageData.GetData(outInfo, 0)
        output.SetDimensions(self.Dimensions)
        output.SetSpacing(abs(xs[1] - xs[0]), abs(ys[1] - ys[0]), abs(ts[1] - ts[0]))
        output.SetOrigin(xs[0], ys[0], ts[0])
        output = dsa.WrapDataObject(output)

        if self.Vectorize:
            start = timer()
            v = self.double_gyre_vectorized(x, y, t)
            end = timer()
            print('Vectorized time: {}s'.format(end - start))
        else:
            start = timer()
            v = self.double_gyre_loop(x, y, t)
            end = timer()
            print('Loop time: {}s'.format(end - start))

        output.PointData.append(v.reshape((-1, 3), order='F'), 'v')
        output.VTKObject.GetPointData().SetActiveScalars('v')

        return 1
		