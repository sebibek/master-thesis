from vtk.util.vtkAlgorithm import *
from vtk.vtkCommonDataModel import vtkPolyData
from vtk.numpy_interface import dataset_adapter as dsa
from paraview.pythonalgorithm import smproxy, smhint, smdomain, smproperty


@smhint.xml('''<ShowInMenu category="My Submenu" />''')
@smproperty.input(name='Input', port_index=0)
@smproxy.filter(label='My Example Filter')
@smdomain.datatype(dataTypes=['vtkPolyData'])
class MyExampleFilter(VTKPythonAlgorithmBase):
    def __init__(self):
        VTKPythonAlgorithmBase.__init__(self, nInputPorts=1, nOutputPorts=1, outputType='vtkPolyData')
 
    def FillInputPortInformation(self, port, info):
        info.Set(self.INPUT_REQUIRED_DATA_TYPE(), 'vtkPolyData')
        return 1
 
    def RequestData(self, request, inInfo, outInfo):
        inp = vtkPolyData.GetData(inInfo[0], 0)
        output = vtkPolyData.GetData(outInfo, 0)
        if inp is None or output is None:
            return 0
        # copy input to output
        output.ShallowCopy(inp)
        # wrap input for numpy interface
        inp = dsa.WrapDataObject(inp)
        my_array = inp.PointData['Array']
        print(my_array)
        return 1

 
if __name__ == '__main__':
    from paraview.pythonalgorithm import get_plugin_xmls
    from xml.dom.minidom import parseString
    for xml in get_plugin_xmls(globals()):
        dom = parseString(xml)
    print(dom.toprettyxml(' ', '\n'))
	