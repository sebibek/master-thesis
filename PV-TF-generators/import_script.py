import vtk
import numpy as np
from vtk.util.vtkAlgorithm import * # use (standalone) vtk module for sake of completeness
from vtk.vtkCommonDataModel import vtkImageData
from vtk.numpy_interface import dataset_adapter as dsa
from timeit import default_timer as timer
from vtk.util import numpy_support