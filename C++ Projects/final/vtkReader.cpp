// INLINE DEFINES

#ifdef WINDOWS
#include <direct.h>
#define GetCurrentDir _getcwd
#else
#include <unistd.h>
#define GetCurrentDir getcwd
#endif

#define MAXBUFSIZE  ((int) 1e8)
#define _USE_MATH_DEFINES


// INCLUDES (IMPORTS)
//#include <SFML/Graphics.hpp>
#include <iostream>
//#include "muParser.h"
#include <sstream>
#include <fstream>
#include <vector>
//#include <unistd.h>
#include <math.h>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
//#include <array>
#include <experimental/filesystem>
#include <string>
//#include<sys/io.h>
//#include <algorithm> 
//#include <cctype>
//#include <locale>
//#include <functional>
#include <numeric>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/tee.hpp>
#include <limits>
#include <ctime>
#include <omp.h>
#include <chrono>
typedef std::chrono::high_resolution_clock Clock;

#include <vtkSmartPointer.h>
#include <vtkProperty.h>
//#include <vtkDataSetMapper.h>
//#include <vtkImageActor.h>
//#include <vtkImageViewer2.h>
#include <vtkXMLImageDataReader.h>
//#include <vtkRenderWindow.h>
//#include <vtkRenderWindowInteractor.h>
//#include <vtkRenderer.h>
//#include <vtkXMLImageDataWriter.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkImageData.h>
#include <vtkPointData.h>
//#include <vtkDoubleArray.h>
#include <vtkCutter.h>
#include <vtkPlane.h>


// NAMESPACE IMPORTS
using namespace std;
using namespace Eigen;

typedef boost::iostreams::tee_device<std::ostream, std::ofstream> Tee;
typedef boost::iostreams::stream<Tee> TeeStream;
typedef std::numeric_limits< double > dbl;

double buffer[MAXBUFSIZE];

//definition of pi
const double pi = M_PI;

// PROTOTYPES
int width;
int height;
int steps; // use n steps for angular resolution - radres
double radres;

template <typename T>
T clip(const T& n, const T& lower, const T& upper); // template clip function

//length and width of window
int wSize = 701;

// SLICE EXTRACTOR
int main(int argc, char* argv[])
{
	// set slice (w) index
	int slice = 8;

	// set file name "brain.vti"
	std::string inputFilename = "brain.vti";

	// Read the file
	vtkSmartPointer<vtkXMLImageDataReader> reader =
		vtkSmartPointer<vtkXMLImageDataReader>::New();
	reader->SetFileName(inputFilename.c_str());
	reader->Update();
	cout << "size: " << reader->GetOutput()->GetScalarSize() << endl;
	
	// Create a plane to cut,here it cuts in the XY direction (xz normal=(1,0,0);XY =(0,0,1),YZ =(0,1,0)
	vtkSmartPointer<vtkPlane> plane =
		vtkSmartPointer<vtkPlane>::New();
	plane->SetOrigin(0, 0, slice);
	plane->SetNormal(0, 0, 1);

	// Create cutter to SLICE a 2D PLANE
	vtkSmartPointer<vtkCutter> cutter =
		vtkSmartPointer<vtkCutter>::New();
	cutter->SetCutFunction(plane);
	cutter->SetInputConnection(reader->GetOutputPort());
	cutter->Update();

	// get polyData of SLICE
	vtkSmartPointer<vtkPolyData> polyData = cutter->GetOutput(); //vtkSmartPointer<vtkImageData>::New();

	// create tensor array of slice and crop down to 2x2 matrices
	vtkSmartPointer < vtkDataArray> tensors = vtkDataArray::SafeDownCast(polyData->GetPointData()->GetArray("tensors"));// cutter->get()->GetScalars();
	cout << "array size: " << tensors->GetSize()/9 << endl;
	cout << "array width: " << sqrt(tensors->GetSize()/9) << endl;
	double tensor[9];
	int dim = 128* 128;
	Eigen::MatrixXd matrix(2, 2);
	std::vector<MatrixXd> matrixList;
	for (int i = 0; i < dim; i++)
	{
		int pointId = i;
		tensors->GetTuple(pointId, tensor);
		// CROP //
		matrix.row(0) << tensor[0], tensor[1]; // tensor[3]
		matrix.row(1) << tensor[3], tensor[4]; // tensor[5]
					  // tensor[6], tensor[7], // tensor[8]
		// APPEND //
		matrixList.push_back(matrix);
	}

	cout << "list size: " << matrixList.size() << endl;

	vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();

	writer->SetFileName("test.vtk");
	//writer->SetDataModeToAscii();
#if VTK_MAJOR_VERSION <= 5
	//writer->SetInputConnection(imageData->GetProducerPort());
#else
	writer->SetInputData(polyData);
#endif
	writer->Write();

	std::system("PaUsE");

	return EXIT_SUCCESS;
}