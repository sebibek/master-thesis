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
//#include <algorithm> 
//#include <cctype>
//#include <locale>
//#include <functional>
#include <numeric>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/tee.hpp>
#include <limits>
#include <ctime>

// VTK Includes
//#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkDoubleArray.h>

//thrust
#include <thrust/device_vector.h>
#include <thrust/transform.h>
#include <thrust/copy.h>
#include <thrust/fill.h>
#include <thrust/functional.h>
#include <thrust/transform_reduce.h>


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

// GENERIC FUNCTION DEFINITIONS

// THRUST's parallelized GPU methods
struct saxpy_functor
{
	const double a;

	saxpy_functor(double _a) : a(_a) {}

	__host__ //__device__
		double operator()(const double& x, const double& y) const {
		return a*x + y; // performs y = a*x + y (SAXPY-OP)
	}
};

struct scale_functor
{
    const double a;

    scale_functor(double _a) : a(_a) {}

    __host__ //__device__
        double operator()(const double& x) const { 
            return a * x; // performs y = a*x (scale-OP)
        }
};

struct abs_functor
{
    abs_functor() {}

    __host__ //__device__
        double operator()(const double& x) const { 
            return abs(x); // performs y = a*x (scale-OP)
        }
};

void saxpy_fast(double a, thrust::host_vector<double>& X, const thrust::host_vector<double>::iterator& dstStart)
{
	// Y <- A * X + Y multiply (scale) vector and add another --> NEEDED
	thrust::transform(X.begin(), X.end(), dstStart, dstStart, saxpy_functor(a));
}


void scale_fast(double a, thrust::host_vector<double>& X)
{
	// Y <- A * X + Y multiply (scale) vector and add another --> NEEDED
	thrust::transform(X.begin(), X.end(), X.begin(), scale_functor(a));
}

thrust::host_vector<double> scale_fast(double a, const thrust::host_vector<double>::iterator& start,const thrust::host_vector<double>::iterator& end)
{
	thrust::host_vector<double> res(end-start, 0.0);
	// Y <- A * X + Y multiply (scale) vector and add another --> NEEDED
	thrust::transform(start, end, res.begin(), scale_functor(a));
	return res;
}

double multsum(const thrust::host_vector<double>::iterator& start,const thrust::host_vector<double>::iterator& end,const thrust::host_vector<double>::iterator& glyphStart)
{
	thrust::host_vector<double> res(end-start, 0.0);
	// Y <- A * X + Y multiply (scale) vector and add another --> NEEDED
	thrust::transform(start, end, glyphStart, res.begin(), thrust::multiplies<double>());
	return thrust::reduce(res.begin(), res.end());
}

double accAbs(const thrust::host_vector<double>& X)
{
	// Y <- A * X + Y multiply (scale) vector and add another --> NEEDED	

	//thrust::device_vector<double> res(X.size(),0.0);
	return thrust::transform_reduce(X.begin(), X.end(), abs_functor(), 0.0, thrust::plus<double>());
	//return thrust::reduce(res.begin(), res.end());
}

// THRUST 


// get current working directory to assign matrix.txt path
std::string GetCurrentWorkingDir(void) {
	char buffer[FILENAME_MAX];
	GetCurrentDir(buffer, FILENAME_MAX);
	std::string current_working_dir(buffer);
	return current_working_dir;
}


// OPTION (cmd args) DEFINITIONS - getopt

// create options for getopt(.c)
std::string functionString = "2.1"; // Prototyping of functionString for strFunction (symbolic string arg parsed by muparser in cpp functional ptr rep)!
std::string workDir;
std::string record_folder = "frames";//directory to write images to, must exist
bool fullscreen = false; //fullscreen flag
bool total_anisotropy = false;
int record_frameskip = 0; // --> disables recording //recording frameskip
int ctrLimit = 0;
//Pair lightSrcPos{ height / 2, width / 2 };
double intensity = 2.1; // --> initial intensity val
double thresh = 0.001;

// parse files
void parse_file(char* filename)
{
	std::ifstream f(filename);
	std::string line;

	//default function attributes
	//sf::Color color = sf::Color::Cyan;
	std::string func_literal = "100*cos(theta)";
	//Pair position{ height / 2, width / 2 };

	if (f.is_open()) {

		while (std::getline(f, line)) {

			//ignore comments
			if (line[0] == '#')
				continue;
			//write function to vector
			/*else if (line == "end") {
				funcs.push_back(func_literal);
				//reset default values
				func_literal = "100*cos(theta)";
				positions.push_back(position);
			}*/
			
			//parse other statements
			else {
				//grab keyword
				std::string::size_type pos;
				pos = line.find(' ', 0);
				std::string tag = line.substr(0, pos);

				//get function literal
				if (tag == "function")
					func_literal = line.substr(pos + 1);
				//check for fullscreen
				else if (line == "fullscreen") {
					fullscreen = true;
				}
				//line color
				/*else if (tag == "color") {
					std::stringstream s;
					s << line.substr(pos + 1);
					int r = 0, g = 0, b = 0;
					s >> r >> g >> b;
					color = sf::Color(r, g, b);
				}*/
				// enter light src position
				/*else if (tag == "pos") 
				{
					std::stringstream s;
					s << line.substr(pos + 1);
					std::string str;
					s >> str;
					std::istringstream(str) >> position;
				}*/
				// thresh
				else if (tag == "thresh") {
					std::stringstream s;
					s << line.substr(pos + 1);
					s >> thresh;
				}
				// total anisotropy permission
				else if (tag == "total_anisotropy") {
					std::stringstream s;
					s << line.substr(pos + 1);
					s >> total_anisotropy;
				}
				// steps
				else if (tag == "steps") {
					std::stringstream s;
					s << line.substr(pos + 1);
					s >> steps;
					radres = 2 * pi / steps;
				}
				// steps
				else if (tag == "ctrLimit") {
					std::stringstream s;
					s << line.substr(pos + 1);
					s >> ctrLimit;
				}
				//window/graph size
				else if (tag == "window_size") {
					std::stringstream s;
					s << line.substr(pos + 1);
					s >> wSize;
				}
				else if (tag == "record_frameskip") {
					std::stringstream s;
					s << line.substr(pos + 1);
					s >> record_frameskip >> record_folder;
				}
			}
		}

		f.close();
	}
	else
		std::cerr << filename << " is not a valid filename.\n";
}

void parse_options(int argc, char* argv[]) {

	int c;
	std::string frameskip_opt = "-1";
	std::string size_opt = "-1";
	std::string directory_opt = "";
	std::string light_opt = "-1";
	std::string intensity_opt = "-1";

	//	extern char *optarg;
	//	extern int optind, optopt;

	//using getopts 
	while ((c = getopt(argc, argv, "fs:d:l:i:r:")) != -1) // --> classic getopt call: argc: # of args; argv: array of args(strings); "": optstring (registered options, followed by colon when taking args itself) 
	{
		int f = -1, s = -1;

		switch (c)
		{
		// correct option use
		/*case 'l': {
			light_opt.assign(optarg);
			std::istringstream(light_opt) >> lightSrcPos;
			positions.push_back(lightSrcPos);
			funcs.push_back(std::to_string(intensity));
			break;
		}*/
		case 'i': {
			intensity_opt.assign(optarg);
			std::istringstream(intensity_opt) >> intensity;
			break;
		}
		case 'r':
			frameskip_opt.assign(optarg);
			std::istringstream(frameskip_opt) >> f;
			if (f <= -1) {
				std::cerr << "invalid argument \'" << frameskip_opt << "\' to option -r\n";
				optind--;
			}
			record_frameskip = f > 0 ? f : 0;
			break;
		case 'f': {
			fullscreen = true;
			break; }
		case 's':
			size_opt.assign(optarg);
			std::istringstream(size_opt) >> s;
			if (s <= 0) {
				std::cerr << "invalid argument \'" << size_opt << "\' to option -s\n";
				optind--;
			}
			wSize = s > 0 ? s : wSize;
			break;
		case 'd':
			directory_opt.assign(optarg);
			record_folder = directory_opt;
			break;
		case ':': // missing option use..
			switch (optopt) {
			case 's':
				std::cout << "option -s requires argument, using default size 700\n";
				break;
			case 'r':
				std::cerr << "using default frameskip of 0 for option -r\n";
				break;
			case 'd':
				std::cerr << "option -d requires argument, disabling recording.\n";
				record_frameskip = -1;
				break;
			case 'l':
				std::cerr << "option -l requires argument, using default position (center-point)\n";
				//lightSrcPos = { width / 2, width / 2 };
				break;
			case 'p':
				std::cerr << "option -p requires argument, using default position (center-point)\n";
				//lightSrcPos = { width / 2, width / 2 };
				break;
			}
			break;
		case '?':
			std::cerr << "Unrecognized option: '-" << optopt << "\n";
			break;
		}
	}
	int optmem = optind; // set up option index memory variable
	for (int i = optind; i < argc; i++)
		parse_file(argv[i]); // parse "filename.txt" passed as LAST cmd line arg
	if (optind == optmem) // if optind did not change.., use standard config.txt in workDir for configuration
	{
		std::string str = workDir + "config.txt";
		char* cstr = new char[str.length() + 1];
		strncpy(cstr, str.c_str(), str.length() + 1);
		parse_file(cstr); // parse config.txt in workDir
		free(cstr);
	}
}


// OPERATOR DEFINITIONS //
//template <typename T> // element-wise plus for std::vector
//std::vector<T> operator+(const std::vector<T>& a, const std::vector<T>& b)
//{
//	std::vector<T> result(b.size());
//
//	for (int i = 0; i < b.size(); i++)
//		result.at(i) = a.at(i) + b.at(i);
//
//	return result;
//}

// vector PLUS(+) and MINUS(-)
template <typename T>
thrust::host_vector<double> operator+(const thrust::host_vector<double>& a, const thrust::host_vector<double>& b)
{
	assert(a.size() == b.size());

	thrust::host_vector<double> result(a.size());;

	//thrust::transform(a.begin(), a.end(), b.begin(), thrust::back_inserter(result), thrust::plus<T>());
	thrust::transform(a.begin(), a.end(), b.begin(), result.begin(), thrust::plus<double>());
	return result;
}

//template <typename T>
thrust::host_vector<double> operator-(const thrust::host_vector<double>& a, const thrust::host_vector<double>& b)
{
	assert(a.size() == b.size());

	thrust::host_vector<double> result(a.size());

	//thrust::transform(a.begin(), a.end(), b.begin(), thrust::back_inserter(result), thrust::minus<T>());
	thrust::transform(a.begin(), a.end(), b.begin(), result.begin(), thrust::minus<double>());
	return result;
}

//template <typename T> // element-wise multiplication for number, std::vector
//std::vector<T> operator*(const T& a, const std::vector<T>& b)
//{
//	std::vector<T> result(b.size());
//
//	for (int i = 0; i < b.size(); i++)
//		result.at(i) = a * b.at(i);
//
//	return result;
//}

//template <typename T> // element-wise minus for vector of vectors
std::vector<thrust::host_vector<double>> operator-(const std::vector<thrust::host_vector<double>>& a, const std::vector<thrust::host_vector<double>>& b)
{
	std::vector<thrust::host_vector<double>> result(b.size());

	for (int i = 0; i < b.size(); i++)
		result.at(i) = a.at(i) - b.at(i);

	return result;
}



//template <typename T> // element-wise plus for std::vector
//std::vector<std::vector<T>> operator*(const T a, const std::vector<std::vector<T>>& b)
//{
//	std::vector<std::vector<T>> result(b.size());
//
//	for (int i = 0; i < b.size(); i++)
//		result.at(i) = a * b.at(i);
//
//	return result;
//}

// SPECIFIC (MATHEMATICAL-PROGRAM USE) FUNCTIONS

// template function for evaluating string functions via ptr 
//template <typename T>
//T strFunction(T theta)
//{
//	// define function parser (muparser)
//	T t = theta; // ->theta: set initial theta
//	mu::Parser parser;
//	parser.DefineConst("pi", pi);
//	parser.DefineVar("theta", &t);
//	parser.SetExpr(functionString);
//	parser.DefineFun(_T("clip"), clipNeg, false);
//	parser.DefineFun(_T("delta"), deltaFunction, false);
//
//	T y = parser.Eval(); // evaluate parser expression
//
//	return y;
//}

template <typename T>
T clip(const T& n, const T& lower, const T& upper) // template clip function
{
	return std::max(lower, std::min(n, upper));
}

string trim(const string& str)
{
	// strip left (begin)
	size_t first = str.find_first_not_of(' ');
	if (string::npos == first) // if empty.., return empty str
		return str;
	// strip right (end)
	size_t last = str.find_last_not_of(' ');
	return str.substr(first, (last - first + 1)); // return cropped substring
}

MatrixXd readMatrix(std::string filepath, int* colsCount, int* rowsCount)
{
	int cols = 0, rows = 0;
	//double buff[MAXBUFSIZE];
	ifstream infile(filepath);

	while (!infile.eof())
	{
		string line;
		getline(infile, line);

		if (line.empty())
			break;

		int temp_cols = 0;
		stringstream stream(trim(line)); // parse stripped (trimmed) line w. stringstream
		while (!stream.eof())
			stream >> buffer[cols*rows + temp_cols++];

		if (temp_cols == 0) // if empty line
			continue;

		if (cols == 0) // set col count from first line
			cols = temp_cols;

		rows++;
	}

	infile.close();
	*colsCount = cols;
	*rowsCount = rows;
	MatrixXd result(rows, cols);
	for (int j = 0; j < rows; j++) // use (j, i)-index for data matrices, use (i, j) for mathematical matrices (w. apt Transpose/Transformations etc.)
		for (int i = 0; i < cols; i++)
			result(j, i) = buffer[cols*j + i];

	return result;
};

void computeGlyphs(thrust::host_vector<double>& glyphBuffer, std::vector<std::vector<bool>>& signMap, std::vector<std::vector<double>>& glyphParameters)
{
	int cols = 0; // create ptr to cols of txt tensor field
	int rows = 0; // create ptr to rows of txt tensor field
	std::string workDir = GetCurrentWorkingDir();
	MatrixXd m = readMatrix(workDir + "/matrix.txt", &cols, &rows);
	int dim = (rows / 2) * (cols / 2); // determine # of dimensions of grid for buffer (string/coefficient etc..) vectors

	// block the read-out matrices into a list of MatrixXd types (hardcoded 2x2 read-out)
	std::vector<MatrixXd> matrixList(dim, MatrixXd::Ones(2, 2));

	// REMEMBER: Grid size is half the amount of numbers for 2x2 matrices!
	for (int i = 0; i < rows / 2; i++)
		for (int j = 0; j < cols / 2; j++)
		{
			MatrixXd a = m.block<2, 2>(2 * i, 2 * j);
			matrixList.at(j + i * (cols / 2)) = a;
			//cout << "Matrix a: " << a << std::endl;
		}

	// compute the SVD (singular value decomposition)of all matrices in matrixList into svdList
	std::vector<JacobiSVD<MatrixXd>> svdList(dim, JacobiSVD<MatrixXd>(matrixList.at(0), ComputeThinU | ComputeThinV));
	for (int i = 0; i < dim; i++)
	{
		// SVD-based
		JacobiSVD<MatrixXd> svd(matrixList.at(i), ComputeThinU | ComputeThinV);
		svdList.at(i) = svd;
	}

	// define parameters
	double radres = (2 * pi) / steps;

	std::complex<double> sigma1(0, 0);
	std::complex<double> sigma2(0, 0);
	std::vector<bool> signs(2, false);

	thrust::host_vector<double>::iterator glyphStart = glyphBuffer.begin();
	thrust::host_vector<double>::iterator glyphEnd = glyphBuffer.begin();
	std::advance(glyphEnd, steps);
	// iterate through the matrixList/svdList (grid) and construct (scaled) ellipses in polar form (function) from the repsective singular values/vectors
	for (int i = 0; i < matrixList.size(); i++)
	{
		double y1 = svdList.at(i).matrixU().col(0)[1]; // use x - coordinate of both semi-axes -- Get LEFT U-vector
		double x1 = svdList.at(i).matrixU().col(0)[0]; // use x - coordinate of both semi-axes
		double y2 = svdList.at(i).matrixU().col(1)[1]; // use x - coordinate of both semi-axes -- Get RIGHT U-vector
		double x2 = svdList.at(i).matrixU().col(1)[0]; // use x - coordinate of both semi-axes
		double xx = matrixList.at(i).row(0)[0]; // "sigma_xx"
		double xy = matrixList.at(i).row(0)[1]; // "sigma_xy"
		double yx = matrixList.at(i).row(1)[1]; // "sigma_yx"
		double yy = matrixList.at(i).row(1)[1]; // "sigma_yy"
		double deg1 = atan2(y1, x1) * 180.0 / M_PI; // use vector atan2 to get rotational angle (phase) of both basis vectors in [-180°,180°]
		double deg2 = atan2(y2, x2) * 180.0 / M_PI; // use vector atan2 to get rotational angle (phase) of both basis vectors [-180°,180°]

		glyphParameters.at(i).at(2) = deg1;

		// calculate principal stresses w. formula.. https://vergleichsspannung.de/vergleichsspannungen/normalspannungshypothese-nh/herleitung-der-hauptspannungen/
		sigma1 = std::complex<double>(0.5*(xx + yy), 0) + 0.5*sqrt(std::complex<double>((xx - yy)*(xx - yy) + 4 * xy*yx, 0));
		sigma2 = std::complex<double>(0.5*(xx + yy), 0) - 0.5*sqrt(std::complex<double>((xx - yy)*(xx - yy) + 4 * xy*yx, 0));

		// crop sign of real part.. https://www.cg.tuwien.ac.at/research/vis/seminar9596/2-topo/evinter.html (rotational (complex) part already present in transformation)
		if (abs(sigma1 - sigma2) < 0) // check order of principal stresses sigma to match corresponding singular values/vectors
		{
			signs.at(0) = sigma2.real() < 0 ? true : false;
			signs.at(1) = sigma1.real() < 0 ? true : false;
		}
		else
		{
			signs.at(0) = sigma1.real() < 0 ? true : false;
			signs.at(1) = sigma2.real() < 0 ? true : false;
		}

		signMap.at(i) = signs; // assign singular value signs in sign map in decreasing order at position i
		// shift (normalize) degs from [-180°,180°] into the interval [0°,360°] - "circular value permutation"
		deg1 = deg1 < 0 ? 360 + deg1 : deg1;
		deg2 = deg2 < 0 ? 360 + deg2 : deg2;

		// singular values, decreasing order, corresponding singular vector order, scale ellipses axes in corresponding directions..
		double sv1 = svdList.at(i).singularValues()[0];
		double sv2 = svdList.at(i).singularValues()[1];
		double dot = sv2 * sv1;

		double sum = 0.0;
		if (sv1 == 0 || sv2 == 0 || sv1 / sv2 > 20.0) // if total anisotropy, needed to hit (match) the indices corresponding to glyph orientation
		{
			glyphStart[round(deg1*steps / 360)] = sv1;
			glyphStart[static_cast<int>(round(deg1*steps / 360 + steps / 2)) % steps] = sv1;
			sum += 2 * sv1;
		}
		else
		{
			for (int j = 0; j < steps; j++) // sample ellipse equation for all steps
			{
				double val = dot / sqrt(sv2*sv2*cos(j*radres - deg1 * (M_PI / 180.0))*cos(j*radres - deg1 * (M_PI / 180.0)) + sv1 * sv1*sin(j*radres - deg1 * (M_PI / 180.0))*sin(j*radres - deg1 * (M_PI / 180.0))); //--> ellipse equation, evaluate for tMean (sum2)
				sum += val;
				glyphStart[j] = val;
			}
		}
		double rMean = sum / steps; // compute rMean from cartesian (rectangular) energy-based integral as opposed to the polar integral relevant to the geometrical (triangular/circular) area
		// write glyphParameters
		glyphParameters.at(i).at(0) = 1.0 / rMean * sv1;
		glyphParameters.at(i).at(1) = 1.0 / rMean * sv2;

		// multiply respective cosine cone by valsum*=radres, because of energy normalization to cosine_sum (pre-computed in constructor)
		thrust::transform(glyphStart, glyphEnd, thrust::make_constant_iterator(1.0 / rMean), glyphStart, thrust::multiplies<double>());
		//std::transform(glyphStart, glyphEnd, glyphStart, std::bind(std::multiplies<double>(), std::placeholders::_1, 1.0 / rMean));
		//glyphBuffer.at(i) = 1.0 / rMean * glyphBuffer.at(i);

		thrust::advance(glyphStart, steps);
		thrust::advance(glyphEnd, steps);
	}
}

int fast_mod(const int input, const int ceil) {
	// apply the modulo operator only when needed
	// (i.e. when the input is greater than the ceiling)
	return input >= ceil ? input % ceil : input;
	// NB: the assumption here is that the numbers are positive
}

typedef thrust::host_vector<int>::iterator IntIterator;
typedef thrust::host_vector<double>::iterator DoubleIterator;;

typedef thrust::tuple<DoubleIterator, DoubleIterator, DoubleIterator,DoubleIterator, DoubleIterator, DoubleIterator,DoubleIterator, DoubleIterator> Double8IteratorTuple; // -->8 dbl iterators
typedef thrust::zip_iterator<Double8IteratorTuple> Double8Iterator; // -> zipped


// We'll use a 8-tuple to store our 8d vector type
typedef thrust::tuple<double,double,double,double,double,double,double,double> Double8Tuple; // --> 1 single element of zip iterator for DEREFERENCE!


// define zip iterator types

/*
typedef thrust::tuple<Double8Iterator, Double8Iterator, Double8Iterator, Double8Iterator, Double8Iterator> zip5IteratorTuple;
typedef thrust::constant_iterator<zip5IteratorTuple> zip5Iterator;

typedef thrust::tuple<DoubleIterator, DoubleIterator, DoubleIterator> zip3IteratorTuple;
typedef thrust::constant_iterator<zip3IteratorTuple> zip3Iterator;


typedef thrust::tuple<IntIterator, zip3Iterator, DoubleIterator, zip5Iterator> zip4IteratorTuple; // -->8 dbl iterators
typedef thrust::zip_iterator<zip4IteratorTuple> zip4Iterator;
*/

typedef thrust::tuple<int, DoubleIterator, Double8IteratorTuple, DoubleIterator, double> zip5Tuple; // --> 1 single element of zip iterator for DEREFERENCE!
typedef thrust::tuple<thrust::counting_iterator<int>, thrust::constant_iterator<DoubleIterator>, thrust::constant_iterator<Double8IteratorTuple>, thrust::constant_iterator<DoubleIterator>, DoubleIterator> zip5IteratorTuple; // -->8 dbl iterators
typedef thrust::zip_iterator<zip5IteratorTuple> zip5Iterator;


/*typedef thrust::tuple<DoubleIterator, DoubleIterator, DoubleIterator, DoubleIterator, DoubleIterator, DoubleIterator, DoubleIterator, DoubleIterator> Double8IteratorTuple;
typedef hrust::zip_iterator<ConstDouble8IteratorTuple> Double8IteratorTuple;*/

typedef thrust::tuple<thrust::constant_iterator<double>, thrust::constant_iterator<double>, thrust::constant_iterator<double>, thrust::constant_iterator<double>, thrust::constant_iterator<double>, thrust::constant_iterator<double>, thrust::constant_iterator<double>, thrust::constant_iterator<double>> ConstDouble8IteratorTuple; // -->8 dbl iterators

typedef thrust::zip_iterator<ConstDouble8IteratorTuple> ConstDouble8Iterator;
/*struct reduceFunctor
{
    __host__ __device__
        Double8 operator()(Double8Iterator firstOut, Double8Iterator lastOut)
        {
        	return thrust::make_tuple(thrust::reduce(firstOut.get_iterator_tuple().get<0>(), firstOut.get_iterator_tuple().get<0>()),thrust::reduce(firstOut.get_iterator_tuple().get<1>(), firstOut.get_iterator_tuple().get<1>()));
        	//return thrust::make_tuple(thrust::reduce(firstOut[0],lastOut[0]),thrust::reduce(firstOut[1],lastOut[1]),thrust::reduce(firstOut[2],lastOut[2]),thrust::reduce(firstOut[3],lastOut[3]),thrust::reduce(firstOut[4],lastOut[4]),thrust::reduce(firstOut[5],lastOut[5]),thrust::reduce(firstOut[6],lastOut[6]),thrust::reduce(firstOut[7],lastOut[7]));
        }
};*/

struct elementMult : public thrust::binary_function<Double8Tuple,Double8Tuple,double>
{
    __host__ //__device__
        Double8Tuple operator()(const Double8Tuple& a, const Double8Tuple& b)
        {
            return thrust::make_tuple(thrust::get<0>(a) * thrust::get<0>(b), thrust::get<1>(a) * thrust::get<1>(b), thrust::get<2>(a) * thrust::get<2>(b), thrust::get<3>(a) * thrust::get<3>(b), thrust::get<4>(a) * thrust::get<4>(b), thrust::get<5>(a) * thrust::get<5>(b), thrust::get<6>(a) * thrust::get<6>(b), thrust::get<7>(a) * thrust::get<7>(b));   
        }
};

struct elementPlus : public thrust::binary_function<Double8Tuple,Double8Tuple,double>
{
    __host__ //__device__
        Double8Tuple operator()(const Double8Tuple& a, const Double8Tuple& b)
        {
            return thrust::make_tuple(thrust::get<0>(a) + thrust::get<0>(b), thrust::get<1>(a) + thrust::get<1>(b), thrust::get<2>(a) + thrust::get<2>(b), thrust::get<3>(a) + thrust::get<3>(b), thrust::get<4>(a) + thrust::get<4>(b), thrust::get<5>(a) + thrust::get<5>(b), thrust::get<6>(a) + thrust::get<6>(b), thrust::get<7>(a) + thrust::get<7>(b));   
        }
};

struct elementAcc : public thrust::unary_function<Double8Tuple,double>
{
    __host__ //__device__
        double operator()(const Double8Tuple& a)
        {
            return thrust::get<0>(a) + thrust::get<1>(a) + thrust::get<2>(a) + thrust::get<3>(a) + thrust::get<4>(a) + thrust::get<5>(a) + thrust::get<6>(a) + thrust::get<7>(a);
        }
};

/*struct propagateCell : public thrust::binary_function<zipTuple7,zipTuple7,double>
{
    __host__ __device__
        Double8 operator()(const zipTuple7& a, const zipTuple7& b)
        {
            return thrust::make_tuple(thrust::get<0>(a) + thrust::get<0>(b), thrust::get<1>(a) + thrust::get<1>(b), thrust::get<2>(a) + thrust::get<2>(b), thrust::get<3>(a) + thrust::get<3>(b), thrust::get<4>(a) + thrust::get<4>(b), thrust::get<5>(a) + thrust::get<5>(b), thrust::get<6>(a) + thrust::get<6>(b), thrust::get<7>(a) + thrust::get<7>(b));   
        }
};*/

struct propagateCell
{ 
    int width;
    int steps;
    double radres;
    thrust::host_vector<double>* initArray;
    std::vector<int> deltaIndexSTL{ 1, 1 - width, -width, -1 - width, -1, -1 + width, width, 1 + width };
    //thrust::host_vector<double> sums = thrust::host_vector<double>(8, 0.0);
   	std::vector<Double8Iterator>* sums;
	std::vector<Double8Iterator>* outVector;
    
    DoubleIterator start;
    DoubleIterator end;
    DoubleIterator readGlyphStart;
    DoubleIterator readGlyphEnd;

    DoubleIterator dst0;
    DoubleIterator dst1;
    DoubleIterator dst2;
    DoubleIterator dst3;
    DoubleIterator dst4;
    DoubleIterator dst5;
    DoubleIterator dst6;
    DoubleIterator dst7;

	Double8Iterator firstReadGlyph;
	Double8Iterator lastReadGlyph; 
    Double8Iterator firstOut;
	Double8Iterator lastOut;
	Double8Iterator destinations;
	Double8Iterator sumIter;

	Double8Iterator* firstWeights;
	Double8Iterator* firstCosine;
	Double8Iterator* lastCosine;

	Double8Iterator firstDst;


public:
    propagateCell(int _width, int _steps, double _radres, Double8Iterator* _firstWeights, Double8Iterator* _firstCosine, Double8Iterator* _lastCosine, thrust::host_vector<double>* _initArray, std::vector<Double8Iterator>* _sums, std::vector<Double8Iterator>* _outVector)
    {
		width = _width;
		steps = _steps;
		radres = _radres;
		initArray = _initArray;
		//outVector = std::vector<thrust::host_vector<double>>(8, *initArray);
		firstWeights = _firstWeights;
		firstCosine = _firstCosine;
		lastCosine = _lastCosine;

		sums = _sums;
		outVector = _outVector;

		// make zip iterators of output iterators
		//firstOut = thrust::make_zip_iterator(thrust::make_tuple(outVector[0].begin(), outVector[1].begin(), outVector[2].begin(), outVector[3].begin(), outVector[4].begin(), outVector[5].begin(), outVector[6].begin(), outVector[7].begin()));// = thrust::get<8>(t);
		//lastOut = thrust::make_zip_iterator(thrust::make_tuple(outVector[0].end(), outVector[1].end(), outVector[2].end(), outVector[3].end(), outVector[4].end(), outVector[5].end(), outVector[6].end(), outVector[7].end())); //= thrust::get<9>(t);
		// make constant iterator for sums as scaling parameters // define zip iterator to write parallel sums for 8 branches
		//sumIter = thrust::make_constant_iterator(thrust::make_zip_iterator(thrust::make_tuple(sums.begin(), sums.begin()+1, sums.begin()+2, sums.begin()+3, sums.begin()+4, sums.begin()+5, sums.begin()+6, sums.begin()+7)));
		//sumIter = thrust::make_zip_iterator(thrust::make_tuple(thrust::make_constant_iterator(sums[0]), thrust::make_constant_iterator(sums[1]), thrust::make_constant_iterator(sums[2]), thrust::make_constant_iterator(sums[3]), thrust::make_constant_iterator(sums[4]), thrust::make_constant_iterator(sums[5]), thrust::make_constant_iterator(sums[6]), thrust::make_constant_iterator(sums[7])));
    }

   template <typename Tuple>
    __host__ //__device__ 
   
    void operator()(Tuple t)
    {

    	int index = thrust::get<0>(t);
    	
    	if (index/width == 0 || index % width == 0 || index / width == height - 1 || index % width == width-1)
    		return;

		// define iterators for accessing current sample
		start = thrust::get<1>(t);
		thrust::advance(start, index * steps);
		end = start;
		thrust::advance(end, steps);

		// check for trivial NULL-sample, if true, skip step (cell)
		if (thrust::equal(start, end, initArray->begin()))
			return;
		
		//outVector = thrust::get<5>(t)[index];
		//sums = thrust::get<6>(t)[index];

		firstOut = outVector[0][index];
		lastOut = firstOut;
		thrust::advance(lastOut, steps);
		sumIter = sums[0][index];

		firstReadGlyph = thrust::get<3>(t);
		thrust::advance(readGlyphStart, index * steps);
		lastReadGlyph = firstReadGlyph;
		thrust::advance(lastReadGlyph, steps); 
		
		readGlyphStart = thrust::get<3>(firstReadGlyph[0]);
		readGlyphEnd = readGlyphStart; // thrust::get<3>(lastReadGlyph[0]);
		thrust::advance(readGlyphEnd, steps);

		// compute iMean from cartesian (rectangular) energy-based integral as opposed to the polar integral relevant to the geometrical (triangular/circular) area
		double iMean = thrust::reduce(start, end) / steps; // -->tinc(dt) is a constant that can be drawn out of the integral
		// compute mean(T) from cartesian (rectangular) energy-based integral as opposed to the polar integral relevant to the geometrical (triangular/circular) area	
		double tMean = thrust::get<4>(t); // -->tinc(dt) is a constant that can be drawn out of the integral
		// compute mean(T*I) from cartesian (rectangular) energy-based integral as opposed to the polar integral relevant to the geometrical (triangular/circular) area
		double tiMean = thrust::reduce(readGlyphStart, readGlyphEnd) / steps; // -->tinc(dt) is a constant that can be drawn out of the integral
		// compute correction factor (scaling to mean=1, subsequent scaling to mean(I)), which follows energy conservation principles
		double cFactor = tiMean > 0.0 ? tMean * iMean / tiMean : 1.0;

		// prepare readGlyphC for whole cell
		thrust::transform(readGlyphStart, readGlyphEnd, thrust::make_constant_iterator(cFactor), readGlyphStart, thrust::multiplies<double>());

		//firstReadGlyph = thrust::make_zip_iterator(thrust::make_tuple(readGlyphStart, readGlyphStart, readGlyphStart, readGlyphStart, readGlyphStart, readGlyphStart, readGlyphStart, readGlyphStart));
		//lastReadGlyph = firstReadGlyph;//thrust::make_zip_iterator(thrust::make_tuple(readGlyphEnd, readGlyphEnd, readGlyphEnd, readGlyphEnd, readGlyphEnd, readGlyphEnd, readGlyphEnd, readGlyphEnd));
		
		
		// multiply weights to readGlyph parallel in 8 different branches (neighbors) contained in firstOut-->outVector
		thrust::transform(firstReadGlyph, lastReadGlyph, *firstWeights, firstOut, elementMult()); // evtl. TODO: Opt for smaller multiplication ranges: 8 hardcoded commands or functor constructor overload w. lowerIndex, uppperIndex

		// sum up 8 branches and write to sumVector to single double sum in 8 branches
		/*std::vector<double> stlSums{thrust::reduce(outVector[0].begin(), outVector[0].end()), thrust::reduce(outVector[1].begin(), outVector[1].end()), thrust::reduce(outVector[2].begin(), outVector[2].end()), thrust::reduce(outVector[3].begin(), outVector[3].end()), thrust::reduce(outVector[4].begin(), outVector[4].end()), thrust::reduce(outVector[5].begin(), outVector[5].end()), thrust::reduce(outVector[6].begin(), outVector[6].end()), thrust::reduce(outVector[7].begin(), outVector[7].end())};
		thrust::host_vector<double> sums = stlSums;*/
		/*sums[0] = thrust::reduce(outVector[0].begin(), outVector[0].end())*radres;
		sums[1] = thrust::reduce(outVector[1].begin() + lowerIndex[1], outVector[1].begin() + upperIndex[1])*radres;
		sums[2] = thrust::reduce(outVector[2].begin() + lowerIndex[2], outVector[2].begin() + upperIndex[2])*radres;
		sums[3] = thrust::reduce(outVector[3].begin() + lowerIndex[3], outVector[3].begin() + upperIndex[3])*radres;
		sums[4] = thrust::reduce(outVector[4].begin() + lowerIndex[4], outVector[4].begin() + upperIndex[4])*radres;
		sums[5] = thrust::reduce(outVector[5].begin() + lowerIndex[5], outVector[5].begin() + upperIndex[5])*radres;
		sums[6] = thrust::reduce(outVector[6].begin() + lowerIndex[6], outVector[6].begin() + upperIndex[6])*radres;
		sums[7] = thrust::reduce(outVector[7].begin() + lowerIndex[7], outVector[7].begin() + upperIndex[7])*radres;*/
		
		/*sums[0] = thrust::reduce(&thrust::get<0>(firstOut[0]), &thrust::get<0>(lastOut[0]))*radres;
		sums[1] = thrust::reduce(&thrust::get<1>(firstOut[0]), &thrust::get<1>(lastOut[0]))*radres;
		sums[2] = thrust::reduce(&thrust::get<2>(firstOut[0]), &thrust::get<2>(lastOut[0]))*radres;
		sums[3] = thrust::reduce(&thrust::get<3>(firstOut[0]), &thrust::get<3>(lastOut[0]))*radres;
		sums[4] = thrust::reduce(&thrust::get<4>(firstOut[0]), &thrust::get<4>(lastOut[0]))*radres;
		sums[5] = thrust::reduce(&thrust::get<5>(firstOut[0]), &thrust::get<5>(lastOut[0]))*radres;
		sums[6] = thrust::reduce(&thrust::get<6>(firstOut[0]), &thrust::get<6>(lastOut[0]))*radres;
		sums[7] = thrust::reduce(&thrust::get<7>(firstOut[0]), &thrust::get<7>(lastOut[0]))*radres;*/

		*thrust::get<0>(sumIter[0]) = thrust::reduce(&thrust::get<0>(firstOut[0]), &thrust::get<0>(lastOut[0]))*radres;
		*thrust::get<1>(sumIter[0]) = thrust::reduce(&thrust::get<1>(firstOut[0]), &thrust::get<1>(lastOut[0]))*radres;
		*thrust::get<2>(sumIter[0]) = thrust::reduce(&thrust::get<2>(firstOut[0]), &thrust::get<2>(lastOut[0]))*radres;
		*thrust::get<3>(sumIter[0]) = thrust::reduce(&thrust::get<3>(firstOut[0]), &thrust::get<3>(lastOut[0]))*radres;
		*thrust::get<4>(sumIter[0]) = thrust::reduce(&thrust::get<4>(firstOut[0]), &thrust::get<4>(lastOut[0]))*radres;
		*thrust::get<5>(sumIter[0]) = thrust::reduce(&thrust::get<5>(firstOut[0]), &thrust::get<5>(lastOut[0]))*radres;
		*thrust::get<6>(sumIter[0]) = thrust::reduce(&thrust::get<6>(firstOut[0]), &thrust::get<6>(lastOut[0]))*radres;
		*thrust::get<7>(sumIter[0]) = thrust::reduce(&thrust::get<7>(firstOut[0]), &thrust::get<7>(lastOut[0]))*radres;


		//thrust::get<5>(t) = thrust::reduce(sums.begin(), sums.end(), 0.0); // sum up sums to obtain mean energy for convergence criterion

		//cout << "sums: " << thrust::get<5>(t) << endl;
		//sumIter = thrust::make_zip_iterator(thrust::make_tuple(thrust::make_constant_iterator(sums[0]), thrust::make_constant_iterator(sums[1]), thrust::make_constant_iterator(sums[2]), thrust::make_constant_iterator(sums[3]), thrust::make_constant_iterator(sums[4]), thrust::make_constant_iterator(sums[5]), thrust::make_constant_iterator(sums[6]), thrust::make_constant_iterator(sums[7])));
		
		// make dst iterators
		DoubleIterator dst0 = thrust::get<0>(thrust::get<2>(t));
		thrust::advance(dst0, (index + deltaIndexSTL[0])*steps);
		DoubleIterator dst1 = thrust::get<1>(thrust::get<2>(t));
		thrust::advance(dst1, (index + deltaIndexSTL[1])*steps);
		DoubleIterator dst2 = thrust::get<2>(thrust::get<2>(t));
		thrust::advance(dst2, (index + deltaIndexSTL[2])*steps);
		DoubleIterator dst3 = thrust::get<3>(thrust::get<2>(t));
		thrust::advance(dst3, (index + deltaIndexSTL[3])*steps);
		DoubleIterator dst4 = thrust::get<4>(thrust::get<2>(t));
		thrust::advance(dst4, (index + deltaIndexSTL[4])*steps);
		DoubleIterator dst5 = thrust::get<5>(thrust::get<2>(t));
		thrust::advance(dst5, (index + deltaIndexSTL[5])*steps);
		DoubleIterator dst6 = thrust::get<6>(thrust::get<2>(t));
		thrust::advance(dst6, (index + deltaIndexSTL[6])*steps);
		DoubleIterator dst7 = thrust::get<7>(thrust::get<2>(t));
		thrust::advance(dst7, (index + deltaIndexSTL[7])*steps);

		// make dst zip iterator
		firstDst = thrust::make_zip_iterator(thrust::make_tuple(dst0, dst1, dst2, dst3, dst4, dst5, dst6, dst7)); // TODO: prepare firstDst as host::vector for every cell index, ALSO prepare start,end, readGlyph

		// VAR 1
		// scale respective cosines lobes w. constant iterator(ptr) to dbl sums in 8 branches
		thrust::transform(*firstCosine, *lastCosine, thrust::make_constant_iterator(sumIter), firstOut, elementMult());
		// add up contribution for respective neighbor in 8 branches
		thrust::transform(firstOut, lastOut, firstDst, firstDst, elementPlus());
    }
};

class propagator
{
	// set principal arture angles in rads
	double alpha = 36.8699 * M_PI / 180;
	double beta = 26.5651 * M_PI / 180;
	// define principal central directions array (12 central directions in 2D --> 30 in 3D, ufff..)
	//std::array<double, 12> centralDirections{ 3.0 / 2 * M_PI + 1.0172232666228471, 0, 0.5535748055013016, 1.0172232666228471, M_PI / 2, M_PI / 2 + 0.5535748055013016, M_PI / 2 + 1.0172232666228471, M_PI, M_PI + 0.5535748055013016, M_PI + 1.0172232666228471, 3.0 / 2 * M_PI, 3.0 / 2 * M_PI + 0.5535748055013016 };
	//std::array<double, 12> apertureAngles{ beta, alpha, beta, beta, alpha, beta, beta, alpha, beta, beta, alpha, beta }; // define aperture angles array in order

	// define function parser (muparser)
	int shiftIndex = steps / 4;
	int betaIndex = (beta) / radres;
	int centralIndex = (alpha / 2) / radres;
	int dim = width * height;
	int propCtr = 0;
	
	double meanA = 0.0;
	double cosine_sum = 0.0;
	// create member vectors (arrays) for storing the sampled directions theta
	thrust::host_vector<double> sampleBufferInit;
	thrust::host_vector<double> sampleBufferA;
	thrust::host_vector<double> sampleBufferB;
	thrust::host_vector<double> readGlyph;
	thrust::host_vector<double>* glyphBuffer;
	
	std::vector<thrust::host_vector<double>> cosines;
	
	thrust::host_vector<int> lowerIndex;
	thrust::host_vector<int> upperIndex;


	thrust::host_vector<double> tMeans;

	std::vector<thrust::host_vector<double>> sumsMem;
	
	std::vector<Double8Iterator> sums;
	std::vector<Double8Iterator> outVector;

	// create sample vector (dynamic)
	thrust::host_vector<double> initArray;

	// create deltaIndex Map to access current neighbor cell for direction k[0..7]
	std::vector<int> deltaIndexSTL{ 1, 1 - width, -width, -1 - width, -1, -1 + width, width, 1 + width };
	thrust::host_vector<int> deltaIndex;

	thrust::host_vector<double> lightSrc;
	std::vector<thrust::host_vector<double>> lightSrcs;

	std::vector<thrust::host_vector<double>> weights;
	//std::vector<thrust::host_vector<double>> diagWeights;

	thrust::host_vector<double> means;

	thrust::host_vector<double>::iterator start;
	thrust::host_vector<double>::iterator end;// = std::next(start, steps);
	thrust::host_vector<double>::iterator glyphStart;
	thrust::host_vector<double>::iterator glyphEnd;
	thrust::host_vector<double>::iterator readGlyphStart;
	thrust::host_vector<double>::iterator readGlyphEnd;
	thrust::host_vector<double>::iterator dstStart;

	thrust::host_vector<double> sampleBuffer0;
	thrust::host_vector<double> sampleBuffer1;
	thrust::host_vector<double> sampleBuffer2;
	thrust::host_vector<double> sampleBuffer3;
	thrust::host_vector<double> sampleBuffer4;
	thrust::host_vector<double> sampleBuffer5;
	thrust::host_vector<double> sampleBuffer6;
	thrust::host_vector<double> sampleBuffer7;

	std::vector<thrust::host_vector<double>> out0;
	std::vector<thrust::host_vector<double>> out1;
	std::vector<thrust::host_vector<double>> out2;
	std::vector<thrust::host_vector<double>> out3;
	std::vector<thrust::host_vector<double>> out4;
	std::vector<thrust::host_vector<double>> out5;
	std::vector<thrust::host_vector<double>> out6;
	std::vector<thrust::host_vector<double>> out7;

	thrust::host_vector<Double8Iterator> cellDestinations;

	zip5Iterator zipFirst;
	zip5Iterator zipLast;

	Double8Iterator firstWeights;
	Double8Iterator firstOut;

	Double8Iterator lastOut;
	Double8Iterator firstCosine;
	Double8Iterator lastCosine;

	Double8IteratorTuple bufferTuple;
	Double8Iterator destinationsFirst; 

	Double8Iterator destinationsLast; 


public:

	propagator(const int dim, thrust::host_vector<double>* ellipseArray)
	{
		// assign ptrs to member vectors
		sampleBufferInit = thrust::host_vector<double>(dim*steps, 0.0);
		readGlyph = sampleBufferInit;
		sampleBufferA = sampleBufferInit;
		sampleBufferB = sampleBufferInit;
		glyphBuffer = ellipseArray;

		// initialize member samples w. 0
		initArray = thrust::host_vector<double>(steps, 0.0);
		lightSrc = initArray;
		outVector = std::vector<Double8Iterator>(dim);
		sums = std::vector<Double8Iterator>(dim);

		sumsMem = std::vector<thrust::host_vector<double>>(dim, thrust::host_vector<double>(8,0.0));
		
		// create output sample buffer for each cell for each of the 8 branches (directions) for intermediate results
		out0 = std::vector<thrust::host_vector<double>>(dim, initArray);
		out1 = std::vector<thrust::host_vector<double>>(dim, initArray);
		out2 = std::vector<thrust::host_vector<double>>(dim, initArray);
		out3 = std::vector<thrust::host_vector<double>>(dim, initArray);
		out4 = std::vector<thrust::host_vector<double>>(dim, initArray);
		out5 = std::vector<thrust::host_vector<double>>(dim, initArray);
		out6 = std::vector<thrust::host_vector<double>>(dim, initArray);
		out7 = std::vector<thrust::host_vector<double>>(dim, initArray);

		// create whole buffer for each of the 8 branches (directions) for intensity accumulation
		sampleBuffer0 = thrust::host_vector<double>(dim*steps, 0.0);
		sampleBuffer1 = sampleBuffer0;
		sampleBuffer2 = sampleBuffer0;
		sampleBuffer3 = sampleBuffer0;
		sampleBuffer4 = sampleBuffer0;
		sampleBuffer5 = sampleBuffer0;
		sampleBuffer6 = sampleBuffer0;
		sampleBuffer7 = sampleBuffer0;

		tMeans = thrust::host_vector<double>(dim, 0.0);
		means = thrust::host_vector<double>(dim, 0.0);

		deltaIndex = deltaIndexSTL;

		for (int k = 0; k < 8; k++) // for each node..
		{
			int midIndex = (k * pi / 4) / radres;

			std::vector<double> cosK(steps, 0.0);
			for (int j = midIndex - shiftIndex; j <= midIndex + shiftIndex; j++) // for each step (along edge)..
			{
				int deltaJ = j - midIndex;
				int j_index = j < 0 ? j+steps : j % steps;

				double res = clip(cos((j_index - midIndex) * radres), 0.0, 1.0);
				cosine_sum += res * radres;
				cosK.at(j_index) = res;
			}
			cosines.push_back(cosK);
		}
		cosine_sum = cosine_sum / 8.0;
		for (int k = 0; k < 8; k++) // for each node..
			thrust::transform(cosines.at(k).begin(), cosines.at(k).end(), thrust::make_constant_iterator(1.0/cosine_sum), cosines.at(k).begin(), thrust::multiplies<double>());
			//std::transform(cosines.at(k).begin(), cosines.at(k).end(), cosines.at(k).begin(), std::bind(std::multiplies<double>(), std::placeholders::_1, 1.0/cosine_sum));

		cout.precision(dbl::max_digits10);
		cout << "cosine_sum: " << cosine_sum << endl;

		// construct a light src vector (delta functions for each sampled direction - normalized to total area (energy) of unit circle 1.0)
		for (int j = 0; j < steps; j++)
		{
			lightSrc[j] = steps;
			lightSrcs.push_back(lightSrc);
			lightSrc = initArray;
		}
		
		weights = std::vector<thrust::host_vector<double>>(8, initArray);
		// iterate through central directions array to distribute (spread) energy (intensity) to the cell neighbors
		for (int k = 0; k < 8; k++) // for each adjacent edge...
		{
			int midIndex = k * shiftIndex / 2;
			int index = fast_mod(k, 2) == 0 ? shiftIndex / 2 : betaIndex;

			lowerIndex.push_back(midIndex - index);
			upperIndex.push_back(midIndex + index + 1);
			// TODO: thrust OP multiplies 3 vectors --> first scale one w. cFactor (constant) then call thrust OP w. read,glyph and prepared weightVector (3 args)
			// --> construct thrust vector w. following loop by running 1 time in constructor RUN k loop in 8 sequential transform and mult calls (also iterator vectors needed)
			for (int t = midIndex - index; t <= midIndex + index; t++) // for each step (along edge)..
			{
				int deltaJ = t - midIndex;
				int t_index = t < 0 ? t + steps : t % steps; // cyclic value permutation in case i exceeds the full circle degree 2pi

				double val = 1.0;

				// split overlapping diagonal cones w.r.t to their relative angular area (obtained from face neighbors)..
				if ((abs(deltaJ) > centralIndex) && fast_mod(k, 2) == 0) // for alphas, use edge overlap > centralIndex
					if (abs(deltaJ) == shiftIndex / 2)
						val = 0.5*0.3622909908722584*val;
					else
						val = 0.3622909908722584*val;
				else if (fast_mod(k, 2) != 0) // for betas (diagonals), use static edge overlap-
					val = 0.6377090091277417*val;

				/*if (fast_mod(k, 2) == 0)
					faceWeights.at(k)[t_index] = val;
				else
					diagWeights.at(k)[t_index] = val;*/
				weights.at(k)[t_index] = val;
			}
		}

		cellDestinations = thrust::host_vector<Double8Iterator>(dim);
		// 1 propagation cycle
		for(int j = 1; j < width -1; j++)
			for (int i = 1; i < width-1; i++) // for each node..
			{
				int index = i + j * width; // compute 1D grid index
				glyphStart = std::next(glyphBuffer->begin(), index *steps);
				glyphEnd = std::next(glyphStart, steps);

				cellDestinations[index] = thrust::make_zip_iterator(thrust::make_tuple(sampleBuffer0.begin() + (index + deltaIndexSTL[0]) * steps, sampleBuffer1.begin() + (index + deltaIndexSTL[1]) * steps, sampleBuffer2.begin() + (index + deltaIndexSTL[2]) * steps, sampleBuffer3.begin() + (index + deltaIndexSTL[3]) * steps, sampleBuffer4.begin() + (index + deltaIndexSTL[4]) * steps, sampleBuffer5.begin() + (index + deltaIndexSTL[5]) * steps, sampleBuffer6.begin() + (index + deltaIndexSTL[6]) * steps, sampleBuffer7.begin() + (index + deltaIndexSTL[7]) * steps));

				outVector[index] = thrust::make_zip_iterator(thrust::make_tuple(out0[index].begin(), out1[index].begin(), out2[index].begin(), out3[index].begin(), out4[index].begin(), out5[index].begin(), out6[index].begin(), out7[index].begin())); 

				sums[index] = thrust::make_zip_iterator(thrust::make_tuple(sumsMem[index].begin(), sumsMem[index].begin()+1, sumsMem[index].begin()+2, sumsMem[index].begin()+3, sumsMem[index].begin()+4, sumsMem[index].begin()+5, sumsMem[index].begin()+6, sumsMem[index].begin()+7));
				//sums[index] = thrust::make_zip_iterator(thrust::make_tuple(thrust::make_constant_iterator(sumsMem[index].begin()), thrust::make_constant_iterator(sumsMem[index].begin()+1), thrust::make_constant_iterator(sumsMem[index].begin()+2), thrust::make_constant_iterator(sumsMem[index].begin()+3), thrust::make_constant_iterator(sumsMem[index].begin()+4), thrust::make_constant_iterator(sumsMem[index].begin()+5), thrust::make_constant_iterator(sumsMem[index].begin()+6), thrust::make_constant_iterator(sumsMem[index].begin())+7));
				
				tMeans[index] = thrust::reduce(glyphStart, glyphEnd) / steps;
			}

		// make zip iterators of weights
		firstWeights = thrust::make_zip_iterator(thrust::make_tuple(weights[0].begin(), weights[1].begin(), weights[2].begin(), weights[3].begin(), weights[4].begin(), weights[5].begin(), weights[6].begin(), weights[7].begin()));
		//auto lastWeights = boost::compute::make_zip_iterator(boost::make_tuple(weights[0).end(), weights[1).end(), weights[2).end(), weights[3).end(), weights[4).end(), weights[5).end(), weights[6).end(), weights[7).end()));
		// make zip iterators of destination iterators
		//firstOut = thrust::make_zip_iterator(thrust::make_tuple(outVector[0].begin(), outVector[1].begin(), outVector[2].begin(), outVector[3].begin(), outVector[4].begin(), outVector[5].begin(), outVector[6].begin(), outVector[7].begin()));
		//lastOut = thrust::make_zip_iterator(thrust::make_tuple(outVector[0].end(), outVector[1].end(), outVector[2].end(), outVector[3].end(), outVector[4].end(), outVector[5].end(), outVector[6].end(), outVector[7].end()));
		
		// make cosine zip iterators			
		firstCosine = thrust::make_zip_iterator(thrust::make_tuple(cosines[0].begin(), cosines[1].begin(), cosines[2].begin(), cosines[3].begin(), cosines[4].begin(), cosines[5].begin(), cosines[6].begin(), cosines[7].begin()));
		lastCosine = thrust::make_zip_iterator(thrust::make_tuple(cosines[0].end(), cosines[1].end(), cosines[2].end(), cosines[3].end(), cosines[4].end(), cosines[5].end(), cosines[6].end(), cosines[7].end()));

		bufferTuple = thrust::make_tuple(sampleBuffer0.begin(), sampleBuffer1.begin(), sampleBuffer2.begin(), sampleBuffer3.begin(), sampleBuffer4.begin(), sampleBuffer5.begin(), sampleBuffer6.begin(), sampleBuffer7.begin());

		destinationsFirst = thrust::make_zip_iterator(bufferTuple);
		destinationsLast = thrust::make_zip_iterator(thrust::make_tuple(sampleBuffer0.end(), sampleBuffer1.end(), sampleBuffer2.end(), sampleBuffer3.end(), sampleBuffer4.end(), sampleBuffer5.end(), sampleBuffer6.end(), sampleBuffer7.end()));
		
		Double8Iterator firstReadGlyph = thrust::make_zip_iterator(thrust::make_tuple(readGlyphStart, readGlyphStart, readGlyphStart, readGlyphStart, readGlyphStart, readGlyphStart, readGlyphStart, readGlyphStart));

		thrust::counting_iterator<int> itFirst(0);
		thrust::counting_iterator<int> itLast = itFirst + dim;

		zipFirst = thrust::make_zip_iterator(thrust::make_tuple(itFirst, thrust::make_constant_iterator(sampleBufferA.begin()), thrust::make_constant_iterator(bufferTuple), thrust::make_constant_iterator(firstReadGlyph), tMeans.begin()));
		zipLast = thrust::make_zip_iterator(thrust::make_tuple(itLast, thrust::make_constant_iterator(sampleBufferA.begin()), thrust::make_constant_iterator(bufferTuple), thrust::make_constant_iterator(firstReadGlyph), tMeans.end()));

		
	}

	void propagate()
	{
		thrust::transform(sampleBufferA.begin(), sampleBufferA.end(), glyphBuffer->begin(), readGlyph.begin(), thrust::multiplies<double>()); // perform read*glyph element-wise via thrust transform method
		
		cout << "Here 1" << endl;
		//zipFirst = thrust::make_zip_iterator(thrust::make_tuple(itFirst, sampleBufferA.begin(), sampleBufferB.begin(), readGlyph.begin(), tMeans.begin(), thrust::make_constant_iterator(zipper5)));
		//zipLast = thrust::make_zip_iterator(thrust::make_tuple(itLast, zipper3, tMeans.end(), thrust::make_constant_iterator(zipper5)));

		//typedef thrust::tuple<int,double,double,double,double, Double8Iterator, Double8Iterator, Double8Iterator, Double8Iterator, Double8Iterator> zipTuple10;

		//thrust::constant_iterator<DoubleIterator> iter = thrust::make_constant_iterator(sampleBufferA.begin());
		

		//propagateCell propagate(width,steps,radres);
		// 1 propagation cycle
		thrust::for_each(zipFirst, zipLast, propagateCell(width, steps, radres, &firstWeights, &firstCosine, &lastCosine, &initArray, &sums, &outVector)); // &cellDestinations
		cout << "Here 2" << endl;

		// add up 8 individual sample buffers
		thrust::transform(destinationsFirst, destinationsLast, sampleBufferB.begin(), elementAcc());
		cout << "Here 3" << endl;

		thrust::fill(sampleBuffer0.begin(), sampleBuffer0.end(), 0.0);
		thrust::fill(sampleBuffer1.begin(), sampleBuffer1.end(), 0.0);
		thrust::fill(sampleBuffer2.begin(), sampleBuffer2.end(), 0.0);
		thrust::fill(sampleBuffer3.begin(), sampleBuffer3.end(), 0.0);
		thrust::fill(sampleBuffer4.begin(), sampleBuffer4.end(), 0.0);
		thrust::fill(sampleBuffer5.begin(), sampleBuffer5.end(), 0.0);
		thrust::fill(sampleBuffer6.begin(), sampleBuffer6.end(), 0.0);
		thrust::fill(sampleBuffer7.begin(), sampleBuffer7.end(), 0.0);
		
		meanA = thrust::reduce(sampleBufferB.begin(), sampleBufferB.end(), 0.0)*radres;
		cout << "Here 4" << endl;
		
		/*if(propCtr < 10)
			cout << "dst sum: " << thrust::reduce(sampleBufferB.begin(), sampleBufferB.end(), 0.0)*radres << endl;
		*/
		

		
	}
	thrust::host_vector<double> propagateDist(int i, int j, int t)
	{
		// DUAL BUFFER PROPAGATION //
		thrust::fill(sampleBufferA.begin(), sampleBufferA.end(), 0.0);
		thrust::fill(sampleBufferB.begin(), sampleBufferB.end(), 0.0);
		double meanMem = 0.0;
		bool finished = false;
		//lightSrc = lightSrcs.at(t);
		//thrust::copy(lightSrc.begin(), lightSrc.end(), sampleBufferA.begin() + steps*(j*width + i)); // WORKS
		int index = (j*width + i)*steps + t; // compute 1D index
		sampleBufferA[index] = steps;
		int ctr = 0;
		// loop over nodes in grid and propagate until error to previous light distribution minimal <thresh
		while (!finished) // perform one single light propagation pass (iteration)
		{
			meanA = 0.0;
			this->propagate(); // propagate until finished..
			//meanA *= (1.0 / radres) / (steps*sampleBufferA.size());
			thrust::swap(sampleBufferA, sampleBufferB);
			//sampleBufferA = sampleBufferB;
			sampleBufferA[index] = steps;
			//thrust::copy(lightSrc.begin(), lightSrc.end(), sampleBufferA.begin() + steps*(j*width + i)); // get pre-computed light src for current direction t

			/*if(t==0 && ctr < 50)
			{
			cout << "meanMem: " << meanMem << endl;
			cout << "meanA: " << meanA << endl;
			}*/
			if (abs(meanA - meanMem) < thresh)
				finished = true;
			meanMem = meanA;

			thrust::fill(sampleBufferB.begin(), sampleBufferB.end(), 0.0);
			//sampleBufferB = sampleBufferInit;
			ctr++;
			/*if(ctr < 50)
			{
			cout << "meanMem: " << meanMem << endl;
			cout << "meanA: " << meanA << endl;
			propCtr++;
			}*/
		}
		
		

		sampleBufferA[index] = 0.0; // remove light src as trivial difference --> if commented, symmetric fields yield NULL (homogeneous) FTLE fields
		//thrust::copy(initArray.begin(), initArray.end(), sampleBufferA.begin() + steps*(j*width + i)); // get pre-computed light src for current direction t

		return sampleBufferA;
	}
};

//template <typename T>
double acc2(std::vector<thrust::host_vector<double>> vec)
{
	/*double sum = std::accumulate(m.begin(), m.end(), 0, [](auto lhs, const auto& rhs)
	{
		return std::accumulate(rhs.begin(), rhs.end(), lhs);
	});
	return sum;*/
	double sum = 0.0;
	for (int i = 0; i < vec.size(); i++)
		sum += accAbs(vec.at(i));
	
	return sum;
}

template<typename Iter_T>
double vectorNorm(Iter_T first, Iter_T last)
{
	return sqrt(inner_product(first, last, first, 0.0));
}

int main(int argc, char* argv[])
{	
	// parse input option file
	parse_options(argc, argv);

	// 2D Grid START //
	int cols = 0; // create cols of txt tensor field
	int rows = 0; // create rows of txt tensor field
	cout << "before matrix read" << endl;
	workDir = GetCurrentWorkingDir();
	MatrixXd m = readMatrix(workDir + "/matrix.txt", &cols, &rows); // call countMatrix to determine rows/cols count #
	width = cols / 2; // determine width of grid for correct indexing
	height = rows / 2;
	cout << "width|height|steps: " << width << "|" << height << "|" << steps << endl;
	const int dim = width * height; // determine # of dimensions of grid for buffer (string/coefficient etc..) vectors

	thrust::host_vector<double> glyphBuffer(dim*steps,0.0);
	
	std::vector<double> initGradient(3, 0.0); // construct 3D Gradient
	std::vector<std::vector<double>> deltaBuffer(dim*steps, initGradient); // initialize #steps 2D-planes w. empty glyphBuffer

	std::vector<std::vector<double>> glyphParameters(dim, std::vector<double>(3, 0.0));
	std::vector<std::vector<bool>> signMap(dim, std::vector<bool>(2, false)); // create a signMap relating normal force signs to singular values
	
	cout << "before compute glyphs" << endl;
	// compute Eigenframes/Superquadrics/Ellipses/Glyphs by calling computeGlyphs w. respective args
	computeGlyphs(glyphBuffer, signMap, glyphParameters);
	double meanA = 0.0; // set up mean variables for threshold comparison as STOP criterion..

	// create propagator object (managing propagation, reprojection, correction, central directions, apertureAngles and more...)
	propagator prop(dim, &glyphBuffer);

	// DELTA (Gradient) COMPUTATION START //
	std::vector<double> gradient(3, 0.0); // dim3: x,y,theta

	thrust::host_vector<double> sampleBufferLeft;
	thrust::host_vector<double> sampleBufferRight;
	cout << "before constructing gradient vector.." << endl;
	double duration; float total = 0.0;
	std::clock_t startTotal = std::clock();
	for (int t = 0; t < steps; t++)
	{
		cout << "before computing gradients for t: " << t << endl;
		double start = std::clock();

		for (int j = 0; j < height; j++)
			for (int i = 0; i < width; i++)
			{
				if (i == 0 || i == width - 1 || j == 0 || j == height - 1)
					continue;

				sampleBufferLeft = prop.propagateDist(i - 1, j, t); // propagate current lower distribution vector
				sampleBufferRight = prop.propagateDist(i + 1, j, t); // propagate current upper distribution vector
				// X-1D central differences.. VARIANT 1: bin by bin - spatial+directional distribution of energies
				meanA = accAbs(sampleBufferRight - sampleBufferLeft);
				gradient.at(0) = meanA / 2.0;

				sampleBufferLeft = prop.propagateDist(i, (j - 1), t); // propagate current distribution vector
				sampleBufferRight = prop.propagateDist(i, (j + 1), t); // propagate current distribution vector
				// Y-1D central differences..
				meanA = accAbs(sampleBufferRight - sampleBufferLeft);
				gradient.at(1) = meanA / 2.0;

				sampleBufferLeft = prop.propagateDist(i, j, (t == 0 ? (steps - 1) : (t - 1))); // propagate current distribution vector
				sampleBufferRight = prop.propagateDist(i, j, (t + 1) % steps); // propagate current distribution vector
				// t-1D central differences..
				meanA = accAbs(sampleBufferRight - sampleBufferLeft);
				gradient.at(2) = meanA / 2.0;//(2.0*radres) for normalization

				// X-1D central differences.. VARIANT 2: cell by cell - spatial distribution of energies (chose because of redundancy for same spatial distribution but differing directional distribution)
				//meanA = 0.0;
				//for(int k = 0; k < dim; k++)
				//	meanA += abs(std::accumulate(distBuffer.at(j*width + i + 1 + t * dim).at(k).begin(), distBuffer.at(j*width + i + 1 + t * dim).at(k).end(),0.0) - std::accumulate(distBuffer.at(j*width + i - 1 + t * dim).at(k).begin(), distBuffer.at(j*width + i - 1 + t * dim).at(k).end(),0.0));
				//gradient.at(0) = meanA / 2.0;
				//// Y-1D central differences..
				//meanA = 0.0;
				//for (int k = 0; k < dim; k++)
				//	meanA += abs(std::accumulate(distBuffer.at((j + 1)*width + i + t * dim).at(k).begin(), distBuffer.at((j + 1)*width + i + t * dim).at(k).end(),0.0) - std::accumulate(distBuffer.at((j - 1)*width + i + t * dim).at(k).begin(), distBuffer.at((j - 1)*width + i + t * dim).at(k).end(),0.0));
				//gradient.at(1) = meanA / 2.0;
				//// t-1D central differences..
				//meanA = 0.0;
				//for (int k = 0; k < dim; k++)
				//	meanA += abs(std::accumulate(distBuffer.at(j*width + i + (t + 1) % steps * dim).at(k).begin(), distBuffer.at(j*width + i + (t + 1) % steps * dim).at(k).end(),0.0) - std::accumulate(distBuffer.at(j *width + i + (t == 0 ? (steps - 1) : (t - 1)) * dim).at(k).begin(), distBuffer.at(j *width + i + (t == 0 ? (steps - 1) : (t - 1)) * dim).at(k).end(),0.0));
				//gradient.at(2) = meanA / 2.0;

				deltaBuffer.at(j*width + i + t * dim) = gradient;
				//cout << "HERE: 1 cycle passed" << endl;
			}
		duration = ((std::clock() - start)*1000.0 / (double)CLOCKS_PER_SEC);
		cout << "timer: " << duration << " ms" << endl;
		total += duration;
	}
	// DELTA (Gradient) COMPUTATION END //

	cout << "..after propagation TOTAL, total timer:" << total << " ms" << endl;

	glyphBuffer.clear();

	glyphBuffer.shrink_to_fit();
	
	// vector norm (Gradient) COMPUTATION START //
	std::vector<double> scalarNorm(width*height*steps, 0.0); // construct 3D Gradient Norm Vector

	vtkSmartPointer<vtkImageData> imageData = vtkSmartPointer<vtkImageData>::New();

	imageData->SetDimensions(width, height, steps);

	vtkSmartPointer<vtkDoubleArray> energy = vtkSmartPointer<vtkDoubleArray>::New();

	energy->SetNumberOfComponents(1);
	energy->SetNumberOfTuples(width * height * steps);

	cout << "before computing gradient (vector) norm.." << endl;
	startTotal = std::clock();
	int ctr = 0;
	for (int t = 0; t < steps; t++)
		for (int j = 0; j < height; j++)
			for (int i = 0; i < width; i++)
			{
				double res = vectorNorm(deltaBuffer.at(j*width + i + t * dim).begin(), deltaBuffer.at(j*width + i + t * dim).end());
				scalarNorm.at(j*width + i + t * dim) = res;
				energy->SetValue(ctr, res);
				ctr++;
			}
	
	// vector norm (Gradient) COMPUTATION END //
	
	duration = ((std::clock() - startTotal)*1000.0 / (double)CLOCKS_PER_SEC);
	cout << "..after, timer: " << duration << " ms" << endl;

	// VTK OUTPUT START //

	imageData->GetPointData()->AddArray(energy);
	energy->SetName("Energy");

	vtkSmartPointer<vtkXMLImageDataWriter> writer = vtkSmartPointer<vtkXMLImageDataWriter>::New();

	writer->SetFileName("test.vti");
	//writer->SetDataModeToAscii();
#if VTK_MAJOR_VERSION <= 5
	writer->SetInputConnection(imageData->GetProducerPort());
#else
	writer->SetInputData(imageData);
#endif
	writer->Write();

	std::system("PaUsE");

	return EXIT_SUCCESS;
}
