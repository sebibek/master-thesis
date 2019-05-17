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

// ZIP2 Iterator
typedef thrust::tuple<int, double> zip2Tuple; // --> 1 single element of zip iterator for DEREFERENCE!
typedef thrust::tuple<thrust::counting_iterator<int>, DoubleIterator> zip2IteratorTuple; // -->8 dbl iterators
typedef thrust::zip_iterator<zip2IteratorTuple> zip2Iterator;


typedef thrust::tuple<thrust::constant_iterator<double>, thrust::constant_iterator<double>, thrust::constant_iterator<double>, thrust::constant_iterator<double>, thrust::constant_iterator<double>, thrust::constant_iterator<double>, thrust::constant_iterator<double>, thrust::constant_iterator<double>> ConstDouble8IteratorTuple; // -->8 dbl iterators

typedef thrust::zip_iterator<ConstDouble8IteratorTuple> ConstDouble8Iterator;

struct elementMult : public thrust::binary_function<Double8Tuple,Double8Tuple,double>
{
    __host__ __device__
        Double8Tuple operator()(const Double8Tuple& a, const Double8Tuple& b)
        {
            return thrust::make_tuple(thrust::get<0>(a) * thrust::get<0>(b), thrust::get<1>(a) * thrust::get<1>(b), thrust::get<2>(a) * thrust::get<2>(b), thrust::get<3>(a) * thrust::get<3>(b), thrust::get<4>(a) * thrust::get<4>(b), thrust::get<5>(a) * thrust::get<5>(b), thrust::get<6>(a) * thrust::get<6>(b), thrust::get<7>(a) * thrust::get<7>(b));   
        }
};

struct elementPlus : public thrust::binary_function<Double8Tuple,Double8Tuple,double>
{
    __host__ __device__
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

struct propagateCell
{
	//width, steps, radres, weightsPtr, cosinesPtr, destination, initArrayD
    int width;
    int steps;
    double radres;
    double* initArrayD;
	int* deltaIndex;// { 1, 1 - width, -width, -1 - width, -1, -1 + width, width, 1 + width };
	double sums[8]; // 8 individual sums
	double* outputs; // 8 individual outputs
    
	double* weightsPtr;
	double* cosinesPtr;
	double** destinations; 
	double* read;
	double* readGlyph; #

	int shift[8];

public:
    propagateCell(int _width, int _steps, double _radres, double* _weightsPtr, double* _cosinesPtr, double** _destinations, double* _read, double* _readGlyph)
    {
		width = _width;
		steps = _steps;
		radres = _radres;
		outputs = new double[8 * steps];
		destinations = _destinations;
		cosinesPtr = _cosinesPtr;
		weightsPtr = _weightsPtr;
		read = _read;
		readGlyph = _readGlyph;
		deltaIndex[0] = 1;
		deltaIndex[1] = 1 - width;
		deltaIndex[2] = -width;
		deltaIndex[3] = -1 - width;
		deltaIndex[4] = -1;
		deltaIndex[5] = -1 + width;
		deltaIndex[6] = width;
		deltaIndex[7] = 1 + width;

    }

   template <typename Tuple>
    __host__ __device__ 
   
    void operator()(Tuple& t)
    {
		// get current cell index
    	int index = thrust::get<0>(t);
    	
    	if (index/width == 0 || index % width == 0 || index / width == height - 1 || index % width == width-1)
    		return;

		int offset = index * steps;
		// run for loop for equality check
		for (int t = 0; t < steps; t++) // for each step..
		{
			if (read[offset + t] != 0) // break on value..
				break;
			if (t == steps - 1)
				return; // return on last null value
		}

		double iMean = 0.0;
		double tiMean = 0.0;
		for (int t = 0; t < steps; t++) // for each step..
		{
			// compute iMean from cartesian (rectangular) energy-based integral as opposed to the polar integral relevant to the geometrical (triangular/circular) area
			iMean += read[offset+t]; // -->tinc(dt) is a constant that can be drawn out of the integral
			// compute mean(T*I) from cartesian (rectangular) energy-based integral as opposed to the polar integral relevant to the geometrical (triangular/circular) area
			tiMean = readGlyph[offset+t]; // -->tinc(dt) is a constant that can be drawn out of the integral
		}
		iMean = iMean / steps;
		tiMean = iMean / steps;
		// compute mean(T) from cartesian (rectangular) energy-based integral as opposed to the polar integral relevant to the geometrical (triangular/circular) area	
		double tMean = thrust::get<1>(t); // -->tinc(dt) is a constant that can be drawn out of the integral
		// compute correction factor (scaling to mean=1, subsequent scaling to mean(I)), which follows energy conservation principles
		double cFactor = tiMean > 0.0 ? tMean * iMean / tiMean : 1.0; // = tMean * iMean / tiMean for uncommented equality check

		// iterate through central directions array to distribute (spread) energy (intensity) to the cell neighbors
		for (int k = 0; k < 8; k++) // for each adjacent edge...
		{
			int midIndex = k * shiftIndex / 2;
			int shift = k % 2 == 0 ? shiftIndex / 2 : betaIndex;

			double val_sum = 0.0;

			for (int t = midIndex - shift; t <= midIndex + shift; t++) // for each step (along edge)..
			{
				int deltaJ = t - midIndex;
				int t_index = t < 0 ? t + steps : t % steps; // cyclic value permutation in case i exceeds the full circle degree 2pi

				double val = cFactor * readGlyph[offset + t_index]*weightsPtr[k*steps+t_index];

				val_sum += val; // val*radres
			}

			val_sum *= radres; // convert luminous intensity to power [W/sr->W]
			int delta = index + deltaIndex[k]; // compute index from deltaIndexMap (stores relative neighbor indices for all 8 directions)
			int offsetDst = delta * steps;

			for (int t = 0; t < steps; t++) // for each step..
				destinations[k][offsetDst + t] += cosines[k*steps + t] * val_sum;
		}
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
	thrust::device_vector<double> sampleBufferA;
	thrust::device_vector<double> sampleBufferB;
	thrust::device_vector<double> readGlyph;
	thrust::device_vector<double> glyphBuffer;
	
	thrust::device_vector<double> cosines;
	
	thrust::host_vector<int> lowerIndex;
	thrust::host_vector<int> upperIndex;

	// make glyph means vector
	thrust::device_vector<double> tMeans;

	// create init vector (dynamic)
	thrust::host_vector<double> initArray;
	thrust::device_vector<double> initArrayD;

	// create deltaIndex Map to access current neighbor cell for direction k[0..7]
	std::vector<int> deltaIndexSTL{ 1, 1 - width, -width, -1 - width, -1, -1 + width, width, 1 + width };
	thrust::host_vector<int> deltaIndex;

	thrust::host_vector<double> lightSrc;
	std::vector<thrust::host_vector<double>> lightSrcs;

	thrust::host_vector<double> weights;

	// 8 sample buffers for each direction k
	thrust::host_vector<double> sampleBuffer0;
	thrust::host_vector<double> sampleBuffer1;
	thrust::host_vector<double> sampleBuffer2;
	thrust::host_vector<double> sampleBuffer3;
	thrust::host_vector<double> sampleBuffer4;
	thrust::host_vector<double> sampleBuffer5;
	thrust::host_vector<double> sampleBuffer6;
	thrust::host_vector<double> sampleBuffer7;

	// make zip2 iterators for directly accessed values
	zip2Iterator zipFirst;
	zip2Iterator zipLast;

	// make counting iterators
	thrust::counting_iterator<int> itFirst(0);
	thrust::counting_iterator<int> itLast;

	double* weightsPtr;
	double* cosinesPtr;// = thrust::raw_pointer_cast(dev_ptr);
	double* read;// = thrust::raw_pointer_cast(dev_ptr);
	double* readGlyph;// = thrust::raw_pointer_cast(dev_ptr);

	double** destinations;


public:

	propagator(const int dim, thrust::host_vector<double> ellipseArray)
	{
		// assign ptrs to member vectors
		sampleBufferInit = thrust::host_vector<double>(dim*steps, 0.0);
		readGlyph = sampleBufferInit;
		sampleBufferA = sampleBufferInit;
		sampleBufferB = sampleBufferInit;
		glyphBuffer = ellipseArray;
		read = thrust::raw_pointer_cast(sampleBufferA.data());
		readGlyph = thrust::raw_pointer_cast(readGlyph.data())

		// initialize member samples w. 0
		initArray = thrust::host_vector<double>(steps, 0.0);
		initArrayD = initArray;// thrust::host_vector<double>(steps, 0.0);
		lightSrc = initArray;
		
		// create whole buffer for each of the 8 branches (directions)
		sampleBuffer0 = sampleBufferInit;
		sampleBuffer1 = sampleBufferInit;
		sampleBuffer2 = sampleBufferInit;
		sampleBuffer3 = sampleBufferInit;
		sampleBuffer4 = sampleBufferInit;
		sampleBuffer5 = sampleBufferInit;
		sampleBuffer6 = sampleBufferInit;
		sampleBuffer7 = sampleBufferInit;

		// cast raw ptrs for destination buffers
		destinations[0] = thrust::raw_pointer_cast(sampleBuffer0.data());
		destinations[1] = thrust::raw_pointer_cast(sampleBuffer1.data());
		destinations[2] = thrust::raw_pointer_cast(sampleBuffer2.data());
		destinations[3] = thrust::raw_pointer_cast(sampleBuffer3.data());
		destinations[4] = thrust::raw_pointer_cast(sampleBuffer4.data());
		destinations[5] = thrust::raw_pointer_cast(sampleBuffer5.data());
		destinations[6] = thrust::raw_pointer_cast(sampleBuffer6.data());
		destinations[7] = thrust::raw_pointer_cast(sampleBuffer7.data());

		tMeans = thrust::host_vector<double>(dim, 0.0);

		deltaIndex = deltaIndexSTL;

		cosines = thrust::host_vector<double>(8*steps, 0.0);

		cosinesPtr = thrust::raw_pointer_cast(weights.data());
		for (int k = 0; k < 8; k++) // for each node/neighbor..
		{
			int midIndex = (k * pi / 4) / radres;
			int offset = k * steps;
			for (int t = midIndex - shiftIndex; t <= midIndex + shiftIndex; t++) // for each step (along edge)..
			{
				int deltaJ = t - midIndex;
				int t_index = t < 0 ? t+steps : t % steps;

				double res = clip(cos((t_index - midIndex) * radres), 0.0, 1.0);
				cosine_sum += res * radres;
				cosines[t_index + offset] = res;
			}
		}
		cosine_sum = cosine_sum / 8.0;
		for (int k = 0; k < 8; k++) // for each node..
			thrust::transform(cosines.begin(), cosines.end(), thrust::make_constant_iterator(1.0/cosine_sum), cosines.at(k).begin(), thrust::multiplies<double>());
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
		
		weights = thrust::host_vector<double>(8 * steps, 0.0);

		weightsPtr = thrust::raw_pointer_cast(weights.data());
		// iterate through central directions array to distribute (spread) energy (intensity) to the cell neighbors
		for (int k = 0; k < 8; k++) // for each adjacent edge...
		{
			int midIndex = k * shiftIndex / 2;
			int index = fast_mod(k, 2) == 0 ? shiftIndex / 2 : betaIndex;
			int offset = k * steps;
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

				weights[t_index + offset] = val;
			}
		}

		// 1 propagation cycle
		for(int j = 1; j < width -1; j++)
			for (int i = 1; i < width-1; i++) // for each node..
			{
				int index = i + j * width; // compute 1D grid index
				glyphStart = std::next(glyphBuffer.begin(), index *steps);
				glyphEnd = std::next(glyphStart, steps);

				tMeans[index] = thrust::reduce(glyphStart, glyphEnd) / steps;
			}

		// assign counting iterator limit
		itLast = itFirst + dim;

		// make zip iterator of values directly accessed !
		zipFirst = thrust::make_zip_iterator(thrust::make_tuple(itFirst, tMeans.begin()));
		zipLast = thrust::make_zip_iterator(thrust::make_tuple(itLast, tMeans.end()));
		
	}

	void propagate()
	{
		thrust::transform(sampleBufferA.begin(), sampleBufferA.end(), glyphBuffer.begin(), readGlyph.begin(), thrust::multiplies<double>()); // perform read*glyph element-wise via thrust transform method
		
		// 1 propagation cycle
		thrust::for_each(zipFirst, zipLast, propagateCell(width, steps, radres, weightsPtr, cosinesPtr, destinations, read, readGlyph));

		// add up 8 individual sample buffers
		thrust::transform(destinationsFirst, destinationsLast, sampleBufferB.begin(), elementAcc());

		thrust::fill(sampleBuffer0.begin(), sampleBuffer0.end(), 0.0);
		thrust::fill(sampleBuffer1.begin(), sampleBuffer1.end(), 0.0);
		thrust::fill(sampleBuffer2.begin(), sampleBuffer2.end(), 0.0);
		thrust::fill(sampleBuffer3.begin(), sampleBuffer3.end(), 0.0);
		thrust::fill(sampleBuffer4.begin(), sampleBuffer4.end(), 0.0);
		thrust::fill(sampleBuffer5.begin(), sampleBuffer5.end(), 0.0);
		thrust::fill(sampleBuffer6.begin(), sampleBuffer6.end(), 0.0);
		thrust::fill(sampleBuffer7.begin(), sampleBuffer7.end(), 0.0);
		
		meanA = thrust::reduce(sampleBufferB.begin(), sampleBufferB.end(), 0.0)*radres;
		
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

			if (abs(meanA - meanMem) < thresh)
				finished = true;
			meanMem = meanA;

			thrust::fill(sampleBufferB.begin(), sampleBufferB.end(), 0.0);
			//sampleBufferB = sampleBufferInit;
			//ctr++;
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
