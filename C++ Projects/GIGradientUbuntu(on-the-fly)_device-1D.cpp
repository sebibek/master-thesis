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

	__host__ __device__
		double operator()(const double& x, const double& y) const {
		return a*x + y; // performs y = a*x + y (SAXPY-OP)
	}
};

struct scale_functor
{
    const double a;

    scale_functor(double _a) : a(_a) {}

    __host__ __device__
        double operator()(const double& x) const { 
            return a * x; // performs y = a*x (scale-OP)
        }
};

struct multsum_functor
{
    multsum_functor() {}

    __host__ __device__
        double operator()(const double& x, const double& y) const { 
            return y * x; // performs y = a*x (scale-OP)
        }
};

struct mult_functor
{
    mult_functor() {}

    __host__ __device__
        double operator()(const double& x, const double& y) const { 
            return y * x; // performs y = a*x (scale-OP)
        }
};

struct abs_functor
{
    abs_functor() {}

    __host__ __device__
        double operator()(const double& x) const { 
            return abs(x); // performs y = a*x (scale-OP)
        }
};

void saxpy_fast(double a, thrust::device_vector<double>& X, const thrust::device_vector<double>::iterator& dstStart)
{
	// Y <- A * X + Y multiply (scale) vector and add another --> NEEDED
	thrust::transform(X.begin(), X.end(), dstStart, dstStart, saxpy_functor(a));
}


void scale_fast(double a, thrust::device_vector<double>& X)
{
	// Y <- A * X + Y multiply (scale) vector and add another --> NEEDED
	thrust::transform(X.begin(), X.end(), X.begin(), scale_functor(a));
}

thrust::device_vector<double> scale_fast(double a, const thrust::device_vector<double>::iterator& start,const thrust::device_vector<double>::iterator& end)
{
	thrust::device_vector<double> res(end-start, 0.0);
	// Y <- A * X + Y multiply (scale) vector and add another --> NEEDED
	thrust::transform(start, end, res.begin(), scale_functor(a));
	return res;
}

double multsum(const thrust::device_vector<double>::iterator& start,const thrust::device_vector<double>::iterator& end,const thrust::device_vector<double>::iterator& glyphStart)
{
	// Y <- A * X + Y multiply (scale) vector and add another --> NEEDED
	thrust::device_vector<double> res(end-start, 0.0);
	thrust::transform(start, end, glyphStart, res.begin(), multsum_functor());;
	return thrust::reduce(res.begin(), res.end());
}

void mult(const thrust::device_vector<double>::iterator& start,const thrust::device_vector<double>::iterator& end,const thrust::device_vector<double>::iterator& glyphStart)
{
	// Y <- A * X + Y multiply (scale) vector and add another --> NEEDED
	thrust::transform(start, end, glyphStart, start, mult_functor());
}

void mult(thrust::device_vector<double>& X, const thrust::device_vector<double>::iterator& start)
{
	//thrust::host_vector<double> res(X.size(), 0.0);
	// Y <- A * X + Y multiply (scale) vector and add another --> NEEDED
	thrust::transform(X.begin(), X.end(), start, X.begin(), mult_functor());
	//return res;
}

double accAbs(const thrust::device_vector<double>& X)
{
	// Y <- A * X + Y multiply (scale) vector and add another --> NEEDED	

	//thrust::device_vector<double> res(X.size(),0.0);
	return thrust::transform_reduce(X.begin(), X.end(), abs_functor(), 0.0, thrust::plus<double>());
	//return thrust::reduce(res.begin(), res.end());
}


// THRUST 

// THRUST 


// get current working directory to assign matrix.txt path
std::string GetCurrentWorkingDir(void) {
	char buff[FILENAME_MAX];
	GetCurrentDir(buff, FILENAME_MAX);
	std::string current_working_dir(buff);
	return current_working_dir;
}

// CLASS DEFINITIONS //

class neighborhood
{
	bool neighborR = true;
	bool neighborT = true;
	bool neighborL = true;
	bool neighborB = true;

public:

	neighborhood() // standard constructor
	{}
	//constructor
	neighborhood(int j, int i) // check the neighbourhood of (j,i) for missing neighbors..
	{
		// check if frame exceeded, outside FOV, offscreen..
		if (i == 0)
			neighborL = false; // left
		else if (i == width - 1)
			neighborR = false; // right
		if (j == 0)
			neighborT = false; // top
		else if (j == height - 1)
			neighborB = false; // bottom
		// leave order.. functional!
	}

	void change(int j, int i)
	{
		// RE-INITIALIZE
		neighborR = true;
		neighborT = true;
		neighborL = true;
		neighborB = true;
		// check if frame exceeded, outside FOV, offscreen..
		if (i == 0)
			neighborL = false; // left
		else if (i == width - 1)
			neighborR = false; // right
		if (j == 0)
			neighborT = false; // top
		else if (j == height - 1)
			neighborB = false; // bottom
		// leave order.. functional!
	}
	
	bool getR() { return neighborR; }
	bool getT() { return neighborT; }
	bool getL() { return neighborL; }
	bool getB() { return neighborB; }
};

class Pair
{
public:
	int jIndex = width / 2;
	int iIndex = width / 2;

	Pair() {}
	Pair(int j, int i)
	{
		jIndex = j;
		iIndex = i;
	}

	friend istringstream& operator>>(istringstream& stream, Pair& pair)
	{
		std::string token;
		int i = 0;
		while (getline(stream, token, ','))
		{
			if (i == 0)
				pair.jIndex = std::stoi(token);
			else
				pair.iIndex = std::stoi(token);
			i++;
		}


		return stream;
	}
};

// OPTION (cmd args) DEFINITIONS - getopt

// create options for getopt(.c)
std::string functionString = "2.1"; // Prototyping of functionString for strFunction (symbolic string arg parsed by muparser in cpp functional ptr rep)!
std::string workDir;
std::string record_folder = "frames";//directory to write images to, must exist
bool fullscreen = false; //fullscreen flag
bool total_anisotropy = false;
int record_frameskip = 0; // --> disables recording //recording frameskip
int ctrLimit = 0;
Pair lightSrcPos{ height / 2, width / 2 };
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
	Pair position{ height / 2, width / 2 };

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
				lightSrcPos = { width / 2, width / 2 };
				break;
			case 'p':
				std::cerr << "option -p requires argument, using default position (center-point)\n";
				lightSrcPos = { width / 2, width / 2 };
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
thrust::device_vector<double> operator-(const thrust::device_vector<double>& a, const thrust::device_vector<double>& b)
{
	assert(a.size() == b.size());

	thrust::device_vector<double> result(a.size());

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

	// REMEMBER: Grid size is half the amount of numbers for 2x2 matrices! CROPPING
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
	std::vector<bool> signs(3, false);
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
		if (sv1 == 0 || sv2 == 0 || sv1 / sv2 > 20.0) // for strong (total anisotropy)
		{
			// define bidirectional delta impulses
			glyphBuffer[i*steps+round(deg1*steps / 360)] = sv1;
			glyphBuffer[i*steps+static_cast<int>(round(deg1*steps / 360 + steps / 2)) % steps] = sv1;
			sum += 2 * sv1;
		}
		else // default case:
		{
			for (int j = 0; j < steps; j++) // sample ellipse equation for all steps
			{
				double val = dot / sqrt(sv2*sv2*cos(j*radres - deg1 * (M_PI / 180.0))*cos(j*radres - deg1 * (M_PI / 180.0)) + sv1 * sv1*sin(j*radres - deg1 * (M_PI / 180.0))*sin(j*radres - deg1 * (M_PI / 180.0))); //--> ellipse equation, evaluate for tMean (sum2)
				sum += val;
				glyphBuffer[i*steps + j] = val;
			}
		}
		double rMean = sum / steps; // compute rMean from cartesian (rectangular) energy-based integral as opposed to the polar integral relevant to the geometrical (triangular/circular) area
		// write glyphParameters
		glyphParameters.at(i).at(0) = 1.0 / rMean * sv1;
		glyphParameters.at(i).at(1) = 1.0 / rMean * sv2;

		// multiply respective cosine cone by valsum*=radres, because of energy normalization to cosine_sum (pre-computed in constructor)
		std::transform(glyphBuffer.begin() + i*steps, glyphBuffer.begin() +(i+1)*steps, glyphBuffer.begin() + i*steps, std::bind(std::multiplies<double>(), std::placeholders::_1, 1.0/rMean));
		//glyphBuffer.at(i) = 1.0 / rMean * glyphBuffer.at(i);
	}
}

int fast_mod(const int input, const int ceil) {
	// apply the modulo operator only when needed
	// (i.e. when the input is greater than the ceiling)
	return input >= ceil ? input % ceil : input;
	// NB: the assumption here is that the numbers are positive
}

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
	
	double meanA = 0.0;
	double cosine_sum = 0.0;
	// create member vectors (arrays) for storing the sampled directions theta
	thrust::host_vector<double> sampleBufferInit;
	thrust::device_vector<double> sampleBufferA;
	thrust::device_vector<double> sampleBufferB;
	thrust::device_vector<double>* glyphBuffer;
	
	std::vector<thrust::device_vector<double>> cosines;
	
	// create sample vector (dynamic)
	thrust::host_vector<double> initArray;
	thrust::device_vector<double> init;


	neighborhood hood; 
	bool flag = false;

	// create deltaIndex Map to access current neighbor cell for direction k[0..7]
	std::vector<int> deltaIndexSTL{ 1, 1 - width, -width, -1 - width, -1, -1 + width, width, 1 + width };
	thrust::host_vector<int> deltaIndex;

	thrust::host_vector<double> lightSrc;
	std::vector<thrust::host_vector<double>> lightSrcs;

	std::vector<thrust::device_vector<double>> weights;


public:
	propagator(const int dim, thrust::device_vector<double>* ellipseArray)
	{
		// initialize member samples w. 0
		initArray = thrust::host_vector<double>(steps, 0.0);
		init = initArray;
		lightSrc = initArray;
		deltaIndex = deltaIndexSTL;

		// assign ptrs to member vectors
		sampleBufferInit = thrust::host_vector<double>(dim*steps, 0.0);
		sampleBufferA = sampleBufferInit;
		sampleBufferB = sampleBufferInit;
		glyphBuffer = ellipseArray;

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
			scale_fast(1.0 / cosine_sum, cosines.at(k));

		cout.precision(dbl::max_digits10);
		cout << "cosine_sum: " << cosine_sum << endl;

		// construct a light src vector (delta functions for each sampled direction - normalized to total area (energy) of unit circle 1.0)
		for (int j = 0; j < steps; j++)
		{
			lightSrc[j] = steps;
			lightSrcs.push_back(lightSrc);
			lightSrc = initArray;
		}

		weights = std::vector<thrust::device_vector<double>>(8, initArray);
		// iterate through central directions array to distribute (spread) energy (intensity) to the cell neighbors
		for (int k = 0; k < 8; k++) // for each adjacent edge...
		{
			int midIndex = k * shiftIndex / 2;
			int index = fast_mod(k, 2) == 0 ? shiftIndex / 2 : betaIndex;

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
	}
	
	void propagate()
	{
		// 1 propagation cycle
		for(int j = 1; j < width -1; j++)
			for (int i = 1; i < width-1; i++) // for each node..
			{
				// define iterators for accessing current sample
				thrust::device_vector<double>::iterator start = sampleBufferA.begin() + (i + j * width) * steps;
				thrust::device_vector<double>::iterator end = sampleBufferA.begin() + (i + j * width + 1) * steps;

				//thrust::copy(start, end, read.begin()); // copy current sample to HOST (CPU) Vector for averaging (normalization)
				/*if (read == initArray)
					continue;*/
				//read = thrust::host_vector<double>(start, end); // copy current sample to HOST (CPU) Vector for averaging (normalization)
				if (thrust::equal(start,end,init.begin())) // test current sample on host
					continue;

				thrust::device_vector<double>::iterator glyphStart = glyphBuffer->begin() + (i + j * width) * steps;
				thrust::device_vector<double>::iterator glyphEnd = glyphBuffer->begin() + (i + j * width + 1) * steps;
			
				//glyph = thrust::host_vector<double>(glyphStart, glyphEnd); // copy current glyph to HOST (CPU) Vector for averaging (normalization)
				//glyph = glyphBuffer->at(i); // copy current glyph to HOST (CPU) Vector for averaging (normalization)

				// calculate means of I(phi) and T(phi)
				double sum1 = thrust::reduce(start, end); // THRUSTs accumulate analog starting w. sum = 0.0
				double sum2 = thrust::reduce(glyphStart, glyphEnd); // THRUSTs accumulate analog starting w. sum = 0.0
				mult(start,end,glyphStart);
				double sum3 = thrust::reduce(start, end);//thrust::transform(read.begin(), read.end(), glyph.begin(), out.begin(), thrust::multiplies<double>()); // perform read*glyph element-wise via trust transform method

				// compute iMean from cartesian (rectangular) energy-based integral as opposed to the polar integral relevant to the geometrical (triangular/circular) area
				double iMean = sum1 / steps; // -->tinc(dt) is a constant that can be drawn out of the integral
				// compute mean(T) from cartesian (rectangular) energy-based integral as opposed to the polar integral relevant to the geometrical (triangular/circular) area	
				double tMean = sum2 / steps; // -->tinc(dt) is a constant that can be drawn out of the integral
				// compute mean(T*I) from cartesian (rectangular) energy-based integral as opposed to the polar integral relevant to the geometrical (triangular/circular) area
				double tiMean = sum3 / steps; // -->tinc(dt) is a constant that can be drawn out of the integral
				// compute correction factor (scaling to mean=1, subsequent scaling to mean(I)), which follows energy conservation principles
				double cFactor = tiMean>0.0? tMean * iMean / tiMean : 1.0;

				// prepare readGlyphC for whole cell
				//readGlyphC = thrust::host_vector<double>(start,end); // crop subset of sample buffer (current sample): REMEMBER, always use constructor to extract subset of vectors!!!
				scale_fast(cFactor, start, end); // multiply cFactor*read*glyph -> readGlyphC yields cFactor*read*glyph for all steps

				double val_sum = 0.0;
				int delta = 0;
				int index = i + j*width;

				#pragma unroll // --> bringt nicht viel, aber stört auch nicht
				for (int k = 0; k < 8; k++) // for each adjacent edge...
				{
					delta = index + deltaIndex[k]; // compute index from deltaIndexMap (stores relative neighbor indices for all 8 directions)

					val_sum = multsum(start, end, weights.at(k).begin());

					thrust::device_vector<double>::iterator dstStart = sampleBufferB.begin() + (delta)*steps;
					saxpy_fast(val_sum *= radres, cosines.at(k), dstStart);
					meanA += val_sum;
				}		
		}
	}
	thrust::device_vector<double> propagateDist(int i, int j, int t)
	{
		// DUAL BUFFER PROPAGATION //
		thrust::fill(sampleBufferA.begin(), sampleBufferA.end(), 0.0);
		thrust::fill(sampleBufferB.begin(), sampleBufferB.end(), 0.0);
		double meanMem = 0.0;
		bool finished = false;
		//cout << "HERE1: " << endl;

		//lightSrc = lightSrcs.at(t);
		//thrust::copy(lightSrc.begin(), lightSrc.end(), sampleBufferA.begin() + steps*(j*width + i)); // WORKS
		int index = (j*width + i)*steps + t; // compute 1D index
		sampleBufferA[index] = steps;
		int ctr = 0;
		//cout << "HERE2: " << endl;
		// loop over nodes in grid and propagate until error to previous light distribution minimal <thresh
		while (!finished) // perform one single light propagation pass (iteration)
		{
			meanA = 0.0;

			this->propagate(); // propagate until finished..
			//	cout << "HERE2: " << endl;
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

		//cout << "ctr: " << ctr << endl;

		sampleBufferA[index] = 0.0; // remove light src as trivial difference --> if commented, symmetric fields yield NULL (homogeneous) FTLE fields
		//thrust::copy(initArray.begin(), initArray.end(), sampleBufferA.begin() + steps*(j*width + i)); // get pre-computed light src for current direction t
		//thrust::copy(sampleBufferA.begin(), sampleBufferA.end(), sampleBufferOut.begin());
		return sampleBufferA;//thrust::host_vector(sampleBufferA.begin(),sampleBufferA.end());
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

	thrust::host_vector<double> glyphBufferHost(dim*steps,0.0);
	
	std::vector<double> initGradient(3, 0.0); // construct 3D Gradient
	std::vector<std::vector<double>> deltaBuffer(dim*steps, initGradient); // initialize #steps 2D-planes w. empty glyphBuffer

	std::vector<std::vector<double>> glyphParameters(dim, std::vector<double>(3, 0.0));
	std::vector<std::vector<bool>> signMap(dim, std::vector<bool>(2, false)); // create a signMap relating normal force signs to singular values
	
	cout << "before compute glyphs" << endl;
	// compute Eigenframes/Superquadrics/Ellipses/Glyphs by calling computeGlyphs w. respective args
	computeGlyphs(glyphBufferHost, signMap, glyphParameters);
	thrust::device_vector<double> glyphBuffer(dim*steps);
	glyphBuffer = glyphBufferHost;
	double meanA = 0.0; // set up mean variables for threshold comparison as STOP criterion..

	// create propagator object (managing propagation, reprojection, correction, central directions, apertureAngles and more...)
	propagator prop(dim, &glyphBuffer);

	// DELTA (Gradient) COMPUTATION START //
	std::vector<double> gradient(3, 0.0); // dim3: x,y,theta

	thrust::device_vector<double> sampleBufferLeft;
	thrust::device_vector<double> sampleBufferRight;
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
				//cout << "HERE1" << endl;

				sampleBufferLeft = prop.propagateDist(i - 1, j, t); // propagate current lower distribution vector
				//cout << "HERE2" << endl;
				sampleBufferRight = prop.propagateDist(i + 1, j, t); // propagate current upper distribution vector
				// X-1D central differences.. VARIANT 1: bin by bin - spatial+directional distribution of energies
				meanA = accAbs(sampleBufferRight - sampleBufferLeft);
				gradient.at(0) = meanA / 2.0;

				

				sampleBufferLeft = prop.propagateDist(i, (j - 1), t); // propagate current distribution vector
				sampleBufferRight = prop.propagateDist(i, (j + 1), t); // propagate current distribution vector
				// Y-1D central differences..
				meanA = accAbs(sampleBufferRight - sampleBufferLeft);
				gradient.at(1) = meanA / 2.0;

				//cout << "HERE3" << endl;

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
