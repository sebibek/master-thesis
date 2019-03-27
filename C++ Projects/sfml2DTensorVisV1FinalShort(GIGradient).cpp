// INLINE DEFINES

#ifdef WINDOWS
#include <direct.h>
#define GetCurrentDir _getcwd
#else
#include <unistd.h>
#define GetCurrentDir getcwd
#endif

#define MAXBUFSIZE  ((int) 1e5)
#define _USE_MATH_DEFINES

// INCLUDES (IMPORTS)
#include <SFML/Graphics.hpp>
#include <iostream>
#include "muParser.h"
#include <sstream>
#include <fstream>
#include <vector>
#include <unistd.h>
#include <math.h>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <array>
#include <filesystem>
#include <string>
#include <algorithm> 
#include <cctype>
#include <locale>
#include <functional>
#include <numeric>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/tee.hpp>
#include <limits>
#include <ctime>
// VTK Includes
#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkXMLImageDataWriter.h>
#include <vtkImageData.h>
#include <vtkPointData.h>
#include <vtkDoubleArray.h>

// NAMESPACE IMPORTS
using namespace std;
using namespace Eigen;
using namespace mu;

typedef boost::iostreams::tee_device<std::ostream, std::ofstream> Tee;
typedef boost::iostreams::stream<Tee> TeeStream;
typedef std::numeric_limits< double > dbl;

//definition of pi
const double pi = M_PI;

// PROTOTYPES
int width;
int height;
int steps = 360; // use n steps for angular resolution - radres
double radres = 2 * pi / steps;

template <typename T>
T clip(const T& n, const T& lower, const T& upper); // template clip function

//length and width of window
int wSize = 701;
sf::Vector2i windowsize;

// GENERIC FUNCTION DEFINITIONS

// -> muparser custom clip function: template functions need to be static callback (call after) functions, which can be passed as arguments
static mu::value_type clipNeg(mu::value_type v) { return clip(v, 0.0, 1.0); }
static mu::value_type deltaFunction(mu::value_type v)
{
	if (abs(v) < radres/2.0 || v == radres/2.0) // use margin related to angular resolution radres..
		return steps; // normalization: should integrate total energy of unit circle into one single direction
	else
		return 0.0;
}
//writes frames to file 
//directory must already exist -- parse file methods
void record(sf::RenderTexture& in, int frameskip, std::string directory)
{
	// create frame index
	static int frame = 0;

	if (frameskip > 0 && frame % frameskip == 0) {
		in.display();
		sf::Texture outTx = in.getTexture();
		sf::Image outImg = outTx.copyToImage();
		std::stringstream ugh;
		ugh << directory << "/" << frame << ".png";
		outImg.saveToFile(ugh.str());
	}
	frame++;
}

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
void parse_file(char* filename, std::vector<std::string>& funcs, std::vector<Pair>& positions)
{
	std::ifstream f(filename);
	std::string line;

	//default function attributes
	sf::Color color = sf::Color::Cyan;
	std::string func_literal = "100*cos(theta)";
	Pair position{ height / 2, width / 2 };

	if (f.is_open()) {

		while (std::getline(f, line)) {

			//ignore comments
			if (line[0] == '#')
				continue;
			//write function to vector
			else if (line == "end") {
				funcs.push_back(func_literal);
				//reset default values
				func_literal = "100*cos(theta)";
				positions.push_back(position);
			}
			
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
				else if (tag == "color") {
					std::stringstream s;
					s << line.substr(pos + 1);
					int r = 0, g = 0, b = 0;
					s >> r >> g >> b;
					color = sf::Color(r, g, b);
				}
				// enter light src position
				else if (tag == "pos") 
				{
					std::stringstream s;
					s << line.substr(pos + 1);
					std::string str;
					s >> str;
					std::istringstream(str) >> position;
				}
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

void parse_options(int argc, char* argv[], std::vector<std::string>& funcs, std::vector<Pair>& positions) {

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
		case 'l': {
			light_opt.assign(optarg);
			std::istringstream(light_opt) >> lightSrcPos;
			positions.push_back(lightSrcPos);
			funcs.push_back(std::to_string(intensity));
			break;
		}
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
		parse_file(argv[i], funcs, positions); // parse "filename.txt" passed as LAST cmd line arg
	if (optind == optmem) // if optind did not change.., use standard config.txt in workDir for configuration
	{
		std::string str = workDir + "/config.txt";
		char* cstr = new char[str.length() + 1];
		strcpy_s(cstr, str.length() + 1, str.c_str());
		parse_file(cstr, funcs, positions); // parse config.txt in workDir
		free(cstr);
	}
}

// OPERATOR DEFINITIONS //
template <typename T> // element-wise plus for std::vector
std::vector<T> operator+(const std::vector<T>& a, const std::vector<T>& b)
{
	std::vector<T> result(b.size());

	for (int i = 0; i < b.size(); i++)
		result.at(i) = a.at(i) + b.at(i);

	return result;
}

template <typename T> // element-wise plus for std::vector
std::vector<T> operator*(const T a, const std::vector<T>& b)
{
	std::vector<T> result(b.size());

	for (int i = 0; i < b.size(); i++)
		result.at(i) = a * b.at(i);

	return result;
}

template <typename T> // element-wise plus for std::vector
std::vector<std::vector<T>> operator-(const std::vector<std::vector<T>>& a, const std::vector<std::vector<T>>& b)
{
	std::vector<std::vector<T>> result(b.size());

	for (int i = 0; i < b.size(); i++)
		result.at(i) = a.at(i) + (-1.0)*b.at(i);

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
template <typename T>
T strFunction(T theta)
{
	// define function parser (muparser)
	T t = theta; // ->theta: set initial theta
	mu::Parser parser;
	parser.DefineConst("pi", pi);
	parser.DefineVar("theta", &t);
	parser.SetExpr(functionString);
	parser.DefineFun(_T("clip"), clipNeg, false);
	parser.DefineFun(_T("delta"), deltaFunction, false);

	T y = parser.Eval(); // evaluate parser expression

	return y;
}

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
	double buff[MAXBUFSIZE];
	ifstream infile(filepath);

	while (!infile.eof())
	{
		string line;
		getline(infile, line);

		int temp_cols = 0;
		stringstream stream(trim(line)); // parse stripped (trimmed) line w. stringstream
		while (!stream.eof())
			stream >> buff[cols*rows + temp_cols++];

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
			result(j, i) = buff[cols*j + i];

	return result;
};

void computeGlyphs(std::vector<std::vector<double>>& glyphBuffer, std::vector<std::vector<bool>>& signMap)
{
	int cols = 0; // create ptr to cols of txt tensor field
	int rows = 0; // create ptr to rows of txt tensor field
	std::string workDir = GetCurrentWorkingDir();
	MatrixXd m = readMatrix(workDir + "/matrix.txt", &cols, &rows);
	int dim = rows / 2 * cols / 2; // determine # of dimensions of grid for buffer (string/coefficient etc..) vectors

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

	std::vector<bool> signs(2, false);
	// iterate through the matrixList/svdList (grid) and construct (scaled) ellipses in polar form (function) from the repsective singular values/vectors
	for (int i = 0; i < matrixList.size(); i++)
	{
		double y1 = svdList.at(i).matrixU().col(0)[1]; // use x - coordinate of both semi-axes 
		double x1 = svdList.at(i).matrixU().col(0)[0]; // use x - coordinate of both semi-axes "sigma_xx"
		double y2 = svdList.at(i).matrixU().col(1)[1]; // use x - coordinate of both semi-axes 
		double x2 = svdList.at(i).matrixU().col(1)[0]; // use x - coordinate of both semi-axes "sigma_xx"
		double xx = matrixList.at(i).row(0)[0]; // "sigma_xx"
		double yy = matrixList.at(i).row(1)[1]; // "sigma_yy"
		double deg1 = atan2(y1, x1) * 180.0 / M_PI; // use vector atan2 to get rotational angle (phase) of both basis vectors in [-180°,180°]
		double deg2 = atan2(y2, x2) * 180.0 / M_PI; // use vector atan2 to get rotational angle (phase) of both basis vectors [-180°,180°]

		if (abs(xx) - abs(yy) < 0) // check for corresponding order of singular values, correct if necessary..
		{
			signs.at(0) = yy < 0 ? true : false;
			signs.at(1) = xx < 0 ? true : false;
		}
		else
		{
			signs.at(0) = xx < 0 ? true : false;
			signs.at(1) = yy < 0 ? true : false;
		}

		signMap.at(i) = signs;
		// shift (normalize) degs from [-180°,180°] into the interval [0°,360°] - "circular value permutation"
		deg1 = deg1 < 0 ? 360 + deg1 : deg1;
		deg2 = deg2 < 0 ? 360 + deg2 : deg2;

		// singular values, decreasing order, corresponding singular vector order, scale ellipses axes in corresponding directions..
		double sv1 = svdList.at(i).singularValues()[0];
		double sv2 = svdList.at(i).singularValues()[1];
		double dot = sv2 * sv1;

		double sum = 0.0;
		if (sv1 == 0 || sv2 == 0 || sv1 / sv2 > 20.0)
		{
			glyphBuffer.at(i).at(round(deg1*steps / 360)) = sv1;
			glyphBuffer.at(i).at(static_cast<int>(round(deg1*steps / 360 + steps / 2)) % steps) = sv1;
			sum += 2 * sv1;
		}
		else
		{

			for (int j = 0; j < steps; j++) // sample ellipse equation for all steps
			{
				double val = dot / sqrt(sv2*sv2*cos(j*radres - deg1 * (M_PI / 180.0))*cos(j*radres - deg1 * (M_PI / 180.0)) + sv1 * sv1*sin(j*radres - deg1 * (M_PI / 180.0))*sin(j*radres - deg1 * (M_PI / 180.0))); //--> ellipse equation, evaluate for tMean (sum2)
				sum += val;
				glyphBuffer.at(i).at(j) = val;
			}
		}
		double rMean = sum / steps; // compute rMean from cartesian (rectangular) energy-based integral as opposed to the polar integral relevant to the geometrical (triangular/circular) area

		glyphBuffer.at(i) = 1.0 / rMean * glyphBuffer.at(i);
	}
}

std::vector<double> sample(double(*f)(double x), std::vector<std::vector<double>>& sampleArray, int jIndex, int iIndex)
{
	std::vector<double> sample(steps, 0.0);
	for (int i = 0; i < steps; i++)
		sample.at(i) = f(i*radres);
	return sample;
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
	
	double* meanA;
	double cosine_sum = 0.0;
	// create member vectors (arrays) for storing the sampled directions theta
	std::vector<std::vector<double>>* sampleBufferA;
	std::vector<std::vector<double>>* sampleBufferB;
	std::vector<std::vector<double>>* glyphBuffer;
	
	std::vector<std::vector<double>> cosines;
	
	// create sample vector (dynamic)
	std::vector<double> read;
	std::vector<double> glyph;
	std::vector<double> out;

	neighborhood hood; 
	bool flag = false;

	std::vector<double> initArray;
public:
	propagator(const int dim, double* mean, std::vector<std::vector<double>>* sampleBuffA, std::vector<std::vector<double>>* sampleBuffB, std::vector<std::vector<double>>* ellipseArray)
	{
		// assign ptrs to member vectors
		sampleBufferA = sampleBuffA;
		sampleBufferB = sampleBuffB;
		glyphBuffer = ellipseArray;
		meanA = mean;

		// initialize member samples w. 0
		read = std::vector<double>(steps, 0.0);
		glyph = std::vector<double>(steps, 0.0);
		out = std::vector<double>(steps, 0.0);
		initArray = std::vector<double>(steps, 0.0);

		double energy_sum = 0.0;
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
		cout.precision(dbl::max_digits10);
		cout << "cosine_sum: " << cosine_sum << endl;
	}

	void propagate()
	{
		// 1 propagation cycle
		for (int i = 0; i < width*height; i++) // for each node..
		{
			read = sampleBufferA->at(i);
			if (read == initArray)
				continue;
			glyph = glyphBuffer->at(i);

			flag = false;
			if (i / width == 0 || i % width == 0 || i / width == height - 1 || i % width == width-1)
				{hood.change(i / width, i%width); flag = true;}
			
			// calculate mean and variance.. of I(phi)
			double sum1 = 0.0;
			double sum2 = 0.0;
			double sum3 = 0.0;

			for (int j = 0; j < steps; j++)
			{
				sum1 += read.at(j); // evaluate for iMean (sum1)
				sum2 += glyphBuffer->at(i).at(j);
				sum3 += read.at(j)*glyphBuffer->at(i).at(j);

				if (read.at(j) < 0)
					cout << "WARNING: Negative values in SAMPLE in grid encountered, UNEXPECTED Behavior can not be excluded (negative energies->mirroring by 180 deg in polar plot) !!!" << endl;
				if (glyphBuffer->at(i).at(j) < 0)
					cout << "WARNING: Negative values in GLYPHS (Tensor Ellipse Eq.) encountered, UNEXPECTED Behavior can not be excluded (negative energies->mirroring by 180 deg in polar plot) !!!" << endl;
			}

			// compute iMean from cartesian (rectangular) energy-based integral as opposed to the polar integral relevant to the geometrical (triangular/circular) area
			double iMean = sum1 / steps; // -->tinc(dt) is a constant that can be drawn out of the integral
			// compute mean(T) from cartesian (rectangular) energy-based integral as opposed to the polar integral relevant to the geometrical (triangular/circular) area	
			double tMean = sum2 / steps; // -->tinc(dt) is a constant that can be drawn out of the integral
			// compute mean(T*I) from cartesian (rectangular) energy-based integral as opposed to the polar integral relevant to the geometrical (triangular/circular) area
			double tiMean = sum3 / steps; // -->tinc(dt) is a constant that can be drawn out of the integral
			// compute correction factor (scaling to mean=1, subsequent scaling to mean(I)), which follows energy conservation principles
			double cFactor = 1.0;
			if(tiMean)
				cFactor = tMean * iMean / tiMean;

			// iterate through central directions array to distribute (spread) energy (intensity) to the cell neighbors
			for (int k = 0; k < 8; k++) // for each adjacent edge...
			{
				// empty (reset) sample, upper and lower for each edge
				std::fill(out.begin(), out.end(), 0);

				if (flag) // if position on grid borders..
				{
					// check the neighborhood for missing (or already processed) neighbors, if missing, skip step..continue
					if (k == 0 && !hood.getR())
						continue;
					if (k == 1 && (!hood.getT() || !hood.getR()))
						continue;
					if (k == 2 && !hood.getT())
						continue;
					if (k == 3 && (!hood.getT() || !hood.getL()))
						continue;
					if (k == 4 && !hood.getL())
						continue;
					if (k == 5 && (!hood.getB() || !hood.getL()))
						continue;
					if (k == 6 && !hood.getB())
						continue;
					if (k == 7 && (!hood.getB() || !hood.getR()))
						continue;
				}

				int midIndex = (k * pi / 4) / radres;
				int index = betaIndex;
				if (k % 2 == 0)
				{
					if (total_anisotropy)
						index = centralIndex;
					else
						index = shiftIndex / 2;
				}

				//double energy_sum = 0.0;
				double val_sum = 0.0;

				for (int j = midIndex - index; j <= midIndex + index; j++) // for each step (along edge)..
				{
					int deltaJ = j - midIndex;
					int j_index = j < 0 ? j+steps : j%steps; // cyclic value permutation in case i exceeds the full circle degree 2pi

					double val = cFactor * read.at(j_index)*glyph.at(j_index);

					if (!total_anisotropy)
					{
						// split overlapping diagonal cones w.r.t to their relative angular area (obtained from face neighbors)..
						if ((abs(deltaJ) > centralIndex) && k % 2 == 0) // for alphas, use edge overlap > centralIndex
							if (abs(deltaJ) == shiftIndex / 2)
								val = 0.5*0.3622909908722584*val;
							else
								val = 0.3622909908722584*val;
						else if (k % 2 != 0) // for betas (diagonals), use static edge overlap-
							val = 0.6377090091277417*val;
					}
					
					
					val_sum += val * radres;
					//for (int l = j - shiftIndex; l < j + shiftIndex; l++) // for each step (along edge)..
					//{
					//	int deltaL = l - midIndex;
					//	if (abs(deltaL) >= shiftIndex)
					//		continue;
					//	int l_index = l % steps;
					//	if (l < 0)
					//		l_index = l + steps; // cyclic (circular) value permutation in case i exceeds the full circle degree 2pi
					//	double res = val * clip(cos((j_index - l_index) * radres), 0.0, 1.0);// *clip(cos(round(offset / radres) - dirIndex * pi / 2), 0.0, 1.0); // integrate over angle in polar
					//	out.at(l_index) += res;// *clip(cos(round(offset / radres) - dirIndex * pi / 2), 0.0, 1.0); // integrate over angle in polar
					//	energy_sum += res * radres;
					//}
				}

				out = (val_sum / cosine_sum) * cosines.at(k);

				*meanA += val_sum;

				switch (k) // propagate correspondent to each edge dir w.r.t forward edges
				{
				case 0:	sampleBufferB->at(i + 1) = sampleBufferB->at(i + 1) + out; break;
				case 1:	sampleBufferB->at(i + 1 - width) = sampleBufferB->at(i + 1 - width) + out; break;
				case 2:	sampleBufferB->at(i - width) = sampleBufferB->at(i - width) + out; break;
				case 3: sampleBufferB->at(i - 1 - width) = sampleBufferB->at(i - 1 - width) + out; break;
				case 4: sampleBufferB->at(i - 1) = sampleBufferB->at(i - 1) + out; break;
				case 5: sampleBufferB->at(i - 1 + width) = sampleBufferB->at(i - 1 + width) + out; break;
				case 6: sampleBufferB->at(i + width) = sampleBufferB->at(i + width) + out; break;
				case 7: sampleBufferB->at(i + 1 + width) = sampleBufferB->at(i + 1 + width) + out; break;
		
				default: break;
				}
			}
		
		}
	}
};

//template <typename T>
double acc(std::vector<double>& vec)
{
	double sum = 0.0;
	for (int i = 0; i < vec.size(); i++)
		sum += abs(vec.at(i));
	return sum;
	/*double sum = std::accumulate(m.begin(), m.end(), 0, [](auto lhs, const auto& rhs)
	{
		return std::accumulate(rhs.begin(), rhs.end(), lhs);
	});
	return sum;*/
}

//template <typename T>
double acc2(std::vector<std::vector<double>>& vec)
{
	double sum = 0.0;
	for (int i = 0; i < vec.size(); i++)
		sum += acc(vec.at(i));
	
	return sum;
}

template<typename Iter_T>
double vectorNorm(Iter_T first, Iter_T last)
{
	return sqrt(inner_product(first, last, first, 0.0));
}

int main(int argc, char* argv[])
{
	// 2D Grid START //
	int cols = 0; // create cols of txt tensor field
	int rows = 0; // create rows of txt tensor field
	cout << "before matrix read" << endl;
	workDir = GetCurrentWorkingDir();
	MatrixXd m = readMatrix(workDir + "/matrix.txt", &cols, &rows); // call countMatrix to determine rows/cols count #
	width = cols / 2; // determine width of grid for correct indexing
	height = rows / 2;
	const int dim = width * height; // determine # of dimensions of grid for buffer (string/coefficient etc..) vectors
	lightSrcPos = { height / 2, width / 2 }; // initialize light src position option w. center point

	// create vector for light src's symbolic user-input functions in convenient string format
	std::vector<std::string> userFunctions;
	std::vector<Pair> userPositions;
	// parse input option file
	parse_options(argc, argv, userFunctions, userPositions);

	const std::vector<double> initArray(steps, 0.0);

	// define dual buffers for propagation
	std::vector<std::vector<double>> sampleBufferA((width*height), initArray);
	std::vector<std::vector<double>> sampleBufferB((width*height), initArray);
	std::vector<std::vector<double>> sampleBufferMem((width*height), initArray);
	std::vector<std::vector<double>> sampleBufferInit((width*height), initArray);
	std::vector<std::vector<double>> glyphBuffer((width*height), initArray);

	// create distribution buffer to gather light distributions for all positions y,x and directions t
	std::vector<std::vector<std::vector<double>>> distBuffer(width*height*steps, glyphBuffer); // initialize #steps 2D-planes w. empty glyphBuffer
	
	// std::vector<std::vector<double>> dist = distBuffer.at(0) - distBuffer.at(0); // Exemplary OPERATOR USE TODO!

	std::vector<double> initGradient(3, 0.0); // construct 3D Gradient
	std::vector<std::vector<double>> deltaBuffer(width*height*steps, initGradient); // initialize #steps 2D-planes w. empty glyphBuffer

	std::vector<std::vector<double>> lightSrcs;
	cout << "before compute glyphs" << endl;
	
	std::vector<std::vector<bool>> signMap((width)*(height), std::vector<bool>(2, false)); // create a signMap relating normal force signs to singular values
	// compute Eigenframes/Superquadrics/Ellipses/Glyphs by calling computeGlyphs w. respective args
	computeGlyphs(glyphBuffer, signMap);
	
	// get defined light srcs positions and intensities...
	//if(userFunctions.size()) // if entered by user, use cmd AND config
	//	for (int i = 0; i < userFunctions.size(); i++)
	//	{
	//		functionString = userFunctions.at(i);
	//		lightSrcPos = userPositions.at(i);
	//		lightSrcs.push_back(sample(strFunction, sampleBufferA, lightSrcPos.jIndex, lightSrcPos.iIndex)); // sample the light profile w. muParser
	//	}
	//else // if entered by default value, use default value in center position
	//{
	//	functionString = std::to_string(intensity);
	//	lightSrcs.push_back(sample(strFunction, sampleBufferA, lightSrcPos.jIndex, lightSrcPos.iIndex)); // sample the light profile w. muParser
	//	userPositions.push_back(Pair(lightSrcPos.jIndex, lightSrcPos.iIndex));
	//}
	//for (int i = 0; i < lightSrcs.size(); i++) // enter the light src's profiles into the grid positions in userPositions
	//	sampleBufferA.at(userPositions.at(i).jIndex*width + userPositions.at(i).iIndex) = lightSrcs.at(i); // initialize grid (sampleBufferA) w. "light src list" lightSrcs

	double meanA = 0.0; // set up mean variables for threshold comparison as STOP criterion..
	double meanMem = 0.0;
	// create propagator object (managing propagation, reprojection, correction, central directions, apertureAngles and more...)
	propagator prop(dim, &meanA, &sampleBufferA, &sampleBufferB, &glyphBuffer);

	std::vector<double> sample = initArray;
	/*std::ofstream file("cout.txt");
	Tee tee(cout, file);
	TeeStream both(tee);
	both.precision(dbl::max_digits10);*/

	// construct a light src vector (delta functions for each sampled direction - normalized to total area (energy) of unit circle 1.0)
	for (int j = 0; j < steps; j++)
	{
		sample.at(j) = steps;
		lightSrcs.push_back(sample);
		sample = initArray;
	}

	int ctr = 0;
	int tCtr = 0;
	double duration;
	double total = 0.0;
	bool finished = false;
	std::clock_t startTotal = std::clock();
	
	std::vector<double> lightSrc;
	for (int t = 0; t < steps; t++)
	{
		lightSrc = lightSrcs.at(t);
		cout << "before propagating t: " << t << endl;
		std::clock_t start = std::clock();
		for (int j = 0; j < height; j++)
			for (int i = 0; i < width; i++)
			{
				//tCtr = 0;
				//lightSrcPos = Pair(j, i);

				// DUAL BUFFER PROPAGATION //
				sampleBufferA.at(j*width + i) = lightSrc; // get pre-computed light src for current direction t
				meanMem = 0.0;
				//ctr = 0;
				finished = false;

				// loop over nodes in grid and propagate until error to previous light distribution minimal <thresh
				while (!finished) // perform one single light propagation pass (iteration)
				{
					meanA = 0.0;
					prop.propagate(); // propagate until finished..
					//meanA *= (1.0 / radres) / (steps*sampleBufferA.size());
					swap(sampleBufferA, sampleBufferB);
					//sampleBufferA = sampleBufferB;
					sampleBufferA.at(j*width + i) = lightSrc;

					if (abs(meanA - meanMem) < thresh)
						finished = true;
					meanMem = meanA;

					//if (ctr == ctrLimit)//6
					//	break;

					//ctr++;
					sampleBufferB = sampleBufferInit;
					//std::fill(sampleBufferB.begin(), sampleBufferB.end(), initArray);

				}

				/*if (ctr <= 3)
				tCtr++;*/
				sampleBufferA.at(j*width + i) = initArray; //remove light src to prevent trivial differences at light src positions 
				distBuffer.at(j*width + i + t * dim) = sampleBufferA;
				//std::fill(sampleBufferA.begin(), sampleBufferA.end(), initArray);
				sampleBufferA = sampleBufferInit;
			}
		duration = ((std::clock() - start)*1000.0 / (double)CLOCKS_PER_SEC);
		cout << "timer: " << duration << " ms" << endl;
		cout << "..after propagation, (last) ctr:" << ctr << endl;
		cout << "..trivial ctr:" << tCtr << endl;
	}
	total = ((std::clock() - startTotal)*1000.0 / (double)CLOCKS_PER_SEC);
	cout << "..after propagation TOTAL, total timer:" << total << " ms" << endl;
	// PROPAGATION SCHEME END //
	sampleBufferB.clear();
	sampleBufferMem.clear();
	sampleBufferInit.clear();
	glyphBuffer.clear();
	
	sampleBufferB.shrink_to_fit();
	sampleBufferMem.shrink_to_fit();
	sampleBufferInit.shrink_to_fit();
	glyphBuffer.shrink_to_fit();
	// DELTA (Gradient) COMPUTATION START //
	std::vector<double> gradient(3, 0.0); // dim3: x,y,theta

	cout << "before constructing gradient vector.." << endl;
	startTotal = std::clock();
	for (int t = 0; t < steps; t++)
		for (int j = 0; j < height; j++)
			for (int i = 0; i < width; i++)
			{
				if (i == 0 || i == width - 1 || j == 0 || j == height - 1)
					continue;
				// X-1D central differences.. VARIANT 1: bin by bin
				sampleBufferA = distBuffer.at(j*width + i + 1 + t * dim) - distBuffer.at(j*width + i - 1 + t * dim);
				meanA = acc2(sampleBufferA);
				gradient.at(0) = meanA /2.0;
				// Y-1D central differences..
				sampleBufferA = distBuffer.at((j+1)*width + i + t * dim) - distBuffer.at((j-1)*width + i + t * dim);
				meanA = acc2(sampleBufferA);
				gradient.at(1) = meanA / 2.0;
				// t-1D central differences..
				sampleBufferA = distBuffer.at(j*width + i + (t+1)%steps * dim) - distBuffer.at(j *width + i + (t == 0 ? (steps-1) : (t-1)) * dim);
				meanA = acc2(sampleBufferA);
				gradient.at(2) = meanA / 2.0;
				
				// X-1D central differences.. VARIANT 2: cell by cell
				//meanA = 0.0;
				//for(int k = 0; k < dim; k++)
				//	meanA += abs(std::accumulate(distBuffer.at(j*width + i + 1 + t * dim).at(k).begin(), distBuffer.at(j*width + i + 1 + t * dim).at(k).end(),0.0) - std::accumulate(distBuffer.at(j*width + i - 1 + t * dim).at(k).begin(), distBuffer.at(j*width + i - 1 + t * dim).at(k).end(),0.0));
				//
				//gradient.at(0) = meanA / 2.0;
				//// Y-1D central differences..
				//meanA = 0.0;
				//for (int k = 0; k < dim; k++)
				//	meanA += abs(std::accumulate(distBuffer.at((j + 1)*width + i + t * dim).at(k).begin(), distBuffer.at((j + 1)*width + i + t * dim).at(k).end(),0.0) - std::accumulate(distBuffer.at((j - 1)*width + i + t * dim).at(k).begin(), distBuffer.at((j - 1)*width + i + t * dim).at(k).end(),0.0));
				//
				//gradient.at(1) = meanA / 2.0;

				//// t-1D central differences..
				//meanA = 0.0;
				//for (int k = 0; k < dim; k++)
				//	meanA += abs(std::accumulate(distBuffer.at(j*width + i + (t + 1) % steps * dim).at(k).begin(), distBuffer.at(j*width + i + (t + 1) % steps * dim).at(k).end(),0.0) - std::accumulate(distBuffer.at(j *width + i + (t == 0 ? (steps - 1) : (t - 1)) * dim).at(k).begin(), distBuffer.at(j *width + i + (t == 0 ? (steps - 1) : (t - 1)) * dim).at(k).end(),0.0));

				//gradient.at(2) = meanA / 2.0;

				deltaBuffer.at(j*width + i + t * dim) = gradient;
			}

	// DELTA (Gradient) COMPUTATION END //
	duration = ((std::clock() - startTotal)*1000.0 / (double)CLOCKS_PER_SEC);
	cout << "..after, timer: " << duration << " ms" << endl;
	
	sampleBufferA.clear();
	distBuffer.clear();
	
	// vector norm (Gradient) COMPUTATION START //
	std::vector<double> scalarNorm(width*height*steps, 0.0); // construct 3D Gradient Norm Vector

	vtkSmartPointer<vtkImageData> imageData = vtkSmartPointer<vtkImageData>::New();

	imageData->SetDimensions(width, height, steps);

	//vtkSmartPointer<vtkDoubleArray> director =	vtkSmartPointer<vtkDoubleArray>::New();

	//director->SetNumberOfComponents(3);
	//director->SetNumberOfTuples(width * height * steps);

	vtkSmartPointer<vtkDoubleArray> energy = vtkSmartPointer<vtkDoubleArray>::New();

	energy->SetNumberOfComponents(1);
	energy->SetNumberOfTuples(width * height * steps);

	cout << "before computing gradient (vector) norm.." << endl;
	startTotal = std::clock();
	ctr = 0;
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