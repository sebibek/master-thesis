// INLINE DEFINES

#ifdef WINDOWS
#include <direct.h>
#define GetCurrentDir _getcwd
#else
#include <unistd.h>
#define GetCurrentDir getcwd
#endif

#define MAXBUFSIZE  ((int) 1e3)
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
#include "exprtk.hpp"
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <array>
#include <filesystem>
#include <string>
//#include <algorithm>
//#include <functional>
//#include <cmath>

// NAMESPACE IMPORTS
using namespace std;
using namespace Eigen;
using namespace mu;

// PROTOTYPES
std::string functionString; // Prototyping of functionString for strFunction (symbolic string arg parsed by muparser in cpp functional ptr rep)!
int width;
int height;
template <typename T>
T clip(const T& n, const T& lower, const T& upper); // template clip function

//length and width of window
int wSize = 700;
sf::Vector2i windowsize;

//definition of pi
const double pi = M_PI;

// GENERIC FUNCTION DEFINITIONS

// -> muparser custom clip function: template functions need to be static callback (call after) functions, which can be passed as arguments
static value_type clipNeg(value_type v) { return clip(v, 0.0, 1.0); }

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

void defineConvexEllipse(sf::ConvexShape* ellipse, double radius_x, double radius_y, unsigned short quality)
{
	for (int i = 0; i < quality; ++i)
	{
		double rad = (360 / quality * i) / (360 / M_PI / 2);
		double x = cos(rad)*radius_x;
		double y = sin(rad)*radius_y;

		ellipse[0].setPoint(i, sf::Vector2f(x, y));
	}
}

// CLASS DEFINITIONS //

class neighborhood
{
	bool neighborR = true;
	bool neighborT = true;
	bool neighborL = true;
	bool neighborB = true;

public:

	//constructor
	neighborhood(int j, int i, const int length) // check the neighbourhood of (j,i) for missing neighbors..
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

	bool getR() { return neighborR; }
	bool getT() { return neighborT; }
	bool getL() { return neighborL; }
	bool getB() { return neighborB; }
};
class polarPlot
{
	// muparser object for parsing mathemathical function
	mu::Parser p;
	//initialize graphics stuff
	sf::Image graph;
	sf::Texture graphTx;
	sf::Sprite graphSpr;

	// draw speed
	int speed;
	int line_width;
	// theta
	double t;
	// radius
	double r;
	double tmax;
	double tinc;
	double rotation;
	double radres;
	sf::Color color;
	std::vector<double> sample;
	bool discretePlot = false;

public:

	//constructor #1 string parser - muParser
	polarPlot(std::string tag,
		int speed_in,
		int line_width_in,
		double tmax_in,
		double tinc_in,
		double rotation_in,
		sf::Color color_in)

		: speed(speed_in),
		line_width(line_width_in),
		t(0),
		tmax(tmax_in),
		tinc(tinc_in),
		rotation(rotation_in),
		color(color_in)
	{
		p.DefineVar("theta", &t);
		p.SetExpr(tag);
	}

	//constructor #2 array parser
	polarPlot(std::vector<double> Sample,
		double res,
		int speed_in,
		int line_width_in,
		double tmax_in,
		double tinc_in,
		double rotation_in,
		sf::Color color_in)

		: speed(speed_in),
		line_width(line_width_in),
		t(0),
		tmax(tmax_in),
		tinc(tinc_in),
		rotation(rotation_in),
		sample(Sample),
		radres(res),
		color(color_in)
	{
		discretePlot = true;
	}

	//initialize muparser and graphics 
	void init() {
		p.DefineVar("theta", &t);
		p.DefineConst("pi", pi);
		graph.create(wSize, wSize, sf::Color::Transparent);
		graphTx.loadFromImage(graph);
		graphSpr.setTexture(graphTx);

		//center on screen
		graphSpr.setOrigin(wSize / 2, wSize / 2);
		graphSpr.setPosition(windowsize.x / 2, windowsize.y / 2);
	}

	//iterate function draw : TODO
	sf::Sprite& update() {
		graphTx.update(graph);
		graphSpr.rotate(rotation * 180.0 / pi);
		return graphSpr;
	}
	//plots point to pixel(s)
	void plot(int index = 0, int mode = 0) {

		//scale
		double scale;
		if (mode == 1)
			scale = 20;
		else
			scale = 50;
		//convert polar to cartesian
		double x = scale * r * cos(t);
		double y = -1 * scale * r * sin(t);

		if (x + 1 + wSize / 2 > wSize || x - 1 + wSize / 2 < 0 || y + 1 + wSize / 2 > wSize || y - 1 + wSize / 2 < 0)
			return;

		int offsetX = (index % width)*(wSize / width) + (wSize / (2 * width)); // check rest (modulo) for x-offset
		int offsetY = (index / width)*(wSize / width) + (wSize / (2 * width)); // check division for y offset

		int xIndex = static_cast<int>(std::round(x + offsetX));
		int yIndex = static_cast<int>(std::round(y + offsetY));

		if (xIndex < 0)
			xIndex = 0;
		if (xIndex >= wSize)
			xIndex = wSize - 1;
		if (yIndex < 0)
			yIndex = 0;
		if (yIndex >= wSize)
			yIndex = wSize - 1;

		graph.setPixel((xIndex), (yIndex), color);
		//graph.setPixel((x + wSize / 2), (y + wSize / 2), color);

		// set 4 -neighborhood
		if (line_width == 1 || line_width == 2) {
			graph.setPixel((xIndex), (yIndex + 1), color);
			graph.setPixel((xIndex), (yIndex - 1), color);
			graph.setPixel((xIndex + 1), (yIndex), color);
			graph.setPixel((xIndex - 1), (yIndex), color);
		}

		// set diagonal neighborhood
		if (line_width == 2) {
			graph.setPixel((xIndex + 1), (yIndex + 1), color);
			graph.setPixel((xIndex + 1), (yIndex - 1), color);
			graph.setPixel((xIndex - 1), (yIndex - 1), color);
			graph.setPixel((xIndex - 1), (yIndex + 1), color);
		}
	}

	//plots draw_speed points per frame
	void animation(int index = 0, int mode = 0)
	{


		for (int i = 0; (i < speed || speed == 0) && t < tmax; i++)
		{
			if (!discretePlot) // use muParser for evaluation 
				r = p.Eval();
			else // nearest neighbor interpolation w. samples in sampleVector - ONLY FOR DUAL GRID PLOT! - use animation overload for cell-based plot
				r = sample.at(static_cast<int>(std::round(t / radres)) % sample.size()); // cyclic value permutation

			// PROPAGATION ATTENUATION TEST //
			/*if (t == tinc && index/width == 3 && mode == 0)
			{
				cout << "r: " << r << endl;
				cout << "attenuation: " << r / 1.0 << endl;
			}*/
			plot(index, mode);
			t += tinc;
		}
	}
	//plots draw_speed points per frame
	void animation(int index, int mode, std::vector<std::vector<double>>& sampleArray)
	{
		for (int i = 0; (i < speed || speed == 0) && t < tmax; i++)
		{
			if (!discretePlot) // use muParser for evaluation 
				r = p.Eval();
			else // nearest neighbor interpolation w. samples in sampleVector - ONLY FOR DUAL GRID PLOT! - use animation overload for cell-based plot
				r = sampleArray.at(index).at(static_cast<int>(std::round(t / radres)) % sample.size());

			// PROPAGATION ATTENUATION TEST //
			if (t == tinc && index / width == width/2 && mode == 0)
			{
				cout << "r: " << r << endl;
				cout << "attenuation: " << r / 1.0 << endl;
			}

			if(t > 45*tinc && t < 47*tinc && index/width == width/2 - 1 && index%width == width/2 + 1 )
			{
				cout << "r(diag): " << r << endl;
				cout << "attenuation: " << r / 1.0 << endl;
			}
			plot(index, mode);
			t += tinc;
		}
	}
};

class myPair
{
public:
	int jIndex = width / 2;
	int iIndex = width / 2;

	myPair() {}
	myPair(int j, int i)
	{
		jIndex = j;
		iIndex = i;
	}

	friend istringstream& operator>>(istringstream& stream, myPair& pair)
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
myPair lightSrcPos{ width / 2, width / 2 };
bool parallel = false;
bool backprop = false;
bool reverse_interpolation = false;
bool fullscreen = false; //fullscreen flag
std::string record_folder = "frames";//directory to write images to, must exist
int record_frameskip = 0; // --> disables recording //recording frameskip 

void parse_options(int argc, char* argv[]) {

	int c;
	std::string frameskip_opt = "-1";
	std::string size_opt = "-1";
	std::string directory_opt = "";
	std::string light_opt = "-1";
	std::string parallel_opt = "-1";
	std::string backprop_opt = "-1";
	std::string reverse_interpolation_opt = "-1";
	//	extern char *optarg;
	//	extern int optind, optopt;

	//using getopts 
	while ((c = getopt(argc, argv, "fs:a:d:l:p:b:r:")) != -1) // --> classic getopt call: argc: # of args; argv: array of args(strings); "": optstring (registered options, followed by colon when taking args itself) 
	{
		int f = -1, s = -1;

		switch (c)
		{
			// correct option use
		case 'l': {
			light_opt.assign(optarg);
			std::istringstream(light_opt) >> lightSrcPos;
			break;
		}
		case 'p': {
			parallel_opt.assign(optarg);
			std::istringstream(parallel_opt) >> parallel;
			break;
		}
		case 'b': {
			backprop_opt.assign(optarg);
			std::istringstream(backprop_opt) >> backprop;
			break;
		}
		case 'r': {
			reverse_interpolation_opt.assign(optarg);
			std::istringstream(reverse_interpolation_opt) >> reverse_interpolation;
			break;
		}
		case 'f': {
			fullscreen = true;
			break; }
		case 'a':
			frameskip_opt.assign(optarg);
			std::istringstream(frameskip_opt) >> f;
			if (f <= -1) {
				std::cerr << "invalid argument \'" << frameskip_opt << "\' to option -r\n";
				optind--;
			}
			record_frameskip = f > 0 ? f : 0;
			break;
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
}

// OPERATOR DEFINITIONS //
template <typename T> // element-wise plus for std::vector
std::array<T, 21> operator*(const T a, const std::array<T, 21>& b)
{
	std::array<T, 21> result;

	for (int i = 0; i < b.size(); i++)
		result.at(i) = a * b.at(i);

	return result;
}

template <typename T> // element-wise plus for std::array
std::array<std::array<T, 21>, 2> operator*(const T a, const std::array<std::array<T, 21>, 2> &b)
{
	std::array<std::array<T, 21>, 2> result;

	for (int i = 0; i < b.size(); i++)
		result.at(i) = a * b.at(i);

	return result;
}

template <typename T> // element-wise plus for std::vector
std::array<T, 21> operator+(const std::array<T, 21>& a, const std::array<T, 21>& b)
{
	assert(a.size() == b.size());

	std::array<T, 21> result;

	for (int i = 0; i < a.size(); i++)
		result.at(i) = a.at(i) + b.at(i);

	return result;
}

template <typename T> // element-wise plus for std::array
std::array<std::array<T, 21>, 2> operator+(const std::array<std::array<T, 21>, 2>& a, const std::array<std::array<T, 21>, 2> &b)
{
	assert(a.size() == b.size());

	std::array<std::array<T, 21>, 2> result;

	for (int i = 0; i < a.size(); i++)
		result.at(i) = a.at(i) + b.at(i);

	return result;
}

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

// SPECIFIC (MATHEMATICAL-PROGRAM USE) FUNCTIONS

// template function for evaluating string functions via ptr 
template <typename T>
T strFunction(T theta)
{
	typedef exprtk::symbol_table<T> symbol_table_t;
	typedef exprtk::expression<T>     expression_t;
	typedef exprtk::parser<T>             parser_t;

	std::string expression_string = functionString;

	symbol_table_t symbol_table;
	symbol_table.add_variable("theta", theta);
	symbol_table.add_constants();
	symbol_table.add_function(_T("clip"), clipNeg);

	expression_t expression;
	expression.register_symbol_table(symbol_table);

	parser_t parser;
	//parser.DefineFun(_T("clip"), clipNeg, false);
	parser.compile(expression_string, expression);
	T y = 0;

	y = expression.value();

	return y;
}

template <typename T>
T clip(const T& n, const T& lower, const T& upper) // template clip function
{
	return std::max(lower, std::min(n, upper));
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
		stringstream stream(line);
		while (!stream.eof())
			stream >> buff[cols*rows + temp_cols++];

		if (temp_cols == 0)
			continue;

		if (cols == 0)
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

void computeGlyphs(std::vector<std::string>& functionStrEllipse, std::vector<std::tuple<double, double, double>>& ellipseArray)
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
		}

	// compute the SVD (singular value decomposition)of all matrices in matrixList into svdList
	std::vector<JacobiSVD<MatrixXd>> svdList(dim, JacobiSVD<MatrixXd>(matrixList.at(0), ComputeThinU | ComputeThinV));
	for (int i = 0; i < dim; i++)
	{
		// SVD-based
		JacobiSVD<MatrixXd> svd(matrixList.at(i), ComputeThinU | ComputeThinV);
		svdList.at(i) = svd;
	}

	// define parser parameters
	double t = 0; // ->theta: set initial theta
	double tinc = 2 * pi / 72; // set theta increment tinc
	int steps = (2 * pi) / tinc; // set # of steps

	// define function parser (muparser)
	mu::Parser parser;
	parser.DefineConst("pi", pi);
	parser.DefineVar("theta", &t);

	// iterate through the matrixList/svdList and construct (scaled) ellipses in polar form (function) from the repsective singular values/vectors
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

		// shift (normalize) degs from [-180°,180°] into the interval [0°,360°] - "circular value permutation"
		if (deg1 < 0)
			deg1 = 360 + deg1;
		if (deg2 < 0)
			deg2 = 360 + deg2;

		double deg = 0;
		if (deg1 >= 270 && deg2 <= 90) // caveat: if basis vector u1 in 4th quadrant and u2 in 1st - continuity/consistency constraint
			deg = atan(y1 / x1) * 180.0 / M_PI; // use u1 as lagging vector
		else if (deg2 >= 270 && deg1 <= 90) // caveat: if basis vector u2 in 4th quadrant and u2 in 1st - continuity/consistency constraint
			deg = atan(y2 / x2) * 180.0 / M_PI;
		else if (deg2 > deg1) // if u2 is leading vector 
			deg = atan(y1 / x1) * 180.0 / M_PI; // use u1 as lagging vector
		else // if u1 is leading vector
			deg = atan(y2 / x2) * 180.0 / M_PI; // use u2 as lagging vector
		// --> u1 and u2 form basis vectors obtained from SVD, whereas the phase(leading vector) - phase(lagging vector) = 90 --> determine u1,u2

		double sv1 = svdList.at(i).singularValues()[0];
		double sv2 = svdList.at(i).singularValues()[1];

		double dot = sv2 * sv1;
		double a = 0; // x-scale
		double b = 0; // y-scale

		// --> eventually, use std::swap to swap ptrs to sv1,sv2 dynamically
		if (((yy < 0 && xx < 0) || (yy >= 0 && xx >= 0))) // if not mirroring at one single axis..
		{
			if (abs(yy) < abs(xx)) // if anisotropic x-scaling, scale x
			{
				a = sv1;
				b = sv2;
			}
			else // if anisotropic y-scaling, scale y
			{
				a = sv2;
				b = sv1;
			}
		}
		else // if mirroring at one single axis.. swap!
		{
			if (abs(yy) < abs(xx)) // if anisotropic x-scaling, scale x
			{
				a = sv2;
				b = sv1;
			}
			else // if anisotropic y-scaling, scale y
			{
				a = sv1;
				b = sv2;
			}
		}

		// assign ellipse parameters in ellipseArray
		ellipseArray.at(i) = { sv1, sv2, deg*(M_PI / 180.0) };

		// define products as string modules to form (construct) ellipse string for current transmission profile T(w)
		std::string aSquaredString = std::to_string(a*a);
		std::string bSquaredString = std::to_string(b*b);
		std::string dotString = std::to_string(dot);
		std::string radString = std::to_string(deg*(M_PI / 180.0));

		std::string ellipse = dotString + "/sqrt(" + bSquaredString + "*cos(theta-" + radString + ")*cos(theta-" + radString + ")+" + aSquaredString + "*sin(theta-" + radString + ")*sin(theta-" + radString + "))";

		// define parser parameters
		t = 0; // ->theta: set initial theta
		parser.SetExpr(ellipse);

		double sum = 0.0;
		// sum r(theta) over (around) the ellipse to obtain mean(r)
		for (int i = 0; i < steps; i++)
		{
			sum += parser.Eval(); t += tinc;
		}

		double rMean = sum / steps; // compute rMean from cartesian (rectangular) energy-based integral as opposed to the polar integral relevant to the geometrical (triangular/circular) area

		if (rMean != 0) // only w. rMean != 0, else neglect normalization to prevent division by zero (NANs)...
		{
			// scale string-modules of ellipse string to rMean(Max)
			aSquaredString = std::to_string(a*a*(1 / rMean)*(1 / rMean));
			bSquaredString = std::to_string(b*b*(1 / rMean)*(1 / rMean));
			dotString = std::to_string(dot*(1 / rMean)*(1 / rMean));
		}

		ellipse = dotString + "/sqrt(" + bSquaredString + "*cos(theta-" + radString + ")*cos(theta-" + radString + ")+" + aSquaredString + "*sin(theta-" + radString + ")*sin(theta-" + radString + "))"; // insert normalization factor right into polar ellipse function at pos. 0

		// TODO: use ellipses normalized to mean(e)=1 for testing!, use ellipses normalized to mean(e)=max(mean(e)) for real fields!
		functionStrEllipse.at(i) = ellipse; // transmission profile T(w)
	}
}

void sample(double(*f)(double x), std::vector<std::vector<double>>& sampleArray, double radres, int steps, int jIndex, int iIndex)
{
	std::vector<double> sample(steps, 0.0);
	for (int i = 0; i < steps; i++)
		sample.at(i) = f(i*radres);
	sampleArray.at(iIndex + jIndex * (width)) = sample;
}

class propagator
{
	// set principal arture angles in rads
	double alpha = 36.8699 * M_PI / 180;
	double beta = 26.5651 * M_PI / 180;
	// define principal central directions array (12 central directions in 2D --> 30 in 3D, ufff..)
	std::array<double, 12> centralDirections{ 3.0 / 2 * M_PI + 1.0172232666228471, 0, 0.5535748055013016, 1.0172232666228471, M_PI / 2, M_PI / 2 + 0.5535748055013016, M_PI / 2 + 1.0172232666228471, M_PI, M_PI + 0.5535748055013016, M_PI + 1.0172232666228471, 3.0 / 2 * M_PI, 3.0 / 2 * M_PI + 0.5535748055013016 };
	std::array<double, 12> apertureAngles{ beta, alpha, beta, beta, alpha, beta, beta, alpha, beta, beta, alpha, beta }; // define aperture angles array in order

	// define function parser (muparser)
	double radres = 0.01745; // set radres for discretized theta (phi) directions for sampling
	int steps = (2 * pi) / radres;
	int stepsQuarter = steps / 4;

	// create member vectors (arrays) for storing the symbolic strings and fourier (CH) coefficients (external)
	std::vector<std::vector<double>>* sampleBufferA;
	std::vector<std::vector<double>>* sampleBufferB;
	std::vector<std::tuple<double, double, double>>* ellipseArray;
	
	// create sample vector (dynamic)
	std::vector<double> read;
	std::vector<double> out;

public:
	propagator(const int dim, std::vector<std::string>* fStr, std::vector<std::vector<double>>* sampleBuffA, std::vector<std::vector<double>>* sampleBuffB, std::vector<std::string>* fStrEllipse, std::vector<std::tuple<double, double, double>>* ellipseArr)
	{
		// assign ptrs to member vectors
		sampleBufferA = sampleBuffA;
		sampleBufferB = sampleBuffB;
		ellipseArray = ellipseArr;

		// initialize member samples w. 0
		read = std::vector<double>(steps, 0.0);
		out = std::vector<double>(steps, 0.0);
	}

	void propagate(bool parallel = false)
	{
		// 1 propagation cycle
		for (int i = 0; i < width*height; i++) // for each node..
		{
			neighborhood hood(i / width, i%width, width);
			
			std::fill(read.begin(), read.end(), 0);
			read = sampleBufferA->at(i);
		
			// iterate through central directions array to distribute (spread) energy (intensity) to the cell neighbors
			for (int k = 0; k < 8; k++) // for each adjacent edge...
			{
				// empty (reset) sample, upper and lower for each edge
				std::fill(out.begin(), out.end(), 0);

					// check the neighborhood for missing (or already processed) neighbors, if missing, skip step..continue
				if (!hood.getR() && k == 0)
					continue;
				if ((!hood.getT() || !hood.getR()) && k == 1)
					continue;
				if (!hood.getT() && k == 2)
					continue;
				if ((!hood.getT() || !hood.getL()) && k == 3)
					continue;
				if (!hood.getL() && k == 4)
					continue;
				if ((!hood.getB() || !hood.getL()) && k == 5)
					continue;
				if (!hood.getB() && k == 6)
					continue;
				if ((!hood.getB() || !hood.getR()) && k == 7)
					continue;

				int betaIndex = (beta)/radres;
				int shiftIndex = pi/2/radres;
				int midIndex = (k * pi / 4) / radres;
				int centralIndex = (alpha / 2) / radres;
				int index = betaIndex;
				if (k % 2 == 0)
					index =shiftIndex/2;

				double sum = 0.0;
				double valsum = 0.0;

				for (int j = midIndex - index; j < midIndex + index; j++) // for each step (along edge)..
				{
					int deltaJ = j - midIndex;
					int j_index = j%steps;

					if (j < 0)
						j_index = j + steps; // cyclic value permutation in case i exceeds the full circle degree 2pi
					double val = read.at(j_index);
					
					if ((abs(deltaJ) > centralIndex) && k % 2 == 0)
						val = 0.3712*val;
					else if (k%2 != 0)
						val = 0.6288*val;
					/*if ((abs(deltaJ) > centralIndex) || (k % 2 != 0))
						val = 0.5*read.at(j_index);*/
					
					valsum += val * radres;
					for (int l = j - shiftIndex; l < j + shiftIndex; l++) // for each step (along edge)..
					{
						int deltaL = l - midIndex;
						if (abs(deltaL) >= shiftIndex)
							continue;
						int l_index = l % steps;
						if (l < 0)
							l_index = l + steps; // cyclic value permutation in case i exceeds the full circle degree 2pi

						double res = val * clip(cos((j_index - l_index) * radres), 0.0, 1.0);// *clip(cos(round(offset / radres) - dirIndex * pi / 2), 0.0, 1.0); // integrate over angle in
						//if (k%2 != 0)
						//	res = val * pow(clip(cos((j_index - l_index) * radres), 0.0, 1.0),1.4);// *clip(cos(round(offset / radres) - dirIndex * pi / 2), 0.0, 1.0); // integrate over angle in

						out.at(l_index) += res;// *clip(cos(round(offset / radres) - dirIndex * pi / 2), 0.0, 1.0); // integrate over angle in
						sum += res * radres;
					}
				}

				if(sum> 0)
					out = (valsum / sum) * out;
						

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

int main(int argc, char* argv[])
{
	// 2D Grid START //
	int cols = 0; // create cols of txt tensor field
	int rows = 0; // create rows of txt tensor field
	cout << "before matrix read" << endl;
	std::string workDir = GetCurrentWorkingDir();
	MatrixXd m = readMatrix(workDir + "/matrix.txt", &cols, &rows); // call countMatrix to determine rows/cols count #
	const int dim = rows / 2 * cols / 2; // determine # of dimensions of grid for buffer (string/coefficient etc..) vectors
	width = cols / 2; // determine width of grid for correct indexing
	height = rows / 2;
	lightSrcPos = { width / 2, width / 2 }; // initialize light src position option w. center point

	// parse input option file
	parse_options(argc, argv);

	// create functionStr vectors to efficiently store the respective function strings for plotting..
	std::vector<std::string> functionStr(dim, "0.0"); // initialize w. 0.0 to prevent parser errors
	std::vector<std::string> functionStrEllipse(dim, "0.0"); // initialize w. 0.0 to prevent parser errors

	double tinc = 2 * pi / 72;
	double radres = 0.01745;
	int steps = 2 * pi / radres;

	// define dual buffers for propagation
	std::vector<double> initArray(steps, 0.0);
	std::vector<std::vector<double>> sampleBufferA((width)*(height), initArray);
	std::vector<std::vector<double>> sampleBufferB((width)*(height), initArray);
	std::vector<std::vector<double>> sampleBufferMem((width)*(height), initArray);
	//
	std::tuple<double, double, double> initTuple = { 0.0,0.0,0.0 }; // -->triple: lambda1, lambda2, deg
	std::vector<std::tuple<double,double,double>> ellipseArray(dim, initTuple); // ..use it to initialize the coefficient array w. dim elements

	//cout << "before compute glyphs" << endl;
	
	// compute Eigenframes/Superquadrics/Ellipses/Glyphs by calling computeGlyphs w. respective args
	computeGlyphs(functionStrEllipse, ellipseArray);
	
	// define excitation (stimulus) polar functions(4 neighbors/directions) normalized to area 1
	std::string circle = "2.1"; // circle w. 100% relative intensity - capturing the spread [0..1] normalized to strongest light src in field
	
	// write light src in coefficientArray
	functionString = circle; // set circular (isotroic) profile for light src
	int jIndex = lightSrcPos.jIndex; int iIndex = lightSrcPos.iIndex; // create light src position indices
	sample(strFunction, sampleBufferA, radres, steps, jIndex, iIndex); // sample the light profile w. muParser
	//parallel = true;
	
	// DUAL BUFFER PROPAGATION //
	// create propagator object (managing propagation, reprojection, correction, central directions, apertureAngles and more...)
	propagator prop(dim, &functionStr, &sampleBufferA, &sampleBufferB, &functionStrEllipse, &ellipseArray);

	double thresh = 0.000001;
	bool finished = false;
	// loop over nodes in grid and propagate until error to previous light distribution minimal <thresh
	int ctr = 0;
	for(int i = 0; i < 10; i++)
	//while (!finished)
	{
		prop.propagate(parallel); // propagate until finished..
		sampleBufferA = sampleBufferB;
		sample(strFunction, sampleBufferA, radres, steps, jIndex, iIndex); // overwrite (stamp) light src in BufferA
		if (sampleBufferA == sampleBufferMem)
			finished = true;
		sampleBufferMem = sampleBufferA;
		ctr++;
		std::fill(sampleBufferB.begin(), sampleBufferB.end(), std::vector<double>(steps, 0.0));
	}

	//maintainFunctionStrings(&functionStr, &coefficientArray);
	cout << "..after propagation, ctr:" << ctr << endl;
	// PROPAGATION SCHEME END //

	// TESTS START //

	//// set options
	//double tmax = 2 * pi;
	//double theta = 0.0;
	//double tInc = pi/4; // set theta increment
	//double t = 0.0;
	//double tinc = 2 * pi / 72;
	//double r = 3.0;
	//int steps = 72;

	//mu::Parser parser;
	//// parser definitions
	//parser.DefineConst("pi", pi);
	//parser.DefineVar("theta", &t);
	//
	//// create threshhold for aborting the fourier series expansion
	//double thresh = 0.0001;
	//
	//double grandSum = 0.0;
	//double grandError = 0.0;
	//// walk along circle in 16 pre-defined (heuristic)
	//for (theta = 0; theta < tmax; theta+=tInc)
	//{
	//	//theta = circularAngles.at(i);
	//	// compute continous x/y-position on grid
	//	double x = r * cos(theta);
	//	double y = r * sin(theta);

	//	// snap to nearest grid cell center
	//	int deltaX = round(x);
	//	int deltaY = round(y);

	//	int xIndex = (width - 1) / 2 + deltaX;
	//	int yIndex = (width - 1) / 2 + deltaY;

	//	// set current fString to current functionString
	//	std::string fString = functionStr.at(xIndex + yIndex * width);

	//	// set parser expression to current intensity approximation fString
	//	parser.SetExpr(fString); // set the current parser expression to fString (Fourier-Series (CH) approximation of intensity profile I(w) at current position
	//	double sum = 0.0;

	//	// sum over (around) the intensity profile I(w) to obtain mean(I) and T(w) to obtain mean(T) and mean(T*I)=mean(exp)
	//	for (t = 0.0; t < tmax; t+=tinc)
	//		sum += parser.Eval(); // evaluate fString for iMean (sum)

	//	// compute iMean from cartesian (rectangular) energy-based integral as opposed to the polar integral relevant to the geometrical (triangular/circular) area
	//	double iMean = sum / steps; // -->tinc(dt) is a constant that can be drawn out of the integral
	//	
	//	grandSum += iMean;
	//}

	//double grandMean = grandSum; // compute grand mean from individual means w. equal sample proportions and corresponding weightings

	//cout << "Mean (Sum) Intensity in center: " << circle << endl;
	//cout << "Mean (Sum) Intensity on circle (line-integral): " << grandMean << endl;
	//	
	//// ellipse parameters
	//unsigned short quality = 70;	
	//sf::ConvexShape ellipse;
	//ellipse.setPointCount(quality);
	//
	//defineConvexEllipse(&ellipse, 100 * 3.0, 100 * 3.0, quality);
	//ellipse.setPosition(wSize/2, wSize / 2); // set ellipse position
	//ellipse.setFillColor(sf::Color::Transparent);
	//ellipse.setOutlineColor(sf::Color::Green);
	//ellipse.setOutlineThickness(1);

	// POLAR GRAPHER START //
	//vector to store functions
	std::vector<polarPlot> funcs;
	std::vector<polarPlot> funcsEllipses;	

	// set options
	int Alpha = 200; // set opacity
	int drawSpeed = 0; // set instant drawing
	int lineWidth = 0; // set line width = 1px
	double thetaMax = 2 * pi;
	double thetaInc = (2*pi) / 360; // set theta increment
	double rotSpeed = 0; // set rotation speed
	sf::Color red = sf::Color(255, 0, 0);
	sf::Color green = sf::Color(0, 255, 0, Alpha);

	// add pre-computed function as string
	for (int i = 0; i < dim; i++)
	{
		funcs.push_back(polarPlot(sampleBufferA.at(i), radres, drawSpeed, lineWidth, thetaMax, thetaInc, rotSpeed, red));
		funcsEllipses.push_back(polarPlot(functionStrEllipse.at(i), drawSpeed, lineWidth, thetaMax, thetaInc, rotSpeed, green));
	}

	//declare renderwindow
	sf::RenderWindow window;
	if (fullscreen == false)
		window.create(sf::VideoMode(wSize, wSize, 32), "Polar Grapher", sf::Style::Default);
	else
		window.create(sf::VideoMode::getDesktopMode(), "Polar Grapher", sf::Style::Fullscreen);
	window.setFramerateLimit(60);
	windowsize.x = window.getSize().x;
	windowsize.y = window.getSize().y;

	//initialize function graphics
	for (int i = 0; i < funcs.size(); i++)
		{funcs.at(i).init(); funcsEllipses.at(i).init();}

	//initialize rendertexture for output image
	sf::RenderTexture out;
	out.create(wSize, wSize);

	// create flags for blending tensors(ellipses)/light distributions(polar plots)
	bool showOverlay{ true };
	bool showBase{ true };

	//main loop
	while (window.isOpen())
	{
		sf::Event event;
		// query window poll events
		while (window.pollEvent(event))
		{
			if (event.type == sf::Event::Closed || event.type == sf::Event::KeyPressed && event.key.code == sf::Keyboard::Escape)
				window.close();
			if (event.type == sf::Event::KeyPressed)
			{
				if (event.key.code == sf::Keyboard::Num1) // press 1 to toggle base/clip image
					showBase = !showBase;
				else if (event.key.code == sf::Keyboard::Num2) // press 2 to toggle overlay image
					showOverlay = !showOverlay;
			}
		}

		// query opened winow, to break from loop to prevent from overwriting the framebuffer w. black (0,0,0)
		if (!window.isOpen())
			break;

		window.clear(sf::Color::Black);
		out.clear(sf::Color::Black);		

		// draw polar function as graph sprite
		for (int i = 0; i < dim; i++) 
		{
			double offsetX = (i % width)*(wSize / width) + (wSize / (2 * width)); // check rest (modulo) for x-offset
			double offsetY = (i / width)*(wSize / width) + (wSize / (2 * width)); // check division for y offset
			// ellipse parameters
			unsigned short quality = 70;
			// define ellipses as convex shapes (with semi-axes u scaled to the specific singular values
			sf::ConvexShape ellipse;
			ellipse.setPointCount(quality);
			defineConvexEllipse(&ellipse, 2.0, 2.0, quality);
			ellipse.setPosition(offsetX, offsetY); // set ellipse position

			funcs.at(i).animation(i, 0, sampleBufferA);
			funcsEllipses.at(i).animation(i, 1);
			sf::Sprite spr = funcs.at(i).update(); // draw sprites[Kobolde/Elfen] (composited bitmaps/images - general term for objects drawn in the framebuffer)
			sf::Sprite sprE = funcsEllipses.at(i).update(); // draw sprites[Kobolde/Elfen] (composited bitmaps/images - general term for objects drawn in the framebuffer)
			spr.setPosition(wSize / 2, wSize / 2);
			window.draw(ellipse);
			if (showBase)
				{window.draw(spr); out.draw(spr);}
			if (showOverlay)
				{window.draw(sprE); out.draw(sprE);}
		}
		// window.draw(ellipse);  // TEST //
		// record(out, record_frameskip, record_folder); // use frameskip to define recording frameskip
		window.display(); // update window texture
	}

	//write to image file
	out.display(); // update output texture
	sf::Texture outTx = out.getTexture();
	sf::Image outImg = outTx.copyToImage();
	outImg.saveToFile("out.png");

	// POLAR GRAPHER END //

	system("PaUsE");

	return EXIT_SUCCESS;
}