#define _USE_MATH_DEFINES
#include <SFML/Graphics.hpp>
#include <iostream>
#include "muParser.h"
#include <cmath>
#include <sstream>
#include <fstream>
#include <vector>
#include <unistd.h>
#include <math.h>
#include "exprtk.hpp"
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <functional>
#include <array>
#include <filesystem>
#include <sstream>
#include <string>

#ifdef WINDOWS
#include <direct.h>
#define GetCurrentDir _getcwd
#else
#include <unistd.h>
#define GetCurrentDir getcwd
#endif

#define MAXBUFSIZE  ((int) 1e3)

using namespace std;
using namespace Eigen;
using namespace mu;

// PROTOTYPES
std::string functionString; // Prototyping of functionString for strFunction!
int width;
int height;

//length and width of window
int wSize = 700;
sf::Vector2i windowsize;

//recording frameskip 
int record_frameskip = 0; // --> disables recording

//directory to write images to, must exist
std::string record_folder = "frames";

//fullscreen flag
bool fullscreen = false;

//definition of pi
const double pi = M_PI;

// get current working directory to assign matrix.txt path
std::string GetCurrentWorkingDir(void) {
	char buff[FILENAME_MAX];
	GetCurrentDir(buff, FILENAME_MAX);
	std::string current_working_dir(buff);
	return current_working_dir;
}

void defineConvexEllipse(sf::ConvexShape* ellipse, double radius_x, double radius_y, unsigned short quality)
{
	for (unsigned short i = 0; i < quality; ++i)
	{
		double rad = (360 / quality * i) / (360 / M_PI / 2);
		double x = cos(rad)*radius_x;
		double y = sin(rad)*radius_y;

		ellipse[0].setPoint(i, sf::Vector2f(x, y));
	}
}


// OPERATORS
template <typename T>
T clip(const T& n, const T& lower, const T& upper) // template clip function
{
	return std::max(lower, std::min(n, upper));
}

// -> muparser custom clip function: template functions need to be static callback (call after) functions, which can be passed as arguments
static value_type clipNeg(value_type v) { return clip(v, 0.0, 1.0); }

//template <typename T1, typename T2>
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

myPair lightSrcPos;

// OPERATOR DEFINITIONS //

void parse_options(int argc, char* argv[]) {

	int c;
	std::string frameskip_opt = "-1";
	std::string size_opt = "-1";
	std::string directory_opt = "";
	std::string light_opt = "-1";
	//	extern char *optarg;
	//	extern int optind, optopt;

	//using getopts 
	while ((c = getopt(argc, argv, "fs:r:d:l:")) != -1) // --> classic getopt call: argc: # of args; argv: array of args(strings); "": optstring (registered options, followed by colon when taking args itself) 
	{
		int f = -1, s = -1;
		// create pair for light src position (j,i)
		lightSrcPos = { width / 2, width / 2 };

		switch (c) 
		{
		// correct option use
		case 'l': {
			light_opt.assign(optarg);
			std::istringstream(light_opt) >> lightSrcPos;
			break;
		}
		case 'f': {
			fullscreen = true;
			break; }
		case 'r':
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
			}
			break;
		case '?':
			std::cerr << "Unrecognized option: '-" << optopt << "\n";
			break;
		}
	}
}


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

// CLASS DEFINITIONS //
class neighborhood
{
	bool neighborR = true;
	bool neighborT = true;
	bool neighborL = true;
	bool neighborB = true;

	std::vector<bool> propagationCue = { true,true,true,true };

public:

	//constructor
	neighborhood(int j, int i, const int length, std::vector<bool>& processMap) // check the neighbourhood of (j,i) for missing neighbors..
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

		// check neighbors in process(ed) map, for neighbors already processed ..
		if (neighborR) // if neighbor existent..
			if (processMap.at(i + 1 + j * width) == true)
				neighborR = false; // right
		if (neighborT)
			if (processMap.at(i + (j - 1) * width) == true)
				neighborT = false; // top
		if (neighborL)
			if (processMap.at(i - 1 + j * width) == true)
				neighborL = false; // left
		if (neighborB)
			if (processMap.at(i + (j + 1) * width) == true)
				neighborB = false; // bottom
		// ordered again..
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
	sf::Color color;

public:

	//constructor
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
		double y = -1*scale * r * sin(t);


		if (x + 1 + wSize / 2 > wSize || x - 1 + wSize / 2 < 0 || y + 1 + wSize / 2 > wSize || y - 1 + wSize / 2 < 0)
			return;

		int offsetX = (index % width)*(wSize / width) + (wSize / (2 * width)); // check rest (modulo) for x-offset
		int offsetY = (index / width)*(wSize / width) + (wSize / (2 * width)); // check division for y offset

		int xIndex = round(x + offsetX);
		int yIndex = round(y + offsetY);

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
		for (int i = 0; (i < speed || speed == 0) && t <= tmax; i++)
		{
			t += tinc;
			r = p.Eval();
			// PROPAGATION ATTENUATION TEST //
			/*if (t == tinc && index/width == 3 && mode == 0)
			{
				cout << "r: " << r << endl;
				cout << "attenuation: " << r / 1.0 << endl;
			}*/
			plot(index, mode);
		}
	}
};

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



// mathematical functions
double circleFunction(double x)
{
	double os = 1.0 / sqrt(M_PI);
	return os;
}

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

double integral(double(*f)(double x), double a, double b, int r)
{
	double step = (b - a) / r;  // width of each small rectangle
	double area = 0.0;  // signed area
	for (int i = 0; i < r; i++)
		area += f(a + (i + 0.5) * step) * step; // sum up each small rectangle

	return area;
}

double aIntegral(double(*f)(double x), double a, double b, int r, int n)
{
	double step = (b - a) / r;  // width of each small rectangle
	double area = 0.0;  // signed area

	for (int i = 0; i < r; i++)
		area += f(a + (i + 0.5) * step)*cos(n*(a + (i + 0.5) * step)) * step; // sum up each small rectangle

	return area;
}

double bIntegral(double(*f)(double x), double a, double b, int r, int n)
{
	double step = (b - a) / r;  // width of each small rectangle
	double area = 0.0;  // signed area
	for (int i = 0; i < r; i++)
		area += f(a + (i + 0.5) * step)*sin(n*(a + (i + 0.5) * step)) * step; // sum up each small rectangle

	return area;
}

std::array<std::array<double, 21>, 2> calcCoeff(double(*f)(double x)) // for 2Pi periodic functions
{
	double a0 = 1 / (2 * M_PI) * integral(f, -M_PI, M_PI, 50);
	
	std::array<double, 21> an; // use static array for better performance
	an.fill(0.0);
	std::array<double, 21> bn; // use static array for better performance
	bn.fill(0.0);
	std::array<std::array<double, 21>, 2> coeff; // init packed (zipped) coefficient vector
	coeff.fill(an);
	double thresh = 0.0001;

	an.at(0) = a0; // push a0 in both lists on frequency 0 (offset)
	bn.at(0) = a0; // push a0 in both lists on frequency 0 (offset)
	for (int i = 1; i < an.size(); i++)
	{
		an.at(i) = 1 / M_PI * aIntegral(f, -M_PI, M_PI, 50, i); // calculate Integrals for the coeff
		bn.at(i) = 1 / M_PI * bIntegral(f, -M_PI, M_PI, 50, i); // calculate Integrals for the coeff

		/*cout << "an: " << an.at(i) << endl;
		cout << "bn: " << bn.at(i) << endl;*/
		if ((bn.at(i) < thresh && an.at(i) < thresh && an.at(i - 1) < thresh && bn.at(i - 1) < thresh)) // if coeff < thresh or #>20, abort (finish) series expansion!
			break;
	}
	coeff.front() = an;
	coeff.back() = bn;

	return coeff;
}

void countMatrix(std::string filepath, int* colsCount, int* rowsCount)
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
};

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

class propagator
{
	// define cosine lobes (4 neighbors/directions)
	std::string cosLobeX = "clip(cos(theta))"; // cosine lobe in +x-direction (in polar plot)
	std::string cosLobeY = "clip(cos(theta-pi/2))"; // cosine lobe in +y-direction (in polar plot)
	std::string cosLobeXneg = "clip(cos(theta+pi))"; // cosine lobe in -x-direction (in polar plot)
	std::string cosLobeYneg = "clip(cos(theta+pi/2))"; // cosine lobe in -y-direction (in polar plot)
	// std::array<std::string, 12> centralLobes{ cosLobeYneg, cosLobeX, cosLobeY, cosLobeX, cosLobeY, cosLobeXneg, cosLobeY, cosLobeXneg, cosLobeYneg, cosLobeXneg, cosLobeYneg, cosLobeX };

	// create arrays for storing pre-computed fourier (CH) coefficients for cosine lobes
	std::array<std::array<double, 21>, 2> cosLobeXCoeff;
	std::array<std::array<double, 21>, 2> cosLobeYCoeff;
	std::array<std::array<double, 21>, 2> cosLobeXnegCoeff;
	std::array<std::array<double, 21>, 2> cosLobeYnegCoeff;

	// set principal arture angles in rads
	double alpha = 36.8699 * M_PI / 180;
	double beta = 26.5651 * M_PI / 180;	
	// define principal central directions array (12 central directions in 2D --> 30 in 3D, ufff..)
	std::array<double, 12> centralDirections{ 3.0 / 2 * M_PI + 1.0172232666228471, 0, 0.5535748055013016, 1.0172232666228471, M_PI / 2, M_PI / 2 + 0.5535748055013016, M_PI / 2 + 1.0172232666228471, M_PI, M_PI + 0.5535748055013016, M_PI + 1.0172232666228471, 3.0 / 2 * M_PI, 3.0 / 2 * M_PI + 0.5535748055013016 };
	std::array<double, 12> apertureAngles{ beta, alpha, beta, beta, alpha, beta, beta, alpha, beta, beta, alpha, beta }; // define aperture angles array in order
	std::array<std::array<std::array<double, 21>, 2>*, 12> centralLobesCoeff{ &cosLobeYnegCoeff, &cosLobeXCoeff, &cosLobeYCoeff, &cosLobeXCoeff, &cosLobeYCoeff, &cosLobeXnegCoeff, &cosLobeYCoeff, &cosLobeXnegCoeff, &cosLobeYnegCoeff, &cosLobeXnegCoeff, &cosLobeYnegCoeff, &cosLobeXCoeff };

	// define member strings for parsing (internal)
	std::string fString = "0.0";
	std::string fStringEllipse = "0.0";

	// create member vectors (arrays) for storing the symbolic strings and fourier (CH) coefficients (external)
	std::vector<std::string>* functionStr; // initialize w. 0.0 to prevent parser errors
	std::vector<std::string>* functionStrEllipse; // initialize w. 0.0 to prevent parser errors
	std::vector<std::array<std::array<double, 21>, 2>>* coefficientArray;
	std::vector<std::tuple<double, double, double>>* ellipseArray;

	// define function parser (muparser)
	mu::Parser parser;
	double t = 0.0;
	double tinc = 2*pi / 72; // set theta increment tinc
	int steps = (2 * pi) / tinc;
	
	// create threshhold for aborting the fourier series expansion
	double thresh = 0.0001;
	int width = 0;

	// create fourier (CH) coefficient vectors
	std::array<double, 21> an;
	std::array<double, 21> bn;// (21, 0.0);

	// create process(ed) map for cells already processed
	std::vector<bool> processMap; // create a binary process(ed) map

public:
	propagator(const int dim, int Width, std::vector<std::string>* fStr, std::vector<std::array<std::array<double, 21>, 2>>* coeffArray, std::vector<std::string>* fStrEllipse, std::vector<std::tuple<double, double, double>>* ellipseArr) : processMap(dim, false)
	{
		// parser definitions
		parser.DefineConst("pi", pi);
		parser.DefineVar("theta", &t);
		parser.DefineFun(_T("clip"), clipNeg, false);

		// initialize fourier (CH) coefficient arrays w. 0.0
		an.fill(0.0);
		bn.fill(0.0);

		// assign ptrs to member vectors
		functionStr = fStr;
		functionStrEllipse = fStrEllipse;
		coefficientArray = coeffArray;
		ellipseArray = ellipseArr;

		// initialize pre-computed fourier coefficient arrays (using wolfram mathematica)
		cosLobeXCoeff.front() = { 0.31831, 0.5, 0.212207, 0., -0.0424413, 0., 0.0181891, 0., -0.0101051, 0., 0.0064305, 0., -0.00445189, 0., 0.00326472, 0., -0.00249655, 0., 0.00197096, 0., -0.00159554 };
		cosLobeXCoeff.back().fill(0); // initialize bn array w. 0 --> no sine coefficients
		cosLobeYCoeff.front() = { 0.31831, 0., -0.212207, 0., -0.0424413 , 0., -0.0181891 , 0., -0.0101051, 0., -0.0064305 , 0., -0.00445189, 0., -0.00326472, 0., -0.00249655, 0., -0.00197096 , 0., -0.00159554 };
		cosLobeYCoeff.back() = { 0., 0.5, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. };
		cosLobeXnegCoeff.front() = { 0.31831, -0.5, 0.212207, 0., -0.0424413, 0., 0.0181891, 0., -0.0101051, 0., 0.0064305, 0., -0.00445189, 0., 0.00326472, 0., -0.00249655, 0., 0.00197096, 0., -0.00159554 };
		cosLobeXnegCoeff.back().fill(0); // initialize bn array w. 0 --> no sine coefficients
		cosLobeYnegCoeff.front() = { 0.31831, 0., -0.212207, 0., -0.0424413 , 0., -0.0181891 , 0., -0.0101051, 0., -0.0064305 , 0., -0.00445189, 0., -0.00326472, 0., -0.00249655, 0., -0.00197096 , 0., -0.00159554 };
		cosLobeYnegCoeff.back() = { 0., -0.5, 0., 0., 0., 0., 0., 0., 0., 0., 0., 0. };


		// assign width
		width = Width;
	}

	void propagate(int jIndex, int iIndex)
	{
		// create neighborhood, for tagging non-existent neighbors..
		neighborhood hood(jIndex, iIndex, width, processMap);

		// assign member coefficient arrays for evaluating current light profile... get from current position (jIndex,iIndex)
		an = coefficientArray->at(iIndex + jIndex * width).front(); // ..cosine coefficient vector
		bn = coefficientArray->at(iIndex + jIndex * width).back(); // ..sine coefficient vector		 

		// calculate mean and variance.. of I(phi)
		double sum1 = 0.0;
		double sum2 = 0.0;
		double sum3 = 0.0;

		// assign singular values and rotational deg(ree)
		double sv1 = std::get<0>(ellipseArray->at(iIndex + jIndex * width));
		double sv2 = std::get<1>(ellipseArray->at(iIndex + jIndex * width));
		double deg = std::get<2>(ellipseArray->at(iIndex + jIndex * width));

		// sum over (around) the intensity profile I(w) to obtain mean(I) and T(w) to obtain mean(T) and mean(T*I)=mean(exp)
		for (int i = 0; i < steps; i++)
		{
			double sum = an.at(0); // initialize sum w. offset a0
			// evaluate fourier coefficient array for magnitude at current angular position i*tinc
			for (int j = 1; j < an.size(); j++)
				sum += an.at(j)*cos(j*i*tinc) + bn.at(j)*sin(j*i*tinc);
			sum1 += sum; // evaluate for iMean (sum1)

			if(sv1 != 0 && sv2 != 0)
				sum2 += sv1 * sv2 / sqrt(sv2*sv2*cos(i*tinc - deg)*cos(i*tinc - deg) + sv1 * sv1*sin(i*tinc - deg)*sin(i*tinc - deg)); //--> ellipse equation, evaluate for tMean (sum2)
			if (sv1 != 0 && sv2 != 0)
				sum3 += sum * sv1 * sv2 / sqrt(sv2*sv2*cos(i*tinc - deg)*cos(i*tinc - deg) + sv1 * sv1*sin(i*tinc - deg)*sin(i*tinc - deg)); // evaluate T*I for tiMean (sum2)
		}
		
		// compute iMean from cartesian (rectangular) energy-based integral as opposed to the polar integral relevant to the geometrical (triangular/circular) area
		double iMean = sum1 / steps; // -->tinc(dt) is a constant that can be drawn out of the integral
		// compute mean(T) from cartesian (rectangular) energy-based integral as opposed to the polar integral relevant to the geometrical (triangular/circular) area	
		double tMean = sum2 / steps; // -->tinc(dt) is a constant that can be drawn out of the integral
		// compute mean(T*I) from cartesian (rectangular) energy-based integral as opposed to the polar integral relevant to the geometrical (triangular/circular) area
		double tiMean = sum3 / steps; // -->tinc(dt) is a constant that can be drawn out of the integral
		// compute correction factor (scaling to mean=1, subsequent scaling to mean(I)), which follows energy conservation principles
		double cFactor = tMean*iMean / tiMean;

		// cout << "cFactor: " << cFactor << endl;

		// define intensity [W/rad] and area [W]
		double area = 0.0;

		// iterate through central directions array to distribute (spread) energy (intensity) to the cell neighbors
		for (int k = 0; k < centralDirections.size(); k++) // k - (Strahl-)keulenindex
		{
			int nIndex = k / 3; // create index to capture the 4 neighbors (in 2D)

			// check the neighborhood for missing (or already processed) neighbors, if missing, skip step..continue
			if (!hood.getR() && nIndex == 0)
				continue;
			if (!hood.getT() && nIndex == 1)
				continue;
			if (!hood.getL() && nIndex == 2)
				continue;
			if (!hood.getB() && nIndex == 3)
				continue;

			// set theta variable to current central direction shifted by half an apertureAngle
			double offset = centralDirections.at(k) - apertureAngles.at(k) / 2.0;
			double area = 0.0; // -->reset energy variables
		
			// compute lSteps to determine # of steps for current apertureAngle
			int lSteps = apertureAngles.at(k) / tinc;

			// integrate over the profile T(w)*I(w) to obtain total intensity received by respective face (of current neighbor nIndex)
			for (long long i = 0; i < lSteps; i++) // --> carry out the integral computation as discrete (weghted) sum
			{
				double intensity = an.at(0); // initialize sum w. offset a0
				// evaluate fourier coefficient array for magnitude at current angular position i*tinc
				for (long long j = 1; j < an.size(); j++)
					intensity += (an.at(j)*cos(j*i*tinc + offset) + bn.at(j)*sin(j*i*tinc + offset))*sv1 * sv2 / sqrt(sv2*sv2*cos(j*i*tinc + offset - deg)*cos(j*i*tinc + offset - deg) + sv1 * sv1*sin(j*i*tinc + offset - deg)*sin(j*i*tinc + offset - deg));

				area += intensity; // integrate over angle in cart. coordinates (int(I(w),0,2Pi) to obtain total luminous flux (power) received by adjacent cell faces
			}

			// calculate magnitude w. corresponding correction factor and weighting by tinc, which is is a constant dt drawn out of the discrete sum
			double mag = cFactor*tinc*area / 2.0; // since int(cos(w)) (w. negative values clamped to zero) yields = 2 --> normalize to energy 1 (spread energy area over the cosine)

			// ReProjection of Flow Lobes //

			switch (nIndex) // switch nIndex to cope with changing neighbor locations..., accumulate fourier coefficients in coefficientArray
			{
			case 0: coefficientArray->at(iIndex + 1 + jIndex * width) = coefficientArray->at(iIndex + 1 + jIndex * width) + mag * *centralLobesCoeff.at(k); break; // right neighbor
			case 1: coefficientArray->at(iIndex + (jIndex - 1) * width) = coefficientArray->at(iIndex + (jIndex - 1) * width) + mag * *centralLobesCoeff.at(k); break; // top neighbor
			case 2: coefficientArray->at(iIndex - 1 + jIndex * width) = coefficientArray->at(iIndex - 1 + jIndex * width) + mag * *centralLobesCoeff.at(k); break;// left neighbor
			case 3: coefficientArray->at(iIndex + (jIndex + 1) * width) = coefficientArray->at(iIndex + (jIndex + 1) * width) + mag * *centralLobesCoeff.at(k); break;// bottom neighbor
			}
		}
		// set flag in process(ed) map for cell already processed
		processMap.at(iIndex + jIndex * width) = true;
	}
};

void maintainFunctionStrings(std::vector<std::string>* fStr, std::vector<std::array<std::array<double, 21>, 2>>* coeffArray)
{
	for (int i = 0; i < fStr->size(); i++)
	{
		double thresh = 0.0001;
		// assign member coefficient arrays for evaluating current light profile... get from current position (jIndex,iIndex)
		std::array<double, 21> an = coeffArray->at(i).front(); // ..cosine coefficient vector
		std::array<double, 21> bn = coeffArray->at(i).back(); // ..sine coefficient vector

		// create fString w. a0 (offset)
		std::string fString = std::to_string(an.front());

		// set up Fourier Series expansion as str expression for the evaluation of the fourier serier (polar function/curve) for each grid cell (j,i)
		for (int i = 1; i < an.size(); i++)
		{
			if (i > 1 && (bn.at(i) < thresh && an.at(i) < thresh && an.at(i - 1) < thresh && bn.at(i - 1) < thresh)) // if threshold exceeded..
				break; // abort Fourier series expansion to obtain shorter functionStrings (approximations)
			fString += "+" + std::to_string(an.at(i)) + "*cos(" + std::to_string(i) + "*theta)" + "+" + std::to_string(bn.at(i)) + "*sin(" + std::to_string(i) + "*theta)";
		}

		// update current functionStr w. fString to maintain the polar functions in gríd
		fStr->at(i) = fString;
	}
}

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
	double tinc = 2*pi / 72; // set theta increment tinc
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
			{sum += parser.Eval(); t += tinc;}

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

	// create fourier (CH) coefficients array (grid)
	std::array<double, 21> init; init.fill(0.0); // create a std::vector w. n elements initialized to 0..
	std::array<std::array<double, 21>, 2> initVector; // create an initArray to represent packed (zipped) fourier coeff. vectors
	initVector.fill(init);
	std::vector<std::array<std::array<double, 21>, 2>> coefficientArray(dim, initVector); // ..use it to initialize the coefficient array w. dim elements
	
	std::tuple<double, double, double> initTuple = { 0.0,0.0,0.0 }; // -->triple: lambda1, lambda2, deg
	std::vector<std::tuple<double,double,double>> ellipseArray(dim, initTuple); // ..use it to initialize the coefficient array w. dim elements

	cout << "before compute glyphs" << endl;
	// compute Eigenframes/Superquadrics/Ellipses/Glyphs by calling computeGlyphs w. respective args
	computeGlyphs(functionStrEllipse, ellipseArray);
	
	// define excitation (stimulus) polar functions(4 neighbors/directions) normalized to area 1
	std::string circle = "1.0"; // circle w. 100% relative intensity - capturing the spread [0..1] normalized to strongest light src in field
	
	// write light src in coefficientArray
	functionString = circle; // set circular (isotroic) profile for light src
	int j = lightSrcPos.jIndex; int i = lightSrcPos.iIndex; // create light src position indices
	coefficientArray.at(i + j * width) = calcCoeff(strFunction);// form coefficient arrays for evaluating current light profile... get from current position (jIndex,iIndex)

	// compute distances to center point
	int deltaJ = abs(width / 2 - j);
	int deltaI = abs(width / 2 - i);
	// create minRadius to determine propRadius (propagation radius: necessary to reach all cells!)
	int minRadius = width - 1; // use larger grid radius to reach outer cells in diagonal propagation - adjust for varying light src position
	// create propagation radius for propagation scheme - deltas added to minRadius!
	int propRadius;
	if (cols%2 == 0) // if even grid
		propRadius = deltaI + deltaJ + minRadius + 1;
	else // ..odd grid
		propRadius = deltaI + deltaJ + minRadius;

	// create propagator object (managing propagation, reprojection, correction, central directions, apertureAngles and more...)
	propagator prop(dim, width, &functionStr, &coefficientArray, &functionStrEllipse, &ellipseArray);

	// center 4- neighborhood - start from light src, walk along diagonals to obtain orthogonal propagation and to cope with ray effects (discretization)!!!
	prop.propagate(j, i); // propagate from the light src position
	int jIndex = j;
	int iIndex = i;

	cout << "before propagation..propagating.." << endl;

	// PROPAGATION SCHEME START //
	for (int p = 0; p <= propRadius; p++) // ..walk along diagonals to encircle point light for unbiased radial propagation
	{
		if (p != 1) // remember: step 0 is the first right step, step 1 already there...
		{
			jIndex += 0; iIndex += 1; // 0 degrees... -> step right once each iteration, except in step 1
			if (iIndex >= 0 && jIndex >= 0 && jIndex < height && iIndex < width) // if inside frame (window)..
				prop.propagate(jIndex, iIndex); // .. propagate from this cell
		}

		for (int l = 0; l < p - 1; l++) // perform p steps of incrementation/decrementation
		{
			jIndex -= 1; iIndex += 1; // + 45 degrees
			if (iIndex >= 0 && jIndex >= 0 && jIndex < height && iIndex < width) // if inside frame (window)..
				prop.propagate(jIndex, iIndex); // .. propagate from this cell
		}

		for (int l = 0; l < p; l++) // perform p steps of incrementation/decrementation
		{
			jIndex -= 1; iIndex -= 1; // + 135 degrees
			if (iIndex >= 0 && jIndex >= 0 && jIndex < height && iIndex < width) // if inside frame (window)..
				prop.propagate(jIndex, iIndex); // .. propagate from this cell
		}

		for (int l = 0; l < p; l++) // perform p steps of incrementation/decrementation
		{
			jIndex += 1; iIndex -= 1; // - 135 degrees
			if (iIndex >= 0 && jIndex >= 0 && jIndex < height && iIndex < width) // if inside frame (window)..
				prop.propagate(jIndex, iIndex); // .. propagate from this cell
		}

		for (int l = 0; l < p; l++) // perform p steps of incrementation/decrementation
		{
			jIndex += 1; iIndex += 1; // - 45 degrees
			if (iIndex >= 0 && jIndex >= 0 && jIndex < height && iIndex < width) // if inside frame (window)..
				prop.propagate(jIndex, iIndex); // .. propagate from this cell
		}
	}

	maintainFunctionStrings(&functionStr, &coefficientArray);
	cout << "..after propagation" << endl;
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
	//// walk along circle in 16 pre-defined steps (heuristic)
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
	double thetaInc = pi / 300; // set theta increment
	double rotSpeed = 0; // set rotation speed
	sf::Color red = sf::Color(255, 0, 0);
	sf::Color green = sf::Color(0, 255, 0, Alpha);

	// add pre-computed function as string
	for (int i = 0; i < width*height; i++)
	{
		funcs.push_back(polarPlot(functionStr.at(i), drawSpeed, lineWidth, thetaMax, thetaInc, rotSpeed, red));
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
	for (unsigned i = 0; i < funcs.size(); i++)
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
		for (unsigned i = 0; i < funcs.size(); i++) 
		{
			funcs.at(i).animation(i, 0);
			funcsEllipses.at(i).animation(i, 1);
			sf::Sprite spr = funcs.at(i).update(); // draw sprites[Kobolde/Elfen] (composited bitmaps/images - general term for objects drawn in the framebuffer)
			sf::Sprite sprE = funcsEllipses.at(i).update(); // draw sprites[Kobolde/Elfen] (composited bitmaps/images - general term for objects drawn in the framebuffer)
			spr.setPosition(wSize / 2, wSize / 2);
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