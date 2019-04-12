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
//#include "exprtk.hpp"
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <functional>
#include <boost/iostreams/stream.hpp>
#include <boost/iostreams/tee.hpp>
typedef boost::iostreams::tee_device<std::ostream, std::ofstream> Tee;
typedef boost::iostreams::stream<Tee> TeeStream;

#define MAXBUFSIZE  ((int) 1e3)
#ifdef WINDOWS
#include <direct.h>
#define GetCurrentDir _getcwd
#else
#include <unistd.h>
#define GetCurrentDir getcwd
#endif

using namespace std;
using namespace Eigen;

//length and width of window
int wSize = 700;
std::string functionString = "2.1"; // Prototyping of functionString for strFunction (symbolic string arg parsed by muparser in cpp functional ptr rep)!
sf::Vector2i windowsize;
int width;
int height;
int steps; // use n steps for angular resolution - radres
double radres;

//definition of pi
const double pi = M_PI;

template <typename T>
T muPlus(T& a, T& b)
{
	T res = a + b;
	return res;
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

	//constructor
	neighborhood(int j, int i, const int length, std::vector<bool>& processMap) // check the neighbourhood of (j,i) for missing neighbors..
	{
		// check if frame exceeded, outside FOV, offscreen..
		if (i == 0)
			neighborL = false; // left
		else if (i == length - 1)
			neighborR = false; // right
		if (j == 0)
			neighborT = false; // top
		else if (j == length - 1)
			neighborB = false; // bottom
		// leave order.. functional!

		// check neighbors in process(ed) map, for neighbors already processed ..
		if (neighborR) // if neighbor existent..
			if (processMap.at(i + 1 + j * length) == true)
				neighborR = false; // right
		if (neighborT)
			if (processMap.at(i + (j - 1) * length) == true)
				neighborT = false; // top
		if (neighborL)
			if (processMap.at(i - 1 + j * length) == true)
				neighborL = false; // left
		if (neighborB)
			if (processMap.at(i + (j + 1) * length) == true)
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
	int steps;
	int ctr = 0;
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

	//constructor #1 string parser - muParser: tinc_in passed as tinc used for parsing
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

	//constructor #2 array parser w. array-inherited resolution radres
	polarPlot(std::vector<std::vector<double>> SampleBuffer,
		int speed_in,
		int line_width_in,
		int steps_in,
		double tmax_in,
		double res,
		double rotation_in,
		sf::Color color_in)

		: speed(speed_in),
		line_width(line_width_in),
		steps(steps_in),
		t(0),
		tmax(tmax_in),
		radres(res),
		rotation(rotation_in),
		color(color_in)
	{
		discretePlot = true;
	}

	////constructor #3 array parser w. 2 different resolutions, array-inherited resolution radres and plot resolution tinc_in, usually coarser resolutions for performance reasons
	//polarPlot(std::vector<double> Sample,
	//	int speed_in,
	//	int line_width_in,
	//	int steps_in,
	//	double tmax_in,
	//	double tinc_in,
	//	double res,
	//	double rotation_in,
	//	sf::Color color_in)

	//	: speed(speed_in),
	//	line_width(line_width_in),
	//	steps(steps_in),
	//	t(0),
	//	tmax(tmax_in),
	//	tinc(tinc_in),
	//	sample(Sample),
	//	rotation(rotation_in),
	//	radres(res),
	//	color(color_in)
	//{
	//	discretePlot = true;
	//}

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
	void plot(int index = 0, int mode = 0)
	{
		// line width
		if (r >= 1.0)
			line_width = line_width + 1;

		//scale
		double scale;
		if (mode == 1)
			scale = 10;// / 360.0;
		else
			scale = 50;
		//convert polar to cartesian
		double x = scale * r * cos(t);
		double y = -1 * scale * r * sin(t);

		if (x + 1 + wSize / 2 > wSize || x - 1 + wSize / 2 < 0 || y + 1 + wSize / 2 > wSize || y - 1 + wSize / 2 < 0)
		{
			line_width -= 1; return;
		}

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
		if ((line_width == 1 || line_width == 2) && xIndex > 0 && xIndex < wSize - 1 && yIndex > 0 && yIndex < wSize - 1) {
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

		if (r >= 1.0)
			line_width = line_width - 1;

	}

	// STRING PARSER OVERLOAD of animation - use w. string/mixed resolution constructor overload
	void animation(int index = 0, int mode = 0)
	{

		for (int i = 0; (i < speed || speed == 0) && t < tmax; i++)
		{
			if (!discretePlot) // use muParser for evaluation 
				r = p.Eval();
			else // nearest neighbor interpolation w. samples in sampleVector - ONLY FOR DUAL GRID PLOT! - use animation overload for cell-based plot
				r = sample.at(static_cast<int>(std::round(t / radres)) % sample.size()); // cyclic value permutation

			plot(index, mode);
			t += tinc;
		}
	}
	// DISCRETE SAMPING OVERLOAD for animation... faster and no string parsing necessary - use for inherited resolution constructor overload
	void animation(int index, int mode, std::vector<std::vector<double>>& sampleArray, TeeStream* both = NULL)
	{
		t = 0.0; // reset theta variable for each cell to prevent from long-term shifting

		for (int i = 0; (i < speed || speed == 0) && i < steps; i++)
		{
			r = sampleArray.at(index).at(i);

			// PROPAGATION ATTENUATION TEST //
			if (i == 0 && index / width == height / 2 && mode == 0 && ctr == 0)
			{

				/**both << "r: " << r << endl;
				*both << "attenuation: " << r / 1.0 << endl;
				ctr++;*/
			}

			/*if(t > 45*tinc && t < 47*tinc && index/width == width/2 - 1 && index%width == width/2 + 1 )
			{
				cout << "r(diag): " << r << endl;
				cout << "attenuation: " << r / 1.0 << endl;
			}*/

			plot(index, mode);
			t += radres;
		}
	}

};

//writes frames to file 
//directory must already exist -- parse file methods
void record(sf::RenderTexture& in, int frameskip, std::string directory)
{

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


// create options for getopt(.c)
Pair lightSrcPos{ height / 2, width / 2 };
bool fullscreen = false; //fullscreen flag
std::string record_folder = "frames";//directory to write images to, must exist
int record_frameskip = 0; // --> disables recording //recording frameskip
double intensity = 2.1; // --> initial intensity val
std::string workDir;
double thresh = 0.001;
bool total_anisotropy = false;
int ctrLimit = 0;

// parse files
void parse_file(char* filename, std::vector<std::string>& funcs, std::vector<Pair>& positions) {

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
// mathematical functions
double circleFunction(double x)
{
	double os = 1.0/sqrt(M_PI);
	return os;
}

//// template function for evaluating string functions via ptr 
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

std::vector<std::vector<double>> calcCoeff(double(*f)(double x)) // for 2Pi periodic functions
{
	double a0 = 1 /(2*M_PI) * integral(f, 0, 2*M_PI, 50);
	cout << "a0: " << a0 << endl;
	bool finish = false;
	std::vector<double> an(21, 0.0);
	std::vector<double> bn(21, 0.0);
	std::vector<std::vector<double>> coeff(2, an); // init packed (zipped) coefficient vector
	double thresh = 0.0001;

	an.at(0) = a0; // push a0 in both lists on frequency 0 (offset)
	bn.at(0) = a0; // push a0 in both lists on frequency 0 (offset)
	for (int i = 1; i < an.size(); i++)
	{
		an.at(i) = 1/M_PI * aIntegral(f, 0, 2*M_PI, 50, i); // calculate Integrals for the coeff
		bn.at(i) = 1/M_PI * bIntegral(f, 0, 2*M_PI, 50, i); // calculate Integrals for the coeff
	
		cout << "an: " << an.at(i) << endl;
		cout << "bn: " << bn.at(i) << endl;
		if ((bn.back() < thresh && an.back() < thresh && an.at(an.size() - 2) < thresh && bn.at(bn.size() - 2) < thresh)) // if coeff < thresh or #>20, abort (finish) series expansion!
			break;
	}
	cout << "an_size(=bn_size): " << an.size() << endl;
	coeff.front() = an;
	coeff.back() = bn;

	return coeff;
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

void computeGlyphs(std::vector<std::vector<double>>& glyphBuffer, std::vector<std::vector<bool>>& signMap, std::vector<std::vector<double>>& glyphParameters)
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
		glyphParameters.at(i).at(0) = 1.0 / rMean * sv1;
		glyphParameters.at(i).at(1) = 1.0 / rMean * sv2;

		std::transform(glyphBuffer.at(i).begin(), glyphBuffer.at(i).end(), glyphBuffer.at(i).begin(), std::bind(std::multiplies<double>(), std::placeholders::_1, 1.0 / rMean));
	}
}

void defineConvexEllipse(sf::ConvexShape* ellipse, double radius_x, double radius_y, unsigned short quality, double rot = 0.0)
{
	for (int i = 0; i < quality; ++i)
	{
		double rad = (360 / quality * i) / (360 / M_PI / 2);
		double x = cos(rad)*radius_x;
		double y = sin(rad)*radius_y;

		ellipse[0].setPoint(i, sf::Vector2f(x, y));
	}
}

int main(int argc, char* argv[])
{
	// SVD/EIGENVALUE DECOMPOSITION START //

	// 2D Grid START //
	int cols = 0; // create cols of txt tensor field
	int rows = 0; // create rows of txt tensor field
	cout << "before matrix read" << endl;
	workDir = GetCurrentWorkingDir();
	MatrixXd m = readMatrix(workDir + "/matrix.txt", &cols, &rows); // call countMatrix to determine rows/cols count #
	width = cols / 2; // determine width of grid for correct indexing
	height = rows / 2;
	const int dim = width * height; // determine # of dimensions of grid for buffer (string/coefficient etc..) vectors
	
	// create vector for light src's symbolic user-input functions in convenient string format
	std::vector<std::string> userFunctions;
	std::vector<Pair> userPositions;
	// parse input option file
	parse_options(argc, argv, userFunctions, userPositions);

	std::vector<double> initArray(steps, 0.0);
	std::vector<std::vector<double>> glyphBuffer((width)*(height), initArray);
	std::vector<std::vector<double>> glyphParameters(width*height, std::vector<double>(3, 0.0));
	std::vector<std::vector<bool>> signMap(width*height, std::vector<bool>(2, false)); // create a signMap relating normal force signs to singular values (sign_xx, sign_yy, dominant_shear)
	// compute Eigenframes/Superquadrics/Ellipses/Glyphs by calling computeGlyphs w. respective args
	computeGlyphs(glyphBuffer, signMap, glyphParameters);

	// POLAR GRAPHER START //
		//vector to store functions

		// set options
	int Alpha = 200; // set opacity ("Deckung")
	int drawSpeed = 0; // set instant drawing
	int lineWidth = 0; // set line width = 1px
	double thetaMax = 2 * pi;
	double thetaInc = (2 * pi) / steps; // set theta increment
	double rotSpeed = 0.0; // set rotation speed
	sf::Color red = sf::Color(255, 0, 0);
	sf::Color green = sf::Color(0, 255, 0, Alpha);


	// add pre-computed function as string
	/*for (int i = 0; i < dim; i++)
	{
		funcs.push_back(polarPlot(sampleBufferA.at(i), drawSpeed, lineWidth, steps, thetaMax, radres, rotSpeed, red));
		funcsEllipses.push_back(polarPlot(glyphBuffer.at(i), drawSpeed, lineWidth, steps, thetaMax, radres, rotSpeed, green));
	}*/

	polarPlot glyphs(glyphBuffer, drawSpeed, lineWidth, steps, thetaMax, radres, rotSpeed, green);

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
	//for (int i = 0; i < funcs.size(); i++)
	glyphs.init();

	//initialize rendertexture for output image
	sf::RenderTexture out;
	out.create(wSize, wSize);
	sf::RenderTexture outAll;
	outAll.create(wSize, wSize);

	// create flags for blending tensors(ellipses)/light distributions(polar plots)
	bool showOverlay{ true };
	bool showBase{ true };
	// ellipse parameters
	unsigned short quality = 90;
	
	int cellWidth = (wSize / (width)); // check rest (modulo) for x-offset
	int cellHeight = (wSize / (width)); // check division for y offset


	// define ellipses: white origin dots
	sf::ConvexShape tensorFieldLinePos;
	tensorFieldLinePos.setPointCount(quality);
	
	tensorFieldLinePos.setOrigin(0,0);
	tensorFieldLinePos.setPosition(3*wSize / 4, wSize / 2);

	sf::ConvexShape tensorFieldLineNeg;
	tensorFieldLineNeg.setPointCount(quality);
	defineConvexEllipse(&tensorFieldLineNeg, 3.0, 3.0, quality);

	tensorFieldLineNeg.setOrigin(0,0);
	tensorFieldLineNeg.setPosition(3 * wSize / 4, wSize / 2);

	double degPos = glyphParameters.at(width/2+height/2*width).at(2);
	double degNeg = glyphParameters.at(width / 2 + height / 2 * width).at(2);
	double xPos = 3 * wSize / 4; double yPos = wSize / 2 - 20;
	double xNeg = 3 * wSize / 4; double yNeg = wSize / 2 - 20;
	for (int i = 0; i < 500; i++)
	{
		/*if (xPos >= wSize || yPos >= wSize || xPos < 0 || yPos < 0 || xNeg >= wSize || yNeg >= wSize || xNeg < 0 || yNeg < 0)
			continue;*/

		// TensorFieldLinePos

		if (xPos / cellWidth >= width-1 || yPos / cellWidth >= height-1 || xPos < 1 || yPos < 1 || xNeg / cellWidth >= width || yNeg / cellWidth >= height-1 || xNeg < 1 || yNeg < 1)
			continue;
		
		double alphaX = abs(xPos - round(xPos));
		std::vector<double> xPosInterpolantFloor = alphaX * glyphParameters.at(ceil((xPos-cellWidth/2)/cellWidth) + floor((yPos - cellWidth / 2) /cellWidth) * width) + (1 - alphaX)*glyphParameters.at(floor((xPos - cellWidth / 2) / cellWidth) + floor((yPos - cellWidth / 2) / cellWidth) * width);
		std::vector<double> xPosInterpolantCeil = alphaX * glyphParameters.at(ceil((xPos - cellWidth / 2) / cellWidth) + ceil((yPos - cellWidth / 2) / cellWidth) * width) + (1 - alphaX)*glyphParameters.at(floor((xPos - cellWidth / 2) / cellWidth) + ceil((yPos - cellWidth / 2) / cellWidth) * width);
		
		double alphaYCeil = yPos - floor(yPos);
		std::vector<double> interpolant = alphaYCeil * xPosInterpolantCeil + (1 - alphaYCeil)*xPosInterpolantFloor;

		defineConvexEllipse(&tensorFieldLinePos, interpolant.at(0), interpolant.at(1), quality, interpolant.at(2));
		degPos = interpolant.at(2);
		
		//convert polar to cartesian
		double dx = 1.0 * cos(degPos * pi/180.0);
		double dy = -1.0 * sin(degPos * pi / 180.0);
		xPos += dx; yPos += dy;

		tensorFieldLinePos.setPosition(round(xPos),round(yPos));
		window.draw(tensorFieldLinePos);
		out.draw(tensorFieldLinePos);

		// TensorFieldLineNeg
		alphaX = abs(xNeg - round(xNeg));
		xPosInterpolantFloor = alphaX * glyphParameters.at(ceil((xNeg - cellWidth / 2) / cellWidth) + floor((yNeg - cellWidth / 2) / cellWidth) * width) + (1 - alphaX)*glyphParameters.at(floor((xNeg - cellWidth / 2) / cellWidth) + floor((yNeg - cellWidth / 2) / cellWidth) * width);
		xPosInterpolantCeil = alphaX * glyphParameters.at(ceil((xNeg - cellWidth / 2) / cellWidth) + ceil((yNeg - cellWidth / 2) / cellWidth) * width) + (1 - alphaX)*glyphParameters.at(floor((xNeg - cellWidth / 2) / cellWidth) + ceil((yNeg - cellWidth / 2) / cellWidth) * width);

		alphaYCeil = yNeg - floor(yNeg);
		interpolant = alphaYCeil * xPosInterpolantCeil + (1 - alphaYCeil)*xPosInterpolantFloor;

		defineConvexEllipse(&tensorFieldLineNeg, interpolant.at(0), interpolant.at(1), quality, interpolant.at(2));
		degNeg = interpolant.at(2);
		
		//convert polar to cartesian
		dx = 1.0 * cos(degNeg*pi / 180.0);
		dy = -1.0 * sin(degNeg*pi / 180.0);
		xNeg -= dx; yNeg -= dy;
		tensorFieldLineNeg.setPosition(round(xNeg), round(yNeg));
		window.draw(tensorFieldLineNeg);
		out.draw(tensorFieldLineNeg);

	}

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

		/*window.clear(sf::Color::Black);
		outAll.clear(sf::Color::Black);
		out.clear(sf::Color::Black);*/
	

		// define ellipses: white origin dots
		sf::ConvexShape ellipse;
		ellipse.setPointCount(quality);
		defineConvexEllipse(&ellipse, 1.5, 1.5, quality);

		// draw polar function as graph sprite
		for (int i = 0; i < dim; i++)
		{
			int offsetX = (i % width)*(wSize / width) + (wSize / (2 * width)); // check rest (modulo) for x-offset
			int offsetY = (i / width)*(wSize / width) + (wSize / (2 * width)); // check division for y offset

			ellipse.setPosition(offsetX, offsetY); // set ellipse position

			// call animation to continously draw sampled arrays in polar space
			glyphs.animation(i, 1, glyphBuffer); // mode 1 for suppressed test outputs
			sf::Sprite sprE = glyphs.update(); // draw sprites[Kobolde/Elfen] (composited bitmaps/images - general term for objects drawn in the framebuffer)
			//spr.setPosition(wSize / 2, wSize / 2);
			window.draw(ellipse);
			out.draw(ellipse);
			outAll.draw(ellipse);
			if (showOverlay)
				window.draw(sprE);
			outAll.draw(sprE);
		}
		// window.draw(ellipse);  // TEST //
		window.display(); // update window texture
	}
	window.display(); // update window texture

	//write to image file
	out.display(); // update output texture
	outAll.display(); // update output texture
	sf::Texture outTx = out.getTexture();
	sf::Texture outTxAll = outAll.getTexture();
	sf::Image outImg = outTx.copyToImage();
	sf::Image outImgAll = outTxAll.copyToImage();

	outImgAll.saveToFile("all-" + std::to_string(lightSrcPos.jIndex) + "," + std::to_string(lightSrcPos.iIndex) + ".png");
	outImg.saveToFile("intensity-" + std::to_string(lightSrcPos.jIndex) + "," + std::to_string(lightSrcPos.iIndex) + ".png");

	// POLAR GRAPHER END //

	system("PaUsE");

	return EXIT_SUCCESS;
}