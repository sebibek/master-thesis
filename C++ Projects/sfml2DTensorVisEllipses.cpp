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

#define MAXBUFSIZE  ((int) 1e3)

using namespace std;
using namespace Eigen;

//length and width of window
int wSize = 700;
sf::Vector2i windowsize;

//recording frameskip 
int record_frameskip = -1;

//directory to write images to, must exist
std::string record_folder = "frames";

//fullscreen flag
bool fullscreen = false;

//definition of pi
const double pi = M_PI;

std::string functionString; // Prototyping

template <typename T>
T muPlus(T& a, T& b)
{
	T res = a + b;
	return res;
}

// OPERATOR DEFINITIONS //
template <typename T>
std::vector<T> operator+(const std::vector<T>& a, const std::vector<T>& b)
{
	assert(a.size() == b.size());
	
	std::vector<T> result(a.size());

	for(int i = 0; i < a.size(); i++)
		result[i] = a[i] + b[i];
	
	return result;
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

	//iterate function draw
	sf::Sprite& update() {
		graphTx.update(graph);
		graphSpr.rotate(rotation * 180.0/ pi);
		return graphSpr;
	}
	//plots point to pixel(s)
	void plot(int index = 0, int length = 0) {

		//scale
		double scale = 40;
		//convert polar to cartesian
		double x = scale*r * cos(t);
		double y = -1*scale*r * sin(t); // CAVEAT: y-coordinates mirrored at x-axis because of inverse array order

		
		if (x + 1 + wSize / 2 > wSize || x - 1 + wSize / 2 < 0 || y + 1 + wSize / 2 > wSize || y - 1 + wSize / 2 < 0)
			return;

		int offsetX = (index % length)*(wSize/length) + (wSize/(2*length)); // check rest (modulo) for x-offset
		int offsetY = (index / length)*(wSize / length) + (wSize / (2 * length)); // check division for y offset
		
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
			graph.setPixel((xIndex + 1), (yIndex + 1 ), color);
			graph.setPixel((xIndex + 1), (yIndex - 1), color);
			graph.setPixel((xIndex - 1), (yIndex - 1), color);
			graph.setPixel((xIndex - 1), (yIndex + 1), color);
		}

		/*if (line_width == 1 || line_width == 2) {
			graph.setPixel((x + wSize / 2), (y + 1 + wSize / 2), color);
			graph.setPixel((x + wSize / 2), (y - 1 + wSize / 2), color);
			graph.setPixel((x + 1 + wSize / 2), (y + wSize / 2), color);
			graph.setPixel((x - 1 + wSize / 2), (y + wSize / 2), color);
		}

		if (line_width == 2) {
			graph.setPixel((x + 1 + wSize / 2), (y + 1 + wSize / 2), color);
			graph.setPixel((x + 1 + wSize / 2), (y - 1 + wSize / 2), color);
			graph.setPixel((x - 1 + wSize / 2), (y - 1 + wSize / 2), color);
			graph.setPixel((x - 1 + wSize / 2), (y + 1 + wSize / 2), color);
		}*/
	}

	//plots draw_speed points per frame
	void animation(int index = 0, int length = 0) {
		// p.DefineVar("theta", &t);
		for (int i = 0; (i < speed || speed == 0) && t <= tmax; i++) {
			t += tinc;
			r = p.Eval();
			plot(index, length);
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

// file parsers - commented //
//void parse_file(char* filename, std::vector<polarPlot>& funcs) {
//
//	std::ifstream f(filename);
//
//	std::string line;
//
//	//default function attributes
//	sf::Color color = sf::Color::Cyan;
//	std::string func_literal = "100*cos(theta)";
//	int speed = 100;
//	int line_width = 0;
//	double tmax = 50 * pi;
//	double tinc = pi / 3000;
//	double rotation = 0;
//
//	if (f.is_open()) {
//
//		while (std::getline(f, line)) {
//
//			//ignore comments
//			if (line[0] == '#')
//				continue;
//			//write function to vector
//			else if (line == "end") {
//				funcs.push_back(polarPlot(func_literal, speed, line_width, tmax, tinc, rotation, color));
//				//reset default values
//				func_literal = "100*cos(theta)";
//				speed = 100;
//				line_width = 0;
//				color = sf::Color::Cyan;
//				tmax = 50 * pi;
//				tinc = pi / 3000;
//				rotation = 0;
//			}
//			//check for fullscreen
//			else if (line == "-f") {
//				fullscreen = true;
//			}
//			//parse other statements
//			else {
//				//grab keyword
//				std::string::size_type pos;
//				pos = line.find(' ', 0);
//				std::string tag = line.substr(0, pos);
//
//				//get function literal
//				if (tag == "function")
//					func_literal = line.substr(pos + 1);
//				//draw speed
//				else if (tag == "draw_speed") {
//					std::stringstream s;
//					s << line.substr(pos + 1);
//					s >> speed;
//				}
//				//line color
//				else if (tag == "color") {
//					std::stringstream s;
//					s << line.substr(pos + 1);
//					int r = 0, g = 0, b = 0;
//					s >> r >> g >> b;
//					color = sf::Color(r, g, b);
//				}
//				//window/graph size
//				else if (tag == "-s") {
//					std::stringstream s;
//					s << line.substr(pos + 1);
//					s >> wSize;
//				}
//				//theta max
//				else if (tag == "theta_max") {
//					mu::Parser p;
//					p.DefineConst("pi", pi);
//					p.SetExpr(line.substr(pos + 1));
//					tmax = p.Eval();
//				}
//				//theta increment
//				else if (tag == "theta_increment") {
//					mu::Parser p;
//					p.DefineConst("pi", pi);
//					p.SetExpr(line.substr(pos + 1));
//					tinc = p.Eval();
//				}
//				//line width
//				else if (tag == "line_width") {
//					std::stringstream s;
//					s << line.substr(pos + 1);
//					s >> line_width;
//				}
//				else if (tag == "rotation_speed") {
//					mu::Parser p;
//					p.DefineConst("pi", pi);
//					p.SetExpr(line.substr(pos + 1));
//					rotation = p.Eval();
//				}
//				else if (tag == "-r") {
//					std::stringstream s;
//					s << line.substr(pos + 1);
//					s >> record_frameskip >> record_folder;
//				}
//			}
//		}
//
//		f.close();
//	}
//	else
//		std::cerr << filename << " is not a valid filename.\n";
//}
//
//void parse_options(int argc, char* argv[], std::vector<polarPlot>& funcs) {
//
//	int c;
//	std::string frameskip_opt = "-1";
//	std::string size_opt = "-1";
//	std::string directory_opt = "";
//	//	extern char *optarg;
//	//	extern int optind, optopt;
//
//		//using getopts 
//	while ((c = getopt(argc, argv, "fs:r:d:")) != -1) {
//
//		int f = -1, s = -1;
//
//		switch (c) {
//
//		case 'f': {
//			fullscreen = true;
//			break; }
//		case 'r':
//			frameskip_opt.assign(optarg);
//			std::istringstream(frameskip_opt) >> f;
//			if (f <= -1) {
//				std::cerr << "invalid argument \'" << frameskip_opt << "\' to option -r\n";
//				optind--;
//			}
//			record_frameskip = f > 0 ? f : 0;
//			break;
//		case 's':
//			size_opt.assign(optarg);
//			std::istringstream(size_opt) >> s;
//			if (s <= 0) {
//				std::cerr << "invalid argument \'" << size_opt << "\' to option -s\n";
//				optind--;
//			}
//			wSize = s > 0 ? s : wSize;
//			break;
//		case 'd':
//			directory_opt.assign(optarg);
//			record_folder = directory_opt;
//			break;
//		case ':':
//			switch (optopt) {
//			case 's':
//				std::cout << "option -s requires argument, using default size 700\n";
//				break;
//			case 'r':
//				std::cerr << "using default frameskip of 0 for option -r\n";
//				break;
//			case 'd':
//				std::cerr << "option -d requires argument, disabling recording.\n";
//				record_frameskip = -1;
//				break;
//			}
//			break;
//		case '?':
//			std::cerr << "Unrecognized option: '-" << optopt << "\n";
//			break;
//		}
//
//
//	}
//	for (int i = optind; i < argc; i++)
//		parse_file(argv[i], funcs);
//}

// mathematical functions
double circleFunction(double x)
{
	double os = 1.0/sqrt(M_PI);
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

	expression_t expression;
	expression.register_symbol_table(symbol_table);

	parser_t parser;
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
	for (int i = 0; i < rows; i++)
		for (int j = 0; j < cols; j++)
			result(i, j) = buff[cols*i + j];

	return result;
};

int main(int argc, char* argv[])
{
	// SVD/EIGENVALUE DECOMPOSITION START //

	// 2D Grid START //
	int* cols = new int(0); // create ptr to cols
	int* rows = new int(0); // create ptr to rows
	MatrixXd m = readMatrix("../matrix.txt", cols, rows);
	const int length = *cols/2;

	std::cout << "Here is the matrix m:" << std::endl << m << std::endl;
	JacobiSVD<MatrixXd> svd(m, ComputeThinU | ComputeThinV);
	std::cout << "Its singular values are:" << std::endl << svd.singularValues() << std::endl;
	std::cout << "Its left singular vectors are the columns of the thin U matrix:" << std::endl << svd.matrixU() << std::endl;
	std::cout << "Its right singular vectors are the columns of the thin V matrix:" << std::endl << svd.matrixV() << std::endl;

	std::cout << "Blocks from left to the right" << endl;
	for (int i = 0; i < *rows / 2; i++) // use rows/2 and cols/2 because of 2D matrices!!!
		for (int j = 0; j < *cols / 2; j++)
			std::cout << m.block<2, 2>(2 * i, 2 * j) << endl << endl;

	// block the read-out matrices into a list of MatrixXd types (hardcoded 2x2 read-out)
	int dim = *rows / 2 * *cols / 2;
	std::vector<MatrixXd> matrixList(dim, MatrixXd::Ones(2, 2));
	for (int i = 0; i < *rows / 2; i++)
		for (int j = 0; j < *cols / 2; j++)
		{
			MatrixXd a = m.block<2, 2>(2 * i, 2 * j);
			matrixList.at(j + i * (*cols / 2)) = (a);
		}

	// compute the SVD of all matrices in matrixList into svdList
	std::vector<JacobiSVD<MatrixXd>> svdList(dim, JacobiSVD<MatrixXd>(matrixList.at(0), ComputeThinU | ComputeThinV));
	for (int i = 0; i < dim; i++)
	{
		// SVD-based
		JacobiSVD<MatrixXd> svd(matrixList.at(i), ComputeThinU | ComputeThinV);
		svdList.at(i) = svd;
		MatrixXd a = matrixList.at(i);

		std::cout << "Here is the matrix A:" << std::endl << a << std::endl;
		std::cout << "Its singular values are:" << std::endl << svdList.at(i).singularValues() << std::endl;
		std::cout << "Its left singular vectors are the columns of the thin U matrix:" << std::endl << svdList.at(i).matrixU() << std::endl;
		std::cout << "Its right singular vectors are the columns of the thin V matrix:" << std::endl << svdList.at(i).matrixV() << std::endl;


		// Eigenvector-based
		/*EigenSolver<MatrixXd> es(a, true);
		std::cout << "Here is the matrix A:" << std::endl << a << std::endl;
		std::cout << "Its eigen values are:" << std::endl << es.eigenvalues() << std::endl;
		std::cout << "Its first eigenvector is:" << std::endl << es.eigenvectors().col(0).real() << std::endl;
		std::cout << "Its second eigenvector is:" << std::endl << es.eigenvectors().col(1).real() << std::endl;*/
	}

	std::string* functionStr = new std::string[length*length];

	for (int j = 0; j < *rows / 2; j++)
	{
		for (int i = 0; i < *cols / 2; i++)
		{
			double y1 = svdList.at(i + j * (*cols / 2)).matrixU().col(0)[1]; // use x - coordinate of both semi-axes 
			double x1 = svdList.at(i + j * (*cols / 2)).matrixU().col(0)[0]; // use x - coordinate of both semi-axes "sigma_xx"
			double y2 = svdList.at(i + j * (*cols / 2)).matrixU().col(1)[1]; // use x - coordinate of both semi-axes 
			double x2 = svdList.at(i + j * (*cols / 2)).matrixU().col(1)[0]; // use x - coordinate of both semi-axes "sigma_xx"
			double xx = matrixList.at(i + j * (*cols / 2)).row(0)[0]; // "sigma_xx"
			double yy = matrixList.at(i + j * (*cols / 2)).row(1)[1]; // "sigma_yy"
			double deg1 = atan2(y1, x1) * 180.0 / M_PI; // use vector atan2 to get rotational angle (phase) of both basis vectors in [-180°,180°]
			double deg2 = atan2(y2, x2) * 180.0 / M_PI; // use vector atan2 to get rotational angle (phase) of both basis vectors [-180°,180°]

			// shift (normalize) degs from [-180°,180°] into the interval [0°,360°] - "circular value permutation"
			if (deg1 < 0)
				deg1 = 360 + deg1;
			if (deg2 < 0)
				deg2 = 360 + deg2; // deg1 = 270, deg2 = 0

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

			std::cout << "deg: " << deg << std::endl;

			double sv1 = svdList.at(i + j * (*cols / 2)).singularValues()[0];
			double sv2 = svdList.at(i + j * (*cols / 2)).singularValues()[1];
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

			std::string aSquaredString = std::to_string(a*a);
			std::string bSquaredString = std::to_string(b*b);
			std::string dotString = std::to_string(dot);
			std::string radString = std::to_string(deg*(M_PI/180.0));

			std::string ellipse = dotString + "/sqrt(" + bSquaredString + "*cos(theta-" + radString + ")*cos(theta-" + radString + ")+" + aSquaredString + "*sin(theta-" + radString + ")*sin(theta-" + radString + "))";
			
			functionStr[i + j * (*cols / 2)] = ellipse;
		}
	}

	// POLAR GRAPHER START //
	//vector to store functions
	std::vector<polarPlot> funcs;

	// set options
	//std::string functionStr = "200 * sin(4 * theta / 9)";
	int drawSpeed = 0; // set instant drawing
	int lineWidth = 1;
	double thetaMax = 2 * pi;
	double thetaInc = pi / 300; // set theta increment
	double rotSpeed = 0; // set rotation speed
	sf::Color red = sf::Color(255, 0, 0);
	
	// add pre-computed function as string
	for (int i = 0; i < length*length; i++)
		funcs.push_back(polarPlot(functionStr[i], drawSpeed, lineWidth, thetaMax, thetaInc, rotSpeed, red));

	//parse input file
	//parse_options(argc, argv, funcs);

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
		funcs[i].init();

	//initialize rendertexture for output image
	sf::RenderTexture out;
	out.create(wSize, wSize);

	//main loop
	while (window.isOpen())
	{
		sf::Event event;
		while (window.pollEvent(event))
		{
			if (event.type == sf::Event::Closed)
				window.close();
			if (event.key.code == sf::Keyboard::Q)
				window.close();
		}

		window.clear(sf::Color::Black);
		out.clear(sf::Color::Black);
		for (unsigned i = 0; i < funcs.size(); i++) {
			funcs[i].animation(i, length);
			sf::Sprite spr = funcs[i].update();
			window.draw(spr);
			spr.setPosition(wSize / 2, wSize / 2);
			out.draw(spr);
		}
		window.display();

		record(out, record_frameskip, record_folder);
	}

	//write to image file
	out.display();
	sf::Texture outTx = out.getTexture();
	sf::Image outImg = outTx.copyToImage();
	outImg.saveToFile("out.png");
	delete[] functionStr;
	system("pause");
	
	return 0;
}