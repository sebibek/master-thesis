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

// define function string
std::string functionStr = "cos(theta) + sin(3*theta)";

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

std::list<std::list<double>> calcCoeff(double(*f)(double x)) // for 2Pi periodic functions
{
	double a0 = 1 /(2*M_PI) * integral(f, -M_PI, M_PI, 50);
	std::cout << "a0: " << a0 << endl;
	bool finish = false;
	std::list<std::list<double>> coeff;
	std::list<double> an;
	std::list<double> bn;
	double thresh = 0.1;

	int n = 1;
	while(!finish)
	{
		an.push_back(1 / M_PI * aIntegral(f, -M_PI, M_PI, 50, n)); // calculate Integrals for the coeff
		bn.push_back(1 / M_PI * bIntegral(f, -M_PI, M_PI, 50, n)); // calculate Integrals for the coeff
	
		n++;
		std::cout << "an: " << an.back() << endl;
		std::cout << "bn: " << bn.back() << endl;
		if (bn.back() < thresh && an.back() < thresh && n > 10) // if coeff < thresh, abort!
			finish = true;
	}
	std::cout << "n: " << n << endl;
	std::cout << "an_size: " << an.size() << endl;
	coeff.push_back(an);
	coeff.push_back(bn);

	return coeff;
}

void defineConvexEllipse(sf::ConvexShape* ellipse, double radius_x, double radius_y, unsigned short quality)
{
	for (unsigned short i = 0; i < quality; ++i)
	{
		float rad = (360 / quality * i) / (360 / M_PI / 2);
		float x = cos(rad)*radius_x;
		float y = sin(rad)*radius_y;

		ellipse[0].setPoint(i, sf::Vector2f(x, y));
	}
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
	//declare renderwindow
	sf::RenderWindow window;
	if (fullscreen == false)
		window.create(sf::VideoMode(wSize, wSize, 32), "Polar Grapher", sf::Style::Default);
	else
		window.create(sf::VideoMode::getDesktopMode(), "Polar Grapher", sf::Style::Fullscreen);
	window.setFramerateLimit(60);
	windowsize.x = window.getSize().x;
	windowsize.y = window.getSize().y;

	////initialize rendertexture for output image
	sf::RenderTexture out;
	out.create(wSize, wSize);

	// 2D Grid START //
	int* cols = new int(0); // create ptr to cols
	int* rows = new int(0); // create ptr to rows
	MatrixXd m = readMatrix("../matrix.txt", cols, rows);

	std::cout << "Here is the matrix m:" << std::endl << m << std::endl;
	JacobiSVD<MatrixXd> svd(m, ComputeThinU | ComputeThinV);
	std::cout << "Its singular values are:" << std::endl << svd.singularValues() << std::endl;
	std::cout << "Its left singular vectors are the columns of the thin U matrix:" << std::endl << svd.matrixU() << std::endl;
	std::cout << "Its right singular vectors are the columns of the thin V matrix:" << std::endl << svd.matrixV() << std::endl;

	std::cout << "Blocks from left to the right" << endl;
	for (int i = 0; i < *rows / 2; i++)
		for (int j = 0; j < *cols / 2; j++)
			std::cout << m.block<2, 2>(2 * i, 2 * j) << endl << endl;

	// block the read-out matrices into a list of MatrixXd types (hardcoded 2x2 read-out)
	int dim = *rows / 2 * *cols / 2;
	std::vector<MatrixXd> matrixList(dim, MatrixXd::Ones(2,2));
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
		JacobiSVD<MatrixXd> svd(matrixList.at(i), ComputeThinU | ComputeThinV);
		svdList.at(i) = svd;
		MatrixXd a = matrixList.at(i);

		std::cout << "Here is the matrix A:" << std::endl << a << std::endl;
		std::cout << "Its singular values are:" << std::endl << svdList.at(i).singularValues() << std::endl;
		std::cout << "Its left singular vectors are the columns of the thin U matrix:" << std::endl << svdList.at(i).matrixU() << std::endl;
		std::cout << "Its right singular vectors are the columns of the thin V matrix:" << std::endl << svdList.at(i).matrixV() << std::endl;

		/*EigenSolver<MatrixXd> es(a, true);
		std::cout << "Here is the matrix A:" << std::endl << a << std::endl;
		std::cout << "Its eigen values are:" << std::endl << es.eigenvalues() << std::endl;
		std::cout << "Its first eigenvector is:" << std::endl << es.eigenvectors().col(0).real() << std::endl;
		std::cout << "Its second eigenvector is:" << std::endl << es.eigenvectors().col(1).real() << std::endl;*/
	}

	// ellipse parameters
	unsigned short quality = 70;
	// define ellipses as convex shapes (with semi-axes u scaled to the specific singular values
	std::vector<sf::ConvexShape> ellipseList;

	for (int j = 0; j < *rows / 2; j++)
	{
		for (int i = 0; i < *cols / 2; i++)
		{
			sf::ConvexShape ellipse;
			ellipse.setPointCount(quality);
			double y1 = svdList.at(i + j * (*cols / 2)).matrixU().col(0)[1]; // use x - coordinate of both semi-axes 
			double x1 = svdList.at(i + j * (*cols / 2)).matrixU().col(0)[0]; // use x - coordinate of both semi-axes "sigma_xx"
			double y2 = svdList.at(i + j * (*cols / 2)).matrixU().col(1)[1]; // use x - coordinate of both semi-axes 
			double x2 = svdList.at(i + j * (*cols / 2)).matrixU().col(1)[0]; // use x - coordinate of both semi-axes "sigma_xx"
			double xx = matrixList.at(i + j * (*cols / 2)).row(0)[0]; // "sigma_xx"
			double yy = matrixList.at(i + j * (*cols / 2)).row(1)[1]; // "sigma_yy"
			double deg1 = atan2(y1,x1) * 180.0 / M_PI; // use vector atan2 to get rotational angle (phase) of both basis vectors in [-180°,180°]
			double deg2 = atan2(y2,x2) * 180.0 / M_PI; // use vector atan2 to get rotational angle (phase) of both basis vectors [-180°,180°]

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

			if ((yy < 0 && xx < 0) || (yy >= 0 && xx >= 0)) // if not mirroring at one single axis..
			{
				if (abs(yy) < abs(xx)) // if anisotropic x-scaling, scale x
					defineConvexEllipse(&ellipse, 50 * svdList.at(i + j * (*cols / 2)).singularValues()[0], 50 * svdList.at(i + j * (*cols / 2)).singularValues()[1], quality);
				else // if anisotropic y-scaling, scale y
					defineConvexEllipse(&ellipse, 50 * svdList.at(i + j * (*cols / 2)).singularValues()[1], 50 * svdList.at(i + j * (*cols / 2)).singularValues()[0], quality);
			}
			else // if mirroring at one single axis.. swap!
			{
				if (abs(yy) < abs(xx)) // if anisotropic x-scaling, scale x
					defineConvexEllipse(&ellipse, 50 * svdList.at(i + j * (*cols / 2)).singularValues()[1], 50 * svdList.at(i + j * (*cols / 2)).singularValues()[0], quality);
				else // if anisotropic y-scaling, scale y
					defineConvexEllipse(&ellipse, 50 * svdList.at(i + j * (*cols / 2)).singularValues()[0], 50 * svdList.at(i + j * (*cols / 2)).singularValues()[1], quality);
			}

			ellipse.setPosition(130 * i + 90, 130 * j + 90); // set ellipse position
			std::cout << "deg: " << deg << std::endl;
			ellipse.setRotation(-1*deg); // use vector (quadrant) atan2 to determine rotation of ellipse
			ellipseList.push_back(ellipse);
		}
	}

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
		for (int i = 0; i < ellipseList.size(); i++)
		{
			window.draw(ellipseList.at(i));
		}
		
		window.display();

	}

	system("pause");
	return 0;
}