#include <iostream>
#include <vector>
#include <stdexcept>
#include <Eigen/Eigen>
using namespace Eigen;

//using namespace std;


class my_matrix {
	std::vector<std::vector<double>> m;
	int width;
	int height;
public:
	my_matrix(unsigned int x, unsigned int y)
	{
		m.resize(x, std::vector<double>(y, 0.0));
		width = x;
		height = y;
	}
	int getWidth() {return width;}
	int getHeight() {return height;}
	class matrix_row {
		std::vector<double>& row;
	public:
		matrix_row(std::vector<double>& r) : row(r) // row gets assigned with adress r
		{
		}
		double& operator[](unsigned int y) // inner operator
		{
			return row.at(y); // This is where the doubles come from !!!
		}
	};
	matrix_row operator[](unsigned int x) // outer operator
	{
		return matrix_row(m.at(x));
	}
};


struct complex { double r, i; bool real; };
template <typename T>
struct pair { T p1, p2; };

pair<complex> midnightFormula(double a, double b, double c)
{
	pair<complex> result = { 0 };

	if (a < 0.000001)    // ==0
	{
		if (b > 0.000001)  // !=0
			result.p1.r = result.p2.r = -c / b;
		else
			if (c > 0.00001) throw std::exception("no solutions");
		return result;
	}

	double delta = b * b - 4 * a*c;
	if (delta >= 0)
	{
		result.p1.r = (-b - sqrt(delta)) / 2 / a;
		result.p2.r = (-b + sqrt(delta)) / 2 / a;
		result.p1.real = result.p2.real = true;
	}
	else
	{
		result.p1.r = result.p2.r = -b / 2 / a;
		result.p1.i = sqrt(-delta) / 2 / a;
		result.p2.i = -sqrt(-delta) / 2 / a;
	}

	return result;
}

pair<complex> Eigenvalues(my_matrix matrix)
{
	double lambda1 = 0.0;
	double lambda2 = 0.0;

	double a = 1;
	double b = -1 * (matrix[0][0] + matrix[1][1]);
	double c = matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];

	pair<complex> result = midnightFormula(a, b, c);

	return result;
}

std::vector<double> Eigenvectors(my_matrix matrix, pair<complex> eigenvalues)
{
	int a = 1;
	return std::vector<double>(0);
}


int main()
{
	// Example usage
	//int width = 2; // --> i/x
	//int height = 2; // --> j/y
	//my_matrix matrix(width, height);
	//matrix[0][0] = 2;
	//matrix[0][1] = 1;
	//matrix[1][0] = 1;
	//matrix[1][1] = 2;
	//for (int j = 0; j < matrix.getHeight(); j++)
	//{
	//	for (int i = 0; i < matrix.getWidth(); i++)
	//		std::cout << "mm[" << j << "][" << i << "]: " << matrix[j][i] << "	";
	//	std::cout << std::endl;
	//}
	//pair<complex> res = Eigenvalues(matrix);
	//std::cout << "p1: " << res.p1.r << "   " << "p2: " << res.p2.r << std::endl;
	Vector2d a(1, 1);
	MatrixXd m = Eigen::MatrixXd::Ones(2,2);
	m.row(0) << 0.5, 1;
	m.row(1) << -0.5, 1;
	std::cout << "Here is the matrix m:" << std::endl << m << std::endl;
	JacobiSVD<MatrixXd> svd(m, ComputeThinU | ComputeThinV);
	std::cout << "Its singular values are:" << std::endl << svd.singularValues() << std::endl;
	std::cout << "Its left singular vectors are the columns of the thin U matrix:" << std::endl << svd.matrixU() << std::endl;
	std::cout << "Its right singular vectors are the columns of the thin V matrix:" << std::endl << svd.matrixV() << std::endl;
	Vector2d rhs(1, 0);
	std::cout << "Now consider this rhs vector:" << std::endl << rhs << std::endl;
	std::cout << "A least-squares solution of m*x = rhs is:" << std::endl << svd.solve(rhs) << std::endl;

	system("pause");
	return 0;
}