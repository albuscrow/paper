/*
 * 可以生成无邻接关系的茶壶obj
 */
#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

struct Point3
{
	double x, y, z;
	Point3(double xx, double yy, double zz)
	{
		x = xx;
		y = yy;
		z = zz;
	}
	Point3 operator+=(const Point3 &p)
	{
		x += p.x;
		y += p.y;
		z += p.z;
		return *this;
	}
};

const Point3 operator-(const Point3 &p0, const Point3 &p1)
{
	return Point3(p0.x - p1.x, p0.y - p1.y, p0.z - p1.z);
}

const Point3 operator*(double v, const Point3 &p)
{
	return Point3(p.x * v, p.y * v, p.z * v);
}

const Point3 operator*(const Point3 &p, double v)
{
	return Point3(p.x * v, p.y * v, p.z * v);
}

ostream &operator<<(ostream &os, const Point3 &p)
{
	return os << p.x << " " << p.y << " " << p.z;
}

const Point3 cross(const Point3 &v0, const Point3 &v1)
{
	double x = v0.y * v1.z - v0.z * v1.y;
	double y = v0.z * v1.x - v0.x * v1.z;
	double z = v0.x * v1.y - v0.y * v1.x;

	return Point3(x, y, z);
}

/*---------------------------------------------------------------------*/

int factorial[] = {1, 1, 2, 6};
double factorial_1[] = {1, 1, 0.5, 0.166666666666666666};

double power(double a, int x)
{
	double result = 1;
	for (int i = 0; i < x; ++i)
		result *= a;
	return result;
}

double B(double t, int n, int i)
{
	
	return factorial[n] * factorial_1[i] * factorial_1[n - i]
			* power(t, i) * power(1 - t, n - i);
}

/*---------------------------------------------------------------------*/

int main()
{
//#define COUNTER_CLOCK
	ifstream fin("/home/cym/program/OBJ/utah_teapot/t.bpt");	// 不需要define cc
	//ifstream fin("/home/cym/program/OBJ/utah_teapot/teapotCGA.bpt");

	//ifstream fin("/home/cym/program/OBJ/utah_teapot/t2.bpt");	// 需要define cc
	//ifstream fin("/home/cym/program/OBJ/utah_teapot/teapotrim.bpt");

	vector<vector<Point3> > ctrl_points;
	int patch_number;
	fin >> patch_number;
	for (int p = 0; p < patch_number; ++p)
	{
		int u_deg, v_deg;
		fin >> u_deg >> v_deg;
		double vx, vy, vz;
		vector<Point3> temp_vector;
		for (int i = 0; i < u_deg + 1; ++i)
		{
			for (int j = 0; j < v_deg + 1; ++j)
			{
				fin >> vx >> vy >> vz;
				temp_vector.push_back(Point3(vx, vy, vz));
			}
		}
		ctrl_points.push_back(temp_vector);
	}
	fin.close();
	
	ofstream fout("teapot.obj");
	//const int u_seg = 3, v_seg = 3;
	const int u_seg = 5, v_seg = 5;
	//const int u_seg = 30, v_seg = 30;
	//const int u_seg = 100, v_seg = 100;
	//int base = 0;
	for (int p = 0; p < patch_number; ++p)
	{
		for (int i = 0; i < u_seg + 1; ++i)
		{
			double u = static_cast<double>(i) / u_seg;
			for (int j = 0; j < v_seg + 1; ++j)
			{
				double v = static_cast<double>(j) / v_seg;
				// 计算(u, v)点的值
				Point3 result(0, 0, 0);
				for (int ii = 0; ii < 4; ++ii)
				{
					Point3 temp_ctrlpoint(0, 0, 0);
					for (int jj = 0; jj < 4; ++jj)
					{
						temp_ctrlpoint += B(v, 3, jj) * ctrl_points[p][ii * 4 + jj];
					}
					result += B(u, 3, ii) * temp_ctrlpoint;
				}
				fout << "v " << result << endl;

				// 计算(u, v)点u方向导矢的值
				Point3 result_u(0, 0, 0);
				for (int jj = 0; jj < 4; ++jj)
				{
					Point3 temp_ctrlpoint(0, 0, 0);
					for (int ii = 0; ii < 3; ++ii)
					{
						temp_ctrlpoint += B(u, 2, ii)
				* (ctrl_points[p][ii * 4 + jj] - ctrl_points[p][(ii + 1) * 4 + jj]);
					}
					result_u += B(v, 3, jj) * temp_ctrlpoint;
				}

				// 计算(u, v)点v方向导矢的值
				Point3 result_v(0, 0, 0);
				for (int ii = 0; ii < 4; ++ii)
				{
					Point3 temp_ctrlpoint(0, 0, 0);
					for (int jj = 0; jj < 3; ++jj)
					{
						temp_ctrlpoint += B(v, 2, jj)
				* (ctrl_points[p][ii * 4 + jj] - ctrl_points[p][ii * 4 + jj + 1]);
					}
					result_v += B(u, 3, ii) * temp_ctrlpoint;
				}
				//Point3 normal = cross(result_u, result_v);
				Point3 normal = cross(result_v, result_u);
				if (p / 4 == 5 && i == 0)	// 最顶端，四边形退化到三角形
					fout << "vn 0 0 1" << endl;
				else if (p / 4 == 7 && i == 0)	// 最底端，四边形退化到三角形
					fout << "vn 0 0 -1" << endl;
				else
					fout << "vn " << normal << endl;
			}
		}
		int base = p * ((u_seg + 1) * (v_seg + 1)) + 1;
		for (int i = 0; i < u_seg; ++i)
		{
			for (int j = 0; j < v_seg; ++j)
			{
				int upperleft = base + i * (v_seg + 1) + j;
				int upperright = base + i * (v_seg + 1) + j + 1;
				int lowerleft = base + (i + 1) * (v_seg + 1) + j;
				int lowerright = base + (i + 1) * (v_seg + 1) + j + 1;
				if ((p / 4 == 5 || p / 4 == 7) && i == 0)	// 四边形退化到三角形
				{
#ifdef COUNTER_CLOCK
					fout << "f " << upperright << "//" << upperright << " "
								 << lowerleft << "//" << lowerleft << " "
								 << lowerright << "//" << lowerright << endl;
#else
					fout << "f " << upperright << "//" << upperright << " "
								 << lowerright << "//" << lowerright << " "
								 << lowerleft << "//" << lowerleft << endl;
#endif
				}
				else
				{
#ifdef COUNTER_CLOCK
					fout << "f " << upperright << "//" << upperright << " "
								 << upperleft << "//" << upperleft << " "
								 << lowerleft << "//" << lowerleft << endl;
					fout << "f " << upperright << "//" << upperright << " "
								 << lowerleft << "//" << lowerleft << " "
								 << lowerright << "//" << lowerright << endl;
#else
					fout << "f " << upperleft << "//" << upperleft << " "
								 << upperright << "//" << upperright << " "
								 << lowerleft << "//" << lowerleft << endl;
					fout << "f " << upperright << "//" << upperright << " "
								 << lowerright << "//" << lowerright << " "
								 << lowerleft << "//" << lowerleft << endl;
#endif
				}
			}
		}
	}

	fout.close();
}
