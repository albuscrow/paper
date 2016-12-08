/*
 * 可以生成有邻接关系的茶壶obj
 */
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

using namespace std;

struct Face
{
	int v0, v1, v2;
	int vt0, vt1, vt2;
	int n0, n1, n2;
	Face(int vv0, int vv1, int vv2,
		 int tt0, int tt1, int tt2,
		 int nn0, int nn1, int nn2)
	{
		v0 = vv0; v1 = vv1; v2 = vv2;
		vt0 = tt0; vt1 = tt1; vt2 = tt2;
		n0 = nn0; n1 = nn1; n2 = nn2;
	}
};

ostream &operator<<(ostream &os, const Face &f)
{
	return os << f.v0 << "/" << f.vt0 << "/" << f.n0 << " "
			  << f.v1 << "/" << f.vt1 << "/" << f.n1 << " "
			  << f.v2 << "/" << f.vt2 << "/" << f.n2;
	//return os << f.v0 << "//" << f.n0 << " "
			  //<< f.v1 << "//" << f.n1 << " "
			  //<< f.v2 << "//" << f.n2;
}

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
	void normalize()
	{
		double length = sqrt(x * x + y * y + z * z);
		x /= length;
		y /= length;
		z /= length;
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

double dist(const Point3 &v0, const Point3 &v1)
{
	double dx = v0.x - v1.x;
	double dy = v0.y - v1.y;
	double dz = v0.z - v1.z;
	return sqrt(dx * dx + dy * dy + dz * dz);
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

	const double ZERO = 0.000001;

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

	for (vector<vector<Point3> >::size_type i = 0; i < ctrl_points.size(); ++i)
	{
		for (vector<Point3>::size_type j = 0; j < ctrl_points[i].size(); ++j)
		{
			cout << ctrl_points[i][j].x << ", "
				 << ctrl_points[i][j].y << ", "
				 << ctrl_points[i][j].z << ", \n";
		}
		cout << endl;
	}
	
	ofstream fout("teapot.obj");

	//const int u_seg = 3, v_seg = 3;

	#define XIAOCHU_CHONGFU
	//const int u_seg = 5, v_seg = 5;
	//const int u_seg = 1, v_seg = 1;

	const int u_seg = 10, v_seg = 10;
	//const int u_seg = 30, v_seg = 30;
	//const int u_seg = 100, v_seg = 100;
	//const int u_seg = 400, v_seg = 400;

	vector<Point3> vertex_vector;	// 全部顶点
	vector<Point3> texture_vector;	// 全部顶点
	vector<Point3> normal_vector;	// 全部法向
	vector<Face> face_vector;		// 全部面片
	vector<pair<Point3, int> > border_vertex_vector;	// 全部边界顶点及其对应编号
	vector<pair<Point3, int> > border_normal_vector;	// 全部边界法向及其对应编号
	for (int p = 0; p < patch_number; ++p)
	{
		cout << "p = " << p << endl;
		int table_v[u_seg + 1][v_seg + 1];	// 从顶点在面片内部的二维编号到全局编号的映射表
		int table_t[u_seg + 1][v_seg + 1];	// 从纹理在面片内部的二维编号到全局编号的映射表
		int table_n[u_seg + 1][v_seg + 1];	// 从法向在面片内部的二维编号到全局编号的映射表
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
				normal.normalize();
				Point3 texture((double)i / u_seg, (double)j / v_seg, p + 0.4);
				if (p < 4 && i == 0)					// 最顶端，四边形退化到三角形
				{
					normal = Point3(0, 0, 1);
					//texture.x = texture.y = 0.0;
				}
				else if (p < 8 && i == 0)				// 最底端，四边形退化到三角形
				{
					normal = Point3(0, 0, -1);
					//texture.x = texture.y = 0.0;
				}

#ifdef XIAOCHU_CHONGFU
				/*-------------------------------------------*/
				if (i == 0 || i == u_seg || j == 0 || j == v_seg)	// 仅有在边界的时候才需要检测重复
				{
					// 找重复顶点
					int existing_index = -1;
					//for (vector<Point3>::size_type ii = 0; ii < vertex_vector.size(); ++ii)
					for (vector<Point3>::size_type ii = 0; ii < border_vertex_vector.size(); ++ii)
					{
						//if (dist(result, vertex_vector[ii]) < ZERO) // 顶点之前出现过
						if (dist(result, border_vertex_vector[ii].first) < ZERO) // 顶点之前出现过
						{
							//table_v[i][j] = existing_index = ii;
							table_v[i][j] = existing_index = border_vertex_vector[ii].second;
							break;
						}
					}
					if (existing_index == -1)	// 顶点之前未出现过
					{
						vertex_vector.push_back(result);
						table_v[i][j] = vertex_vector.size() - 1;
						border_vertex_vector.push_back(pair<Point3, int>(result, table_v[i][j]));
					}

					// 找重复法向
					existing_index = -1;
					//for (vector<Point3>::size_type ii = 0; ii < normal_vector.size(); ++ii)
					for (vector<Point3>::size_type ii = 0; ii < border_normal_vector.size(); ++ii)
					{
						//if (dist(normal, normal_vector[ii]) < ZERO)	// 法向之前出现过
						if (dist(normal, border_normal_vector[ii].first) < ZERO)	// 法向之前出现过
						{
							//table_n[i][j] = existing_index = ii;
							table_n[i][j] = existing_index = border_normal_vector[ii].second;
							break;
						}
					}
					if (existing_index == -1)	// 法向之前未出现过
					{
						normal_vector.push_back(normal);
						table_n[i][j] = normal_vector.size() - 1;
						border_normal_vector.push_back(pair<Point3, int>(normal, table_n[i][j]));
					}

					texture_vector.push_back(texture);
					table_t[i][j] = texture_vector.size() - 1;
				}
				else
#endif
				{
					vertex_vector.push_back(result);
					table_v[i][j] = vertex_vector.size() - 1;
					normal_vector.push_back(normal);
					table_n[i][j] = normal_vector.size() - 1;
					texture_vector.push_back(texture);
					table_t[i][j] = texture_vector.size() - 1;
				}
			}
		}

		for (int i = 0; i < u_seg; ++i)
		{
			for (int j = 0; j < v_seg; ++j)
			{
				int upperleft_v = table_v[i][j] + 1;
				int upperright_v = table_v[i][j + 1] + 1;
				int lowerleft_v = table_v[i + 1][j] + 1;
				int lowerright_v = table_v[i + 1][j + 1] + 1;

				int upperleft_n = table_n[i][j] + 1;
				int upperright_n = table_n[i][j + 1] + 1;
				int lowerleft_n = table_n[i + 1][j] + 1;
				int lowerright_n = table_n[i + 1][j + 1] + 1;

				int upperleft_t = table_t[i][j] + 1;
				int upperright_t = table_t[i][j + 1] + 1;
				int lowerleft_t = table_t[i + 1][j] + 1;
				int lowerright_t = table_t[i + 1][j + 1] + 1;
				//if ((p / 4 == 5 || p / 4 == 7) && i == 0)	// 四边形退化到三角形
				if (p < 8 && i == 0)	// 四边形退化到三角形
				{
#ifdef COUNTER_CLOCK
					fout << "f " << upperright << "//" << upperright << " "
								 << lowerleft << "//" << lowerleft << " "
								 << lowerright << "//" << lowerright << endl;
#else
					//face_vector.push_back(Face(upperright_v, lowerright_v, lowerleft_v,
											   //upperright_t, lowerright_t, lowerleft_t,
											   //upperright_n, lowerright_n, lowerleft_n));
					face_vector.push_back(Face(upperleft_v, lowerright_v, lowerleft_v,
											   upperleft_t, lowerright_t, lowerleft_t,
											   upperleft_n, lowerright_n, lowerleft_n));
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
					face_vector.push_back(Face(upperleft_v, lowerright_v, lowerleft_v,
											   upperleft_t, lowerright_t, lowerleft_t,
											   upperleft_n, lowerright_n, lowerleft_n));
					face_vector.push_back(Face(upperleft_v, upperright_v, lowerright_v,
											   upperleft_t, upperright_t, lowerright_t,
											   upperleft_n, upperright_n, lowerright_n));
					//face_vector.push_back(Face(upperleft_v, upperright_v, lowerleft_v,
											   //upperleft_t, upperright_t, lowerleft_t,
											   //upperleft_n, upperright_n, lowerleft_n));
					//face_vector.push_back(Face(upperright_v, lowerright_v, lowerleft_v,
											   //upperright_t, lowerright_t, lowerleft_t,
											   //upperright_n, lowerright_n, lowerleft_n));
#endif
				}
			}
		}
	}

	for (vector<Point3>::iterator it = vertex_vector.begin(); it != vertex_vector.end(); ++it)
	{
		fout << "v " << *it << endl;
	}
	for (vector<Point3>::iterator it = texture_vector.begin(); it != texture_vector.end(); ++it)
	{
		fout << "vt " << *it << endl;
	}
	for (vector<Point3>::iterator it = normal_vector.begin(); it != normal_vector.end(); ++it)
		fout << "vn " << *it << endl;
	for (vector<Face>::iterator it = face_vector.begin(); it != face_vector.end(); ++it)
		fout << "f " << *it << endl;

	fout.close();
}
