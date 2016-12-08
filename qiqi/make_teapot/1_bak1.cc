/*
 * 直接把控制网当作顶点绘制
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
};

//#define COUNTER_CLOCK
#ifdef COUNTER_CLOCK
int table[18 * 3] = 
{
	0, 4, 1,
	1, 4, 5,
	1, 5, 2,
	2, 5, 6,
	2, 6, 3,
	3, 6, 7,

	4, 8, 5,
	5, 8, 9,
	5, 9, 6,
	6, 9, 10,
	6, 10, 7,
	7, 10, 11,

	8, 12, 9,
	9, 12, 13,
	9, 13, 10,
	10, 13, 14,
	10, 14, 11,
	11, 14, 15,
};
#else
int table_inv[18 * 3] = 
{
	0, 1, 4,
	1, 5, 4,
	1, 2, 5,
	2, 6, 5,
	2, 3, 6,
	3, 7, 6,

	4, 5, 8,
	5, 9, 8,
	5, 6, 9,
	6, 10, 9,
	6, 7, 10,
	7, 11, 10,

	8, 9, 12,
	9, 13, 12,
	9, 10, 13,
	10, 14, 13,
	10, 11, 14,
	11, 15, 14,
};
#endif

ostream &operator<<(ostream &os, const Point3 &p)
{
	return os << p.x << " " << p.y << " " << p.z << endl;
}

int main()
{
	ifstream fin("/home/cym/program/OBJ/utah_teapot/t.bpt");
	//ifstream fin("/home/cym/program/OBJ/utah_teapot/teapotCGA.bpt");
	//ifstream fin("/home/cym/program/OBJ/utah_teapot/t2.bpt");
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
	
	//for (int p = 0; p < patch_number; ++p)
	//{
		//for (vector<Point3>::size_type i = 0; i < ctrl_points[i].size(); ++i)
		//{
			//cout << i << " " << ctrl_points[p][i] << endl;
		//}
		//cout << endl;
	//}

	ofstream fout("teapot.obj");

	for (int p = 0; p < patch_number; ++p)
	{
		for (int i = 0; i < 3 + 1; ++i)
		{
			for (int j = 0; j < 3 + 1; ++j)
			{
				fout << "v " << ctrl_points[p][i * 4 + j];
			}
		}
		int base = p * 16 + 1;
		for (int i = 0; i < 18; ++i)
		{
#ifdef COUNTER_CLOCK
			fout << "f " << base + table[i * 3] << " " <<
							base + table[i * 3 + 1] << " " <<
							base + table[i * 3 + 2] << endl;
#else
			fout << "f " << base + table_inv[i * 3] << " " <<
							base + table_inv[i * 3 + 1] << " " <<
							base + table_inv[i * 3 + 2] << endl;
#endif
		}
		//for (int i = 0; i < 4; ++i)
		//{
			//for (int j = 0; j < 4; ++j)
			//{
				
			//}
		//}
	}

	fout.close();
}
