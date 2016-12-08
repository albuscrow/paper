#include <iostream>
#include <fstream>

using namespace std;

int main()
{
	ifstream fin("/home/cym/program/OBJ/utah_teapot/teacup.bpt");
	ofstream fout("teapot.obj");

	int patch_number;
	fin >> patch_number;
	string line;
	while(getline(fin, line))
	{
		
	}
}
