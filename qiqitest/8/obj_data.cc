#include "obj_data.h"
#include <cmath>
#include <fstream>
#include <sstream>
#include <iostream>

namespace objdata
{
	using std::string;
	using std::vector;
	using std::ifstream;
	using std::istringstream;

	const Vector3 cross(const Vector3 &v0, const Vector3 &v1)
	{
		double x = v0.component1_ * v1.component2_ - v0.component2_ * v1.component1_;
		double y = v0.component2_ * v1.component0_ - v0.component0_ * v1.component2_;
		double z = v0.component0_ * v1.component1_ - v0.component1_ * v1.component0_;

		return Vector3(x, y, z);
	}

	Vector3 &Vector3::operator+=(const Vector3 &v)
	{
		component0_ += v.component0_;
		component1_ += v.component1_;
		component2_ += v.component2_;
		return *this;
	}

	void Vector3::normalize()
	{
		double length = norm();
		component0_ /= length;
		component1_ /= length;
		component2_ /= length;
	}

	/**************** Face类的成员函数 **************/
	Face::Face()
	{
		m_eFaceCase = V_0_0;
		m_nMtlIdx = -1;
	}

	Face::Face(const Face& face)
	{
		m_eFaceCase = face.m_eFaceCase;
		m_nMtlIdx = face.m_nMtlIdx;

		vertexCoordIndex = face.vertexCoordIndex;
		textureCoordIndex = face.textureCoordIndex;
		normalCoordIndex = face.normalCoordIndex;
	}

	/**************** MtlTex类的成员函数 **************/
	MtlRange::MtlRange(const string &mtlName, int begin)
	{
		m_sMtlName = mtlName;
		m_nBegin = begin;
	}

	MtlRange::MtlRange(const MtlRange &mtlRange)
	{
		m_sMtlName = mtlRange.m_sMtlName;
		m_nBegin = mtlRange.m_nBegin;
		m_nEnd = mtlRange.m_nEnd;
	}

	/**************** MtlTex类的成员函数 **************/
	MtlTex::MtlTex(const string mtlTexName)
	{
		m_sMtlTexName = mtlTexName;
		m_sTexFileName = string("");
		m_nTexBindingIdx = 2048;
		for (int i = 0; i < 3; ++i)
		{
			//m_fKa[i] = 0.5;
			//m_fKd[i] = 0.5;
			//m_fKs[i] = 0.9;
			m_fKa[i] = 0.2;
			m_fKd[i] = 1.0;
			m_fKs[i] = 0.8;
		}
	}

	MtlTex::MtlTex(const MtlTex &mtlTex)
	{
		m_sMtlTexName = mtlTex.m_sMtlTexName;
		for (int i = 0; i < 3; ++i)
		{
			m_fKa[i] = mtlTex.m_fKa[i];
			m_fKd[i] = mtlTex.m_fKd[i];
			m_fKs[i] = mtlTex.m_fKs[i];
		}
		m_sTexFileName = mtlTex.m_sTexFileName;
		m_nTexBindingIdx = mtlTex.m_nTexBindingIdx;
	}

	/**************** ObjData类的成员函数 **************/
	void ObjData::parseVertex(istringstream &iss)			// 解析顶点坐标
	{
		double x, y, z;

		iss >> x >> y >> z;
		vertexCoordList.push_back(VertexCoord(x, y, z));
	}

	void ObjData::parseTexture(istringstream &iss)			// 解析纹理坐标
	{
		double u, v, w;

		iss >> u >> v >> w;
		textureCoordList.push_back(TextureCoord(u, v, w));
	}

	void ObjData::parseNormal(istringstream &iss)			// 解析法向坐标
	{
		double i, j, k;

		iss >> i >> j >> k;
		//normalCoordList.push_back(NormalCoord(i, j, k));
		double length = sqrt(i * i + j * j + k * k);
		normalCoordList.push_back(NormalCoord(i / length, j / length, k / length));
	}

	void ObjData::parseFace(istringstream &iss, const string &str)	// 解析面信息
	{
		Face face;
		string useless;

		string::size_type firstSlashPos = str.find('/');
		if (firstSlashPos == string::npos)
			face.m_eFaceCase = V_0_0;		// 找不到'/'，说明仅有顶点编号
		else if (str[firstSlashPos + 1] == '/')
			face.m_eFaceCase = V_0_N;		// 两个'/'连续，说明是形如"3//4"的顶点
		else
			face.m_eFaceCase = V_T_N;		// 两个'/'不连续，说明是形如"1/2/3"的顶点

		int vIndex, vtIndex, vnIndex;
		vIndex = 0;
		vtIndex = 0;
		vnIndex = 0;
		char slash0, slash1;

		switch (face.m_eFaceCase)
		{
			case V_0_0:							// 顶点信息中不包含纹理坐标和法向坐标
				//while (iss >> vIndex >> slash0 >> slash1)
				while (iss >> vIndex)
				{
					if (vIndex > 0)				// 索引为正
						vIndex--;				// 索引从1开始，数组从零开始，需要调整一下
					else
						vIndex += vertexCoordList.size();	// 负数索引需要往上数相应的行
					face.vertexCoordIndex.push_back(vIndex);
				}
				break;
			case V_0_N:							// 顶点信息中不包含纹理坐标，包含法向坐标
				while(iss >> vIndex >> slash0 >> slash1 >> vnIndex)
				{
					if (vIndex > 0)
						vIndex--;
					else
						vIndex += vertexCoordList.size();
					if (vnIndex > 0)
						vnIndex--;
					else
						vnIndex += normalCoordList.size();
					face.vertexCoordIndex.push_back(vIndex);
					face.normalCoordIndex.push_back(vnIndex);
				}
				break;		
			case V_T_N:							// 顶点信息中包含纹理坐标和法向坐标
				while(iss >> vIndex >> slash0 >> vtIndex >> slash1 >> vnIndex)
				{
					if (vIndex > 0)
						vIndex--;
					else
						vIndex += vertexCoordList.size();
					if (vtIndex > 0)
						vtIndex--;
					else
						vtIndex += textureCoordList.size();
					if (vnIndex > 0)
						vnIndex--;
					else
						vnIndex += normalCoordList.size();
					face.vertexCoordIndex.push_back(vIndex);
					face.textureCoordIndex.push_back(vtIndex);
					face.normalCoordIndex.push_back(vnIndex);
				}
				break;
			default:
				break;
		}
		faceList.push_back(face);
	}

	string ObjData::parseMtlLib(istringstream &iss) const	// 解析mtl文件信息
	{
		string mtlName;

		iss >> mtlName;
		return string(mtlName);
	}

	void ObjData::parseUseMtl(istringstream &iss, vector<MtlRange> &mtlRangeList)	// 解析使用mtl的信息
	{
		string result;

		iss >> result;

		int size = mtlRangeList.size();
		if (size != 0)
			mtlRangeList[size - 1].m_nEnd = faceList.size() - 1;
		mtlRangeList.push_back(MtlRange(result, faceList.size()));
	}

	void ObjData::parseNewMtl(istringstream &iss)			// 解析新的mtl信息
	{
		string name;

		iss >> name;
		mtlTexList.push_back(MtlTex(name));
	}

	void ObjData::parseTexFile(istringstream &iss, const string &filePath)	// 解析纹理文件的信息
	{
		string texName;
		static unsigned int texBindingIdx = 0;

		iss >> texName;
		if (texName.size())
		{
			mtlTexList[mtlTexList.size() - 1].m_sTexFileName = filePath + texName;
			mtlTexList[mtlTexList.size() - 1].m_nTexBindingIdx = texBindingIdx++;
		}
	}

	void ObjData::parseKa(istringstream &iss)				// 解析ka信息
	{
		float ka[3];
		iss >> ka[0] >> ka[1] >> ka[2];
		for (int i = 0; i < 3; ++i)
			mtlTexList[mtlTexList.size() - 1].m_fKa[i] = ka[i];
	}

	void ObjData::parseKd(istringstream &iss)				// 解析kd信息
	{
		float kd[3];
		iss >> kd[0] >> kd[1] >> kd[2];
		for (int i = 0; i < 3; ++i)
			mtlTexList[mtlTexList.size() - 1].m_fKd[i] = kd[i];
	}

	void ObjData::parseKs(istringstream &iss)				// 解析ks信息
	{
		float ks[3];
		iss >> ks[0] >> ks[1] >> ks[2];
		for (int i = 0; i < 3; ++i)
			mtlTexList[mtlTexList.size() - 1].m_fKs[i] = ks[i];
	}

	void ObjData::readObj(const char *filePathName)			// 读取obj文件中的数据
	{
		/**************************** 初始化相关参数 ********************************/
		vertexCoordList.clear();
		textureCoordList.clear();
		normalCoordList.clear();
		faceList.clear();
		mtlTexList.clear();

		/***************************** 解析 obj 文件 ********************************/
		string mtlName;
		vector<MtlRange> mtlRangeList;
		ifstream inFile(filePathName);
		if (inFile != 0)
		{
			string line;
			while(getline(inFile, line))
			{
				istringstream iss(line);
				string primitive;
				iss >> primitive;
				if (primitive == "v")				// 这一行是顶点坐标
					parseVertex(iss);
				else if (primitive == "vt")			// 这一行是纹理坐标
					parseTexture(iss);
				else if (primitive == "vn")			// 这一行是法向坐标
					parseNormal(iss);
				else if (primitive == "f")			// 这一行是面信息
					parseFace(iss, line);
				else if (primitive == "mtllib")		// 这一行是mtl文件信息
					mtlName = parseMtlLib(iss);
				else if (primitive == "usemtl")		// 这一行是使用材质信息
					parseUseMtl(iss, mtlRangeList);
			}
			int size = mtlRangeList.size();
			if (size != 0)
				mtlRangeList[size - 1].m_nEnd = faceList.size() - 1;
		}
		std::cout << "mtlRangeListSize = " << mtlRangeList.size() << std::endl;
		inFile.close();
		inFile.clear();

		/***************************** 解析 mtl 文件 ********************************/
		string filePath(filePathName);
		int lastSlashPos = filePath.rfind('/');
		filePath = filePath.substr(0, lastSlashPos + 1);

		inFile.open((filePath + mtlName).c_str());
		if (inFile != 0)
		{
			string line;
			while(getline(inFile, line))
			{
				istringstream iss(line);
				string primitive;
				iss >> primitive;
				if (primitive == "newmtl")
					parseNewMtl(iss);
				else if (primitive == "map_Kd")
					parseTexFile(iss, filePath);
				else if (primitive == "Ka")
					parseKa(iss);
				else if (primitive == "Kd")
					parseKd(iss);
				else if (primitive == "Ks")
					parseKs(iss);
			}
			inFile.close();
			// 设置所有面片的材质编号
			for (std::vector<MtlRange>::size_type i = 0; i < mtlRangeList.size(); ++i)
			{
				string name = mtlRangeList[i].m_sMtlName;
				int begin = mtlRangeList[i].m_nBegin;
				int end = mtlRangeList[i].m_nEnd;
				int mtlIdx = -1;
				for (std::vector<MtlTex>::size_type j = 0; j < mtlTexList.size(); ++j)
				{
					if (mtlTexList[j].m_sMtlTexName == name)
					{
						mtlIdx = j;
						break;
					}
				}
				for (int j = begin; j <= end; ++j)
				{
					faceList[j].m_nMtlIdx = mtlIdx;
				}
			}
		}

		/************************************* 打印 ************************************/
		/*for (std::vector<MtlTex>::size_type i = 0; i < mtlTexList.size(); ++i)
		{
			cout << mtlTexList[i].m_sMtlTexName << endl;
			cout << "ka = " << mtlTexList[i].m_fKa[0] << " "
							<< mtlTexList[i].m_fKa[1] << " "
							<< mtlTexList[i].m_fKa[2] << " " << endl;
			cout << "kd = " << mtlTexList[i].m_fKd[0] << " "
							<< mtlTexList[i].m_fKd[1] << " "
							<< mtlTexList[i].m_fKd[2] << " " << endl;
			cout << "ks = " << mtlTexList[i].m_fKs[0] << " "
							<< mtlTexList[i].m_fKs[1] << " "
							<< mtlTexList[i].m_fKs[2] << " " << endl;
			cout << "texFileName = " << mtlTexList[i].m_sTexFileName << endl << endl;
		}
		cout << "*****************************************" << endl;

		cout << "mtlRangeList.size() = " << mtlRangeList.size() << endl;
		for (std::vector<MtlRange>::size_type i = 0; i < mtlRangeList.size(); ++i)
		{
			cout << mtlRangeList[i].m_sMtlName << ": ("
				 << mtlRangeList[i].m_nBegin << ", " << mtlRangeList[i].m_nEnd << ")" << endl;
		}*/

		/*FILE *fptr = fopen("ObjDataOutput.txt", "w");

		for (vector<Face>::size_type i = 0; i < faceList.size(); ++i)
		{
			fprintf(fptr, "%d\n", faceList[i].m_nMtlIdx);
		}

		fclose(fptr);*/
	}
}
