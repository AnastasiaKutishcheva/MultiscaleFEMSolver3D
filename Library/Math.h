#pragma once
#include <vector>
#include <functional>
//#include "DenseMatrix.h"
#include "Tensor.h"

#define M_PI 3.14159265358979323846
#define mu0 4*M_PI*1E-7
#define eps0 8.85*1E-12

namespace math
{
	double SolveLengthVector(Point<double>  A, Point<double> B);
	void ParserStringToVectorInt(char* A, std::vector<int>& B, const char* separator);
	template <typename T> std::vector<T> GetConfluence(std::vector<T>& A, std::vector<T>& B);

	struct SimpleGrid
	{
		std::vector<std::vector<int>> nvtr;
		std::vector<int> nvkat;
		std::vector<Point<double>> xyz;
		std::vector<int> vertex_map;

		bool read_from_zero;

		bool ReadFromSalomeDat(char *in_file, int dimension)
		{
			FILE *fin;
			fopen_s(&fin, in_file, "r");
			if (fin == NULL)
			{
				return false;
			}

			int elements_count, vertexes_count;
			fscanf_s(fin, "%d %d", &vertexes_count, &elements_count);

			this->xyz.resize(vertexes_count);
			for (int i = 0; i < vertexes_count; i++)
			{
				int _tmp;
				fscanf_s(fin, "%d %lf %lf %lf", &_tmp, &this->xyz[i].x, &this->xyz[i].y, &this->xyz[i].z);
			}

			for (int i = 0; i < elements_count; i++)
			{
				int type;
				std::vector<int> _tmp(dimension+1);
				fscanf_s(fin, "%d %d", &_tmp[0], &type);
				switch (type)
				{
				case 102: //1D elements
					fscanf_s(fin, "%d %d", &_tmp[0], &_tmp[1]);
					_tmp[0]--;
					_tmp[1]--;
					this->nvtr.push_back(_tmp);
					break;
				case 203: //2D elements
					fscanf_s(fin, "%d %d %d", &_tmp[0], &_tmp[1], &_tmp[2]);
					_tmp[0]--;
					_tmp[1]--;
					_tmp[2]--;
					this->nvtr.push_back(_tmp);
					break;
				case 304: //3D elements
					fscanf_s(fin, "%d %d %d %d", &_tmp[0], &_tmp[1], &_tmp[2], &_tmp[3]);
					_tmp[0]--;
					_tmp[1]--;
					_tmp[2]--;
					_tmp[3]--;
					this->nvtr.push_back(_tmp);
					break;
				default:
					break;
				}
			}
			this->nvkat.resize(this->nvtr.size());

			fclose(fin);
			return true;
		}
		bool ReadFromSalomeDat(char* in_file, int dimension, int map_size)
		{
			FILE* fin;
			fopen_s(&fin, in_file, "r");
			if (fin == NULL)
			{
				return false;
			}

			int elements_count, vertexes_count;
			fscanf_s(fin, "%d %d", &vertexes_count, &elements_count);
			
			this->vertex_map.resize(vertexes_count);

			this->xyz.resize(vertexes_count);
			for (int i = 0; i < vertexes_count; i++)
			{
				int _tmp;
				fscanf_s(fin, "%d %lf %lf %lf", &_tmp, &this->xyz[i].x, &this->xyz[i].y, &this->xyz[i].z);
				vertex_map[i] = _tmp - 1;
			}

			for (int i = 0; i < elements_count; i++)
			{
				int type;
				std::vector<int> _tmp(dimension + 1);
				fscanf_s(fin, "%d %d", &_tmp[0], &type);
				switch (type)
				{
				case 102: //1D elements
					fscanf_s(fin, "%d %d", &_tmp[0], &_tmp[1]);
					_tmp[0]--;
					_tmp[1]--;
					this->nvtr.push_back(_tmp);
					break;
				case 203: //2D elements
					fscanf_s(fin, "%d %d %d", &_tmp[0], &_tmp[1], &_tmp[2]);
					_tmp[0]--;
					_tmp[1]--;
					_tmp[2]--;
					this->nvtr.push_back(_tmp);
					break;
				case 304: //3D elements
					fscanf_s(fin, "%d %d %d %d", &_tmp[0], &_tmp[1], &_tmp[2], &_tmp[3]);
					_tmp[0]--;
					_tmp[1]--;
					_tmp[2]--;
					_tmp[3]--;
					this->nvtr.push_back(_tmp);
					break;
				default:
					break;
				}
			}
			this->nvkat.resize(this->nvtr.size());

			fclose(fin);
			return true;
		}
		bool ReadFromNVTR(char *directory, int vertex_count_in_element)
		{
			char name_nvtr[1000], name_nvkat[1000], name_xyz[1000], name_param[1000];
			sprintf_s(name_nvtr, sizeof(name_nvtr), "%s/nvtr.txt", directory);
			sprintf_s(name_nvkat, sizeof(name_nvkat), "%s/nvkat.txt", directory);
			sprintf_s(name_xyz, sizeof(name_xyz), "%s/xyz.txt", directory);
			sprintf_s(name_param, sizeof(name_param), "%s/param.txt", directory);
			
			FILE *f_nvtr, *f_xyz, *f_nvkat, *f_param;
			fopen_s(&f_nvtr, name_nvtr, "r");
			if (f_nvtr == NULL)
			{
				return false;
			}
			fopen_s(&f_param, name_param, "r");
			if (f_param == NULL)
			{
				return false;
			}
			fopen_s(&f_xyz, name_xyz, "r");
			if (f_xyz == NULL)
			{
				return false;
			}
			fopen_s(&f_nvkat, name_nvkat, "r");

			int elements_count, vertexes_count;
			fscanf_s(f_param, "%d %d", &vertexes_count, &elements_count);

			this->xyz.resize(vertexes_count);
			for (int i = 0; i < vertexes_count; i++)
			{
				Point<double> x;
				fscanf_s(f_xyz, "%lf %lf %lf", &x.x, &x.y, &x.z);
				this->xyz[i] =x;
			}
			vertexes_count = (int)this->xyz.size();

			bool find_0 = false;
			for (int i = 0; i < elements_count; i++)
			{
				std::vector<int> _tmp;

				if (vertex_count_in_element == -1) //polyhedron;
				{
					int n_nodes;
					fscanf_s(f_nvtr, "%d", &n_nodes);
					_tmp.push_back(n_nodes);
					for (int n = 0; n < n_nodes; n++)
					{
						int _n;
						fscanf_s(f_nvtr, "%d", &_n);
						_tmp.push_back(_n);
					}
					int num_faces;
					fscanf_s(f_nvtr, "%d", &num_faces);
					_tmp.push_back(num_faces);
					for (int f = 0; f < num_faces; f++)
					{
						int n_in_face;
						fscanf_s(f_nvtr, "%d", &n_in_face);
						_tmp.push_back(n_in_face);
						for (int nf = 0; nf < n_in_face; nf++)
						{
							int _n;
							fscanf_s(f_nvtr, "%d", &_n);
							_tmp.push_back(_n);
						}
					}
				}
				else {
					_tmp.resize(vertex_count_in_element);
					for (int ii = 0; ii < vertex_count_in_element; ii++)
					{
						fscanf_s(f_nvtr, "%d", &_tmp[ii]);
						//if (_tmp[ii] == 0) find_0 = true;
					}
				}

				if (feof(f_nvtr) == 0)
				{
					for (int ii = 0; ii < _tmp.size(); ii++)
					{
						if (_tmp[ii] == 0) find_0 = true;
					}
					this->nvtr.push_back(_tmp);
				}
				else break;
			}
			elements_count = (int)this->nvtr.size();
			this->read_from_zero = true;
			if (!find_0)
			{
				this->read_from_zero = false;
				for (int i = 0; i < elements_count; i++)
				{
					if(vertex_count_in_element == -1)
					{
						int ii = 0;
						int n_nodes = this->nvtr[i][ii];
						for (int n = 0, ii = 1; n < n_nodes; n++)
						{
							this->nvtr[i][ii]--;
							ii++;
						}
						int faces = this->nvtr[i][ii];
						ii++;
						for (int f = 0; f < faces; f++)
						{
							int n_faces = this->nvtr[i][ii];
							ii++;
							for (int nf = 0; nf < n_faces; nf++)
							{
								this->nvtr[i][ii]--;
								ii++;
							}
						}
					}
					else {
						for (int ii = 0; ii < vertex_count_in_element; ii++)
						{
							this->nvtr[i][ii] --;
						}
					}
				}
			}

			this->nvkat.resize(this->nvtr.size());
			
			if (f_nvkat != NULL)
			{
				for (int i = 0; i < elements_count; i++)
				{
					fscanf_s(f_nvkat, "%d", &this->nvkat[i]);
				}
				fclose(f_nvkat);
			}

			fclose(f_xyz);
			fclose(f_nvtr);
			return true;
		}
		void ReadFacesBoundaryNVTR(char* directory, std::vector<std::vector<std::vector<int>>>& boundary_faces)
		{
			char boundary_file[1000];
			sprintf_s(boundary_file, sizeof(boundary_file), "%s/kraev2_faces.txt", directory);
			FILE* fin;
			fopen_s(&fin, boundary_file, "r");
			int N;
			fscanf_s(fin, "%d", &N);
			std::vector<std::vector<int>> tmp_vect2;
			std::vector<int> tmp_vect1;
			for (int i = 0; i < N; i++)
			{
				int parent_elem0, parent_elem1, type;
				double value;
				std::vector<int> id_face(3);
				fscanf_s(fin, "%d %d %d %d %d %d", &parent_elem0, &parent_elem1, &id_face[0], &id_face[1], &id_face[2], &type);
				parent_elem0 -= 2;
				parent_elem1 -= 2;
				id_face[0]--;
				id_face[1]--;
				id_face[2]--;
				type--;

				while (type >= boundary_faces.size())
				{
					boundary_faces.push_back(tmp_vect2);
				}
				boundary_faces[type].push_back(tmp_vect1);
				boundary_faces[type][boundary_faces[type].size() - 1].push_back(parent_elem0 < 0 ? parent_elem1 : parent_elem0);
				boundary_faces[type][boundary_faces[type].size() - 1].push_back(id_face[0]);
				boundary_faces[type][boundary_faces[type].size() - 1].push_back(id_face[1]);
				boundary_faces[type][boundary_faces[type].size() - 1].push_back(id_face[2]);
			}
			fclose(fin);
		}
		void ReadFromMSH(char* fileMSH, double scale, int vertex_count_in_element, int position_of_material, std::vector<std::vector<std::vector<int>>>& boundary_faces)
		{
			char separator[3] = { ' ', '\t', '\0' };
			char Str[1000];
			std::vector<int> Str_int;
			bool modified_file = false;

			FILE* fin;

			fopen_s(&fin, fileMSH, "r");
			fgets(Str, sizeof(Str), fin);
			fgets(Str, sizeof(Str), fin);
			ParserStringToVectorInt(Str, Str_int, separator);
			if (Str_int.size() == 4)
			{
				printf_s("It is modified file\n");
				modified_file = true;
			}
			fgets(Str, sizeof(Str), fin);
			fgets(Str, sizeof(Str), fin);

			int num_node;
			fscanf_s(fin, "%d", &num_node);

			xyz.resize(num_node);
			for (int i = 0; i < num_node; i++)
			{
				int tmp;
				fscanf_s(fin, "%d %lf %lf %lf", &tmp, &xyz[i].x, &xyz[i].y, &xyz[i].z);
				xyz[i].x *= scale;
				xyz[i].y *= scale;
				xyz[i].z *= scale;
			}

			fgets(Str, sizeof(Str), fin);
			fgets(Str, sizeof(Str), fin);
			fgets(Str, sizeof(Str), fin);

			int FullNum;
			fgets(Str, sizeof(Str), fin);
			ParserStringToVectorInt(Str, Str_int, separator);
			//fscanf_s(fin, "%d", &FullNum);
			FullNum = Str_int[0];
			printf_s("FullNum = %d\n", FullNum);

			//читаем ребра
			fgets(Str, sizeof(Str), fin);
			ParserStringToVectorInt(Str, Str_int, separator);
			int edge_size = 0;
			while (Str_int.size() == 8 && Str_int[1] == 1)
			{
				edge_size++;
				fgets(Str, sizeof(Str), fin);
				ParserStringToVectorInt(Str, Str_int, separator);
			}
			//читаем границы в элементы-точки
			//fgets(Str, sizeof(Str), fin);
			//Parser_String_to_VectInt(Str, Str_int, separator);
			int face_size = 0;
			printf_s("Str = %s", Str);
			printf_s("Str_int = {");
			for (int i = 0; i < Str_int.size(); i++)
			{
				printf_s("%d; ", Str_int[i]);
			}
			printf_s("}\n");

			while (Str_int.size() == 9 && Str_int[1] == 2)
			{
				face_size++;
				int offset = 6;
				std::vector<std::vector<int>> tmp_vect2;
				std::vector<int> tmp_vect1;

				while (Str_int[3] >= boundary_faces.size())
				{
					boundary_faces.push_back(tmp_vect2);
				}
				boundary_faces[Str_int[3]].push_back(tmp_vect1);

				if (modified_file)
				{
					boundary_faces[Str_int[3]][boundary_faces[Str_int[3]].size() - 1].push_back(Str_int[4] - 1);
				}
				else {
					boundary_faces[Str_int[3]][boundary_faces[Str_int[3]].size() - 1].push_back(-1);
				}

				for (int ii = 0; ii < 3; ii++)
					boundary_faces[Str_int[3]][boundary_faces[Str_int[3]].size()-1].push_back(Str_int[offset + ii] - 1);
				
				fgets(Str, sizeof(Str), fin);
				math::ParserStringToVectorInt(Str, Str_int, separator);
			}
			for (int i = 0; i < boundary_faces.size(); i++)
			{
				printf_s("Boundary[%d] = %d; ", i, boundary_faces[i].size());
			}
			printf_s("\n");

			nvtr.reserve(FullNum - edge_size - face_size);
			nvkat.reserve(FullNum - edge_size - face_size);
			bool flag = false;
			{
				int i = 0;
				//fgets(Str, sizeof(Str), fin);
				while (Str_int.size() == 10 /*&& Str_int[0] != FullNum*/)
				{
					int offset = 6;
					//if (Str_int.size() == 8) offset = 5;

					std::vector<int> nums(vertex_count_in_element);
					for (int j = 0; j < vertex_count_in_element; j++)
					{
						nums[j] = Str_int[offset + j] - 1;
					}

					int id_dom = Str_int[position_of_material - 1];//Str_int[3]-1;
					nvtr.push_back(nums);
					nvkat.push_back(id_dom);

					fgets(Str, sizeof(Str), fin);
					math::ParserStringToVectorInt(Str, Str_int, separator);
					i++;
				}
			};
			fclose(fin);

			if (!modified_file)
			{
				printf_s("\n");

				for (int id_type = 0; id_type < boundary_faces.size(); id_type++)
				{
					for (int id_triang = 0; id_triang < boundary_faces[id_type].size(); id_triang++)
					{
						printf_s("Create boundary topology %d/%d (triangle %d/%d)\r", id_type, boundary_faces.size(), id_triang, boundary_faces[id_type].size());
						int test_vertex = -1;
						int base_elem = -1;
						for (int id_elem = 0; id_elem < nvtr.size(); id_elem++)
						{
							auto _tmp = math::GetConfluence(nvtr[id_elem], boundary_faces[id_type][id_triang]);
							if (_tmp.size() == boundary_faces[id_type][id_triang].size() - 1)
							{
								base_elem = id_elem;
								break;
							}
						}
						boundary_faces[id_type][id_triang][0] = base_elem;
					}
				}

				fopen_s(&fin, fileMSH, "w");
				fprintf_s(fin, "$MeshFormat\n");
				fprintf_s(fin, "2 0 8 1\n");
				fprintf_s(fin, "$EndMeshFormat\n");
				fprintf_s(fin, "$Nodes\n");
				fprintf_s(fin, "\t%d\n", xyz.size());
				for (int i = 0; i < xyz.size(); i++)
				{
					fprintf_s(fin, "\t\t%d\t%.16lf\t%.16lf\t%.16lf\n", i + 1, xyz[i].x, xyz[i].y, xyz[i].z);
				}
				fprintf_s(fin, "$EndNodes\n");
				fprintf_s(fin, "$Elements\n");
				fprintf_s(fin, "\t%d\n", FullNum);
				int i = 0;
				for (int b = 0; b < boundary_faces.size(); b++)
					for (int j = 0; j < boundary_faces[b].size(); j++)
					{
						fprintf_s(fin, "\t\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
							i + 1, 2, 3, b, boundary_faces[b][j][0]+1, 0,
							boundary_faces[b][j][1]+1, boundary_faces[b][j][2]+1, boundary_faces[b][j][3]+1);
						i++;
					}
				for (int j = 0; j < nvtr.size(); j++)
				{
					fprintf_s(fin, "\t\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
						i + 1, 4, 3, nvkat[j], nvkat[j], 0,
						nvtr[j][0]+1, nvtr[j][1]+1, nvtr[j][2]+1, nvtr[j][3]+1);
					i++;
				}
				fprintf_s(fin, "$EndElements\n");
				fclose(fin);
				
			}
		}
		void WriteMSH(char* fileMSH, int position_of_material, std::vector<std::vector<std::vector<int>>>& boundary_faces)
		{
			FILE* fin;

			fopen_s(&fin, fileMSH, "w");
			fprintf_s(fin, "$MeshFormat\n");
			fprintf_s(fin, "2 0 8 1\n");
			fprintf_s(fin, "$EndMeshFormat\n");
			fprintf_s(fin, "$Nodes\n");
			fprintf_s(fin, "\t%d\n", xyz.size());
			for (int i = 0; i < xyz.size(); i++)
			{
				fprintf_s(fin, "\t\t%d\t%.16lf\t%.16lf\t%.16lf\n", i + 1, xyz[i].x, xyz[i].y, xyz[i].z);
			}
			fprintf_s(fin, "$EndNodes\n");
			fprintf_s(fin, "$Elements\n");
			fprintf_s(fin, "\t%d\n", boundary_faces.size() + nvtr.size());
			int i = 0;
			for (int b = 0; b < boundary_faces.size(); b++)
				for (int j = 0; j < boundary_faces[b].size(); j++)
				{
					fprintf_s(fin, "\t\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
						i + 1, 2, 3, b, boundary_faces[b][j][0] + 1, 0,
						boundary_faces[b][j][1] + 1, boundary_faces[b][j][2] + 1, boundary_faces[b][j][3] + 1);
					i++;
				}
			for (int j = 0; j < nvtr.size(); j++)
			{
				fprintf_s(fin, "\t\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n",
					i + 1, 4, 3, nvkat[j]+1, nvkat[j]+1, 0,
					nvtr[j][0] + 1, nvtr[j][1] + 1, nvtr[j][2] + 1, nvtr[j][3] + 1);
				i++;
			}
			fprintf_s(fin, "$EndElements\n");
			fclose(fin);
		}

		int TransferIdSelfIntoGlobal(int local_id)
		{
			return this->vertex_map[local_id];
		}
		int TransferIdGlobalIntoSelf(int global_id)
		{
			for (int i = 0; i < vertex_map.size(); i++)
			{
				if (vertex_map[i] == global_id)
					return i;
			}
			return -1;
		}

		void printTecPlot3D(FILE* fdat, std::vector<std::vector<double>>& value, std::vector<std::vector<char>> name_value, char* name_zone)
		{
			int DEL_domain = 10;
			fprintf_s(fdat, "TITLE     = \"numerical\"\n");
			fprintf_s(fdat, "VARIABLES = \"x\"\n \"y\"\n \"z\"\n");
			for (int i = 0; i < value.size(); i++)
			{
				fprintf_s(fdat, " \"");
				for (int j = 0; j < value[i].size(); j++)
				{
					if (name_value[i][j] != '\0')
					{
						fprintf_s(fdat, "%c", name_value[i][j]);
					}
					else
					{
						break;
					}
				}
				fprintf_s(fdat, "\"\n");
			}

			int num_elem = this->nvtr.size();
			fprintf_s(fdat, "ZONE T=\"%s\"\n", name_zone);
			fprintf_s(fdat, " N=%d,  E=%d, F=FEBLOCK ET=", this->xyz.size(), num_elem);
			switch (this->nvtr[0].size())
			{
			case 4: fprintf_s(fdat, "Tetrahedron "); break;
			case 3: fprintf_s(fdat, "Triangle "); break;
			default:
				break;
			}
			fprintf_s(fdat, "\n");
			fprintf_s(fdat, " VARLOCATION=(NODAL NODAL NODAL");
			for (int i = 0; i < value.size(); i++)
			{
				fprintf_s(fdat, " CELLCENTERED");
			}
			fprintf_s(fdat, ")\n");

			for (int i = 0; i < this->xyz.size(); i++)
				fprintf_s(fdat, "%.10e\n", this->xyz[i].x);
			fprintf_s(fdat, "\n");
			for (int i = 0; i < this->xyz.size(); i++)
				fprintf_s(fdat, "%.10e\n", this->xyz[i].y);
			fprintf_s(fdat, "\n");
			for (int i = 0; i < this->xyz.size(); i++)
				fprintf_s(fdat, "%.10e\n", this->xyz[i].z);
			fprintf_s(fdat, "\n");

			for (int i = 0; i < value.size(); i++)
			{
				for (int j = 0; j < this->nvtr.size(); j++)
				{
					fprintf_s(fdat, "%.10lf\n", value[i][j]);
				}
				fprintf_s(fdat, "\n");
			}

			for (int i = 0; i < this->nvtr.size(); i++)
			{
				for (int j = 0; j < this->nvtr[i].size(); j++)
					fprintf_s(fdat, "%d ", this->nvtr[i][j] + 1);

				fprintf_s(fdat, "\n");
			}

			fclose(fdat);
		}
	};

	//------------->algebra
	bool IsEqual(double A, double B)
	{
		if (abs(A - B) < 1e-10) return true;
		return false;
	}
	bool IsEqual(Point<double> A, Point<double> B)
	{
		if (SolveLengthVector(A, B) < 1e-10) return true;
		return false;
	}
	template <typename T>
	bool IsEqual(T A, T B)
	{
		if (A != B) return false;
		return true;
	}
	bool IsEqual(int A, int B)
	{
		if (A == B) return true;
		return false;
	}

	double GetRound(double X, int num)
	{
		//num = 8;
		double result = X;

		/*double p_d = pow(10, num);
		long_int p_i = pow_long(10, num);
		long_int p_inv = pow_long(10, num);
		p_inv.inv = true;
		result = p_inv * round(p_i*X);*/

		char text[1000];
		switch (num)
		{
		case 0: sprintf_s(text, "%.0lf", X); break;
		case 1: sprintf_s(text, "%.1lf", X); break;
		case 2: sprintf_s(text, "%.2lf", X); break;
		case 3: sprintf_s(text, "%.3lf", X); break;
		case 4: sprintf_s(text, "%.4lf", X); break;
		case 5: sprintf_s(text, "%.5lf", X); break;
		case 6: sprintf_s(text, "%.6lf", X); break;
		case 7: sprintf_s(text, "%.7lf", X); break;
		case 8: sprintf_s(text, "%.8lf", X); break;
		case 9: sprintf_s(text, "%.9lf", X); break;
		case 10: sprintf_s(text, "%.10lf", X); break;
		case 11: sprintf_s(text, "%.11lf", X); break;
		case 12: sprintf_s(text, "%.12lf", X); break;
		case 13: sprintf_s(text, "%.13lf", X); break;
		case 14: sprintf_s(text, "%.14lf", X); break;
		default: sprintf_s(text, "%.0lf", X); break;
		}
		sscanf_s(text, "%lf", &result);

		return result;
	}

	int GetSignum(double X)
	{
		/*X = round(X, 5);*/
		if (IsEqual(X, 0.0)) return 0;
		if (X < 0.0) return -1;
		if (X > 0.0) return 1;
		return 0;
	}
	double SolvePow(double X, int P)
	{
		if (P > 0)
			return pow(X, P);
		if (P == 0)
			return 0;
		if (P < 0)
			return 1. / pow(X, P);
	}

	void BisectionMethod_iteration(double A, double B, std::function<double(double)>& F, std::vector<double>& S, double _eps)
	{
		if (GetSignum(F(A)) != GetSignum(F(B)))
		{
			if (abs(B - A) > _eps)
			{
				BisectionMethod_iteration(A, (A + B) / 2., F, S, _eps);
				BisectionMethod_iteration((A + B) / 2., B, F, S, _eps);
			}
			else {
				S.push_back((A + B) / 2.);
			}
		}
	}
	void SolveEquationViaBisectionMethod(double A, double B, std::function<double(double)> F, std::vector<double>& S)
	{
		if (A < B)
		{
			BisectionMethod_iteration(A, B, F, S, 1E-6);
		}
		else {
			BisectionMethod_iteration(B, A, F, S, 1E-6);
		}
	}

	//-------Matrix_Vector operations
	//result[][]=A[][]*B[][]
	template<typename T>
	void MultMatrixTMatrix(std::vector<std::vector<T>>& A_T, std::vector<std::vector<T>>& B, std::vector<std::vector<T>>& result)
	{
		try
		{
			result.resize(A_T[0].size());
			for (int i = 0; i < result.size(); i++)
			{
				result[i].resize(B[0].size());
				for (int j = 0; j < result[i].size(); j++)
				{
					T tmp = 0;
					for (int k = 0; k < A_T.size(); k++)
					{
						tmp += A_T[k][i] * B[k][j];
					}
					result[i][j] = tmp;
				}
			}
		}
		catch (const std::exception&)
		{
			printf_s("Error: Library/Math.h/void MultMatrixMatrix(std::vector<std::vector<T>> &A, std::vector<std::vector<T>> &B, std::vector<std::vector<T>> &result)\n");
			if (A_T[0].size() != B.size())
			{
				printf("\t\tError in mult matrix (%d,%d)x(%d,%d)\n", (int)A_T.size(), (int)A_T[0].size(), (int)B.size(), (int)B[0].size());
			}
		}
	}
	template<typename T>
	void MultMatrixMatrix(std::vector<std::vector<T>>& A, std::vector<std::vector<T>>& B, std::vector<std::vector<T>>& result)
	{
		try
		{
			result.resize(A.size());
			for (int i = 0; i < result.size(); i++)
			{
				result[i].resize(B[0].size());
				for (int j = 0; j < result[i].size(); j++)
				{
					T tmp = 0;
					for (int k = 0; k < A[0].size(); k++)
					{
						tmp += A[i][k] * B[k][j];
					}
					result[i][j] = tmp;
				}
			}
		}
		catch (const std::exception&)
		{
			printf_s("Error: Library/Math.h/void MultMatrixMatrix(std::vector<std::vector<T>> &A, std::vector<std::vector<T>> &B, std::vector<std::vector<T>> &result)\n");
			if (A[0].size() != B.size())
			{
				printf("\t\tError in mult matrix (%d,%d)x(%d,%d)\n", (int)A.size(), (int)A[0].size(), (int)B.size(), (int)B[0].size());
			}
		}
	}
	template<typename T>
	void MakeSummMatrixMatrix(std::vector<std::vector<T>>& A, std::vector<std::vector<T>>& B, std::vector<std::vector<T>>& result)
	{
		try
		{
			result.resize(A.size());
			for (int i = 0; i < result.size(); i++)
			{
				result[i].resize(B[0].size());
				for (int j = 0; j < result[i].size(); j++)
				{
					result[i][j] = A[i][j] + B[i][j];
				}
			}
		}
		catch (const std::exception&)
		{
			printf_s("Error: Library/Math.h/void MultMatrixMatrix(std::vector<std::vector<T>> &A, std::vector<std::vector<T>> &B, std::vector<std::vector<T>> &result)\n");
			if (A[0].size() != B.size())
			{
				printf("\t\tError in mult matrix (%d,%d)x(%d,%d)\n", (int)A.size(), (int)A[0].size(), (int)B.size(), (int)B[0].size());
			}
		}
	}
	template<typename T>
	void MultMatrixScalar(std::vector<std::vector<T>>& A, T B, std::vector<std::vector<T>>& result)
	{
		try
		{
			result.resize(A.size());
			for (int i = 0; i < result.size(); i++)
			{
				result[i].resize(A.size());
				for (int j = 0; j < result[i].size(); j++)
				{
					result[i][j] = A[i][j] * B;
				}
			}
		}
		catch (const std::exception&)
		{
			printf_s("Error: Library/Math.h/void MultMatrixMatrix(std::vector<std::vector<T>> &A, std::vector<std::vector<T>> &B, std::vector<std::vector<T>> &result)\n");
			if (A[0].size() != B.size())
			{
				printf("\t\tError in mult matrix (%d,%d)x(%d,%d)\n", (int)A.size(), (int)A[0].size(), (int)B.size(), (int)B[0].size());
			}
		}
	}
	template<typename T1, typename T2>
	void MultiplicationMatrixVector(std::vector<std::vector<T1>>& A, std::vector<T2>& B, std::vector<T2>& result)
	{
		try
		{
			result.resize(A.size());
			for (int i = 0; i < result.size(); i++)
			{
				T2 tmp = 0;
				for (int j = 0; j < A[i].size(); j++)
				{
					tmp += A[i][j] * B[j];
				}
				result[i] = tmp;
			}
		}
		catch (const std::exception&)
		{
			printf_s("Error: Library/Math.h/void MultiplicationMatrixVector(std::vector<std::vector<T1>> &A, std::vector<T2> &B, std::vector<T2> &result)\n");

		}
	}
	template<typename T1, typename T2>
	void MultiplicationTMatrixVector(std::vector<std::vector<T1>>& A_T, std::vector<T2>& B, std::vector<T2>& result)
	{
		try
		{
			result.resize(A_T[0].size());
			for (int i = 0; i < B.size(); i++)
			{
				for (int j = 0; j < A_T[i].size(); j++)
				{
					result[j] += A_T[i][j] * B[i];
				}
			}
		}
		catch (const std::exception&)
		{
			printf_s("Error: Library/Math.h/void MultiplicationMatrixVector(std::vector<std::vector<T1>> &A, std::vector<T2> &B, std::vector<T2> &result)\n");

		}
	}
	template<typename T>
	void ResizeVector(std::vector<std::vector<T>>& A, int size_I, int size_J)
	{
		A.resize(size_I);
		for (int i = 0; i < size_I; i++)
		{
			A[i].resize(size_J);
		}
	}
	template<typename T1, typename T1k, typename T2, typename T2k, typename TResult>
	void MakeSummVectors(std::vector<T1>& A, T1k koefA, std::vector<T2>& B, T2k koefB, std::vector<TResult>& result)
	{
		try
		{
			result.resize(A.size());
			for (int i = 0; i < result.size(); i++)
			{
				result[i] = A[i] * koefA + B[i] * koefB;
			}
		}
		catch (const std::exception&)
		{
			printf_s("Error: Library/Math.h/void MakeSummVectors(std::vector<T1> &A, T1k koefA, std::vector<T2> &B, T2k koefB, std::vector<TResult> &Result)\n");
		}
	}
	template<typename T, typename T2>
	void MultiplicationVector(std::vector<T>& A, T2 koef, std::vector<T>& result)
	{
		try
		{
			if (result.size() != A.size()) result.resize(A.size());
			for (int i = 0; i < A.size(); i++)
			{
				result[i] = A[i] * koef;
			}
		}
		catch (const std::exception&)
		{
			printf_s("Error: Library/Math.h/void MultiplicationVector(std::vector<T> &A, T koef, std::vector<T> &result)\n");
		}
	}
	template<typename T1, typename T2>
	double MakeInnerProduct(std::vector<T1>& A, std::vector<T2>& B)
	{
		try
		{
			double result = 0.0;
			double tmp = 0;
//#pragma omp parallel shared(result) private(tmp)
//#pragma omp for reduction(+:result) nowait
			for (int i = 0; i < A.size(); i++)
			{
				tmp = A[i] * B[i];
				result += tmp;
			}

			return result;
		}
		catch (const std::exception&)
		{
			printf_s("Error: Library/Math.h/double MakeInnerProduct(std::vector<T1> &A, std::vector<T2> &B)\n");
		}
	}
	template<typename T1>
	double MakeNorma_inH1(std::vector<T1>& A)
	{
		try
		{
			//return sqrt(MakeInnerProduct(A, A));
			return sqrt(MakeInnerProduct(A, A));
		}
		catch (const std::exception&)
		{
			printf_s("Error: Library/Math.h/double MakeInnerProduct(std::vector<T1> &A, std::vector<T2> &B)\n");
		}
	}

	//---------->vectors|matrix
	
	//быстрая сортировка вектора А
	//low - начало сортировки
	//high - конец сортировки
	template <typename T>
	void MakeQuickSort(std::vector <T>& A, int low, int high)
	{
		if (high >= low)
		{
			int i = low;
			int j = high;
			T x = A[(low + high) / 2];
			do {
				while (A[i] < x) ++i;
				while (A[j] > x) --j;
				if (i <= j) {
					T temp = A[i];
					A[i] = A[j];
					A[j] = temp;
					i++; j--;
				}
			} while (i <= j);

			if (low < j) MakeQuickSort(A, low, j);
			if (i < high) MakeQuickSort(A, i, high);
		}
	};
	template <typename T>
	void MakeQuickSort(std::vector <T>& A)
	{
		MakeQuickSort(A, 0, (int)A.size() - 1);
	};
	//удаление повтоений в отсортированном векторе А
	//результат в В
	//eps - с какой точностью считать числа равными
	template <typename T_old, typename T_new>
	void MakeRemovalOfDuplication(std::vector <T_old>& A, std::vector <T_new>& B)
	{
		if (A.size() != 0)
		{
			int k = 0, j = 0;

			B.push_back(T_new());
			B[k] = A[0];
			k++;

			for (int i = 1; i < A.size(); i++)
			{
				if (!IsEqual(A[j], A[i])) //не повтор
				{
					j = i;
					B.push_back(T_new());
					B[k] = A[i];
					k++;
				}
			}
		}
	}

	template <class T>
	void MakeCopyVector_A_into_B(std::vector<T>& A, std::vector<T>& B)
	{
		if (A.size() != 0)
		{
			B.resize(A.size());
			for (int i = 0; i < A.size(); i++)
			{
				B[i] = A[i];
			}
		}
	}
	void MakeCopyVector_A_into_B(std::vector<double>& A, std::vector<Point<double>>& B)
	{
		if (A.size() != 0)
		{
			B.resize(A.size() / 3);
			for (int i = 0; i < B.size(); i++)
			{
				B[i].x = A[i * 3 + 0];
				B[i].y = A[i * 3 + 1];
				B[i].z = A[i * 3 + 2];
			}
		}
	}

	template <typename T1, typename T2>
	void InitializationVector(std::vector<T1>& A, T2 b)
	{
		for (int i = 0; i < A.size(); i++)
		{
			A[i] = b;
		}
	}
	template <typename T1, typename T2>
	void InitializationVector(std::vector<std::vector<T1>>& A, T2 b)
	{
		for (int i = 0; i < A.size(); i++)
		{
			for (int j = 0; j < A[i].size(); j++)
			{
				A[i][j] = b;
			}
		}
	}

	template <typename T>
	std::vector<T> GetConfluence(std::vector<T>& A, std::vector<T>& B)
	{
		std::vector<T> tmp;
		for (int i = 0; i < A.size(); i++)
			for (int j = 0; j < B.size(); j++)
				if (A[i] == B[j])
					tmp.push_back(A[i]);
		return tmp;
	}

	template <typename T>
	T GetPositionInSortVector(std::vector<T>& A, T& a)
	{
		
		int position = -1;
		for (int i = 0; i < A.size() && a >= A[i]; i++)
		{
			if (a == A[i])
			{
				position = i;
			}
		}
		return position;
	}



	//тензорное произведение точек
	template<typename T>
	std::vector<std::vector<T>> GetTensorMult(Point<T>& A, Point<T>& B)
	{
		std::vector<std::vector<T>> res(3);
		res[0].resize(3);
		res[1].resize(3);
		res[2].resize(3);
		res[0][0] = A.x*B.x;
		res[0][1] = A.x*B.y;
		res[0][2] = A.x*B.z;
		res[1][0] = A.y*B.x;
		res[1][1] = A.y*B.y;
		res[1][2] = A.y*B.z;
		res[2][0] = A.z*B.x;
		res[2][1] = A.z*B.y;
		res[2][2] = A.z*B.z;

		return res;
	}
	template<typename T>
	void GetTensorMult(Point<T>& A, Point<T>& B, Tensor2Rank3D &res)
	{
		res.val[0][0] = A.x * B.x;
		res.val[0][1] = A.x * B.y;
		res.val[0][2] = A.x * B.z;
		res.val[1][0] = A.y * B.x;
		res.val[1][1] = A.y * B.y;
		res.val[1][2] = A.y * B.z;
		res.val[2][0] = A.z * B.x;
		res.val[2][1] = A.z * B.y;
		res.val[2][2] = A.z * B.z;
	}
	template<typename T>
	void GetTensorMult(Point<T>& A, Point<T>& B, std::vector<std::vector<T>> &res)
	{
		res.resize(3);
		res[0].resize(3);
		res[1].resize(3);
		res[2].resize(3);
		res[0][0] = A.x*B.x;
		res[0][1] = A.x*B.y;
		res[0][2] = A.x*B.z;
		res[1][0] = A.y*B.x;
		res[1][1] = A.y*B.y;
		res[1][2] = A.y*B.z;
		res[2][0] = A.z*B.x;
		res[2][1] = A.z*B.y;
		res[2][2] = A.z*B.z;

		return;
	}

	double GetDeterminantForMatrix3x3(double A[3][3])
	{
		return A[0][0] * A[1][1] * A[2][2] + A[0][1] * A[1][2] * A[2][0] + A[1][0] * A[0][2] * A[2][1]
			- A[0][2] * A[1][1] * A[2][0] - A[0][0] * A[2][1] * A[1][2] - A[0][1] * A[1][0] * A[2][2];
	}
	double GetDeterminantForMatrix3x3(std::vector<std::vector<double>> &A)
	{
		return A[0][0] * A[1][1] * A[2][2] + A[0][1] * A[1][2] * A[2][0] + A[1][0] * A[0][2] * A[2][1]
			- A[0][2] * A[1][1] * A[2][0] - A[0][0] * A[2][1] * A[1][2] - A[0][1] * A[1][0] * A[2][2];
	}
	double GetDeterminantForMatrix3x3(Tensor2Rank3D &A)
	{
		return A.val[0][0] * A.val[1][1] * A.val[2][2] + A.val[0][1] * A.val[1][2] * A.val[2][0] + A.val[1][0] * A.val[0][2] * A.val[2][1]
			- A.val[0][2] * A.val[1][1] * A.val[2][0] - A.val[0][0] * A.val[2][1] * A.val[1][2] - A.val[0][1] * A.val[1][0] * A.val[2][2];
	}
	

	std::vector< std::vector <double> > SolveInverseMatrix3x3(std::vector< std::vector <double> > &A)
	{
		double detA = GetDeterminantForMatrix3x3(A);
		//printf( "detA = %.5e\n", detA );
		std::vector< std::vector <double> > res(3);
		res[0].resize(3);
		res[1].resize(3);
		res[2].resize(3);

		res[0][0] = (A[1][1] * A[2][2] - A[1][2] * A[2][1]) / (A[2][0] * A[0][1] * A[1][2] - A[2][0] * A[0][2] * A[1][1] - A[1][0] * A[0][1] * A[2][2] + A[1][0] * A[0][2] * A[2][1] + A[0][0] * A[1][1] * A[2][2] - A[0][0] * A[1][2] * A[2][1]);
		res[0][1] = -(A[0][1] * A[2][2] - A[0][2] * A[2][1]) / (A[2][0] * A[0][1] * A[1][2] - A[2][0] * A[0][2] * A[1][1] - A[1][0] * A[0][1] * A[2][2] + A[1][0] * A[0][2] * A[2][1] + A[0][0] * A[1][1] * A[2][2] - A[0][0] * A[1][2] * A[2][1]);
		res[0][2] = (A[0][1] * A[1][2] - A[0][2] * A[1][1]) / (A[2][0] * A[0][1] * A[1][2] - A[2][0] * A[0][2] * A[1][1] - A[1][0] * A[0][1] * A[2][2] + A[1][0] * A[0][2] * A[2][1] + A[0][0] * A[1][1] * A[2][2] - A[0][0] * A[1][2] * A[2][1]);

		res[1][0] = -(-A[2][0] * A[1][2] + A[1][0] * A[2][2]) / (A[2][0] * A[0][1] * A[1][2] - A[2][0] * A[0][2] * A[1][1] - A[1][0] * A[0][1] * A[2][2] + A[1][0] * A[0][2] * A[2][1] + A[0][0] * A[1][1] * A[2][2] - A[0][0] * A[1][2] * A[2][1]);
		res[1][1] = (-A[2][0] * A[0][2] + A[0][0] * A[2][2]) / (A[2][0] * A[0][1] * A[1][2] - A[2][0] * A[0][2] * A[1][1] - A[1][0] * A[0][1] * A[2][2] + A[1][0] * A[0][2] * A[2][1] + A[0][0] * A[1][1] * A[2][2] - A[0][0] * A[1][2] * A[2][1]);
		res[1][2] = -(-A[1][0] * A[0][2] + A[0][0] * A[1][2]) / (A[2][0] * A[0][1] * A[1][2] - A[2][0] * A[0][2] * A[1][1] - A[1][0] * A[0][1] * A[2][2] + A[1][0] * A[0][2] * A[2][1] + A[0][0] * A[1][1] * A[2][2] - A[0][0] * A[1][2] * A[2][1]);

		res[2][0] = (-A[2][0] * A[1][1] + A[1][0] * A[2][1]) / (A[2][0] * A[0][1] * A[1][2] - A[2][0] * A[0][2] * A[1][1] - A[1][0] * A[0][1] * A[2][2] + A[1][0] * A[0][2] * A[2][1] + A[0][0] * A[1][1] * A[2][2] - A[0][0] * A[1][2] * A[2][1]);
		res[2][1] = -(-A[2][0] * A[0][1] + A[0][0] * A[2][1]) / (A[2][0] * A[0][1] * A[1][2] - A[2][0] * A[0][2] * A[1][1] - A[1][0] * A[0][1] * A[2][2] + A[1][0] * A[0][2] * A[2][1] + A[0][0] * A[1][1] * A[2][2] - A[0][0] * A[1][2] * A[2][1]);
		res[2][2] = (-A[1][0] * A[0][1] + A[0][0] * A[1][1]) / (A[2][0] * A[0][1] * A[1][2] - A[2][0] * A[0][2] * A[1][1] - A[1][0] * A[0][1] * A[2][2] + A[1][0] * A[0][2] * A[2][1] + A[0][0] * A[1][1] * A[2][2] - A[0][0] * A[1][2] * A[2][1]);

		/*for( int i = 0; i < 3; i++ )
			for( int j = 0; j < 3; j++ )
				res[i][j] /= detA;*/

		return res;
	}
	Tensor2Rank3D SolveInverseMatrix3x3(Tensor2Rank3D &A)
	{
		double detA = GetDeterminantForMatrix3x3(A);
		//printf( "detA = %.5e\n", detA );
		Tensor2Rank3D res;

		res.val[0][0] = (A.val[1][1] * A.val[2][2] - A.val[1][2] * A.val[2][1]) / (A.val[2][0] * A.val[0][1] * A.val[1][2] - A.val[2][0] * A.val[0][2] * A.val[1][1] - A.val[1][0] * A.val[0][1] * A.val[2][2] + A.val[1][0] * A.val[0][2] * A.val[2][1] + A.val[0][0] * A.val[1][1] * A.val[2][2] - A.val[0][0] * A.val[1][2] * A.val[2][1]);
		res.val[0][1] = -(A.val[0][1] * A.val[2][2] - A.val[0][2] * A.val[2][1]) / (A.val[2][0] * A.val[0][1] * A.val[1][2] - A.val[2][0] * A.val[0][2] * A.val[1][1] - A.val[1][0] * A.val[0][1] * A.val[2][2] + A.val[1][0] * A.val[0][2] * A.val[2][1] + A.val[0][0] * A.val[1][1] * A.val[2][2] - A.val[0][0] * A.val[1][2] * A.val[2][1]);
		res.val[0][2] = (A.val[0][1] * A.val[1][2] - A.val[0][2] * A.val[1][1]) / (A.val[2][0] * A.val[0][1] * A.val[1][2] - A.val[2][0] * A.val[0][2] * A.val[1][1] - A.val[1][0] * A.val[0][1] * A.val[2][2] + A.val[1][0] * A.val[0][2] * A.val[2][1] + A.val[0][0] * A.val[1][1] * A.val[2][2] - A.val[0][0] * A.val[1][2] * A.val[2][1]);
		   
		res.val[1][0] = -(-A.val[2][0] * A.val[1][2] + A.val[1][0] * A.val[2][2]) / (A.val[2][0] * A.val[0][1] * A.val[1][2] - A.val[2][0] * A.val[0][2] * A.val[1][1] - A.val[1][0] * A.val[0][1] * A.val[2][2] + A.val[1][0] * A.val[0][2] * A.val[2][1] + A.val[0][0] * A.val[1][1] * A.val[2][2] - A.val[0][0] * A.val[1][2] * A.val[2][1]);
		res.val[1][1] = (-A.val[2][0] * A.val[0][2] + A.val[0][0] * A.val[2][2]) / (A.val[2][0] * A.val[0][1] * A.val[1][2] - A.val[2][0] * A.val[0][2] * A.val[1][1] - A.val[1][0] * A.val[0][1] * A.val[2][2] + A.val[1][0] * A.val[0][2] * A.val[2][1] + A.val[0][0] * A.val[1][1] * A.val[2][2] - A.val[0][0] * A.val[1][2] * A.val[2][1]);
		res.val[1][2] = -(-A.val[1][0] * A.val[0][2] + A.val[0][0] * A.val[1][2]) / (A.val[2][0] * A.val[0][1] * A.val[1][2] - A.val[2][0] * A.val[0][2] * A.val[1][1] - A.val[1][0] * A.val[0][1] * A.val[2][2] + A.val[1][0] * A.val[0][2] * A.val[2][1] + A.val[0][0] * A.val[1][1] * A.val[2][2] - A.val[0][0] * A.val[1][2] * A.val[2][1]);
		   
		res.val[2][0] = (-A.val[2][0] * A.val[1][1] + A.val[1][0] * A.val[2][1]) / (A.val[2][0] * A.val[0][1] * A.val[1][2] - A.val[2][0] * A.val[0][2] * A.val[1][1] - A.val[1][0] * A.val[0][1] * A.val[2][2] + A.val[1][0] * A.val[0][2] * A.val[2][1] + A.val[0][0] * A.val[1][1] * A.val[2][2] - A.val[0][0] * A.val[1][2] * A.val[2][1]);
		res.val[2][1] = -(-A.val[2][0] * A.val[0][1] + A.val[0][0] * A.val[2][1]) / (A.val[2][0] * A.val[0][1] * A.val[1][2] - A.val[2][0] * A.val[0][2] * A.val[1][1] - A.val[1][0] * A.val[0][1] * A.val[2][2] + A.val[1][0] * A.val[0][2] * A.val[2][1] + A.val[0][0] * A.val[1][1] * A.val[2][2] - A.val[0][0] * A.val[1][2] * A.val[2][1]);
		res.val[2][2] = (-A.val[1][0] * A.val[0][1] + A.val[0][0] * A.val[1][1]) / (A.val[2][0] * A.val[0][1] * A.val[1][2] - A.val[2][0] * A.val[0][2] * A.val[1][1] - A.val[1][0] * A.val[0][1] * A.val[2][2] + A.val[1][0] * A.val[0][2] * A.val[2][1] + A.val[0][0] * A.val[1][1] * A.val[2][2] - A.val[0][0] * A.val[1][2] * A.val[2][1]);

		/*for( int i = 0; i < 3; i++ )
			for( int j = 0; j < 3; j++ )
				res[i][j] /= detA;*/

		return res;
	}

	
	//Copy ONLY!!! CSSD_matrix<Tensor2Rank3D, Point<double>> into CSSD_matrix<double, double>
	template <typename T1, typename T2>
	void MakeCopyMatrix_A_into_B(T1 &A, T2 &B)
	{
		//crate new columns id
		std::vector<std::vector<int>> id_column_for_B_up(A.id_column_for_A_up.size()*3);
		std::vector<std::vector<int>> id_column_for_B_down(A.id_column_for_A_down.size()*3);
		for (int id_row_A = 0; id_row_A < A.id_column_for_A_up.size(); id_row_A++)
		{
			for (int ii = 0; ii < 3; ii++)
			{
				//поправка на диагональ
				id_column_for_B_up[id_row_A * 3 + ii].resize(2 - ii + A.id_column_for_A_up[id_row_A].size()*3);
				for (int jj = 0; jj < 2 - ii; jj++)
				{
					id_column_for_B_up[id_row_A * 3 + ii][jj] = id_row_A * 3 + ii + jj + 1;
				}
				//остальные элементы строк верхнего тругольника
				for (int id_column_A = 0; id_column_A < A.id_column_for_A_up[id_row_A].size(); id_column_A++)
				{
					for (int jj = 0; jj < 3; jj++)
					{
						id_column_for_B_up[id_row_A * 3 + ii][id_column_A * 3 + jj + 2 - ii] = A.id_column_for_A_up[id_row_A][id_column_A]*3+jj;
					}
				}

				id_column_for_B_down[id_row_A * 3 + ii].resize(ii + A.id_column_for_A_down[id_row_A].size()*3);
				//остальные элементы строк нижнего тругольника
				for (int id_column_A = 0; id_column_A < A.id_column_for_A_down[id_row_A].size(); id_column_A++)
				{
					for (int jj = 0; jj < 3; jj++)
					{
						id_column_for_B_down[id_row_A * 3 + ii][id_column_A * 3 + jj] = A.id_column_for_A_down[id_row_A][id_column_A]*3+jj;
					}
				}
				//поправка на диагональ
				for (int jj = 0; jj < ii; jj++)
				{
					int id_j = id_column_for_B_down[id_row_A * 3 + ii].size() - ii + jj;
					id_column_for_B_down[id_row_A * 3 + ii][id_j] = id_row_A * 3 + jj;
				}
			}
		}

		B.Initialization(id_column_for_B_up, id_column_for_B_down);

		for (int id_row_A = 0; id_row_A < A.Diag.size(); id_row_A++)
		{
			for (int ii = 0; ii < 3; ii++)
			{
				B.Diag[id_row_A * 3 + ii] = A.Diag[id_row_A].val[ii][ii];
			}
			B.X[id_row_A * 3 + 0] = A.X[id_row_A].x;
			B.X[id_row_A * 3 + 1] = A.X[id_row_A].y;
			B.X[id_row_A * 3 + 2] = A.X[id_row_A].z;
			B.F[id_row_A * 3 + 0] = A.F[id_row_A].x;
			B.F[id_row_A * 3 + 1] = A.F[id_row_A].y;
			B.F[id_row_A * 3 + 2] = A.F[id_row_A].z;

			for (int ii = 0; ii < 3; ii++)
			{
				//поправка на диагональ
				for (int jj = 0; jj < 2 - ii; jj++)
				{
					B.A_up[id_row_A * 3 + ii][jj] = A.Diag[id_row_A].val[ii][ii+jj+1];
				}
				//остальные элементы строк верхнего тругольника
				for (int id_column_A = 0; id_column_A < A.id_column_for_A_up[id_row_A].size(); id_column_A++)
				{
					for (int jj = 0; jj < 3; jj++)
					{
						B.A_up[id_row_A * 3 + ii][id_column_A * 3 + jj + 2 - ii] = A.A_up[id_row_A][id_column_A].val[ii][jj];
					}
				}

				//остальные элементы строк нижнего тругольника
				for (int id_column_A = 0; id_column_A < A.id_column_for_A_down[id_row_A].size(); id_column_A++)
				{
					for (int jj = 0; jj < 3; jj++)
					{
						B.A_down[id_row_A * 3 + ii][id_column_A * 3 + jj] = A.A_down[id_row_A][id_column_A].val[ii][jj];
					}
				}
				//поправка на диагональ
				for (int jj = 0; jj < ii; jj++)
				{
					int id_j = id_column_for_B_down[id_row_A * 3 + ii].size() - ii + jj;
					id_column_for_B_down[id_row_A * 3 + ii][id_j] = id_row_A * 3 + jj;
					B.A_down[id_row_A * 3 + ii][id_j] = A.Diag[id_row_A].val[ii][jj];
				}
			}
		}
	}

	//-------------> geometry

	//A*p0+B*p1+C*p2+D=0
	void GetPlaneEquation(std::vector<Point<double>> P, double &A, double &B, double &C, double &D)
	{
		double Matr[3][3];
		Matr[0][0] = 1; Matr[0][1] = P[0].y; Matr[0][2] = P[0].z;
		Matr[1][0] = 1; Matr[1][1] = P[1].y; Matr[1][2] = P[1].z;
		Matr[2][0] = 1; Matr[2][1] = P[2].y; Matr[2][2] = P[2].z;
		A = GetDeterminantForMatrix3x3(Matr);
		Matr[0][0] = P[0].x; Matr[0][1] = 1; Matr[0][2] = P[0].z;
		Matr[1][0] = P[1].x; Matr[1][1] = 1; Matr[1][2] = P[1].z;
		Matr[2][0] = P[2].x; Matr[2][1] = 1; Matr[2][2] = P[2].z;
		B = GetDeterminantForMatrix3x3(Matr);
		Matr[0][0] = P[0].x; Matr[0][1] = P[0].y; Matr[0][2] = 1;
		Matr[1][0] = P[1].x; Matr[1][1] = P[1].y; Matr[1][2] = 1;
		Matr[2][0] = P[2].x; Matr[2][1] = P[2].y; Matr[2][2] = 1;
		C = GetDeterminantForMatrix3x3(Matr);
		Matr[0][0] = P[0].x; Matr[0][1] = P[0].y; Matr[0][2] = P[0].z;
		Matr[1][0] = P[1].x; Matr[1][1] = P[1].y; Matr[1][2] = P[1].z;
		Matr[2][0] = P[2].x; Matr[2][1] = P[2].y; Matr[2][2] = P[2].z;
		D = -GetDeterminantForMatrix3x3(Matr);
		return;
	}
	void GetPlaneEquation(Point<double> P[3], double &A, double &B, double &C, double &D)
	{
		double Matr[3][3];
		Matr[0][0] = 1; Matr[0][1] = P[0].y; Matr[0][2] = P[0].z;
		Matr[1][0] = 1; Matr[1][1] = P[1].y; Matr[1][2] = P[1].z;
		Matr[2][0] = 1; Matr[2][1] = P[2].y; Matr[2][2] = P[2].z;
		A = GetDeterminantForMatrix3x3(Matr);
		Matr[0][0] = P[0].x; Matr[0][1] = 1; Matr[0][2] = P[0].z;
		Matr[1][0] = P[1].x; Matr[1][1] = 1; Matr[1][2] = P[1].z;
		Matr[2][0] = P[2].x; Matr[2][1] = 1; Matr[2][2] = P[2].z;
		B = GetDeterminantForMatrix3x3(Matr);
		Matr[0][0] = P[0].x; Matr[0][1] = P[0].y; Matr[0][2] = 1;
		Matr[1][0] = P[1].x; Matr[1][1] = P[1].y; Matr[1][2] = 1;
		Matr[2][0] = P[2].x; Matr[2][1] = P[2].y; Matr[2][2] = 1;
		C = GetDeterminantForMatrix3x3(Matr);
		Matr[0][0] = P[0].x; Matr[0][1] = P[0].y; Matr[0][2] = P[0].z;
		Matr[1][0] = P[1].x; Matr[1][1] = P[1].y; Matr[1][2] = P[1].z;
		Matr[2][0] = P[2].x; Matr[2][1] = P[2].y; Matr[2][2] = P[2].z;
		D = -GetDeterminantForMatrix3x3(Matr);
		return;
	}
	//A*p0+B*p1+C*p2+D=0
	void GetPlaneEquation(Point<double> P0, Point<double> P1, Point<double> P2, double &A, double &B, double &C, double &D)
	{
		Point<double> P[3];
		P[0] = P0;
		P[1] = P1;
		P[2] = P2;
		GetPlaneEquation(P, A, B, C, D);
	}
	Point<double> GetIntersectionPointOfLineAndPlane(double A, double B, double C, double D, Point<double> line0, Point<double> line1)
	{
		Point<double> delta = line1 - line0;
		double t = (-D - A * line0.x - B * line0.y - C * line0.z) / (A*delta.x + B * delta.y + C * delta.z);
		return line0 + delta * t;
	}
	Point<double> GetIntersectionPointOfLineAndPlane(Point<double> P0, Point<double> P1, Point<double> P2, Point<double> line0, Point<double> line1)
	{
		double A, B, C, D;
		GetPlaneEquation(P0, P1, P2, A, B, C, D);
		return GetIntersectionPointOfLineAndPlane(A, B, C, D, line0, line1);
	}
	double SolveLengthVector(Point<double>  A, Point<double> B)
	{
		return sqrt((A.x - B.x)*(A.x - B.x) + (A.y - B.y)*(A.y - B.y) + (A.z - B.z)*(A.z - B.z));
	}
	double SolveLengthVector(Point<double>  A)
	{
		return sqrt((A.x)*(A.x) + (A.y)*(A.y) + (A.z)*(A.z));
	}
	Point<double> GetIntersectionPointOfLineAndLine(Point<double> line0_left, Point<double> line0_right, Point<double> line1_left, Point<double> line1_right)
	{
		Point<double> test_point;
		Point<double> n_line0 = line0_right - line0_left;
		Point<double> n_line1 = line1_right - line1_left;
		int i = 0;
		do
		{
			double rx = (rand() * 1.0 / RAND_MAX);
			double ry = (rand() * 1.0 / RAND_MAX);
			double rz = (rand() * 1.0 / RAND_MAX);
			test_point.x = line0_left.x + rx * (n_line0.x + rx * 2);
			test_point.y = line0_left.y + ry * (n_line0.y + ry * 2);
			test_point.z = line0_left.z + rz * (n_line0.z + rz * 2);

			Point<double> lambda;
			lambda = test_point - line0_left;
			lambda.x /= n_line0.x;
			lambda.y /= n_line0.y;
			lambda.z /= n_line0.z;
			if (IsEqual(n_line0.x, 0.0)) lambda.x = 0;
			if (IsEqual(n_line0.y, 0.0)) lambda.y = 0;
			if (IsEqual(n_line0.z, 0.0)) lambda.z = 0;

			if (!IsEqual(lambda.x, lambda.y) || !IsEqual(lambda.y, lambda.z) || !IsEqual(lambda.x, lambda.z))
			{
				Point<double> n0; double D0;
				Point<double> n1; double D1;
				GetPlaneEquation(test_point, line0_left, line0_right, n0.x, n0.y, n0.z, D0);
				GetPlaneEquation(line0_left, line1_left, line1_right, n1.x, n1.y, n1.z, D1);
				n0 /= SolveLengthVector(n0);
				n1 /= SolveLengthVector(n1);
				if (!IsEqual(abs(n0.x), abs(n1.x)) || !IsEqual(abs(n0.y), abs(n1.y)) || !IsEqual(abs(n0.z), abs(n1.z)))
				{
					break;
				}
			}
			i++;
		} while (true);

		Point<double> result = GetIntersectionPointOfLineAndPlane(test_point, line0_left, line0_right, line1_left, line1_right);
		return result;
	}
	bool IsLineContainThePoint(Point<double> line0, Point<double> line1, Point<double> test)
	{
		double len_left = SolveLengthVector(test, line0);
		double len_right = SolveLengthVector(test, line1);
		double len_full = SolveLengthVector(line0, line1);
		if (IsEqual(len_left + len_right, len_full))
		{
			return true;
		}
		return false;
	}
	double SolveSquareTriangle(Point<double> P1, Point<double> P2, Point<double> P3)
	{
		double S;
		double edge01_s = SolveLengthVector(P1, P2);
		double edge02_s = SolveLengthVector(P1, P3);
		double edge12_s = SolveLengthVector(P2, P3);
		double p_s = (edge01_s + edge02_s + edge12_s) / 2.;
		S = sqrt(p_s*(p_s - edge01_s)*(p_s - edge02_s)*(p_s - edge12_s));

		return S;
	}
	template<typename DenseMatrix_d_d>
	Point<double> MakeNormalForLine(Point<double> Line[2], Point<double> test_point, DenseMatrix_d_d&Matrix)
	{
		Point<double> Plane[3], _X;
		Plane[0] = Line[0];
		Plane[1] = Line[1];
		Plane[2] = test_point;
		double A, B, C, D;
		GetPlaneEquation(Plane, A, B, C, D);
		Point<double> AB = Line[1] - Line[0];
		Point<double> AD = test_point - Line[0];

		//DenseMatrix<double, double> Matrix;
		Matrix.SetSize(3);
		Matrix.SetElementOfA(0, 0, A);
		Matrix.SetElementOfA(0, 1, B);
		Matrix.SetElementOfA(0, 2, C);

		Matrix.SetElementOfA(1, 0, AB.x);
		Matrix.SetElementOfA(1, 1, AB.y);
		Matrix.SetElementOfA(1, 2, AB.z);

		Matrix.SetElementOfA(2, 0, AD.x);
		Matrix.SetElementOfA(2, 1, AD.y);
		Matrix.SetElementOfA(2, 2, AD.z);

		Matrix.SetF(0, 0);
		Matrix.SetF(1, 0);
		Matrix.SetF(2, -1);

		Matrix.Gauss();

		_X = Point<double>(Matrix.GetElementOfX(0), Matrix.GetElementOfX(1), Matrix.GetElementOfX(2));
		_X = _X / SolveLengthVector(_X);
		return _X;
	}
	//angle = [rad]
	double SolveAngleBetweenV1AndV2(Point<double> v1, Point<double> v2)
	{
		double cos_v1v2 = (v1*v2) / (SolveLengthVector(v1)*SolveLengthVector(v2));
		double angle = acos(cos_v1v2);
		return angle;
	}
	//Old->New - self_basis
	void MakeNewBasisOnPlane(Point<double> plane[3], std::vector<std::vector<double>> &self_basis, std::vector<std::vector<double>> &self_basis_revers)
	{
		struct Plane {
		public:
			Point<double> A, B, C;
			Plane()
			{}
			Plane(Point<double> A, Point<double> B, Point<double> C)
			{
				this->A = A;
				this->B = B;
				this->C = C;
			}
			bool IsContain(Point<double> X)
			{
				Point<double> tmp;
				tmp.x = math::GetRound(X.x, 6);
				tmp.y = math::GetRound(X.y, 6);
				tmp.z = math::GetRound(X.z, 6);
				double matr[3][3] = { { tmp.x - A.x, tmp.y - A.y, tmp.z - A.z },
									  { B.x - A.x, B.y - A.y, B.z - A.z },
									  { C.x - A.x, C.y - A.y, C.z - A.z } };
				double d = abs(math::GetDeterminantForMatrix3x3(matr));
				if (d <= 1E-7)
					return true;
				return false;
			}
			double SolveDistanceToPoint(Point<double> X)
			{
				double a = (B.y - A.y)*(C.z - A.z) - (C.y - A.y)*(B.z - A.z);
				double b = (B.z - A.z)*(C.x - A.x) - (C.z - A.z)*(B.x - A.x);
				double c = (B.x - A.x)*(C.y - A.y) - (C.x - A.x)*(B.y - A.y);
				double d = -A.x*(B.y - A.y)*(C.z - A.z) - A.y*(C.x - A.x)*(B.z - A.z) - A.z*(B.x - A.x)*(C.y - A.y)
					+ A.x*(C.y - A.y)*(B.z - A.z) + A.y*(B.x - A.x)*(C.z - A.z) + A.z*(C.x - A.x)*(B.y - A.y);

				return abs(a*X.x + b * X.y + c * X.z + d) / sqrt(a*a + b * b + c * c);
			}
			void GetCoefficients(double &a, double &b, double &c, double &d)
			{
				a = (B.y - A.y)*(C.z - A.z) - (C.y - A.y)*(B.z - A.z);
				b = (B.z - A.z)*(C.x - A.x) - (C.z - A.z)*(B.x - A.x);
				c = (B.x - A.x)*(C.y - A.y) - (C.x - A.x)*(B.y - A.y);
				d = -A.x*(B.y - A.y)*(C.z - A.z) - A.y*(C.x - A.x)*(B.z - A.z) - A.z*(B.x - A.x)*(C.y - A.y)
					+ A.x*(C.y - A.y)*(B.z - A.z) + A.y*(B.x - A.x)*(C.z - A.z) + A.z*(C.x - A.x)*(B.y - A.y);
			}
		};

		self_basis.resize(3);
		for (int i = 0; i < 3; i++)
			self_basis[i].resize(3);

		Point<double> X, Y, Z; //??? ???? ????? ???? (?? ????????? ?????? ????????? ? (0,0,0))

					   //?????? ???? X
		X = plane[1] - plane[0];
		X /= SolveLengthVector(X, Point<double>(0, 0, 0)); //??????????

		double t = SolveLengthVector(X, Point<double>(0, 0, 0));

		//?????? ???? Y ?? ??????????? ??????????????? ? ? ? ?????????????? ? ????????? ??????
		Plane tr(plane[0], plane[1], plane[2]);
		double _kx, _ky, _kz, _kd, kx, ky, kz, kd;
		tr.GetCoefficients(_kx, _ky, _kz, _kd);
		kx = _kx / sqrt(_kx*_kx + _ky * _ky + _kz * _kz);
		ky = _ky / sqrt(_kx*_kx + _ky * _ky + _kz * _kz);
		kz = _kz / sqrt(_kx*_kx + _ky * _ky + _kz * _kz);
		kd = _kd / sqrt(_kx*_kx + _ky * _ky + _kz * _kz);

		t = tr.SolveDistanceToPoint(Point<double>(0, 0, 0));

		if (X.z != 0)
		{
			if ((ky - X.y*kz / X.z) != 0)
			{
				Y.x = 1.0;
				Y.y = Y.x*(-kx + X.x*kz / X.z) / (ky - X.y*kz / X.z);
				Y.z = (-X.y*Y.y - X.x*Y.x) / X.z;
			}
			else {
				if ((kx - X.x*kz / X.z) != 0)
				{
					Y.y = 1.0;
					Y.x = Y.y*(-ky + X.y*kz / X.z) / (kx - X.x*kz / X.z);
					Y.z = (-X.y*Y.y - X.x*Y.x) / X.z;
				}
				else goto next_if1;
			}
		}
		else {
		next_if1:
			if (X.y != 0)
			{
				if ((kz - X.z*ky / X.y) != 0)
				{
					Y.x = 1.0;
					Y.z = Y.x*(-kx + X.x*ky / X.y) / (kz - X.z*ky / X.y);
					Y.y = (-X.z*Y.z - X.x*Y.x) / X.y;
				}
				else {
					if ((kx - X.x*ky / X.y) != 0)
					{
						Y.z = 1.0;
						Y.x = Y.z*(-kz + X.z*ky / X.y) / (kx - X.x*ky / X.y);
						Y.y = (-X.z*Y.z - X.x*Y.x) / X.y;
					}
					else goto next_if2;
				}
			}
			else {
			next_if2:
				if ((ky - X.y*kx / X.x) != 0)
				{
					Y.z = 1.0;
					Y.y = Y.z*(-kz + X.z*kx / X.x) / (ky - X.y*kx / X.x);
					Y.x = (-X.z*Y.z - X.y*Y.y) / X.x;
				}
				else {
					if ((kz - X.z*kx / X.x) != 0)
					{
						Y.y = 1.0;
						Y.z = Y.y*(-ky + X.y*kx / X.x) / (kz - X.z*kx / X.x);
						Y.x = (-X.z*Y.z - X.y*Y.y) / X.x;
					}
					else goto next_if2;
				}
			}
		}
		Y /= SolveLengthVector(Y, Point<double>(0, 0, 0)); //??????????

		t = SolveLengthVector(Y, Point<double>(0, 0, 0));

		//?????? ???? Z ?? ??????????? ??????????????? ? ? ? Y
		if (X.z != 0)
		{
			Z.x = 1.0;
			Z.y = (X.x*Y.z / X.z - Y.x)*Z.x / (Y.y - X.y*Y.z / X.z);
			Z.z = (-X.x*Z.x - X.y*Z.y) / X.z;

			if (Z.z != Z.z || Z.x != Z.x || Z.y != Z.y || IsEqual((Y.y - X.y*Y.z / X.z), 0.0))
			{
				Z.y = 1.0;
				Z.x = (X.y*Y.z / X.z - Y.y)*Z.y / (Y.x - X.x*Y.z / X.z);
				Z.z = (-X.x*Z.x - X.y*Z.y) / X.z;
			}
		}
		else {
			if (X.y != 0)
			{
				Z.x = 1.0;
				Z.z = (X.x*Y.y / X.y - Y.x)*Z.x / (Y.z - X.z*Y.y / X.y);
				Z.y = (-X.x*Z.x - X.z*Z.z) / X.y;

				if (Z.z != Z.z || Z.x != Z.x || Z.y != Z.y || IsEqual((Y.z - X.z*Y.y / X.y), 0.0))
				{
					Z.z = 1.0;
					Z.x = (X.z*Y.y / X.y - Y.z)*Z.z / (Y.x - X.x*Y.y / X.y);
					Z.y = (-X.x*Z.x - X.z*Z.z) / X.y;
				}
			}
			else {
				Z.z = 1.0;
				Z.y = (X.z*Y.x / X.x - Y.z) / (Y.y - X.y*Y.x / X.x);
				Z.x = (-X.z - X.y*Z.y) / X.x;

				if (Z.z != Z.z || Z.x != Z.x || Z.y != Z.y || IsEqual((Y.y - X.y*Y.x / X.x), 0.0))
				{
					Z.y = 1.0;
					Z.z = (X.y*Y.x / X.x - Y.y) / (Y.z - X.z*Y.x / X.x);
					Z.x = (-X.z - X.y*Z.y) / X.x;
				}
			}
		}
		Z /= SolveLengthVector(Z, Point<double>(0, 0, 0)); //??????????


		t = X * Y;
		t = X * Z;
		t = Z * Y;

		//????????
		if (!(X*Y <= 1E-6 && X*Z <= 1e-6 && Y*Z <= 1e-6))
		{
			printf_s("ERROR in self basis 2D!!\n");
			//Sleep(10000);
		}

		self_basis[0][0] = X.x;
		self_basis[0][1] = Y.x;
		self_basis[0][2] = Z.x;
		self_basis[1][0] = X.y;
		self_basis[1][1] = Y.y;
		self_basis[1][2] = Z.y;
		self_basis[2][0] = X.z;
		self_basis[2][1] = Y.z;
		self_basis[2][2] = Z.z;

		self_basis_revers.resize(3);
		self_basis_revers[0].resize(3);
		self_basis_revers[1].resize(3);
		self_basis_revers[2].resize(3);

		self_basis_revers[0][0] = (self_basis[1][1] * self_basis[2][2] - self_basis[1][2] * self_basis[2][1]) / (self_basis[2][0] * self_basis[0][1] * self_basis[1][2] - self_basis[2][0] * self_basis[0][2] * self_basis[1][1] - self_basis[1][0] * self_basis[0][1] * self_basis[2][2] + self_basis[1][0] * self_basis[0][2] * self_basis[2][1] + self_basis[0][0] * self_basis[1][1] * self_basis[2][2] - self_basis[0][0] * self_basis[1][2] * self_basis[2][1]);
		self_basis_revers[0][1] = -(self_basis[0][1] * self_basis[2][2] - self_basis[0][2] * self_basis[2][1]) / (self_basis[2][0] * self_basis[0][1] * self_basis[1][2] - self_basis[2][0] * self_basis[0][2] * self_basis[1][1] - self_basis[1][0] * self_basis[0][1] * self_basis[2][2] + self_basis[1][0] * self_basis[0][2] * self_basis[2][1] + self_basis[0][0] * self_basis[1][1] * self_basis[2][2] - self_basis[0][0] * self_basis[1][2] * self_basis[2][1]);
		self_basis_revers[0][2] = (self_basis[0][1] * self_basis[1][2] - self_basis[0][2] * self_basis[1][1]) / (self_basis[2][0] * self_basis[0][1] * self_basis[1][2] - self_basis[2][0] * self_basis[0][2] * self_basis[1][1] - self_basis[1][0] * self_basis[0][1] * self_basis[2][2] + self_basis[1][0] * self_basis[0][2] * self_basis[2][1] + self_basis[0][0] * self_basis[1][1] * self_basis[2][2] - self_basis[0][0] * self_basis[1][2] * self_basis[2][1]);

		self_basis_revers[1][0] = -(-self_basis[2][0] * self_basis[1][2] + self_basis[1][0] * self_basis[2][2]) / (self_basis[2][0] * self_basis[0][1] * self_basis[1][2] - self_basis[2][0] * self_basis[0][2] * self_basis[1][1] - self_basis[1][0] * self_basis[0][1] * self_basis[2][2] + self_basis[1][0] * self_basis[0][2] * self_basis[2][1] + self_basis[0][0] * self_basis[1][1] * self_basis[2][2] - self_basis[0][0] * self_basis[1][2] * self_basis[2][1]);
		self_basis_revers[1][1] = (-self_basis[2][0] * self_basis[0][2] + self_basis[0][0] * self_basis[2][2]) / (self_basis[2][0] * self_basis[0][1] * self_basis[1][2] - self_basis[2][0] * self_basis[0][2] * self_basis[1][1] - self_basis[1][0] * self_basis[0][1] * self_basis[2][2] + self_basis[1][0] * self_basis[0][2] * self_basis[2][1] + self_basis[0][0] * self_basis[1][1] * self_basis[2][2] - self_basis[0][0] * self_basis[1][2] * self_basis[2][1]);
		self_basis_revers[1][2] = -(-self_basis[1][0] * self_basis[0][2] + self_basis[0][0] * self_basis[1][2]) / (self_basis[2][0] * self_basis[0][1] * self_basis[1][2] - self_basis[2][0] * self_basis[0][2] * self_basis[1][1] - self_basis[1][0] * self_basis[0][1] * self_basis[2][2] + self_basis[1][0] * self_basis[0][2] * self_basis[2][1] + self_basis[0][0] * self_basis[1][1] * self_basis[2][2] - self_basis[0][0] * self_basis[1][2] * self_basis[2][1]);

		self_basis_revers[2][0] = (-self_basis[2][0] * self_basis[1][1] + self_basis[1][0] * self_basis[2][1]) / (self_basis[2][0] * self_basis[0][1] * self_basis[1][2] - self_basis[2][0] * self_basis[0][2] * self_basis[1][1] - self_basis[1][0] * self_basis[0][1] * self_basis[2][2] + self_basis[1][0] * self_basis[0][2] * self_basis[2][1] + self_basis[0][0] * self_basis[1][1] * self_basis[2][2] - self_basis[0][0] * self_basis[1][2] * self_basis[2][1]);
		self_basis_revers[2][1] = -(-self_basis[2][0] * self_basis[0][1] + self_basis[0][0] * self_basis[2][1]) / (self_basis[2][0] * self_basis[0][1] * self_basis[1][2] - self_basis[2][0] * self_basis[0][2] * self_basis[1][1] - self_basis[1][0] * self_basis[0][1] * self_basis[2][2] + self_basis[1][0] * self_basis[0][2] * self_basis[2][1] + self_basis[0][0] * self_basis[1][1] * self_basis[2][2] - self_basis[0][0] * self_basis[1][2] * self_basis[2][1]);
		self_basis_revers[2][2] = (-self_basis[1][0] * self_basis[0][1] + self_basis[0][0] * self_basis[1][1]) / (self_basis[2][0] * self_basis[0][1] * self_basis[1][2] - self_basis[2][0] * self_basis[0][2] * self_basis[1][1] - self_basis[1][0] * self_basis[0][1] * self_basis[2][2] + self_basis[1][0] * self_basis[0][2] * self_basis[2][1] + self_basis[0][0] * self_basis[1][1] * self_basis[2][2] - self_basis[0][0] * self_basis[1][2] * self_basis[2][1]);
	}
	Point<double> MakeProgectionOfPointIntoLine(Point<double> Line[2], Point<double> X)
	{
		Point<double> _X;
		std::vector<std::vector<double>> self_basis; //Old->New
		std::vector<std::vector<double>> self_basis_revers; //New->Old

		Point<double> Plane[3];
		Plane[0] = Line[0];
		Plane[1] = Line[1];
		Plane[2] = X;
		MakeNewBasisOnPlane(Plane, self_basis, self_basis_revers);

		Point<double> NewLine[2], NewX, New_X;
		//NewLine[0].x = self_basis_revers[0][0] * Line[0].x + self_basis_revers[0][1] * Line[0].y + self_basis_revers[0][2] * Line[0].z;
		NewLine[0].y = self_basis_revers[1][0] * Line[0].x + self_basis_revers[1][1] * Line[0].y + self_basis_revers[1][2] * Line[0].z;
		NewLine[0].z = self_basis_revers[2][0] * Line[0].x + self_basis_revers[2][1] * Line[0].y + self_basis_revers[2][2] * Line[0].z;

		//NewLine[1].x = self_basis_revers[0][0] * Line[1].x + self_basis_revers[0][1] * Line[1].y + self_basis_revers[0][2] * Line[1].z;
		//NewLine[1].y = self_basis_revers[1][0] * Line[1].x + self_basis_revers[1][1] * Line[1].y + self_basis_revers[1][2] * Line[1].z;
		//NewLine[1].z = self_basis_revers[2][0] * Line[1].x + self_basis_revers[2][1] * Line[1].y + self_basis_revers[2][2] * Line[1].z;

		NewX.x = self_basis_revers[0][0] * X.x + self_basis_revers[0][1] * X.y + self_basis_revers[0][2] * X.z;
		//NewX.y = self_basis_revers[1][0] * X.x + self_basis_revers[1][1] * X.y + self_basis_revers[1][2] * X.z;
		//NewX.z = self_basis_revers[2][0] * X.x + self_basis_revers[2][1] * X.y + self_basis_revers[2][2] * X.z;

		New_X.x = NewX.x;
		New_X.y = NewLine[0].y;
		New_X.z = NewLine[0].z;

		_X.x = self_basis[0][0] * New_X.x + self_basis[0][1] * New_X.y + self_basis[0][2] * New_X.z;
		_X.y = self_basis[1][0] * New_X.x + self_basis[1][1] * New_X.y + self_basis[1][2] * New_X.z;
		_X.z = self_basis[2][0] * New_X.x + self_basis[2][1] * New_X.y + self_basis[2][2] * New_X.z;

		return _X;
	}
	Point<double> MakeProgectionOfPointIntoLine(Point<double> Line1, Point<double> Line2, Point<double> X)
	{
		Point<double> _X;
		std::vector<std::vector<double>> self_basis; //Old->New
		std::vector<std::vector<double>> self_basis_revers; //New->Old

		Point<double> Plane[3];
		Plane[0] = Line1;
		Plane[1] = Line2;
		Plane[2] = X;
		MakeNewBasisOnPlane(Plane, self_basis, self_basis_revers);

		Point<double> NewLine[2], NewX, New_X;
		//NewLine[0].x = self_basis_revers[0][0] * Line[0].x + self_basis_revers[0][1] * Line[0].y + self_basis_revers[0][2] * Line[0].z;
		NewLine[0].y = self_basis_revers[1][0] * Line1.x + self_basis_revers[1][1] * Line1.y + self_basis_revers[1][2] * Line1.z;
		NewLine[0].z = self_basis_revers[2][0] * Line1.x + self_basis_revers[2][1] * Line1.y + self_basis_revers[2][2] * Line1.z;

		//NewLine[1].x = self_basis_revers[0][0] * Line[1].x + self_basis_revers[0][1] * Line[1].y + self_basis_revers[0][2] * Line[1].z;
		//NewLine[1].y = self_basis_revers[1][0] * Line[1].x + self_basis_revers[1][1] * Line[1].y + self_basis_revers[1][2] * Line[1].z;
		//NewLine[1].z = self_basis_revers[2][0] * Line[1].x + self_basis_revers[2][1] * Line[1].y + self_basis_revers[2][2] * Line[1].z;

		NewX.x = self_basis_revers[0][0] * X.x + self_basis_revers[0][1] * X.y + self_basis_revers[0][2] * X.z;
		//NewX.y = self_basis_revers[1][0] * X.x + self_basis_revers[1][1] * X.y + self_basis_revers[1][2] * X.z;
		//NewX.z = self_basis_revers[2][0] * X.x + self_basis_revers[2][1] * X.y + self_basis_revers[2][2] * X.z;

		New_X.x = NewX.x;
		New_X.y = NewLine[0].y;
		New_X.z = NewLine[0].z;

		_X.x = self_basis[0][0] * New_X.x + self_basis[0][1] * New_X.y + self_basis[0][2] * New_X.z;
		_X.y = self_basis[1][0] * New_X.x + self_basis[1][1] * New_X.y + self_basis[1][2] * New_X.z;
		_X.z = self_basis[2][0] * New_X.x + self_basis[2][1] * New_X.y + self_basis[2][2] * New_X.z;

		return _X;
	}
	double SolveLengthPointAndLine(Point<double> Line1, Point<double> Line2, Point<double> X)
	{
		Point<double> O = MakeProgectionOfPointIntoLine(Line1, Line2, X);
		return SolveLengthVector(O, X);
	}
	template<typename DenseMatrix_d_d>
	Point<double> MakeProgectionOfPointIntoPlane(Point<double> P[3], Point<double> X, DenseMatrix_d_d & M)
	{
		Point<double> _X;
		double A, B, C, D;
		GetPlaneEquation(P, A, B, C, D);
		M.SetSize(3);
		/*Matrix.set_A(0, 0, B);
		Matrix.set_A(0, 1, -A);
		Matrix.set_A(0, 2, 0);

		Matrix.set_A(1, 0, -C);
		Matrix.set_A(1, 1, C);
		Matrix.set_A(1, 2, -(B-A));

		Matrix.set_A(2, 0, A);
		Matrix.set_A(2, 1, B);
		Matrix.set_A(2, 2, C);

		Matrix.set_b(0, X.x*B - X.y*A);
		Matrix.set_b(1, X.y*C - X.x*C - X.z*(B-A));
		Matrix.set_b(2, -D);*/

		Point<double> a = P[0];
		Point<double> b = P[1];
		Point<double> c = P[2];
		Point<double> s = X;
		M.SetElementOfA(0, 0, (b.y - a.y) * (c.z - a.z) - (b.z - a.z) * (c.y - a.y));
		M.SetElementOfA(0, 1, (b.z - a.z) * (c.x - a.x) - (b.x - a.x) * (c.z - a.z));
		M.SetElementOfA(0, 2, (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x));

		M.SetElementOfA(1, 0, b.x - a.x);
		M.SetElementOfA(1, 1, b.y - a.y);
		M.SetElementOfA(1, 2, b.z - a.z);

		M.SetElementOfA(2, 0, c.x - a.x);
		M.SetElementOfA(2, 1, c.y - a.y);
		M.SetElementOfA(2, 2, c.z - a.z);

		M.SetF(0, a.x * M.GetElementOfA(0, 0) + a.y * M.GetElementOfA(0, 1) + a.z * a.x * M.GetElementOfA(0, 2));
		M.SetF(1, s.x * M.GetElementOfA(1, 0) + s.y * M.GetElementOfA(1, 1) + s.z * a.x * M.GetElementOfA(1, 2));
		M.SetF(2, s.x * M.GetElementOfA(2, 0) + s.y * M.GetElementOfA(2, 1) + s.z * a.x * M.GetElementOfA(2, 2));

		M.Gauss();

		_X = Point<double>(M.GetElementOfX(0), M.GetElementOfX(1), M.GetElementOfX(2));
		return _X;
	}
	//newT=A*T;
	Point<double> MakeTransferPoint(Point<double> &T, Tensor2Rank3D &A)
	{
		Point<double> _newT;
		_newT.x = T.x*A.val[0][0] + T.y*A.val[0][1] + T.z*A.val[0][2];
		_newT.y = T.x*A.val[1][0] + T.y*A.val[1][1] + T.z*A.val[1][2];
		_newT.z = T.x*A.val[2][0] + T.y*A.val[2][1] + T.z*A.val[2][2];
		return _newT;
	}
	//newT=A*T+A*Offset;
	Point<double> MakeTransferPoint(Point<double> &T, Tensor2Rank3D &A, Point<double> &Offset)
	{
			Point<double> _newT, _newOffset;
			_newT = MakeTransferPoint(T, A);
			_newOffset = MakeTransferPoint(Offset, A);
			return _newT + _newOffset;
	}
	//T in (xyz), new T in (xyz)*, (xyz)*=A(xyz): newT[ij]=A[ik]T[kj]
	Tensor2Rank3D MakeTransferTensor(Tensor2Rank3D T, Tensor2Rank3D A)
	{
		Tensor2Rank3D newT;

		for (int I = 0; I < 3; I++)
		{
			for (int i = 0; i < 3; i++)
			{
				for (int j = 0; j < 3; j++)
				{
					newT.val[i][I] += A.val[i][j] * T.val[j][I];
				}
			}
		}
		return newT;
	}
	Tensor2Rank3D MakeTransferTensor(Tensor2Rank3D T, Tensor2Rank3D A, Point<double> &Offset)
	{
		Tensor2Rank3D newT;
		std::vector<double> Off(3);
		Off[0] = Offset.x;
		Off[1] = Offset.y;
		Off[2] = Offset.z;

		for (int I = 0; I < 3; I++)
		{
			for (int i = 0; i < 3; i++)
			{
				for (int j = 0; j < 3; j++)
				{
					newT.val[i][I] += A.val[i][j] * (T.val[j][I] + Off[j]);
				}
			}
		}
		return newT;
	}
	std::vector<std::vector<std::vector<double>>> MakeTransferTensor(std::vector<std::vector<std::vector<double>>> T, Tensor2Rank3D A)
	{
		std::vector<std::vector<std::vector<double>>> newT(3);
		for (int i = 0; i < newT.size(); i++)
		{
			newT[i].resize(3);
			for (int j = 0; j < newT[i].size(); j++)
			{
				newT[i][j].resize(3);
			}
		}

		for (int k = 0; k < 3; k++)
		{
			for (int I = 0; I < 3; I++)
			{
				for (int i = 0; i < 3; i++)
				{
					for (int j = 0; j < 3; j++)
					{
						newT[i][I][k] += A.val[i][j] * T[j][I][k];
					}
				}
			}
		}
		return newT;
	}

	//------------->input|output
	
	//A is string with end '\0' 
	void ParserStringToVectorFloat(char *A, std::vector<float> &B, const char *separator)
	{
		std::vector<float> v_el;
		std::vector<float>(v_el).swap(B);

		printf_s("%s", separator);

		char b[21];
		int curr_ib = 0;
		printf_s("strlen(%s) = %d, strlen(separator) = %d\r", A, strlen(A), strlen(separator));
		for (int ia = 0; ia < strlen(A); ia++)
		{
			bool flag = false;
			for (int js = 0; js < strlen(separator); js++)
			{
				if (separator[js] == A[ia])
				{
					flag = true;
					b[curr_ib] = '\0';

					if (curr_ib != 0)
					{
						B.push_back(atof(b));
					}
					curr_ib = 0;
					break;
				}
			}

			if (!flag)
			{
				b[curr_ib] = A[ia];
				curr_ib++;

				if (ia == strlen(A) - 1)
				{
					b[curr_ib + 1] = '\0';
					B.push_back(atof(b));
				}
			}
		}
	}
	void ParserStringToVectorInt(char *A, std::vector<int> &B, const char *separator)
	{
		std::vector<int> v_el;
		std::vector<int>(v_el).swap(B);

		printf_s("%s", separator);

		char b[21];
		int curr_ib = 0;
		printf_s("strlen(%s) = %d, strlen(separator) = %d\r", A, strlen(A), strlen(separator));
		for (int ia = 0; ia < strlen(A); ia++)
		{
			bool flag = false;
			for (int js = 0; js < strlen(separator); js++)
			{
				if (separator[js] == A[ia])
				{
					flag = true;
					b[curr_ib] = '\0';

					if (curr_ib != 0)
					{
						B.push_back(atoi(b));
					}
					curr_ib = 0;
					break;
				}
			}

			if (!flag)
			{
				b[curr_ib] = A[ia];
				curr_ib++;

				if (ia == strlen(A) - 1)
				{
					b[curr_ib + 1] = '\0';
					B.push_back(atoi(b));
				}
			}
		}
	}
	int ParserCharToInt(char A)
	{
		char B[2];
		B[0] = A;
		B[1] = '\0';
		return atoi(B);
	}
	wchar_t* MakeConvertationFromCharToWchar(char* c, size_t max)
	{
		wchar_t* w = new wchar_t[max];
		mbstowcs_s(&max, w, max, c, max);
		return w;
	}

	void ReadNonEmptyLine(FILE *fin, char *result_line)
	{
		char stop_sumbols[100] = "\n#";
		while (feof(fin) == 0)
		{
			char line[1000];
			fgets(line, 1000, fin);
			int k = 0;
			bool find_stop = false;
			for (int i = 0; i < strlen(line) && find_stop == false; i++)
			{
				for (int j = 0; j < strlen(stop_sumbols); j++)
				{
					if (line[i] == stop_sumbols[j])
					{
						find_stop = true;
						break;
					}
				}

				if (!find_stop && line[i] != ' ' && line[i] != '\t')
				{
					result_line[k] = line[i];
					k++;
				}
			}
			if (k != 0)
			{
				result_line[k] = '\0';
				break;
			}
		}
		return;
	}
	void ReadNonEmptyLine_forNumbers(FILE *fin, char *result_line)
	{
		char stop_sumbols[100] = "\n#";
		while (feof(fin) == 0)
		{
			char line[1000];
			fgets(line, 1000, fin);
			int k = 0;
			bool find_stop = false;
			for (int i = 0; i < strlen(line) && !find_stop; i++)
			{
				for (int j = 0; j < strlen(stop_sumbols); j++)
				{
					if (line[i] == stop_sumbols[j])
					{
						find_stop = true;
						break;
					}
				}

				if (!find_stop && line[i] != '\t')
				{
					result_line[k] = line[i];
					k++;
				}
			}
			if (k != 0)
			{
				result_line[k] = '\0';
				break;
			}
		}

		printf_s("%s\n", result_line);
		return;
	}


	Point<double> SolveSqrt(Point<double> a)
	{
		return Point<double>(sqrt(a.x), sqrt(a.y), sqrt(a.z));
	}
	double SolveSqrt(double a)
	{
		return sqrt(a);
	}

}