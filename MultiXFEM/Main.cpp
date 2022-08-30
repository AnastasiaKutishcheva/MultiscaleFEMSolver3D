#pragma once
#include "MultiXFEM.h"
#include <Windows.h>
#include <iostream>

//void Solve_StatiñElasticProblem_withCracks()
//{
//	bool is_print_logFile = false; //input
//								   //char* problem_directory = { "box_200x200x200/800el" }; //input
//	char problem_directory[1000] = { "box_200x200x200/90el" }; //input
//															   //char* problem_directory = { "box_200x200x200/28el" }; //input
//															   //char problem_directory[1000] = { "box_200x200x200/32kel" }; //input
//															   //char problem_directory[1000] = { "box_200x200x200/246kel" }; //input
//	char cracks_directory[1000]; //input
//	sprintf_s(cracks_directory, sizeof(cracks_directory), "%s/Crack_100x200_Z=50", problem_directory);
//	sprintf_s(cracks_directory, sizeof(cracks_directory), "%s/2Crack_100x200_Z=50_Z=150", problem_directory);
//
//	math::SimpleGrid geo_grid; //input
//	std::vector<std::vector<int>> first_boundary;//input
//	std::vector<std::vector<int>> second_boundary;//input
//	MultiXFEM::Grid solver_grid; //output
//	std::vector<Point<double>> Solution; //output
//
//										 //geo_grid.ReadFromSalomeDat("box_200x200x200/800el/Mesh_1.dat", 3);
//	char name_grid[1000];
//	sprintf_s(name_grid, sizeof(name_grid), "%s/Mesh_1.dat", problem_directory);
//	geo_grid.ReadFromSalomeDat(name_grid, 3);
//	for (int i = 0; i < geo_grid.xyz.size(); i++)
//	{
//		if (math::IsEqual(geo_grid.xyz[i].z, 200.0))
//		{
//			std::vector<int> vertex(2);
//			vertex[0] = 1; //type boundary
//			vertex[1] = i;
//			first_boundary.push_back(vertex);
//		}
//		if (math::IsEqual(geo_grid.xyz[i].z, 0.0))
//		{
//			std::vector<int> vertex(2);
//			vertex[0] = 0; //type boundary
//			vertex[1] = i;
//			first_boundary.push_back(vertex);
//		}
//	}
//
//	solver_grid.AddDomain();
//	auto domain = solver_grid.GetDomain(0);
//	domain->forMech.SetE(300 * 1e+9);
//	domain->forMech.SetV(0.27);
//
//	MultiXFEM::MultiXFEM_forElasticDeformation(
//		is_print_logFile,
//		cracks_directory, //input
//		geo_grid, //input
//		first_boundary,//input
//		second_boundary,//input
//		problem_directory,
//		solver_grid, //output
//		Solution //output
//	);
//}

void Solve_StatiñElasticProblem_withCracks()
{
	MultiXFEM::Grid solver_grid; //output
	std::vector<Point<double>> Solution; //output


	char properties_file[1000] = { "E:/Box/100x50x200/XFEM/50k/param_for_solver.txt" };
	//char properties_file[1000] = { "D:/Documents/universiti/MySOFT/TEST/box_200x200x200_crañks/800el_disc/base_properties.txt" };
	//char properties_file[1000] = { "base_properties.txt" };
	printf_s("The properties file contains the information:\n==========================================\n");
	printf_s("\tDo you need to print to a log file? 0-false/1-true\n");
	printf_s("\tinput: Mesh directory: box_200x200x200/90el\n");
	printf_s("\tinput: Crack directory: box_200x200x200/90el/Crack_100x200_Z=50\n");
	printf_s("\t\t<Step size for crack prolongation>\n");
	printf_s("\t\t<Points for crack smoothing>\n");
	printf_s("\t\t<Critical residual value for SLAE solution>\n");
	printf_s("\t<properties: Boundary values>\n");
	printf_s("\t<properties: Materials>\n");
	printf_s("\t<output: Result directory: box_200x200x200/90el/result>\n==========================================\n");


	printf_s("Enter the name of the properties file: ");
	scanf_s("%s", &properties_file);

	FILE *f_properties;
	fopen_s(&f_properties, properties_file, "r");
	if (f_properties == NULL)
	{
		printf_s("\nError in properties file\n");
	}
	bool is_print_logFile = false;
	char mesh_directory[1000];
	char cracks_directory[1000];
	char base_result_directory[1000];


	int _flag;
	fscanf_s(f_properties, "%d", &_flag);
	if (_flag == 1) is_print_logFile = true;

	math::ReadNonEmptyLine(f_properties, mesh_directory);
	math::SimpleGrid geo_grid; //input
	geo_grid.ReadFromNVTR(mesh_directory, 4);

	math::ReadNonEmptyLine(f_properties, cracks_directory);

	double step_size;
	{
		char _line[1000];
		math::ReadNonEmptyLine_forNumbers(f_properties, _line);
		std::vector<float> val;
		math::ParserStringToVectorFloat(_line, val, " ");
		step_size = val[0];
	}
	int crack_smoothing_coefficient;
	{
		char _line[1000];
		math::ReadNonEmptyLine_forNumbers(f_properties, _line);
		std::vector<int> val;
		math::ParserStringToVectorInt(_line, val, " ");
		crack_smoothing_coefficient = val[0];
	}
	double critical_residual;
	{
		char _line[1000];
		math::ReadNonEmptyLine_forNumbers(f_properties, _line);
		std::vector<float> val;
		math::ParserStringToVectorFloat(_line, val, " ");
		critical_residual = val[0];
	}

	//read Boundary values and vertexes
	struct _Dirichlet {
		//Point<double> value;
		//Point<bool> is_condition;
		std::vector<int> id_vertexes;
		std::function<Point<double>(Point<bool>&)> value;
	};
	std::vector<_Dirichlet> first_boundaries;
	{
		int Nb;
		{
			char _line[1000];
			math::ReadNonEmptyLine_forNumbers(f_properties, _line);
			std::vector<int> val;
			math::ParserStringToVectorInt(_line, val, " ");
			Nb = val[0];
		}
		first_boundaries.resize(Nb);
		char line[1000];
		for (int i = 0; i < Nb; i++)
		{
			Point<double> value;
			Point<bool> is_condition;

			math::ReadNonEmptyLine_forNumbers(f_properties, line);
			int jj = 0;
			int curr_i = 0;
			char _tmp_line[1000];
			for (int j = 0; j < 1000; j++)
			{
				switch (line[j])
				{
				case '\0': j = 1000; break;
				case '\n': j = 1000; break;
				case '#': j = 1000; break;
				case '*':
				{
					switch (jj)
					{
					case 0: is_condition.x = false; break;
					case 1: is_condition.y = false; break;
					case 2: is_condition.z = false; break;
					default:
						break;
					}
					jj++;
					break;
				}
				case ' ':
				{
					if (curr_i > 0)
					{
						std::vector<float> _val;
						_tmp_line[curr_i] = '\0';
						math::ParserStringToVectorFloat(_tmp_line, _val, " *");
						curr_i = 0;

						switch (jj)
						{
						case 0: is_condition.x = true; value.x = _val[0]; break;
						case 1: is_condition.y = true; value.y = _val[0]; break;
						case 2: is_condition.z = true; value.z = _val[0]; break;
						default:
							break;
						}
						jj++;
					}
					break;
				}
				case '\t': break;
				default:
				{
					_tmp_line[curr_i] = line[j];
					curr_i++;
				}
				break;
				}
			}

			first_boundaries[i].value = [value, is_condition](Point<bool> &is_take) -> Point<double>
			{
				is_take = is_condition;
				return value;
			};
		}

		if (Nb != 0)
		{
			char name_boundary_file[1000];
			sprintf_s(name_boundary_file, sizeof(name_boundary_file), "%s/kraev1.txt", mesh_directory);
			FILE *fbv;
			fopen_s(&fbv, name_boundary_file, "r");
			int N;
			fscanf_s(fbv, "%d", &N);
			for (int i = 0; i < N; i++)
			{
				int _type, _vertex;
				fscanf_s(fbv, "%d %d", &_vertex, &_type);
				first_boundaries[_type].id_vertexes.push_back(_vertex);
			}
			if (geo_grid.read_from_zero == false)
			{
				for (int t = 0; t < first_boundaries.size(); t++)
				{
					for (int i = 0; i < first_boundaries[t].id_vertexes.size(); i++)
					{
						first_boundaries[t].id_vertexes[i]--;
					}
				}
			}
			fclose(fbv);
		}
	}

	struct _Neumann {
		/*double value;
		Point<double> vector;*/
		std::function<Point<double>(Point<double>)> value;
		std::vector<std::function<Point<double>(Point<double>)>> values;
		std::vector<std::vector<int>> id_vertexes_as_triangle;
		std::vector<int> id_base_element;
	};
	std::vector<_Neumann> second_boundary;
	{
		int Nb;
		{
			char _line[1000];
			math::ReadNonEmptyLine_forNumbers(f_properties, _line);
			std::vector<int> val;
			math::ParserStringToVectorInt(_line, val, " ");
			Nb = val[0];
		}
		second_boundary.resize(Nb);
		char line[1000];
		std::vector<bool> is_individual_values(Nb);
		std::vector<double> individual_values(Nb);
		for (int i = 0; i < Nb; i++)
		{
			is_individual_values[i] = false;
			double value;
			Point<double> vector;

			math::ReadNonEmptyLine_forNumbers(f_properties, line);
			int jj = 0;
			int curr_i = 0;
			char _tmp_line[1000];
			for (int j = 0; j < 1000; j++)
			{
				if ((curr_i > 0) && (line[j] == '\n' || line[j] == '\0' || line[j] == '#'))
				{
					std::vector<float> val;
					_tmp_line[curr_i] = '\0';
					math::ParserStringToVectorFloat(_tmp_line, val, " ");
					curr_i = 0;

					value = val[0];
					vector.x = val[1];
					vector.y = val[2];
					vector.z = val[3];

					individual_values[i] = value;
					if (math::IsEqual(math::SolveLengthVector(vector), 0.0))
					{
						is_individual_values[i] = true;
					}

					break;
				}
				else
				{
					if (line[j] != '\t')
					{
						_tmp_line[curr_i] = line[j];
						curr_i++;
					}
				}
			}

			if (is_individual_values[i] == false)
			{
				second_boundary[i].value = [value, vector](Point<double> X) -> Point<double>
				{
					Point <double> res;
					res.x = vector.x *value;
					res.y = vector.y *value;
					res.z = vector.z *value;
					return res;
				};
			}
		}

		if (Nb != 0)
		{
			std::vector<std::vector<int>> id_vertexes;
			id_vertexes.resize(Nb);
			char name_boundary_file[1000];
			sprintf_s(name_boundary_file, sizeof(name_boundary_file), "%s/kraev2.txt", mesh_directory);
			FILE *fbv;
			fopen_s(&fbv, name_boundary_file, "r");
			int N;
			fscanf_s(fbv, "%d", &N);
			for (int i = 0; i < N; i++)
			{
				int _type, _vertex;
				fscanf_s(fbv, "%d %d", &_vertex, &_type);
				id_vertexes[_type].push_back(_vertex);
			}
			if (geo_grid.read_from_zero == false)
			{
				for (int t = 0; t < id_vertexes.size(); t++)
				{
					for (int i = 0; i < id_vertexes[t].size(); i++)
					{
						id_vertexes[t][i]--;
					}
				}
			}
			fclose(fbv);

			for (int id_type = 0; id_type < Nb; id_type++)
			{
				for (int id_element = 0; id_element < geo_grid.nvtr.size(); id_element++)
				{
					auto _tmp = math::GetConfluence(geo_grid.nvtr[id_element], id_vertexes[id_type]);
					if (_tmp.size() == 3)
					{
						//second_boundary[id_type].id_vertexes_as_triangle.push_back(_tmp);
						//second_boundary[id_type].id_base_element.push_back(id_element);

						int test_vertex = -1;
						for (int t = 0; t < geo_grid.nvtr[id_element].size(); t++)
						{
							bool find_id = false;
							for (int tt = 0; tt < _tmp.size(); tt++)
							{
								if (geo_grid.nvtr[id_element][t] == _tmp[tt])
								{
									find_id = true;
									break;
								}
							}
							if (!find_id)
							{
								test_vertex = geo_grid.nvtr[id_element][t];
								break;
							}
						}

						double A, B, C, D;
						Point<double> test_vector = geo_grid.xyz[test_vertex];
						std::vector<Point<double>> vertexes(3);
						vertexes[0] = geo_grid.xyz[_tmp[0]];
						vertexes[1] = geo_grid.xyz[_tmp[1]];
						vertexes[2] = geo_grid.xyz[_tmp[2]];
						bool reverse = false;
						math::GetPlaneEquation(vertexes, A, B, C, D);

						if (Point<double>(A, B, C)*test_vector > 0)
						{
							int t = _tmp[1];
							_tmp[1] = _tmp[2];
							_tmp[2] = t;

							vertexes[0] = geo_grid.xyz[_tmp[0]];
							vertexes[1] = geo_grid.xyz[_tmp[1]];
							vertexes[2] = geo_grid.xyz[_tmp[2]];
						}
						Point<double> normal;
						double d;
						math::GetPlaneEquation(vertexes, normal.x, normal.y, normal.z, d);
						normal /= math::SolveLengthVector(normal);

						second_boundary[id_type].id_vertexes_as_triangle.push_back(_tmp);
						second_boundary[id_type].id_base_element.push_back(id_element);

						if (is_individual_values[id_type] == true)
						{
							double _val = individual_values[id_type];
							std::function<Point<double>(Point<double>)> curr_val = [_val, normal](Point<double> X) -> Point<double>
							{
								Point <double> res;
								res.x = normal.x * _val;
								res.y = normal.y * _val;
								res.z = normal.z * _val;
								return res;
							};
							second_boundary[id_type].values.push_back(curr_val);
						}
					}
				}
			}
		}
	}

	//read Materials
	struct _material
	{
		double _E;
		double _v;

		_material(double _E, double _v)
		{
			this->_E = _E;
			this->_v = _v;
		}
	};
	std::vector<_material> _materials;
	int N_domain;
	{
		char _line[1000];
		math::ReadNonEmptyLine_forNumbers(f_properties, _line);
		std::vector<int> val;
		math::ParserStringToVectorInt(_line, val, " ");
		N_domain = val[0];
	}
	for (int id_domain = 0; id_domain < N_domain; id_domain++)
	{
		double _E, _v;
		{
			char _line[1000];
			math::ReadNonEmptyLine_forNumbers(f_properties, _line);
			std::vector<float> val;
			math::ParserStringToVectorFloat(_line, val, " ");
			_E = val[0];
			_v = val[1];
		}
		_materials.push_back(_material(_E, _v));
	}

	math::ReadNonEmptyLine(f_properties, base_result_directory);
	fclose(f_properties);

	CreateDirectory((LPCTSTR)base_result_directory, NULL);


	bool res = false;
	{
		printf_s("\n=============================*===============================\n");
		printf_s("================= Start solution================\n");

		for (int i = 0; i < _materials.size(); i++)
		{
			solver_grid.AddDomain();
			auto domain = solver_grid.GetDomain(i);
			domain->forMech.SetE(_materials[i]._E);
			domain->forMech.SetV(_materials[i]._v);
		}

		MultiXFEM::MultiXFEM_forElasticDeformation(
			is_print_logFile,
			critical_residual,
			cracks_directory, //input
			step_size,
			crack_smoothing_coefficient,
			geo_grid, //input
			first_boundaries, //input
			second_boundary, //input
			base_result_directory, //output
			solver_grid, //output
			Solution //output
		);

		//output solution
		printf_s("Print the mech result into .dat file... \n");
		{
			FILE *fout_tech;
			char name_u_tech[5000];
			sprintf_s(name_u_tech, "%s/U.dat", base_result_directory);
			fopen_s(&fout_tech, name_u_tech, "w");
			char name_in_file[1000];
			sprintf_s(name_in_file, "Elastic");
			std::vector<std::vector<char>> name_value(6);
			char name_v_tmp[6][100];
			sprintf_s(name_v_tmp[0], "sigma_xx");
			sprintf_s(name_v_tmp[1], "sigma_yy");
			sprintf_s(name_v_tmp[2], "sigma_zz");
			sprintf_s(name_v_tmp[3], "sigma_inv");
			sprintf_s(name_v_tmp[4], "eps_inv");
			sprintf_s(name_v_tmp[5], "material");
			for (int i = 0; i < name_value.size(); i++)
			{
				name_value[i].resize(100);
				for (int j = 0; j < name_value[i].size(); j++)
				{
					name_value[i][j] = name_v_tmp[i][j];
				}
			}
			std::vector<std::vector<double>> value(3 * 2);
			value[0].resize(solver_grid.GetElementsCount());
			value[1].resize(solver_grid.GetElementsCount());
			value[2].resize(solver_grid.GetElementsCount());
			value[3].resize(solver_grid.GetElementsCount());
			value[4].resize(solver_grid.GetElementsCount());
			value[5].resize(solver_grid.GetElementsCount());
			double sigma_inv_max = 0;
			int elem_sigma_max = 0;
			for (int i = 0; i < solver_grid.GetElementsCount(); i++)
			{
				auto element = solver_grid.GetElement(i);

				Point<double> Centr = element->GetWeightCentr();

				Point<double> U = solver_grid.GetSolutionInPoint(i, Centr, Solution);
				Point<Point<double>> dU = solver_grid.GetDerevativeFromSolutionInPoint(i, Centr, Solution);

				double eps[6] = { dU.x.x, dU.y.y, dU.z.z, dU.x.y + dU.y.x, dU.y.z + dU.z.y, dU.x.z + dU.z.x };
				double v = solver_grid.GetDomain(0)->forMech.v;
				double E = solver_grid.GetDomain(0)->forMech.GetE(0);
				double a = v / (1 - v);
				double b = (1 - 2 * v) / (2 * (1 - v));
				//double k = E * (1 - v) / ((1 + v)*(1 - v));
				double k = E * (1. - v) / (1. + v) / (1. - 2. * v);
				double sigma[6] = { k*(eps[0] * 1 + eps[1] * a + eps[2] * a), k*(eps[0] * a + eps[1] * 1 + eps[2] * a), k*(eps[0] * a + eps[1] * a + eps[2] * 1),
					k*(2 * eps[3] * b), k*(2 * eps[4] * b), k*(2 * eps[5] * b) };
				double sigma_inv = sqrt((sigma[0] - sigma[1])*(sigma[0] - sigma[1]) + (sigma[1] - sigma[2])*(sigma[1] - sigma[2]) + (sigma[2] - sigma[0])*(sigma[2] - sigma[0]) +
					6 * (sigma[3] * sigma[3] + sigma[4] * sigma[4] + sigma[5] * sigma[5])) / sqrt(2.);
				double eps_inv = sqrt((eps[0] - eps[1])*(eps[0] - eps[1]) + (eps[1] - eps[2])*(eps[1] - eps[2]) + (eps[2] - eps[0])*(eps[2] - eps[0]) +
					3 * (eps[3] * eps[3] + eps[4] * eps[4] + eps[5] * eps[5]) / 2.) * sqrt(2.) / 3;

				value[0][i] = sigma[0];
				value[1][i] = sigma[1];
				value[2][i] = sigma[2];

				value[3][i] = eps[0];
				value[4][i] = eps[1];
				value[5][i] = eps[2];

				value[3][i] = U.x;
				value[4][i] = U.y;
				value[5][i] = element->GetIdDomain();// U.z;

				/*value[3][i] = sigma_inv;
				value[4][i] = eps_inv;
				value[5][i] = 0;*/

				if (sigma_inv > sigma_inv_max)
				{
					sigma_inv_max = sigma_inv;
					elem_sigma_max = i;
				}

				/*if (base_grid.elem[i].get_id_domain() == 1)
				{
				value[0][i] = 0;
				value[1][i] = 0;
				value[2][i] = 0;

				value[3][i] = 0;
				value[4][i] = 0;
				value[5][i] = 0;
				}*/


				/*value[0][i] = eps[0];
				value[1][i] = eps[1];
				value[2][i] = eps[2];*/
			}
			for (int id_DOF = 0; id_DOF < Solution.size(); id_DOF++)
			{
				solver_grid.MoveTheVertex(solver_grid.accordance_DOF_and_vertex[id_DOF], Solution[id_DOF]);
			}
			solver_grid.printTecPlot3D(fout_tech, value, name_value, name_in_file);
			for (int id_DOF = 0; id_DOF < Solution.size(); id_DOF++)
			{
				solver_grid.MoveTheVertex(solver_grid.accordance_DOF_and_vertex[id_DOF], Solution[id_DOF] * (-1));
			}
			fclose(fout_tech);
		}
		if(false){
			auto OutputIntoLine = [&](int line_size, Point<double> base_x, Point<bool> variated,
				std::vector<Point<double>> &Line, std::vector<std::vector<double>> &value)
			{
				Point<double> max, min;
				max = solver_grid.GetMaxCoordinate();
				min = solver_grid.GetMinCoordinate();
				/*Line.resize(line_size + 1);
				value.resize(6);
				for (int i = 0; i < value.size(); i++)
					value[i].resize(Line.size());*/

				if (variated.x)
				{
					double step = (max.x - min.x) / line_size;
					for (int i = 0; i < line_size + 1; i++)
					{
						Line[i] = Point<double>(min.x + step * i, base_x.y, base_x.z);
					}
				}
				if (variated.y)
				{
					double step = (max.y - min.y) / line_size;
					for (int i = 0; i < line_size + 1; i++)
					{
						Line[i] = Point<double>(base_x.x, min.y + step * i, base_x.z);
					}
				}
				if (variated.z)
				{
					double step = (max.z - min.z) / line_size;
					for (int i = 0; i < line_size + 1; i++)
					{
						Line[i] = Point<double>(base_x.x, base_x.y, min.z + step * i);
					}
				}

				for (int i = 0; i < value[0].size(); i++)
				{
					int base_element = solver_grid.GetElementID(Line[i]);
					auto element = solver_grid.GetElement(base_element);

					Point<double> Centr = Line[i];//= element->GetWeightCentr();

					Point<double> U = solver_grid.GetSolutionInPoint(base_element, Centr, Solution);
					Point<Point<double>> dU = solver_grid.GetDerevativeFromSolutionInPoint(base_element, Centr, Solution);

					double eps[6] = { dU.x.x, dU.y.y, dU.z.z, dU.x.y + dU.y.x, dU.y.z + dU.z.y, dU.x.z + dU.z.x };
					double v = solver_grid.GetDomain(0)->forMech.v;
					double E = solver_grid.GetDomain(0)->forMech.GetE(0);
					double a = v / (1 - v);
					double b = (1 - 2 * v) / (2 * (1 - v));
					//double k = E * (1 - v) / (1 - v * v);
					double k = E * (1. - v) / (1. + v) / (1. - 2. * v);
					double sigma[6] = { k*(eps[0] * 1 + eps[1] * a + eps[2] * a), k*(eps[0] * a + eps[1] * 1 + eps[2] * a), k*(eps[0] * a + eps[1] * a + eps[2] * 1),
						k*(2 * eps[3] * b), k*(2 * eps[4] * b), k*(2 * eps[5] * b) };
					double sigma_inv = sqrt((sigma[0] - sigma[1])*(sigma[0] - sigma[1]) + (sigma[1] - sigma[2])*(sigma[1] - sigma[2]) + (sigma[2] - sigma[0])*(sigma[2] - sigma[0]) +
						6 * (sigma[3] * sigma[3] + sigma[4] * sigma[4] + sigma[5] * sigma[5])) / sqrt(2.);
					double eps_inv = sqrt((eps[0] - eps[1])*(eps[0] - eps[1]) + (eps[1] - eps[2])*(eps[1] - eps[2]) + (eps[2] - eps[0])*(eps[2] - eps[0]) +
						3 * (eps[3] * eps[3] + eps[4] * eps[4] + eps[5] * eps[5]) / 2.) * sqrt(2.) / 3;

					value[0][i] = sigma[0];
					value[1][i] = sigma[1];
					value[2][i] = sigma[2];

					value[3][i] = eps[0];
					value[4][i] = eps[1];
					value[5][i] = eps[2];

					value[3][i] = U.x;
					value[4][i] = U.y;
					value[5][i] = U.z;
				}

				return;
			};

			std::vector<double> variant_x(6);
			variant_x[0] = 50;
			variant_x[1] = 71;
			variant_x[2] = 72;
			variant_x[3] = 75;
			variant_x[4] = 80;
			variant_x[5] = 95;
			for (int i = 0; i < variant_x.size(); i++)
			{
				std::vector<Point<double>> line;
				std::vector<std::vector<double>> value;
				Point<double> base_point(variant_x[i], 25, 0);
				int line_size = 400;

				line.resize(line_size + 1);
				value.resize(6);
				for (int i = 0; i < value.size(); i++)
					value[i].resize(line.size());


				OutputIntoLine(
					400, 
					base_point, 
					Point<bool>(false, false, true),
					line, value);

				FILE *fout;
				char name_f[5000];
				sprintf_s(name_f, "%s/Result in line (%.2lf, %.2lf, %.2lf).txt", base_result_directory,
					base_point.x, base_point.y, base_point.z);
				printf_s("%s/Result in line (%.2lf, %.2lf, %.2lf).txt\n", base_result_directory,
					base_point.x, base_point.y, base_point.z);
				fopen_s(&fout, name_f, "w");
				fprintf_s(fout, "X Y Z sigma_xx sigma_yy sigma_zz Ux Uy Uz\n");
				for (int j = 0; j < line.size(); j++)
				{
					fprintf_s(fout, "%.4e %.4e %.4e ", line[j].x, line[j].y, line[j].z);
					for (int k = 0; k < value.size(); k++)
					{
						fprintf_s(fout, "%.12e ", value[k][j]);
					}
					fprintf_s(fout, "\n");
				}
				fclose(fout);
			}
		}
		printf_s("\t complite\n");

		sprintf_s(cracks_directory, sizeof(cracks_directory), "%s", base_result_directory);
		solver_grid.~Grid();
	}
}

void Solve_QuasiStatiñElasticProblem_withCracks()
{
	MultiXFEM::Grid solver_grid; //output
	std::vector<Point<double>> Solution; //output


	clock_t t_after = clock();
	double start = omp_get_wtime();


	//char properties_file[1000] = { "E:/+cyl/800el/param.txt" };
	//char properties_file[1000] = { "D:/Documents/universiti/MySOFT/TEST/box_200x200x200_crañks/800el_disc/base_properties.txt" };
	//char properties_file[1000] = { "E:/27_hydrofracture/area_200x200x200/340kel/Crack_30x200_Z=100_Y0grad/param.txt" };
	char properties_file[1000] = { "D:/Babenko_2022/1cracks/XFEM/50/param.txt" };
	//char properties_file[1000] = { "base_properties.txt" };
	printf_s("The properties file contains the information:\n==========================================\n");
	printf_s("\tDo you need to print to a log file? 0-false/1-true\n");
	printf_s("\tinput: Mesh directory: box_200x200x200/90el\n");
	printf_s("\tinput: Crack directory: box_200x200x200/90el/Crack_100x200_Z=50\n");
	printf_s("\tproperties: Iterative process:\n");
	printf_s("\t\t<Start iteration> <End iteration>\n");
	printf_s("\t\t<Step size for crack prolongation>\n");
	printf_s("\t\t<Critical residual value for SLAE solution>\n");
	printf_s("\t<properties: Boundary values>\n");
	printf_s("\t<properties: Materials>\n");
	printf_s("\t<output: Result directory: box_200x200x200/90el/result>\n==========================================\n");
	

	printf_s("Enter the name of the properties file: ");
	//scanf_s("%s", &properties_file);

	FILE *f_properties;
	fopen_s(&f_properties, properties_file, "r");
	if (f_properties == NULL)
	{
		printf_s("\nError in properties file\n");
	}
	bool is_print_logFile = false;
	char mesh_directory[1000];
	char cracks_directory[1000];
	char base_result_directory[1000];


	int _flag;
	fscanf_s(f_properties, "%d", &_flag);
	if (_flag == 1) is_print_logFile = true;

	math::ReadNonEmptyLine(f_properties, mesh_directory);
	math::SimpleGrid geo_grid; //input
	geo_grid.ReadFromNVTR(mesh_directory, 4);

	math::ReadNonEmptyLine(f_properties, cracks_directory);

	int star_iteration, end_iteration;
	{
		char _line[1000];
		math::ReadNonEmptyLine_forNumbers(f_properties, _line);
		std::vector<int> val;
		math::ParserStringToVectorInt(_line, val, " ");
		star_iteration = val[0];
		end_iteration = val[1];
	}
	double step_size;
	{
		char _line[1000];
		math::ReadNonEmptyLine_forNumbers(f_properties, _line);
		std::vector<float> val;
		math::ParserStringToVectorFloat(_line, val, " ");
		step_size = val[0];
	}
	int crack_smoothing_coefficient;
	{
		char _line[1000];
		math::ReadNonEmptyLine_forNumbers(f_properties, _line);
		std::vector<int> val;
		math::ParserStringToVectorInt(_line, val, " ");
		crack_smoothing_coefficient = val[0];
	}
	double critical_residual;
	{
		char _line[1000];
		math::ReadNonEmptyLine_forNumbers(f_properties, _line);
		std::vector<float> val;
		math::ParserStringToVectorFloat(_line, val, " ");
		critical_residual = val[0];
	}
	double pressure_of_fluid;
	{
		char _line[1000];
		math::ReadNonEmptyLine_forNumbers(f_properties, _line);
		std::vector<float> val;
		math::ParserStringToVectorFloat(_line, val, " ");
		pressure_of_fluid = val[0];
	}

	//read Boundary values and vertexes
	struct _Dirichlet {
		//Point<double> value;
		//Point<bool> is_condition;
		std::vector<int> id_vertexes;
		std::function<Point<double>(Point<bool>&)> value;
	};
	std::vector<_Dirichlet> first_boundaries;
	{
		int Nb;
		{
			char _line[1000];
			math::ReadNonEmptyLine_forNumbers(f_properties, _line);
			std::vector<int> val;
			math::ParserStringToVectorInt(_line, val, " ");
			Nb = val[0];
		}
		first_boundaries.resize(Nb);
		char line[1000];
		for (int i = 0; i < Nb; i++)
		{
			Point<double> value;
			Point<bool> is_condition;

			math::ReadNonEmptyLine_forNumbers(f_properties, line);
			int jj = 0;
			int curr_i = 0;
			char _tmp_line[1000];
			for (int j = 0; j < 1000; j++)
			{
				switch (line[j])
				{
				case '\0': j = 1000; break;
				case '\n': j = 1000; break;
				case '#': j = 1000; break;
				case '*':
				{
					switch (jj)
					{
					case 0: is_condition.x = false; break;
					case 1: is_condition.y = false; break;
					case 2: is_condition.z = false; break;
					default:
						break;
					}
					jj++;
					break;
				}
				case ' ':
				{
					if (curr_i > 0)
					{
						std::vector<float> _val;
						_tmp_line[curr_i] = '\0';
						math::ParserStringToVectorFloat(_tmp_line, _val, " *");
						curr_i = 0;

						switch (jj)
						{
						case 0: is_condition.x = true; value.x = _val[0]; break;
						case 1: is_condition.y = true; value.y = _val[0]; break;
						case 2: is_condition.z = true; value.z = _val[0]; break;
						default:
							break;
						}
						jj++;
					}
					break;
				}
				case '\t': break;
				default:
				{
					_tmp_line[curr_i] = line[j];
					curr_i++;
				}
				break;
				}
			}

			first_boundaries[i].value = [value, is_condition](Point<bool> &is_take) -> Point<double>
			{
				is_take = is_condition;
				return value;
			};
		}

		if (Nb != 0)
		{
			char name_boundary_file[1000];
			sprintf_s(name_boundary_file, sizeof(name_boundary_file), "%s/kraev1.txt", mesh_directory);
			FILE *fbv;
			fopen_s(&fbv, name_boundary_file, "r");
			int N;
			fscanf_s(fbv, "%d", &N);
			for (int i = 0; i < N; i++)
			{
				int _type, _vertex;
				fscanf_s(fbv, "%d %d", &_vertex, &_type);
				first_boundaries[_type].id_vertexes.push_back(_vertex);
			}
			if (geo_grid.read_from_zero == false)
			{
				for (int t = 0; t < first_boundaries.size(); t++)
				{
					for (int i = 0; i < first_boundaries[t].id_vertexes.size(); i++)
					{
						first_boundaries[t].id_vertexes[i]--;
					}
				}
			}
			fclose(fbv);
		}
	}

	struct _Neumann {
		/*double value;
		Point<double> vector;*/
		std::function<Point<double>(Point<double>)> value;
		std::vector<std::function<Point<double>(Point<double>)>> values;
		std::vector<std::vector<int>> id_vertexes_as_triangle;
		std::vector<int> id_base_element;
	};
	std::vector<_Neumann> second_boundary;
	{
		int Nb;
		{
			char _line[1000];
			math::ReadNonEmptyLine_forNumbers(f_properties, _line);
			std::vector<int> val;
			math::ParserStringToVectorInt(_line, val, " ");
			Nb = val[0];
		}
		second_boundary.resize(Nb);
		char line[1000];
		std::vector<bool> is_individual_values(Nb);
		std::vector<double> individual_values(Nb);
		for (int i = 0; i < Nb; i++)
		{
			is_individual_values[i] = false;
			double value;
			Point<double> vector;

			math::ReadNonEmptyLine_forNumbers(f_properties, line);
			int jj = 0;
			int curr_i = 0;
			char _tmp_line[1000];
			for (int j = 0; j < 1000; j++)
			{
				if ((curr_i > 0) &&(line[j] == '\n' || line[j] == '\0' || line[j] == '#'))
				{
					std::vector<float> val;
					_tmp_line[curr_i] = '\0';
					math::ParserStringToVectorFloat(_tmp_line, val, " ");
					curr_i = 0;

					value = val[0];
					vector.x = val[1];
					vector.y = val[2];
					vector.z = val[3];

					individual_values[i] = value;
					if (math::IsEqual(math::SolveLengthVector(vector), 0.0))
					{
						is_individual_values[i] = true;
					}

					break;
				}
				else
				{
					if (line[j] != '\t')
					{
						_tmp_line[curr_i] = line[j];
						curr_i++;
					}
				}
			}

			if (is_individual_values[i] == false)
			{
				second_boundary[i].value = [value, vector](Point<double> X) -> Point<double>
				{
					Point <double> res;
					res.x = vector.x *value;
					res.y = vector.y *value;
					res.z = vector.z *value;
					return res;
				};
			}
		}

		if (Nb != 0)
		{
			std::vector<std::vector<int>> id_vertexes;
			id_vertexes.resize(Nb);
			char name_boundary_file[1000];
			sprintf_s(name_boundary_file, sizeof(name_boundary_file), "%s/kraev2.txt", mesh_directory);
			FILE *fbv;
			fopen_s(&fbv, name_boundary_file, "r");
			int N;
			fscanf_s(fbv, "%d", &N);
			for (int i = 0; i < N; i++)
			{
				int _type, _vertex;
				fscanf_s(fbv, "%d %d", &_vertex, &_type);
				id_vertexes[_type].push_back(_vertex);
			}
			if (geo_grid.read_from_zero == false)
			{
				for (int t = 0; t < id_vertexes.size(); t++)
				{
					for (int i = 0; i < id_vertexes[t].size(); i++)
					{
						id_vertexes[t][i]--;
					}
				}
			}
			fclose(fbv);

			for (int id_type = 0; id_type < Nb; id_type++)
			{
				for (int id_element = 0; id_element < geo_grid.nvtr.size(); id_element++)
				{
					auto _tmp = math::GetConfluence(geo_grid.nvtr[id_element], id_vertexes[id_type]);
					if (_tmp.size() == 3)
					{
						//second_boundary[id_type].id_vertexes_as_triangle.push_back(_tmp);
						//second_boundary[id_type].id_base_element.push_back(id_element);

						int test_vertex=-1;
						for (int t = 0; t < geo_grid.nvtr[id_element].size(); t++)
						{
							bool find_id = false;
							for (int tt = 0; tt < _tmp.size(); tt++)
							{
								if (geo_grid.nvtr[id_element][t] == _tmp[tt])
								{
									find_id = true;
									break;
								}
							}
							if (!find_id)
							{
								test_vertex = geo_grid.nvtr[id_element][t];
								break;
							}
						}

						double A, B, C, D;
						Point<double> test_vector = geo_grid.xyz[test_vertex];
						std::vector<Point<double>> vertexes(3);
						vertexes[0] = geo_grid.xyz[_tmp[0]];
						vertexes[1] = geo_grid.xyz[_tmp[1]];
						vertexes[2] = geo_grid.xyz[_tmp[2]];
						bool reverse = false;
						math::GetPlaneEquation(vertexes, A, B, C, D);

						if (Point<double>(A, B, C)*test_vector > 0)
						{
							int t = _tmp[1];
							_tmp[1] = _tmp[2];
							_tmp[2] = t;

							vertexes[0] = geo_grid.xyz[_tmp[0]];
							vertexes[1] = geo_grid.xyz[_tmp[1]];
							vertexes[2] = geo_grid.xyz[_tmp[2]];
						}
						Point<double> normal;
						double d;
						math::GetPlaneEquation(vertexes, normal.x, normal.y, normal.z, d);
						normal /= math::SolveLengthVector(normal);

						second_boundary[id_type].id_vertexes_as_triangle.push_back(_tmp);
						second_boundary[id_type].id_base_element.push_back(id_element);

						if (is_individual_values[id_type] == true)
						{
							double _val = individual_values[id_type];
							std::function<Point<double>(Point<double>)> curr_val = [_val, normal](Point<double> X) -> Point<double>
							{
								Point <double> res;
								res.x = normal.x * _val;
								res.y = normal.y * _val;
								res.z = normal.z * _val;
								return res;
							};
							second_boundary[id_type].values.push_back(curr_val);
						}
					}
				}
			}
		}
	}

	//read Materials
	struct _material
	{
		double _E;
		double _v;
		double _K_Ic;
		double _pressure_of_fluid;

		_material(double _E, double _v, double _K_Ic, double _pressure_of_fluid)
		{
			this->_E = _E;
			this->_v = _v;
			this->_K_Ic = _K_Ic;
			this->_pressure_of_fluid = _pressure_of_fluid;

		}
	};
	std::vector<_material> _materials;
	int N_domain;
	{
		char _line[1000];
		math::ReadNonEmptyLine_forNumbers(f_properties, _line);
		std::vector<int> val;
		math::ParserStringToVectorInt(_line, val, " ");
		N_domain = val[0];
	}
	for (int id_domain = 0; id_domain < N_domain; id_domain++)
	{
		double _E, _v, _K_Ic;
		{
			char _line[1000];
			math::ReadNonEmptyLine_forNumbers(f_properties, _line);
			std::vector<float> val;
			math::ParserStringToVectorFloat(_line, val, " ");
			_E = val[0];
			_v = val[1];
			_K_Ic = val[2];
		}
		_materials.push_back(_material(_E, _v, _K_Ic, pressure_of_fluid));
	}

	math::ReadNonEmptyLine(f_properties, base_result_directory);
	fclose(f_properties);

	CreateDirectory((LPCTSTR)base_result_directory, NULL);

	if (star_iteration != 0)
	{
		sprintf_s(cracks_directory, sizeof(cracks_directory), "%s/STEP_%d", base_result_directory, star_iteration);
	}

	bool res = false;
	for (int id_STEP = star_iteration; id_STEP < end_iteration; id_STEP++)
	{
		printf_s("\n=============================*===============================\n");
		printf_s("================= Start solution of %d STEP ================\n", id_STEP);
		//Create new directories
		char result_directory[1000];
		sprintf_s(result_directory, sizeof(result_directory), "%s/STEP_%d", base_result_directory, id_STEP);
		CreateDirectory((LPCTSTR)result_directory, NULL);

		for (int i = 0; i < _materials.size(); i++)
		{
			solver_grid.AddDomain();
			auto domain = solver_grid.GetDomain(i);
			domain->forMech.SetE(_materials[i]._E);
			domain->forMech.SetV(_materials[i]._v);
			domain->forMech.SetK_Ic(_materials[i]._K_Ic);
			domain->forMech.SetPressureOfFluid(_materials[i]._pressure_of_fluid);
		}

		MultiXFEM::MultiXFEM_forElasticDeformation(
			is_print_logFile,
			critical_residual,
			cracks_directory, //input
			step_size,
			crack_smoothing_coefficient,
			geo_grid, //input
			first_boundaries, //input
			second_boundary, //input
			result_directory, //output
			solver_grid, //output
			Solution //output
		);

		double end = omp_get_wtime();
		printf("start = %.16g\nend = %.16g\ndiff = %.16g\n", start, end, end - start);

		FILE* ftime;
		char ftime_name[1000];
		sprintf_s(ftime_name, "%s/time_result.dat", base_result_directory);
		fopen_s(&ftime, ftime_name, "w");
		fprintf_s(ftime, "start = %.16g\nend = %.16g\ndiff = %.16g\n", start, end, end - start);
		fclose(ftime);

		//output solution
		printf_s("Print the mech result into .dat file... ");
		{ //via full 3D mesh
			FILE *fout_tech;
			char name_u_tech[5000];
			sprintf_s(name_u_tech, "%s/U.dat", result_directory);
			fopen_s(&fout_tech, name_u_tech, "w");
			char name_in_file[1000];
			sprintf_s(name_in_file, "Elastic");
			std::vector<std::vector<char>> name_value(6);
			char name_v_tmp[6][100];
			sprintf_s(name_v_tmp[0], "sigma_xx");
			sprintf_s(name_v_tmp[1], "sigma_yy");
			sprintf_s(name_v_tmp[2], "sigma_zz");
			sprintf_s(name_v_tmp[3], "VonMises_stress");
			sprintf_s(name_v_tmp[4], "eps_inv");
			sprintf_s(name_v_tmp[5], "material");
			for (int i = 0; i < name_value.size(); i++)
			{
				name_value[i].resize(100);
				for (int j = 0; j < name_value[i].size(); j++)
				{
					name_value[i][j] = name_v_tmp[i][j];
				}
			}
			std::vector<std::vector<double>> value(3 * 2);
			value[0].resize(solver_grid.GetElementsCount());
			value[1].resize(solver_grid.GetElementsCount());
			value[2].resize(solver_grid.GetElementsCount());
			value[3].resize(solver_grid.GetElementsCount());
			value[4].resize(solver_grid.GetElementsCount());
			value[5].resize(solver_grid.GetElementsCount());
			double sigma_inv_max = 0;
			int elem_sigma_max = 0;
			for (int i = 0; i < solver_grid.GetElementsCount(); i++)
			{
				auto element = solver_grid.GetElement(i);

				Point<double> Centr = element->GetWeightCentr();

				Point<double> U = solver_grid.GetSolutionInPoint(i, Centr, Solution);
				Point<Point<double>> dU = solver_grid.GetDerevativeFromSolutionInPoint(i, Centr, Solution);

				double eps[6] = { dU.x.x, dU.y.y, dU.z.z, dU.x.y + dU.y.x, dU.y.z + dU.z.y, dU.x.z + dU.z.x };
				double v = solver_grid.GetDomain(0)->forMech.v;
				double E = solver_grid.GetDomain(0)->forMech.GetE(0);
				double a = v / (1 - v);
				double b = (1 - 2 * v) / (2 * (1 - v));
				double k = E * (1 - v) / ((1 + v)*(1 - v));
				double sigma[6] = { k*(eps[0] * 1 + eps[1] * a + eps[2] * a), k*(eps[0] * a + eps[1] * 1 + eps[2] * a), k*(eps[0] * a + eps[1] * a + eps[2] * 1),
					k*(2 * eps[3] * b), k*(2 * eps[4] * b), k*(2 * eps[5] * b) };
				double sigma_inv = sqrt((sigma[0] - sigma[1])*(sigma[0] - sigma[1]) + (sigma[1] - sigma[2])*(sigma[1] - sigma[2]) + (sigma[2] - sigma[0])*(sigma[2] - sigma[0]) +
					6 * (sigma[3] * sigma[3] + sigma[4] * sigma[4] + sigma[5] * sigma[5])) / sqrt(2.);
				double eps_inv = sqrt((eps[0] - eps[1])*(eps[0] - eps[1]) + (eps[1] - eps[2])*(eps[1] - eps[2]) + (eps[2] - eps[0])*(eps[2] - eps[0]) +
					3 * (eps[3] * eps[3] + eps[4] * eps[4] + eps[5] * eps[5]) / 2.) * sqrt(2.) / 3;

				double mises_stress = sqrt((sigma[0] - sigma[1]) * (sigma[0] - sigma[1])
					+ (sigma[1] - sigma[2]) * (sigma[1] - sigma[2])
					+ (sigma[0] - sigma[2]) * (sigma[0] - sigma[2])
					+ 6 * (sigma[3] * sigma[3] + sigma[4] * sigma[4] + sigma[5] * sigma[5]) / 2.0);

				value[0][i] = sigma[0];
				value[1][i] = sigma[1];
				value[2][i] = sigma[2];

				value[3][i] = eps[0];
				value[4][i] = eps[1];
				value[5][i] = eps[2];

				value[3][i] = mises_stress;
				value[4][i] = U.y;
				value[5][i] = element->GetDOFsCount() <= 4 ? 0 : 1;

				/*value[0][i] = 0;
				value[1][i] = 0;
				value[2][i] = 0;

				value[3][i] = 0;
				value[4][i] = 0;*/
				value[5][i] = element->GetDOFsCount();
				
			}
			/*for (int id_DOF = 0; id_DOF < Solution.size(); id_DOF++)
			{
				solver_grid.MoveTheVertex(solver_grid.accordance_DOF_and_vertex[id_DOF], Solution[id_DOF]);
			}*/
			solver_grid.printTecPlot3D(fout_tech, value, name_value, name_in_file);
			/*for (int id_DOF = 0; id_DOF < Solution.size(); id_DOF++)
			{
				solver_grid.MoveTheVertex(solver_grid.accordance_DOF_and_vertex[id_DOF], Solution[id_DOF]*(-1));
			}*/
			fclose(fout_tech);
		}
		if (true) //via face 2D mesh
		{
			math::SimpleGrid out_grid;
			char name_grid[5000];
			sprintf_s(name_grid, "%s/Face_OYZ_forOut.dat", mesh_directory);
			out_grid.ReadFromSalomeDat(name_grid, 2);

			FILE* fout_tech;
			FILE* fout_res;
			char name_u_tech[5000];
			char name_res[5000];
			sprintf_s(name_u_tech, "%s/U_in_face.dat", result_directory);
			sprintf_s(name_res, "%s/u_z_in_face.txt", result_directory);
			fopen_s(&fout_tech, name_u_tech, "w");
			fopen_s(&fout_res, name_res, "w");
			
			char name_in_file[1000];
			sprintf_s(name_in_file, "Elastic");
			std::vector<std::vector<char>> name_value(6);
			char name_v_tmp[6][100];
			sprintf_s(name_v_tmp[0], "sigma_xx");
			sprintf_s(name_v_tmp[1], "sigma_yy");
			sprintf_s(name_v_tmp[2], "sigma_zz");
			sprintf_s(name_v_tmp[3], "VonMises_stress");
			sprintf_s(name_v_tmp[4], "eps_inv");
			sprintf_s(name_v_tmp[5], "material");
			for (int i = 0; i < name_value.size(); i++)
			{
				name_value[i].resize(100);
				for (int j = 0; j < name_value[i].size(); j++)
				{
					name_value[i][j] = name_v_tmp[i][j];
				}
			}
			std::vector<std::vector<double>> value(3 * 2);
			value[0].resize(out_grid.nvtr.size());
			value[1].resize(out_grid.nvtr.size());
			value[2].resize(out_grid.nvtr.size());
			value[3].resize(out_grid.nvtr.size());
			value[4].resize(out_grid.nvtr.size());
			value[5].resize(out_grid.nvtr.size());
			double sigma_inv_max = 0;
			int elem_sigma_max = 0;
			for (int i = 0; i < out_grid.nvtr.size(); i++)
			{

				Point<double> Centr;
				for (int ii = 0; ii < out_grid.nvtr[i].size(); ii++)
				{
					Centr += out_grid.xyz[out_grid.nvtr[i][ii]] / 3.0;
				}

				double len;
				int element_id = abs(solver_grid.GetNearestElementID(Centr, len));
				auto element = solver_grid.GetElement(element_id);
				

				Point<double> U = solver_grid.GetSolutionInPoint(element_id, Centr, Solution);
				Point<Point<double>> dU = solver_grid.GetDerevativeFromSolutionInPoint(element_id, Centr, Solution);

				double eps[6] = { dU.x.x, dU.y.y, dU.z.z, dU.x.y + dU.y.x, dU.y.z + dU.z.y, dU.x.z + dU.z.x };
				double v = solver_grid.GetDomain(0)->forMech.v;
				double E = solver_grid.GetDomain(0)->forMech.GetE(0);
				double a = v / (1 - v);
				double b = (1 - 2 * v) / (2 * (1 - v));
				double k = E * (1 - v) / ((1 + v) * (1 - v));
				double sigma[6] = { k * (eps[0] * 1 + eps[1] * a + eps[2] * a), k * (eps[0] * a + eps[1] * 1 + eps[2] * a), k * (eps[0] * a + eps[1] * a + eps[2] * 1),
					k * (2 * eps[3] * b), k * (2 * eps[4] * b), k * (2 * eps[5] * b) };
				double sigma_inv = sqrt((sigma[0] - sigma[1]) * (sigma[0] - sigma[1]) + (sigma[1] - sigma[2]) * (sigma[1] - sigma[2]) + (sigma[2] - sigma[0]) * (sigma[2] - sigma[0]) +
					6 * (sigma[3] * sigma[3] + sigma[4] * sigma[4] + sigma[5] * sigma[5])) / sqrt(2.);
				double eps_inv = sqrt((eps[0] - eps[1]) * (eps[0] - eps[1]) + (eps[1] - eps[2]) * (eps[1] - eps[2]) + (eps[2] - eps[0]) * (eps[2] - eps[0]) +
					3 * (eps[3] * eps[3] + eps[4] * eps[4] + eps[5] * eps[5]) / 2.) * sqrt(2.) / 3;

				double mises_stress = sqrt((sigma[0] - sigma[1]) * (sigma[0] - sigma[1])
					+ (sigma[1] - sigma[2]) * (sigma[1] - sigma[2])
					+ (sigma[0] - sigma[2]) * (sigma[0] - sigma[2])
					+ 6 * (sigma[3] * sigma[3] + sigma[4] * sigma[4] + sigma[5] * sigma[5]) / 2.0);

				value[0][i] = sigma[0];
				value[1][i] = sigma[1];
				value[2][i] = sigma[2];

				value[3][i] = eps[0];
				value[4][i] = eps[1];
				value[5][i] = eps[2];

				value[3][i] = mises_stress;
				value[4][i] = U.y;
				value[5][i] = element->GetDOFsCount() <= 4 ? 0 : 1;

				fprintf_s(fout_res, "%.15e\n", U.z);

				/*value[0][i] = 0;
				value[1][i] = 0;
				value[2][i] = 0;

				value[3][i] = 0;
				value[4][i] = 0;*/
				value[5][i] = element->GetDOFsCount();

			}

			out_grid.printTecPlot3D(fout_tech, value, name_value, name_in_file);
			fclose(fout_res);
		}
		printf_s("\t complite\n");

		//MultiXFEM::CrackPropagation_3D(solver_grid, Solution, result_directory, id_STEP);
		//MultiXFEM::CrackPropagation_3D_Cherepanov(solver_grid, Solution, result_directory, id_STEP);
		//MultiXFEM::CrackPropagation_3D_Cherepanov_v2(solver_grid, Solution, result_directory, id_STEP);

		sprintf_s(cracks_directory, sizeof(cracks_directory), "%s", result_directory);
		solver_grid.~Grid();
	}
}

int main()
{
	//Solve_StatiñElasticProblem_withCracks();
	
	Solve_QuasiStatiñElasticProblem_withCracks();

	int a;
	scanf_s("%d", &a);
}