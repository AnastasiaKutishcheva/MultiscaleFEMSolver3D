#pragma once
#include "MsFEM.h"
#include <Windows.h>
#include <iostream>

void CreateCubGridAsPoly_self_format()
{
	Point<double> size(200,200,200);
	Point<int> N_elem(10, 10, 10);
	Point<double> h;
	h.x = size.x / N_elem.x;
	h.y = size.y / N_elem.y;
	h.z = size.z / N_elem.z;
	FILE* fout;

	fopen_s(&fout, "param.txt", "w");
	fprintf_s(fout, "%d %d", (N_elem.x + 1) * (N_elem.y + 1) * (N_elem.z + 1), N_elem.x * N_elem.y * N_elem.z);
	fclose(fout);

	fopen_s(&fout, "xyz.txt", "w");
	for (int i = 0; i < N_elem.x + 1; i++)
	{
		for (int j = 0; j < N_elem.y + 1; j++)
		{
			for (int k = 0; k < N_elem.z + 1; k++)
			{
				fprintf_s(fout, "%lf %lf %lf\n", h.x * i, h.y * j, h.z * k);
			}
		}
	}
	fclose(fout);

	fopen_s(&fout, "kraev1.txt", "w");
	fprintf_s(fout, "%d\n", (N_elem.x + 1) * (N_elem.y + 1) * 2);
	for (int i = 0; i < N_elem.x + 1; i++)
	{
		for (int j = 0; j < N_elem.y + 1; j++)
		{
			fprintf_s(fout, "%d %d\n", i + j * (N_elem.x+1), 0);
		}
	}
	for (int i = 0; i < N_elem.x + 1; i++)
	{
		for (int j = 0; j < N_elem.y + 1; j++)
		{
			fprintf_s(fout, "%d %d\n", i + j * (N_elem.x + 1) + N_elem.z * (N_elem.x + 1) * (N_elem.y + 1), 1);
		}
	}
	fclose(fout);

	fopen_s(&fout, "nvtr.txt", "w");
	for (int i = 0; i < N_elem.x; i++)
	{
		for (int j = 0; j < N_elem.y; j++)
		{
			for (int k = 0; k < N_elem.z; k++)
			{
				fprintf_s(fout, "%d\n", 8);
				fprintf_s(fout, "%d ", i + j * (N_elem.x + 1) + k * (N_elem.x + 1) * (N_elem.y + 1));
				fprintf_s(fout, "%d ", i + 1 + j * (N_elem.x + 1) + k * (N_elem.x + 1) * (N_elem.y + 1));
				fprintf_s(fout, "%d ", i + (j + 1) * (N_elem.x + 1) + k * (N_elem.x + 1) * (N_elem.y + 1));
				fprintf_s(fout, "%d ", i + 1 + (j + 1) * (N_elem.x + 1) + k * (N_elem.x + 1) * (N_elem.y + 1));
				fprintf_s(fout, "%d ", i + j * (N_elem.x + 1) + (k+1) * (N_elem.x + 1) * (N_elem.y + 1));
				fprintf_s(fout, "%d ", i + 1 + j * (N_elem.x + 1) + (k+1) * (N_elem.x + 1) * (N_elem.y + 1));
				fprintf_s(fout, "%d ", i + (j + 1) * (N_elem.x + 1) + (k+1) * (N_elem.x + 1) * (N_elem.y + 1));
				fprintf_s(fout, "%d\n", i + 1 + (j + 1) * (N_elem.x + 1) + (k+1) * (N_elem.x + 1) * (N_elem.y + 1));
				fprintf_s(fout, "%d\n", 6);
				fprintf_s(fout, "%d %d %d %d %d\n", 4, 0, 1, 3, 2);
				fprintf_s(fout, "%d %d %d %d %d\n", 4, 4, 5, 7, 6);
				fprintf_s(fout, "%d %d %d %d %d\n", 4, 0, 2, 6, 4);
				fprintf_s(fout, "%d %d %d %d %d\n", 4, 2, 3, 7, 6);
				fprintf_s(fout, "%d %d %d %d %d\n", 4, 1, 3, 7, 5);
				fprintf_s(fout, "%d %d %d %d %d\n", 4, 0, 1, 5, 4);
			}
		}
	}
	fclose(fout);
}
void CreateCubGridAsPoly_dat_format()
{
	Point<double> size(200, 200, 200);
	Point<int> N_elem(10, 10, 10);

	printf_s("Enter size(x,y,z): ");
	scanf_s("%lf %lf %lf", &size.x, &size.y, &size.z);
	printf_s("Enter number of elements(x,y,z): ");
	scanf_s("%d %d %d", &N_elem.x, &N_elem.y, &N_elem.z);

	Point<double> h;
	h.x = size.x / N_elem.x;
	h.y = size.y / N_elem.y;
	h.z = size.z / N_elem.z;
	
	
	
	FILE* fout;

	fopen_s(&fout, "Volume_poly.dat", "w");
	fprintf_s(fout, "%d %d\n", (N_elem.x + 1) * (N_elem.y + 1) * (N_elem.z + 1), N_elem.x * N_elem.y * N_elem.z);

	int num = 1;
	for (int i = 0; i < N_elem.x + 1; i++)
	{
		for (int j = 0; j < N_elem.y + 1; j++)
		{
			for (int k = 0; k < N_elem.z + 1; k++)
			{
				fprintf_s(fout, "%d %lf %lf %lf\n", num, h.x * i, h.y * j, h.z * k);
				num++;
			}
		}
	}

	for (int i = 0; i < N_elem.x; i++)
	{
		for (int j = 0; j < N_elem.y; j++)
		{
			for (int k = 0; k < N_elem.z; k++)
			{
				fprintf_s(fout, "%d ", 8);
				fprintf_s(fout, "%d ", i + j * (N_elem.x + 1) + k * (N_elem.x + 1) * (N_elem.y + 1));
				fprintf_s(fout, "%d ", i + 1 + j * (N_elem.x + 1) + k * (N_elem.x + 1) * (N_elem.y + 1));
				fprintf_s(fout, "%d ", i + (j + 1) * (N_elem.x + 1) + k * (N_elem.x + 1) * (N_elem.y + 1));
				fprintf_s(fout, "%d ", i + 1 + (j + 1) * (N_elem.x + 1) + k * (N_elem.x + 1) * (N_elem.y + 1));
				fprintf_s(fout, "%d ", i + j * (N_elem.x + 1) + (k + 1) * (N_elem.x + 1) * (N_elem.y + 1));
				fprintf_s(fout, "%d ", i + 1 + j * (N_elem.x + 1) + (k + 1) * (N_elem.x + 1) * (N_elem.y + 1));
				fprintf_s(fout, "%d ", i + (j + 1) * (N_elem.x + 1) + (k + 1) * (N_elem.x + 1) * (N_elem.y + 1));
				fprintf_s(fout, "%d ", i + 1 + (j + 1) * (N_elem.x + 1) + (k + 1) * (N_elem.x + 1) * (N_elem.y + 1));
				fprintf_s(fout, "%d ", 6);
				fprintf_s(fout, "%d %d %d %d %d ", 4, 0, 1, 3, 2);
				fprintf_s(fout, "%d %d %d %d %d ", 4, 4, 5, 7, 6);
				fprintf_s(fout, "%d %d %d %d %d ", 4, 0, 2, 6, 4);
				fprintf_s(fout, "%d %d %d %d %d ", 4, 2, 3, 7, 6);
				fprintf_s(fout, "%d %d %d %d %d ", 4, 1, 3, 7, 5);
				fprintf_s(fout, "%d %d %d %d %d\n", 4, 0, 1, 5, 4);
			}
		}
	}
	fclose(fout);
}

void Solve_ElasticDeformationProblem()
{
	MsFEM::Grid_forMech solver_grid; //output
	std::vector<Point<double>> Solution; //output
	CSSD_Matrix<Tensor2Rank3D, Point<double>> global_SLAE;

	//char properties_file[1000] = { "E:/+cyl/800el/param.txt" };
	//char properties_file[1000] = { "E:/Box/100x10x200/BoxWithCracks/67k/param.txt" };
	//char properties_file[1000] = { "E:/Box/100x10x200/BoxWithCracks/638k/param.txt" };
	//char properties_file[1000] = { "E:/Box/100x50x200/param.txt" };
	//char properties_file[1000] = { "D:/P-ref_MsFEM/2_el/param_for_solver.txt" };
	//char properties_file[1000] = { "D:/P-ref_MsFEM/Sample_for_convergense/PolyMesh_1Order/H2_h1/param_for_solver.txt" };
	//char properties_file[1000] = { "D:/P-ref_MsFEM/Sample2/macro_H1_micro_h1/param_for_solver.txt" };
	char properties_file[1000] = { "./param_for_solver.txt" };
	//char properties_file[1000] = { "base_properties.txt" };
	printf_s("The properties file contains the information:\n==========================================\n");
	printf_s("\tDo you need to print to a log file? 0-false/1-true\n");
	printf_s("\tinput: Mesh directory: box_200x200x200/90el\n");
	printf_s("\tproperties: Iterative process:\n");
	printf_s("\t\t<Critical residual value for SLAE solution>\n");
	printf_s("\t<properties: Boundary values>\n");
	printf_s("\t<properties: Materials>\n");
	printf_s("\t<output: Result directory: box_200x200x200/90el/result>\n==========================================\n");


	printf_s("Enter the name of the properties file: ");
	//scanf_s("%s", &properties_file);

	FILE* f_properties;
	fopen_s(&f_properties, properties_file, "r");
	if (f_properties == NULL)
	{
		printf_s("\nError in properties file\n");
	}
	bool is_print_logFile = false;
	char mesh_directory[1000];
	char base_result_directory[1000];


	int _flag;
	fscanf_s(f_properties, "%d", &_flag);
	if (_flag == 1) is_print_logFile = true;

	math::ReadNonEmptyLine(f_properties, mesh_directory);
	math::SimpleGrid geo_grid; //input
	geo_grid.ReadFromNVTR(mesh_directory, -1); //-1 is special flag for polyhedral grid 
	char fine_mesh_dir[1000];
	sprintf_s(fine_mesh_dir, sizeof(fine_mesh_dir), "%s/FinalGrids", mesh_directory);

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

			first_boundaries[i].value = [value, is_condition](Point<bool>& is_take) -> Point<double>
			{
				is_take = is_condition;
				return value;
			};
		}

		if (Nb != 0)
		{
			char name_boundary_file[1000];
			sprintf_s(name_boundary_file, sizeof(name_boundary_file), "%s/kraev1.txt", mesh_directory);
			FILE* fbv;
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
					res.x = vector.x * value;
					res.y = vector.y * value;
					res.z = vector.z * value;
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
			FILE* fbv;
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

						if (Point<double>(A, B, C) * test_vector > 0)
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
		for (int i = 0; i < _materials.size(); i++)
		{
			solver_grid.AddDomain();
			auto domain = solver_grid.GetDomain(i);
			domain->forMech.SetE(_materials[i]._E);
			domain->forMech.SetV(_materials[i]._v);
		}

		MsFEM::MsFEM_forElasticDeformation_OrderBF2_forPoly(
			is_print_logFile,
			critical_residual,
			fine_mesh_dir,
			geo_grid, //input
			first_boundaries, //input
			second_boundary, //input
			base_result_directory, //output
			solver_grid, //output
			Solution, //output
			global_SLAE
		);

		//output solution
		printf_s("Print the mech result into .dat file... ");
		for(int id_elem = 0; id_elem < solver_grid.GetElementsCount(); id_elem++)
		{
			auto macro_element = solver_grid.GetElement(id_elem);

			FILE* fout_tech;
			char name_u_tech[5000];
			sprintf_s(name_u_tech, "%s/U_elem%d.dat", base_result_directory, id_elem);
			fopen_s(&fout_tech, name_u_tech, "w");
			char name_in_file[1000];
			sprintf_s(name_in_file, "Elastic_elem%d", id_elem);
			std::vector<std::vector<char>> name_value(6);
			char name_v_tmp[6][100];
			sprintf_s(name_v_tmp[0], "sigma_xx");
			sprintf_s(name_v_tmp[1], "sigma_yy");
			sprintf_s(name_v_tmp[2], "sigma_zz");
			sprintf_s(name_v_tmp[3], "Ux");
			sprintf_s(name_v_tmp[4], "Uy");
			sprintf_s(name_v_tmp[5], "Uz");
			for (int i = 0; i < name_value.size(); i++)
			{
				name_value[i].resize(100);
				for (int j = 0; j < name_value[i].size(); j++)
				{
					name_value[i][j] = name_v_tmp[i][j];
				}
			}
			std::vector<std::vector<double>> value(3 * 2);
			value[0].resize( macro_element->self_grid.GetElementsCount());
			value[1].resize( macro_element->self_grid.GetElementsCount());
			value[2].resize( macro_element->self_grid.GetElementsCount());
			value[3].resize( macro_element->self_grid.GetElementsCount());
			value[4].resize( macro_element->self_grid.GetElementsCount());
			value[5].resize( macro_element->self_grid.GetElementsCount());
			double sigma_inv_max = 0;
			int elem_sigma_max = 0;
			for (int i = 0; i < macro_element->self_grid.GetElementsCount(); i++)
			{
				auto element = macro_element->self_grid.GetElement(i);

				Point<double> Centr_micro = element->GetWeightCentr();

				Point<double> U_macro;
				Point<Point<double>> dU_macro;

				std::vector<Point<double>> bf_macro(macro_element->GetDOFsCount());
				std::vector<Point<Point<double>>> d_bf_macro(macro_element->GetDOFsCount());
				for (int j = 0; j < bf_macro.size(); j++)
				{
					bf_macro[j] = macro_element->self_grid.GetSolutionInPoint(i, Centr_micro, macro_element->self_basis_functions[j]);
					d_bf_macro[j] = macro_element->self_grid.GetDerevativeFromSolutionInPoint(i, Centr_micro, macro_element->self_basis_functions[j]);

					Point<double> val;
					val.x = Solution[macro_element->GetDOFInLocalID(j)].x;
					val.y = Solution[macro_element->GetDOFInLocalID(j)].y;
					val.z = Solution[macro_element->GetDOFInLocalID(j)].z;

					U_macro.x += bf_macro[j].x * Solution[macro_element->GetDOFInLocalID(j)].x;
					U_macro.y += bf_macro[j].y * Solution[macro_element->GetDOFInLocalID(j)].y;
					U_macro.z += bf_macro[j].z * Solution[macro_element->GetDOFInLocalID(j)].z;

					dU_macro.x += d_bf_macro[j].x * Solution[macro_element->GetDOFInLocalID(j)].x;
					dU_macro.y += d_bf_macro[j].y * Solution[macro_element->GetDOFInLocalID(j)].y;
					dU_macro.z += d_bf_macro[j].z * Solution[macro_element->GetDOFInLocalID(j)].z;
				}

				double len;
				int id_micro = macro_element->self_grid.GetNearestElementID(Centr_micro, len);
				double eps[6] = { dU_macro.x.x, dU_macro.y.y, dU_macro.z.z, dU_macro.x.y + dU_macro.y.x, dU_macro.y.z + dU_macro.z.y, dU_macro.x.z + dU_macro.z.x };
				double v = solver_grid.GetDomain(macro_element->self_grid.GetElement(id_micro)->GetIdDomain())->forMech.v;
				double E = solver_grid.GetDomain(macro_element->self_grid.GetElement(id_micro)->GetIdDomain())->forMech.GetE(0);
				double a = v / (1 - v);
				double b = (1 - 2 * v) / (2 * (1 - v));
				//double k = E * (1 - v) / ((1 + v)*(1 - v));
				double k = E * (1. - v) / (1. + v) / (1. - 2. * v);
				double sigma[6] = { k * (eps[0] * 1 + eps[1] * a + eps[2] * a), k * (eps[0] * a + eps[1] * 1 + eps[2] * a), k * (eps[0] * a + eps[1] * a + eps[2] * 1),
					k * (2 * eps[3] * b), k * (2 * eps[4] * b), k * (2 * eps[5] * b) };
				double sigma_inv = sqrt((sigma[0] - sigma[1]) * (sigma[0] - sigma[1]) + (sigma[1] - sigma[2]) * (sigma[1] - sigma[2]) + (sigma[2] - sigma[0]) * (sigma[2] - sigma[0]) +
					6 * (sigma[3] * sigma[3] + sigma[4] * sigma[4] + sigma[5] * sigma[5])) / sqrt(2.);
				double eps_inv = sqrt((eps[0] - eps[1]) * (eps[0] - eps[1]) + (eps[1] - eps[2]) * (eps[1] - eps[2]) + (eps[2] - eps[0]) * (eps[2] - eps[0]) +
					3 * (eps[3] * eps[3] + eps[4] * eps[4] + eps[5] * eps[5]) / 2.) * sqrt(2.) / 3;

				value[0][i] = macro_element->self_grid.GetElement(id_micro)->GetIdDomain();
				value[1][i] = sigma[1];
				value[2][i] = sigma[2];

				value[3][i] = eps[0];
				value[4][i] = eps[1];
				value[5][i] = eps[2];

				value[3][i] = U_macro.x;
				value[4][i] = U_macro.y;
				value[5][i] = U_macro.z;

				if (sigma_inv > sigma_inv_max)
				{
					sigma_inv_max = sigma_inv;
					elem_sigma_max = i;
				}
			}
			/*for (int id_DOF = 0; id_DOF < Solution.size(); id_DOF++)
			{
				solver_grid.MoveTheVertex(solver_grid.accordance_DOF_and_vertex[id_DOF], Solution[id_DOF]);
			}*/

			macro_element->self_grid.printTecPlot3D_DiffDomains(fout_tech, value, name_value, name_in_file);

			/*for (int id_DOF = 0; id_DOF < Solution.size(); id_DOF++)
			{
				solver_grid.MoveTheVertex(solver_grid.accordance_DOF_and_vertex[id_DOF], Solution[id_DOF] * (-1));
			}*/
			
			fclose(fout_tech);
		}
		if (false) {
			auto OutputIntoLine = [&](int line_size, Point<double> base_x, Point<bool> variated,
				std::vector<Point<double>>& Line, std::vector<std::vector<double>>& value)
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
					double sigma[6] = { k * (eps[0] * 1 + eps[1] * a + eps[2] * a), k * (eps[0] * a + eps[1] * 1 + eps[2] * a), k * (eps[0] * a + eps[1] * a + eps[2] * 1),
						k * (2 * eps[3] * b), k * (2 * eps[4] * b), k * (2 * eps[5] * b) };
					double sigma_inv = sqrt((sigma[0] - sigma[1]) * (sigma[0] - sigma[1]) + (sigma[1] - sigma[2]) * (sigma[1] - sigma[2]) + (sigma[2] - sigma[0]) * (sigma[2] - sigma[0]) +
						6 * (sigma[3] * sigma[3] + sigma[4] * sigma[4] + sigma[5] * sigma[5])) / sqrt(2.);
					double eps_inv = sqrt((eps[0] - eps[1]) * (eps[0] - eps[1]) + (eps[1] - eps[2]) * (eps[1] - eps[2]) + (eps[2] - eps[0]) * (eps[2] - eps[0]) +
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

				FILE* fout;
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
		//Output into plane 
		if (true)
		{
			math::SimpleGrid geo_grid;
			char mesh_plane_directory[1000];
			sprintf_s(mesh_plane_directory, sizeof(mesh_plane_directory), "%s/Mesh_for_print_XOZ.dat", mesh_directory);
			geo_grid.ReadFromSalomeDat(mesh_plane_directory, 2);
			
			FILE* fout_tech;
			char name_u_tech[5000];
			sprintf_s(name_u_tech, "%s/Plane_XOZ.dat", base_result_directory);
			fopen_s(&fout_tech, name_u_tech, "w");
			char name_in_file[1000];
			sprintf_s(name_in_file, "Elastic_plane_XOZ");
			std::vector<std::vector<char>> name_value(6);
			char name_v_tmp[6][100];
			sprintf_s(name_v_tmp[0], "sigma_xx");
			sprintf_s(name_v_tmp[1], "sigma_yy");
			sprintf_s(name_v_tmp[2], "sigma_zz");
			sprintf_s(name_v_tmp[3], "Ux");
			sprintf_s(name_v_tmp[4], "Uy");
			sprintf_s(name_v_tmp[5], "Uz");
			for (int i = 0; i < name_value.size(); i++)
			{
				name_value[i].resize(100);
				for (int j = 0; j < name_value[i].size(); j++)
				{
					name_value[i][j] = name_v_tmp[i][j];
				}
			}
			std::vector<std::vector<double>> value(3 * 2);
			value[0].resize(geo_grid.nvtr.size());
			value[1].resize(geo_grid.nvtr.size());
			value[2].resize(geo_grid.nvtr.size());
			value[3].resize(geo_grid.nvtr.size());
			value[4].resize(geo_grid.nvtr.size());
			value[5].resize(geo_grid.nvtr.size());
			double sigma_inv_max = 0;
			int elem_sigma_max = 0;
			for (int id_elem = 0; id_elem < geo_grid.nvtr.size(); id_elem++)
			{
				Point<double> Target_point;
				for (int i = 0; i < geo_grid.nvtr[id_elem].size(); i++)
				{
					Target_point += geo_grid.xyz[geo_grid.nvtr[id_elem][i]];
				}
				Target_point /= geo_grid.nvtr[id_elem].size();

				Point<double> U_macro;
				Point<Point<double>> dU_macro;

				double len;
				int id_macro = solver_grid.GetNearestElementID(Target_point, len);
				auto macro_element = solver_grid.GetElement(id_macro);
				int id_micro = macro_element->self_grid.GetNearestElementID(Target_point, len);

				std::vector<Point<double>> bf_macro(macro_element->GetDOFsCount());
				std::vector<Point<Point<double>>> d_bf_macro(macro_element->GetDOFsCount());
				for (int j = 0; j < bf_macro.size(); j++)
				{
					bf_macro[j] = macro_element->self_grid.GetSolutionInPoint(id_micro, Target_point, macro_element->self_basis_functions[j]);
					d_bf_macro[j] = macro_element->self_grid.GetDerevativeFromSolutionInPoint(id_micro, Target_point, macro_element->self_basis_functions[j]);

					Point<double> val;
					val.x = Solution[macro_element->GetDOFInLocalID(j)].x;
					val.y = Solution[macro_element->GetDOFInLocalID(j)].y;
					val.z = Solution[macro_element->GetDOFInLocalID(j)].z;

					U_macro.x += bf_macro[j].x * val.x;
					U_macro.y += bf_macro[j].y * val.y;
					U_macro.z += bf_macro[j].z * val.z;

					dU_macro.x += d_bf_macro[j].x * val.x;
					dU_macro.y += d_bf_macro[j].y * val.y;
					dU_macro.z += d_bf_macro[j].z * val.z;
				}


				double eps[6] = { dU_macro.x.x, dU_macro.y.y, dU_macro.z.z, dU_macro.x.y + dU_macro.y.x, dU_macro.y.z + dU_macro.z.y, dU_macro.x.z + dU_macro.z.x };
				double v = solver_grid.GetDomain(macro_element->self_grid.GetElement(id_micro)->GetIdDomain())->forMech.v;
				double E = solver_grid.GetDomain(macro_element->self_grid.GetElement(id_micro)->GetIdDomain())->forMech.GetE(0);
				double a = v / (1 - v);
				double b = (1 - 2 * v) / (2 * (1 - v));
				//double k = E * (1 - v) / ((1 + v)*(1 - v));
				double k = E * (1. - v) / (1. + v) / (1. - 2. * v);
				double sigma[6] = { k * (eps[0] * 1 + eps[1] * a + eps[2] * a), k * (eps[0] * a + eps[1] * 1 + eps[2] * a), k * (eps[0] * a + eps[1] * a + eps[2] * 1),
					k * (2 * eps[3] * b), k * (2 * eps[4] * b), k * (2 * eps[5] * b) };
				double sigma_inv = sqrt((sigma[0] - sigma[1]) * (sigma[0] - sigma[1]) + (sigma[1] - sigma[2]) * (sigma[1] - sigma[2]) + (sigma[2] - sigma[0]) * (sigma[2] - sigma[0]) +
					6 * (sigma[3] * sigma[3] + sigma[4] * sigma[4] + sigma[5] * sigma[5])) / sqrt(2.);
				double eps_inv = sqrt((eps[0] - eps[1]) * (eps[0] - eps[1]) + (eps[1] - eps[2]) * (eps[1] - eps[2]) + (eps[2] - eps[0]) * (eps[2] - eps[0]) +
					3 * (eps[3] * eps[3] + eps[4] * eps[4] + eps[5] * eps[5]) / 2.) * sqrt(2.) / 3;

				value[0][id_elem] = macro_element->self_grid.GetElement(id_micro)->GetIdDomain();
				value[1][id_elem] = sigma[1];
				value[2][id_elem] = sigma[2];

				value[3][id_elem] = eps[0];
				value[4][id_elem] = eps[1];
				value[5][id_elem] = eps[2];

				value[3][id_elem] = U_macro.x;
				value[4][id_elem] = U_macro.y;
				value[5][id_elem] = U_macro.z;
			}

			geo_grid.printTecPlot3D(fout_tech, value, name_value, name_in_file);

			FILE* fout_U;
			char name_out_U[1000];
			sprintf_s(name_out_U, "%s/Solution_in_plane.dat", base_result_directory);
			fopen_s(&fout_U, name_out_U, "w");
			fprintf_s(fout_U, "%d\n", value[0].size());
			for (int i = 0; i < value[0].size(); i++)
			{
				fprintf_s(fout_U, "%.8e %.8e %.8e\n", value[3][i], value[4][i], value[5][i]);
			}
			fclose(fout_U);
		}
		printf_s("\t complite\n");

		FILE* fout_U;
		char name_out_U[1000];
		sprintf_s(name_out_U, "%s/Solution.dat", base_result_directory);
		fopen_s(&fout_U, name_out_U, "w");
		for (int i = 0; i < Solution.size(); i++)
		{
			fprintf_s(fout_U, "%.8e %.8e %.8e\n", Solution[i].x, Solution[i].y, Solution[i].z);
		}
		fclose(fout_U);

		solver_grid.~Grid_forMech();
	}
}

int main()
{
	//CreateCubGridAsPoly_dat_format();
	
	Solve_ElasticDeformationProblem();

	int a;
	scanf_s("%d", &a);
}