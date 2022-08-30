#include "../FEM/FEM.h"
#include <Windows.h>
#include <iostream>

void EffectiveElastisity_MSH(char properties_file[1000])
{
	FEM::Grid_forMech solver_grid; //output
	std::vector<Point<double>> Solution; //output

	//char properties_file[1000] = { "E:/+cyl/800el/param.txt" };
	//char properties_file[1000] = { "E:/Box/100x10x200/BoxWithCracks/67k/param.txt" };
	//char properties_file[1000] = { "E:/Box/100x10x200/BoxWithCracks/638k/param.txt" };
	//char properties_file[1000] = { "E:/Box/100x50x200/param.txt" };
	//char properties_file[1000] = { "E:/Box/200x200x200/200k_fem/param_for_solver.txt" };
	//char properties_file[1000] = { "D:/SeismicGroup/sand_init_mech/m1_mech_porosity/param_for_solver.txt" };
	//char properties_file[1000] = { "D:/SeismicGroup/simpl_model/param_for_solver.txt" };
	//char properties_file[1000] = { "D:/SeismicGroup/simpl_model/10spheres/param_for_solver.txt" };
	//char properties_file[1000] = { "base_properties.txt" };
	printf_s("The properties file contains the information:\n==========================================\n");
	printf_s("\tDo you need to print to a log file? 0-false/1-true\n");
	printf_s("\tinput: Mesh directory: box_200x200x200/90el\n");
	printf_s("\tproperties: Iterative process:\n");
	printf_s("\t\t<Critical residual value for SLAE solution>\n");
	printf_s("\t<properties: Boundary values>\n");
	printf_s("\t<properties: Materials>\n");
	printf_s("\t<output: Result directory: box_200x200x200/90el/result>\n==========================================\n");


	//printf_s("Enter the name of the properties file: ");
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
	std::vector<std::vector<std::vector<int>>> boundary_faces;
	geo_grid.ReadFromMSH(mesh_directory, 1, 4, 4, boundary_faces);

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
		/*Point<double> value_const;
		Point<bool> is_condition;*/
		std::vector<int> id_vertexes;
		std::function<Point<double>(Point<bool>&, int)> value;
	};
	std::vector<_Dirichlet> first_boundaries;
	{
		int Nb;
		std::vector<int> id_boundaries;
		{
			char _line[1000];
			math::ReadNonEmptyLine_forNumbers(f_properties, _line);
			std::vector<int> val;
			math::ParserStringToVectorInt(_line, val, " ");
			Nb = val[0];
		}

		first_boundaries.resize(Nb);
		id_boundaries.resize(Nb);
		char line[1000];
		for (int i = 0; i < Nb; i++)
		{
			Point<double> value;
			Point<bool> is_condition;

			math::ReadNonEmptyLine_forNumbers(f_properties, line);

			if (line[0] != '#')
			{
				id_boundaries[i] = math::ParserCharToInt(line[0]);
				int jj = 0;
				int curr_i = 0;
				char _tmp_line[1000];
				for (int j = 1; j < 1000; j++)
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
			}

			first_boundaries[i].value = [value, is_condition](Point<bool>& is_take, int id)->Point<double>
			{
				is_take = is_condition;
				return value;
			};
		}

		if (Nb != 0)
		{
			for (int id_type = 0; id_type < Nb; id_type++)
			{
				std::vector<int> tmp_vert;
				for (int i = 0; i < boundary_faces[id_boundaries[id_type]].size(); i++)
				{
					for (int j = 1; j < boundary_faces[id_boundaries[id_type]][i].size(); j++)
					{
						tmp_vert.push_back(boundary_faces[id_boundaries[id_type]][i][j]);
					}
				}
				math::MakeQuickSort(tmp_vert);
				math::MakeRemovalOfDuplication(tmp_vert, first_boundaries[id_type].id_vertexes);
			}
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
		std::vector<int> id_boundaries(Nb);
		char line[1000];
		std::vector<bool> is_individual_values(Nb);
		std::vector<double> individual_values(Nb);
		for (int i = 0; i < Nb; i++)
		{
			is_individual_values[i] = false;
			double value;
			Point<double> vector;

			math::ReadNonEmptyLine_forNumbers(f_properties, line);
			id_boundaries[i] = math::ParserCharToInt(line[0]);
			int jj = 0;
			int curr_i = 0;
			char _tmp_line[1000];
			for (int j = 1; j < 1000; j++)
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
			for (int id_type = 0; id_type < Nb; id_type++)
			{
				second_boundary[id_type].id_vertexes_as_triangle.resize(boundary_faces[id_boundaries[id_type]].size());
				for (int id_triang = 0; id_triang < boundary_faces[id_boundaries[id_type]].size(); id_triang++)
				{
					printf_s("Create boundary condition %d (triangle %d/%d)\r", id_boundaries[id_type], id_triang, boundary_faces[id_boundaries[id_type]].size());
					int test_vertex = -1;
					int base_elem = boundary_faces[id_boundaries[id_type]][id_triang][0];

					for (int i = 1; i < boundary_faces[id_boundaries[id_type]][id_triang].size(); i++)
						second_boundary[id_type].id_vertexes_as_triangle[id_triang].push_back(boundary_faces[id_boundaries[id_type]][id_triang][i]);
					second_boundary[id_type].id_base_element.push_back(base_elem);

					//по собственным нормалям
					if (is_individual_values[id_type] == true)
					{
						double A, B, C, D;
						Point<double> test_vector = geo_grid.xyz[test_vertex];
						std::vector<Point<double>> vertexes(3);
						vertexes[0] = geo_grid.xyz[second_boundary[id_type].id_vertexes_as_triangle[id_triang][0]];
						vertexes[1] = geo_grid.xyz[second_boundary[id_type].id_vertexes_as_triangle[id_triang][1]];
						vertexes[2] = geo_grid.xyz[second_boundary[id_type].id_vertexes_as_triangle[id_triang][2]];
						bool reverse = false;
						math::GetPlaneEquation(vertexes, A, B, C, D);

						if (Point<double>(A, B, C) * test_vector > 0)
						{
							int t = second_boundary[id_type].id_vertexes_as_triangle[id_triang][1];
							second_boundary[id_type].id_vertexes_as_triangle[id_triang][1] = second_boundary[id_type].id_vertexes_as_triangle[id_triang][2];
							second_boundary[id_type].id_vertexes_as_triangle[id_triang][2] = t;

							vertexes[0] = geo_grid.xyz[second_boundary[id_type].id_vertexes_as_triangle[id_triang][0]];
							vertexes[1] = geo_grid.xyz[second_boundary[id_type].id_vertexes_as_triangle[id_triang][1]];
							vertexes[2] = geo_grid.xyz[second_boundary[id_type].id_vertexes_as_triangle[id_triang][2]];
						}
						Point<double> normal;
						double d;
						math::GetPlaneEquation(vertexes, normal.x, normal.y, normal.z, d);
						normal /= math::SolveLengthVector(normal);

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

		std::vector<double> volumes_before(_materials.size()), volumes_after(_materials.size());
		double obj_volume_before = 0, obj_volume_after = 0;
		for (int j = 0; j < geo_grid.nvkat.size(); j++)
		{
			double vol = geometry::Tetrahedron::GetVolume(
				geo_grid.xyz[geo_grid.nvtr[j][0]],
				geo_grid.xyz[geo_grid.nvtr[j][1]],
				geo_grid.xyz[geo_grid.nvtr[j][2]],
				geo_grid.xyz[geo_grid.nvtr[j][3]]);
			obj_volume_before += vol;
			volumes_before[geo_grid.nvkat[j]] += vol;
		}

		FEM::FEM_forElasticDeformation(
			is_print_logFile,
			critical_residual,
			geo_grid, //input
			first_boundaries, //input
			second_boundary, //input
			base_result_directory, //output
			solver_grid, //output
			Solution //output
		);

		//изменение объема от нагрузки
		if (true)
		{
			double full_volume_by_part_before = 0;// solver_grid.GetFullVolumeByPart();

			for (int id_DOF = 0; id_DOF < Solution.size(); id_DOF++)
			{
				solver_grid.MoveTheVertex(solver_grid.accordance_DOF_and_vertex[id_DOF], Solution[id_DOF]);
			}
			for (int j = 0; j < solver_grid.GetElementsCount(); j++)
			{
				double vol = solver_grid.GetElement(j)->GetVolume();
				obj_volume_after += vol;
				volumes_after[solver_grid.GetElement(j)->GetIdDomain()] += vol;
			}

			double full_volume_by_part_after = 0;// solver_grid.GetFullVolumeByPart();

			for (int id_DOF = 0; id_DOF < Solution.size(); id_DOF++)
			{
				solver_grid.MoveTheVertex(solver_grid.accordance_DOF_and_vertex[id_DOF], Solution[id_DOF] * (-1));
			}

			FILE* fout_volume;
			char name_vol[5000];
			sprintf_s(name_vol, "%s/Volume_result.txt", base_result_directory);
			fopen_s(&fout_volume, name_vol, "w");

			fprintf_s(fout_volume, "\n----exterior forces (%s)----\n", properties_file);
			for (int i = 0; i < first_boundaries.size(); i++)
			{
				Point<bool> is_use;
				auto value = first_boundaries[i].value(is_use, 0);
				fprintf_s(fout_volume, "First condition number is %d: value(", i);
				if (is_use.x) fprintf_s(fout_volume, "%.2e, ", value.x);
				else fprintf_s(fout_volume, "*, ", value.x);
				if (is_use.y) fprintf_s(fout_volume, "%.2e, ", value.y);
				else fprintf_s(fout_volume, "*, ", value.y);
				if (is_use.z) fprintf_s(fout_volume, "%.2e)\n", value.z);
				else fprintf_s(fout_volume, "*)\n", value.z);
			}
			for (int i = 0; i < second_boundary.size(); i++)
			{
				auto value = second_boundary[i].value(Point<double>(0, 0, 0));
				fprintf_s(fout_volume, "Second condition is %d: value(%.2e, %.2e, %.2e)\n", i, value.x, value.y, value.z);
			}

			fprintf_s(fout_volume, "\n----result before deformations----\n");
			fprintf_s(fout_volume, "Full volume %.6e (by part %.6e)\n", obj_volume_before, full_volume_by_part_before);
			for (int i = 0; i < volumes_before.size(); i++)
			{
				if (volumes_before[i] / obj_volume_before * 100 > 0.01)
					fprintf_s(fout_volume, "%d domain: %.6e (%.6lf%%)\n", i, volumes_before[i], volumes_before[i] / obj_volume_before * 100);
				else
					fprintf_s(fout_volume, "%d domain: %.6e (%.6e%%)\n", i, volumes_before[i], volumes_before[i] / obj_volume_before * 100);
			}

			fprintf_s(fout_volume, "\n----result after deformations----\n");
			fprintf_s(fout_volume, "Full volume %.6e (by part %.6e)\n", obj_volume_after, full_volume_by_part_after);
			for (int i = 0; i < volumes_after.size(); i++)
			{
				if (volumes_after[i] / obj_volume_after * 100 > 0.01)
					fprintf_s(fout_volume, "%d domain: %.6e (%.6lf%%)\n", i, volumes_after[i], volumes_after[i] / obj_volume_after * 100);
				else
					fprintf_s(fout_volume, "%d domain: %.6e (%.6e%%)\n", i, volumes_after[i], volumes_after[i] / obj_volume_after * 100);
			}
			fclose(fout_volume);
		}

		//output solution
		printf_s("Print the mech result into .dat file... ");
		//with deformations
		if (true) {
			FILE* fout_tech;
			char name_u_tech[5000];
			sprintf_s(name_u_tech, "%s/U_deformations.dat", base_result_directory);
			fopen_s(&fout_tech, name_u_tech, "w");
			char name_in_file[1000];
			sprintf_s(name_in_file, "Elastic");
			std::vector<std::vector<char>> name_value(6);
			char name_v_tmp[6][100];
			sprintf_s(name_v_tmp[0], "sigma_mises");
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
				double v = solver_grid.GetDomain(element->GetIdDomain())->forMech.v;
				double E = solver_grid.GetDomain(element->GetIdDomain())->forMech.GetE(0);
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

				if (E > 0)
				{
					value[0][i] = sigma_inv;
					value[1][i] = sigma[1];
					value[2][i] = sigma[2];

					value[3][i] = eps[0];
					value[4][i] = eps[1];
					value[5][i] = eps[2];
				}
				else {
					value[0][i] = 0;
					value[1][i] = 0;
					value[2][i] = 0;
				}

				if (value[0][i] != value[0][i]) value[0][i] = 0;
				if (value[1][i] != value[1][i]) value[1][i] = 0;
				if (value[2][i] != value[2][i]) value[2][i] = 0;

				value[3][i] = U.x;
				value[4][i] = U.y;
				value[5][i] = U.z;

				if (sigma_inv > sigma_inv_max)
				{
					sigma_inv_max = sigma_inv;
					elem_sigma_max = i;
				}
			}
			for (int id_DOF = 0; id_DOF < Solution.size(); id_DOF++)
			{
				solver_grid.MoveTheVertex(solver_grid.accordance_DOF_and_vertex[id_DOF], Solution[id_DOF]);
			}
			solver_grid.printTecPlot3D_DiffDomains(fout_tech, value, name_value, name_in_file);
			for (int id_DOF = 0; id_DOF < Solution.size(); id_DOF++)
			{
				solver_grid.MoveTheVertex(solver_grid.accordance_DOF_and_vertex[id_DOF], Solution[id_DOF] * (-1));
			}
			fclose(fout_tech);
		}
		//without deformations
		if (true) {
			FILE* fout_tech;
			char name_u_tech[5000];
			sprintf_s(name_u_tech, "%s/U_non_deformation.dat", base_result_directory);
			fopen_s(&fout_tech, name_u_tech, "w");
			char name_in_file[1000];
			sprintf_s(name_in_file, "Elastic");
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
				double v = solver_grid.GetDomain(element->GetIdDomain())->forMech.v;
				double E = solver_grid.GetDomain(element->GetIdDomain())->forMech.GetE(0);
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

				if (E > 0)
				{
					value[0][i] = sigma[0];
					value[1][i] = sigma[1];
					value[2][i] = sigma[2];

					value[3][i] = eps[0];
					value[4][i] = eps[1];
					value[5][i] = eps[2];
				}
				else {
					value[0][i] = 0;
					value[1][i] = 0;
					value[2][i] = 0;
				}

				if (value[0][i] != value[0][i]) value[0][i] = 0;
				if (value[1][i] != value[1][i]) value[1][i] = 0;
				if (value[2][i] != value[2][i]) value[2][i] = 0;


				value[3][i] = U.x;
				value[4][i] = U.y;
				value[5][i] = U.z;

				if (sigma_inv > sigma_inv_max)
				{
					sigma_inv_max = sigma_inv;
					elem_sigma_max = i;
				}
			}

			solver_grid.printTecPlot3D_DiffDomains(fout_tech, value, name_value, name_in_file);
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
		printf_s("\t complite\n");

		solver_grid.~Grid_forMech();
	}
}

void main()
{
	EffectiveElastisity_MSH();
}