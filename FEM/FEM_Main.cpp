#pragma once
#include "FEM.h"
#include <Windows.h>
#include <iostream>

double EffectiveResistance_new(FEM::Grid_forScal &grid, std::vector<double> &U, double &norm_Jsqrt, double &rpho_eff_summ, double &GLOBAL_SIGMA_hmm)
{
	//вычисляем сразу ||J|| и ||gardU||
	double norm_J = 0;
	norm_Jsqrt = 0;
	double norm_gradU = 0;

	GLOBAL_SIGMA_hmm = 0;
	double norm_gradUsqrt_on_points = 0;
	double norm_Jsqrt_on_points = 0;

	double V_global = 0;
	for (int elem_coarse = 0; elem_coarse < grid.GetElementsCount(); elem_coarse++)
	{
		V_global += grid.GetElement(elem_coarse)->GetVolume();
	}

	std::vector<double> macro_norm_J_master;
	std::vector<double> macro_norm_Jsqrt_master;
	std::vector<double> macro_norm_gradU_master;
	std::vector<double> macro_norm_Jsqrt_on_points_master;
	std::vector<double> macro_norm_gradUsqrt_on_points_master;
	for (int elem_coarse = 0; elem_coarse < grid.GetElementsCount(); elem_coarse++)
	{
		auto element = grid.GetElement(elem_coarse);
		element->SetIntegrationLaw(5);

		std::function<double(Point<double>)> norm_J_fem = [&](Point<double> X) ->double
		{
			Point<double> J = grid.GetDerevativeFromSolutionInPoint(elem_coarse, X, U);
			double sigma;
			sigma = grid.GetDomain(element->GetIdDomain())->forThermal.lambda;
			return sigma * sigma*(J.x*J.x + J.y*J.y + J.z*J.z);
		};
		std::function<double(Point<double>)> norm_Jsqrt_fem = [&](Point<double> X) ->double
		{
			Point<double> J = grid.GetDerevativeFromSolutionInPoint(elem_coarse, X, U);
			double sigma;
			sigma = grid.GetDomain(element->GetIdDomain())->forThermal.lambda;
			return sigma * sqrt(J.x*J.x + J.y*J.y + J.z*J.z);
		};
		std::function<double(Point<double>)> norm_gradU_fem = [&](Point<double> X) ->double
		{
			Point<double> gradU = grid.GetDerevativeFromSolutionInPoint(elem_coarse, X, U);
			return (gradU.x*gradU.x + gradU.y*gradU.y + gradU.z*gradU.z);
		};
		std::function<double(Point<double>)> norm_gradUsqrt_fem = [&](Point<double> X) ->double
		{
			Point<double> gradU = grid.GetDerevativeFromSolutionInPoint(elem_coarse, X, U);
			return sqrt(gradU.x*gradU.x + gradU.y*gradU.y + gradU.z*gradU.z);
		};

		std::function<double(Point<double>)> func = [&](Point<double> X) ->double
		{
			return 1;
		};

		double tmpJ = element->SolveIntegral(norm_J_fem);
		double tmpU = element->SolveIntegral(norm_gradU_fem);
		double tmpJsqrt = element->SolveIntegral(norm_Jsqrt_fem);
		double integr = element->SolveIntegral(func);

		double tmpJsqrt_on_points = norm_Jsqrt_fem(element->GetWeightCentr());//grid.elem[elem_coarse].integrate(1, norm_Jsqrt_fem);
		double tmpUsqrt_on_points = norm_gradUsqrt_fem(element->GetWeightCentr());//grid.elem[elem_coarse].integrate(1, norm_gradUsqrt_fem);

		norm_Jsqrt += tmpJsqrt;
		norm_J += tmpJ;
		norm_gradU += tmpU;
		norm_Jsqrt_on_points += tmpJsqrt_on_points;
		norm_gradUsqrt_on_points += tmpUsqrt_on_points;
	}

	norm_J = sqrt(abs(norm_J));
	printf("\n----------->\nnorm_J = %.15lf\n<-----------\n", norm_J);
	norm_gradU = sqrt(abs(norm_gradU));
	printf("\n----------->\nnorm_gradU = %.15lf\n<-----------\n", norm_gradU);

	printf("\n----------->\nI = %.15lf\n<-----------\n", norm_Jsqrt);

	rpho_eff_summ = norm_gradUsqrt_on_points / norm_Jsqrt_on_points;

	return (norm_J/norm_gradU);
}

void Solve_TermalProblem_inTime()
{
	FEM::Grid_forScal solver_grid; //output
	std::vector<double> Solution; //output

	printf_s("\n====================TermalProblem_inTime======================\n");
	char properties_file[1000] = { "D:/Documents/universiti/MySOFT/TEST/box_200x200x200_craсks/800el_1_trmaal/p.txt" };
	//char properties_file[1000] = { "base_properties.txt" };
	printf_s("The properties file contains the information:\n==========================================\n");
	printf_s("\tDo you need to print to a log file? 0-false/1-true\n");
	printf_s("\tinput: Mesh directory: box_200x200x200/90el\n");
	printf_s("\tproperties: Iterative process:\n");
	printf_s("\t\t<Start iteration> <End iteration>\n");
	printf_s("\t\t<Step size for Time>\n");
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
	
	int start_iteration, end_iteration;
	{
		char _line[1000];
		math::ReadNonEmptyLine_forNumbers(f_properties, _line);
		std::vector<int> val;
		math::ParserStringToVectorInt(_line, val, " ");
		start_iteration = val[0];
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
		std::function<double(Point<double>X)> value;
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
			double value;

			math::ReadNonEmptyLine_forNumbers(f_properties, line);
			std::vector<float> val;
			math::ParserStringToVectorFloat(line, val, " ");
			value = val[0];

			first_boundaries[i].value = [value](Point<double> X) -> double
			{
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
		//double value;
		//Point<double> vector;
		std::function<Point<double>(Point<double>)> value;
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
		for (int i = 0; i < Nb; i++)
		{
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

			second_boundary[i].value = [value, vector](Point<double> X) -> Point<double>
			{
				Point <double> res;
				res.x = vector.x *value;
				res.y = vector.y *value;
				res.z = vector.z *value;
				return res;
			};
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
						second_boundary[id_type].id_vertexes_as_triangle.push_back(_tmp);
						second_boundary[id_type].id_base_element.push_back(id_element);
					}
				}
			}
		}
	}

	//read Materials
	struct _material
	{
		double _c;
		double _lambda;
		double _rpho;

		_material(double _c, double _lambda, double _rpho)
		{
			this->_c = _c;
			this->_lambda = _lambda;
			this->_rpho = _rpho;
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
		double _c;
		double _lambda;
		double _rpho;
		{
			char _line[1000];
			math::ReadNonEmptyLine_forNumbers(f_properties, _line);
			std::vector<float> val;
			math::ParserStringToVectorFloat(_line, val, " ");
			_lambda = val[0];
			_c = val[1];
			_rpho = val[2];
		}
		_materials.push_back(_material(_c, _lambda, _rpho));
	}

	math::ReadNonEmptyLine(f_properties, base_result_directory);
	fclose(f_properties);

	CreateDirectory((LPCTSTR)base_result_directory, NULL);

	if (start_iteration != 0)
	{
		sprintf_s(cracks_directory, sizeof(cracks_directory), "%s/STEP_%d", base_result_directory, start_iteration);
	}

	bool res = false;
	double TIME_h = step_size;
	std::vector<double> T_prev(geo_grid.xyz.size());
	Solution.resize(geo_grid.xyz.size());
	math::InitializationVector(T_prev, 5);
	math::InitializationVector(Solution, 5);
	for (int id_STEP = start_iteration; id_STEP < end_iteration; id_STEP++)
	{
		printf_s("\n=============================*===============================\n");
		printf_s("================= Start solution of %d STEP (time = %.2lf) ================\n", id_STEP, TIME_h*id_STEP);
		//Create new directories
		char result_directory[1000];
		sprintf_s(result_directory, sizeof(result_directory), "%s/STEP_%d", base_result_directory, id_STEP);
		CreateDirectory((LPCTSTR)result_directory, NULL);

		for (int i = 0; i < _materials.size(); i++)
		{
			solver_grid.AddDomain();
			auto domain = solver_grid.GetDomain(i);
			domain->forThermal.c = _materials[i]._c;
			domain->forThermal.lambda = _materials[i]._lambda;
			domain->forThermal.rpho = _materials[i]._rpho;
		}

		std::function<double(int, Point<double>)> StiffnessCoef = [&](int elem, Point<double> X)->double
		{
			//return (20E-7)*(20E-7);

			return solver_grid.GetDomain(solver_grid.GetElement(elem)->GetIdDomain())->forThermal.lambda;
		};
		std::function<double(int, Point<double>)> MassCoef = [&](int elem, Point<double> X)->double
		{
			//return 1/TIME_h;

			auto domain = (solver_grid.GetDomain(solver_grid.GetElement(elem)->GetIdDomain()));
			return domain->forThermal.c * domain->forThermal.rpho / TIME_h;
		};
		std::function<double(int, Point<double>)> F = [&](int elem, Point<double> X)->double
		{
			double _T = solver_grid.GetSolutionInPoint(elem, X, T_prev);

			//return +1 * _T / TIME_h;


			auto domain = (solver_grid.GetDomain(solver_grid.GetElement(elem)->GetIdDomain()));
			return +1 * _T* domain->forThermal.c* domain->forThermal.rpho / TIME_h;
		};

		FEM::FEM_forBaseEliptic(
			is_print_logFile,
			critical_residual,
			geo_grid, //input
			first_boundaries, //input
			second_boundary, //input
			StiffnessCoef,
			MassCoef,
			F,
			result_directory, //output
			solver_grid, //output
			Solution //output
		);

		printf("Solve effective resistans... \n");
		{
			double R = 0, R_new = 0, norm_Jsqrt = 0, rpho_eff_summ, GLOBAL_SIGMA_hmm;
			R_new = EffectiveResistance_new(solver_grid, Solution, norm_Jsqrt, rpho_eff_summ, GLOBAL_SIGMA_hmm);
			printf("\n\n----------------------------\n----------------------------\nR_new = %lf // %.6e\n----------------------------\n----------------------------\n\n", R_new, R_new);
			FILE *fout;
			FILE *fout_new;
			char name_out[5000];
			sprintf_s(name_out, "%s/Effectiev properties.txt", result_directory);
			fopen_s(&fout, name_out, "a");

			fprintf_s(fout, "\n\n----------------------------\n");
			for (int i = 0; i < solver_grid.GetDomainsCount(); i++)
			{
				fprintf_s(fout, "lambda_%d = %.5e\n", i, solver_grid.GetDomain(i)->forThermal.lambda);
				fprintf_s(fout, "c_%d = %.5e\n", i, solver_grid.GetDomain(i)->forThermal.c);
				fprintf_s(fout, "rpho_%d = %.5e\n", i, solver_grid.GetDomain(i)->forThermal.rpho);
			}
			fprintf_s(fout, "lambda_eff_new = %.15lf\n", R_new);

			double V = 0;
			std::vector<double> Vin(solver_grid.GetDomainsCount());
			int qwe = 0;
			double conc = 0.0;

			for (int i = 0; i < solver_grid.GetElementsCount(); i++)
			{
				auto element = solver_grid.GetElement(i);
				Vin[element->GetIdDomain()] += element->GetVolume();
				V += element->GetVolume();
			}
			for (int i = 0; i < Vin.size(); i++)
			{
				printf("concentration(id%d) = %.5lf %%\n", i, (Vin[i] / V) * 100.);
				fprintf_s(fout, "concentration(id%d) = %.5lf %%\n", i, (Vin[i] / V) * 100.);
			}
			fclose(fout);
		}

		//output solution
		printf_s("Print the mech result into .dat file... ");
		{
			FILE *fout_tech;
			char name_u_tech[5000];
			sprintf_s(name_u_tech, "%s/T.dat", result_directory);
			fopen_s(&fout_tech, name_u_tech, "w");
			char name_in_file[1000];
			sprintf_s(name_in_file, "Termal_in_time_step%d", id_STEP);
			std::vector<std::vector<char>> name_value(5);
			char name_v_tmp[5][100];
			sprintf_s(name_v_tmp[0], "T");
			sprintf_s(name_v_tmp[1], "gradT_x");
			sprintf_s(name_v_tmp[2], "gradT_y");
			sprintf_s(name_v_tmp[3], "gradT_z");
			sprintf_s(name_v_tmp[4], "material");
			for (int i = 0; i < name_value.size(); i++)
			{
				name_value[i].resize(100);
				for (int j = 0; j < name_value[i].size(); j++)
				{
					name_value[i][j] = name_v_tmp[i][j];
				}
			}
			std::vector<std::vector<double>> value(5);
			value[0].resize(solver_grid.GetElementsCount());
			value[1].resize(solver_grid.GetElementsCount());
			value[2].resize(solver_grid.GetElementsCount());
			value[3].resize(solver_grid.GetElementsCount());
			value[4].resize(solver_grid.GetElementsCount());
			double sigma_inv_max = 0;
			int elem_sigma_max = 0;
			for (int i = 0; i < solver_grid.GetElementsCount(); i++)
			{
				auto element = solver_grid.GetElement(i);

				Point<double> Centr = element->GetWeightCentr();

				double U = solver_grid.GetSolutionInPoint(i, Centr, Solution);
				Point<double> dU = solver_grid.GetDerevativeFromSolutionInPoint(i, Centr, Solution);

				value[0][i] = U;
				value[1][i] = dU.x;
				value[2][i] = dU.y;
				value[3][i] = dU.z;
				value[4][i] = element->GetIdDomain();
			}
			solver_grid.printTecPlot3D(fout_tech, value, name_value, name_in_file);

			fclose(fout_tech);
		}
		printf_s("\t complite\n");

		math::MakeCopyVector_A_into_B(Solution, T_prev);
		solver_grid.~Grid_forScal();
	}
}

void Solve_ElasticDeformationProblem()
{
	FEM::Grid_forMech solver_grid; //output
	std::vector<Point<double>> Solution; //output


	clock_t t_after = clock();
	double start = omp_get_wtime();

	//char properties_file[1000] = { "E:/+cyl/800el/param.txt" };
	//char properties_file[1000] = { "E:/Box/100x10x200/BoxWithCracks/67k/param.txt" };
	//char properties_file[1000] = { "E:/Box/100x10x200/BoxWithCracks/638k/param.txt" };
	//char properties_file[1000] = { "E:/Box/100x50x200/param.txt" };
	//char properties_file[1000] = { "E:/Box/200x200x200/200k_fem/param_for_solver.txt" };
	//char properties_file[1000] = { "D:/Babenko_2022/1cracks/FEM/75/param_for_solver.txt" };
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

	FILE *f_properties;
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

	{
		char _line[1000];
		math::ReadNonEmptyLine_forNumbers(f_properties, _line);
		std::vector<int> val;
		math::ParserStringToVectorInt(_line, val, " ");
		math::NUM_THREADS = val[0];
	}

	math::ReadNonEmptyLine(f_properties, mesh_directory);
	math::SimpleGrid geo_grid; //input
	geo_grid.ReadFromNVTR(mesh_directory, 4);

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
		std::function<Point<double>(Point<bool>&, int)> value;
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

			first_boundaries[i].value = [value, is_condition](Point<bool> &is_take, int id) -> Point<double>
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
					if(id_element%100)
					printf_s("create boundary elem[%d-%d]\r", id_element, geo_grid.nvtr.size());
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
		for (int i = 0; i < _materials.size(); i++)
		{
			solver_grid.AddDomain();
			auto domain = solver_grid.GetDomain(i);
			domain->forMech.SetE(_materials[i]._E);
			domain->forMech.SetV(_materials[i]._v);
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

		double end = omp_get_wtime();
		printf("start = %.16g\nend = %.16g\ndiff = %.16g\n", start, end, end - start);

		FILE* ftime;
		char ftime_name[1000];
		sprintf_s(ftime_name, "%s/time_result.txt", base_result_directory);
		fopen_s(&ftime, ftime_name, "w");
		fprintf_s(ftime, "start = %.16g\nend = %.16g\ndiff = %.16g\n", start, end, end - start);
		fclose(ftime);

		FILE* fres;
		char fres_name[1000];
		sprintf_s(fres_name, "%s/U_result.dat", base_result_directory);
		fopen_s(&fres, fres_name, "w");
		for (int i = 0; i < Solution.size(); i++)
			fprintf_s(fres, "%.15e %.15e %.15e\n", Solution[i].x, Solution[i].y, Solution[i].z);
		fclose(fres);

		//output solution
		printf_s("Print the mech result into .dat file... ");
		//with deformations
		if(true){
			FILE *fout_tech;
			char name_u_tech[5000];
			sprintf_s(name_u_tech, "%s/U_deformations.dat", base_result_directory);
			fopen_s(&fout_tech, name_u_tech, "w");
			char name_in_file[1000];
			sprintf_s(name_in_file, "Elastic");
			std::vector<std::vector<char>> name_value(6);
			char name_v_tmp[6][100];
			sprintf_s(name_v_tmp[0], "sigma_xx");
			sprintf_s(name_v_tmp[1], "sigma_yy");
			sprintf_s(name_v_tmp[2], "sigma_zz");
			sprintf_s(name_v_tmp[3], "VonMises_stress");
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
			solver_grid.printTecPlot3D(fout_tech, value, name_value, name_in_file);
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
				double v = solver_grid.GetDomain(0)->forMech.v;
				double E = solver_grid.GetDomain(0)->forMech.GetE(0);
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

				value[0][i] = sigma[0];
				value[1][i] = sigma[1];
				value[2][i] = sigma[2];

				value[3][i] = eps[0];
				value[4][i] = eps[1];
				value[5][i] = eps[2];

				value[3][i] = U.x;
				value[4][i] = U.y;
				value[5][i] = U.z;

				if (sigma_inv > sigma_inv_max)
				{
					sigma_inv_max = sigma_inv;
					elem_sigma_max = i;
				}
			}
			
			solver_grid.printTecPlot3D(fout_tech, value, name_value, name_in_file);
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

		solver_grid.~Grid_forMech();
	}
}

void Solve_ElastodynamicsProblem()
{
	FEM::Grid_forMech solver_grid; //output
	std::vector<Point<double>> Solution; //output


	//char properties_file[1000] = { "E:/+cyl/800el/param.txt" };
	//char properties_file[1000] = { "E:/Box/100x10x200/BoxWithCracks/67k/param.txt" };
	//char properties_file[1000] = { "E:/Box/100x10x200/BoxWithCracks/638k/param.txt" };
	//char properties_file[1000] = { "E:/Box/100x50x200/param.txt" };
	//char properties_file[1000] = { "E:/Box/200x200x200/200k_fem/param_for_solver.txt" };
	//char properties_file[1000] = { "./param_for_solver.txt" };
	char properties_file[1000] = { "D:/Elastodynamic/homocyl/7k/param_for_solver.txt" };
	//char properties_file[1000] = { "base_properties.txt" };
	printf_s("The properties file contains the information:\n==========================================\n");
	printf_s("\tDo you need to print to a log file? 0-false/1-true\n");
	printf_s("\tinput: Mesh directory: box_200x200x200/90el\n");
	printf_s("\tproperties: Iterative process:\n");
	printf_s("\t\t<Start iteration> <End iteration>\n");
	printf_s("\t\t<Step size for Time>\n");
	printf_s("\t\t<Critical residual value for SLAE solution>\n");
	printf_s("\t<properties: Boundary values>\n");
	printf_s("\t<properties: Materials (E, v, rpho)>\n");
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
	math::ReadNonEmptyLine(f_properties, base_result_directory);
	math::SimpleGrid geo_grid; //input
	std::vector<std::vector<std::vector<int>>> boundary_faces;
	geo_grid.ReadFromNVTR(mesh_directory, 4);
	geo_grid.ReadFacesBoundaryNVTR(mesh_directory, boundary_faces);

	int start_iteration, end_iteration;
	{
		char _line[1000];
		math::ReadNonEmptyLine_forNumbers(f_properties, _line);
		std::vector<int> val;
		math::ParserStringToVectorInt(_line, val, " ");
		start_iteration = val[0];
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

			//
			/*first_boundaries[0].id_vertexes.resize(1);
			first_boundaries[0].id_vertexes[0] = 8;*/
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
		double _rpho;

		_material(double _E, double _v, double _rpho)
		{
			this->_E = _E;
			this->_v = _v;
			this->_rpho = _rpho;
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
		double _E, _v, _rpho;
		{
			char _line[1000];
			math::ReadNonEmptyLine_forNumbers(f_properties, _line);
			std::vector<float> val;
			math::ParserStringToVectorFloat(_line, val, " ");
			_E = val[0];
			_v = val[1];
			_rpho = val[2];
		}
		_materials.push_back(_material(_E, _v, _rpho));
	}

	fclose(f_properties);

	CreateDirectory((LPCTSTR)base_result_directory, NULL);

	bool res = false;
	double TIME_h = step_size;
	std::vector<Point<double>> U_prev(geo_grid.xyz.size()), U_prevprev(geo_grid.xyz.size());
	Solution.resize(geo_grid.xyz.size());
	math::InitializationVector(U_prev, 0);
	math::InitializationVector(U_prevprev, 0);
	math::InitializationVector(Solution, 1.0e-8);
	for (int id_STEP = start_iteration; id_STEP < end_iteration; id_STEP++)
	{
		printf_s("\n=============================*===============================\n");
		printf_s("================= Start solution of %d STEP (time = %.2lf) ================\n", id_STEP, TIME_h * id_STEP);
		double TIME_curr = TIME_h * id_STEP;
		//Create new directories
		char result_directory[1000];
		sprintf_s(result_directory, sizeof(result_directory), "%s/STEP_%d_t=%.2e", base_result_directory, id_STEP, TIME_curr);
		CreateDirectory((LPCTSTR)result_directory, NULL);

		for (int i = 0; i < _materials.size(); i++)
		{
			solver_grid.AddDomain();
			auto domain = solver_grid.GetDomain(i);
			domain->forThermal.rpho = _materials[i]._rpho;
			domain->forMech.SetE(_materials[i]._E);
			domain->forMech.SetV(_materials[i]._v);
		}

		if (id_STEP >= start_iteration)
		{
			{
				Point<double> power(0, 0, -1e+6);
				double tay_plus = 0.0001; //половина длины импульса
				double tay_0 = 0.0001; //середина импульса
				double degree = 2; //степень функции
				std::function<Point<double>(Point<double>)> new_sourse_value = [TIME_curr, power, tay_0, degree, tay_plus](Point<double> X) -> Point<double>
				{
					//double sourse = exp(-pow(TIME_curr - tay_0, degree) / pow(tay_plus, degree));
					//return power;
					//return Point<double>(power.x * sourse, power.y * sourse, power.z * sourse);
					double sourse = TIME_curr < 14e-7 ? 10e+6 : 0;
					return Point<double>(0 * sourse, 0 * sourse, -1* sourse);
				};
				if(second_boundary.size() > 0)
					second_boundary[0].value = new_sourse_value;

				std::function<Point<double>(Point<bool>&, int)> new_sourse_value_D = [TIME_curr](Point<bool>& is_take, int id)->Point<double>
				{
					//double sourse = exp(-pow(TIME_curr - tay_0, degree) / pow(tay_plus, degree));
					//return power;
					//return Point<double>(power.x * sourse, power.y * sourse, power.z * sourse);
					is_take.x = true;
					is_take.y = true;
					is_take.z = true;
					double sourse = TIME_curr < 14e-7 ? 1e-8 : 0;
					return Point<double>(0 * sourse, 0 * sourse, -1 * sourse);
				};
				if (first_boundaries.size() > 1)
					first_boundaries[1].value = new_sourse_value_D;
			}

			std::function<std::vector<std::vector<double>>(int, Point<double>)> StiffnessCoef = [&](int elem, Point<double> X)->std::vector<std::vector<double>>
			{
				std::vector<std::vector<double>> D = solver_grid.GetDomain(solver_grid.GetElement(elem)->GetIdDomain())->forMech.GetD(3);
				for (int i = 0; i < D.size(); i++)
					for (int j = 0; j < D[i].size(); j++)
						D[i][j] *= 1;
				return D;
			};
			std::function<double(int, Point<double>)> MassCoef = [&](int elem, Point<double> X)->double
			{
				//return 0;

				auto domain = (solver_grid.GetDomain(solver_grid.GetElement(elem)->GetIdDomain()));
				return (1) * domain->forThermal.rpho / TIME_h / TIME_h;
			};
			std::function<Point<double>(int, Point<double>)> F = [&](int elem, Point<double> X)->Point<double>
			{
				Point<double> U_prev_in_X = solver_grid.GetSolutionInPoint(elem, X, U_prev);
				Point<double> U_prevprev_in_X = solver_grid.GetSolutionInPoint(elem, X, U_prevprev);

				return ( U_prevprev_in_X *(-1) + U_prev_in_X * 2)
					* solver_grid.GetDomain(solver_grid.GetElement(elem)->GetIdDomain())->forThermal.rpho
					/ TIME_h / TIME_h;
			};

			FEM::FEM_forElasticDeformation(
				is_print_logFile,
				critical_residual,
				geo_grid, //input
				first_boundaries, //input
				second_boundary, //input
				StiffnessCoef,
				MassCoef,
				F,
				result_directory, //output
				solver_grid, //output
				Solution //output
			);

			math::MakeCopyVector_A_into_B(U_prev, U_prevprev);
			math::MakeCopyVector_A_into_B(Solution, U_prev);
		}
		else {
			std::function<std::vector<std::vector<double>>(int, Point<double>)> StiffnessCoef = [&](int elem, Point<double> X)->std::vector<std::vector<double>>
			{
				return solver_grid.GetDomain(solver_grid.GetElement(elem)->GetIdDomain())->forMech.GetD(3);
			};
			std::function<double(int, Point<double>)> MassCoef = [&](int elem, Point<double> X)->double
			{
				return 0;
			};
			std::function<Point<double>(int, Point<double>)> F = [&](int elem, Point<double> X)->Point<double>
			{
				return Point<double>(0,0,0);
			};

			FEM::FEM_forElasticDeformation(
				is_print_logFile,
				critical_residual,
				geo_grid, //input
				first_boundaries, //input
				second_boundary, //input
				StiffnessCoef,
				MassCoef,
				F,
				result_directory, //output
				solver_grid, //output
				Solution //output
			);
		}

		//output solution
		printf_s("Print the mech result into .dat file... ");
		//without deformations
		if (true) {
			FILE* fout_tech;
			char name_u_tech[5000];
			sprintf_s(name_u_tech, "%s/U_non_deformation.dat", result_directory);
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

				value[0][i] = sigma[0];
				value[1][i] = sigma[1];
				value[2][i] = sigma[2];

				value[3][i] = eps[0];
				value[4][i] = eps[1];
				value[5][i] = eps[2];

				value[3][i] = U.x;
				value[4][i] = U.y;
				value[5][i] = U.z;

				if (sigma_inv > sigma_inv_max)
				{
					sigma_inv_max = sigma_inv;
					elem_sigma_max = i;
				}
			}
			value[3].resize(solver_grid.GetDOFsCount());
			value[4].resize(solver_grid.GetDOFsCount());
			value[5].resize(solver_grid.GetDOFsCount());
			for (int i = 0; i < value[3].size(); i++)
			{
				value[3][i] = Solution[i].x;
				value[4][i] = Solution[i].y;
				value[5][i] = Solution[i].z;

			}
			//solver_grid.MoveCoordinates(Solution);

			solver_grid.printTecPlot3D_DiffDomains(fout_tech, value, name_value, name_in_file);
			fclose(fout_tech);
		}
		printf_s("\t complite\n");

		solver_grid.~Grid_forMech();
	}
}
void Solve_ElastodynamicsProblem_Explicit()
{
	FEM::Grid_forMech solver_grid; //output
	//FEM::Grid_forMech_Order2 solver_grid;


	//char properties_file[1000] = { "E:/+cyl/800el/param.txt" };
	//char properties_file[1000] = { "E:/Box/100x10x200/BoxWithCracks/67k/param.txt" };
	//char properties_file[1000] = { "E:/Box/100x10x200/BoxWithCracks/638k/param.txt" };
	//char properties_file[1000] = { "E:/Box/100x50x200/param.txt" };
	//char properties_file[1000] = { "E:/Box/200x200x200/200k_fem/param_for_solver.txt" };
	//char properties_file[1000] = { "./param_for_solver.txt" };
	char properties_file[1000] = { "D:/Elastodynamic/homocyl/73k/param_for_solver.txt" };
	//char properties_file[1000] = { "base_properties.txt" };
	printf_s("The properties file contains the information:\n==========================================\n");
	printf_s("\tDo you need to print to a log file? 0-false/1-true\n");
	printf_s("\tinput: Mesh directory: box_200x200x200/90el\n");
	printf_s("\tproperties: Iterative process:\n");
	printf_s("\t\t<Start iteration> <End iteration>\n");
	printf_s("\t\t<Step size for Time>\n");
	printf_s("\t\t<Critical residual value for SLAE solution>\n");
	printf_s("\t<properties: Boundary values>\n");
	printf_s("\t<properties: Materials (E, v, rpho)>\n");
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
	math::ReadNonEmptyLine(f_properties, base_result_directory);
	math::SimpleGrid geo_grid; //input
	std::vector<std::vector<std::vector<int>>> boundary_faces;
	geo_grid.ReadFromNVTR(mesh_directory, 4);
	geo_grid.ReadFacesBoundaryNVTR(mesh_directory, boundary_faces);

	int start_iteration, end_iteration;
	{
		char _line[1000];
		math::ReadNonEmptyLine_forNumbers(f_properties, _line);
		std::vector<int> val;
		math::ParserStringToVectorInt(_line, val, " ");
		start_iteration = val[0];
		end_iteration = val[1];
	}
	double step_size;
	int step_for_out;
	{
		char _line[1000];
		math::ReadNonEmptyLine_forNumbers(f_properties, _line);
		std::vector<float> val;
		math::ParserStringToVectorFloat(_line, val, " ");
		step_size = val[0];
		std::vector<int> val2;
		math::ReadNonEmptyLine_forNumbers(f_properties, _line);
		math::ParserStringToVectorInt(_line, val2, " ");
		step_for_out = (int)val2[0];
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

			//
			/*first_boundaries[0].id_vertexes.resize(1);
			first_boundaries[0].id_vertexes[0] = 8;*/
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
		double _rpho;

		_material(double _E, double _v, double _rpho)
		{
			this->_E = _E;
			this->_v = _v;
			this->_rpho = _rpho;
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
		double _E, _v, _rpho;
		{
			char _line[1000];
			math::ReadNonEmptyLine_forNumbers(f_properties, _line);
			std::vector<float> val;
			math::ParserStringToVectorFloat(_line, val, " ");
			_E = val[0];
			_v = val[1];
			_rpho = val[2];
		}
		_materials.push_back(_material(_E, _v, _rpho));
	}

	fclose(f_properties);

	CreateDirectory((LPCTSTR)base_result_directory, NULL);

	//Create base matrixes
	CSSD_Matrix<Tensor2Rank3D, Point<double>> MassMatrix, MassMatrix_base, StiffnessMatrix, StiffnessMatrix_base;
	{
		printf_s("================= Create base matrix ================\n");
		printf("Initialization of grid...\t");
		solver_grid.Initialization(geo_grid, first_boundaries, second_boundary);
		for (int i = 0; i < _materials.size(); i++)
		{
			solver_grid.AddDomain();
			auto domain = solver_grid.GetDomain(i);
			domain->forThermal.rpho = _materials[i]._rpho;
			domain->forMech.SetE(_materials[i]._E);
			domain->forMech.SetV(_materials[i]._v);
		}
		printf_s("complite\n");
		printf("Creation the SLAE portrait...\t");
		solver_grid.CreationPortrait(MassMatrix_base);
		StiffnessMatrix_base.Initialization(MassMatrix_base.id_column_for_A_up, MassMatrix_base.id_column_for_A_down);
		printf_s("\t\tcomplite\n");

		std::function<std::vector<std::vector<double>>(int, Point<double>)> StiffnessCoef = [&](int elem, Point<double> X)->std::vector<std::vector<double>>
		{
			std::vector<std::vector<double>> D = solver_grid.GetDomain(solver_grid.GetElement(elem)->GetIdDomain())->forMech.GetD(3);
			/*for (int i = 0; i < D.size(); i++)
				for (int j = 0; j < D[i].size(); j++)
					D[i][j] *= 1;*/
			return D;
		};
		std::function<double(int, Point<double>)> MassCoef = [&](int elem, Point<double> X)->double
		{
			auto domain = (solver_grid.GetDomain(solver_grid.GetElement(elem)->GetIdDomain()));
			return domain->forThermal.rpho;
		};
		std::function<Point<double>(int, Point<double>)> VolumeForсe = [&](int elem, Point<double> X)->Point<double>
		{
			return Point<double>(0, 0, 0);
		};
		//SLAE assembling
		{
			printf("STIFFNESS assembling...\n");
			std::vector<DenseMatrix<Tensor2Rank3D, Point<double>>> local_SLAE_stiffness(solver_grid.GetElementsCount());
			omp_set_num_threads(math::NUM_THREADS);
#pragma omp parallel for
			for (int id_elem = 0; id_elem < solver_grid.GetElementsCount(); id_elem++)
			{
				if (id_elem % 1000 == 0)
					printf("Solve element[%d]\r", id_elem);
				auto element = solver_grid.GetElement(id_elem);

				std::function<std::vector<std::vector<double>>(Point<double>)> D = [&](Point<double> X) {return StiffnessCoef(id_elem, X); };
				std::function<double(Point<double>)> M = [&](Point<double> X) {return 0/*MassCoef(id_elem, X)*/; };
				std::function<Point<double>(Point<double>)> F = [&](Point<double> X) {return VolumeForсe(id_elem, X); };

				element->SolveLocalMatrix(local_SLAE_stiffness[id_elem], D, M, F);
			}
			for (int id_elem = 0; id_elem < solver_grid.GetElementsCount(); id_elem++)
			{
				if (id_elem % 1000 == 0)
					printf("Add the local matrix of element[%d]\r", id_elem);
				StiffnessMatrix_base.SummPartOfMatrix(local_SLAE_stiffness[id_elem], *solver_grid.GetElementDOFs(id_elem));
			}
			printf_s("                                                                                    \r");
			printf_s("\t\tcomplite\n");

			printf("MASS assembling...\n");
			std::vector<DenseMatrix<Tensor2Rank3D, Point<double>>> local_SLAE_mass(solver_grid.GetElementsCount());
			omp_set_num_threads(math::NUM_THREADS);
#pragma omp parallel for
			for (int id_elem = 0; id_elem < solver_grid.GetElementsCount(); id_elem++)
			{
				if (id_elem % 1000 == 0)
					printf("Solve element[%d]\r", id_elem);
				auto element = solver_grid.GetElement(id_elem);

				std::function<std::vector<std::vector<double>>(Point<double>)> D = [&](Point<double> X) 
				{
					std::vector <std::vector<double>> dd(6);
					for (int i = 0; i < dd.size(); i++)
					{
						dd[i].resize(6);
						math::InitializationVector(dd[i], 0);
					}
					return dd; 
				};
				std::function<double(Point<double>)> M = [&](Point<double> X) {return MassCoef(id_elem, X); };
				std::function<Point<double>(Point<double>)> F = [&](Point<double> X) {return Point<double>(0,0,0); };

				element->SolveLocalMatrix(local_SLAE_mass[id_elem], D, M, F);
			}
			for (int id_elem = 0; id_elem < solver_grid.GetElementsCount(); id_elem++)
			{
				if (id_elem % 1000 == 0)
					printf("Add the local matrix of element[%d]\r", id_elem);
				MassMatrix_base.SummPartOfMatrix(local_SLAE_mass[id_elem], *solver_grid.GetElementDOFs(id_elem));
			}
			printf_s("                                                                                    \r");
			printf_s("\t\tcomplite\n");
		}
	}

	bool res = false;
	double TIME_h = step_size;

	std::vector<Point<double>> U_prev(solver_grid.GetDOFsCount()), U_prevprev(solver_grid.GetDOFsCount()), U_curr(solver_grid.GetDOFsCount());
	math::InitializationVector(U_curr, 0);
	math::InitializationVector(U_prev, 0);
	math::InitializationVector(U_prevprev, 0);
	
	for (int id_STEP = start_iteration; id_STEP < end_iteration; id_STEP++)
	{
		printf_s("\n================= Start solution of %d STEP (time = %.2lf) ================\n", id_STEP, TIME_h * id_STEP);
		double TIME_curr = TIME_h * id_STEP;
		bool is_print_result = false;
		char result_directory[1000];
		if (id_STEP % step_for_out == 0)
		{
			is_print_result = true;
			sprintf_s(result_directory, sizeof(result_directory), "%s/STEP_%d_t=%.2e", base_result_directory, id_STEP, TIME_curr);
			CreateDirectory((LPCTSTR)result_directory, NULL);
		}
		

		StiffnessMatrix.SetMatrix(StiffnessMatrix_base);
		MassMatrix.SetMatrix(MassMatrix_base);

		{
			Point<double> power(0, 0, -1e+6);
			double tay_plus = 0.0001; //половина длины импульса
			double tay_0 = 0.0001; //середина импульса
			double degree = 2; //степень функции
			std::function<Point<double>(Point<double>)> new_sourse_value = [TIME_curr, power, tay_0, degree, tay_plus](Point<double> X) -> Point<double>
			{
				//double sourse = exp(-pow(TIME_curr - tay_0, degree) / pow(tay_plus, degree));
				//return power;
				//return Point<double>(power.x * sourse, power.y * sourse, power.z * sourse);
				double sourse = TIME_curr < 14e-7 ? 10e+6 : 0;
				return Point<double>(0 * sourse, 0 * sourse, -1 * sourse);
			};
			if (second_boundary.size() > 0)
			{
				second_boundary[0].value = new_sourse_value;

				for (int i = 0; i < solver_grid.boundary_faces.size(); i++)
				{
					solver_grid.boundary_faces[i].boundary_value = second_boundary[solver_grid.boundary_faces[i].id_type].value;
				}
			}

			std::function<Point<double>(Point<bool>&, int)> new_sourse_value_D = [TIME_curr](Point<bool>& is_take, int id)->Point<double>
			{
				//double sourse = exp(-pow(TIME_curr - tay_0, degree) / pow(tay_plus, degree));
				//return power;
				//return Point<double>(power.x * sourse, power.y * sourse, power.z * sourse);
				is_take.x = true;
				is_take.y = true;
				is_take.z = true;
				double sourse = TIME_curr < 14e-7 ? 1e-8 : 0;
				return Point<double>(0 * sourse, 0 * sourse, -1 * sourse);
			};
			if (false && first_boundaries.size() > 1)
			{
				first_boundaries[1].value = new_sourse_value_D;

				for (int i = 0; i < solver_grid.boundary_vertexes.size(); i++)
				{
					if (solver_grid.boundary_vertexes[i].id_type == 1)
					{
						solver_grid.boundary_vertexes[i].boundary_value = first_boundaries[1].value;
					}
				}
			}
		}

		//Boundary condition
		if (false)
		{
			//first condition
			printf("First boundary conditions...");
			//Stiffness
			for (int id_vertex = 0; id_vertex < solver_grid.boundary_vertexes.size(); id_vertex++)
			{
				Tensor2Rank3D diag, non_diag;
				Point<double> X, F;
				diag.InitializationAsI();
				non_diag.InitializationAs0();

				auto boundary = &(solver_grid.boundary_vertexes[id_vertex]);

				Point<bool> is_take;
				Point<double> boundary_value = boundary->boundary_value(is_take, id_vertex);
				///-->>>>>
				/*boundary_value.x = 0;
				boundary_value.y = 0;
				boundary_value.z = 20*(*solver_grid.GetPtrCoordinateViaID(id_vertex)).z - 20*200*0.5;*/
				///-->>>>>
				int global_id = boundary->GetDOFInLocalID(0);

				auto enter_boundary = [&](int position) {
					switch (position)
					{
					case 0:
						StiffnessMatrix.X[global_id].x = boundary_value.x;
						StiffnessMatrix.F[global_id].x = boundary_value.x;
						if (StiffnessMatrix.F[global_id].x != StiffnessMatrix.F[global_id].x)
						{
							printf_s("F[%d] = %.2e\n", global_id, StiffnessMatrix.F[global_id].x);
						}
						break;
					case 1:
						StiffnessMatrix.X[global_id].y = boundary_value.y;
						StiffnessMatrix.F[global_id].y = boundary_value.y;
						if (StiffnessMatrix.F[global_id].y != StiffnessMatrix.F[global_id].y)
						{
							printf_s("F[%d] = %.2e\n", global_id, StiffnessMatrix.F[global_id].y);
						}
						break;
					case 2:
						StiffnessMatrix.X[global_id].z = boundary_value.z;
						StiffnessMatrix.F[global_id].z = boundary_value.z;
						if (StiffnessMatrix.F[global_id].z != StiffnessMatrix.F[global_id].z)
						{
							printf_s("F[%d] = %.2e\n", global_id, StiffnessMatrix.F[global_id].z);
						}
						break;
					default:
						break;
					}

					StiffnessMatrix.Diag[global_id].val[position][0] = 0.0;
					StiffnessMatrix.Diag[global_id].val[position][1] = 0.0;
					StiffnessMatrix.Diag[global_id].val[position][2] = 0.0;
					StiffnessMatrix.Diag[global_id].val[position][position] = 1.0;
					for (int j = 0; j < StiffnessMatrix.A_down[global_id].size(); j++)
					{
						StiffnessMatrix.A_down[global_id][j].val[position][0] = 0.0;
						StiffnessMatrix.A_down[global_id][j].val[position][1] = 0.0;
						StiffnessMatrix.A_down[global_id][j].val[position][2] = 0.0;
					}
					for (int j = 0; j < StiffnessMatrix.A_up[global_id].size(); j++)
					{
						StiffnessMatrix.A_up[global_id][j].val[position][0] = 0.0;
						StiffnessMatrix.A_up[global_id][j].val[position][1] = 0.0;
						StiffnessMatrix.A_up[global_id][j].val[position][2] = 0.0;
					}
				};
				if (is_take.x == true) enter_boundary(0);
				if (is_take.y == true) enter_boundary(1);
				if (is_take.z == true) enter_boundary(2);
			}
			//Mass
			for (int id_vertex = 0; id_vertex < solver_grid.boundary_vertexes.size(); id_vertex++)
			{
				Tensor2Rank3D diag, non_diag;
				Point<double> X, F;
				diag.InitializationAsI();
				non_diag.InitializationAs0();

				auto boundary = &(solver_grid.boundary_vertexes[id_vertex]);

				Point<bool> is_take;
				Point<double> boundary_value = boundary->boundary_value(is_take, id_vertex);
				///-->>>>>
				/*boundary_value.x = 0;
				boundary_value.y = 0;
				boundary_value.z = 20*(*solver_grid.GetPtrCoordinateViaID(id_vertex)).z - 20*200*0.5;*/
				///-->>>>>
				int global_id = boundary->GetDOFInLocalID(0);

				auto enter_boundary = [&](int position) {
					switch (position)
					{
					case 0:
						MassMatrix.X[global_id].x = boundary_value.x;
						MassMatrix.F[global_id].x = boundary_value.x;
						if (MassMatrix.F[global_id].x != MassMatrix.F[global_id].x)
						{
							printf_s("F[%d] = %.2e\n", global_id, MassMatrix.F[global_id].x);
						}
						break;
					case 1:
						MassMatrix.X[global_id].y = boundary_value.y;
						MassMatrix.F[global_id].y = boundary_value.y;
						if (MassMatrix.F[global_id].y != MassMatrix.F[global_id].y)
						{
							printf_s("F[%d] = %.2e\n", global_id, MassMatrix.F[global_id].y);
						}
						break;
					case 2:
						MassMatrix.X[global_id].z = boundary_value.z;
						MassMatrix.F[global_id].z = boundary_value.z;
						if (MassMatrix.F[global_id].z != MassMatrix.F[global_id].z)
						{
							printf_s("F[%d] = %.2e\n", global_id, MassMatrix.F[global_id].z);
						}
						break;
					default:
						break;
					}

					MassMatrix.Diag[global_id].val[position][0] = 0.0;
					MassMatrix.Diag[global_id].val[position][1] = 0.0;
					MassMatrix.Diag[global_id].val[position][2] = 0.0;
					MassMatrix.Diag[global_id].val[position][position] = 1.0;
					for (int j = 0; j < MassMatrix.A_down[global_id].size(); j++)
					{
						MassMatrix.A_down[global_id][j].val[position][0] = 0.0;
						MassMatrix.A_down[global_id][j].val[position][1] = 0.0;
						MassMatrix.A_down[global_id][j].val[position][2] = 0.0;
					}
					for (int j = 0; j < MassMatrix.A_up[global_id].size(); j++)
					{
						MassMatrix.A_up[global_id][j].val[position][0] = 0.0;
						MassMatrix.A_up[global_id][j].val[position][1] = 0.0;
						MassMatrix.A_up[global_id][j].val[position][2] = 0.0;
					}
				};
				if (is_take.x == true) enter_boundary(0);
				if (is_take.y == true) enter_boundary(1);
				if (is_take.z == true) enter_boundary(2);
			}

		}
		//second condition
		for (int id_triangle = 0; id_triangle < solver_grid.boundary_faces.size(); id_triangle++)
		{
			std::vector<Point<double>> local_vector_SLAE;
			solver_grid.boundary_faces[id_triangle].SolveLocalBoundaryVector(local_vector_SLAE, solver_grid.boundary_faces[id_triangle].boundary_value);
			StiffnessMatrix.SummPartOfVector(local_vector_SLAE, *solver_grid.boundary_faces[id_triangle].GetElementDOFs());
		};

		//solve (F-Stiffness*U_prev)
		StiffnessMatrix.MultiplicationMatrixVector(U_prev, MassMatrix.F);
		math::MakeSummVectors(StiffnessMatrix.F, Point<double>(1.0, 1.0, 1.0), MassMatrix.F, Point<double>(-1.0, -1.0, -1.0), MassMatrix.F);

		//костыль для "обращения" матрицы
		//результат (Mass^-1)*(F-Stiffness*U_prev) в MassMatrix.X
		{
			CSSD_Matrix<double, double> MassMatrix_doubleSLAE;
			math::MakeCopyMatrix_A_into_B(MassMatrix, MassMatrix_doubleSLAE);

			if (id_STEP == 0 && false)
			{
				for (int i = 0; i < MassMatrix_doubleSLAE.GetMatrixSize() / 3; i++)
				{
					MassMatrix_doubleSLAE.X[i * 3 + 0] = 1e-10;
					MassMatrix_doubleSLAE.X[i * 3 + 1] = 1e-10;
					MassMatrix_doubleSLAE.X[i * 3 + 2] = 1e-10;
				}
			}
			else {
				for (int i = 0; i < MassMatrix_doubleSLAE.GetMatrixSize() / 3; i++)
				{
					MassMatrix_doubleSLAE.X[i * 3 + 0] = U_prev[i].x;
					MassMatrix_doubleSLAE.X[i * 3 + 1] = U_prev[i].y;
					MassMatrix_doubleSLAE.X[i * 3 + 2] = U_prev[i].z;
				}
			}
						
			printf("Soluting SLAY... (%d)\n", MassMatrix_doubleSLAE.GetMatrixSize());
			int MaxSize = MassMatrix_doubleSLAE.GetMatrixSize();
			MaxSize = MaxSize / 10 < 10 ? 10 : MaxSize / 10;
			std::vector<double> best_solution;
			math::MakeCopyVector_A_into_B(MassMatrix_doubleSLAE.X, best_solution);
			double current_residual = 1, best_residual = 1e+25;
			int MAX_STEPS = 15;
			int ii = 1;
			CSSD_Matrix<double, double> Predcondor;
			//Predcondor.PrecondorSSOR_summetric(0.75, MassMatrix_doubleSLAE);
			for (int i = 0; i <= MAX_STEPS; i++)
			{
				double needed_residual = pow(10., -1 * (ii + 1));
				//current_residual /= 2.;
				printf_s("//---> I = %d/%d (full size %d) - needed residual %.2e\n", i, MAX_STEPS, MassMatrix_doubleSLAE.GetMatrixSize(), needed_residual);
				

				//current_residual = abs(newSLAE.BiCG_Stab(MaxSize, needed_residual));
				current_residual = abs(MassMatrix_doubleSLAE.MSG(MaxSize, needed_residual));
				//current_residual = abs(MassMatrix_doubleSLAE.MSG_Preconditioning(MaxSize, needed_residual, Predcondor));
				//current_residual = abs(newSLAE.BiCG_Stab_Preconditioning(MaxSize, needed_residual, Predcondor));

				if (current_residual < needed_residual)
				{
					i = 0;
					ii++;
					if (current_residual > best_residual + (1e-10))
					{
						//math::MakeCopyVector_A_into_B(best_solution, newSLAE.X);
						//printf_s("//---> BEST residual %.2e\n", best_residual);
						break;
					}
				}
				if (current_residual <= critical_residual)
				{
					best_residual = current_residual;
					break;
				}
				if (current_residual < best_residual)
				{
					best_residual = current_residual;
					math::MakeCopyVector_A_into_B(MassMatrix_doubleSLAE.X, best_solution);
				}
			}
			math::MakeCopyVector_A_into_B(best_solution, MassMatrix_doubleSLAE.X);
			printf_s("//---> BEST residual %.2e\n", best_residual);
			math::MakeCopyVector_A_into_B(MassMatrix_doubleSLAE.X, MassMatrix.X);

			if (best_residual >= 1)
			{
				printf_s("There are problems in the solution\n");
				int a;
				scanf_s("%d", &a);
			}
		}

		//solve U_curr = ht^2 * (Mass^-1)*(F-Stiffness*U_prev) + 2 * U_prev - U_prevprev
		for (int i = 0; i < MassMatrix.X.size(); i++)
		{
			U_curr[i] = MassMatrix.X[i] * TIME_h * TIME_h + U_prev[i] * 2 - U_prevprev[i];
		}
		//поправка на первые краевые
		for (int id_vertex = 0; id_vertex < solver_grid.boundary_vertexes.size(); id_vertex++)
		{
			auto boundary = &(solver_grid.boundary_vertexes[id_vertex]);

			Point<bool> is_take;
			Point<double> boundary_value = boundary->boundary_value(is_take, id_vertex);
			if (boundary_value.z > 0)
			{
				printf("%d\n", id_vertex);
			}
			///-->>>>>
			/*boundary_value.x = 0;
			boundary_value.y = 0;
			boundary_value.z = 20*(*solver_grid.GetPtrCoordinateViaID(id_vertex)).z - 20*200*0.5;*/
			///-->>>>>
			int global_id = boundary->GetDOFInLocalID(0);

			auto enter_boundary = [&](int position) {
				switch (position)
				{
				case 0:
					U_curr[global_id].x = boundary_value.x;
					break;
				case 1:
					U_curr[global_id].y = boundary_value.y;
					break;
				case 2:
					U_curr[global_id].z = boundary_value.z;
				default:
					break;
				}
			};
			if (is_take.x == true) enter_boundary(0);
			if (is_take.y == true) enter_boundary(1);
			if (is_take.z == true) enter_boundary(2);
		}

		math::MakeCopyVector_A_into_B(U_prev, U_prevprev);
		math::MakeCopyVector_A_into_B(U_curr, U_prev);

		//output solution
		printf_s("Print the mech result into .dat file... ");
		//without deformations
		if (is_print_result) {
			FILE* fout_tech;
			char name_u_tech[5000];
			sprintf_s(name_u_tech, "%s/U_non_deformation.dat", result_directory);
			fopen_s(&fout_tech, name_u_tech, "w");
			char name_in_file[1000];
			sprintf_s(name_in_file, "Elastic");
			std::vector<std::vector<char>> name_value(6);
			char name_v_tmp[6][100];
			sprintf_s(name_v_tmp[0], "Mises_stress");
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

				Point<double> U = solver_grid.GetSolutionInPoint(i, Centr, U_curr);
				Point<Point<double>> dU = solver_grid.GetDerevativeFromSolutionInPoint(i, Centr, U_curr);

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

				value[0][i] = sigma_inv;
				value[1][i] = sigma[1];
				value[2][i] = sigma[2];

				value[3][i] = eps[0];
				value[4][i] = eps[1];
				value[5][i] = eps[2];

				value[3][i] = U.x;
				value[4][i] = U.y;
				value[5][i] = U.z;

				if (sigma_inv > sigma_inv_max)
				{
					sigma_inv_max = sigma_inv;
					elem_sigma_max = i;
				}
			}
			value[3].resize(solver_grid.GetDOFsCount());
			value[4].resize(solver_grid.GetDOFsCount());
			value[5].resize(solver_grid.GetDOFsCount());
			for (int i = 0; i < value[3].size(); i++)
			{
				value[3][i] = U_curr[i].x;
				value[4][i] = U_curr[i].y;
				value[5][i] = U_curr[i].z;

			}
			//solver_grid.MoveCoordinates(Solution);
			for (int i = 0; i < solver_grid.GetVertexCount(); i++)
			{
				solver_grid.MoveTheVertex(i, U_curr[i]);
			}

			solver_grid.printTecPlot3D_DiffDomains(fout_tech, value, name_value, name_in_file);

			for (int i = 0; i < solver_grid.GetVertexCount(); i++)
			{
				solver_grid.MoveTheVertex(i, U_curr[i] * (-1) );
			}

			fclose(fout_tech);
		}
		printf_s("\t complite\n");

		//solver_grid.~Grid_forMech();
	}
}
void Solve_ElastodynamicsProblem_Explicit_v2()
{
	FEM::Grid_forMech solver_grid; //output
	std::vector<Point<double>> Solution; //output

	bool is_STATIONARY = true;


	//char properties_file[1000] = { "E:/+cyl/800el/param.txt" };
	//char properties_file[1000] = { "E:/Box/100x10x200/BoxWithCracks/67k/param.txt" };
	//char properties_file[1000] = { "E:/Box/100x10x200/BoxWithCracks/638k/param.txt" };
	//char properties_file[1000] = { "E:/Box/100x50x200/param.txt" };
	//char properties_file[1000] = { "E:/Box/200x200x200/200k_fem/param_for_solver.txt" };
	//char properties_file[1000] = { "./param_for_solver.txt" };
	char properties_file[1000] = { "D:/Elastodynamic/homocyl/67k/param_for_solver.txt" };
	//char properties_file[1000] = { "base_properties.txt" };
	printf_s("The properties file contains the information:\n==========================================\n");
	printf_s("\tDo you need to print to a log file? 0-false/1-true\n");
	printf_s("\tinput: Mesh directory: box_200x200x200/90el\n");
	printf_s("\tproperties: Iterative process:\n");
	printf_s("\t\t<Start iteration> <End iteration>\n");
	printf_s("\t\t<Step size for Time>\n");
	printf_s("\t\t<Critical residual value for SLAE solution>\n");
	printf_s("\t<properties: Boundary values>\n");
	printf_s("\t<properties: Materials (E, v, rpho)>\n");
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
	math::ReadNonEmptyLine(f_properties, base_result_directory);
	math::SimpleGrid geo_grid; //input
	std::vector<std::vector<std::vector<int>>> boundary_faces;
	geo_grid.ReadFromNVTR(mesh_directory, 4);
	geo_grid.ReadFacesBoundaryNVTR(mesh_directory, boundary_faces);

	int start_iteration, end_iteration;
	{
		char _line[1000];
		math::ReadNonEmptyLine_forNumbers(f_properties, _line);
		std::vector<int> val;
		math::ParserStringToVectorInt(_line, val, " ");
		start_iteration = val[0];
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

			//
			/*first_boundaries[0].id_vertexes.resize(1);
			first_boundaries[0].id_vertexes[0] = 8;*/
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
		double _rpho;

		_material(double _E, double _v, double _rpho)
		{
			this->_E = _E;
			this->_v = _v;
			this->_rpho = _rpho;
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
		double _E, _v, _rpho;
		{
			char _line[1000];
			math::ReadNonEmptyLine_forNumbers(f_properties, _line);
			std::vector<float> val;
			math::ParserStringToVectorFloat(_line, val, " ");
			_E = val[0];
			_v = val[1];
			_rpho = val[2];
		}
		_materials.push_back(_material(_E, _v, _rpho));
	}

	fclose(f_properties);

	CreateDirectory((LPCTSTR)base_result_directory, NULL);
	bool res = false;
	int step_coef = 1;
	double TIME_h = step_size / step_coef;
	if (is_STATIONARY) {
		end_iteration = 1;
		start_iteration = 0;
	}
	std::vector<Point<double>> U_prev(geo_grid.xyz.size()), U_prevprev(geo_grid.xyz.size()), U_curr(geo_grid.xyz.size());
	math::InitializationVector(U_curr, 0);
	math::InitializationVector(U_prev, 0);
	math::InitializationVector(U_prevprev, 0);

	//Create base matrixes
	CSSD_Matrix<Tensor2Rank3D, Point<double>> StiffnessMatrix, StiffnessMatrix_base;
	{
		printf_s("================= Create base matrix ================\n");
		printf("Initialization of grid...\t");
		solver_grid.Initialization(geo_grid, first_boundaries, second_boundary);
		for (int i = 0; i < _materials.size(); i++)
		{
			solver_grid.AddDomain();
			auto domain = solver_grid.GetDomain(i);
			domain->forThermal.rpho = _materials[i]._rpho;
			domain->forMech.SetE(_materials[i]._E);
			domain->forMech.SetV(_materials[i]._v);
		}
		printf_s("complite\n");
		printf("Creation the SLAE portrait...\t");
		solver_grid.CreationPortrait(StiffnessMatrix_base);
		StiffnessMatrix.Initialization(StiffnessMatrix_base.id_column_for_A_up, StiffnessMatrix_base.id_column_for_A_down);
		printf_s("\t\tcomplite\n");

		std::function<std::vector<std::vector<double>>(int, Point<double>)> StiffnessCoef = [&](int elem, Point<double> X)->std::vector<std::vector<double>>
		{
			std::vector<std::vector<double>> D = solver_grid.GetDomain(solver_grid.GetElement(elem)->GetIdDomain())->forMech.GetD(3);
			/*for (int i = 0; i < D.size(); i++)
				for (int j = 0; j < D[i].size(); j++)
					D[i][j] *= 1;*/
			return D;
		};
		std::function<double(int, Point<double>)> MassCoef = [&](int elem, Point<double> X)->double
		{
			if (is_STATIONARY)
			{
				return 0;
			}
			else {
				auto domain = (solver_grid.GetDomain(solver_grid.GetElement(elem)->GetIdDomain()));
				return domain->forThermal.rpho / (TIME_h * TIME_h);
			}
		};
		std::function<Point<double>(int, Point<double>)> VolumeForсe = [&](int elem, Point<double> X)->Point<double>
		{
			return Point<double>(0, 0, 0);
		};
		//SLAE assembling
		{
			printf("Matrix assembling...\n");
			std::vector<DenseMatrix<Tensor2Rank3D, Point<double>>> local_SLAE_stiffness(solver_grid.GetElementsCount());
			for (int id_elem = 0; id_elem < solver_grid.GetElementsCount(); id_elem++)
			{
				if (id_elem % 1000 == 0)
					printf("Solve element[%d]\r", id_elem);
				auto element = solver_grid.GetElement(id_elem);

				std::function<std::vector<std::vector<double>>(Point<double>)> D = [&](Point<double> X) {return StiffnessCoef(id_elem, X); };
				std::function<double(Point<double>)> M = [&](Point<double> X) {return MassCoef(id_elem, X); };
				std::function<Point<double>(Point<double>)> F = [&](Point<double> X) {return VolumeForсe(id_elem, X); };

				element->SolveLocalMatrix(local_SLAE_stiffness[id_elem], D, M, F);
			}
			for (int id_elem = 0; id_elem < solver_grid.GetElementsCount(); id_elem++)
			{
				if (id_elem % 1000 == 0)
					printf("Add the local matrix of element[%d]\r", id_elem);
				StiffnessMatrix_base.SummPartOfMatrix(local_SLAE_stiffness[id_elem], *solver_grid.GetElementDOFs(id_elem));
			}
			printf_s("                                                                                    \r");
			printf_s("\t\tcomplite\n");
		}
	}
	
	for (int id_STEP = start_iteration; id_STEP < end_iteration * step_coef; id_STEP++)
	{
		printf_s("\n================= Start solution of %d STEP (time = %.2lf) ================\n", id_STEP, TIME_h * id_STEP);
		double TIME_curr = TIME_h * id_STEP;
		bool is_print_result = false;
		char result_directory[1000];
		if (id_STEP % step_coef == 0)
		{
			is_print_result = true;
			sprintf_s(result_directory, sizeof(result_directory), "%s/STEP_%d_t=%.2e", base_result_directory, id_STEP, TIME_curr);
			CreateDirectory((LPCTSTR)result_directory, NULL);
		}

		StiffnessMatrix.SetMatrix(StiffnessMatrix_base);

		//переопределяем источник
		if(true) {
			Point<double> power(0, 0, -1e+6);
			double tay_plus = 0.0001; //половина длины импульса
			double tay_0 = 0.0001; //середина импульса
			double degree = 2; //степень функции
			std::function<Point<double>(Point<double>)> new_sourse_value = [TIME_curr, power, tay_0, degree, tay_plus](Point<double> X) -> Point<double>
			{
				//double sourse = exp(-pow(TIME_curr - tay_0, degree) / pow(tay_plus, degree));
				//return power;
				//return Point<double>(power.x * sourse, power.y * sourse, power.z * sourse);
				double sourse = TIME_curr < 14e-7 ? 10e+6 : 0;
				return Point<double>(0 * sourse, 0 * sourse, -1 * sourse);
			};
			if (second_boundary.size() > 0)
			{
				second_boundary[0].value = new_sourse_value;

				for (int i = 0; i < solver_grid.boundary_faces.size(); i++)
				{
					solver_grid.boundary_faces[i].boundary_value = second_boundary[solver_grid.boundary_faces[i].id_type].value;
				}
			}

			std::function<Point<double>(Point<bool>&, int)> new_sourse_value_D = [TIME_curr](Point<bool>& is_take, int id)->Point<double>
			{
				//double sourse = exp(-pow(TIME_curr - tay_0, degree) / pow(tay_plus, degree));
				//return power;
				//return Point<double>(power.x * sourse, power.y * sourse, power.z * sourse);
				is_take.x = true;
				is_take.y = true;
				is_take.z = true;
				double sourse = TIME_curr < 14e-7 ? 1e-8 : 0;
				return Point<double>(0 * sourse, 0 * sourse, -1 * sourse);
			};
			if (true && first_boundaries.size() > 1)
			{
				first_boundaries[1].value = new_sourse_value_D;

				for (int i = 0; i < solver_grid.boundary_vertexes.size(); i++)
				{
					if (solver_grid.boundary_vertexes[i].id_type == 1)
					{
						solver_grid.boundary_vertexes[i].boundary_value = first_boundaries[1].value;
					}
				}
			}
		}

		//second boyndary condition
		for (int id_triangle = 0; id_triangle < solver_grid.boundary_faces.size(); id_triangle++)
		{
			std::vector<Point<double>> local_vector_SLAE;
			solver_grid.boundary_faces[id_triangle].SolveLocalBoundaryVector(local_vector_SLAE, solver_grid.boundary_faces[id_triangle].boundary_value);
			StiffnessMatrix.SummPartOfVector(local_vector_SLAE, *solver_grid.boundary_faces[id_triangle].GetElementDOFs());
		};

		//additional right side
		if (true) {
			printf("NEW right side assembling...\n");
			std::vector<std::vector<Point<double>>> local_SLAE_vector(solver_grid.GetElementsCount());
			for (int id_elem = 0; id_elem < solver_grid.GetElementsCount(); id_elem++)
			{
				if (id_elem % 1000 == 0)
					printf("Solve element[%d]\r", id_elem);
				auto element = solver_grid.GetElement(id_elem);

				std::function<Point<double>(Point<double>)> F = [&](Point<double> X) -> Point<double>
				{
					if (is_STATIONARY)
					{
						return Point<double>(0, 0, 0);
					}
					else {
						Point<double> U_prev_in_X = solver_grid.GetSolutionInPoint(id_elem, X, U_prev);
						Point<double> U_prevprev_in_X = solver_grid.GetSolutionInPoint(id_elem, X, U_prevprev);

						return (U_prevprev_in_X * (-1) + U_prev_in_X * 2)
							* solver_grid.GetDomain(solver_grid.GetElement(id_elem)->GetIdDomain())->forThermal.rpho
							/ TIME_h / TIME_h;
					}
				};

				element->SolveRightSide(local_SLAE_vector[id_elem], F);
			}
			for (int id_elem = 0; id_elem < solver_grid.GetElementsCount(); id_elem++)
			{
				if (id_elem % 1000 == 0)
					printf("Add the local matrix of element[%d]\r", id_elem);
				StiffnessMatrix.SummPartOfVector(local_SLAE_vector[id_elem], *solver_grid.GetElementDOFs(id_elem));
			}
			printf_s("                                                                                    \r");
			printf_s("\t\tcomplite\n");
		}

		//добавляем первые краевые и решаем СЛАУ
		{
			CSSD_Matrix<double, double> StiffnessMatrix_doubleSLAE;
			math::MakeCopyMatrix_A_into_B(StiffnessMatrix, StiffnessMatrix_doubleSLAE);


			printf("First boundary conditions...");
			for (int id_vertex = 0; id_vertex < solver_grid.boundary_vertexes.size(); id_vertex++)
			{
				auto boundary = &(solver_grid.boundary_vertexes[id_vertex]);

				Point<bool> is_take;
				Point<double> boundary_value = boundary->boundary_value(is_take, solver_grid.boundary_vertexes[id_vertex].GetIdNode(0));
				int global_id = boundary->GetDOFInLocalID(0);

				if (global_id < StiffnessMatrix.GetMatrixSize())
				{
					auto enter_value = [&StiffnessMatrix_doubleSLAE](int value_id, double value) ->void
					{
						StiffnessMatrix_doubleSLAE.X[value_id] = value;
						StiffnessMatrix_doubleSLAE.F[value_id] = value;
						StiffnessMatrix_doubleSLAE.Diag[value_id] = 1;
						//Обнуляем строку
						for (int i = 0; i < StiffnessMatrix_doubleSLAE.A_down[value_id].size(); i++)
						{
							StiffnessMatrix_doubleSLAE.A_down[value_id][i] = 0;
						}
						for (int i = 0; i < StiffnessMatrix_doubleSLAE.A_up[value_id].size(); i++)
						{
							StiffnessMatrix_doubleSLAE.A_up[value_id][i] = 0;
						}
						//обнуляем столбец
						for (int i = 0; i < StiffnessMatrix_doubleSLAE.GetMatrixSize(); i++)
						{
							//верхний треугольник
							if (i < value_id)
							{
								for (int jj = 0; jj < StiffnessMatrix_doubleSLAE.A_up[i].size(); jj++)
								{
									if (StiffnessMatrix_doubleSLAE.id_column_for_A_up[i][jj] == value_id)
									{
										StiffnessMatrix_doubleSLAE.F[i] -= StiffnessMatrix_doubleSLAE.A_up[i][jj] * StiffnessMatrix_doubleSLAE.F[value_id];
										StiffnessMatrix_doubleSLAE.A_up[i][jj] = 0;
									}
								}
							}
							//нижний треугольник
							if (i > value_id)
							{
								for (int jj = 0; jj < StiffnessMatrix_doubleSLAE.A_down[i].size(); jj++)
								{
									if (StiffnessMatrix_doubleSLAE.id_column_for_A_down[i][jj] == value_id)
									{
										StiffnessMatrix_doubleSLAE.F[i] -= StiffnessMatrix_doubleSLAE.A_down[i][jj] * StiffnessMatrix_doubleSLAE.F[value_id];
										StiffnessMatrix_doubleSLAE.A_down[i][jj] = 0;
									}
								}
							}
						}
						return;
					};

					if (is_take.x) enter_value(global_id * 3 + 0, boundary_value.x);
					if (is_take.y) enter_value(global_id * 3 + 1, boundary_value.y);
					if (is_take.z)
					{
						//if (boundary_value.z > 0) boundary_value.z = 1;
						enter_value(global_id * 3 + 2, boundary_value.z);
					}
				}

			}

			if (id_STEP == 0)
			{
				for (int i = 0; i < StiffnessMatrix_doubleSLAE.GetMatrixSize() / 3; i++)
				{
					StiffnessMatrix_doubleSLAE.X[i * 3 + 0] = 1e-10;
					StiffnessMatrix_doubleSLAE.X[i * 3 + 1] = 1e-10;
					StiffnessMatrix_doubleSLAE.X[i * 3 + 2] = 1e-10;
				}
			}
			else {
				for (int i = 0; i < StiffnessMatrix_doubleSLAE.GetMatrixSize() / 3; i++)
				{
					StiffnessMatrix_doubleSLAE.X[i * 3 + 0] = U_prev[i].x;
					StiffnessMatrix_doubleSLAE.X[i * 3 + 1] = U_prev[i].y;
					StiffnessMatrix_doubleSLAE.X[i * 3 + 2] = U_prev[i].z;
				}
			}

			printf("Soluting SLAY... (%d)\n", StiffnessMatrix_doubleSLAE.GetMatrixSize());
			int MaxSize = StiffnessMatrix_doubleSLAE.GetMatrixSize();
			MaxSize = MaxSize / 10 < 10 ? 10 : MaxSize / 10;
			std::vector<double> best_solution;
			math::MakeCopyVector_A_into_B(StiffnessMatrix_doubleSLAE.X, best_solution);
			double current_residual = 1, best_residual = 1e+25;
			int MAX_STEPS = 15;
			int ii = 1;
			CSSD_Matrix<double, double> Predcondor;
			Predcondor.PrecondorSSOR_summetric(0.75, StiffnessMatrix_doubleSLAE);
			for (int i = 0; i <= MAX_STEPS; i++)
			{
				double needed_residual = pow(10., -1 * (ii + 1));
				//current_residual /= 2.;
				printf_s("//---> I = %d/%d (full size %d) - needed residual %.2e\n", i, MAX_STEPS, StiffnessMatrix_doubleSLAE.GetMatrixSize(), needed_residual);


				//current_residual = abs(newSLAE.BiCG_Stab(MaxSize, needed_residual));
				//current_residual = abs(StiffnessMatrix_doubleSLAE.MSG(MaxSize, needed_residual));
				//current_residual = abs(StiffnessMatrix_doubleSLAE.MSG_Preconditioning(MaxSize, needed_residual, Predcondor));
				//current_residual = abs(newSLAE.BiCG_Stab_Preconditioning(MaxSize, needed_residual, Predcondor));

				if (current_residual < needed_residual)
				{
					i = 0;
					ii++;
					if (current_residual > best_residual + (1e-10))
					{
						//math::MakeCopyVector_A_into_B(best_solution, newSLAE.X);
						//printf_s("//---> BEST residual %.2e\n", best_residual);
						break;
					}
				}
				if (current_residual <= critical_residual)
				{
					best_residual = current_residual;
					break;
				}
				if (current_residual < best_residual)
				{
					best_residual = current_residual;
					math::MakeCopyVector_A_into_B(StiffnessMatrix_doubleSLAE.X, best_solution);
				}
			}
			math::MakeCopyVector_A_into_B(best_solution, StiffnessMatrix_doubleSLAE.X);
			printf_s("//---> BEST residual %.2e\n", best_residual);
			math::MakeCopyVector_A_into_B(StiffnessMatrix_doubleSLAE.X, StiffnessMatrix.X);
			math::MakeCopyVector_A_into_B(StiffnessMatrix_doubleSLAE.X, U_curr);

			if (best_residual >= 1)
			{
				printf_s("There are problems in the solution\n");
				int a;
				scanf_s("%d", &a);
			}
		}

		math::MakeCopyVector_A_into_B(U_prev, U_prevprev);
		math::MakeCopyVector_A_into_B(U_curr, U_prev);

		//output solution
		printf_s("Print the mech result into .dat file... ");
		//without deformations
		if (is_print_result) {
			FILE* fout_tech;
			char name_u_tech[5000];
			sprintf_s(name_u_tech, "%s/U_non_deformation.dat", result_directory);
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

				Point<double> U = solver_grid.GetSolutionInPoint(i, Centr, U_curr);
				Point<Point<double>> dU = solver_grid.GetDerevativeFromSolutionInPoint(i, Centr, U_curr);

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

				value[0][i] = sigma[0];
				value[1][i] = sigma[1];
				value[2][i] = sigma[2];

				value[3][i] = eps[0];
				value[4][i] = eps[1];
				value[5][i] = eps[2];

				value[3][i] = U.x;
				value[4][i] = U.y;
				value[5][i] = U.z;

				if (sigma_inv > sigma_inv_max)
				{
					sigma_inv_max = sigma_inv;
					elem_sigma_max = i;
				}
			}
			value[3].resize(solver_grid.GetDOFsCount());
			value[4].resize(solver_grid.GetDOFsCount());
			value[5].resize(solver_grid.GetDOFsCount());
			for (int i = 0; i < value[3].size(); i++)
			{
				value[3][i] = U_curr[i].x;
				value[4][i] = U_curr[i].y;
				value[5][i] = U_curr[i].z;

			}
			//solver_grid.MoveCoordinates(Solution);

			solver_grid.printTecPlot3D_DiffDomains(fout_tech, value, name_value, name_in_file);
			fclose(fout_tech);
		}
		printf_s("\t complite\n");

		//solver_grid.~Grid_forMech();
	}
}
void Solve_ElastodynamicsProblem_SelfDeform()
{
	FEM::Grid_forMech solver_grid; //output
	std::vector<Point<double>> Solution; //output

	bool is_STATIONARY = false;


	//char properties_file[1000] = { "E:/+cyl/800el/param.txt" };
	//char properties_file[1000] = { "E:/Box/100x10x200/BoxWithCracks/67k/param.txt" };
	//char properties_file[1000] = { "E:/Box/100x10x200/BoxWithCracks/638k/param.txt" };
	//char properties_file[1000] = { "E:/Box/100x50x200/param.txt" };
	//char properties_file[1000] = { "E:/Box/200x200x200/200k_fem/param_for_solver.txt" };
	//char properties_file[1000] = { "./param_for_solver.txt" };
	char properties_file[1000] = { "D:/Elastodynamic/homocyl/67k/param_for_solver.txt" };
	//char properties_file[1000] = { "base_properties.txt" };
	printf_s("The properties file contains the information:\n==========================================\n");
	printf_s("\tDo you need to print to a log file? 0-false/1-true\n");
	printf_s("\tinput: Mesh directory: box_200x200x200/90el\n");
	printf_s("\tproperties: Iterative process:\n");
	printf_s("\t\t<Start iteration> <End iteration>\n");
	printf_s("\t\t<Step size for Time>\n");
	printf_s("\t\t<Critical residual value for SLAE solution>\n");
	printf_s("\t<properties: Boundary values>\n");
	printf_s("\t<properties: Materials (E, v, rpho)>\n");
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
	math::ReadNonEmptyLine(f_properties, base_result_directory);
	math::SimpleGrid geo_grid; //input
	std::vector<std::vector<std::vector<int>>> boundary_faces;
	geo_grid.ReadFromNVTR(mesh_directory, 4);
	geo_grid.ReadFacesBoundaryNVTR(mesh_directory, boundary_faces);

	int start_iteration, end_iteration;
	{
		char _line[1000];
		math::ReadNonEmptyLine_forNumbers(f_properties, _line);
		std::vector<int> val;
		math::ParserStringToVectorInt(_line, val, " ");
		start_iteration = val[0];
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

			//
			/*first_boundaries[0].id_vertexes.resize(1);
			first_boundaries[0].id_vertexes[0] = 8;*/
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
		double _rpho;

		_material(double _E, double _v, double _rpho)
		{
			this->_E = _E;
			this->_v = _v;
			this->_rpho = _rpho;
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
		double _E, _v, _rpho;
		{
			char _line[1000];
			math::ReadNonEmptyLine_forNumbers(f_properties, _line);
			std::vector<float> val;
			math::ParserStringToVectorFloat(_line, val, " ");
			_E = val[0];
			_v = val[1];
			_rpho = val[2];
		}
		_materials.push_back(_material(_E, _v, _rpho));
	}

	fclose(f_properties);

	CreateDirectory((LPCTSTR)base_result_directory, NULL);
	bool res = false;
	int step_coef = 1;
	double TIME_h = step_size / step_coef*10;
	if (is_STATIONARY) {
		end_iteration = 1;
		start_iteration = 0;
	}
	std::vector<Point<double>> U_prev(geo_grid.xyz.size()), U_prevprev(geo_grid.xyz.size()), U_curr(geo_grid.xyz.size());
	math::InitializationVector(U_curr, 0);
	math::InitializationVector(U_prev, 0);
	math::InitializationVector(U_prevprev, 0);

	CSSD_Matrix<Tensor2Rank3D, Point<double>> StiffnessMatrix;
	printf("Initialization of grid...\t");
	solver_grid.Initialization(geo_grid, first_boundaries, second_boundary);
	for (int i = 0; i < _materials.size(); i++)
	{
		solver_grid.AddDomain();
		auto domain = solver_grid.GetDomain(i);
		domain->forThermal.rpho = _materials[i]._rpho;
		domain->forMech.SetE(_materials[i]._E);
		domain->forMech.SetV(_materials[i]._v);
	}
	printf_s("complite\n");
	printf("Creation the SLAE portrait...\t");
	solver_grid.CreationPortrait(StiffnessMatrix);
	printf_s("\t\tcomplite\n");

	for (int id_STEP = start_iteration; id_STEP < end_iteration * step_coef; id_STEP++)
	{
		printf_s("\n================= Start solution of %d STEP (time = %.2e) ================\n", id_STEP, TIME_h * id_STEP);
		double TIME_curr = TIME_h * id_STEP;
		bool is_print_result = false;
		char result_directory[1000];
		if (id_STEP % step_coef == 0)
		{
			is_print_result = true;
			sprintf_s(result_directory, sizeof(result_directory), "%s/STEP_%d_t=%.2e", base_result_directory, id_STEP, TIME_curr);
			CreateDirectory((LPCTSTR)result_directory, NULL);
		}

		//переопределяем источник
		if (true) {
			Point<double> power(0, 0, -1e+6 / (2*M_PI*0.004 * 0.004));
			double tay_plus = 0.0001; //половина длины импульса
			double tay_0 = 0.0001; //середина импульса
			double degree = 2; //степень функции
			std::function<Point<double>(Point<double>)> new_sourse_value = [TIME_curr, power, tay_0, degree, tay_plus](Point<double> X) -> Point<double>
			{
				//double sourse = exp(-pow(TIME_curr - tay_0, degree) / pow(tay_plus, degree));
				//return power;
				//return Point<double>(power.x * sourse, power.y * sourse, power.z * sourse);
				double sourse = TIME_curr < 14e-7 ? 10e+6 : 0;
				return Point<double>(0 * sourse, 0 * sourse, -1 * sourse / (2 * M_PI * 0.004 * 0.004));
			};
			if (second_boundary.size() > 0)
			{
				second_boundary[0].value = new_sourse_value;

				for (int i = 0; i < solver_grid.boundary_faces.size(); i++)
				{
					solver_grid.boundary_faces[i].boundary_value = second_boundary[solver_grid.boundary_faces[i].id_type].value;
				}
			}

			std::function<Point<double>(Point<bool>&, int)> new_sourse_value_D = [TIME_curr](Point<bool>& is_take, int id)->Point<double>
			{
				//double sourse = exp(-pow(TIME_curr - tay_0, degree) / pow(tay_plus, degree));
				//return power;
				//return Point<double>(power.x * sourse, power.y * sourse, power.z * sourse);
				is_take.x = true;
				is_take.y = true;
				is_take.z = true;
				double sourse = TIME_curr < 14e-7 ? 1e-8 : 0;
				return Point<double>(0 * sourse, 0 * sourse, -1 * sourse);
			};
			if (false && first_boundaries.size() > 1)
			{
				first_boundaries[1].value = new_sourse_value_D;

				for (int i = 0; i < solver_grid.boundary_vertexes.size(); i++)
				{
					if (solver_grid.boundary_vertexes[i].id_type == 1)
					{
						solver_grid.boundary_vertexes[i].boundary_value = first_boundaries[1].value;
					}
				}
			}
		}
		
		//printf_s("================= Create matrix ================\n");
		{
			std::function<std::vector<std::vector<double>>(int, Point<double>)> StiffnessCoef = [&](int elem, Point<double> X)->std::vector<std::vector<double>>
			{
				std::vector<std::vector<double>> D = solver_grid.GetDomain(solver_grid.GetElement(elem)->GetIdDomain())->forMech.GetD(3);
				/*for (int i = 0; i < D.size(); i++)
					for (int j = 0; j < D[i].size(); j++)
						D[i][j] *= 1;*/
				return D;
			};
			std::function<double(int, Point<double>)> MassCoef = [&](int elem, Point<double> X)->double
			{
				if (is_STATIONARY)
				{
					return 0;
				}
				else {
					auto domain = (solver_grid.GetDomain(solver_grid.GetElement(elem)->GetIdDomain()));
					return domain->forThermal.rpho / (TIME_h * TIME_h);
				}
			};
			std::function<Point<double>(int, Point<double>)> VolumeForсe = [&](int elem, Point<double> X)->Point<double>
			{
				if (is_STATIONARY)
				{
					return Point<double>(0, 0, 0);
				}
				else {
					Point<double> U_prev_in_X = solver_grid.GetSolutionInPoint(elem, X, U_prev);
					Point<double> U_prevprev_in_X = solver_grid.GetSolutionInPoint(elem, X, U_prevprev);

					return (U_prevprev_in_X * (-1) + U_prev_in_X * 2)
						* solver_grid.GetDomain(solver_grid.GetElement(elem)->GetIdDomain())->forThermal.rpho
						/ TIME_h / TIME_h;
				}
			};
			//SLAE assembling
			{
				printf("Matrix assembling...\n");
				StiffnessMatrix.ClearVariables();
				std::vector<DenseMatrix<Tensor2Rank3D, Point<double>>> local_SLAE_stiffness(solver_grid.GetElementsCount());
				for (int id_elem = 0; id_elem < solver_grid.GetElementsCount(); id_elem++)
				{
					if (id_elem % 1000 == 0)
						printf("Solve element[%d]\r", id_elem);
					auto element = solver_grid.GetElement(id_elem);

					std::function<std::vector<std::vector<double>>(Point<double>)> D = [&](Point<double> X) {return StiffnessCoef(id_elem, X); };
					std::function<double(Point<double>)> M = [&](Point<double> X) {return MassCoef(id_elem, X); };
					std::function<Point<double>(Point<double>)> F = [&](Point<double> X) {return VolumeForсe(id_elem, X); };

					element->SolveLocalMatrix(local_SLAE_stiffness[id_elem], D, M, F);
				}
				for (int id_elem = 0; id_elem < solver_grid.GetElementsCount(); id_elem++)
				{
					if (id_elem % 1000 == 0)
						printf("Add the local matrix of element[%d]\r", id_elem);
					StiffnessMatrix.SummPartOfMatrix(local_SLAE_stiffness[id_elem], *solver_grid.GetElementDOFs(id_elem));
				}
				printf_s("                                                                                    \r");
				printf_s("\t\tcomplite\n");
			}
		}

		//second boyndary condition
		for (int id_triangle = 0; id_triangle < solver_grid.boundary_faces.size(); id_triangle++)
		{
			std::vector<Point<double>> local_vector_SLAE;
			solver_grid.boundary_faces[id_triangle].SolveLocalBoundaryVector(local_vector_SLAE, solver_grid.boundary_faces[id_triangle].boundary_value);
			StiffnessMatrix.SummPartOfVector(local_vector_SLAE, *solver_grid.boundary_faces[id_triangle].GetElementDOFs());
		};

		//добавляем первые краевые и решаем СЛАУ
		{
			CSSD_Matrix<double, double> StiffnessMatrix_doubleSLAE;
			math::MakeCopyMatrix_A_into_B(StiffnessMatrix, StiffnessMatrix_doubleSLAE);


			printf("First boundary conditions...");
			for (int id_vertex = 0; id_vertex < solver_grid.boundary_vertexes.size(); id_vertex++)
			{
				auto boundary = &(solver_grid.boundary_vertexes[id_vertex]);

				Point<bool> is_take;
				Point<double> boundary_value = boundary->boundary_value(is_take, solver_grid.boundary_vertexes[id_vertex].GetIdNode(0));
				int global_id = boundary->GetDOFInLocalID(0);

				if (global_id < StiffnessMatrix.GetMatrixSize())
				{
					auto enter_value = [&StiffnessMatrix_doubleSLAE](int value_id, double value) ->void
					{
						StiffnessMatrix_doubleSLAE.X[value_id] = value;
						StiffnessMatrix_doubleSLAE.F[value_id] = value;
						StiffnessMatrix_doubleSLAE.Diag[value_id] = 1;
						//Обнуляем строку
						for (int i = 0; i < StiffnessMatrix_doubleSLAE.A_down[value_id].size(); i++)
						{
							StiffnessMatrix_doubleSLAE.A_down[value_id][i] = 0;
						}
						for (int i = 0; i < StiffnessMatrix_doubleSLAE.A_up[value_id].size(); i++)
						{
							StiffnessMatrix_doubleSLAE.A_up[value_id][i] = 0;
						}
						//обнуляем столбец
						for (int i = 0; i < StiffnessMatrix_doubleSLAE.GetMatrixSize(); i++)
						{
							//верхний треугольник
							if (i < value_id)
							{
								for (int jj = 0; jj < StiffnessMatrix_doubleSLAE.A_up[i].size(); jj++)
								{
									if (StiffnessMatrix_doubleSLAE.id_column_for_A_up[i][jj] == value_id)
									{
										StiffnessMatrix_doubleSLAE.F[i] -= StiffnessMatrix_doubleSLAE.A_up[i][jj] * StiffnessMatrix_doubleSLAE.F[value_id];
										StiffnessMatrix_doubleSLAE.A_up[i][jj] = 0;
									}
								}
							}
							//нижний треугольник
							if (i > value_id)
							{
								for (int jj = 0; jj < StiffnessMatrix_doubleSLAE.A_down[i].size(); jj++)
								{
									if (StiffnessMatrix_doubleSLAE.id_column_for_A_down[i][jj] == value_id)
									{
										StiffnessMatrix_doubleSLAE.F[i] -= StiffnessMatrix_doubleSLAE.A_down[i][jj] * StiffnessMatrix_doubleSLAE.F[value_id];
										StiffnessMatrix_doubleSLAE.A_down[i][jj] = 0;
									}
								}
							}
						}
						return;
					};

					if (is_take.x) enter_value(global_id * 3 + 0, boundary_value.x);
					if (is_take.y) enter_value(global_id * 3 + 1, boundary_value.y);
					if (is_take.z)
					{
						//if (boundary_value.z > 0) boundary_value.z = 1;
						enter_value(global_id * 3 + 2, boundary_value.z);
					}
				}

			}

			if (id_STEP == 0)
			{
				for (int i = 0; i < StiffnessMatrix_doubleSLAE.GetMatrixSize() / 3; i++)
				{
					StiffnessMatrix_doubleSLAE.X[i * 3 + 0] = 1e-10;
					StiffnessMatrix_doubleSLAE.X[i * 3 + 1] = 1e-10;
					StiffnessMatrix_doubleSLAE.X[i * 3 + 2] = 1e-10;
				}
			}
			else {
				for (int i = 0; i < StiffnessMatrix_doubleSLAE.GetMatrixSize() / 3; i++)
				{
					StiffnessMatrix_doubleSLAE.X[i * 3 + 0] = U_prev[i].x;
					StiffnessMatrix_doubleSLAE.X[i * 3 + 1] = U_prev[i].y;
					StiffnessMatrix_doubleSLAE.X[i * 3 + 2] = U_prev[i].z;
				}
			}

			printf("Soluting SLAY... (%d)\n", StiffnessMatrix_doubleSLAE.GetMatrixSize());
			int MaxSize = StiffnessMatrix_doubleSLAE.GetMatrixSize();
			MaxSize = MaxSize / 10 < 10 ? 10 : MaxSize / 10;
			std::vector<double> best_solution;
			math::MakeCopyVector_A_into_B(StiffnessMatrix_doubleSLAE.X, best_solution);
			double current_residual = 1, best_residual = 1e+25;
			int MAX_STEPS = 15;
			int ii = 1;
			CSSD_Matrix<double, double> Predcondor;
			Predcondor.PrecondorSSOR_summetric(0.75, StiffnessMatrix_doubleSLAE);
			for (int i = 0; i <= MAX_STEPS; i++)
			{
				double needed_residual = pow(10., -1 * (ii + 1));
				//current_residual /= 2.;
				printf_s("//---> I = %d/%d (full size %d) - needed residual %.2e\n", i, MAX_STEPS, StiffnessMatrix_doubleSLAE.GetMatrixSize(), needed_residual);


				//current_residual = abs(newSLAE.BiCG_Stab(MaxSize, needed_residual));
				current_residual = abs(StiffnessMatrix_doubleSLAE.MSG(MaxSize, needed_residual));
				//current_residual = abs(StiffnessMatrix_doubleSLAE.MSG_Preconditioning(MaxSize, needed_residual, Predcondor));
				//current_residual = abs(newSLAE.BiCG_Stab_Preconditioning(MaxSize, needed_residual, Predcondor));

				if (current_residual < needed_residual)
				{
					i = 0;
					ii++;
					if (current_residual > best_residual + (1e-10))
					{
						//math::MakeCopyVector_A_into_B(best_solution, newSLAE.X);
						//printf_s("//---> BEST residual %.2e\n", best_residual);
						break;
					}
				}
				if (current_residual <= critical_residual)
				{
					best_residual = current_residual;
					break;
				}
				if (current_residual < best_residual)
				{
					best_residual = current_residual;
					math::MakeCopyVector_A_into_B(StiffnessMatrix_doubleSLAE.X, best_solution);
				}
			}
			math::MakeCopyVector_A_into_B(best_solution, StiffnessMatrix_doubleSLAE.X);
			printf_s("//---> BEST residual %.2e\n", best_residual);
			math::MakeCopyVector_A_into_B(StiffnessMatrix_doubleSLAE.X, StiffnessMatrix.X);
			math::MakeCopyVector_A_into_B(StiffnessMatrix_doubleSLAE.X, U_curr);

			if (best_residual >= 1)
			{
				printf_s("There are problems in the solution\n");
				int a;
				scanf_s("%d", &a);
			}
		}

		math::MakeCopyVector_A_into_B(U_prev, U_prevprev);
		math::MakeCopyVector_A_into_B(U_curr, U_prev);

		//обновляем сетку
		if(true) {
			solver_grid.MoveCoordinates(U_curr);
			for (int i = 0; i < solver_grid.GetElementsCount(); i++)
			{
				solver_grid.GetElement(i)->SolveAlphaMatrix();
			}
		}

		//output solution
		printf_s("Print the mech result into .dat file... ");
		//without deformations
		if (is_print_result) {
			FILE* fout_tech;
			char name_u_tech[5000];
			sprintf_s(name_u_tech, "%s/U_non_deformation.dat", result_directory);
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

				Point<double> U = solver_grid.GetSolutionInPoint(i, Centr, U_curr);
				Point<Point<double>> dU = solver_grid.GetDerevativeFromSolutionInPoint(i, Centr, U_curr);

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

				value[0][i] = sigma[0];
				value[1][i] = sigma[1];
				value[2][i] = sigma[2];

				value[3][i] = eps[0];
				value[4][i] = eps[1];
				value[5][i] = eps[2];

				value[3][i] = U.x;
				value[4][i] = U.y;
				value[5][i] = U.z;

				if (sigma_inv > sigma_inv_max)
				{
					sigma_inv_max = sigma_inv;
					elem_sigma_max = i;
				}
			}
			value[3].resize(solver_grid.GetDOFsCount());
			value[4].resize(solver_grid.GetDOFsCount());
			value[5].resize(solver_grid.GetDOFsCount());
			for (int i = 0; i < value[3].size(); i++)
			{
				value[3][i] = U_curr[i].x;
				value[4][i] = U_curr[i].y;
				value[5][i] = U_curr[i].z;

			}
			//solver_grid.MoveCoordinates(Solution);

			solver_grid.printTecPlot3D_DiffDomains(fout_tech, value, name_value, name_in_file);
			fclose(fout_tech);
		}
		printf_s("\t complite\n");

		//solver_grid.~Grid_forMech();
	}
}
void Solve_ElastodynamicsProblem_SelfDeform_2Order()
{
	FEM::Grid_forMech_Order2 solver_grid; //output
	//FEM::Grid_forMech solver_grid; //output

	bool is_STATIONARY = false;


	//char properties_file[1000] = { "E:/+cyl/800el/param.txt" };
	//char properties_file[1000] = { "E:/Box/100x10x200/BoxWithCracks/67k/param.txt" };
	//char properties_file[1000] = { "E:/Box/100x10x200/BoxWithCracks/638k/param.txt" };
	//char properties_file[1000] = { "E:/Box/100x50x200/param.txt" };
	//char properties_file[1000] = { "E:/Box/200x200x200/200k_fem/param_for_solver.txt" };
	//char properties_file[1000] = { "./param_for_solver.txt" };
	//char properties_file[1000] = { "D:/Elastodynamic/homocyl/74k/param_for_solver.txt" };
	char properties_file[1000] = { "D:/Elastodynamic/pores_cyl/param_for_solver.txt" };
	//char properties_file[1000] = { "base_properties.txt" };
	printf_s("The properties file contains the information:\n==========================================\n");
	printf_s("\tDo you need to print to a log file? 0-false/1-true\n");
	printf_s("\tinput: Mesh directory: box_200x200x200/90el\n");
	printf_s("\tproperties: Iterative process:\n");
	printf_s("\t\t<Start iteration> <End iteration>\n");
	printf_s("\t\t<Step size for Time>\n");
	printf_s("\t\t<Critical residual value for SLAE solution>\n");
	printf_s("\t<properties: Boundary values>\n");
	printf_s("\t<properties: Materials (E, v, rpho)>\n");
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
	math::ReadNonEmptyLine(f_properties, base_result_directory);
	math::SimpleGrid geo_grid; //input
	std::vector<std::vector<std::vector<int>>> boundary_faces;
	geo_grid.ReadFromNVTR(mesh_directory, 4);
	geo_grid.ReadFacesBoundaryNVTR(mesh_directory, boundary_faces);

	int start_iteration, end_iteration;
	{
		char _line[1000];
		math::ReadNonEmptyLine_forNumbers(f_properties, _line);
		std::vector<int> val;
		math::ParserStringToVectorInt(_line, val, " ");
		start_iteration = val[0];
		end_iteration = val[1];
	}
	double step_size;
	int step_for_out;
	{
		char _line[1000];
		math::ReadNonEmptyLine_forNumbers(f_properties, _line);
		std::vector<float> val;
		math::ParserStringToVectorFloat(_line, val, " ");
		step_size = val[0];
		std::vector<int> val2;
		math::ReadNonEmptyLine_forNumbers(f_properties, _line);
		math::ParserStringToVectorInt(_line, val2, " ");
		step_for_out = (int)val2[0];
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
				if (curr_i > 0)
				{
					std::vector<float> _val;
					_tmp_line[curr_i] = '\0';
					math::ParserStringToVectorFloat(_tmp_line, _val, " *");
					printf_s("                                                                                                                  \r");
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

				//
				/*if (id_type == 0)
				{
					for (int i = 0; i < geo_grid.xyz.size(); i++)
					{
						if (math::IsEqual(0.0, geo_grid.xyz[i].z))
						{
							tmp_vert.push_back(i);
						}
					}
				}*/

				math::MakeQuickSort(tmp_vert);
				math::MakeRemovalOfDuplication(tmp_vert, first_boundaries[id_type].id_vertexes);
			}

			//
			/*first_boundaries[0].id_vertexes.resize(1);
			first_boundaries[0].id_vertexes[0] = 8;*/
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
		double _rpho;

		_material(double _E, double _v, double _rpho)
		{
			this->_E = _E;
			this->_v = _v;
			this->_rpho = _rpho;
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
		double _E, _v, _rpho;
		{
			char _line[1000];
			math::ReadNonEmptyLine_forNumbers(f_properties, _line);
			std::vector<float> val;
			math::ParserStringToVectorFloat(_line, val, " ");
			_E = val[0];
			_v = val[1];
			_rpho = val[2];
		}
		_materials.push_back(_material(_E, _v, _rpho));
	}

	fclose(f_properties);

	CSSD_Matrix<Tensor2Rank3D, Point<double>> StiffnessMatrix;
	printf("Initialization of grid...\t");
	solver_grid.Initialization(geo_grid, first_boundaries, second_boundary);
	for (int i = 0; i < _materials.size(); i++)
	{
		solver_grid.AddDomain();
		auto domain = solver_grid.GetDomain(i);
		domain->forThermal.rpho = _materials[i]._rpho;
		domain->forMech.SetE(_materials[i]._E);
		domain->forMech.SetV(_materials[i]._v);
	}
	printf_s("complite\n");
	printf("Creation the SLAE portrait...\t");
	solver_grid.CreationPortrait(StiffnessMatrix);
	printf_s("\t\tcomplite\n");

	CreateDirectory((LPCTSTR)base_result_directory, NULL);
	bool res = false;
	double TIME_h = step_size;
	if (is_STATIONARY) {
		end_iteration = 1;
		start_iteration = 0;
	}
	std::vector<Point<double>> U_prev(solver_grid.GetDOFsCount()), U_prevprev(solver_grid.GetDOFsCount()), U_curr(solver_grid.GetDOFsCount());
	math::InitializationVector(U_curr, 0);
	math::InitializationVector(U_prev, 0);
	math::InitializationVector(U_prevprev, 0);

	for (int id_STEP = start_iteration; id_STEP < end_iteration; id_STEP++)
	{
		if (id_STEP == 1) TIME_h = step_size;
		printf_s("\n================= Start solution of %d STEP (time = %.2e) ================\n", id_STEP, TIME_h * id_STEP);
		double TIME_curr = TIME_h * (id_STEP);
		bool is_print_result = false;
		char result_directory[1000];
		if (id_STEP % step_for_out == 0)
		{
			is_print_result = true;
			sprintf_s(result_directory, sizeof(result_directory), "%s/STEP_%d_t=%.2e", base_result_directory, id_STEP, TIME_curr);
			CreateDirectory((LPCTSTR)result_directory, NULL);
		}

		//переопределяем источник
		if (true) {
			Point<double> power(0, 0, -10e+6 / (2 * M_PI * 0.004 * 0.004));
			double tay_plus = 0.0001; //половина длины импульса
			double tay_0 = 0.0001; //середина импульса
			double degree = 2; //степень функции
			std::function<Point<double>(Point<double>)> new_sourse_value = [TIME_curr, power, tay_0, degree, tay_plus, TIME_h](Point<double> X) -> Point<double>
			{
				//double sourse = exp(-pow(TIME_curr - tay_0, degree) / pow(tay_plus, degree));
				//return power;
				//return Point<double>(power.x * sourse, power.y * sourse, power.z * sourse);
				double sourse = TIME_curr < 14e-7 ? 10e+6 : 0;
				return Point<double>(0 * sourse, 0 * sourse, -1 * sourse /** TIME_h / 14.0e-7*/);
			};
			if (second_boundary.size() > 0)
			{
				second_boundary[0].value = new_sourse_value;

				for (int i = 0; i < solver_grid.boundary_faces.size(); i++)
				{
					solver_grid.boundary_faces[i].boundary_value = second_boundary[solver_grid.boundary_faces[i].id_type].value;
				}
			}

			std::function<Point<double>(Point<bool>&, int)> new_sourse_value_D = [TIME_curr](Point<bool>& is_take, int id)->Point<double>
			{
				//double sourse = exp(-pow(TIME_curr - tay_0, degree) / pow(tay_plus, degree));
				//return power;
				//return Point<double>(power.x * sourse, power.y * sourse, power.z * sourse);
				is_take.x = true;
				is_take.y = true;
				is_take.z = true;
				double sourse = TIME_curr < 14e-7 ? 1e-8 : 0;
				return Point<double>(0 * sourse, 0 * sourse, -1 * sourse);
			};
			if (false && first_boundaries.size() > 1)
			{
				first_boundaries[1].value = new_sourse_value_D;

				for (int i = 0; i < solver_grid.boundary_vertexes.size(); i++)
				{
					if (solver_grid.boundary_vertexes[i].id_type == 1)
					{
						solver_grid.boundary_vertexes[i].boundary_value = first_boundaries[1].value;
					}
				}
			}
		}

		//printf_s("================= Create matrix ================\n");
		{
			std::function<std::vector<std::vector<double>>(int, Point<double>)> StiffnessCoef = [&](int elem, Point<double> X)->std::vector<std::vector<double>>
			{
				std::vector<std::vector<double>> D = solver_grid.GetDomain(solver_grid.GetElement(elem)->GetIdDomain())->forMech.GetD(3);
				/*for (int i = 0; i < D.size(); i++)
					for (int j = 0; j < D[i].size(); j++)
						D[i][j] *= 1;*/
				return D;
			};
			std::function<double(int, Point<double>)> MassCoef = [&](int elem, Point<double> X)->double
			{
				if (is_STATIONARY)
				{
					return 0;
				}
				else {
					auto domain = (solver_grid.GetDomain(solver_grid.GetElement(elem)->GetIdDomain()));
					return domain->forThermal.rpho / (TIME_h * TIME_h);
				}
			};
			std::function<Point<double>(int, Point<double>)> VolumeForсe = [&](int elem, Point<double> X)->Point<double>
			{
				if (is_STATIONARY)
				{
					return Point<double>(0, 1e+10, 0);
					//return Point<double>(0, 0, -1e+10);
				}
				else {
					Point<double> U_prev_in_X = solver_grid.GetSolutionInPoint(elem, X, U_prev);
					Point<double> U_prevprev_in_X = solver_grid.GetSolutionInPoint(elem, X, U_prevprev);

					return (U_prevprev_in_X * (-1) + U_prev_in_X * 2)
						* solver_grid.GetDomain(solver_grid.GetElement(elem)->GetIdDomain())->forThermal.rpho
						/ TIME_h / TIME_h;
				}
			};
			//SLAE assembling
			{
				

				printf("Matrix assembling...\n");
				StiffnessMatrix.ClearVariables();
				std::vector<DenseMatrix<Tensor2Rank3D, Point<double>>> local_SLAE_stiffness(solver_grid.GetElementsCount());
				
				//#pragma omp parallel num_threads(8)
				omp_set_num_threads(math::NUM_THREADS);
#pragma omp parallel for
				/*for (int i = 0; i < 150; i++)
				{
					printf_s("%d\n", i);
				}*/
				for (int id_elem = 0; id_elem < solver_grid.GetElementsCount(); id_elem++)
				{
					if (id_elem % 1000 == 0)
						printf("Solve element[%d]\r", id_elem);
					auto element = solver_grid.GetElement(id_elem);

					std::function<std::vector<std::vector<double>>(Point<double>)> D = [&](Point<double> X) {return StiffnessCoef(id_elem, X); };
					std::function<double(Point<double>)> M = [&](Point<double> X) {return MassCoef(id_elem, X); };
					std::function<Point<double>(Point<double>)> F = [&](Point<double> X) {return VolumeForсe(id_elem, X); };

					element->SolveLocalMatrix(local_SLAE_stiffness[id_elem], D, M, F);
				}

				for (int id_elem = 0; id_elem < solver_grid.GetElementsCount(); id_elem++)
				{
					if (id_elem % 1000 == 0)
						printf("Add the local matrix of element[%d]\r", id_elem);
					StiffnessMatrix.SummPartOfMatrix(local_SLAE_stiffness[id_elem], *solver_grid.GetElementDOFs(id_elem));
				}
				printf_s("                                                                                    \r");
				printf_s("\t\tcomplite\n");			
			}
		}

		//second boyndary condition
		for (int id_triangle = 0; id_triangle < solver_grid.boundary_faces.size(); id_triangle++)
		{
			std::vector<Point<double>> local_vector_SLAE;
			solver_grid.boundary_faces[id_triangle].SolveLocalBoundaryVector(local_vector_SLAE, solver_grid.boundary_faces[id_triangle].boundary_value);
			StiffnessMatrix.SummPartOfVector(local_vector_SLAE, *solver_grid.boundary_faces[id_triangle].GetElementDOFs());
		};

		//добавляем первые краевые и решаем СЛАУ
		{
			CSSD_Matrix<double, double> StiffnessMatrix_doubleSLAE;
			math::MakeCopyMatrix_A_into_B(StiffnessMatrix, StiffnessMatrix_doubleSLAE);

			//double start = omp_get_wtime();

			printf("First boundary conditions...\n");
			if (false) { //old version
				for (int id_vertex = 0; id_vertex < solver_grid.boundary_vertexes.size(); id_vertex++)
				{
					if (id_vertex % 100 == 0)
					{
						printf("\tcurrent %d/%d\r", id_vertex, solver_grid.boundary_vertexes.size());
					}
					auto boundary = &(solver_grid.boundary_vertexes[id_vertex]);

					int id_type = boundary->id_type;
					Point<bool> is_take;
					Point<double> boundary_value = boundary->boundary_value(is_take, solver_grid.boundary_vertexes[id_vertex].GetIdNode(0));
					int global_id = boundary->GetDOFInLocalID(0);

					//printf("global_id = %d\tid_type = %d\tvalue=(%.2e; %.2e; %.2e)\n", global_id, id_type, boundary_value.x, boundary_value.y, boundary_value.z);

					if (global_id < StiffnessMatrix.GetMatrixSize())
					{
						auto enter_value = [&StiffnessMatrix_doubleSLAE](int value_id, double value) ->void
						{
							StiffnessMatrix_doubleSLAE.X[value_id] = value;
							StiffnessMatrix_doubleSLAE.F[value_id] = value;
							StiffnessMatrix_doubleSLAE.Diag[value_id] = 1;
							//Обнуляем строку
							for (int i = 0; i < StiffnessMatrix_doubleSLAE.A_down[value_id].size(); i++)
							{
								StiffnessMatrix_doubleSLAE.A_down[value_id][i] = 0;
							}
							for (int i = 0; i < StiffnessMatrix_doubleSLAE.A_up[value_id].size(); i++)
							{
								StiffnessMatrix_doubleSLAE.A_up[value_id][i] = 0;
							}
							//обнуляем столбец
							for (int i = 0; i < StiffnessMatrix_doubleSLAE.GetMatrixSize(); i++)
							{
								//верхний треугольник
								if (i < value_id)
								{
									int jj = math::GetPositionInSortVector(StiffnessMatrix_doubleSLAE.id_column_for_A_up[i], value_id);
									if (jj >= 0)
									{
										StiffnessMatrix_doubleSLAE.F[i] -= StiffnessMatrix_doubleSLAE.A_up[i][jj] * StiffnessMatrix_doubleSLAE.F[value_id];
										StiffnessMatrix_doubleSLAE.A_up[i][jj] = 0;
									}
									/*for (int jj = 0; jj < StiffnessMatrix_doubleSLAE.A_up[i].size(); jj++)
									{
										if (StiffnessMatrix_doubleSLAE.id_column_for_A_up[i][jj] == value_id)
										{
											StiffnessMatrix_doubleSLAE.F[i] -= StiffnessMatrix_doubleSLAE.A_up[i][jj] * StiffnessMatrix_doubleSLAE.F[value_id];
											StiffnessMatrix_doubleSLAE.A_up[i][jj] = 0;
											break;
										}
									}*/
								}
								//нижний треугольник
								if (i > value_id)
								{
									int jj = math::GetPositionInSortVector(StiffnessMatrix_doubleSLAE.id_column_for_A_down[i], value_id);
									if (jj >= 0)
									{
										StiffnessMatrix_doubleSLAE.F[i] -= StiffnessMatrix_doubleSLAE.A_down[i][jj] * StiffnessMatrix_doubleSLAE.F[value_id];
										StiffnessMatrix_doubleSLAE.A_down[i][jj] = 0;
									}

									/*for (int jj = 0; jj < StiffnessMatrix_doubleSLAE.A_down[i].size(); jj++)
									{
										if (StiffnessMatrix_doubleSLAE.id_column_for_A_down[i][jj] == value_id)
										{
											StiffnessMatrix_doubleSLAE.F[i] -= StiffnessMatrix_doubleSLAE.A_down[i][jj] * StiffnessMatrix_doubleSLAE.F[value_id];
											StiffnessMatrix_doubleSLAE.A_down[i][jj] = 0;
											break;
										}
									}*/
								}
							}
							return;
						};

						if (is_take.x) enter_value(global_id * 3 + 0, boundary_value.x);
						if (is_take.y) enter_value(global_id * 3 + 1, boundary_value.y);
						if (is_take.z)
						{
							//if (boundary_value.z > 0) boundary_value.z = 1;
							enter_value(global_id * 3 + 2, boundary_value.z);
						}
					}

				}
			}
			if (true) { //new version
				//обнуляем строки
				for (int id_vertex = 0; id_vertex < solver_grid.boundary_vertexes.size(); id_vertex++)
				{
					auto boundary = &(solver_grid.boundary_vertexes[id_vertex]);

					int id_type = boundary->id_type;
					Point<bool> is_take;
					Point<double> boundary_value = boundary->boundary_value(is_take, solver_grid.boundary_vertexes[id_vertex].GetIdNode(0));
					int global_id = boundary->GetDOFInLocalID(0);

					/*int id = solver_grid.boundary_vertexes[id_vertex].GetIdNode(0);
					Point<double> node;
					if (id < solver_grid.GetVertexCount())
					{
						node = solver_grid.GetCoordinateViaID(id);
					}

					printf("global_id = %d\tid_type = %d\tvalue=(%.2e; %.2e; %.2e)\tcoord=(%.2e; %.2e; %.2e)\n", 
						global_id, id_type, boundary_value.x, boundary_value.y, boundary_value.z,
						node.x, node.y, node.z);*/

					if (global_id < StiffnessMatrix.GetMatrixSize())
					{
						auto enter_value = [&StiffnessMatrix_doubleSLAE](int value_id, double value) ->void
						{
							StiffnessMatrix_doubleSLAE.X[value_id] = value;
							StiffnessMatrix_doubleSLAE.F[value_id] = value;
							StiffnessMatrix_doubleSLAE.Diag[value_id] = 1;
							//Обнуляем строку
							for (int i = 0; i < StiffnessMatrix_doubleSLAE.A_down[value_id].size(); i++)
							{
								StiffnessMatrix_doubleSLAE.A_down[value_id][i] = 0;
							}
							for (int i = 0; i < StiffnessMatrix_doubleSLAE.A_up[value_id].size(); i++)
							{
								StiffnessMatrix_doubleSLAE.A_up[value_id][i] = 0;
							}
							return;
						};

						if (is_take.x) enter_value(global_id * 3 + 0, boundary_value.x);
						if (is_take.y) enter_value(global_id * 3 + 1, boundary_value.y);
						if (is_take.z) enter_value(global_id * 3 + 2, boundary_value.z);
					}

				}
				//симметризация
				omp_set_num_threads(math::NUM_THREADS);
#pragma omp parallel for 
				for(int id_row = 0; id_row < StiffnessMatrix_doubleSLAE.GetMatrixSize(); id_row++)
				{
					if (id_row % 100 == 0)
					{
						printf("\tcurrent %d/%d\r", id_row, StiffnessMatrix_doubleSLAE.GetMatrixSize());
					}
					int iterator_in_boundary = 0;
					for (int jj = 0; jj < StiffnessMatrix_doubleSLAE.id_column_for_A_up[id_row].size(); jj++)
					{
						int id_column = StiffnessMatrix_doubleSLAE.id_column_for_A_up[id_row][jj];
						for (/*iterator_in_boundary = 0*/; iterator_in_boundary < solver_grid.boundary_vertexes.size(); iterator_in_boundary++)
						{
							if (solver_grid.boundary_vertexes[iterator_in_boundary].GetDOFInLocalID(0) * 3 + 0 == id_column)
							{
								Point<bool> is_take;
								Point<double> boundary_value = solver_grid.boundary_vertexes[iterator_in_boundary].boundary_value(is_take, solver_grid.boundary_vertexes[iterator_in_boundary].GetIdNode(0));

								if (is_take.x)
								{
									StiffnessMatrix_doubleSLAE.F[id_row] -= StiffnessMatrix_doubleSLAE.A_up[id_row][jj] * StiffnessMatrix_doubleSLAE.F[id_column];
									StiffnessMatrix_doubleSLAE.A_up[id_row][jj] = 0;
								}

								break;
							}
							if (solver_grid.boundary_vertexes[iterator_in_boundary].GetDOFInLocalID(0) * 3 + 1 == id_column)
							{
								Point<bool> is_take;
								Point<double> boundary_value = solver_grid.boundary_vertexes[iterator_in_boundary].boundary_value(is_take, solver_grid.boundary_vertexes[iterator_in_boundary].GetIdNode(0));

								if (is_take.y)
								{
									StiffnessMatrix_doubleSLAE.F[id_row] -= StiffnessMatrix_doubleSLAE.A_up[id_row][jj] * StiffnessMatrix_doubleSLAE.F[id_column];
									StiffnessMatrix_doubleSLAE.A_up[id_row][jj] = 0;
								}

								break;
							}
							if (solver_grid.boundary_vertexes[iterator_in_boundary].GetDOFInLocalID(0) * 3 + 2 == id_column)
							{
								Point<bool> is_take;
								Point<double> boundary_value = solver_grid.boundary_vertexes[iterator_in_boundary].boundary_value(is_take, solver_grid.boundary_vertexes[iterator_in_boundary].GetIdNode(0));

								if (is_take.z)
								{
									StiffnessMatrix_doubleSLAE.F[id_row] -= StiffnessMatrix_doubleSLAE.A_up[id_row][jj] * StiffnessMatrix_doubleSLAE.F[id_column];
									StiffnessMatrix_doubleSLAE.A_up[id_row][jj] = 0;
								}

								break;
							}
							if (solver_grid.boundary_vertexes[iterator_in_boundary].GetDOFInLocalID(0) * 3 > id_column)
							{
								//iterator_in_boundary--;
								break;
							}
						}
					}

					iterator_in_boundary = 0;
					for (int jj = 0; jj < StiffnessMatrix_doubleSLAE.id_column_for_A_down[id_row].size(); jj++)
					{
						int id_column = StiffnessMatrix_doubleSLAE.id_column_for_A_down[id_row][jj];
						for (/*iterator_in_boundary = 0*/; iterator_in_boundary < solver_grid.boundary_vertexes.size(); iterator_in_boundary++)
						{
							if (solver_grid.boundary_vertexes[iterator_in_boundary].GetDOFInLocalID(0) * 3 + 0 == id_column)
							{
								Point<bool> is_take;
								Point<double> boundary_value = solver_grid.boundary_vertexes[iterator_in_boundary].boundary_value(is_take, solver_grid.boundary_vertexes[iterator_in_boundary].GetIdNode(0));

								if (is_take.x)
								{
									StiffnessMatrix_doubleSLAE.F[id_row] -= StiffnessMatrix_doubleSLAE.A_down[id_row][jj] * StiffnessMatrix_doubleSLAE.F[id_column];
									StiffnessMatrix_doubleSLAE.A_down[id_row][jj] = 0;
								}

								break;
							}
							if (solver_grid.boundary_vertexes[iterator_in_boundary].GetDOFInLocalID(0) * 3 + 1 == id_column)
							{
								Point<bool> is_take;
								Point<double> boundary_value = solver_grid.boundary_vertexes[iterator_in_boundary].boundary_value(is_take, solver_grid.boundary_vertexes[iterator_in_boundary].GetIdNode(0));

								if (is_take.y)
								{
									StiffnessMatrix_doubleSLAE.F[id_row] -= StiffnessMatrix_doubleSLAE.A_down[id_row][jj] * StiffnessMatrix_doubleSLAE.F[id_column];
									StiffnessMatrix_doubleSLAE.A_down[id_row][jj] = 0;
								}

								break;
							}
							if (solver_grid.boundary_vertexes[iterator_in_boundary].GetDOFInLocalID(0) * 3 + 2 == id_column)
							{
								Point<bool> is_take;
								Point<double> boundary_value = solver_grid.boundary_vertexes[iterator_in_boundary].boundary_value(is_take, solver_grid.boundary_vertexes[iterator_in_boundary].GetIdNode(0));

								if (is_take.z)
								{
									StiffnessMatrix_doubleSLAE.F[id_row] -= StiffnessMatrix_doubleSLAE.A_down[id_row][jj] * StiffnessMatrix_doubleSLAE.F[id_column];
									StiffnessMatrix_doubleSLAE.A_down[id_row][jj] = 0;
								}

								break;
							}
							if (solver_grid.boundary_vertexes[iterator_in_boundary].GetDOFInLocalID(0) * 3 > id_column)
							{
								//iterator_in_boundary--;
								break;
							}
						}
					}
				}
			}
			printf("\t                                                                              \r");
			printf("\t\tcomplite\n");

			/*double end = omp_get_wtime();
			printf_s("\nTime for assembling %lf sec\n", end - start);
			int ttt;
			scanf_s("%d", &ttt);*/

			if (id_STEP == 0)
			{
				for (int i = 0; i < StiffnessMatrix_doubleSLAE.GetMatrixSize() / 3; i++)
				{
					StiffnessMatrix_doubleSLAE.X[i * 3 + 0] = 1e-10;
					StiffnessMatrix_doubleSLAE.X[i * 3 + 1] = 1e-10;
					StiffnessMatrix_doubleSLAE.X[i * 3 + 2] = 1e-10;
				}
			}
			else {
				for (int i = 0; i < StiffnessMatrix_doubleSLAE.GetMatrixSize() / 3; i++)
				{
					StiffnessMatrix_doubleSLAE.X[i * 3 + 0] = U_prev[i].x;
					StiffnessMatrix_doubleSLAE.X[i * 3 + 1] = U_prev[i].y;
					StiffnessMatrix_doubleSLAE.X[i * 3 + 2] = U_prev[i].z;
				}
			}

			
			printf("Soluting SLAY... (%d)\n", StiffnessMatrix_doubleSLAE.GetMatrixSize());
			int MaxSize = StiffnessMatrix_doubleSLAE.GetMatrixSize();
			MaxSize = MaxSize / 10 < 10 ? 10 : MaxSize / 10;
			MaxSize = min(1000, MaxSize);
			std::vector<double> best_solution;
			math::MakeCopyVector_A_into_B(StiffnessMatrix_doubleSLAE.X, best_solution);
			double current_residual = 1, best_residual = 1e+25;
			int MAX_STEPS = 5;
			int ii = 1;
			CSSD_Matrix<double, double> Predcondor;
			//Predcondor.PrecondorSSOR_summetric(0.75, StiffnessMatrix_doubleSLAE);
			for (int i = 0; i <= MAX_STEPS; i++)
			{
				double needed_residual = pow(10., -1 * (ii + 1));
				//current_residual /= 2.;
				printf_s("//---> I = %d/%d (full size %d) - needed residual %.2e\n", i, MAX_STEPS, StiffnessMatrix_doubleSLAE.GetMatrixSize(), needed_residual);

				/*for (int i = 0; i < 10; i++)
				{
					current_residual = abs(StiffnessMatrix_doubleSLAE.BiCG_Stab(MaxSize, needed_residual));
				}
				current_residual = abs(StiffnessMatrix_doubleSLAE.MSG(MaxSize, needed_residual));*/

				for (int i = 0; i < 5 && current_residual > needed_residual; i++)
				{
					current_residual = abs(StiffnessMatrix_doubleSLAE.BiCG_Stab(MaxSize, needed_residual));
					if (current_residual < needed_residual) break;
					current_residual = abs(StiffnessMatrix_doubleSLAE.BiCG_Stab(MaxSize, needed_residual));
					if (current_residual < needed_residual) break;
					current_residual = abs(StiffnessMatrix_doubleSLAE.BiCG_Stab(MaxSize, needed_residual));
					if (current_residual < needed_residual) break;

					current_residual = abs(StiffnessMatrix_doubleSLAE.MSG(100, needed_residual));
				}

				//current_residual = abs(StiffnessMatrix_doubleSLAE.MSG_Preconditioning(MaxSize, needed_residual, Predcondor));
				//current_residual = abs(newSLAE.BiCG_Stab_Preconditioning(MaxSize, needed_residual, Predcondor));

				if (current_residual < needed_residual)
				{
					i = 0;
					ii++;
					if (current_residual > best_residual + (1e-10))
					{
						//math::MakeCopyVector_A_into_B(best_solution, newSLAE.X);
						//printf_s("//---> BEST residual %.2e\n", best_residual);
						break;
					}
				}
				if (current_residual <= critical_residual)
				{
					best_residual = current_residual;
					break;
				}
				if (current_residual < best_residual)
				{
					best_residual = current_residual;
					math::MakeCopyVector_A_into_B(StiffnessMatrix_doubleSLAE.X, best_solution);
				}
			}
			math::MakeCopyVector_A_into_B(best_solution, StiffnessMatrix_doubleSLAE.X);
			printf_s("//---> BEST residual %.2e\n", best_residual);
			math::MakeCopyVector_A_into_B(StiffnessMatrix_doubleSLAE.X, StiffnessMatrix.X);
			math::MakeCopyVector_A_into_B(StiffnessMatrix_doubleSLAE.X, U_curr);

			if (best_residual >= 1)
			{
				printf_s("There are problems in the solution\n");
				int a;
				scanf_s("%d", &a);
			}
		}

		math::MakeCopyVector_A_into_B(U_prev, U_prevprev);
		math::MakeCopyVector_A_into_B(U_curr, U_prev);

		//обновляем сетку
		if (false) {
			solver_grid.MoveCoordinates(U_curr);
			for (int i = 0; i < solver_grid.GetElementsCount(); i++)
			{
				solver_grid.GetElement(i)->SolveAlphaMatrix();
			}
		}

		//output solution
		printf_s("Print the mech result into .dat file... ");
		//without deformations
		if (is_print_result) {
			FILE* fout_tech;
			char name_u_tech[5000];
			sprintf_s(name_u_tech, "%s/U_non_deformation.dat", result_directory);
			fopen_s(&fout_tech, name_u_tech, "w");
			char name_in_file[1000];
			sprintf_s(name_in_file, "Elastic_T=%.2e", TIME_curr);
			std::vector<std::vector<char>> name_value(6);
			char name_v_tmp[6][100];
			sprintf_s(name_v_tmp[0], "Mises_stress");
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

				Point<double> U = solver_grid.GetSolutionInPoint(i, Centr, U_curr);
				Point<Point<double>> dU = solver_grid.GetDerevativeFromSolutionInPoint(i, Centr, U_curr);

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

				double mises_stress = sqrt((sigma[0] - sigma[1]) * (sigma[0] - sigma[1])
					+ (sigma[1] - sigma[2]) * (sigma[1] - sigma[2])
					+ (sigma[0] - sigma[2]) * (sigma[0] - sigma[2])
					+ 6 * (sigma[3] * sigma[3] + sigma[4] * sigma[4] + sigma[5] * sigma[5]) / 2.0);

				value[0][i] = mises_stress;
				value[1][i] = sigma[1];
				value[2][i] = sigma[2];


				value[3][i] = U.x;
				value[4][i] = U.y;
				value[5][i] = U.z;

				if (sigma_inv > sigma_inv_max)
				{
					sigma_inv_max = sigma_inv;
					elem_sigma_max = i;
				}
			}
			/*value[3].resize(solver_grid.GetVertexCount());
			value[4].resize(solver_grid.GetVertexCount());
			value[5].resize(solver_grid.GetVertexCount());
			for (int i = 0; i < value[3].size(); i++)
			{
				value[3][i] = U_curr[i].x;
				value[4][i] = U_curr[i].y;
				value[5][i] = U_curr[i].z;

			}*/
			//solver_grid.MoveCoordinates(U_curr);

			for (int i = 0; i < solver_grid.GetVertexCount(); i++)
			{
				solver_grid.MoveTheVertex(i, U_curr[i]);
			}

			solver_grid.printTecPlot3D_DiffDomains(fout_tech, value, name_value, name_in_file);

			for (int i = 0; i < solver_grid.GetVertexCount(); i++)
			{
				solver_grid.MoveTheVertex(i, U_curr[i] * (-1));
			}

			fclose(fout_tech);
		}
		//via YZ plane
		if (is_print_result) {
			math::SimpleGrid grid_YZ_plane; //input
			char in_file[1000];
			sprintf_s(in_file, "%s/Plane_YZ.dat", mesh_directory);
			grid_YZ_plane.ReadFromSalomeDat(in_file, 2);

			FILE* fout_tech;
			char name_u_tech[5000];
			sprintf_s(name_u_tech, "%s/U_plane_YZ.dat", result_directory);
			fopen_s(&fout_tech, name_u_tech, "w");
			char name_in_file[1000];
			sprintf_s(name_in_file, "Time_%.4e_YZ", TIME_curr);
			std::vector<std::vector<char>> name_value(6);
			char name_v_tmp[6][100];
			sprintf_s(name_v_tmp[0], "Mises_stress");
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
			value[0].resize(grid_YZ_plane.nvtr.size());
			value[1].resize(grid_YZ_plane.nvtr.size());
			value[2].resize(grid_YZ_plane.nvtr.size());
			value[3].resize(grid_YZ_plane.nvtr.size());
			value[4].resize(grid_YZ_plane.nvtr.size());
			value[5].resize(grid_YZ_plane.nvtr.size());
			double sigma_inv_max = 0;
			int elem_sigma_max = 0;
			for (int i = 0; i < grid_YZ_plane.nvtr.size(); i++)
			{
				Point<double> Centr;
				for (int j = 0; j < grid_YZ_plane.nvtr[i].size(); j++)
				{
					Centr += grid_YZ_plane.xyz[grid_YZ_plane.nvtr[i][j]];
				}
				Centr /= grid_YZ_plane.nvtr[i].size();

				double len;
				int id_elem = solver_grid.GetNearestElementID(Centr, len);
				if (id_elem >= 0)
				{
					auto element = solver_grid.GetElement(id_elem);

					Point<double> U = solver_grid.GetSolutionInPoint(id_elem, Centr, U_curr);
					Point<Point<double>> dU = solver_grid.GetDerevativeFromSolutionInPoint(id_elem, Centr, U_curr);

					double eps[6] = { dU.x.x, dU.y.y, dU.z.z, dU.x.y + dU.y.x, dU.y.z + dU.z.y, dU.x.z + dU.z.x };
					double v = solver_grid.GetDomain(element->GetIdDomain())->forMech.v;
					double E = solver_grid.GetDomain(element->GetIdDomain())->forMech.GetE(0);
					double a = v / (1 - v);
					double b = (1 - 2 * v) / (2 * (1 - v));
					//double k = E * (1 - v) / ((1 + v)*(1 - v));
					double k = E * (1. - v) / (1. + v) / (1. - 2. * v);
					double sigma[6] = { k * (eps[0] * 1 + eps[1] * a + eps[2] * a), k * (eps[0] * a + eps[1] * 1 + eps[2] * a), k * (eps[0] * a + eps[1] * a + eps[2] * 1),
						k * (2 * eps[3] * b), k * (2 * eps[4] * b), k * (2 * eps[5] * b) };
					double eps_inv = sqrt((eps[0] - eps[1]) * (eps[0] - eps[1]) + (eps[1] - eps[2]) * (eps[1] - eps[2]) + (eps[2] - eps[0]) * (eps[2] - eps[0]) +
						3 * (eps[3] * eps[3] + eps[4] * eps[4] + eps[5] * eps[5]) / 2.) * sqrt(2.) / 3;

					double mises_stress = sqrt((sigma[0] - sigma[1]) * (sigma[0] - sigma[1])
						+ (sigma[1] - sigma[2]) * (sigma[1] - sigma[2])
						+ (sigma[0] - sigma[2]) * (sigma[0] - sigma[2])
						+ 6 * (sigma[3] * sigma[3] + sigma[4] * sigma[4] + sigma[5] * sigma[5]) / 2.0);

					value[0][i] = mises_stress;
					value[1][i] = sigma[1];
					value[2][i] = sigma[2];

					value[3][i] = U.x;
					value[4][i] = U.y;
					value[5][i] = U.z;

					if (mises_stress > sigma_inv_max)
					{
						sigma_inv_max = mises_stress;
						elem_sigma_max = i;
					}
				}
			}
			/*value[3].resize(solver_grid.GetVertexCount());
			value[4].resize(solver_grid.GetVertexCount());
			value[5].resize(solver_grid.GetVertexCount());
			for (int i = 0; i < value[3].size(); i++)
			{
				value[3][i] = U_curr[i].x;
				value[4][i] = U_curr[i].y;
				value[5][i] = U_curr[i].z;

			}*/
			//solver_grid.MoveCoordinates(Solution);

			grid_YZ_plane.printTecPlot3D(fout_tech, value, name_value, name_in_file);
			fclose(fout_tech);
		}
		printf_s("\t complite\n");

		//solver_grid.~Grid_forMech();
	}
}

void Solve_ElasticDeformationProblem_MSH(char properties_file[1000])
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

					for(int i = 1; i < boundary_faces[id_boundaries[id_type]][id_triang].size(); i++)
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
				if(is_use.x) fprintf_s(fout_volume, "%.2e, ", value.x);
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

void Solve_ViscoElasticDeformationProblem()
{
	FEM::Grid_forMech solver_grid; //output
	std::vector<Point<double>> Solution; //output
	CSSD_Matrix<Tensor2Rank3D, Point<double>> global_SLAE;
	
	char properties_file[1000] = { "E:/Box/200x200x200/200k_fem_nonlin/param_for_solver.txt" };
	//char properties_file[1000] = { "E:/Box/200x200x200/300el/param.txt" };
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

	FILE *f_properties;
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
	geo_grid.ReadFromNVTR(mesh_directory, 4);

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
		std::function<Point<double>(Point<bool>&, int)> value;
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

			first_boundaries[i].value = [value, is_condition](Point<bool> &is_take, int id) -> Point<double>
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
	//zero step (assembling of global matrix and obtain elastic solution)
	{
		for (int i = 0; i < _materials.size(); i++)
		{
			solver_grid.AddDomain();
			auto domain = solver_grid.GetDomain(i);
			domain->forMech.SetE(_materials[i]._E);
			domain->forMech.SetV(_materials[i]._v);
		}

		FEM::FEM_forElasticDeformation(
			is_print_logFile,
			critical_residual,
			geo_grid, //input
			first_boundaries, //input
			second_boundary, //input
			base_result_directory, //output
			solver_grid, //output
			Solution, //output
			global_SLAE //output
		);

		//output solution
		printf_s("Print the mech result into .dat file... ");
		{
			FILE *fout_tech;
			char name_u_tech[5000];
			sprintf_s(name_u_tech, "%s/U_0.dat", base_result_directory);
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
			solver_grid.printTecPlot3D(fout_tech, value, name_value, name_in_file);
			for (int id_DOF = 0; id_DOF < Solution.size(); id_DOF++)
			{
				solver_grid.MoveTheVertex(solver_grid.accordance_DOF_and_vertex[id_DOF], Solution[id_DOF] * (-1));
			}
			fclose(fout_tech);
		}
		printf_s("\t complite\n");
	}

	//iteration process for nonlinear
	double residual = 1;
	int iter = 0;
	do
	{
		iter++;

		//only the right hand is changed
		for (int id_el = 0; id_el < solver_grid.GetElementsCount(); id_el++)
		{
			if(id_el%10 == 0)
				printf("Solve Q for elem[%d/%d]\r", id_el, solver_grid.GetElementsCount());

			auto element = solver_grid.GetElement(id_el);
			element->SetIntegrationLaw(1);
			std::vector<double> W;
			std::vector<Point<double>> X;
			element->GetIntegrationProperties(W, X);

			std::vector<double> Q_new(3 * element->GetDOFsCount());
			
			for(int INTEGR = 0; INTEGR < W.size(); INTEGR++)
			{
				std::vector<std::vector<double>> B;
				math::ResizeVector(B, 6, 3 * element->GetDOFsCount());
				for (int i = 0; i < element->GetDOFsCount(); i++)
				{
					auto dphi_i = *element->GetDerivativeOfBasisFunctionInLocalID(i);
					B[0][i * 3 + 0] = dphi_i(X[INTEGR]).x.x;
					B[1][i * 3 + 1] = dphi_i(X[INTEGR]).y.y;
					B[2][i * 3 + 2] = dphi_i(X[INTEGR]).z.z;

					B[3][i * 3 + 0] = dphi_i(X[INTEGR]).x.y;
					B[3][i * 3 + 1] = dphi_i(X[INTEGR]).y.x;
					B[4][i * 3 + 1] = dphi_i(X[INTEGR]).y.z;
					B[4][i * 3 + 2] = dphi_i(X[INTEGR]).z.y;
					B[5][i * 3 + 0] = dphi_i(X[INTEGR]).x.z;
					B[5][i * 3 + 2] = dphi_i(X[INTEGR]).z.x;
				}

				std::vector<std::vector<double>> D = solver_grid.GetDomain(element->GetIdDomain())->forMech.GetD(3);

				std::vector<std::vector<double>> G_xx, G_yy, G_zz, G_xy, G_yx, G_yz, G_zy, G_xz, G_zx;
				math::ResizeVector(G_xx, 3 * element->GetDOFsCount(), 3 * element->GetDOFsCount());
				math::ResizeVector(G_yy, 3 * element->GetDOFsCount(), 3 * element->GetDOFsCount());
				math::ResizeVector(G_zz, 3 * element->GetDOFsCount(), 3 * element->GetDOFsCount());
				math::ResizeVector(G_xy, 3 * element->GetDOFsCount(), 3 * element->GetDOFsCount());
				math::ResizeVector(G_yz, 3 * element->GetDOFsCount(), 3 * element->GetDOFsCount());
				math::ResizeVector(G_xz, 3 * element->GetDOFsCount(), 3 * element->GetDOFsCount());
				math::ResizeVector(G_yx, 3 * element->GetDOFsCount(), 3 * element->GetDOFsCount());
				math::ResizeVector(G_zy, 3 * element->GetDOFsCount(), 3 * element->GetDOFsCount());
				math::ResizeVector(G_zx, 3 * element->GetDOFsCount(), 3 * element->GetDOFsCount());

				std::vector<double> U(element->GetDOFsCount() * 3);
				for (int i = 0; i < element->GetDOFsCount(); i++)
				{
					auto dphi_i = *element->GetDerivativeOfBasisFunctionInLocalID(i);
					Point<double> dphi_i_dx(dphi_i(X[INTEGR]).x.x, dphi_i(X[INTEGR]).y.x, dphi_i(X[INTEGR]).z.x);
					Point<double> dphi_i_dy(dphi_i(X[INTEGR]).x.y, dphi_i(X[INTEGR]).y.y, dphi_i(X[INTEGR]).z.y);
					Point<double> dphi_i_dz(dphi_i(X[INTEGR]).x.z, dphi_i(X[INTEGR]).y.z, dphi_i(X[INTEGR]).z.z);

					U[i * 3 + 0] = Solution[element->GetDOFInLocalID(i)].x;
					U[i * 3 + 1] = Solution[element->GetDOFInLocalID(i)].y;
					U[i * 3 + 2] = Solution[element->GetDOFInLocalID(i)].z;


					for (int j = 0; j < element->GetDOFsCount(); j++)
					{
						auto dphi_j = *element->GetDerivativeOfBasisFunctionInLocalID(j);
						Point<double> dphi_j_dx(dphi_j(X[INTEGR]).x.x, dphi_j(X[INTEGR]).y.x, dphi_j(X[INTEGR]).z.x);
						Point<double> dphi_j_dy(dphi_j(X[INTEGR]).x.y, dphi_j(X[INTEGR]).y.y, dphi_j(X[INTEGR]).z.y);
						Point<double> dphi_j_dz(dphi_j(X[INTEGR]).x.z, dphi_j(X[INTEGR]).y.z, dphi_j(X[INTEGR]).z.z);
						
						auto Solve_G_pk = [](int i,int j,Point<double> dphi_i_dp, Point<double> dphi_j_dk, std::vector<std::vector<double>> &G_pk) -> void
						{
							std::vector<std::vector<double>> G_ij_pk;
							math::GetTensorMult(dphi_i_dp, dphi_j_dk, G_ij_pk);
							for (int ii = 0; ii < G_ij_pk.size(); ii++)
							{
								for (int jj = 0; jj < G_ij_pk[ii].size(); jj++)
								{
									G_pk[i * 3 + ii][j * 3 + jj] = G_ij_pk[ii][jj];
								}
							}
							return;
						};

						Solve_G_pk(i, j, dphi_i_dx, dphi_j_dx, G_xx);
						Solve_G_pk(i, j, dphi_i_dy, dphi_j_dy, G_yy);
						Solve_G_pk(i, j, dphi_i_dz, dphi_j_dz, G_zz);
						Solve_G_pk(i, j, dphi_i_dx, dphi_j_dy, G_xy);
						Solve_G_pk(i, j, dphi_i_dy, dphi_j_dz, G_yz);
						Solve_G_pk(i, j, dphi_i_dx, dphi_j_dz, G_xz);

						Solve_G_pk(i, j, dphi_i_dy, dphi_j_dx, G_yx);
						Solve_G_pk(i, j, dphi_i_dz, dphi_j_dy, G_zy);
						Solve_G_pk(i, j, dphi_i_dz, dphi_j_dx, G_zx);
					}
				}
				std::vector<double> eps_nonlin(6);
				auto eps_nonlin_pk = [](std::vector<std::vector<double>> &G_pk, std::vector<double> &U) -> double
				{
					double res = 0;
					for (int ii = 0; ii < G_pk.size(); ii++)
					{
						for (int jj = 0; jj < G_pk[0].size(); jj++)
						{
							res += U[ii] * G_pk[ii][jj] * U[jj];
						}
					}
					res /= 2.;
					return res;
				};
				eps_nonlin[0] = eps_nonlin_pk(G_xx, U);
				eps_nonlin[1] = eps_nonlin_pk(G_yy, U);
				eps_nonlin[2] = eps_nonlin_pk(G_zz, U);
				eps_nonlin[3] = 2*eps_nonlin_pk(G_yx, U);
				eps_nonlin[4] = 2*eps_nonlin_pk(G_yz, U);
				eps_nonlin[5] = 2*eps_nonlin_pk(G_xz, U);

				std::vector<std::vector<double>> de_du(6);
				math::ResizeVector(de_du, 6, 3 * element->GetDOFsCount());
				/*math::MultiplicationMatrixVector(G_xx, U, de_du[0]);
				math::MultiplicationMatrixVector(G_yy, U, de_du[1]);
				math::MultiplicationMatrixVector(G_zz, U, de_du[2]);
				math::MultiplicationMatrixVector(G_xy, U, de_du[3]);
				math::MultiplicationVector(de_du[3], 2, de_du[3]);
				math::MultiplicationMatrixVector(G_yz, U, de_du[4]);
				math::MultiplicationVector(de_du[4], 2, de_du[4]);
				math::MultiplicationMatrixVector(G_xz, U, de_du[5]);
				math::MultiplicationVector(de_du[5], 2, de_du[5]);*/
				for (int ii = 0; ii < G_xx.size(); ii++)
				{
					for (int jj = 0; jj < G_xx[0].size(); jj++)
					{
						de_du[0][ii] += G_xx[ii][jj] * U[jj];
						de_du[1][ii] += G_yy[ii][jj] * U[jj];
						de_du[2][ii] += G_zz[ii][jj] * U[jj];
						de_du[3][ii] += (G_xy[ii][jj] + G_yx[ii][jj]) * U[jj];
						de_du[4][ii] += (G_yz[ii][jj] + G_zy[ii][jj]) * U[jj];
						de_du[5][ii] += (G_xz[ii][jj] + G_zx[ii][jj]) * U[jj];
					}
				}


				// 1 integral
				std::vector<std::vector<double>> BT_D;
				math::MultMatrixTMatrix(B, D, BT_D);
				std::vector<double> I1;
				math::MultiplicationMatrixVector(BT_D, eps_nonlin, I1);

				// 2 integral
				std::vector<std::vector<double>> deT_du_D;
				math::MultMatrixTMatrix(de_du, D, deT_du_D);
				std::vector<std::vector<double>> deT_du_D_B;
				math::MultMatrixMatrix(deT_du_D, B, deT_du_D_B);
				std::vector<double> I2;
				math::MultiplicationMatrixVector(deT_du_D_B, U, I2);

				
				// 3 integral
				std::vector<double> I3;
				math::MultiplicationMatrixVector(deT_du_D, eps_nonlin, I3);

				for (int ii = 0; ii < Q_new.size(); ii++)
				{
					Q_new[ii] += W[INTEGR] * (I1[ii] + I2[ii] + I3[ii]);
				}

			}

			for (int _id_DOF = 0; _id_DOF < element->GetDOFsCount(); _id_DOF++)
			{
				int global_DOF = element->GetDOFInLocalID(_id_DOF);
				global_SLAE.F[global_DOF].x -= Q_new[_id_DOF * 3 + 0];
				global_SLAE.F[global_DOF].y -= Q_new[_id_DOF * 3 + 1];
				global_SLAE.F[global_DOF].z -= Q_new[_id_DOF * 3 + 2];
			}
			

			/*for (int _id_DOF = 0; _id_DOF < element->GetDOFsCount(); _id_DOF++)
			{
				std::function<Point<double>(Point<double>)> IntegralSumm_for_I = [&](Point<double> X) -> Point<double>
				{
					Point<double> result;
					auto gradS_phi = [&](int i, Point<double> x, std::vector<std::vector<double>> &gradS_phi_i) ->void
					{
						auto dphi_i = *element->GetDerivativeOfBasisFunctionInLocalID(i);

						gradS_phi_i[0][0] = dphi_i(x).x.x;
						gradS_phi_i[1][1] = dphi_i(x).y.y;
						gradS_phi_i[2][2] = dphi_i(x).z.z;

						if (gradS_phi_i.size() == 6)
						{
							gradS_phi_i[3][0] = gradS_phi_i[1][1];
							gradS_phi_i[5][0] = gradS_phi_i[2][2];

							gradS_phi_i[3][1] = gradS_phi_i[0][0];
							gradS_phi_i[4][1] = gradS_phi_i[2][2];

							gradS_phi_i[4][2] = gradS_phi_i[1][1];
							gradS_phi_i[5][2] = gradS_phi_i[0][0];
						}
						else {
							gradS_phi_i[0][3] = gradS_phi_i[1][1];
							gradS_phi_i[0][5] = gradS_phi_i[2][2];

							gradS_phi_i[1][3] = gradS_phi_i[0][0];
							gradS_phi_i[1][4] = gradS_phi_i[2][2];

							gradS_phi_i[2][4] = gradS_phi_i[1][1];
							gradS_phi_i[2][5] = gradS_phi_i[0][0];
						}
						return;
					};
					std::vector<std::vector<double>> gradST_phi_i, gradS_phi_i;
					math::ResizeVector(gradST_phi_i, 3, 6);
					math::ResizeVector(gradS_phi_i, 6, 3);
					gradS_phi(_id_DOF, X, gradST_phi_i);
					gradS_phi(_id_DOF, X, gradS_phi_i);

					Point<double> U_inX = solver_grid.GetSolutionInPoint(id_el, X, Solution);
					Tensor2Rank3D e_tens = solver_grid.GetStrainTensorFromSolutionInPoint(id_el, X, Solution);
					std::vector<double> eps_nonlin = solver_grid.GetDomain(element->GetIdDomain())->forMech.transfer_tensor_into_vector_EPS(e_tens);

					std::vector<std::vector<double>> D = solver_grid.GetDomain(element->GetIdDomain())->forMech.GetD(3);

					// 1 integral
					std::vector<std::vector<double>> dphi_D;
					math::MultMatrixMatrix(gradST_phi_i, D, dphi_D);
					std::vector<double> I1___(3);
					math::MultiplicationMatrixVector(dphi_D, eps_nonlin, I1___);

					// 2 integral
					std::vector<std::vector<double>> deT_du;
					math::ResizeVector(deT_du, 3 * element->GetDOFsCount(), 6);
					std::vector<std::vector<double>> GIj_pk, GIj_kp;//p,k = {x,y,z}
					math::ResizeVector(GIj_pk, 3, 3);
					math::ResizeVector(GIj_kp, 3, 3);

					auto dphi_I = *element->GetDerivativeOfBasisFunctionInLocalID(_id_DOF);
					Point<double> dphi_I_dx(dphi_I(X).x.x, dphi_I(X).y.x, dphi_I(X).z.x);
					Point<double> dphi_I_dy(dphi_I(X).x.y, dphi_I(X).y.y, dphi_I(X).z.y);
					Point<double> dphi_I_dz(dphi_I(X).x.z, dphi_I(X).y.z, dphi_I(X).z.z);
					for (int j = 0; false && j < element->GetDOFsCount(); j++)
					{
						Point<double> u_j = Solution[element->GetDOFInLocalID(j)];
						Point<double> GIj_pk_u_j;
						auto dphi_j = *element->GetDerivativeOfBasisFunctionInLocalID(j);
						Point<double> dphi_j_dx(dphi_j(X).x.x, dphi_j(X).y.x, dphi_j(X).z.x);
						Point<double> dphi_j_dy(dphi_j(X).x.y, dphi_j(X).y.y, dphi_j(X).z.y);
						Point<double> dphi_j_dz(dphi_j(X).x.z, dphi_j(X).y.z, dphi_j(X).z.z);

						//G_Ij_xx;
						math::GetTensorMult(dphi_I_dx, dphi_j_dx, GIj_pk);
						GIj_pk_u_j.x = GIj_pk[0][0] * u_j.x + GIj_pk[0][1] * u_j.y + GIj_pk[0][2] * u_j.z;
						GIj_pk_u_j.y = GIj_pk[1][0] * u_j.x + GIj_pk[1][1] * u_j.y + GIj_pk[1][2] * u_j.z;
						GIj_pk_u_j.z = GIj_pk[2][0] * u_j.x + GIj_pk[2][1] * u_j.y + GIj_pk[2][2] * u_j.z;
						deT_du[3 * j + 0][0] += GIj_pk_u_j.x;
						deT_du[3 * j + 1][0] += GIj_pk_u_j.y;
						deT_du[3 * j + 2][0] += GIj_pk_u_j.z;

						//G_Ij_yy;
						math::GetTensorMult(dphi_I_dy, dphi_j_dy, GIj_pk);
						GIj_pk_u_j.x = GIj_pk[0][0] * u_j.x + GIj_pk[0][1] * u_j.y + GIj_pk[0][2] * u_j.z;
						GIj_pk_u_j.y = GIj_pk[1][0] * u_j.x + GIj_pk[1][1] * u_j.y + GIj_pk[1][2] * u_j.z;
						GIj_pk_u_j.z = GIj_pk[2][0] * u_j.x + GIj_pk[2][1] * u_j.y + GIj_pk[2][2] * u_j.z;
						deT_du[3 * j + 0][1] += GIj_pk_u_j.x;
						deT_du[3 * j + 1][1] += GIj_pk_u_j.y;
						deT_du[3 * j + 2][1] += GIj_pk_u_j.z;

						//G_Ij_zz;
						math::GetTensorMult(dphi_I_dz, dphi_j_dz, GIj_pk);
						GIj_pk_u_j.x = GIj_pk[0][0] * u_j.x + GIj_pk[0][1] * u_j.y + GIj_pk[0][2] * u_j.z;
						GIj_pk_u_j.y = GIj_pk[1][0] * u_j.x + GIj_pk[1][1] * u_j.y + GIj_pk[1][2] * u_j.z;
						GIj_pk_u_j.z = GIj_pk[2][0] * u_j.x + GIj_pk[2][1] * u_j.y + GIj_pk[2][2] * u_j.z;
						deT_du[3 * j + 0][2] += GIj_pk_u_j.x;
						deT_du[3 * j + 1][2] += GIj_pk_u_j.y;
						deT_du[3 * j + 2][2] += GIj_pk_u_j.z;

						//G_Ij_xy + G_Ij_yx;
						math::GetTensorMult(dphi_I_dx, dphi_j_dy, GIj_pk);
						math::GetTensorMult(dphi_I_dy, dphi_j_dx, GIj_kp);
						GIj_pk_u_j.x = (GIj_pk[0][0] + GIj_kp[0][0]) * u_j.x + (GIj_pk[0][1] + GIj_kp[0][1]) * u_j.y + (GIj_pk[0][2] + GIj_kp[0][2]) * u_j.z;
						GIj_pk_u_j.y = (GIj_pk[1][0] + GIj_kp[1][0]) * u_j.x + (GIj_pk[1][1] + GIj_kp[1][1]) * u_j.y + (GIj_pk[1][2] + GIj_kp[1][2]) * u_j.z;
						GIj_pk_u_j.z = (GIj_pk[2][0] + GIj_kp[2][0]) * u_j.x + (GIj_pk[2][1] + GIj_kp[2][1]) * u_j.y + (GIj_pk[2][2] + GIj_kp[2][2]) * u_j.z;
						deT_du[3 * j + 0][3] += GIj_pk_u_j.x;
						deT_du[3 * j + 1][3] += GIj_pk_u_j.y;
						deT_du[3 * j + 2][3] += GIj_pk_u_j.z;

						//G_Ij_yz + G_Ij_zy;
						math::GetTensorMult(dphi_I_dy, dphi_j_dz, GIj_pk);
						math::GetTensorMult(dphi_I_dz, dphi_j_dy, GIj_kp);
						GIj_pk_u_j.x = (GIj_pk[0][0] + GIj_kp[0][0]) * u_j.x + (GIj_pk[0][1] + GIj_kp[0][1]) * u_j.y + (GIj_pk[0][2] + GIj_kp[0][2]) * u_j.z;
						GIj_pk_u_j.y = (GIj_pk[1][0] + GIj_kp[1][0]) * u_j.x + (GIj_pk[1][1] + GIj_kp[1][1]) * u_j.y + (GIj_pk[1][2] + GIj_kp[1][2]) * u_j.z;
						GIj_pk_u_j.z = (GIj_pk[2][0] + GIj_kp[2][0]) * u_j.x + (GIj_pk[2][1] + GIj_kp[2][1]) * u_j.y + (GIj_pk[2][2] + GIj_kp[2][2]) * u_j.z;
						deT_du[3 * j + 0][4] += GIj_pk_u_j.x;
						deT_du[3 * j + 1][4] += GIj_pk_u_j.y;
						deT_du[3 * j + 2][4] += GIj_pk_u_j.z;

						//G_Ij_xz + G_Ij_zx;
						math::GetTensorMult(dphi_I_dx, dphi_j_dz, GIj_pk);
						math::GetTensorMult(dphi_I_dz, dphi_j_dx, GIj_kp);
						GIj_pk_u_j.x = (GIj_pk[0][0] + GIj_kp[0][0]) * u_j.x + (GIj_pk[0][1] + GIj_kp[0][1]) * u_j.y + (GIj_pk[0][2] + GIj_kp[0][2]) * u_j.z;
						GIj_pk_u_j.y = (GIj_pk[1][0] + GIj_kp[1][0]) * u_j.x + (GIj_pk[1][1] + GIj_kp[1][1]) * u_j.y + (GIj_pk[1][2] + GIj_kp[1][2]) * u_j.z;
						GIj_pk_u_j.z = (GIj_pk[2][0] + GIj_kp[2][0]) * u_j.x + (GIj_pk[2][1] + GIj_kp[2][1]) * u_j.y + (GIj_pk[2][2] + GIj_kp[2][2]) * u_j.z;
						deT_du[3 * j + 0][5] += GIj_pk_u_j.x;
						deT_du[3 * j + 1][5] += GIj_pk_u_j.y;
						deT_du[3 * j + 2][5] += GIj_pk_u_j.z;
					}

					std::vector<std::vector<double>> deT_du_D, deT_du_D_gradS_phi_i;
					math::MultMatrixMatrix(deT_du, D, deT_du_D);
					math::MultMatrixMatrix(deT_du_D, gradS_phi_i, deT_du_D_gradS_phi_i);
					std::vector<double> I2___(3);
					I2___[0] = deT_du_D_gradS_phi_i[0][0] * U_inX.x + deT_du_D_gradS_phi_i[0][1] * U_inX.y + deT_du_D_gradS_phi_i[0][2] * U_inX.z;
					I2___[1] = deT_du_D_gradS_phi_i[1][0] * U_inX.x + deT_du_D_gradS_phi_i[1][1] * U_inX.y + deT_du_D_gradS_phi_i[1][2] * U_inX.z;
					I2___[2] = deT_du_D_gradS_phi_i[2][0] * U_inX.x + deT_du_D_gradS_phi_i[2][1] * U_inX.y + deT_du_D_gradS_phi_i[2][2] * U_inX.z;

					// 3 integral
					std::vector<double> I3___(3);
					math::MultiplicationMatrixVector(deT_du_D, eps_nonlin, I3___);

					result.x = I1___[0] + I2___[0] + I3___[0];
					result.y = I1___[1] + I2___[1] + I3___[1];
					result.z = I1___[2] + I2___[2] + I3___[2];

					return result;
				};
				
				Point<double> integr_res = element->SolveIntegral(IntegralSumm_for_I);

				int global_id = element->GetDOFInLocalID(_id_DOF);
				global_SLAE.F[global_id] -= integr_res;

			}*/
		}

		//first condition
		printf("First boundary conditions...");
		for (int id_vertex = 0; id_vertex < solver_grid.boundary_vertexes.size(); id_vertex++)
		{
			Point<double> X, F;
			auto boundary = &(solver_grid.boundary_vertexes[id_vertex]);

			Point<bool> is_take;
			Point<double> boundary_value = boundary->boundary_value(is_take, id_vertex);
			int global_id = boundary->GetDOFInLocalID(0);

			auto enter_boundary = [&](int position) {
				switch (position)
				{
				case 0:
					global_SLAE.X[global_id].x = boundary_value.x;
					global_SLAE.F[global_id].x = boundary_value.x;
					if (global_SLAE.F[global_id].x != global_SLAE.F[global_id].x)
					{
						printf_s("F[%d] = %.2e\n", global_id, global_SLAE.F[global_id].x);
					}
					break;
				case 1:
					global_SLAE.X[global_id].y = boundary_value.y;
					global_SLAE.F[global_id].y = boundary_value.y;
					if (global_SLAE.F[global_id].y != global_SLAE.F[global_id].y)
					{
						printf_s("F[%d] = %.2e\n", global_id, global_SLAE.F[global_id].y);
					}
					break;
				case 2:
					global_SLAE.X[global_id].z = boundary_value.z;
					global_SLAE.F[global_id].z = boundary_value.z;
					if (global_SLAE.F[global_id].z != global_SLAE.F[global_id].z)
					{
						printf_s("F[%d] = %.2e\n", global_id, global_SLAE.F[global_id].z);
					}
					break;
				default:
					break;
				}
			};
			if (is_take.x == true) enter_boundary(0);
			if (is_take.y == true) enter_boundary(1);
			if (is_take.z == true) enter_boundary(2);
		}

		//obtain the new solution
		if (true)
		{
			printf("Soluting SLAY... (%d)\n", global_SLAE.GetMatrixSize());
			math::MakeCopyVector_A_into_B(Solution, global_SLAE.X);
			//for (int i = 0; i < global_SLAE.X.size(); i++)
			//	global_SLAE.X[i] = 0;

			CSSD_Matrix<double, double> newSLAE;
			//CSSD_Matrix<Tensor2Rank3D, Point<double>> newSLAE;
			math::MakeCopyMatrix_A_into_B(global_SLAE, newSLAE);
			newSLAE.MakeDivideOnSqrtFormSumm();

			int MaxSize = global_SLAE.GetMatrixSize();
			std::vector<double> best_solution;
			math::MakeCopyVector_A_into_B(newSLAE.X, best_solution);
			double current_residual = 1, best_residual = 1e+25;
			int MAX_STEPS = 12;
			for (int i = 0; i <= MAX_STEPS; i++)
			{
				current_residual = pow(10., -1 * (i + 1));
				//current_residual /= 2.;
				printf_s("//---> I = %d/%d (full size %d) - needed residual %.2e\n", i, MAX_STEPS, newSLAE.GetMatrixSize(), current_residual);
				current_residual = abs(newSLAE.BiCG_Stab(MaxSize, current_residual));
				if (current_residual <= critical_residual)
					break;
				//if (current_residual > best_residual + (1e-10))
				//{
				//	math::MakeCopyVector_A_into_B(best_solution, newSLAE.X);
				//	printf_s("//---> BEST residual %.2e\n", best_residual);
				//	break;
				//}
				if (current_residual < best_residual)
				{
					best_residual = current_residual;
					math::MakeCopyVector_A_into_B(newSLAE.X, best_solution);
				}
			}
			math::MakeCopyVector_A_into_B(newSLAE.X, global_SLAE.X);
		}

		double max_prev = 0, max_curr = 0;
		for (int i = 0; i < Solution.size(); i++)
		{
			if (abs(Solution[i].x) > max_prev) max_prev = abs(Solution[i].x);
			if (abs(Solution[i].y) > max_prev) max_prev = abs(Solution[i].y);
			if (abs(Solution[i].z) > max_prev) max_prev = abs(Solution[i].z);

			if (abs(global_SLAE.X[i].x) > max_curr) max_curr = abs(global_SLAE.X[i].x);
			if (abs(global_SLAE.X[i].y) > max_curr) max_curr = abs(global_SLAE.X[i].y);
			if (abs(global_SLAE.X[i].z) > max_curr) max_curr = abs(global_SLAE.X[i].z);
		}
		residual = abs(max_curr - max_prev) / max_curr;
		math::MakeCopyVector_A_into_B(global_SLAE.X, Solution);
		printf("\n---------------------\nIteration[%d]: Residual = %.4e\n=======================\n", iter, residual);

		//output solution
		printf_s("Print the mech result into .dat file... ");
		{
			FILE *fout_tech;
			char name_u_tech[5000];
			sprintf_s(name_u_tech, "%s/U_iter%d.dat", base_result_directory, iter);
			fopen_s(&fout_tech, name_u_tech, "w");
			char name_in_file[1000];
			sprintf_s(name_in_file, "ViscoElastic_iter%d", iter);
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
			solver_grid.printTecPlot3D(fout_tech, value, name_value, name_in_file);
			for (int id_DOF = 0; id_DOF < Solution.size(); id_DOF++)
			{
				solver_grid.MoveTheVertex(solver_grid.accordance_DOF_and_vertex[id_DOF], Solution[id_DOF] * (-1));
			}
			fclose(fout_tech);
		}
		printf_s("\t complite\n\n");

	} while (residual > 1e-4);

	solver_grid.~Grid_forMech();
}
void Solve_ViscoElasticDeformationProblem_viaNewton()
{
	FEM::Grid_forMech solver_grid; //output
	std::vector<Point<double>> Solution; //output
	CSSD_Matrix<Tensor2Rank3D, Point<double>> global_SLAE;

	//char properties_file[1000] = { "E:/Box/200x200x200/200k_fem_nonlin/param_for_solver.txt" };
	char properties_file[1000] = { "E:/Box/beam/1x0.2x0.02/middle_up_force/33k/param_for_solver.txt" };
	//char properties_file[1000] = { "E:/Box/beam/1x0.2x0.02/middle_up_force/295k/param_for_solver.txt" };
	//char properties_file[1000] = { "E:/Box/200x200x200/300el/param.txt" };
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

	FILE *f_properties;
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
	geo_grid.ReadFromNVTR(mesh_directory, 4);

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
		std::function<Point<double>(Point<bool>&, int)> value;
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

			first_boundaries[i].value = [value, is_condition](Point<bool> &is_take, int id) -> Point<double>
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
	//zero step (assembling of global matrix and obtain elastic solution)
	{
		for (int i = 0; i < _materials.size(); i++)
		{
			solver_grid.AddDomain();
			auto domain = solver_grid.GetDomain(i);
			domain->forMech.SetE(_materials[i]._E);
			domain->forMech.SetV(_materials[i]._v);
		}

		FEM::FEM_forElasticDeformation(
			is_print_logFile,
			critical_residual,
			geo_grid, //input
			first_boundaries, //input
			second_boundary, //input
			base_result_directory, //output
			solver_grid, //output
			Solution, //output
			global_SLAE //output
		);

		//output solution
		printf_s("Print the mech result into .dat file... ");
		{
			FILE *fout_tech;
			char name_u_tech[5000];
			sprintf_s(name_u_tech, "%s/U_0.dat", base_result_directory);
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
			solver_grid.printTecPlot3D(fout_tech, value, name_value, name_in_file);
			for (int id_DOF = 0; id_DOF < Solution.size(); id_DOF++)
			{
				solver_grid.MoveTheVertex(solver_grid.accordance_DOF_and_vertex[id_DOF], Solution[id_DOF] * (-1));
			}
			fclose(fout_tech);
		}
		printf_s("\t complite\n");
	}

	//iteration process for nonlinear
	double residual = 1;
	int iter = 0;
	do
	{
		iter++;

		math::InitializationVector(global_SLAE.A_down, 0.0);
		math::InitializationVector(global_SLAE.A_up, 0.0);
		math::InitializationVector(global_SLAE.Diag, 0.0);
		math::InitializationVector(global_SLAE.F, 0.0);
		math::InitializationVector(global_SLAE.X, 0.0);

		//only the additional part is changed
		for (int id_el = 0; id_el < solver_grid.GetElementsCount(); id_el++)
		{
			if (id_el % 10 == 0)
				printf("Solve matrix for elem[%d/%d]\r", id_el, solver_grid.GetElementsCount());

			auto element = solver_grid.GetElement(id_el);
			element->SetIntegrationLaw(1);
			std::vector<double> W;
			std::vector<Point<double>> X;
			element->GetIntegrationProperties(W, X);

			std::vector<std::vector<double>> K;
			std::vector<double> b(3 * element->GetDOFsCount());
			math::ResizeVector(K, 3 * element->GetDOFsCount(), 3 * element->GetDOFsCount());

			for (int INTEGR = 0; INTEGR < W.size(); INTEGR++)
			{
				std::vector<std::vector<double>> D = solver_grid.GetDomain(element->GetIdDomain())->forMech.GetD(3);

				std::vector<std::vector<double>> B0, BL, B, A, G, K_0_L, M, K_s;
				math::ResizeVector(B0, 6, 3 * element->GetDOFsCount());
				math::ResizeVector(BL, 6, 3 * element->GetDOFsCount());
				math::ResizeVector(B, 6, 3 * element->GetDOFsCount());
				math::ResizeVector(A, 6, 9);
				math::ResizeVector(M, 9, 9);
				math::ResizeVector(G, 9, 3 * element->GetDOFsCount());
				std::vector<double> U_prev(element->GetDOFsCount() * 3);
				Point<double> Qx, Qy, Qz;
				std::vector<double> e(6), sigma(6);
				for (int i = 0; i < element->GetDOFsCount(); i++)
				{
					auto dphi_i = *element->GetDerivativeOfBasisFunctionInLocalID(i);
					Point<double> dphi_i_dx(dphi_i(X[INTEGR]).x.x, dphi_i(X[INTEGR]).y.x, dphi_i(X[INTEGR]).z.x);
					Point<double> dphi_i_dy(dphi_i(X[INTEGR]).x.y, dphi_i(X[INTEGR]).y.y, dphi_i(X[INTEGR]).z.y);
					Point<double> dphi_i_dz(dphi_i(X[INTEGR]).x.z, dphi_i(X[INTEGR]).y.z, dphi_i(X[INTEGR]).z.z);
					U_prev[i * 3 + 0] = Solution[element->GetDOFInLocalID(i)].x;
					U_prev[i * 3 + 1] = Solution[element->GetDOFInLocalID(i)].y;
					U_prev[i * 3 + 2] = Solution[element->GetDOFInLocalID(i)].z;
					//solve matr Qx, Qy, Qz
					{
						Qx.x += U_prev[i * 3 + 0] * dphi_i_dx.x;
						Qx.y += U_prev[i * 3 + 1] * dphi_i_dx.y;
						Qx.z += U_prev[i * 3 + 2] * dphi_i_dx.z;
						Qy.x += U_prev[i * 3 + 0] * dphi_i_dy.x;
						Qy.y += U_prev[i * 3 + 1] * dphi_i_dy.y;
						Qy.z += U_prev[i * 3 + 2] * dphi_i_dy.z;
						Qz.x += U_prev[i * 3 + 0] * dphi_i_dz.x;
						Qz.y += U_prev[i * 3 + 1] * dphi_i_dz.y;
						Qz.z += U_prev[i * 3 + 2] * dphi_i_dz.z;

					}
					//solve matr G[9x12]
					{
						G[0 * 3 + 0][i * 3 + 0] = dphi_i_dx.x;
						G[0 * 3 + 1][i * 3 + 1] = dphi_i_dx.y;
						G[0 * 3 + 2][i * 3 + 2] = dphi_i_dx.z;

						G[1 * 3 + 0][i * 3 + 0] = dphi_i_dy.x;
						G[1 * 3 + 1][i * 3 + 1] = dphi_i_dy.y;
						G[1 * 3 + 2][i * 3 + 2] = dphi_i_dy.z;

						G[2 * 3 + 0][i * 3 + 0] = dphi_i_dz.x;
						G[2 * 3 + 1][i * 3 + 1] = dphi_i_dz.y;
						G[2 * 3 + 2][i * 3 + 2] = dphi_i_dz.z;
					}

					B0[0][i * 3 + 0] = dphi_i(X[INTEGR]).x.x;
					B0[1][i * 3 + 1] = dphi_i(X[INTEGR]).y.y;
					B0[2][i * 3 + 2] = dphi_i(X[INTEGR]).z.z;
					B0[3][i * 3 + 0] = dphi_i(X[INTEGR]).x.y;
					B0[3][i * 3 + 1] = dphi_i(X[INTEGR]).y.x;
					B0[4][i * 3 + 1] = dphi_i(X[INTEGR]).y.z;
					B0[4][i * 3 + 2] = dphi_i(X[INTEGR]).z.y;
					B0[5][i * 3 + 0] = dphi_i(X[INTEGR]).x.z;
					B0[5][i * 3 + 2] = dphi_i(X[INTEGR]).z.x;
				}
				//solve matr A[6x9]
				{
					A[0][0 * 3 + 0] = Qx.x;
					A[0][0 * 3 + 1] = Qx.y;
					A[0][0 * 3 + 2] = Qx.z;
					A[1][1 * 3 + 0] = Qy.x;
					A[1][1 * 3 + 1] = Qy.y;
					A[1][1 * 3 + 2] = Qy.z;
					A[2][2 * 3 + 0] = Qz.x;
					A[2][2 * 3 + 1] = Qz.y;
					A[2][2 * 3 + 2] = Qz.z;
					A[3][1 * 3 + 0] = Qz.x;
					A[3][1 * 3 + 1] = Qz.y;
					A[3][1 * 3 + 2] = Qz.z;
					A[3][2 * 3 + 0] = Qy.x;
					A[3][2 * 3 + 1] = Qy.y;
					A[3][2 * 3 + 2] = Qy.z;
					A[4][0 * 3 + 0] = Qz.x;
					A[4][0 * 3 + 1] = Qz.y;
					A[4][0 * 3 + 2] = Qz.z;
					A[4][2 * 3 + 0] = Qx.x;
					A[4][2 * 3 + 1] = Qx.y;
					A[4][2 * 3 + 2] = Qx.z;
					A[5][0 * 3 + 0] = Qy.x;
					A[5][0 * 3 + 1] = Qy.y;
					A[5][0 * 3 + 2] = Qy.z;
					A[5][1 * 3 + 0] = Qx.x;
					A[5][1 * 3 + 1] = Qx.y;
					A[5][1 * 3 + 2] = Qx.z;
				}
				math::MultMatrixMatrix(A, G, BL);
				math::MakeSummMatrixMatrix(B0, BL, B);
				std::vector<std::vector<double>> Bt_D;
				math::MultMatrixTMatrix(B, D, Bt_D);
				math::MultMatrixMatrix(Bt_D, B, K_0_L);

				//solve e=e0+eL, sigma = D*e, M[9x9]
				{
					e[0] = Qx.x + (Qx.x*Qx.x + Qx.y*Qx.y + Qx.z * Qx.z) / 2.;
					e[1] = Qy.y + (Qy.x*Qy.x + Qy.y*Qy.y + Qy.z * Qy.z) / 2.;
					e[2] = Qz.z + (Qz.x*Qz.x + Qz.y*Qz.y + Qz.z * Qz.z) / 2.;
					e[3] = Qz.y + Qy.z + (Qx.x*Qy.x + Qx.y*Qy.y + Qx.z * Qy.z);
					e[4] = Qz.x + Qx.z + (Qz.x*Qy.x + Qz.y*Qy.y + Qz.z * Qy.z);
					e[5] = Qx.y + Qy.y + (Qx.x*Qz.x + Qx.y*Qz.y + Qx.z * Qz.z);
					math::MultiplicationMatrixVector(D, e, sigma);

					for (int k = 0; k < 3; k++)
					{
						M[0 * 3 + k][0 * 3 + k] = sigma[0];
						M[1 * 3 + k][1 * 3 + k] = sigma[1];
						M[2 * 3 + k][2 * 3 + k] = sigma[2];
						M[0 * 3 + k][1 * 3 + k] = sigma[5];
						M[1 * 3 + k][0 * 3 + k] = sigma[5];
						M[0 * 3 + k][2 * 3 + k] = sigma[4];
						M[2 * 3 + k][0 * 3 + k] = sigma[4];
						M[1 * 3 + k][2 * 3 + k] = sigma[3];
						M[2 * 3 + k][1 * 3 + k] = sigma[3];
					}
				}
				std::vector<std::vector<double>> Gt_M;
				math::MultMatrixTMatrix(G, M, Gt_M);
				math::MultMatrixMatrix(Gt_M, G, K_s);
				std::vector<double> Bt_sigma;
				math::MultiplicationTMatrixVector(B, sigma, Bt_sigma);

				for (int ii = 0; ii < K.size(); ii++)
				{
					for (int jj = 0; jj < K.size(); jj++)
					{
						K[ii][jj] += W[INTEGR] * (K_0_L[ii][jj] + K_s[ii][jj]);
					}
					b[ii] += W[INTEGR] * Bt_sigma[ii];
				}
			}

			DenseMatrix<Tensor2Rank3D, Point<double>> local_matix;
			local_matix.SetSize(element->GetDOFsCount());
			for (int I = 0; I < element->GetDOFsCount(); I++)
			{
				for (int J = 0; J < element->GetDOFsCount(); J++)
				{
					Tensor2Rank3D k_IJ;
					for (int ii = 0; ii < 3; ii++)
					{
						for (int jj = 0; jj < 3; jj++)
						{
							k_IJ.val[ii][jj] = K[I * 3 + ii][J * 3 + jj];
						}
					}
					local_matix.A[I][J] = k_IJ;
				}
				local_matix.F[I].x -= b[I * 3 + 0];
				local_matix.F[I].y -= b[I * 3 + 1];
				local_matix.F[I].z -= b[I * 3 + 2];
			}
			global_SLAE.SummPartOfMatrix(local_matix, *solver_grid.GetElementDOFs(id_el));
		}

		//second condition
		for (int id_triangle = 0; id_triangle < solver_grid.boundary_faces.size(); id_triangle++)
		{
			std::vector<Point<double>> local_vector_SLAE;
			solver_grid.boundary_faces[id_triangle].SolveLocalBoundaryVector(local_vector_SLAE, solver_grid.boundary_faces[id_triangle].boundary_value);
			global_SLAE.SummPartOfVector(local_vector_SLAE, *solver_grid.boundary_faces[id_triangle].GetElementDOFs());
		};

		//first condition
		//first condition
		printf("First boundary conditions...");
		for (int id_vertex = 0; id_vertex < solver_grid.boundary_vertexes.size(); id_vertex++)
		{
			Tensor2Rank3D diag, non_diag;
			Point<double> X, F;
			diag.InitializationAsI();
			non_diag.InitializationAs0();

			auto boundary = &(solver_grid.boundary_vertexes[id_vertex]);

			Point<bool> is_take;
			Point<double> boundary_value = boundary->boundary_value(is_take, id_vertex);
			boundary_value = Point<double>(0, 0, 0);
			///-->>>>>
			/*boundary_value.x = 0;
			boundary_value.y = 0;
			boundary_value.z = 20*(*solver_grid.GetPtrCoordinateViaID(id_vertex)).z - 20*200*0.5;*/
			///-->>>>>
			int global_id = boundary->GetDOFInLocalID(0);
			//printf_s("\n%d\r", global_id);

			auto enter_boundary = [&](int position) {
				switch (position)
				{
				case 0:
					global_SLAE.X[global_id].x = boundary_value.x;
					global_SLAE.F[global_id].x = boundary_value.x;
					if (global_SLAE.F[global_id].x != global_SLAE.F[global_id].x)
					{
						printf_s("F[%d] = %.2e\n", global_id, global_SLAE.F[global_id].x);
					}
					break;
				case 1:
					global_SLAE.X[global_id].y = boundary_value.y;
					global_SLAE.F[global_id].y = boundary_value.y;
					if (global_SLAE.F[global_id].y != global_SLAE.F[global_id].y)
					{
						printf_s("F[%d] = %.2e\n", global_id, global_SLAE.F[global_id].y);
					}
					break;
				case 2:
					global_SLAE.X[global_id].z = boundary_value.z;
					global_SLAE.F[global_id].z = boundary_value.z;
					if (global_SLAE.F[global_id].z != global_SLAE.F[global_id].z)
					{
						printf_s("F[%d] = %.2e\n", global_id, global_SLAE.F[global_id].z);
					}
					break;
				default:
					break;
				}

				global_SLAE.Diag[global_id].val[position][0] = 0.0;
				global_SLAE.Diag[global_id].val[position][1] = 0.0;
				global_SLAE.Diag[global_id].val[position][2] = 0.0;
				global_SLAE.Diag[global_id].val[position][position] = 1.0;
				for (int j = 0; j < global_SLAE.A_down[global_id].size(); j++)
				{
					global_SLAE.A_down[global_id][j].val[position][0] = 0.0;
					global_SLAE.A_down[global_id][j].val[position][1] = 0.0;
					global_SLAE.A_down[global_id][j].val[position][2] = 0.0;
				}
				for (int j = 0; j < global_SLAE.A_up[global_id].size(); j++)
				{
					global_SLAE.A_up[global_id][j].val[position][0] = 0.0;
					global_SLAE.A_up[global_id][j].val[position][1] = 0.0;
					global_SLAE.A_up[global_id][j].val[position][2] = 0.0;
				}

				//simmetrisation
				//int position = 0;
				if (false)
				{
					double target_coef_in_firgt_side = 0;
					switch (position)
					{
					case 0:
						target_coef_in_firgt_side = global_SLAE.F[global_id].x;
						global_SLAE.F[global_id].y -= target_coef_in_firgt_side * global_SLAE.Diag[global_id].val[1][position];
						global_SLAE.F[global_id].z -= target_coef_in_firgt_side * global_SLAE.Diag[global_id].val[2][position];
						global_SLAE.Diag[global_id].val[0][position] = 0.0;
						global_SLAE.Diag[global_id].val[1][position] = 0.0;
						break;
					case 1:
						target_coef_in_firgt_side = global_SLAE.F[global_id].y;
						global_SLAE.F[global_id].x -= target_coef_in_firgt_side * global_SLAE.Diag[global_id].val[0][position];
						global_SLAE.F[global_id].z -= target_coef_in_firgt_side * global_SLAE.Diag[global_id].val[2][position];
						global_SLAE.Diag[global_id].val[0][position] = 0.0;
						global_SLAE.Diag[global_id].val[2][position] = 0.0;
						break;
					case 2:
						target_coef_in_firgt_side = global_SLAE.F[global_id].z;
						global_SLAE.F[global_id].x -= target_coef_in_firgt_side * global_SLAE.Diag[global_id].val[0][position];
						global_SLAE.F[global_id].y -= target_coef_in_firgt_side * global_SLAE.Diag[global_id].val[1][position];
						global_SLAE.Diag[global_id].val[0][position] = 0.0;
						global_SLAE.Diag[global_id].val[1][position] = 0.0;
						break;
					}
					for (int i = 0; i < global_id; i++)
					{
						for (int j = 0; j < global_SLAE.id_column_for_A_up[i].size() &&
							global_SLAE.id_column_for_A_up[i][j] <= global_id; j++)
						{
							if (global_id == global_SLAE.id_column_for_A_up[i][j])
							{
								global_SLAE.F[i].x -= target_coef_in_firgt_side * global_SLAE.A_up[i][j].val[0][position];
								global_SLAE.F[i].y -= target_coef_in_firgt_side * global_SLAE.A_up[i][j].val[1][position];
								global_SLAE.F[i].z -= target_coef_in_firgt_side * global_SLAE.A_up[i][j].val[2][position];
								global_SLAE.A_up[i][j].val[0][position] = 0.0;
								global_SLAE.A_up[i][j].val[1][position] = 0.0;
								global_SLAE.A_up[i][j].val[2][position] = 0.0;

								if (global_SLAE.F[i].x != global_SLAE.F[i].x || abs(global_SLAE.F[i].x) > 1e+50)
								{
									printf_s("i = %d; j = %d; A = %.2e; F[%d] = %.2e\n", i, j, global_SLAE.A_up[i][j].val[0][position], global_id, global_SLAE.F[global_id].x);
								}
								if (global_SLAE.F[i].y != global_SLAE.F[i].y || abs(global_SLAE.F[i].y) > 1e+50)
								{
									printf_s("i = %d; j = %d; A = %.2e; F[%d] = %.2e\n", i, j, global_SLAE.A_up[i][j].val[1][position], global_id, global_SLAE.F[global_id].y);
								}
								if (global_SLAE.F[i].z != global_SLAE.F[i].z || abs(global_SLAE.F[i].z) > 1e+50)
								{
									printf_s("i = %d; j = %d; A = %.2e; F[%d] = %.2e\n", i, j, global_SLAE.A_up[i][j].val[2][position], global_id, global_SLAE.F[global_id].z);
								}
							}
						}
					}
					for (int i = global_id + 1; i < global_SLAE.id_column_for_A_down.size(); i++)
					{
						for (int j = 0; j < global_SLAE.id_column_for_A_down[i].size() &&
							global_SLAE.id_column_for_A_down[i][j] <= global_id; j++)
						{
							if (global_id == global_SLAE.id_column_for_A_down[i][j])
							{
								global_SLAE.F[i].x -= target_coef_in_firgt_side * global_SLAE.A_down[i][j].val[0][position];
								global_SLAE.F[i].y -= target_coef_in_firgt_side * global_SLAE.A_down[i][j].val[1][position];
								global_SLAE.F[i].z -= target_coef_in_firgt_side * global_SLAE.A_down[i][j].val[2][position];
								global_SLAE.A_down[i][j].val[0][position] = 0.0;
								global_SLAE.A_down[i][j].val[1][position] = 0.0;
								global_SLAE.A_down[i][j].val[2][position] = 0.0;
								if (global_SLAE.F[i].x != global_SLAE.F[i].x || abs(global_SLAE.F[i].x) > 1e+50)
								{
									printf_s("i = %d; j = %d; A = %.2e; F[%d] = %.2e\n", i, j, global_SLAE.A_up[i][j].val[0][position], global_id, global_SLAE.F[global_id].x);
								}
								if (global_SLAE.F[i].y != global_SLAE.F[i].y || abs(global_SLAE.F[i].y) > 1e+50)
								{
									printf_s("i = %d; j = %d; A = %.2e; F[%d] = %.2e\n", i, j, global_SLAE.A_up[i][j].val[1][position], global_id, global_SLAE.F[global_id].y);
								}
								if (global_SLAE.F[i].z != global_SLAE.F[i].z || abs(global_SLAE.F[i].z) > 1e+50)
								{
									printf_s("i = %d; j = %d; A = %.2e; F[%d] = %.2e\n", i, j, global_SLAE.A_up[i][j].val[2][position], global_id, global_SLAE.F[global_id].z);
								}
							}
						}
					}
				}
			};
			if (is_take.x == true) enter_boundary(0);
			if (is_take.y == true) enter_boundary(1);
			if (is_take.z == true) enter_boundary(2);

		}
		
		//obtain the new additioanl solution
		if (true)
		{
			printf("Soluting SLAY... (%d)\n", global_SLAE.GetMatrixSize());
			math::MakeCopyVector_A_into_B(Solution, global_SLAE.X);
			//for (int i = 0; i < global_SLAE.X.size(); i++)
			//	global_SLAE.X[i] = 0;

			CSSD_Matrix<double, double> newSLAE;
			//CSSD_Matrix<Tensor2Rank3D, Point<double>> newSLAE;
			math::MakeCopyMatrix_A_into_B(global_SLAE, newSLAE);
			newSLAE.MakeDivideOnSqrtFormSumm();

			int MaxSize = global_SLAE.GetMatrixSize();
			std::vector<double> best_solution;
			math::MakeCopyVector_A_into_B(newSLAE.X, best_solution);
			double current_residual = 1, best_residual = 1e+25;
			int MAX_STEPS = 12;
			for (int i = 0; i <= MAX_STEPS; i++)
			{
				current_residual = pow(10., -1 * (i + 1));
				//current_residual /= 2.;
				printf_s("//---> I = %d/%d (full size %d) - needed residual %.2e\n", i, MAX_STEPS, newSLAE.GetMatrixSize(), current_residual);
				current_residual = abs(newSLAE.BiCG_Stab(MaxSize, current_residual));
				if (current_residual <= critical_residual)
					break;
				//if (current_residual > best_residual + (1e-10))
				//{
				//	math::MakeCopyVector_A_into_B(best_solution, newSLAE.X);
				//	printf_s("//---> BEST residual %.2e\n", best_residual);
				//	break;
				//}
				if (current_residual < best_residual)
				{
					best_residual = current_residual;
					math::MakeCopyVector_A_into_B(newSLAE.X, best_solution);
				}
			}
			math::MakeCopyVector_A_into_B(newSLAE.X, global_SLAE.X);
		}

		double summ_add = 0, summ_right = 0;
		for (int i = 0; i < global_SLAE.X.size(); i++)
		{
			summ_add += abs(global_SLAE.X[i].x) + abs(global_SLAE.X[i].y) + abs(global_SLAE.X[i].z);
			summ_right += abs(global_SLAE.F[i].x*global_SLAE.F[i].x) + abs(global_SLAE.F[i].y*global_SLAE.F[i].y) + abs(global_SLAE.F[i].z*global_SLAE.F[i].z);
		}
		residual = sqrt(summ_add) / global_SLAE.X.size();
		double residual_right = sqrt(summ_right) / global_SLAE.X.size();
		math::MakeSummVectors(Solution, 1, global_SLAE.X, 1, Solution);
		printf("\n---------------------\nIteration[%d]: Residual = %.4e; Residual_right = %.4e\n=======================\n", iter, residual, residual_right);

		//output solution
		printf_s("Print the mech result into .dat file... ");
		{
			FILE *fout_tech;
			char name_u_tech[5000];
			sprintf_s(name_u_tech, "%s/U_iter%d.dat", base_result_directory, iter);
			fopen_s(&fout_tech, name_u_tech, "w");
			char name_in_file[1000];
			sprintf_s(name_in_file, "ViscoElastic_iter%d", iter);
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
			solver_grid.printTecPlot3D(fout_tech, value, name_value, name_in_file);
			for (int id_DOF = 0; id_DOF < Solution.size(); id_DOF++)
			{
				solver_grid.MoveTheVertex(solver_grid.accordance_DOF_and_vertex[id_DOF], Solution[id_DOF] * (-1));
			}
			fclose(fout_tech);
		}
		printf_s("\t complite\n\n");

	} while (residual > 1e-4);

	solver_grid.~Grid_forMech();
}

void TransferSelfIntoMSH()
{
	//char msh_directory[1000] = { "D:/SeismicGroup/simpl_model/new_mesh.msh" };
	//char self_directory[1000] = { "D:/SeismicGroup/simpl_model/salome_dat" };
	char msh_directory[1000] = { "./new_mesh.msh" };
	char self_directory[1000] = { "." };
	std::vector<std::vector<std::vector<int>>> boundary_faces;
	math::SimpleGrid geo_grid; //input
	geo_grid.ReadFromNVTR(self_directory, 4);

	char boundary_file[1000];
	sprintf_s(boundary_file, sizeof(boundary_file), "%s/kraev2_faces.txt", self_directory);
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
		parent_elem0--;
		parent_elem1--;
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

	geo_grid.WriteMSH(msh_directory, 4, boundary_faces);

	int a;
	scanf_s("%d", &a);
}
void TransferSelfIntoMSH_v2208()
{
	//char msh_directory[1000] = { "D:/ForEpov/tensor_rotation/mesh_BoxHomog/new_mesh.msh" };
	//char self_directory[1000] = { "D:/ForEpov/tensor_rotation/mesh_BoxHomog" };
	char msh_directory[1000] = { "./new_mesh.msh" };
	char self_directory[1000] = { "." };
	std::vector<std::vector<std::vector<int>>> boundary_faces;
	math::SimpleGrid geo_grid; //input
	geo_grid.ReadFromNVTR(self_directory, 4);

	char boundary_file[1000];
	sprintf_s(boundary_file, sizeof(boundary_file), "%s/kraev2_faces.txt", self_directory);
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
		parent_elem0--;
		parent_elem1--;
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

	geo_grid.WriteMSH_v2208(msh_directory, 4, boundary_faces);

	printf_s("Transfer was complited\n");

	int a;
	scanf_s("%d", &a);
}

void SLAE_testing()
{
	CSSD_Matrix<double, double> A, M;
	int N = 1000;
	A.RandomMatrix(N, 1.0);

	for (int i = 0; i < A.X.size(); i++)
	{
		A.X[i] = 1e-10;
	}

	M.PrecondorSSOR(0.75, A);
	double residual;
	//residual = A.MSG_PreconditioningSSOR(N, 1e-10, M);
	residual = A.MSG(N, 1e-10);
}

int main()
{
	//SLAE_testing();
	//Solve_TermalProblem_inTime();

	Solve_ElasticDeformationProblem();

	//TransferSelfIntoMSH();
	//Solve_ElastodynamicsProblem();
	//Solve_ElastodynamicsProblem_SelfDeform();
	//Solve_ElastodynamicsProblem_SelfDeform_2Order();
	//Solve_ElastodynamicsProblem_Explicit();

	//Solve_ViscoElasticDeformationProblem();

	//Solve_ViscoElasticDeformationProblem_viaNewton();

	//TransferSelfIntoMSH();
	//TransferSelfIntoMSH_v2208();

	//char base_name[1000] = { "./param_for_solver" };
	////char base_name[1000] = {"D:/Lab1104_2019/Mesh1/param_for_solver"};
	//char properties_file[1000];
	//int I = 0;
	//sprintf_s(properties_file, sizeof(properties_file), "%s_%d.txt", base_name, I);
	//FILE* fparam;
	//fopen_s(&fparam, properties_file, "r");
	//while (fparam != NULL)
	//{
	//	fclose(fparam);
	//	Solve_ElasticDeformationProblem_MSH(properties_file);
	//	I++;
	//	sprintf_s(properties_file, sizeof(properties_file), "%s_%d.txt", base_name, I);
	//	fopen_s(&fparam, properties_file, "r");
	//}

	int a;
	scanf_s("%d", &a);
}