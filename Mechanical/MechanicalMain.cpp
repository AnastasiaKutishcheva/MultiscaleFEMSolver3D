#include "../FEM/FEM.h"
#include <Windows.h>
#include <iostream>
#include "../MsFEM/MsFEM.h"

void EffectiveElastisity_MSH(char properties_file[1000])
{
	FEM::Grid_forMech solver_grid; //output
	std::vector<std::vector<Point<double>>> Solution; //output
	std::vector<std::vector<double>> D_eff; //output

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


	double _flag;
	fscanf_s(f_properties, "%lf", &_flag);
	if (_flag > 0) is_print_logFile = true;

	math::ReadNonEmptyLine(f_properties, mesh_directory);
	math::SimpleGrid geo_grid; //input
	std::vector<std::vector<std::vector<int>>> boundary_faces;
	//geo_grid.ReadFromMSH(mesh_directory, 1, 4, 4, boundary_faces); //нумерация граней не важна, но все внешние должны быть внешними
	geo_grid.ReadFromMSH_v2208(mesh_directory, 1, 4, 4, boundary_faces); //нумерация граней не важна, но все внешние должны быть внешними

	//для таяния
	{
		if (_flag > 0 && _flag < 1)
		{
			int target_material = 2;
			int new_material = 3;

			double minZ = geo_grid.xyz[0].z, maxZ = geo_grid.xyz[0].z;
			for (int i = 0; i < geo_grid.xyz.size(); i++)
			{
				if (geo_grid.xyz[i].z < minZ) minZ = geo_grid.xyz[i].z;
				if (geo_grid.xyz[i].z > maxZ) maxZ = geo_grid.xyz[i].z;
			}
			double Zmax_curr = minZ + (maxZ - minZ) * _flag;

			for (int i = 0; i < geo_grid.nvtr.size(); i++)
			{
				if (geo_grid.nvkat[i] == target_material)
				{
					bool isLower = true;
					for (int j = 0; j < geo_grid.nvtr[i].size(); j++)
					{
						if (geo_grid.xyz[geo_grid.nvtr[i][j]].z > Zmax_curr)
						{
							isLower = false;
							break;
						}
					}
					if (isLower)
					{
						geo_grid.nvkat[i] = new_material;
					}
				}
			}
		}
	}

	double critical_residual;
	{
		char _line[1000];
		math::ReadNonEmptyLine_forNumbers(f_properties, _line);
		std::vector<float> val;
		math::ParserStringToVectorFloat(_line, val, " ");
		critical_residual = val[0];
	}

	//Краевые НЕ ЧИТАЕМ из файла параметров
	//read Boundary values and vertexes
	struct _Dirichlet {
		/*Point<double> value_const;
		Point<bool> is_condition;*/
		std::vector<int> id_vertexes;
		std::function<Point<double>(Point<bool>&, int)> value;
	};
	std::vector<_Dirichlet> first_boundaries;
	{
		//все грани получают одну функцию и мы ничего из файла не читаем
		first_boundaries.resize(1);
		first_boundaries[0].value = [](Point<bool>& is_take, int id)->Point<double>
		{
			is_take = Point<bool>(true, true, true);
			return Point<double>(1, 1, 1);
		};

		std::vector<int> tmp_vert;
		for (int id_type = 0; id_type < boundary_faces.size(); id_type++)
		{
			
			for (int i = 0; i < boundary_faces[id_type].size(); i++)
			{
				//читаем не с 0, так как там номер базового тетраэдра
				for (int j = 1; j < boundary_faces[id_type][i].size(); j++)
				{
					tmp_vert.push_back(boundary_faces[id_type][i][j]);
				}
			}
		}
		math::MakeQuickSort(tmp_vert);
		math::MakeRemovalOfDuplication(tmp_vert, first_boundaries[0].id_vertexes);
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
	for (int i = 0; i < _materials.size(); i++)
	{
		solver_grid.AddDomain();
		auto domain = solver_grid.GetDomain(i);
		domain->forMech.SetE(_materials[i]._E);
		domain->forMech.SetV(_materials[i]._v);
	}

	math::ReadNonEmptyLine(f_properties, base_result_directory);
	fclose(f_properties);

	wchar_t _tmp_wc[1000];
	math::Char_To_Wchar_t(base_result_directory, _tmp_wc, 1000);
	CreateDirectory((LPCTSTR)_tmp_wc, NULL);

	//Формируем FEM-сетку и собираем 0-матрицу (без краевых) только один раз
	CSSD_Matrix<Tensor2Rank3D, Point<double>> global_base_SLAE;
	{
		printf("Initialization of grid...\n");
		solver_grid.Initialization(geo_grid, first_boundaries, second_boundary);
		printf_s("complite                   \n\n");

		printf("Creation the SLAE portrait...");
		solver_grid.CreationPortrait(global_base_SLAE);
		printf_s("complite\n\n");

		//SLAE assembling
		if(true) {
			printf("SLAE assembling...\n");
			std::vector<DenseMatrix<Tensor2Rank3D, Point<double>>> local_SLAE(solver_grid.GetElementsCount());
			omp_set_num_threads(math::NUM_THREADS);
#pragma omp parallel for schedule(dynamic)
			for (int id_elem = 0; id_elem < solver_grid.GetElementsCount(); id_elem++)
			{
				if (id_elem % 1000 == 0)
					printf("Solve element[%d]\r", id_elem);
				auto element = solver_grid.GetElement(id_elem);

				if (solver_grid.GetDomain(element->GetIdDomain())->forMech.GetE(0) <= 0)
				{
					//домен не должен идти в расчет
					local_SLAE[id_elem].SetSize(element->GetDOFsCount());
				}
				else {
					std::function< std::vector<std::vector<double>>(Point<double> X) > koefD = [&](Point<double> X)->std::vector<std::vector<double>>
					{
						return (solver_grid.GetDomain(element->GetIdDomain()))->forMech.GetD(3);
					};
					element->SolveLocalMatrix(local_SLAE[id_elem], koefD);
				}
			}
			for (int id_elem = 0; id_elem < solver_grid.GetElementsCount(); id_elem++)
			{
				if (id_elem % 1000 == 0)
					printf("Add the local matrix of element[%d]\r", id_elem);
				global_base_SLAE.SummPartOfMatrix(local_SLAE[id_elem], *solver_grid.GetElementDOFs(id_elem));
			}
			printf_s("                                                                                    \r");
			printf_s("complite\n\n");
		}
	}

	//Просчитываем все нужные варианты задач
	std::function<Point<double>(int, Point<double>, Point<double>, double)> func_boundary = [](int Type, Point<double> H, Point<double> X, double betta) -> Point<double>
	{
		Point<double> U;
		Point<double> M = H / 2.;
		switch (Type + 1)
		{
		case 1:
			U.x = betta * (X.x - M.x);
			U.y = 0;
			U.z = 0;
			break;
		case 6:
			U.x = betta * (X.y - M.y);
			U.y = 0;
			U.z = 0;
			break;
		case 5:
			U.x = betta * (X.z - M.z);
			U.y = 0;
			U.z = 0;
			break;
		case 2:
			U.x = 0;
			U.y = betta * (X.y - M.y);
			U.z = 0;
			break;
		case 4:
			U.x = 0;
			U.y = betta * (X.z - M.z);
			U.z = 0;
			break;
		case 3:
			U.x = 0;
			U.y = 0;
			U.z = betta * (X.z - M.z);
			break;
		default:
			break;
		}

		return U;
	};
	Point<double> sampleLenght, sampleCentr;
	sampleLenght = solver_grid.GetMaxCoordinate() - solver_grid.GetMinCoordinate();
	sampleCentr = (solver_grid.GetMaxCoordinate() + solver_grid.GetMinCoordinate()) / 2.0;
	double BETTA = 0.01 * sampleLenght.x;
	//Solution.resize(6);
	math::ResizeVector(Solution, 6, solver_grid.GetVertexCount());
	for (int STAGE = 0; STAGE < 6; STAGE++)
	{
		printf("========== STAGE = %d ==========\n", STAGE);

		//делаем действительную матрицу и перебрасываем туда 0-матрицу
		CSSD_Matrix<double, double> global_current_SLAE;
		math::MakeCopyMatrix_A_into_B(global_base_SLAE, global_current_SLAE);

		//добавляем первые краевые из func_boundary	
		{
			std::vector<int> use_id;
			printf("First boundary conditions...");
			//работа со строками
			for (int id_vertex = 0; id_vertex < solver_grid.boundary_vertexes.size(); id_vertex++)
			{
				int global_id = solver_grid.boundary_vertexes[id_vertex].GetDOFInLocalID(0);
				Point<double> vertex = solver_grid.GetCoordinateViaID(solver_grid.boundary_vertexes[id_vertex].GetIdNode(0));
				Point<double> boundary_value = func_boundary(STAGE, sampleLenght, vertex, BETTA);

				if (global_id < global_base_SLAE.GetMatrixSize())
				{
					auto enter_value = [&global_current_SLAE, &use_id](int value_id, double value) ->void
					{
						use_id.push_back(value_id);
						global_current_SLAE.X[value_id] = value;
						global_current_SLAE.F[value_id] = value;
						global_current_SLAE.Diag[value_id] = 1;
						//Обнуляем строку
						for (int i = 0; i < global_current_SLAE.A_down[value_id].size(); i++)
						{
							global_current_SLAE.A_down[value_id][i] = 0;
						}
						for (int i = 0; i < global_current_SLAE.A_up[value_id].size(); i++)
						{
							global_current_SLAE.A_up[value_id][i] = 0;
						}
						return;
					};

					enter_value(global_id * 3 + 0, boundary_value.x);
					enter_value(global_id * 3 + 1, boundary_value.y);
					enter_value(global_id * 3 + 2, boundary_value.z);
				}
			}

			//симметризация
			omp_set_num_threads(math::NUM_THREADS);
#pragma omp parallel for schedule(dynamic) 
			for (int id_row = 0; id_row < global_current_SLAE.GetMatrixSize(); id_row++)
			{
				if (id_row % 1000 == 0)
				{
					printf("\tcurrent %d/%d\r", id_row, global_current_SLAE.GetMatrixSize());
				}
				int iterator_in_boundary = 0;
				for (int jj = 0; jj < global_current_SLAE.id_column_for_A_up[id_row].size(); jj++)
				{
					int id_column = global_current_SLAE.id_column_for_A_up[id_row][jj];
					for (; iterator_in_boundary < use_id.size(); iterator_in_boundary++)
					{
						if (use_id[iterator_in_boundary] == id_column)
						{
							global_current_SLAE.F[id_row] -= global_current_SLAE.A_up[id_row][jj] * global_current_SLAE.F[id_column];
							global_current_SLAE.A_up[id_row][jj] = 0;

							break;
						}
						if (use_id[iterator_in_boundary] > id_column)
						{
							break;
						}
					}
				}

				iterator_in_boundary = 0;
				for (int jj = 0; jj < global_current_SLAE.id_column_for_A_down[id_row].size(); jj++)
				{
					int id_column = global_current_SLAE.id_column_for_A_down[id_row][jj];
					for (; iterator_in_boundary < use_id.size(); iterator_in_boundary++)
					{
						if (use_id[iterator_in_boundary] == id_column)
						{
							global_current_SLAE.F[id_row] -= global_current_SLAE.A_down[id_row][jj] * global_current_SLAE.F[id_column];
							global_current_SLAE.A_down[id_row][jj] = 0;

							break;
						}
						if (use_id[iterator_in_boundary] > id_column)
						{
							break;
						}
					}
				}
			}
		}

		//решаем СЛАУ
		if (true)
		{
			printf("Soluting SLAY... (%d)\n", global_current_SLAE.GetMatrixSize());
			int MaxSize = global_current_SLAE.GetMatrixSize();
			MaxSize = MaxSize / 10 < 10 ? 10 : MaxSize / 10;
			std::vector<double> best_solution;
			math::MakeCopyVector_A_into_B(global_current_SLAE.X, best_solution);
			double current_residual = 1, best_residual = 1e+25;
			int MAX_STEPS = 4;
			int ii = 1;
			//newSLAE.MakeDivideOnSqrtFormSumm();
			CSSD_Matrix<double, double> Precond;
			Precond.PrecondorSSOR(0.75, global_current_SLAE);
			for (int i = 0; i <= MAX_STEPS; i++)
			{
				double needed_residual = pow(10., -1 * (ii + 1));
				//current_residual /= 2.;
				printf_s("//---> I = %d/%d (full size %d) - needed residual %.2e\n", i, MAX_STEPS, global_current_SLAE.GetMatrixSize(), needed_residual);

				current_residual = abs(global_current_SLAE.MSG_PreconditioningSSOR(MaxSize, needed_residual, Precond));

				if (current_residual < needed_residual)
				{
					i = 0;
					ii++;
				}
				if (current_residual > best_residual * 100)
				{
					break;
				}
				if (current_residual <= critical_residual)
				{
					best_residual = current_residual;
					break;
				}
				if (current_residual < best_residual)
				{
					best_residual = current_residual;
					math::MakeCopyVector_A_into_B(global_current_SLAE.X, best_solution);
				}
			}
			math::MakeCopyVector_A_into_B(best_solution, global_current_SLAE.X);
			printf_s("//---> BEST residual %.2e\n", best_residual);
			//переписываем лучшее решение в векторный вид
			math::MakeCopyVector_A_into_B(global_current_SLAE.X, Solution[STAGE]);

			if (best_residual >= 1)
			{
				printf_s("There are problems in the solution\n");
				int a;
				scanf_s("%d", &a);
			}
		}

		//вывод результатов
		printf_s("\nPrint the mech result into .dat file...\n ");
		//without deformations
		if (true) {
			FILE* fout_tech;
			char name_u_tech[5000];
			sprintf_s(name_u_tech, "%s/U%d_non_deformation.dat", base_result_directory, STAGE);
			fopen_s(&fout_tech, name_u_tech, "w");
			char name_in_file[1000];
			sprintf_s(name_in_file, "Elastic_%d", STAGE);
			std::vector<std::vector<char>> name_value(5);
			char name_v_tmp[5][100];
			sprintf_s(name_v_tmp[0], "VonMises, MPa");
			sprintf_s(name_v_tmp[1], "Ux, mm");
			sprintf_s(name_v_tmp[2], "Uy, mm");
			sprintf_s(name_v_tmp[3], "Uz, mm");
			sprintf_s(name_v_tmp[4], "Material");
			switch (STAGE)
			{
			case 0: //x
				sprintf_s(name_v_tmp[1], "Ux");
				sprintf_s(name_v_tmp[2], "Stress_xx");
				sprintf_s(name_v_tmp[3], "Strain_xx");
				break;
			case 1: //y
				sprintf_s(name_v_tmp[1], "Uy");
				sprintf_s(name_v_tmp[2], "Stress_yy");
				sprintf_s(name_v_tmp[3], "Strain_yy");
				break;
			case 2: //z
				sprintf_s(name_v_tmp[1], "Uz");
				sprintf_s(name_v_tmp[2], "Stress_zz");
				sprintf_s(name_v_tmp[3], "Strain_zz");
				break;
			case 3: //yz
				sprintf_s(name_v_tmp[1], "Uy");
				sprintf_s(name_v_tmp[2], "Stress_yz");
				sprintf_s(name_v_tmp[3], "Strain_yz");
				break;
			case 4: //xz
				sprintf_s(name_v_tmp[1], "Ux");
				sprintf_s(name_v_tmp[2], "Stress_xz");
				sprintf_s(name_v_tmp[3], "Strain_xz");
				break;
			case 5: //xy
				sprintf_s(name_v_tmp[1], "Ux");
				sprintf_s(name_v_tmp[2], "Stress_xy");
				sprintf_s(name_v_tmp[3], "Strain_xy");
				break;
			}
			for (int i = 0; i < name_value.size(); i++)
			{
				name_value[i].resize(100);
				for (int j = 0; j < name_value[i].size(); j++)
				{
					name_value[i][j] = name_v_tmp[i][j];
				}
			}
			std::vector<std::vector<double>> value(name_value.size());
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

				Point<double> U = solver_grid.GetSolutionInPoint(i, Centr, Solution[STAGE]);
				Point<Point<double>> dU = solver_grid.GetDerevativeFromSolutionInPoint(i, Centr, Solution[STAGE]);
				auto EPS = solver_grid.GetStrainTensorFromSolutionInPoint(i, Centr, dU);
				auto SIG = solver_grid.GetStressTensorFromSolutionInPoint(i, Centr, EPS);
				double mises_stress = solver_grid.GetVonMisesStress(SIG) * 1e-6;
				U *= 1e+3;

				if (solver_grid.GetDomain(element->GetIdDomain())->forMech.GetE(0) > 0)
				{
					value[0][i] = mises_stress;
					value[1][i] = U.x;
					value[2][i] = U.y;
					value[3][i] = U.z;
					value[4][i] = element->GetIdDomain();

					switch (STAGE)
					{
					case 0: //x
						value[1][i] = U.x;
						value[2][i] = SIG.val[0][0];
						value[3][i] = EPS.val[0][0];
						break;
					case 1: //y
						value[1][i] = U.y;
						value[2][i] = SIG.val[1][1];
						value[3][i] = EPS.val[1][1];
						break;
					case 2: //z
						value[1][i] = U.z;
						value[2][i] = SIG.val[2][2];
						value[3][i] = EPS.val[2][2];
						break;
					case 3: //yz
						value[1][i] = U.y;
						value[2][i] = SIG.val[1][2];
						value[3][i] = EPS.val[1][2];
						break;
					case 4: //xz
						value[1][i] = U.x;
						value[2][i] = SIG.val[0][2];
						value[3][i] = EPS.val[0][2];
						break;
					case 5: //xy
						value[1][i] = U.x;
						value[2][i] = SIG.val[0][1];
						value[3][i] = EPS.val[0][1];
						break;
					}
				}
				else {
					value[0][i] = 0;
					value[1][i] = 0;
					value[2][i] = 0;
					value[3][i] = 0;
					value[4][i] = element->GetIdDomain();
				}
			}

			solver_grid.printTecPlot3D_DiffDomains(fout_tech, value, name_value, name_in_file);
			fclose(fout_tech);
		}


		//удаляем ненужную матрицу
		global_current_SLAE.~CSSD_Matrix();
	}

	//считаем эффективный тензор упругости
	math::ResizeVector(D_eff, 6, 6);
	{
		printf("Solve effective D... \n");
		double Volume = solver_grid.GetVolume();
		Volume = 0;
		for (int id_elem = 0; id_elem < solver_grid.GetElementsCount(); id_elem++)
		{
			if (solver_grid.GetElement(id_elem)->GetIdDomain() != 3)
			{
				Volume += solver_grid.GetElement(id_elem)->GetVolume();
			}
		}

		for (int I = 0; I < 6; I++)
		{
			for (int J = I; J < 6; J++)
			{
				printf_s("Solve D_eff[%d][%d]\r", I, J);
				double E_ij = 0;
				for (int id_elem = 0; id_elem < solver_grid.GetElementsCount(); id_elem++)
				{
					auto elem = solver_grid.GetElement(id_elem);

					std::function<double(Point<double>)> func = [I, J, &elem, &id_elem, &solver_grid, &Solution](Point<double> X) -> double
					{
						double res = 0;

						Point<Point<double>> dU_I = solver_grid.GetDerevativeFromSolutionInPoint(id_elem, X, Solution[I]);
						Point<Point<double>> dU_J = solver_grid.GetDerevativeFromSolutionInPoint(id_elem, X, Solution[J]);

						auto EPS_I = solver_grid.GetStrainTensorFromSolutionInPoint(id_elem, X, dU_I);
						auto SIG_J = solver_grid.GetStressTensorFromSolutionInPoint(id_elem, X, dU_J);

						res = solver_grid.GetDomain(elem->GetIdDomain())->forMech.mult_SIGMA_EPS(EPS_I, SIG_J);

						return res;
					};

					if (elem->GetIdDomain() != 3)
					{
						elem->SetIntegrationLaw(1);
						E_ij += elem->SolveIntegral(func);
					}
				}
				D_eff[I][J] = E_ij / (BETTA * BETTA * Volume);
				D_eff[J][I] = E_ij / (BETTA * BETTA * Volume);
				//if (i == j) D_eff[i][j] *= 2;
			}
		}
	}

	//вывод результатов
	printf_s("\nPrint the mech result into .dat file...\n ");
	//without deformations
	for (int STAGE = 0; false && STAGE < 6; STAGE++)
	{
		if (true) {
			FILE* fout_tech;
			char name_u_tech[5000];
			sprintf_s(name_u_tech, "%s/U%d_non_deformation.dat", base_result_directory, STAGE);
			fopen_s(&fout_tech, name_u_tech, "w");
			char name_in_file[1000];
			sprintf_s(name_in_file, "Elastic_%d", STAGE);
			std::vector<std::vector<char>> name_value(5);
			char name_v_tmp[5][100];
			sprintf_s(name_v_tmp[0], "VonMises");
			sprintf_s(name_v_tmp[1], "Ux");
			sprintf_s(name_v_tmp[2], "Uy");
			sprintf_s(name_v_tmp[3], "Uz");
			sprintf_s(name_v_tmp[4], "Material");
			switch (STAGE)
			{
			case 0: //x
				sprintf_s(name_v_tmp[1], "Ux");
				sprintf_s(name_v_tmp[2], "Stress_xx");
				sprintf_s(name_v_tmp[3], "Strain_xx");
				break;
			case 1: //y
				sprintf_s(name_v_tmp[1], "Uy");
				sprintf_s(name_v_tmp[2], "Stress_yy");
				sprintf_s(name_v_tmp[3], "Strain_yy");
				break;
			case 2: //z
				sprintf_s(name_v_tmp[1], "Uz");
				sprintf_s(name_v_tmp[2], "Stress_zz");
				sprintf_s(name_v_tmp[3], "Strain_zz");
				break;
			case 3: //yz
				sprintf_s(name_v_tmp[1], "Uy");
				sprintf_s(name_v_tmp[2], "Stress_yz");
				sprintf_s(name_v_tmp[3], "Strain_yz");
				break;
			case 4: //xz
				sprintf_s(name_v_tmp[1], "Ux");
				sprintf_s(name_v_tmp[2], "Stress_xz");
				sprintf_s(name_v_tmp[3], "Strain_xz");
				break;
			case 5: //xy
				sprintf_s(name_v_tmp[1], "Ux");
				sprintf_s(name_v_tmp[2], "Stress_xy");
				sprintf_s(name_v_tmp[3], "Strain_xy");
				break;
			}
			for (int i = 0; i < name_value.size(); i++)
			{
				name_value[i].resize(100);
				for (int j = 0; j < name_value[i].size(); j++)
				{
					name_value[i][j] = name_v_tmp[i][j];
				}
			}
			std::vector<std::vector<double>> value(name_value.size());
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

				Point<double> U = solver_grid.GetSolutionInPoint(i, Centr, Solution[STAGE]);
				Point<Point<double>> dU = solver_grid.GetDerevativeFromSolutionInPoint(i, Centr, Solution[STAGE]);
				auto EPS = solver_grid.GetStrainTensorFromSolutionInPoint(i, Centr, dU);
				auto SIG = solver_grid.GetStressTensorFromSolutionInPoint(i, Centr, EPS);
				double mises_stress = solver_grid.GetVonMisesStress(SIG);

				if (solver_grid.GetDomain(element->GetIdDomain())->forMech.GetE(0) > 0)
				{
					value[0][i] = mises_stress;
					value[1][i] = U.x;
					value[2][i] = U.y;
					value[3][i] = U.z;
					value[4][i] = element->GetIdDomain();

					switch (STAGE)
					{
					case 0: //x
						value[1][i] = U.x;
						value[2][i] = SIG.val[0][0];
						value[3][i] = EPS.val[0][0];
						break;
					case 1: //y
						value[1][i] = U.y;
						value[2][i] = SIG.val[1][1];
						value[3][i] = EPS.val[1][1];
						break;
					case 2: //z
						value[1][i] = U.z;
						value[2][i] = SIG.val[2][2];
						value[3][i] = EPS.val[2][2];
						break;
					case 3: //yz
						value[1][i] = U.y;
						value[2][i] = SIG.val[1][2];
						value[3][i] = EPS.val[1][2];
						break;
					case 4: //xz
						value[1][i] = U.x;
						value[2][i] = SIG.val[0][2];
						value[3][i] = EPS.val[0][2];
						break;
					case 5: //xy
						value[1][i] = U.x;
						value[2][i] = SIG.val[0][1];
						value[3][i] = EPS.val[0][1];
						break;
					}
				}
				else {
					value[0][i] = 0;
					value[1][i] = 0;
					value[2][i] = 0;
					value[3][i] = 0;
					value[4][i] = element->GetIdDomain();
				}
			}

			solver_grid.printTecPlot3D_DiffDomains(fout_tech, value, name_value, name_in_file);
			fclose(fout_tech);
		}
	}
	//вывод тензора упругости
	{
		FILE* fout_D;
		char name_D[5000];
		sprintf_s(name_D, "%s/D_eff.txt", base_result_directory);
		fopen_s(&fout_D, name_D, "a");
		std::vector<double> Volume_incl(solver_grid.GetDomainsCount());
		double fullVolume = solver_grid.GetVolume();
		for (int i = 0; i < solver_grid.GetElementsCount(); i++)
		{
			Volume_incl[solver_grid.GetElement(i)->GetIdDomain()] += solver_grid.GetElement(i)->GetVolume() / fullVolume * 100;
		}
		for (int d = 0; d < solver_grid.GetDomainsCount(); d++)
		{
			fprintf_s(fout_D, "------ D[%d] (E = %.2e;  v = %.2e) ------\n", d, solver_grid.GetDomain(d)->forMech.GetE(0), solver_grid.GetDomain(d)->forMech.v);
			fprintf_s(fout_D, "--conc = %.5lf--\n", Volume_incl[d]);
			solver_grid.GetDomain(d)->forMech.print_D(fout_D);
			fprintf_s(fout_D, "-------------------\n");
		}
		double v_eff_zz = (2 * D_eff[5][5] - D_eff[2][2]) / (2 * D_eff[5][5] - 2 * D_eff[2][2]);
		double E_eff_zz = 2 * D_eff[5][5] * (1 + v_eff_zz);
		fprintf_s(fout_D, "------ D_eff (E_eff_zz = %.2e;  v_eff_zz = %.2e) ------\n", E_eff_zz, v_eff_zz);
		for (int i = 0; i < D_eff.size(); i++)
		{
			for (int j = 0; j < D_eff[0].size(); j++)
			{
				if (D_eff[i][j] > 0)
				{
					fprintf_s(fout_D, " ");
				}
				fprintf_s(fout_D, "%.4e ", D_eff[i][j]);
			}
			fprintf_s(fout_D, "\n");
		}
		fprintf_s(fout_D, "-------------------\n");
		for (int i = 0; i < D_eff.size(); i++)
		{
			for (int j = 0; j < D_eff[0].size(); j++)
			{
				if (i < 3 && j < 3 || i == j)
				{
					if (D_eff[i][j] > 0)
					{
						fprintf_s(fout_D, " ");
					}
					fprintf_s(fout_D, "%.2e ", D_eff[i][j]);
				}
				else {
					fprintf_s(fout_D, " ");
					fprintf_s(fout_D, "%.2e ", 0.0);
				}
			}
			fprintf_s(fout_D, "\n");
		}
		fprintf_s(fout_D, "-------------------\n");

		fclose(fout_D);
	}

	solver_grid.~Grid_forMech();

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

	{
		char _line[1000];
		math::ReadNonEmptyLine_forNumbers(f_properties, _line);
		std::vector<int> val;
		math::ParserStringToVectorInt(_line, val, " ");
		math::NUM_THREADS = val[0];
	}

	math::SimpleGrid geo_grid; //input
	std::vector<std::vector<std::vector<int>>> boundary_faces;

	int num_in_FE;
	{
		char _line[1000];
		math::ReadNonEmptyLine_forNumbers(f_properties, _line);
		std::vector<int> val_i;
		math::ParserStringToVectorInt(_line, val_i, " ");

		num_in_FE = val_i[0];
	}
	math::ReadNonEmptyLine(f_properties, mesh_directory);
	switch (num_in_FE)
	{
	case 9:
		geo_grid.ReadFromMSH_v2208(mesh_directory, 1, 4, 4, boundary_faces);
		break;
	case 10:
		geo_grid.ReadFromMSH(mesh_directory, 1, 4, 4, boundary_faces);
		break;
	default:
		break;
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

	wchar_t _tmp_wc[1000];
	math::Char_To_Wchar_t(base_result_directory, _tmp_wc, 1000);
	CreateDirectory((LPCTSTR)_tmp_wc, NULL);

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

		//FEM::FEM_forElasticDeformation(
		//	is_print_logFile,
		//	critical_residual,
		//	geo_grid, //input
		//	first_boundaries, //input
		//	second_boundary, //input
		//	base_result_directory, //output
		//	solver_grid, //output
		//	Solution //output
		//);

		FEM::FEM_forElasticDeformation_renumeration(
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
			value[3].resize(solver_grid.GetVertexCount());
			value[4].resize(solver_grid.GetVertexCount());
			value[5].resize(solver_grid.GetVertexCount());
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

					/*value[3][i] = eps[0];
					value[4][i] = eps[1];
					value[5][i] = eps[2];*/
				}
				else {
					value[0][i] = 0;
					value[1][i] = 0;
					value[2][i] = 0;
				}

				if (value[0][i] != value[0][i]) value[0][i] = 0;
				if (value[1][i] != value[1][i]) value[1][i] = 0;
				if (value[2][i] != value[2][i]) value[2][i] = 0;

				/*value[3][i] = U.x;
				value[4][i] = U.y;
				value[5][i] = U.z;*/

				if (sigma_inv > sigma_inv_max)
				{
					sigma_inv_max = sigma_inv;
					elem_sigma_max = i;
				}
			}
			for (int id_DOF = 0; id_DOF < Solution.size(); id_DOF++)
			{
				solver_grid.MoveTheVertex(solver_grid.accordance_DOF_and_vertex[id_DOF], Solution[id_DOF]);

				value[3][id_DOF] = Solution[id_DOF].x;
				value[4][id_DOF] = Solution[id_DOF].y;
				value[5][id_DOF] = Solution[id_DOF].z;
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
			value[3].resize(solver_grid.GetVertexCount());
			value[4].resize(solver_grid.GetVertexCount());
			value[5].resize(solver_grid.GetVertexCount());
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

					/*value[3][i] = eps[0];
					value[4][i] = eps[1];
					value[5][i] = eps[2];*/
				}
				else {
					value[0][i] = 0;
					value[1][i] = 0;
					value[2][i] = 0;
				}

				if (value[0][i] != value[0][i]) value[0][i] = 0;
				if (value[1][i] != value[1][i]) value[1][i] = 0;
				if (value[2][i] != value[2][i]) value[2][i] = 0;


				/*value[3][i] = U.x;
				value[4][i] = U.y;
				value[5][i] = U.z;*/

				if (sigma_inv > sigma_inv_max)
				{
					sigma_inv_max = sigma_inv;
					elem_sigma_max = i;
				}
			}
			for (int id_DOF = 0; id_DOF < Solution.size(); id_DOF++)
			{
				value[3][id_DOF] = Solution[id_DOF].x;
				value[4][id_DOF] = Solution[id_DOF].y;
				value[5][id_DOF] = Solution[id_DOF].z;
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


void ElasticDeformation_MsFEM_Poly(char properties_file[1000])
{
	try {
		MsFEM::Grid_forMech_Poly_Order1 solver_grid; //output
		std::vector<Point<double>> Solution; //output
		CSSD_Matrix<Tensor2Rank3D, Point<double>> global_SLAE;

		clock_t t_after = clock();
		double start = omp_get_wtime();

		//char properties_file[1000] = { "E:/+cyl/800el/param.txt" };
		//char properties_file[1000] = { "E:/Box/100x10x200/BoxWithCracks/67k/param.txt" };
		//char properties_file[1000] = { "E:/Box/100x10x200/BoxWithCracks/638k/param.txt" };
		//char properties_file[1000] = { "E:/Box/100x50x200/param.txt" };
		//char properties_file[1000] = { "D:/P-ref_MsFEM/2_el/param_for_solver.txt" };
		//char properties_file[1000] = { "D:/P-ref_MsFEM/Sample_for_convergense/PolyMesh_1Order/H2_h1/param_for_solver.txt" };
		//char properties_file[1000] = { "D:/P-ref_MsFEM/CubMesh/5x5x5_Order1/param_for_solver.txt" };
		//char properties_file[1000] = { "./param_for_solver.txt" };
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
		char fine_mesh_dir[1000];
		char TEST_mesh_dir[1000];


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
		math::ReadNonEmptyLine(f_properties, fine_mesh_dir);
		math::SimpleGrid geo_grid; //input
		geo_grid.ReadFromNVTR(mesh_directory, -1); //-1 is special flag for polyhedral grid 

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
					case '\t':
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
					default:
					{
						_tmp_line[curr_i] = line[j];
						curr_i++;
					}
					break;
					}
				}

				first_boundaries[i].value = [value, is_condition](Point<bool>& is_take) /*-> Point<double>*/
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
		math::ReadNonEmptyLine(f_properties, TEST_mesh_dir);
		fclose(f_properties);

		wchar_t _tmp_wc[1000];
		math::Char_To_Wchar_t(base_result_directory, _tmp_wc, 1000);
		CreateDirectory((LPCTSTR)_tmp_wc, NULL);

		bool res = false;
		{
			for (int i = 0; i < _materials.size(); i++)
			{
				solver_grid.AddDomain();
				auto domain = solver_grid.GetDomain(i);
				domain->forMech.SetE(_materials[i]._E);
				domain->forMech.SetV(_materials[i]._v);
			}

			MsFEM::MsFEM_forElasticDeformation_OrderBF1_forPoly(
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

			clock_t t_before = clock();
			double end = omp_get_wtime();
			
			int num_micro_nodes = 0, micro_nodes_min = solver_grid.GetElement(0)->self_grid.GetVertexCount(), micro_nodes_max = 0;
			int num_micro_tetr = 0, micro_tetr_min = solver_grid.GetElement(0)->self_grid.GetElementsCount(), micro_tetr_max = 0;
			for (int i = 0; i < solver_grid.GetElementsCount(); i++)
			{
				int nodes = solver_grid.GetElement(i)->self_grid.GetVertexCount();
				num_micro_nodes += nodes;
				if (micro_nodes_min > nodes) micro_nodes_min = nodes;
				if (micro_nodes_max < nodes) micro_nodes_max = nodes;

				int tetr = solver_grid.GetElement(i)->self_grid.GetElementsCount();
				num_micro_tetr += tetr;
				if (micro_tetr_min > tetr) micro_tetr_min = tetr;
				if (micro_tetr_max < tetr) micro_tetr_max = tetr;
			}
			FILE* fout_U;
			char name_out_time[1000];
			sprintf_s(name_out_time, "%s/Time_result.txt", base_result_directory);
			fopen_s(&fout_U, name_out_time, "a");
			fprintf_s(fout_U, "---------------------------------------------\n");
			fprintf_s(fout_U, "TIME = %lf sec\n", (t_before - t_after) * 1.0 / CLK_TCK);
			fprintf_s(fout_U, "start = %.16g\nend = %.16g\ndiff = %.16g\n", start, end, end - start);
			fprintf_s(fout_U, "Micro nodes (min/max/summ): %d\t%d\t%d\n", micro_nodes_min, micro_nodes_max, num_micro_nodes);
			fprintf_s(fout_U, "Micro tetr (min/max/summ): %d\t%d\t%d\n", micro_tetr_min, micro_tetr_max, num_micro_tetr);

			fclose(fout_U);


			printf_s("Compare solution with FEM...\n\t%s\n", TEST_mesh_dir);
			Point<double> result;
			if (true)
			{
				math::SimpleGrid FEM_mesh_simple; //input
				FEM_mesh_simple.ReadFromNVTR(TEST_mesh_dir, 4);

				FEM::Grid_forMech FEM_mesh;
				struct _Dirichlet_fem {
					//Point<double> value;
					//Point<bool> is_condition;
					std::vector<int> id_vertexes;
					std::function<Point<double>(Point<bool>&, int)> value;
				};
				std::vector<_Dirichlet_fem> first_boundaries_fem;
				std::vector<_Neumann> second_boundaries_fem;
				FEM_mesh.Initialization(FEM_mesh_simple, first_boundaries_fem, second_boundaries_fem);
				std::vector<Point<double>> FEM_solution(FEM_mesh.GetDOFsCount());

				FILE* fin_fem;
				char name_in_fem[1000];
				sprintf_s(name_in_fem, "%s/Solution.txt", TEST_mesh_dir);
				fopen_s(&fin_fem, name_in_fem, "r");
				for (int i = 0; i < FEM_mesh.GetDOFsCount(); i++)
				{
					fscanf_s(fin_fem, "%lf %lf %lf", &FEM_solution[i].x, &FEM_solution[i].y, &FEM_solution[i].z);
				}
				fclose(fin_fem);

				std::vector<Point<double>> result_parallel(solver_grid.GetElementsCount());
				std::vector<Point<double>> result_parallel_FEM(solver_grid.GetElementsCount());
				omp_set_num_threads(math::NUM_THREADS);
#pragma omp parallel for schedule(dynamic)
				for (int id_macro = 0; id_macro < solver_grid.GetElementsCount(); id_macro++)
				{
					printf_s("Compare id_macro[%d/%d]                                                       \r", id_macro, solver_grid.GetElementsCount());
					std::function<Point<double>(Point<double>)> integral_value = [&id_macro, &Solution, &solver_grid, &FEM_mesh, &FEM_solution](Point<double> X) -> Point<double>
					{
						Point<double> U_FEHMM = solver_grid.GetSolutionInPoint(id_macro, X, Solution);

						double len;
						int id_fem = FEM_mesh.GetNearestElementID(X, len);
						Point<double> U_FEM = FEM_mesh.GetSolutionInPoint(id_fem, X, FEM_solution);

						Point<double> result;
						result.x = (U_FEHMM.x - U_FEM.x) / U_FEM.x;
						result.y = (U_FEHMM.y - U_FEM.y) / U_FEM.y;
						result.z = (U_FEHMM.z - U_FEM.z) / U_FEM.z;

						return Point<double>(result.x * result.x, result.y * result.y, result.z * result.z);
					};
					//solver_grid.GetElement(id_macro)->SetIntegrationLaw(4);
					//result += solver_grid.GetElement(id_macro)->SolveIntegral(integral_value);

					auto macro_el = solver_grid.GetElement(id_macro);
					for (int id_micro = 0; id_micro < macro_el->self_grid.GetElementsCount(); id_micro++)
					{
						auto micro_el = macro_el->self_grid.GetElement(id_micro);
						micro_el->SetIntegrationLaw(1);

						std::function<Point<double>(Point<double>)> integral_value_viaMicro = [&macro_el, &id_micro, &Solution, &solver_grid, &FEM_mesh, &FEM_solution](Point<double> X) -> Point<double>
						{
							Point<double> U_FEHMM;
							for (int id_bf = 0; id_bf < macro_el->GetDOFsCount(); id_bf++)
							{
								Point<double> bf_X = (*macro_el->GetBasisFunctionInLocalID_viaMicro(id_bf))(X, id_micro);
								U_FEHMM.x += bf_X.x * Solution[macro_el->GetIdNode(id_bf)].x;
								U_FEHMM.y += bf_X.y * Solution[macro_el->GetIdNode(id_bf)].y;
								U_FEHMM.z += bf_X.z * Solution[macro_el->GetIdNode(id_bf)].z;
							}

							double len;
							int id_fem = FEM_mesh.GetNearestElementID(X, len);
							Point<double> U_FEM = FEM_mesh.GetSolutionInPoint(id_fem, X, FEM_solution);

							Point<double> result;
							result.x = (U_FEHMM.x - U_FEM.x);
							result.y = (U_FEHMM.y - U_FEM.y);
							result.z = (U_FEHMM.z - U_FEM.z);

							return Point<double>(result.x * result.x, result.y * result.y, result.z * result.z);
						};
						std::function<Point<double>(Point<double>)> integral_value_viaMicro_FEM = [&macro_el, &id_micro, &Solution, &solver_grid, &FEM_mesh, &FEM_solution](Point<double> X) -> Point<double>
						{
							double len;
							int id_fem = FEM_mesh.GetNearestElementID(X, len);
							Point<double> U_FEM = FEM_mesh.GetSolutionInPoint(id_fem, X, FEM_solution);

							Point<double> result;
							result.x = U_FEM.x;
							result.y = U_FEM.y;
							result.z = U_FEM.z;

							return Point<double>(result.x * result.x, result.y * result.y, result.z * result.z);
						};
						
						result_parallel[id_macro] += micro_el->SolveIntegral(integral_value_viaMicro);
						result_parallel_FEM[id_macro] += micro_el->SolveIntegral(integral_value_viaMicro_FEM);
					}
				}

				Point<double> result_fem;
				for (int i = 0; i < result_parallel.size(); i++)
				{
					result += result_parallel[i];
					result_fem += result_parallel_FEM[i];
				}
				result.x = sqrt(result.x) / (sqrt(result_fem.x) < 1e-10 ? 1 : sqrt(result_fem.x));
				result.y = sqrt(result.y) / (sqrt(result_fem.y) < 1e-10 ? 1 : sqrt(result_fem.y));
				result.z = sqrt(result.z) / (sqrt(result_fem.z) < 1e-10 ? 1 : sqrt(result_fem.z));
			}
			if (false) //analitical solution
			{
				std::vector<Point<double>> result_parallel(solver_grid.GetElementsCount());
				std::vector<Point<double>> true_result_parallel(solver_grid.GetElementsCount());
				omp_set_num_threads(math::NUM_THREADS);
#pragma omp parallel for schedule(dynamic)
				for (int id_macro = 0; id_macro < solver_grid.GetElementsCount(); id_macro++)
				{
					printf_s("Compare id_macro[%d/%d]                                                       \r", id_macro, solver_grid.GetElementsCount());
					std::function<Point<double>(Point<double>)> integral_value = [&id_macro, &Solution, &solver_grid](Point<double> X) -> Point<double>
					{
						Point<double> U_FEHMM = solver_grid.GetSolutionInPoint(id_macro, X, Solution);

						double len;
						Point<double> U_FEM;
						U_FEM.x = 0;
						U_FEM.y = 0;
						Point<double> H = solver_grid.GetMaxCoordinate() - solver_grid.GetMinCoordinate();
						U_FEM.z = 0.01 * (X.z - H.z / 2.0);

						Point<double> result;
						result.x = (U_FEHMM.x - U_FEM.x);
						result.y = (U_FEHMM.y - U_FEM.y);
						result.z = (U_FEHMM.z - U_FEM.z) / U_FEM.z;

						return Point<double>(result.x * result.x, result.y * result.y, result.z * result.z);
					};
					//solver_grid.GetElement(id_macro)->SetIntegrationLaw(4);
					//result += solver_grid.GetElement(id_macro)->SolveIntegral(integral_value);

					auto macro_el = solver_grid.GetElement(id_macro);
					for (int id_micro = 0; id_micro < macro_el->self_grid.GetElementsCount(); id_micro++)
					{
						auto micro_el = macro_el->self_grid.GetElement(id_micro);
						micro_el->SetIntegrationLaw(1);

						std::function<Point<double>(Point<double>)> integral_value_viaMicro = [&macro_el, &id_micro, &Solution, &solver_grid](Point<double> X) -> Point<double>
						{
							Point<double> U_FEHMM;
							for (int id_bf = 0; id_bf < macro_el->GetDOFsCount(); id_bf++)
							{
								Point<double> bf_X = (*macro_el->GetBasisFunctionInLocalID_viaMicro(id_bf))(X, id_micro);
								U_FEHMM.x += bf_X.x * Solution[macro_el->GetIdNode(id_bf)].x;
								U_FEHMM.y += bf_X.y * Solution[macro_el->GetIdNode(id_bf)].y;
								U_FEHMM.z += bf_X.z * Solution[macro_el->GetIdNode(id_bf)].z;
							}

							Point<double> U_FEM;
							U_FEM.x = 0;
							U_FEM.y = 0;
							Point<double> H = solver_grid.GetMaxCoordinate() - solver_grid.GetMinCoordinate();
							U_FEM.z = 0.01 * (X.z - H.z / 2.0);

							Point<double> result;
							result.x = (U_FEHMM.x - U_FEM.x);
							result.y = (U_FEHMM.y - U_FEM.y);
							//result.z = (U_FEHMM.z - U_FEM.z) / U_FEM.z;
							result.z = (U_FEHMM.z - U_FEM.z);

							return Point<double>(result.x * result.x, result.y * result.y, result.z * result.z);
						};
						std::function<Point<double>(Point<double>)> integral_value_viaMicro_TRUE = [&macro_el, &id_micro, &Solution, &solver_grid](Point<double> X) -> Point<double>
						{
							Point<double> U_FEM;
							U_FEM.x = 0;
							U_FEM.y = 0;
							Point<double> H = solver_grid.GetMaxCoordinate() - solver_grid.GetMinCoordinate();
							U_FEM.z = 0.01 * (X.z - H.z / 2.0);

							Point<double> result;
							result.x = U_FEM.x;
							result.y = U_FEM.y;
							//result.z = (U_FEHMM.z - U_FEM.z) / U_FEM.z;
							result.z = U_FEM.z;

							return Point<double>(result.x * result.x, result.y * result.y, result.z * result.z);
						};

						result_parallel[id_macro] += micro_el->SolveIntegral(integral_value_viaMicro);
						true_result_parallel[id_macro] += micro_el->SolveIntegral(integral_value_viaMicro_TRUE);
					}
				}

				Point<double> true_result;
				for (int i = 0; i < result_parallel.size(); i++)
				{
					result += result_parallel[i];
					true_result += true_result_parallel[i];
				}
				result.x = sqrt(result.x) / (sqrt(true_result.x) < 1e-10 ? 1 : sqrt(true_result.x));
				result.y = sqrt(result.y) / (sqrt(true_result.y) < 1e-10 ? 1 : sqrt(true_result.y));
				result.z = sqrt(result.z) / (sqrt(true_result.z) < 1e-10 ? 1 : sqrt(true_result.z));
			}

			

			//FILE* fout_U;
			//char name_out_time[1000];
			//sprintf_s(name_out_time, "%s/Time_result.txt", base_result_directory);
			fopen_s(&fout_U, name_out_time, "a");
			fprintf_s(fout_U, "Error(X) %.5e\n", result.x);
			fprintf_s(fout_U, "Error(Y) %.5e\n", result.y);
			fprintf_s(fout_U, "Error(Z) %.5e\n", result.z);
			fprintf_s(fout_U, "Error(middle) %.5e\n", (result.x + result.y + result.z) / 3.0);
			fprintf_s(fout_U, "%.5e\t%.5e\t%.5e\n", result.x, result.y, result.z);
			fprintf_s(fout_U, "---------------------------------------------\n\n");
			fclose(fout_U);

			//output solution
			printf_s("Print the mech result into .dat file... ");
			omp_set_num_threads(math::NUM_THREADS);
#pragma omp parallel for schedule(dynamic)
			for (int id_elem = 0; id_elem < solver_grid.GetElementsCount(); id_elem++)
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
				value[0].resize(macro_element->self_grid.GetElementsCount());
				value[1].resize(macro_element->self_grid.GetElementsCount());
				value[2].resize(macro_element->self_grid.GetElementsCount());
				value[3].resize(macro_element->self_grid.GetElementsCount());
				value[4].resize(macro_element->self_grid.GetElementsCount());
				value[5].resize(macro_element->self_grid.GetElementsCount());
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

					Point<double> H = solver_grid.GetMaxCoordinate() - solver_grid.GetMinCoordinate();

					value[3][i] = U_macro.x;
					//value[4][i] = 0.01 * (Centr_micro.z - H.z / 2.0);
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
				math::SimpleGrid geo_grid_plane;
				char mesh_plane_directory[1000];
				sprintf_s(mesh_plane_directory, sizeof(mesh_plane_directory), "%s/Mesh_for_print_XOZ.dat", mesh_directory);
				geo_grid_plane.ReadFromSalomeDat(mesh_plane_directory, 2);

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
				value[0].resize(geo_grid_plane.nvtr.size());
				value[1].resize(geo_grid_plane.nvtr.size());
				value[2].resize(geo_grid_plane.nvtr.size());
				value[3].resize(geo_grid_plane.nvtr.size());
				value[4].resize(geo_grid_plane.nvtr.size());
				value[5].resize(geo_grid_plane.nvtr.size());
				double sigma_inv_max = 0;
				int elem_sigma_max = 0;
				for (int id_elem = 0; id_elem < geo_grid_plane.nvtr.size(); id_elem++)
				{
					Point<double> Target_point;
					for (int i = 0; i < geo_grid_plane.nvtr[id_elem].size(); i++)
					{
						Target_point += geo_grid_plane.xyz[geo_grid_plane.nvtr[id_elem][i]];
					}
					Target_point /= geo_grid_plane.nvtr[id_elem].size();

					Point<double> U_macro;
					Point<Point<double>> dU_macro;

					double len;
					int id_macro = solver_grid.GetNearestElementID(Target_point, len);
					if (id_macro >= 0)
					{
						auto macro_element = solver_grid.GetElement(id_macro);
						int id_micro = macro_element->self_grid.GetNearestElementID(Target_point, len);

						printf_s("id_elem[%d] - Macro[%d] - Micro[%d]                                        \r", id_elem, id_macro, id_micro);

						if (id_micro >= 0)
						{
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

							value[3][id_elem] = U_macro.x;
							value[4][id_elem] = U_macro.y;
							value[5][id_elem] = U_macro.z;
						}
					}
				}

				geo_grid_plane.printTecPlot3D(fout_tech, value, name_value, name_in_file);

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

			//FILE* fout_U;
			char name_out_U[1000];
			sprintf_s(name_out_U, "%s/Solution.txt", base_result_directory);
			fopen_s(&fout_U, name_out_U, "w");
			for (int i = 0; i < Solution.size(); i++)
			{
				fprintf_s(fout_U, "%.8e %.8e %.8e\n", Solution[i].x, Solution[i].y, Solution[i].z);
			}
			fclose(fout_U);


			solver_grid.~Grid_forMech_Poly_Order1();
		}

	}
	catch (const std::exception&)
	{
		printf_s("Error: ElasticDeformation_MsFEM_Poly\n");
	}
}

void ElastodynamicsProblem_Explicit(char properties_file[1000])
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
	//char properties_file[1000] = { "D:/Elastodynamic/homocyl/67k/param_for_solver.txt" };
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
	bool res = false;
	double TIME_h = step_size;
	if (is_STATIONARY) {
		end_iteration = 1;
		start_iteration = 0;
	}
	std::vector<Point<double>> U_prev(geo_grid.xyz.size()), U_prevprev(geo_grid.xyz.size()), U_curr(geo_grid.xyz.size());
	math::InitializationVector(U_curr, 0);
	math::InitializationVector(U_prev, 0);
	math::InitializationVector(U_prevprev, 0);

	CSSD_Matrix<Tensor2Rank3D, Point<double>> StiffnessMassMatrix;
	CSSD_Matrix<Tensor2Rank3D, Point<double>> MassMatrix;
	CSSD_Matrix<Tensor2Rank3D, Point<double>> VolumeForceVector;
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
	solver_grid.CreationPortrait(StiffnessMassMatrix);
	MassMatrix.SetMatrix(StiffnessMassMatrix);
	VolumeForceVector.F.resize(StiffnessMassMatrix.F.size());
	printf_s("\t\tcomplite\n");

	//строим базовые матрицы жесткости/массы/вектор силы тяжести
	//матрица массы идет с коэффициентом rpho
	//первые краевые без симметризации
	printf_s("================= Create matrix ================\n");
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
			return solver_grid.GetDomain(solver_grid.GetElement(elem)->GetIdDomain())->forThermal.rpho;
		};
		std::function<Point<double>(int, Point<double>)> VolumeForсe = [&](int elem, Point<double> X)->Point<double>
		{
			return Point<double>(0, 0, 0);
		};

		printf("Matrix assembling...\n");
		std::vector<DenseMatrix<Tensor2Rank3D, Point<double>>> local_SLAE_stiffness(solver_grid.GetElementsCount());
		std::vector<DenseMatrix<Tensor2Rank3D, Point<double>>> local_SLAE_mass(solver_grid.GetElementsCount());
		std::vector< std::vector<Point<double>>> local_force(solver_grid.GetElementsCount());
		for (int id_elem = 0; id_elem < solver_grid.GetElementsCount(); id_elem++)
		{
			if (id_elem % 1000 == 0)
				printf("Solve element[%d]\r", id_elem);
			auto element = solver_grid.GetElement(id_elem);

			std::function<std::vector<std::vector<double>>(Point<double>)> D = [&](Point<double> X) {return StiffnessCoef(id_elem, X); };
			std::function<double(Point<double>)> M = [&](Point<double> X) {return MassCoef(id_elem, X); };
			std::function<Point<double>(Point<double>)> F = [&](Point<double> X) {return VolumeForсe(id_elem, X); };

			element->SolveLocalMatrix(local_SLAE_stiffness[id_elem], D);
			element->SolveMassMatrix(local_SLAE_mass[id_elem], M);
			element->SolveRightSide(local_force[id_elem], F);
		}
		for (int id_elem = 0; id_elem < solver_grid.GetElementsCount(); id_elem++)
		{
			if (id_elem % 1000 == 0)
				printf("Add the local matrix of element[%d]\r", id_elem);
			StiffnessMassMatrix.SummPartOfMatrix(local_SLAE_stiffness[id_elem], *solver_grid.GetElementDOFs(id_elem));
			MassMatrix.SummPartOfMatrix(local_SLAE_mass[id_elem], *solver_grid.GetElementDOFs(id_elem));
			VolumeForceVector.SummPartOfVector(local_force[id_elem], *solver_grid.GetElementDOFs(id_elem));
		}
		printf_s("                                                                                    \r");
		printf_s("\t\tcomplite\n");

		//суммируем массу и жесткость
		{
			for (int i = 0; i < StiffnessMassMatrix.Diag.size(); i++)
			{
				StiffnessMassMatrix.Diag[i] += MassMatrix.Diag[i] * 1.0 / (TIME_h * TIME_h);

				for (int j = 0; j < StiffnessMassMatrix.A_down[i].size(); j++)
				{
					StiffnessMassMatrix.A_down[i][j] += MassMatrix.A_down[i][j] * 1.0 / (TIME_h * TIME_h);
				}

				for (int j = 0; j < StiffnessMassMatrix.A_up[i].size(); j++)
				{
					StiffnessMassMatrix.A_up[i][j] += MassMatrix.A_up[i][j] * 1.0 / (TIME_h * TIME_h);
				}
			}
		}
	}

	for (int id_STEP = start_iteration; id_STEP < end_iteration; id_STEP++)
	{
		printf_s("\n================= Start solution of %d STEP (time = %.2e) ================\n", id_STEP, TIME_h * id_STEP);
		double TIME_curr = TIME_h * id_STEP;
		bool is_print_result = false;
		char result_directory[1000];
		if (id_STEP % step_for_out == 0)
		{
			is_print_result = true;
			sprintf_s(result_directory, sizeof(result_directory), "%s", base_result_directory);
			//sprintf_s(result_directory, sizeof(result_directory), "%s/STEP_%d_t=%.2e", base_result_directory, id_STEP, TIME_curr);
			wchar_t _tmp_wc[1000];
			math::Char_To_Wchar_t(result_directory, _tmp_wc, 1000);
			CreateDirectory((LPCTSTR)_tmp_wc, NULL);
		}

		//переопределяем источник и строим вектор в правую часть (через вторые краевые)
		CSSD_Matrix<Tensor2Rank3D, Point<double>> SourseVector;
		SourseVector.F.resize(StiffnessMassMatrix.F.size());
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

			for (int id_triangle = 0; id_triangle < solver_grid.boundary_faces.size(); id_triangle++)
			{
				std::vector<Point<double>> local_vector_SLAE;
				solver_grid.boundary_faces[id_triangle].SolveLocalBoundaryVector(local_vector_SLAE, solver_grid.boundary_faces[id_triangle].boundary_value);
				SourseVector.SummPartOfVector(local_vector_SLAE, *solver_grid.boundary_faces[id_triangle].GetElementDOFs());
			};
		}

		//считаем правую часть основного уравнения
		//пишем в StiffnessMassMatrix.F
		{
			std::vector<Point<double>> Mu_prev;
			MassMatrix.MultiplicationMatrixVector(U_prev, Mu_prev);

			std::vector<Point<double>> Mu_prevprev;
			MassMatrix.MultiplicationMatrixVector(U_prevprev, Mu_prevprev);

			for (int i = 0; i < StiffnessMassMatrix.F.size(); i++)
			{
				StiffnessMassMatrix.F[i] = Point<double>(0,0,0)
					+ Mu_prev[i] * 2.0 / (TIME_h * TIME_h)
					- Mu_prevprev[i] * 1.0 / (TIME_h * TIME_h)
					//- Ku_prev[i]
					+ VolumeForceVector.F[i]
					+ SourseVector.F[i];
			}
		}

		//решаем СЛАУ
		{
			CSSD_Matrix<double, double> StiffnessMassMatrix_doubleSLAE;
			math::MakeCopyMatrix_A_into_B(StiffnessMassMatrix, StiffnessMassMatrix_doubleSLAE);
						
			printf("First boundary conditions...");
			if (true) { //new version
				//обнуляем строки
				for (int id_vertex = 0; id_vertex < solver_grid.boundary_vertexes.size(); id_vertex++)
				{
					auto boundary = &(solver_grid.boundary_vertexes[id_vertex]);

					int id_type = boundary->id_type;
					Point<bool> is_take;
					Point<double> boundary_value = boundary->boundary_value(is_take, solver_grid.boundary_vertexes[id_vertex].GetIdNode(0));
					int global_id = boundary->GetDOFInLocalID(0);

					if (global_id < MassMatrix.GetMatrixSize())
					{
						auto enter_value = [&StiffnessMassMatrix_doubleSLAE](int value_id, double value) ->void
						{
							StiffnessMassMatrix_doubleSLAE.X[value_id] = value;
							StiffnessMassMatrix_doubleSLAE.F[value_id] = value;
							StiffnessMassMatrix_doubleSLAE.Diag[value_id] = 1;
							//Обнуляем строку
							for (int i = 0; i < StiffnessMassMatrix_doubleSLAE.A_down[value_id].size(); i++)
							{
								StiffnessMassMatrix_doubleSLAE.A_down[value_id][i] = 0;
							}
							for (int i = 0; i < StiffnessMassMatrix_doubleSLAE.A_up[value_id].size(); i++)
							{
								StiffnessMassMatrix_doubleSLAE.A_up[value_id][i] = 0;
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
#pragma omp parallel for schedule(dynamic) 
				for (int id_row = 0; id_row < StiffnessMassMatrix_doubleSLAE.GetMatrixSize(); id_row++)
				{
					if (id_row % 100 == 0)
					{
						printf("\tcurrent %d/%d\r", id_row, StiffnessMassMatrix_doubleSLAE.GetMatrixSize());
					}
					int iterator_in_boundary = 0;
					for (int jj = 0; jj < StiffnessMassMatrix_doubleSLAE.id_column_for_A_up[id_row].size(); jj++)
					{
						int id_column = StiffnessMassMatrix_doubleSLAE.id_column_for_A_up[id_row][jj];
						for (/*iterator_in_boundary = 0*/; iterator_in_boundary < solver_grid.boundary_vertexes.size(); iterator_in_boundary++)
						{
							if (solver_grid.boundary_vertexes[iterator_in_boundary].GetDOFInLocalID(0) * 3 + 0 == id_column)
							{
								Point<bool> is_take;
								Point<double> boundary_value = solver_grid.boundary_vertexes[iterator_in_boundary].boundary_value(is_take, solver_grid.boundary_vertexes[iterator_in_boundary].GetIdNode(0));

								if (is_take.x)
								{
									StiffnessMassMatrix_doubleSLAE.F[id_row] -= StiffnessMassMatrix_doubleSLAE.A_up[id_row][jj] * StiffnessMassMatrix_doubleSLAE.F[id_column];
									StiffnessMassMatrix_doubleSLAE.A_up[id_row][jj] = 0;
								}

								break;
							}
							if (solver_grid.boundary_vertexes[iterator_in_boundary].GetDOFInLocalID(0) * 3 + 1 == id_column)
							{
								Point<bool> is_take;
								Point<double> boundary_value = solver_grid.boundary_vertexes[iterator_in_boundary].boundary_value(is_take, solver_grid.boundary_vertexes[iterator_in_boundary].GetIdNode(0));

								if (is_take.y)
								{
									StiffnessMassMatrix_doubleSLAE.F[id_row] -= StiffnessMassMatrix_doubleSLAE.A_up[id_row][jj] * StiffnessMassMatrix_doubleSLAE.F[id_column];
									StiffnessMassMatrix_doubleSLAE.A_up[id_row][jj] = 0;
								}

								break;
							}
							if (solver_grid.boundary_vertexes[iterator_in_boundary].GetDOFInLocalID(0) * 3 + 2 == id_column)
							{
								Point<bool> is_take;
								Point<double> boundary_value = solver_grid.boundary_vertexes[iterator_in_boundary].boundary_value(is_take, solver_grid.boundary_vertexes[iterator_in_boundary].GetIdNode(0));

								if (is_take.z)
								{
									StiffnessMassMatrix_doubleSLAE.F[id_row] -= StiffnessMassMatrix_doubleSLAE.A_up[id_row][jj] * StiffnessMassMatrix_doubleSLAE.F[id_column];
									StiffnessMassMatrix_doubleSLAE.A_up[id_row][jj] = 0;
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
					for (int jj = 0; jj < StiffnessMassMatrix_doubleSLAE.id_column_for_A_down[id_row].size(); jj++)
					{
						int id_column = StiffnessMassMatrix_doubleSLAE.id_column_for_A_down[id_row][jj];
						for (/*iterator_in_boundary = 0*/; iterator_in_boundary < solver_grid.boundary_vertexes.size(); iterator_in_boundary++)
						{
							if (solver_grid.boundary_vertexes[iterator_in_boundary].GetDOFInLocalID(0) * 3 + 0 == id_column)
							{
								Point<bool> is_take;
								Point<double> boundary_value = solver_grid.boundary_vertexes[iterator_in_boundary].boundary_value(is_take, solver_grid.boundary_vertexes[iterator_in_boundary].GetIdNode(0));

								if (is_take.x)
								{
									StiffnessMassMatrix_doubleSLAE.F[id_row] -= StiffnessMassMatrix_doubleSLAE.A_down[id_row][jj] * StiffnessMassMatrix_doubleSLAE.F[id_column];
									StiffnessMassMatrix_doubleSLAE.A_down[id_row][jj] = 0;
								}

								break;
							}
							if (solver_grid.boundary_vertexes[iterator_in_boundary].GetDOFInLocalID(0) * 3 + 1 == id_column)
							{
								Point<bool> is_take;
								Point<double> boundary_value = solver_grid.boundary_vertexes[iterator_in_boundary].boundary_value(is_take, solver_grid.boundary_vertexes[iterator_in_boundary].GetIdNode(0));

								if (is_take.y)
								{
									StiffnessMassMatrix_doubleSLAE.F[id_row] -= StiffnessMassMatrix_doubleSLAE.A_down[id_row][jj] * StiffnessMassMatrix_doubleSLAE.F[id_column];
									StiffnessMassMatrix_doubleSLAE.A_down[id_row][jj] = 0;
								}

								break;
							}
							if (solver_grid.boundary_vertexes[iterator_in_boundary].GetDOFInLocalID(0) * 3 + 2 == id_column)
							{
								Point<bool> is_take;
								Point<double> boundary_value = solver_grid.boundary_vertexes[iterator_in_boundary].boundary_value(is_take, solver_grid.boundary_vertexes[iterator_in_boundary].GetIdNode(0));

								if (is_take.z)
								{
									StiffnessMassMatrix_doubleSLAE.F[id_row] -= StiffnessMassMatrix_doubleSLAE.A_down[id_row][jj] * StiffnessMassMatrix_doubleSLAE.F[id_column];
									StiffnessMassMatrix_doubleSLAE.A_down[id_row][jj] = 0;
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

			//начальное приближение
			if (id_STEP == 0)
			{
				for (int i = 0; i < StiffnessMassMatrix_doubleSLAE.GetMatrixSize() / 3; i++)
				{
					StiffnessMassMatrix_doubleSLAE.X[i * 3 + 0] = 1e-10;
					StiffnessMassMatrix_doubleSLAE.X[i * 3 + 1] = 1e-10;
					StiffnessMassMatrix_doubleSLAE.X[i * 3 + 2] = 1e-10;
				}
			}
			else {
				for (int i = 0; i < StiffnessMassMatrix_doubleSLAE.GetMatrixSize() / 3; i++)
				{
					StiffnessMassMatrix_doubleSLAE.X[i * 3 + 0] = U_prev[i].x;
					StiffnessMassMatrix_doubleSLAE.X[i * 3 + 1] = U_prev[i].y;
					StiffnessMassMatrix_doubleSLAE.X[i * 3 + 2] = U_prev[i].z;
				}
			}

			printf("Soluting SLAY... (%d)\n", StiffnessMassMatrix_doubleSLAE.GetMatrixSize());
			int MaxSize = StiffnessMassMatrix_doubleSLAE.GetMatrixSize();
			MaxSize = MaxSize / 10 < 1000 ? 1000 : MaxSize / 10;
			std::vector<double> best_solution;
			math::MakeCopyVector_A_into_B(StiffnessMassMatrix_doubleSLAE.X, best_solution);
			double current_residual = 1, best_residual = 1e+25;
			int MAX_STEPS = 4;
			int ii = 1;
			CSSD_Matrix<double, double> Precond;
			Precond.PrecondorSSOR(0.75, StiffnessMassMatrix_doubleSLAE);
			for (int i = 0; i <= MAX_STEPS; i++)
			{
				double needed_residual = pow(10., -1 * (ii + 1));
				//current_residual /= 2.;
				printf_s("//---> I = %d/%d (full size %d) - needed residual %.2e\n", i, MAX_STEPS, StiffnessMassMatrix_doubleSLAE.GetMatrixSize(), needed_residual);

				StiffnessMassMatrix_doubleSLAE.print_logs = true;
				//current_residual = abs(MassMatrix_doubleSLAE.MSG(MaxSize, needed_residual));
				//current_residual = abs(MassMatrix_doubleSLAE.BCG_Stab2(MaxSize, needed_residual));
				current_residual = abs(StiffnessMassMatrix_doubleSLAE.MSG_PreconditioningSSOR(MaxSize, needed_residual, Precond));

				if (current_residual < needed_residual)
				{
					i = 0;
					ii++;
				}
				if (current_residual <= critical_residual)
				{
					best_residual = current_residual;
					break;
				}
				if (current_residual < best_residual)
				{
					best_residual = current_residual;
					math::MakeCopyVector_A_into_B(StiffnessMassMatrix_doubleSLAE.X, best_solution);
				}
			}
			math::MakeCopyVector_A_into_B(best_solution, StiffnessMassMatrix_doubleSLAE.X);
			printf_s("//---> BEST residual %.2e\n", best_residual);
			math::MakeCopyVector_A_into_B(StiffnessMassMatrix_doubleSLAE.X, U_curr);

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
			sprintf_s(name_u_tech, "%s/U_deformation_%d_t%.2e.dat", result_directory, id_STEP, TIME_curr);
			fopen_s(&fout_tech, name_u_tech, "w");
			char name_in_file[1000];
			sprintf_s(name_in_file, "step%d_time%.2e", id_STEP, TIME_curr);
			std::vector<std::vector<char>> name_value(6);
			char name_v_tmp[6][100];
			sprintf_s(name_v_tmp[0], "sigma_Mises");
			sprintf_s(name_v_tmp[1], "eps_inv");
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

				//Point<double> U = solver_grid.GetSolutionInPoint(i, Centr, U_curr);
				Point<Point<double>> dU = solver_grid.GetDerevativeFromSolutionInPoint(i, Centr, U_curr);
				auto Eps = solver_grid.GetStrainTensorFromSolutionInPoint(i, Centr, dU);
				auto Sigma = solver_grid.GetStressTensorFromSolutionInPoint(i, Centr, Eps);
				auto MisesSigma = solver_grid.GetVonMisesStress(Sigma);
				
				auto Eps_inv = sqrt((Eps.val[0][0] - Eps.val[1][1]) * (Eps.val[0][0] - Eps.val[1][1])
					+ (Eps.val[1][1] - Eps.val[2][2]) * (Eps.val[1][1] - Eps.val[2][2])
					+ (Eps.val[0][0] - Eps.val[2][2]) * (Eps.val[0][0] - Eps.val[2][2])
					+ 3 * (Eps.val[0][1] * Eps.val[1][0] + Eps.val[0][2] * Eps.val[2][0] + Eps.val[1][2] * Eps.val[2][1]) / 2.0) * sqrt(2.) / 3.;

				value[0][i] = MisesSigma;
				value[1][i] = Eps_inv;
				value[2][i] = Sigma.val[2][2];
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

			//solver_grid.MoveCoordinates(U_curr);
			solver_grid.printTecPlot3D_DiffDomains(fout_tech, value, name_value, name_in_file);
			//solver_grid.ReMoveCoordinates(U_curr);
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
			sprintf_s(name_u_tech, "%s/U_plane_YZ_s%d_t%.2e.dat", result_directory, id_STEP, TIME_curr);
			fopen_s(&fout_tech, name_u_tech, "w");
			char name_in_file[1000];
			sprintf_s(name_in_file, "Time_%.4e_YZ", TIME_curr);
			std::vector<std::vector<char>> name_value(6);
			char name_v_tmp[6][100];
			sprintf_s(name_v_tmp[0], "sigma_Mises");
			sprintf_s(name_v_tmp[1], "eps_inv");
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
					auto Eps = solver_grid.GetStrainTensorFromSolutionInPoint(id_elem, Centr, dU);
					auto Sigma = solver_grid.GetStressTensorFromSolutionInPoint(id_elem, Centr, Eps);
					auto MisesSigma = solver_grid.GetVonMisesStress(Sigma);

					auto Eps_inv = sqrt((Eps.val[0][0] - Eps.val[1][1]) * (Eps.val[0][0] - Eps.val[1][1])
						+ (Eps.val[1][1] - Eps.val[2][2]) * (Eps.val[1][1] - Eps.val[2][2])
						+ (Eps.val[0][0] - Eps.val[2][2]) * (Eps.val[0][0] - Eps.val[2][2])
						+ 3 * (Eps.val[0][1] * Eps.val[1][0] + Eps.val[0][2] * Eps.val[2][0] + Eps.val[1][2] * Eps.val[2][1]) / 2.0) * sqrt(2.) / 3.;

					value[0][i] = MisesSigma;
					value[1][i] = Eps_inv;
					value[2][i] = Sigma.val[2][2];

					value[3][i] = U.x;
					value[4][i] = U.y;
					value[5][i] = U.z;
				}
			}

			grid_YZ_plane.printTecPlot3D(fout_tech, value, name_value, name_in_file);
			fclose(fout_tech);
		}
		printf_s("\t complite\n");

		//solver_grid.~Grid_forMech();
	}
}
void ElastodynamicsProblem_ExplicitSimple(char properties_file[1000])
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
	//char properties_file[1000] = { "D:/Elastodynamic/homocyl/67k/param_for_solver.txt" };
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


	{
		char _line[1000];
		std::vector<int> val;

		math::ReadNonEmptyLine_forNumbers(f_properties, _line);
		math::ParserStringToVectorInt(_line, val, " "); 
		if (val[0] == 1)
			is_print_logFile = true;


		math::ReadNonEmptyLine_forNumbers(f_properties, _line);
		math::ParserStringToVectorInt(_line, val, " ");
		math::NUM_THREADS = val[0];
	}

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

	wchar_t _tmp_wc[1000];
	math::Char_To_Wchar_t(base_result_directory, _tmp_wc, 1000);
	CreateDirectory((LPCTSTR)_tmp_wc, NULL);

	bool res = false;
	double TIME_h = step_size;
	if (is_STATIONARY) {
		end_iteration = 1;
		start_iteration = 0;
	}
	std::vector<Point<double>> U_prev(geo_grid.xyz.size()), U_prevprev(geo_grid.xyz.size()), U_curr(geo_grid.xyz.size());
	math::InitializationVector(U_curr, 0);
	math::InitializationVector(U_prev, 0);
	math::InitializationVector(U_prevprev, 0);

	CSSD_Matrix<Tensor2Rank3D, Point<double>> StiffnessMatrix;
	CSSD_Matrix<Tensor2Rank3D, Point<double>> MassMatrix;
	CSSD_Matrix<Tensor2Rank3D, Point<double>> VolumeForceVector;
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
	MassMatrix.SetMatrix(StiffnessMatrix);
	VolumeForceVector.F.resize(StiffnessMatrix.F.size());
	printf_s("\t\tcomplite\n");

	//строим базовые матрицы жесткости/массы/вектор силы тяжести
	//матрица массы идет с коэффициентом rpho
	printf_s("================= Create matrix ================\n");
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
			return solver_grid.GetDomain(solver_grid.GetElement(elem)->GetIdDomain())->forThermal.rpho;
		};
		std::function<Point<double>(int, Point<double>)> VolumeForсe = [&](int elem, Point<double> X)->Point<double>
		{
			return Point<double>(0, 0, 0);
		};

		printf("Matrix assembling...\n");
		std::vector<DenseMatrix<Tensor2Rank3D, Point<double>>> local_SLAE_stiffness(solver_grid.GetElementsCount());
		std::vector<DenseMatrix<Tensor2Rank3D, Point<double>>> local_SLAE_mass(solver_grid.GetElementsCount());
		std::vector< std::vector<Point<double>>> local_force(solver_grid.GetElementsCount());
		omp_set_num_threads(math::NUM_THREADS);
#pragma omp parallel for schedule(dynamic)
		for (int id_elem = 0; id_elem < solver_grid.GetElementsCount(); id_elem++)
		{
			if (id_elem % 1000 == 0)
				printf("Solve element[%d]\r", id_elem);
			auto element = solver_grid.GetElement(id_elem);

			std::function<std::vector<std::vector<double>>(Point<double>)> D = [&](Point<double> X) {return StiffnessCoef(id_elem, X); };
			std::function<double(Point<double>)> M = [&](Point<double> X) {return MassCoef(id_elem, X); };
			std::function<Point<double>(Point<double>)> F = [&](Point<double> X) {return VolumeForсe(id_elem, X); };

			element->SolveLocalMatrix(local_SLAE_stiffness[id_elem], D);
			element->SolveMassMatrix(local_SLAE_mass[id_elem], M);
			element->SolveRightSide(local_force[id_elem], F);
		}
		for (int id_elem = 0; id_elem < solver_grid.GetElementsCount(); id_elem++)
		{
			if (id_elem % 1000 == 0)
				printf("Add the local matrix of element[%d]\r", id_elem);
			StiffnessMatrix.SummPartOfMatrix(local_SLAE_stiffness[id_elem], *solver_grid.GetElementDOFs(id_elem));
			MassMatrix.SummPartOfMatrix(local_SLAE_mass[id_elem], *solver_grid.GetElementDOFs(id_elem));
			VolumeForceVector.SummPartOfVector(local_force[id_elem], *solver_grid.GetElementDOFs(id_elem));
		}
		printf_s("                                                                                    \r");
		printf_s("\t\tcomplite\n");
	}

	///-------------------
	Point<double> point_for_out(0, 0, 0.5 - 0.01);
	std::vector<double> Uz_in_point;
	std::vector<double> Time_for_Uz;
	double TIME_L;
	int id_elem_for_out;
	char name_uz[5000];
	{
		double v = solver_grid.GetDomain(0)->forMech.GetLongitudinalWaveVelocity_Vp(solver_grid.GetDomain(0)->forThermal.rpho);
		double L = 0.5;
		TIME_L = L / v;
		double len;
		id_elem_for_out = solver_grid.GetNearestElementID(point_for_out, len);

		wchar_t _tmp_wc[1000];
		math::Char_To_Wchar_t(base_result_directory, _tmp_wc, 1000);
		CreateDirectory((LPCTSTR)_tmp_wc, NULL);

		FILE* fout_Uz;
		sprintf_s(name_uz, "%s/Uz_in_point.txt", base_result_directory);
		fopen_s(&fout_Uz, name_uz, "w");
		fprintf_s(fout_Uz, "Point: (%.2e, %.2e, %.2e)\n", point_for_out.x, point_for_out.y, point_for_out.z);
		fprintf_s(fout_Uz, "Id element: %d\n", id_elem_for_out);
		fprintf_s(fout_Uz, "T0: %.5e sec\n", TIME_L);
		fclose(fout_Uz);

		if (id_elem_for_out < 0) return;
	}
	///-------------------


	for (int id_STEP = start_iteration; id_STEP <= end_iteration; id_STEP++)
	{
		printf_s("\n================= Start solution of %d STEP (time = %.2e) ================\n", id_STEP, TIME_h * id_STEP);
		double TIME_curr = TIME_h * id_STEP+ TIME_h;
		bool is_print_result = false;
		char result_directory[1000];
		if (id_STEP % step_for_out == 0)
		{
			is_print_result = true;
			sprintf_s(result_directory, sizeof(result_directory), "%s", base_result_directory);
			//sprintf_s(result_directory, sizeof(result_directory), "%s/STEP_%d_t=%.2e", base_result_directory, id_STEP, TIME_curr);
			wchar_t _tmp_wc[1000];
			math::Char_To_Wchar_t(result_directory, _tmp_wc, 1000);
			CreateDirectory((LPCTSTR)_tmp_wc, NULL);
		}

		//переопределяем источник и строим вектор в правую часть (через вторые краевые)
		CSSD_Matrix<Tensor2Rank3D, Point<double>> SourseVector;
		SourseVector.F.resize(StiffnessMatrix.F.size());
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
			std::function<Point<double>(Point<double>)> new_sourse_value_sin = [TIME_curr, power, &solver_grid, TIME_h](Point<double> X) -> Point<double>
			{
				//double sourse = exp(-pow(TIME_curr - tay_0, degree) / pow(tay_plus, degree));
				//return power;
				//return Point<double>(power.x * sourse, power.y * sourse, power.z * sourse);
				double v = solver_grid.GetDomain(0)->forMech.GetLongitudinalWaveVelocity_Vp(solver_grid.GetDomain(0)->forThermal.rpho);
				double L = 0.5;
				double t0 = v / L;
				double w = 2 * M_PI / t0 * 10;
				double alpha = 5e+4;

				double sourse = 10e+6 * exp(-1 * alpha * TIME_curr) * sin(w * TIME_curr);
				if (TIME_curr > t0) sourse = 0;
				return Point<double>(0 * sourse, 0 * sourse, -1 * sourse /** TIME_h / 14.0e-7*/);
			};
			if (second_boundary.size() > 0)
			{
				second_boundary[0].value = new_sourse_value_sin;
				//second_boundary[1].value = new_sourse_value_sin;
				//second_boundary[2].value = new_sourse_value;

				for (int i = 0; i < solver_grid.boundary_faces.size(); i++)
				{
					solver_grid.boundary_faces[i].boundary_value = second_boundary[solver_grid.boundary_faces[i].id_type].value;
				}
			}

			for (int id_triangle = 0; id_triangle < solver_grid.boundary_faces.size(); id_triangle++)
			{
				std::vector<Point<double>> local_vector_SLAE;
				solver_grid.boundary_faces[id_triangle].SolveLocalBoundaryVector(local_vector_SLAE, solver_grid.boundary_faces[id_triangle].boundary_value);
				SourseVector.SummPartOfVector(local_vector_SLAE, *solver_grid.boundary_faces[id_triangle].GetElementDOFs());
			};
		}

		//считаем правую часть основного уравнения
		//пишем в MassMatrix.F
		{
			std::vector<Point<double>> Mu_prev;
			MassMatrix.MultiplicationMatrixVector(U_prev, Mu_prev);

			std::vector<Point<double>> Mu_prevprev;
			MassMatrix.MultiplicationMatrixVector(U_prevprev, Mu_prevprev);

			std::vector<Point<double>> Ku_prev;
			StiffnessMatrix.MultiplicationMatrixVector(U_prev, Ku_prev);

			for (int i = 0; i < MassMatrix.F.size(); i++)
			{
				MassMatrix.F[i] = Point<double>(0, 0, 0)
					+ Mu_prev[i] * 2.0 / (TIME_h * TIME_h)
					- Mu_prevprev[i] * 1.0 / (TIME_h * TIME_h)
					- Ku_prev[i]
					+ VolumeForceVector.F[i]
					+ SourseVector.F[i];
				MassMatrix.F[i] /= 1.0 / (TIME_h * TIME_h);
			}
		}

		//решаем СЛАУ
		{
			CSSD_Matrix<double, double> StiffnessMassMatrix_doubleSLAE;
			math::MakeCopyMatrix_A_into_B(MassMatrix, StiffnessMassMatrix_doubleSLAE);

			printf("First boundary conditions...");
			if (true) { 
				//обнуляем строки
				for (int id_vertex = 0; id_vertex < solver_grid.boundary_vertexes.size(); id_vertex++)
				{
					auto boundary = &(solver_grid.boundary_vertexes[id_vertex]);

					int id_type = boundary->id_type;
					Point<bool> is_take;
					Point<double> boundary_value = boundary->boundary_value(is_take, solver_grid.boundary_vertexes[id_vertex].GetIdNode(0));
					int global_id = boundary->GetDOFInLocalID(0);

					if (global_id < MassMatrix.GetMatrixSize())
					{
						auto enter_value = [&StiffnessMassMatrix_doubleSLAE](int value_id, double value) ->void
						{
							StiffnessMassMatrix_doubleSLAE.X[value_id] = value;
							StiffnessMassMatrix_doubleSLAE.F[value_id] = value;
							StiffnessMassMatrix_doubleSLAE.Diag[value_id] = 1;
							//Обнуляем строку
							for (int i = 0; i < StiffnessMassMatrix_doubleSLAE.A_down[value_id].size(); i++)
							{
								StiffnessMassMatrix_doubleSLAE.A_down[value_id][i] = 0;
							}
							for (int i = 0; i < StiffnessMassMatrix_doubleSLAE.A_up[value_id].size(); i++)
							{
								StiffnessMassMatrix_doubleSLAE.A_up[value_id][i] = 0;
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
#pragma omp parallel for schedule(dynamic) 
				for (int id_row = 0; id_row < StiffnessMassMatrix_doubleSLAE.GetMatrixSize(); id_row++)
				{
					if (id_row % 100 == 0)
					{
						printf("\tcurrent %d/%d\r", id_row, StiffnessMassMatrix_doubleSLAE.GetMatrixSize());
					}
					int iterator_in_boundary = 0;
					for (int jj = 0; jj < StiffnessMassMatrix_doubleSLAE.id_column_for_A_up[id_row].size(); jj++)
					{
						int id_column = StiffnessMassMatrix_doubleSLAE.id_column_for_A_up[id_row][jj];
						for (/*iterator_in_boundary = 0*/; iterator_in_boundary < solver_grid.boundary_vertexes.size(); iterator_in_boundary++)
						{
							if (solver_grid.boundary_vertexes[iterator_in_boundary].GetDOFInLocalID(0) * 3 + 0 == id_column)
							{
								Point<bool> is_take;
								Point<double> boundary_value = solver_grid.boundary_vertexes[iterator_in_boundary].boundary_value(is_take, solver_grid.boundary_vertexes[iterator_in_boundary].GetIdNode(0));

								if (is_take.x)
								{
									StiffnessMassMatrix_doubleSLAE.F[id_row] -= StiffnessMassMatrix_doubleSLAE.A_up[id_row][jj] * StiffnessMassMatrix_doubleSLAE.F[id_column];
									StiffnessMassMatrix_doubleSLAE.A_up[id_row][jj] = 0;
								}

								break;
							}
							if (solver_grid.boundary_vertexes[iterator_in_boundary].GetDOFInLocalID(0) * 3 + 1 == id_column)
							{
								Point<bool> is_take;
								Point<double> boundary_value = solver_grid.boundary_vertexes[iterator_in_boundary].boundary_value(is_take, solver_grid.boundary_vertexes[iterator_in_boundary].GetIdNode(0));

								if (is_take.y)
								{
									StiffnessMassMatrix_doubleSLAE.F[id_row] -= StiffnessMassMatrix_doubleSLAE.A_up[id_row][jj] * StiffnessMassMatrix_doubleSLAE.F[id_column];
									StiffnessMassMatrix_doubleSLAE.A_up[id_row][jj] = 0;
								}

								break;
							}
							if (solver_grid.boundary_vertexes[iterator_in_boundary].GetDOFInLocalID(0) * 3 + 2 == id_column)
							{
								Point<bool> is_take;
								Point<double> boundary_value = solver_grid.boundary_vertexes[iterator_in_boundary].boundary_value(is_take, solver_grid.boundary_vertexes[iterator_in_boundary].GetIdNode(0));

								if (is_take.z)
								{
									StiffnessMassMatrix_doubleSLAE.F[id_row] -= StiffnessMassMatrix_doubleSLAE.A_up[id_row][jj] * StiffnessMassMatrix_doubleSLAE.F[id_column];
									StiffnessMassMatrix_doubleSLAE.A_up[id_row][jj] = 0;
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
					for (int jj = 0; jj < StiffnessMassMatrix_doubleSLAE.id_column_for_A_down[id_row].size(); jj++)
					{
						int id_column = StiffnessMassMatrix_doubleSLAE.id_column_for_A_down[id_row][jj];
						for (/*iterator_in_boundary = 0*/; iterator_in_boundary < solver_grid.boundary_vertexes.size(); iterator_in_boundary++)
						{
							if (solver_grid.boundary_vertexes[iterator_in_boundary].GetDOFInLocalID(0) * 3 + 0 == id_column)
							{
								Point<bool> is_take;
								Point<double> boundary_value = solver_grid.boundary_vertexes[iterator_in_boundary].boundary_value(is_take, solver_grid.boundary_vertexes[iterator_in_boundary].GetIdNode(0));

								if (is_take.x)
								{
									StiffnessMassMatrix_doubleSLAE.F[id_row] -= StiffnessMassMatrix_doubleSLAE.A_down[id_row][jj] * StiffnessMassMatrix_doubleSLAE.F[id_column];
									StiffnessMassMatrix_doubleSLAE.A_down[id_row][jj] = 0;
								}

								break;
							}
							if (solver_grid.boundary_vertexes[iterator_in_boundary].GetDOFInLocalID(0) * 3 + 1 == id_column)
							{
								Point<bool> is_take;
								Point<double> boundary_value = solver_grid.boundary_vertexes[iterator_in_boundary].boundary_value(is_take, solver_grid.boundary_vertexes[iterator_in_boundary].GetIdNode(0));

								if (is_take.y)
								{
								StiffnessMassMatrix_doubleSLAE.F[id_row] -= StiffnessMassMatrix_doubleSLAE.A_down[id_row][jj] * StiffnessMassMatrix_doubleSLAE.F[id_column];
								StiffnessMassMatrix_doubleSLAE.A_down[id_row][jj] = 0;
								}

								break;
							}
							if (solver_grid.boundary_vertexes[iterator_in_boundary].GetDOFInLocalID(0) * 3 + 2 == id_column)
							{
								Point<bool> is_take;
								Point<double> boundary_value = solver_grid.boundary_vertexes[iterator_in_boundary].boundary_value(is_take, solver_grid.boundary_vertexes[iterator_in_boundary].GetIdNode(0));

								if (is_take.z)
								{
									StiffnessMassMatrix_doubleSLAE.F[id_row] -= StiffnessMassMatrix_doubleSLAE.A_down[id_row][jj] * StiffnessMassMatrix_doubleSLAE.F[id_column];
									StiffnessMassMatrix_doubleSLAE.A_down[id_row][jj] = 0;
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

			//начальное приближение
			if (id_STEP == 0)
			{
				for (int i = 0; i < StiffnessMassMatrix_doubleSLAE.GetMatrixSize() / 3; i++)
				{
					StiffnessMassMatrix_doubleSLAE.X[i * 3 + 0] = 1e-10;
					StiffnessMassMatrix_doubleSLAE.X[i * 3 + 1] = 1e-10;
					StiffnessMassMatrix_doubleSLAE.X[i * 3 + 2] = 1e-10;
				}
			}
			else {
				for (int i = 0; i < StiffnessMassMatrix_doubleSLAE.GetMatrixSize() / 3; i++)
				{
					StiffnessMassMatrix_doubleSLAE.X[i * 3 + 0] = U_prev[i].x;
					StiffnessMassMatrix_doubleSLAE.X[i * 3 + 1] = U_prev[i].y;
					StiffnessMassMatrix_doubleSLAE.X[i * 3 + 2] = U_prev[i].z;
				}
			}

			printf("Soluting SLAY... (%d)\n", StiffnessMassMatrix_doubleSLAE.GetMatrixSize());
			int MaxSize = StiffnessMassMatrix_doubleSLAE.GetMatrixSize();
			MaxSize = MaxSize / 10 < 1000 ? 1000 : MaxSize / 10;
			std::vector<double> best_solution;
			math::MakeCopyVector_A_into_B(StiffnessMassMatrix_doubleSLAE.X, best_solution);
			double current_residual = 1, best_residual = 1e+25;
			int MAX_STEPS = 4;
			int ii = 1;
			CSSD_Matrix<double, double> Precond;
			Precond.PrecondorSSOR(0.75, StiffnessMassMatrix_doubleSLAE);
			for (int i = 0; i <= MAX_STEPS; i++)
			{
				double needed_residual = critical_residual;// pow(10., -1 * (ii + 1));
				//current_residual /= 2.;
				printf_s("//---> I = %d/%d (full size %d) - needed residual %.2e\n", i, MAX_STEPS, StiffnessMassMatrix_doubleSLAE.GetMatrixSize(), needed_residual);

				StiffnessMassMatrix_doubleSLAE.print_logs = true;
				//current_residual = abs(MassMatrix_doubleSLAE.MSG(MaxSize, needed_residual));
				//current_residual = abs(MassMatrix_doubleSLAE.BCG_Stab2(MaxSize, needed_residual));
				current_residual = abs(StiffnessMassMatrix_doubleSLAE.MSG_PreconditioningSSOR(MaxSize, needed_residual, Precond));

				if (current_residual < needed_residual)
				{
					i = 0;
					ii++;
				}
				if (current_residual <= critical_residual)
				{
					best_residual = current_residual;
					math::MakeCopyVector_A_into_B(StiffnessMassMatrix_doubleSLAE.X, best_solution);
					break;
				}
				if (current_residual < best_residual)
				{
					best_residual = current_residual;
					math::MakeCopyVector_A_into_B(StiffnessMassMatrix_doubleSLAE.X, best_solution);
				}
			}
			math::MakeCopyVector_A_into_B(best_solution, StiffnessMassMatrix_doubleSLAE.X);
			printf_s("//---> BEST residual %.2e\n", best_residual);
			math::MakeCopyVector_A_into_B(StiffnessMassMatrix_doubleSLAE.X, U_curr);

			if (best_residual >= 1)
			{
				printf_s("There are problems in the solution\n");
				int a;
				scanf_s("%d", &a);
			}
		}

		//вывод сейсмограммы
		if(TIME_curr > TIME_L || true) {
			Point<double> U = solver_grid.GetSolutionInPoint(id_elem_for_out, point_for_out, U_curr);
			Uz_in_point.push_back(U.z);
			Time_for_Uz.push_back(TIME_curr);
		}
		
		//output solution
		//without deformations
		if (is_print_result) {
			printf_s("Print the mech result into .dat file... ");

			FILE* fout_tech;
			char name_u_tech[5000];
			sprintf_s(name_u_tech, "%s/U_deformation_%d_t%.2e.dat", result_directory, id_STEP, TIME_curr);
			fopen_s(&fout_tech, name_u_tech, "w");
			char name_in_file[1000];
			sprintf_s(name_in_file, "step%d_time%.2e", id_STEP, TIME_curr);
			std::vector<std::vector<char>> name_value(6);
			char name_v_tmp[6][100];
			sprintf_s(name_v_tmp[0], "sigma_Mises");
			sprintf_s(name_v_tmp[1], "sigma_zz");
			sprintf_s(name_v_tmp[2], "Ux");
			sprintf_s(name_v_tmp[3], "Uy");
			sprintf_s(name_v_tmp[4], "Uz");
			sprintf_s(name_v_tmp[5], "Vz");
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
			double sigma_inv_max = 0;
			int elem_sigma_max = 0;
			for (int i = 0; i < solver_grid.GetElementsCount(); i++)
			{
				auto element = solver_grid.GetElement(i);

				Point<double> Centr = element->GetWeightCentr();

				//Point<double> U = solver_grid.GetSolutionInPoint(i, Centr, U_curr);
				Point<Point<double>> dU = solver_grid.GetDerevativeFromSolutionInPoint(i, Centr, U_curr);
				auto Eps = solver_grid.GetStrainTensorFromSolutionInPoint(i, Centr, dU);
				auto Sigma = solver_grid.GetStressTensorFromSolutionInPoint(i, Centr, Eps);
				auto MisesSigma = solver_grid.GetVonMisesStress(Sigma);

				auto Eps_inv = sqrt((Eps.val[0][0] - Eps.val[1][1]) * (Eps.val[0][0] - Eps.val[1][1])
					+ (Eps.val[1][1] - Eps.val[2][2]) * (Eps.val[1][1] - Eps.val[2][2])
					+ (Eps.val[0][0] - Eps.val[2][2]) * (Eps.val[0][0] - Eps.val[2][2])
					+ 3 * (Eps.val[0][1] * Eps.val[1][0] + Eps.val[0][2] * Eps.val[2][0] + Eps.val[1][2] * Eps.val[2][1]) / 2.0) * sqrt(2.) / 3.;

				value[0][i] = MisesSigma;
				value[1][i] = Sigma.val[2][2];
			}
			value[2].resize(solver_grid.GetDOFsCount());
			value[3].resize(solver_grid.GetDOFsCount());
			value[4].resize(solver_grid.GetDOFsCount());
			value[5].resize(solver_grid.GetDOFsCount());
			for (int i = 0; i < value[3].size(); i++)
			{
				value[2][i] = U_curr[i].x;
				value[3][i] = U_curr[i].y;
				value[4][i] = U_curr[i].z;
				value[5][i] = (U_curr[i].z - U_prev[i].z) / TIME_h;
			}

			//solver_grid.MoveCoordinates(U_curr);
			solver_grid.printTecPlot3D_DiffDomains(fout_tech, value, name_value, name_in_file, TIME_curr);
			//solver_grid.ReMoveCoordinates(U_curr);
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
			sprintf_s(name_u_tech, "%s/U_plane_YZ_s%d_t%.2e.dat", result_directory, id_STEP, TIME_curr);
			fopen_s(&fout_tech, name_u_tech, "w");
			char name_in_file[1000];
			sprintf_s(name_in_file, "Time_%.4e_YZ", TIME_curr);
			std::vector<std::vector<char>> name_value(6);
			char name_v_tmp[6][100];
			sprintf_s(name_v_tmp[0], "sigma_Mises");
			sprintf_s(name_v_tmp[1], "sigma_zz");
			sprintf_s(name_v_tmp[2], "Ux");
			sprintf_s(name_v_tmp[3], "Uy");
			sprintf_s(name_v_tmp[4], "Uz");
			sprintf_s(name_v_tmp[5], "Vz");
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
					Point<double> U_pr = solver_grid.GetSolutionInPoint(id_elem, Centr, U_prev);
					Point<Point<double>> dU = solver_grid.GetDerevativeFromSolutionInPoint(id_elem, Centr, U_curr);
					auto Eps = solver_grid.GetStrainTensorFromSolutionInPoint(id_elem, Centr, dU);
					auto Sigma = solver_grid.GetStressTensorFromSolutionInPoint(id_elem, Centr, Eps);
					auto MisesSigma = solver_grid.GetVonMisesStress(Sigma);

					auto Eps_inv = sqrt((Eps.val[0][0] - Eps.val[1][1]) * (Eps.val[0][0] - Eps.val[1][1])
						+ (Eps.val[1][1] - Eps.val[2][2]) * (Eps.val[1][1] - Eps.val[2][2])
						+ (Eps.val[0][0] - Eps.val[2][2]) * (Eps.val[0][0] - Eps.val[2][2])
						+ 3 * (Eps.val[0][1] * Eps.val[1][0] + Eps.val[0][2] * Eps.val[2][0] + Eps.val[1][2] * Eps.val[2][1]) / 2.0) * sqrt(2.) / 3.;

					value[0][i] = MisesSigma;
					value[1][i] = Sigma.val[2][2];

					value[2][i] = U.x;
					value[3][i] = U.y;
					value[4][i] = U.z;
					value[5][i] = (U.z - U_pr.z) / TIME_h;
				}
			}

			grid_YZ_plane.printTecPlot3D(fout_tech, value, name_value, name_in_file, TIME_curr);
			fclose(fout_tech);
		}
		//Uz in point
		if (is_print_result && TIME_curr > TIME_L)
		{
			FILE* fout_Uz;
			fopen_s(&fout_Uz, name_uz, "a");
			for (int i = 0; i < Uz_in_point.size(); i++)
			{
				fprintf_s(fout_Uz, "%.10e %.10e\n", Time_for_Uz[i], Uz_in_point[i]);
			}
			Time_for_Uz.clear();
			Uz_in_point.clear();
			fclose(fout_Uz);
		}

		printf_s("\t complite\n");

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
	}
}
void ElastodynamicsProblem_ExplicitSimple_fast(char properties_file[1000])
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
	//char properties_file[1000] = { "D:/Elastodynamic/homocyl/67k/param_for_solver.txt" };
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


	{
		char _line[1000];
		std::vector<int> val;

		math::ReadNonEmptyLine_forNumbers(f_properties, _line);
		math::ParserStringToVectorInt(_line, val, " ");
		if (val[0] == 1)
			is_print_logFile = true;


		math::ReadNonEmptyLine_forNumbers(f_properties, _line);
		math::ParserStringToVectorInt(_line, val, " ");
		math::NUM_THREADS = val[0];
	}

	math::ReadNonEmptyLine(f_properties, mesh_directory);
	math::ReadNonEmptyLine(f_properties, base_result_directory);
	math::SimpleGrid geo_grid; //input
	std::vector<std::vector<std::vector<int>>> boundary_faces;
	geo_grid.ReadFromNVTR(mesh_directory, 4);
	geo_grid.ReadFacesBoundaryNVTR(mesh_directory, boundary_faces);
	
	struct Receiver
	{
		Point<double> point;
		std::vector<Point<double>> U_in_point;
		int id_node;
		char file_name[1000];
	};
	std::vector<Receiver> receivers;
	std::vector<double> receivers_times;
	{
		char _line[1000];
		math::ReadNonEmptyLine_forNumbers(f_properties, _line);
		std::vector<double> val;
		math::ParserStringToVectorDouble(_line, val, " ");

		int N = (int)val[0];

		receivers.resize(N);
		for (int i = 0; i < N; i++)
		{
			receivers[i].point.x = val[1 + 0 + i * 3];
			receivers[i].point.y = val[1 + 1 + i * 3];
			receivers[i].point.z = val[1 + 2 + i * 3];
		}
	}
	Point<double> crack_top;
	Point<double> crack_bottom;
	{
		char _line[1000];
		math::ReadNonEmptyLine_forNumbers(f_properties, _line);
		std::vector<double> val;
		math::ParserStringToVectorDouble(_line, val, " ");

		crack_bottom.x = val[0];
		crack_bottom.y = val[1];
		crack_bottom.z = val[2];

		crack_top.x = val[3];
		crack_top.y = val[4];
		crack_top.z = val[5];
	}

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
	struct Gauss_source {
		bool is_used = false;
		int target_boundary=0;
		int degree=0;
		double tay=0;
		double t0=0;
		double power=0;
		Point<double> direction;

		Point<double> value(double Time_curr)
		{
			double arg = (Time_curr - t0) / tay;
			return direction * power * exp(-1.0 * pow(arg, degree));
		}
	} gauss_source;
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
			int id_boundary = math::ParserCharToInt(line[0]);
			
			if (line[0] == '-')
				id_boundary = math::ParserCharToInt(line[1]) * -1;

			id_boundaries[i] = abs(id_boundary);
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

					//special boundary
					if (id_boundary < 0)
					{
						gauss_source.is_used = true;
						gauss_source.power = val[1];
						gauss_source.direction.x = val[2];
						gauss_source.direction.y = val[3];
						gauss_source.direction.z = val[4];
						gauss_source.degree = (int)val[5];
						gauss_source.t0 = val[6];
						gauss_source.tay = val[7];
						gauss_source.target_boundary = i;
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

	{
		N_domain++;
		double _E = 0;
		double _v = 0;
		double _rpho = 0;
		_materials.push_back(_material(_E, _v, _rpho));

		for (int i = 0; i < geo_grid.nvtr.size(); i++)
		{
			Point<double> centr;
			for (int j = 0; j < geo_grid.nvtr[i].size(); j++)
			{
				centr += geo_grid.xyz[geo_grid.nvtr[i][j]] / geo_grid.nvtr[i].size();
			}
			if (crack_bottom.x < centr.x && centr.x < crack_top.x &&
				crack_bottom.y < centr.y && centr.y < crack_top.y &&
				crack_bottom.z < centr.z && centr.z < crack_top.z)
			{
				geo_grid.nvkat[i] = N_domain - 1;
			}
		}
	}

	fclose(f_properties);

	wchar_t _tmp_wc[1000];
	math::Char_To_Wchar_t(base_result_directory, _tmp_wc, 1000);
	CreateDirectory((LPCTSTR)_tmp_wc, NULL);

	bool res = false;
	double TIME_h = step_size;
	if (is_STATIONARY) {
		end_iteration = 1;
		start_iteration = 0;
	}
	std::vector<Point<double>> U_prev(geo_grid.xyz.size()), U_prevprev(geo_grid.xyz.size()), U_curr(geo_grid.xyz.size());
	math::InitializationVector(U_curr, 0);
	math::InitializationVector(U_prev, 0);
	math::InitializationVector(U_prevprev, 0);

	CSSD_Matrix<Tensor2Rank3D, Point<double>> StiffnessMatrix;
	CSSD_Matrix<Tensor2Rank3D, Point<double>> MassMatrix;
	CSSD_Matrix<double, double> MassMatrix_doubleSLAE;
	CSSD_Matrix<double, double> Precond;
	CSSD_Matrix<Tensor2Rank3D, Point<double>> VolumeForceVector;
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
	MassMatrix.SetMatrix(StiffnessMatrix);
	VolumeForceVector.F.resize(StiffnessMatrix.F.size());
	printf_s("\t\tcomplite\n");

	std::vector<int> clear_vertexes;

	//строим базовые матрицы жесткости/массы/вектор силы тяжести
	//матрица массы идет с коэффициентом rpho
	printf_s("================= Create matrix ================\n");
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
			return solver_grid.GetDomain(solver_grid.GetElement(elem)->GetIdDomain())->forThermal.rpho;
		};
		std::function<Point<double>(int, Point<double>)> VolumeForсe = [&](int elem, Point<double> X)->Point<double>
		{
			return Point<double>(0, 0, 0);
		};

		printf("Matrix assembling...\n");
		std::vector<DenseMatrix<Tensor2Rank3D, Point<double>>> local_SLAE_stiffness(solver_grid.GetElementsCount());
		std::vector<DenseMatrix<Tensor2Rank3D, Point<double>>> local_SLAE_mass(solver_grid.GetElementsCount());
		std::vector< std::vector<Point<double>>> local_force(solver_grid.GetElementsCount());
		omp_set_num_threads(math::NUM_THREADS);
#pragma omp parallel for schedule(dynamic)
		for (int id_elem = 0; id_elem < solver_grid.GetElementsCount(); id_elem++)
		{
			if (id_elem % 1000 == 0)
				printf("Solve element[%d]\r", id_elem);
			auto element = solver_grid.GetElement(id_elem);

			if (element->GetIdDomain() != solver_grid.GetDomainsCount() - 1)
			{
				std::function<std::vector<std::vector<double>>(Point<double>)> D = [&](Point<double> X) {return StiffnessCoef(id_elem, X); };
				std::function<double(Point<double>)> M = [&](Point<double> X) {return MassCoef(id_elem, X); };
				std::function<Point<double>(Point<double>)> F = [&](Point<double> X) {return VolumeForсe(id_elem, X); };

				element->SolveLocalMatrix(local_SLAE_stiffness[id_elem], D);
				element->SolveMassMatrix(local_SLAE_mass[id_elem], M);
				element->SolveRightSide(local_force[id_elem], F);
			}
			else
			{
				local_SLAE_mass[id_elem].SetSize(element->GetDOFsCount());
				local_SLAE_stiffness[id_elem].SetSize(element->GetDOFsCount());
				local_force[id_elem].resize(element->GetDOFsCount());
			}
		}
		for (int id_elem = 0; id_elem < solver_grid.GetElementsCount(); id_elem++)
		{
			if (id_elem % 1000 == 0)
				printf("Add the local matrix of element[%d]\r", id_elem);
			StiffnessMatrix.SummPartOfMatrix(local_SLAE_stiffness[id_elem], *solver_grid.GetElementDOFs(id_elem));
			MassMatrix.SummPartOfMatrix(local_SLAE_mass[id_elem], *solver_grid.GetElementDOFs(id_elem));
			VolumeForceVector.SummPartOfVector(local_force[id_elem], *solver_grid.GetElementDOFs(id_elem));
		}
		printf_s("                                                                                    \r");
		printf_s("\t\tcomplite\n");

		//учитываем первые краевые ТОЛЬКО в матрице массы БЕЗ симметризации
		if (false) {
			//обнуляем строки
			for (int id_vertex = 0; id_vertex < solver_grid.boundary_vertexes.size(); id_vertex++)
			{
				auto boundary = &(solver_grid.boundary_vertexes[id_vertex]);

				int id_type = boundary->id_type;
				Point<bool> is_take;
				Point<double> boundary_value = boundary->boundary_value(is_take, solver_grid.boundary_vertexes[id_vertex].GetIdNode(0));
				int global_id = boundary->GetDOFInLocalID(0);

				auto enter_boundary = [global_id, &boundary_value](int position, CSSD_Matrix<Tensor2Rank3D, Point<double>>& global_SLAE) {
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
				};
				if (is_take.x == true) enter_boundary(0, MassMatrix);
				if (is_take.y == true) enter_boundary(1, MassMatrix);
				if (is_take.z == true) enter_boundary(2, MassMatrix);
			}
		}
		if (true) {
			//большое число на диагональ
			for (int id_vertex = 0; id_vertex < solver_grid.boundary_vertexes.size(); id_vertex++)
			{
				auto boundary = &(solver_grid.boundary_vertexes[id_vertex]);

				int id_type = boundary->id_type;
				Point<bool> is_take;
				Point<double> boundary_value = boundary->boundary_value(is_take, solver_grid.boundary_vertexes[id_vertex].GetIdNode(0));
				int global_id = boundary->GetDOFInLocalID(0);

				auto enter_boundary = [global_id, &boundary_value](int position, CSSD_Matrix<Tensor2Rank3D, Point<double>>& global_SLAE) {
					global_SLAE.Diag[global_id].val[position][position] = 1e+30;
				};
				if (is_take.x == true) enter_boundary(0, MassMatrix);
				if (is_take.y == true) enter_boundary(1, MassMatrix);
				if (is_take.z == true) enter_boundary(2, MassMatrix);
			}
		}
		for (int i = 0; i < MassMatrix.Diag.size(); i++)
		{
			if (abs(MassMatrix.Diag[i].val[0][0]) < 1e-12 || abs(MassMatrix.Diag[i].val[1][1]) < 1e-12 || abs(MassMatrix.Diag[i].val[2][2]) < 1e-12)
			{
				MassMatrix.Diag[i].InitializationAsI();
				clear_vertexes.push_back(i);
			}
		}

		//переводим матрицу массы в действительную
		math::MakeCopyMatrix_A_into_B(MassMatrix, MassMatrix_doubleSLAE);
		Precond.PrecondorSSOR(0.75, MassMatrix_doubleSLAE);
	}

	///-------------------
	//Point<double> point_for_out(0, 0, 0.5 - 0.01);
	
	double TIME_L;
	{
		double v = solver_grid.GetDomain(0)->forMech.GetLongitudinalWaveVelocity_Vp(solver_grid.GetDomain(0)->forThermal.rpho);
		double L = 0.5;
		TIME_L = L / v;

		wchar_t _tmp_wc[1000];
		math::Char_To_Wchar_t(base_result_directory, _tmp_wc, 1000);
		CreateDirectory((LPCTSTR)_tmp_wc, NULL);

		for (int i = 0; i < receivers.size(); i++)
		{
			double len;
			int id_elem_for_out = solver_grid.GetNearestElementID(receivers[i].point, len);

			if (id_elem_for_out >= 0)
			{
				double min_len = 1e+50;
				int target_point = -1;
				auto elem = solver_grid.GetElement(id_elem_for_out);
				for (int j = 0; j < elem->GetNodesCount(); j++)
				{
					double _len = math::SolveLengthVector(receivers[i].point, elem->GetNode(j));
					if (_len < min_len)
					{
						min_len = _len;
						target_point = elem->GetIdNode(j);
					}
				}
				receivers[i].id_node = target_point;
			}
			else
			{
				return;
			}

			FILE* fout_Uz;
			sprintf_s(receivers[i].file_name, "%s/Uz_in_point_%d.txt", base_result_directory, i);
			fopen_s(&fout_Uz, receivers[i].file_name, "w");
			fprintf_s(fout_Uz, "Point: (%.2e, %.2e, %.2e)\n", receivers[i].point.x, receivers[i].point.y, receivers[i].point.z);
			fprintf_s(fout_Uz, "Id vertex: %d\n", receivers[i].id_node);
			fprintf_s(fout_Uz, "T0(Z): %.5e sec\n", TIME_L);
			fclose(fout_Uz);
		}
	}
	///-------------------

	double TIME_curr = TIME_h;
	for (int id_STEP = start_iteration; id_STEP <= end_iteration; id_STEP++)
	{
		TIME_curr += TIME_h;
		printf_s("\n================= Start solution of %d STEP (time = %.2e) ================\n", id_STEP, TIME_curr);
		bool is_print_result = false;
		char result_directory[1000];
		char result_directory_XZ[1000];
		char result_directory_XY[1000];
		char result_directory_YZ[1000];
		char result_directory_save[1000];
		if (id_STEP % step_for_out == 0)
		{
			is_print_result = true;
			sprintf_s(result_directory, sizeof(result_directory), "%s", base_result_directory);
			//sprintf_s(result_directory, sizeof(result_directory), "%s/STEP_%d_t=%.2e", base_result_directory, id_STEP, TIME_curr);
			wchar_t _tmp_wc[1000];
			math::Char_To_Wchar_t(result_directory, _tmp_wc, 1000);
			CreateDirectory((LPCTSTR)_tmp_wc, NULL);

			
			sprintf_s(result_directory_XY, sizeof(result_directory_XY), "%s/Plane_XY", base_result_directory);
			math::Char_To_Wchar_t(result_directory_XY, _tmp_wc, 1000);
			CreateDirectory((LPCTSTR)_tmp_wc, NULL);
			sprintf_s(result_directory_XZ, sizeof(result_directory_XZ), "%s/Plane_XZ", base_result_directory);
			math::Char_To_Wchar_t(result_directory_XZ, _tmp_wc, 1000);
			CreateDirectory((LPCTSTR)_tmp_wc, NULL);
			sprintf_s(result_directory_YZ, sizeof(result_directory_YZ), "%s/Plane_YZ", base_result_directory);
			math::Char_To_Wchar_t(result_directory_YZ, _tmp_wc, 1000);
			CreateDirectory((LPCTSTR)_tmp_wc, NULL);
			sprintf_s(result_directory_save, sizeof(result_directory_save), "%s/Saves", base_result_directory);
			math::Char_To_Wchar_t(result_directory_save, _tmp_wc, 1000);
			CreateDirectory((LPCTSTR)_tmp_wc, NULL);
		}

		//переопределяем источник и строим вектор в правую часть (через вторые краевые)
		CSSD_Matrix<Tensor2Rank3D, Point<double>> SourseVector;
		SourseVector.F.resize(StiffnessMatrix.F.size());
		if (true) {
			//Point<double> power(0, 0, -10e+6 / (2 * M_PI * 0.004 * 0.004));
			//double tay_plus = 0.0001; //половина длины импульса
			//double tay_0 = 0.0001; //середина импульса
			//double degree = 2; //степень функции
			//std::function<Point<double>(Point<double>)> new_sourse_value = [TIME_curr, power, tay_0, degree, tay_plus, TIME_h](Point<double> X) -> Point<double>
			//{
			//	//double sourse = exp(-pow(TIME_curr - tay_0, degree) / pow(tay_plus, degree));
			//	//return power;
			//	//return Point<double>(power.x * sourse, power.y * sourse, power.z * sourse);
			//	double sourse = TIME_curr < 14e-7 ? 10e+6 : 0;
			//	return Point<double>(0 * sourse, 0 * sourse, -1 * sourse /** TIME_h / 14.0e-7*/);
			//};
			//std::function<Point<double>(Point<double>)> new_sourse_value_sin = [TIME_curr, power, &solver_grid, TIME_h, &sourse_power](Point<double> X) -> Point<double>
			//{
			//	//double sourse = exp(-pow(TIME_curr - tay_0, degree) / pow(tay_plus, degree));
			//	//return power;
			//	//return Point<double>(power.x * sourse, power.y * sourse, power.z * sourse);
			//	double v = solver_grid.GetDomain(0)->forMech.GetLongitudinalWaveVelocity_Vp(solver_grid.GetDomain(0)->forThermal.rpho);
			//	double L = 0.5;
			//	double t0 = L / v;
			//	double w = 2 * M_PI / t0 * 9.5;
			//	double alpha = 2.5e+4;
			//
			//	double sourse = exp(-1 * alpha * TIME_curr) * sin(w * TIME_curr);
			//	if (TIME_curr > t0) sourse = 0;
			//	return sourse_power * sourse;
			//};
			//std::function<Point<double>(Point<double>)> new_sourse_value_gauss8 = [TIME_curr, power, &solver_grid, TIME_h, &sourse_power](Point<double> X) -> Point<double>
			//{
			//	//double sourse = exp(-pow(TIME_curr - tay_0, degree) / pow(tay_plus, degree));
			//	//return power;
			//	//return Point<double>(power.x * sourse, power.y * sourse, power.z * sourse);
			//	double v = solver_grid.GetDomain(0)->forMech.GetLongitudinalWaveVelocity_Vp(solver_grid.GetDomain(0)->forThermal.rpho);
			//	double L = 0.5;
			//	double t0 = L / v;
			//	double tay = t0/50.0;
			//	double gauss_0 = 0;
			//	double arg = (TIME_curr - gauss_0) / tay;
			//	double sourse = exp(-1.0 * arg * arg * arg * arg * arg * arg * arg * arg);
			//	if (TIME_curr > t0) sourse = 0;
			//	return sourse_power * sourse;
			//};
			//std::function<Point<double>(Point<double>)> round_cond = [TIME_curr, power, &solver_grid, TIME_h, round_power](Point<double> X) -> Point<double>
			//{
			//	//return Point<double>(X.x * round_power.x, X.y * round_power.y, 0);
			//	return round_power;
			//};
			//if (second_boundary.size() > 0)
			//{
			//	//second_boundary[0].value = new_sourse_value_sin;
			//	second_boundary[0].value = new_sourse_value_gauss8;
			//
			//	if(second_boundary.size() == 2)
			//		second_boundary[1].value = round_cond;
			//	//second_boundary[2].value = new_sourse_value;
			//
			//	for (int i = 0; i < solver_grid.boundary_faces.size(); i++)
			//	{
			//		solver_grid.boundary_faces[i].boundary_value = second_boundary[solver_grid.boundary_faces[i].id_type].value;
			//	}
			//}

			std::function<Point<double>(Point<double>)> time_source_condition = [TIME_curr, &gauss_source](Point<double> X) -> Point<double>
			{
				return gauss_source.value(TIME_curr);
			};
			if (gauss_source.is_used)
			{
				second_boundary[gauss_source.target_boundary].value = time_source_condition;
				for (int i = 0; i < solver_grid.boundary_faces.size(); i++)
				{
					solver_grid.boundary_faces[i].boundary_value = second_boundary[solver_grid.boundary_faces[i].id_type].value;
				}
			}

			for (int id_triangle = 0; id_triangle < solver_grid.boundary_faces.size(); id_triangle++)
			{
				std::vector<Point<double>> local_vector_SLAE;
				solver_grid.boundary_faces[id_triangle].SolveLocalBoundaryVector(local_vector_SLAE, solver_grid.boundary_faces[id_triangle].boundary_value);
				SourseVector.SummPartOfVector(local_vector_SLAE, *solver_grid.boundary_faces[id_triangle].GetElementDOFs());
			};
		}

		auto SolveSLAE = [&MassMatrix_doubleSLAE, &solver_grid, &Precond, &clear_vertexes](std::vector<double>& F, std::vector<double>& result, double critical_residual) -> void
		{
			//краевые в правую часть + симметризация 
			if (false) {
				//добавляем значения в правую часть
				for (int id_vertex = 0; id_vertex < solver_grid.boundary_vertexes.size(); id_vertex++)
				{
					auto boundary = &(solver_grid.boundary_vertexes[id_vertex]);

					int id_type = boundary->id_type;
					Point<bool> is_take;
					Point<double> boundary_value = boundary->boundary_value(is_take, solver_grid.boundary_vertexes[id_vertex].GetIdNode(0));
					int global_id = boundary->GetDOFInLocalID(0);

					if (global_id < MassMatrix_doubleSLAE.GetMatrixSize() / 3)
					{
						auto enter_value = [&result, &F](int value_id, double value) ->void
						{
							result[value_id] = value;
							F[value_id] = value;
							return;
						};

						if (is_take.x) enter_value(global_id * 3 + 0, boundary_value.x);
						if (is_take.y) enter_value(global_id * 3 + 1, boundary_value.y);
						if (is_take.z) enter_value(global_id * 3 + 2, boundary_value.z);
					}

				}
				//симметризация
				omp_set_num_threads(math::NUM_THREADS);
#pragma omp parallel for schedule(dynamic) 
				for (int id_row = 0; id_row < MassMatrix_doubleSLAE.GetMatrixSize(); id_row++)
				{
					if (id_row % 100 == 0)
					{
						printf("\tcurrent %d/%d\r", id_row, MassMatrix_doubleSLAE.GetMatrixSize());
					}
					int iterator_in_boundary = 0;
					for (int jj = 0; jj < MassMatrix_doubleSLAE.id_column_for_A_up[id_row].size(); jj++)
					{
						int id_column = MassMatrix_doubleSLAE.id_column_for_A_up[id_row][jj];
						for (/*iterator_in_boundary = 0*/; iterator_in_boundary < solver_grid.boundary_vertexes.size(); iterator_in_boundary++)
						{
							if (solver_grid.boundary_vertexes[iterator_in_boundary].GetDOFInLocalID(0) * 3 + 0 == id_column)
							{
								Point<bool> is_take;
								Point<double> boundary_value = solver_grid.boundary_vertexes[iterator_in_boundary].boundary_value(is_take, solver_grid.boundary_vertexes[iterator_in_boundary].GetIdNode(0));

								if (is_take.x)
								{
									F[id_row] -= MassMatrix_doubleSLAE.A_up[id_row][jj] * F[id_column];
									MassMatrix_doubleSLAE.A_up[id_row][jj] = 0;
								}

								break;
							}
							if (solver_grid.boundary_vertexes[iterator_in_boundary].GetDOFInLocalID(0) * 3 + 1 == id_column)
							{
								Point<bool> is_take;
								Point<double> boundary_value = solver_grid.boundary_vertexes[iterator_in_boundary].boundary_value(is_take, solver_grid.boundary_vertexes[iterator_in_boundary].GetIdNode(0));

								if (is_take.y)
								{
									F[id_row] -= MassMatrix_doubleSLAE.A_up[id_row][jj] * F[id_column];
									MassMatrix_doubleSLAE.A_up[id_row][jj] = 0;
								}

								break;
							}
							if (solver_grid.boundary_vertexes[iterator_in_boundary].GetDOFInLocalID(0) * 3 + 2 == id_column)
							{
								Point<bool> is_take;
								Point<double> boundary_value = solver_grid.boundary_vertexes[iterator_in_boundary].boundary_value(is_take, solver_grid.boundary_vertexes[iterator_in_boundary].GetIdNode(0));

								if (is_take.z)
								{
									F[id_row] -= MassMatrix_doubleSLAE.A_up[id_row][jj] * F[id_column];
									MassMatrix_doubleSLAE.A_up[id_row][jj] = 0;
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
					for (int jj = 0; jj < MassMatrix_doubleSLAE.id_column_for_A_down[id_row].size(); jj++)
					{
						int id_column = MassMatrix_doubleSLAE.id_column_for_A_down[id_row][jj];
						for (/*iterator_in_boundary = 0*/; iterator_in_boundary < solver_grid.boundary_vertexes.size(); iterator_in_boundary++)
						{
							if (solver_grid.boundary_vertexes[iterator_in_boundary].GetDOFInLocalID(0) * 3 + 0 == id_column)
							{
								Point<bool> is_take;
								Point<double> boundary_value = solver_grid.boundary_vertexes[iterator_in_boundary].boundary_value(is_take, solver_grid.boundary_vertexes[iterator_in_boundary].GetIdNode(0));

								if (is_take.x)
								{
									F[id_row] -= MassMatrix_doubleSLAE.A_down[id_row][jj] * F[id_column];
									MassMatrix_doubleSLAE.A_down[id_row][jj] = 0;
								}

								break;
							}
							if (solver_grid.boundary_vertexes[iterator_in_boundary].GetDOFInLocalID(0) * 3 + 1 == id_column)
							{
								Point<bool> is_take;
								Point<double> boundary_value = solver_grid.boundary_vertexes[iterator_in_boundary].boundary_value(is_take, solver_grid.boundary_vertexes[iterator_in_boundary].GetIdNode(0));

								if (is_take.y)
								{
									F[id_row] -= MassMatrix_doubleSLAE.A_down[id_row][jj] * F[id_column];
									MassMatrix_doubleSLAE.A_down[id_row][jj] = 0;
								}

								break;
							}
							if (solver_grid.boundary_vertexes[iterator_in_boundary].GetDOFInLocalID(0) * 3 + 2 == id_column)
							{
								Point<bool> is_take;
								Point<double> boundary_value = solver_grid.boundary_vertexes[iterator_in_boundary].boundary_value(is_take, solver_grid.boundary_vertexes[iterator_in_boundary].GetIdNode(0));

								if (is_take.z)
								{
									F[id_row] -= MassMatrix_doubleSLAE.A_down[id_row][jj] * F[id_column];
									MassMatrix_doubleSLAE.A_down[id_row][jj] = 0;
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
			//краевые в правую часть большим числом 
			if (true) {
				//добавляем значения в правую часть
				for (int id_vertex = 0; id_vertex < solver_grid.boundary_vertexes.size(); id_vertex++)
				{
					auto boundary = &(solver_grid.boundary_vertexes[id_vertex]);

					int id_type = boundary->id_type;
					Point<bool> is_take;
					Point<double> boundary_value = boundary->boundary_value(is_take, solver_grid.boundary_vertexes[id_vertex].GetIdNode(0));
					int global_id = boundary->GetDOFInLocalID(0);

					if (global_id < MassMatrix_doubleSLAE.GetMatrixSize() / 3)
					{
						auto enter_value = [&result, &F, &MassMatrix_doubleSLAE](int value_id, double value) ->void
						{
							result[value_id] = value;
							F[value_id] = value * MassMatrix_doubleSLAE.Diag[value_id];
							return;
						};

						if (is_take.x) enter_value(global_id * 3 + 0, boundary_value.x);
						if (is_take.y) enter_value(global_id * 3 + 1, boundary_value.y);
						if (is_take.z) enter_value(global_id * 3 + 2, boundary_value.z);
					}

				}
				for (int i = 0; i < clear_vertexes.size(); i++)
				{
					result[clear_vertexes[i]] = 0;
					F[clear_vertexes[i]] = 0;
				}
			}

			printf("Soluting SLAY... (%d)\n", MassMatrix_doubleSLAE.GetMatrixSize());
			int MaxSize = MassMatrix_doubleSLAE.GetMatrixSize();
			MaxSize = MaxSize / 10 < 1000 ? 1000 : MaxSize / 10;
			std::vector<double> best_solution;
			math::MakeCopyVector_A_into_B(result, best_solution);
			double current_residual = 1, best_residual = 1e+25;
			int MAX_STEPS = 4;
			int ii = 1;

			for (int i = 0; i <= MAX_STEPS; i++)
			{
				printf_s("//---> I = %d/%d (full size %d) - needed residual %.2e\n", i, MAX_STEPS, MassMatrix_doubleSLAE.GetMatrixSize(), critical_residual);

				MassMatrix_doubleSLAE.print_logs = true;
				current_residual = abs(MassMatrix_doubleSLAE.MSG_PreconditioningSSOR(F, result, MaxSize, critical_residual, Precond));

				if (current_residual <= critical_residual)
				{
					best_residual = current_residual;
					math::MakeCopyVector_A_into_B(result, best_solution);
					break;
				}
				if (current_residual < best_residual)
				{
					best_residual = current_residual;
					math::MakeCopyVector_A_into_B(result, best_solution);
				}
			}
			math::MakeCopyVector_A_into_B(best_solution, result);

			if (best_residual >= 1)
			{
				printf_s("There are problems in the solution\n");
				int a;
				scanf_s("%d", &a);
			}
		};

		std::vector<double> F_slae(U_prev.size() * 3);

		//считаем правую часть основного уравнения
		//пишем в MassMatrix.F
		{
			std::vector<Point<double>> Mu_prev;
			MassMatrix.MultiplicationMatrixVector(U_prev, Mu_prev);

			std::vector<Point<double>> Mu_prevprev;
			MassMatrix.MultiplicationMatrixVector(U_prevprev, Mu_prevprev);

			std::vector<Point<double>> Ku_prev;
			StiffnessMatrix.MultiplicationMatrixVector(U_prev, Ku_prev);

			for (int i = 0; i < MassMatrix.F.size(); i++)
			{
				Point<double> res = Point<double>(0, 0, 0)
					+ Mu_prev[i] * 2.0 / (TIME_h * TIME_h)
					- Mu_prevprev[i] * 1.0 / (TIME_h * TIME_h)
					- Ku_prev[i]
					+ VolumeForceVector.F[i]
					+ SourseVector.F[i];
				MassMatrix.F[i] /= 1.0 / (TIME_h * TIME_h);
				F_slae[i * 3 + 0] = res.x / (1.0 / (TIME_h * TIME_h));
				F_slae[i * 3 + 1] = res.y / (1.0 / (TIME_h * TIME_h));
				F_slae[i * 3 + 2] = res.z / (1.0 / (TIME_h * TIME_h));
			}
		}
		//Решаем СЛАУ
		if (id_STEP == 0)
			math::InitializationVector(MassMatrix_doubleSLAE.X, 1e-10);
		else
			math::MakeCopyVector_A_into_B(U_prev, MassMatrix_doubleSLAE.X);
		SolveSLAE(F_slae, MassMatrix_doubleSLAE.X, 1e-10);
		math::MakeCopyVector_A_into_B(MassMatrix_doubleSLAE.X, U_curr);

		//вывод сейсмограммы
		if (TIME_curr > TIME_L || true) 
		{
			receivers_times.push_back(TIME_curr);
			for (int r = 0; r < receivers.size(); r++)
			{
				receivers[r].U_in_point.push_back(U_curr[receivers[r].id_node]);
			}
		}

		//save reserve solution
		if (id_STEP % (step_for_out * 10) == 0)
		{
			FILE* fout_save;
			char name_u_save[5000];

			sprintf_s(name_u_save, "%s/U_s%d_t%.2e.dat", result_directory_save, id_STEP, TIME_curr);
			fopen_s(&fout_save, name_u_save, "w");
			for (int i = 0; i < U_curr.size(); i++)
			{
				fprintf_s(fout_save, "%.16e %.16e %.16e\n", U_curr[i].x, U_curr[i].y, U_curr[i].z);
			}
			fclose(fout_save);

			sprintf_s(name_u_save, "%s/U_s%d_t%.2e.dat", result_directory_save, id_STEP - 1, TIME_curr - TIME_h);
			fopen_s(&fout_save, name_u_save, "w");
			for (int i = 0; i < U_prev.size(); i++)
			{
				fprintf_s(fout_save, "%.16e %.16e %.16e\n", U_prev[i].x, U_prev[i].y, U_prev[i].z);
			}
			fclose(fout_save);

			sprintf_s(name_u_save, "%s/U_s%d_t%.2e.dat", result_directory_save, id_STEP - 2, TIME_curr - TIME_h * 2);
			fopen_s(&fout_save, name_u_save, "w");
			for (int i = 0; i < U_prevprev.size(); i++)
			{
				fprintf_s(fout_save, "%.16e %.16e %.16e\n", U_prevprev[i].x, U_prevprev[i].y, U_prevprev[i].z);
			}
			fclose(fout_save);
		}

		//output solution
		//without deformations
		if (is_print_result && id_STEP == start_iteration) {
			printf_s("Print the mech result into .dat file... ");

			FILE* fout_tech;
			char name_u_tech[5000];
			sprintf_s(name_u_tech, "%s/U_deformation_%d_t%.2e.dat", result_directory, id_STEP, TIME_curr);
			fopen_s(&fout_tech, name_u_tech, "w");
			char name_in_file[1000];
			sprintf_s(name_in_file, "step%d_time%.2e", id_STEP, TIME_curr);
			std::vector<std::vector<char>> name_value(6);
			char name_v_tmp[6][100];
			sprintf_s(name_v_tmp[0], "sigma_Mises");
			sprintf_s(name_v_tmp[1], "sigma_zz");
			sprintf_s(name_v_tmp[2], "Ux");
			sprintf_s(name_v_tmp[3], "Uy");
			sprintf_s(name_v_tmp[4], "Uz");
			sprintf_s(name_v_tmp[5], "Vz");
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
			double sigma_inv_max = 0;
			int elem_sigma_max = 0;
			for (int i = 0; i < solver_grid.GetElementsCount(); i++)
			{
				auto element = solver_grid.GetElement(i);

				Point<double> Centr = element->GetWeightCentr();

				//Point<double> U = solver_grid.GetSolutionInPoint(i, Centr, U_curr);
				Point<Point<double>> dU = solver_grid.GetDerevativeFromSolutionInPoint(i, Centr, U_curr);
				auto Eps = solver_grid.GetStrainTensorFromSolutionInPoint(i, Centr, dU);
				auto Sigma = solver_grid.GetStressTensorFromSolutionInPoint(i, Centr, Eps);
				auto MisesSigma = solver_grid.GetVonMisesStress(Sigma);

				auto Eps_inv = sqrt((Eps.val[0][0] - Eps.val[1][1]) * (Eps.val[0][0] - Eps.val[1][1])
					+ (Eps.val[1][1] - Eps.val[2][2]) * (Eps.val[1][1] - Eps.val[2][2])
					+ (Eps.val[0][0] - Eps.val[2][2]) * (Eps.val[0][0] - Eps.val[2][2])
					+ 3 * (Eps.val[0][1] * Eps.val[1][0] + Eps.val[0][2] * Eps.val[2][0] + Eps.val[1][2] * Eps.val[2][1]) / 2.0) * sqrt(2.) / 3.;

				value[0][i] = MisesSigma;
				value[1][i] = Sigma.val[2][2];
			}
			value[2].resize(solver_grid.GetDOFsCount());
			value[3].resize(solver_grid.GetDOFsCount());
			value[4].resize(solver_grid.GetDOFsCount());
			value[5].resize(solver_grid.GetDOFsCount());
			for (int i = 0; i < value[3].size(); i++)
			{
				value[2][i] = U_curr[i].x;
				value[3][i] = U_curr[i].y;
				value[4][i] = U_curr[i].z;
				value[5][i] = (U_curr[i].z - U_prev[i].z) / TIME_h;;
			}

			//solver_grid.MoveCoordinates(U_curr);
			solver_grid.printTecPlot3D_DiffDomains(fout_tech, value, name_value, name_in_file, TIME_curr);
			//solver_grid.ReMoveCoordinates(U_curr);
			fclose(fout_tech);
		}
		//via XY plane
		if (is_print_result) {
			math::SimpleGrid grid_XY_plane; //input
			char in_file[1000];
			sprintf_s(in_file, "%s/Plane_XY.dat", mesh_directory);

			FILE* fin_mesh;
			fopen_s(&fin_mesh, in_file, "r");
			if (fin_mesh != NULL)
			{
				fclose(fin_mesh);
				grid_XY_plane.ReadFromSalomeDat(in_file, 2);

				FILE* fout_tech;
				char name_u_tech[5000];
				sprintf_s(name_u_tech, "%s/U_plane_XY_s%d_t%.2e.dat", result_directory_YZ, id_STEP, TIME_curr);
				fopen_s(&fout_tech, name_u_tech, "w");
				char name_in_file[1000];
				sprintf_s(name_in_file, "Time_%.4e_XY", TIME_curr);
				std::vector<std::vector<char>> name_value(6);
				char name_v_tmp[6][100];
				sprintf_s(name_v_tmp[0], "sigma_Mises");
				sprintf_s(name_v_tmp[1], "sigma_xx");
				sprintf_s(name_v_tmp[2], "Ux");
				sprintf_s(name_v_tmp[3], "Uy");
				sprintf_s(name_v_tmp[4], "Uz");
				sprintf_s(name_v_tmp[5], "Vx");
				for (int i = 0; i < name_value.size(); i++)
				{
					name_value[i].resize(100);
					for (int j = 0; j < name_value[i].size(); j++)
					{
						name_value[i][j] = name_v_tmp[i][j];
					}
				}
				std::vector<std::vector<double>> value(3 * 2);
				value[0].resize(grid_XY_plane.nvtr.size());
				value[1].resize(grid_XY_plane.nvtr.size());
				value[2].resize(grid_XY_plane.nvtr.size());
				value[3].resize(grid_XY_plane.nvtr.size());
				value[4].resize(grid_XY_plane.nvtr.size());
				value[5].resize(grid_XY_plane.nvtr.size());
				double sigma_inv_max = 0;
				int elem_sigma_max = 0;
				for (int i = 0; i < grid_XY_plane.nvtr.size(); i++)
				{
					Point<double> Centr;
					for (int j = 0; j < grid_XY_plane.nvtr[i].size(); j++)
					{
						Centr += grid_XY_plane.xyz[grid_XY_plane.nvtr[i][j]];
					}
					Centr /= grid_XY_plane.nvtr[i].size();

					double len;
					int id_elem = solver_grid.GetNearestElementID(Centr, len);
					if (id_elem >= 0)
					{
						auto element = solver_grid.GetElement(id_elem);

						Point<double> U = solver_grid.GetSolutionInPoint(id_elem, Centr, U_curr);
						Point<double> U_prev_ = solver_grid.GetSolutionInPoint(id_elem, Centr, U_prev);
						Point<Point<double>> dU = solver_grid.GetDerevativeFromSolutionInPoint(id_elem, Centr, U_curr);
						auto Eps = solver_grid.GetStrainTensorFromSolutionInPoint(id_elem, Centr, dU);
						auto Sigma = solver_grid.GetStressTensorFromSolutionInPoint(id_elem, Centr, Eps);
						auto MisesSigma = solver_grid.GetVonMisesStress(Sigma);

						auto Eps_inv = sqrt((Eps.val[0][0] - Eps.val[1][1]) * (Eps.val[0][0] - Eps.val[1][1])
							+ (Eps.val[1][1] - Eps.val[2][2]) * (Eps.val[1][1] - Eps.val[2][2])
							+ (Eps.val[0][0] - Eps.val[2][2]) * (Eps.val[0][0] - Eps.val[2][2])
							+ 3 * (Eps.val[0][1] * Eps.val[1][0] + Eps.val[0][2] * Eps.val[2][0] + Eps.val[1][2] * Eps.val[2][1]) / 2.0) * sqrt(2.) / 3.;

						value[0][i] = MisesSigma;
						value[1][i] = Sigma.val[0][0];

						value[2][i] = U.x;
						value[3][i] = U.y;
						value[4][i] = U.z;
						value[5][i] = (U_prev_.x - U.x) / TIME_h;
					}
				}

				grid_XY_plane.printTecPlot3D(fout_tech, value, name_value, name_in_file, TIME_curr);
				fclose(fout_tech);
			}
		}
		//via YZ plane
		if (is_print_result) {
			math::SimpleGrid grid_YZ_plane; //input
			char in_file[1000];
			sprintf_s(in_file, "%s/Plane_YZ.dat", mesh_directory);

			FILE* fin_mesh;
			fopen_s(&fin_mesh, in_file, "r");
			if (fin_mesh != NULL)
			{
				fclose(fin_mesh);
				grid_YZ_plane.ReadFromSalomeDat(in_file, 2);

				FILE* fout_tech;
				char name_u_tech[5000];
				sprintf_s(name_u_tech, "%s/U_plane_YZ_s%d_t%.2e.dat", result_directory_YZ, id_STEP, TIME_curr);
				fopen_s(&fout_tech, name_u_tech, "w");
				char name_in_file[1000];
				sprintf_s(name_in_file, "Time_%.4e_YZ", TIME_curr);
				std::vector<std::vector<char>> name_value(6);
				char name_v_tmp[6][100];
				sprintf_s(name_v_tmp[0], "sigma_Mises");
				sprintf_s(name_v_tmp[1], "sigma_zz");
				sprintf_s(name_v_tmp[2], "Ux");
				sprintf_s(name_v_tmp[3], "Uy");
				sprintf_s(name_v_tmp[4], "Uz");
				sprintf_s(name_v_tmp[5], "Vz");
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
						Point<double> U_prev_ = solver_grid.GetSolutionInPoint(id_elem, Centr, U_prev);
						Point<Point<double>> dU = solver_grid.GetDerevativeFromSolutionInPoint(id_elem, Centr, U_curr);
						auto Eps = solver_grid.GetStrainTensorFromSolutionInPoint(id_elem, Centr, dU);
						auto Sigma = solver_grid.GetStressTensorFromSolutionInPoint(id_elem, Centr, Eps);
						auto MisesSigma = solver_grid.GetVonMisesStress(Sigma);

						auto Eps_inv = sqrt((Eps.val[0][0] - Eps.val[1][1]) * (Eps.val[0][0] - Eps.val[1][1])
							+ (Eps.val[1][1] - Eps.val[2][2]) * (Eps.val[1][1] - Eps.val[2][2])
							+ (Eps.val[0][0] - Eps.val[2][2]) * (Eps.val[0][0] - Eps.val[2][2])
							+ 3 * (Eps.val[0][1] * Eps.val[1][0] + Eps.val[0][2] * Eps.val[2][0] + Eps.val[1][2] * Eps.val[2][1]) / 2.0) * sqrt(2.) / 3.;

						value[0][i] = MisesSigma;
						value[1][i] = Sigma.val[2][2];

						value[2][i] = U.x;
						value[3][i] = U.y;
						value[4][i] = U.z;
						value[5][i] = (U_prev_.z - U.z) / TIME_h;
					}
				}

				grid_YZ_plane.printTecPlot3D(fout_tech, value, name_value, name_in_file, TIME_curr);
				fclose(fout_tech);
			}
		}
		//via XZ plane
		if (is_print_result) {
			math::SimpleGrid grid_XZ_plane; //input
			char in_file[1000];
			sprintf_s(in_file, "%s/Plane_XZ.dat", mesh_directory);

			FILE* fin_mesh;
			fopen_s(&fin_mesh, in_file, "r");
			if (fin_mesh != NULL)
			{
				fclose(fin_mesh);
				grid_XZ_plane.ReadFromSalomeDat(in_file, 2);

				FILE* fout_tech;
				char name_u_tech[5000];
				sprintf_s(name_u_tech, "%s/U_plane_XZ_s%d_t%.2e.dat", result_directory_XZ, id_STEP, TIME_curr);
				fopen_s(&fout_tech, name_u_tech, "w");
				char name_in_file[1000];
				sprintf_s(name_in_file, "Time_%.4e_XZ", TIME_curr);
				std::vector<std::vector<char>> name_value(6);
				char name_v_tmp[6][100];
				sprintf_s(name_v_tmp[0], "sigma_Mises");
				sprintf_s(name_v_tmp[1], "sigma_zz");
				sprintf_s(name_v_tmp[2], "Ux");
				sprintf_s(name_v_tmp[3], "Uy");
				sprintf_s(name_v_tmp[4], "Uz");
				sprintf_s(name_v_tmp[5], "Vz");
				for (int i = 0; i < name_value.size(); i++)
				{
					name_value[i].resize(100);
					for (int j = 0; j < name_value[i].size(); j++)
					{
						name_value[i][j] = name_v_tmp[i][j];
					}
				}
				std::vector<std::vector<double>> value(3 * 2);
				value[0].resize(grid_XZ_plane.nvtr.size());
				value[1].resize(grid_XZ_plane.nvtr.size());
				value[2].resize(grid_XZ_plane.nvtr.size());
				value[3].resize(grid_XZ_plane.nvtr.size());
				value[4].resize(grid_XZ_plane.nvtr.size());
				value[5].resize(grid_XZ_plane.nvtr.size());
				double sigma_inv_max = 0;
				int elem_sigma_max = 0;
				for (int i = 0; i < grid_XZ_plane.nvtr.size(); i++)
				{
					Point<double> Centr;
					for (int j = 0; j < grid_XZ_plane.nvtr[i].size(); j++)
					{
						Centr += grid_XZ_plane.xyz[grid_XZ_plane.nvtr[i][j]];
					}
					Centr /= grid_XZ_plane.nvtr[i].size();

					double len;
					int id_elem = solver_grid.GetNearestElementID(Centr, len);
					if (id_elem >= 0)
					{
						auto element = solver_grid.GetElement(id_elem);

						Point<double> U = solver_grid.GetSolutionInPoint(id_elem, Centr, U_curr);
						Point<double> U_prev_ = solver_grid.GetSolutionInPoint(id_elem, Centr, U_prev);
						Point<Point<double>> dU = solver_grid.GetDerevativeFromSolutionInPoint(id_elem, Centr, U_curr);
						auto Eps = solver_grid.GetStrainTensorFromSolutionInPoint(id_elem, Centr, dU);
						auto Sigma = solver_grid.GetStressTensorFromSolutionInPoint(id_elem, Centr, Eps);
						auto MisesSigma = solver_grid.GetVonMisesStress(Sigma);

						auto Eps_inv = sqrt((Eps.val[0][0] - Eps.val[1][1]) * (Eps.val[0][0] - Eps.val[1][1])
							+ (Eps.val[1][1] - Eps.val[2][2]) * (Eps.val[1][1] - Eps.val[2][2])
							+ (Eps.val[0][0] - Eps.val[2][2]) * (Eps.val[0][0] - Eps.val[2][2])
							+ 3 * (Eps.val[0][1] * Eps.val[1][0] + Eps.val[0][2] * Eps.val[2][0] + Eps.val[1][2] * Eps.val[2][1]) / 2.0) * sqrt(2.) / 3.;

						value[0][i] = MisesSigma;
						value[1][i] = Sigma.val[2][2];

						value[2][i] = U.x;
						value[3][i] = U.y;
						value[4][i] = U.z;
						value[5][i] = (U_prev_.z - U.z) / TIME_h;
					}
				}

				grid_XZ_plane.printTecPlot3D(fout_tech, value, name_value, name_in_file, TIME_curr);
				fclose(fout_tech);
			}
		}
		//U in point
		if (is_print_result /*&& TIME_curr > TIME_L*/)
		{
			for (int r = 0; r < receivers.size(); r++)
			{
				FILE* fout_Uz;
				fopen_s(&fout_Uz, receivers[r].file_name, "a");
				for (int i = 0; i < receivers[r].U_in_point.size(); i++)
				{
					fprintf_s(fout_Uz, "%.10e %.10e %.10e %.10e\n", receivers_times[i], receivers[r].U_in_point[i].x, receivers[r].U_in_point[i].y, receivers[r].U_in_point[i].z);
				}
				receivers[r].U_in_point.clear();
				fclose(fout_Uz);
			}
			receivers_times.clear();
		}

		printf_s("\t complite\n");

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
	}
}
void ElastodynamicsProblem_ExplicitRungeKutta4(char properties_file[1000])
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
	//char properties_file[1000] = { "D:/Elastodynamic/homocyl/67k/param_for_solver.txt" };
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


	{
		char _line[1000];
		std::vector<int> val;

		math::ReadNonEmptyLine_forNumbers(f_properties, _line);
		math::ParserStringToVectorInt(_line, val, " ");
		if (val[0] == 1)
			is_print_logFile = true;


		math::ReadNonEmptyLine_forNumbers(f_properties, _line);
		math::ParserStringToVectorInt(_line, val, " ");
		math::NUM_THREADS = val[0];
	}

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

	wchar_t _tmp_wc[1000];
	math::Char_To_Wchar_t(base_result_directory, _tmp_wc, 1000);
	CreateDirectory((LPCTSTR)_tmp_wc, NULL);

	bool res = false;
	double TIME_h = step_size;
	if (is_STATIONARY) {
		end_iteration = 1;
		start_iteration = 0;
	}
	std::vector<Point<double>> U_prev(geo_grid.xyz.size()), U_curr(geo_grid.xyz.size());
	math::InitializationVector(U_curr, 0);
	math::InitializationVector(U_prev, 0);
	std::vector<Point<double>> dU_dt_prev(geo_grid.xyz.size()), dU_dt_curr(geo_grid.xyz.size());
	math::InitializationVector(dU_dt_curr, 0);
	math::InitializationVector(dU_dt_prev, 0);

	CSSD_Matrix<Tensor2Rank3D, Point<double>> StiffnessMatrix;
	CSSD_Matrix<Tensor2Rank3D, Point<double>> MassMatrix;
	CSSD_Matrix<double, double> MassMatrix_doubleSLAE;
	CSSD_Matrix<double, double> Precond;
	CSSD_Matrix<Tensor2Rank3D, Point<double>> VolumeForceVector;
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
	MassMatrix.SetMatrix(StiffnessMatrix);
	VolumeForceVector.F.resize(StiffnessMatrix.F.size());
	printf_s("\t\tcomplite\n");

	//строим базовые матрицы жесткости/массы/вектор силы тяжести
	//матрица массы идет с коэффициентом rpho
	printf_s("================= Create matrix ================\n");
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
			return solver_grid.GetDomain(solver_grid.GetElement(elem)->GetIdDomain())->forThermal.rpho;
		};
		std::function<Point<double>(int, Point<double>)> VolumeForсe = [&](int elem, Point<double> X)->Point<double>
		{
			return Point<double>(0, 0, 0);
		};

		printf("Matrix assembling...\n");
		std::vector<DenseMatrix<Tensor2Rank3D, Point<double>>> local_SLAE_stiffness(solver_grid.GetElementsCount());
		std::vector<DenseMatrix<Tensor2Rank3D, Point<double>>> local_SLAE_mass(solver_grid.GetElementsCount());
		std::vector< std::vector<Point<double>>> local_force(solver_grid.GetElementsCount());
		omp_set_num_threads(math::NUM_THREADS);
#pragma omp parallel for schedule(dynamic)
		for (int id_elem = 0; id_elem < solver_grid.GetElementsCount(); id_elem++)
		{
			if (id_elem % 1000 == 0)
				printf("Solve element[%d]\r", id_elem);
			auto element = solver_grid.GetElement(id_elem);

			std::function<std::vector<std::vector<double>>(Point<double>)> D = [&](Point<double> X) {return StiffnessCoef(id_elem, X); };
			std::function<double(Point<double>)> M = [&](Point<double> X) {return MassCoef(id_elem, X); };
			std::function<Point<double>(Point<double>)> F = [&](Point<double> X) {return VolumeForсe(id_elem, X); };

			element->SolveLocalMatrix(local_SLAE_stiffness[id_elem], D);
			element->SolveMassMatrix(local_SLAE_mass[id_elem], M);
			element->SolveRightSide(local_force[id_elem], F);
		}
		for (int id_elem = 0; id_elem < solver_grid.GetElementsCount(); id_elem++)
		{
			if (id_elem % 1000 == 0)
				printf("Add the local matrix of element[%d]\r", id_elem);
			StiffnessMatrix.SummPartOfMatrix(local_SLAE_stiffness[id_elem], *solver_grid.GetElementDOFs(id_elem));
			MassMatrix.SummPartOfMatrix(local_SLAE_mass[id_elem], *solver_grid.GetElementDOFs(id_elem));
			VolumeForceVector.SummPartOfVector(local_force[id_elem], *solver_grid.GetElementDOFs(id_elem));
		}
		printf_s("                                                                                    \r");
		printf_s("\t\tcomplite\n");

		//учитываем первые краевые ТОЛЬКО в матрице массы БЕЗ симметризации
		if(false) {
			//обнуляем строки
			for (int id_vertex = 0; id_vertex < solver_grid.boundary_vertexes.size(); id_vertex++)
			{
				auto boundary = &(solver_grid.boundary_vertexes[id_vertex]);

				int id_type = boundary->id_type;
				Point<bool> is_take;
				Point<double> boundary_value = boundary->boundary_value(is_take, solver_grid.boundary_vertexes[id_vertex].GetIdNode(0));
				int global_id = boundary->GetDOFInLocalID(0);

				auto enter_boundary = [global_id, &boundary_value](int position, CSSD_Matrix<Tensor2Rank3D, Point<double>>& global_SLAE) {
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
				};
				if (is_take.x == true) enter_boundary(0, MassMatrix);
				if (is_take.y == true) enter_boundary(1, MassMatrix);
				if (is_take.z == true) enter_boundary(2, MassMatrix);
			}
		}
		if (true) {
			//большое число на диагональ
			for (int id_vertex = 0; id_vertex < solver_grid.boundary_vertexes.size(); id_vertex++)
			{
				auto boundary = &(solver_grid.boundary_vertexes[id_vertex]);

				int id_type = boundary->id_type;
				Point<bool> is_take;
				Point<double> boundary_value = boundary->boundary_value(is_take, solver_grid.boundary_vertexes[id_vertex].GetIdNode(0));
				int global_id = boundary->GetDOFInLocalID(0);

				auto enter_boundary = [global_id, &boundary_value](int position, CSSD_Matrix<Tensor2Rank3D, Point<double>>& global_SLAE) {
					global_SLAE.Diag[global_id].val[position][position] = 1e+30;
				};
				if (is_take.x == true) enter_boundary(0, MassMatrix);
				if (is_take.y == true) enter_boundary(1, MassMatrix);
				if (is_take.z == true) enter_boundary(2, MassMatrix);
			}
		}

		//переводим матрицу массы в действительную
		math::MakeCopyMatrix_A_into_B(MassMatrix, MassMatrix_doubleSLAE);
		Precond.PrecondorSSOR(0.75, MassMatrix_doubleSLAE);
	}

	///-------------------
	Point<double> point_for_out(0, 0, 0.5 - 0.01);
	std::vector<double> Uz_in_point;
	std::vector<double> Time_for_Uz;
	double TIME_L;
	int id_elem_for_out;
	char name_uz[5000];
	{
		double v = solver_grid.GetDomain(0)->forMech.GetLongitudinalWaveVelocity_Vp(solver_grid.GetDomain(0)->forThermal.rpho);
		double L = 0.5;
		TIME_L = L / v;
		double len;
		id_elem_for_out = solver_grid.GetNearestElementID(point_for_out, len);

		wchar_t _tmp_wc[1000];
		math::Char_To_Wchar_t(base_result_directory, _tmp_wc, 1000);
		CreateDirectory((LPCTSTR)_tmp_wc, NULL);

		FILE* fout_Uz;
		sprintf_s(name_uz, "%s/Uz_in_point.txt", base_result_directory);
		fopen_s(&fout_Uz, name_uz, "w");
		fprintf_s(fout_Uz, "Point: (%.2e, %.2e, %.2e)\n", point_for_out.x, point_for_out.y, point_for_out.z);
		fprintf_s(fout_Uz, "Id element: %d\n", id_elem_for_out);
		fprintf_s(fout_Uz, "T0: %.5e sec\n", TIME_L);
		fclose(fout_Uz);

		if (id_elem_for_out < 0) return;
	}
	///-------------------


	for (int id_STEP = start_iteration; id_STEP <= end_iteration; id_STEP++)
	{
		printf_s("\n================= Start solution of %d STEP (time = %.2e) ================\n", id_STEP, TIME_h * id_STEP);
		double TIME_curr = TIME_h * id_STEP + TIME_h;
		bool is_print_result = false;
		char result_directory[1000];
		if (id_STEP % step_for_out == 0)
		{
			is_print_result = true;
			sprintf_s(result_directory, sizeof(result_directory), "%s", base_result_directory);
			//sprintf_s(result_directory, sizeof(result_directory), "%s/STEP_%d_t=%.2e", base_result_directory, id_STEP, TIME_curr);
			wchar_t _tmp_wc[1000];
			math::Char_To_Wchar_t(result_directory, _tmp_wc, 1000);
			CreateDirectory((LPCTSTR)_tmp_wc, NULL);
		}

		//переопределяем источник и строим вектор в правую часть (через вторые краевые)
		CSSD_Matrix<Tensor2Rank3D, Point<double>> SourseVector;
		SourseVector.F.resize(StiffnessMatrix.F.size());
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
			std::function<Point<double>(Point<double>)> new_sourse_value_sin = [TIME_curr, power, &solver_grid, TIME_h](Point<double> X) -> Point<double>
			{
				//double sourse = exp(-pow(TIME_curr - tay_0, degree) / pow(tay_plus, degree));
				//return power;
				//return Point<double>(power.x * sourse, power.y * sourse, power.z * sourse);
				double v = solver_grid.GetDomain(0)->forMech.GetLongitudinalWaveVelocity_Vp(solver_grid.GetDomain(0)->forThermal.rpho);
				double L = 0.5;
				double t0 = v / L;
				double w = 2 * M_PI / t0 * 10;
				double alpha = 5e+4;

				double sourse = 10e+6 * exp(-1 * alpha * TIME_curr) * sin(w * TIME_curr);
				if (TIME_curr > t0) sourse = 0;
				return Point<double>(0 * sourse, 0 * sourse, -1 * sourse /** TIME_h / 14.0e-7*/);
			};
			if (second_boundary.size() > 0)
			{
				second_boundary[0].value = new_sourse_value_sin;
				//second_boundary[1].value = new_sourse_value_sin;
				//second_boundary[2].value = new_sourse_value;

				for (int i = 0; i < solver_grid.boundary_faces.size(); i++)
				{
					solver_grid.boundary_faces[i].boundary_value = second_boundary[solver_grid.boundary_faces[i].id_type].value;
				}
			}

			for (int id_triangle = 0; id_triangle < solver_grid.boundary_faces.size(); id_triangle++)
			{
				std::vector<Point<double>> local_vector_SLAE;
				solver_grid.boundary_faces[id_triangle].SolveLocalBoundaryVector(local_vector_SLAE, solver_grid.boundary_faces[id_triangle].boundary_value);
				SourseVector.SummPartOfVector(local_vector_SLAE, *solver_grid.boundary_faces[id_triangle].GetElementDOFs());
			};
		}

		//строим коэффициенты Рунге-Кутта через решение СЛАУ
		{
			std::vector<double> K[4];
			std::vector<double> L[4];
			std::vector<Point<double>> Ku0_prev;
			std::vector<double> F_slae(U_prev.size() * 3);

			auto SolveSLAE = [&MassMatrix_doubleSLAE, &solver_grid, &Precond](std::vector<double>& F, std::vector<double>& result, double critical_residual) -> void
			{
				//краевые в правую часть + симметризация 
				if(false){
					//добавляем значения в правую часть
					for (int id_vertex = 0; id_vertex < solver_grid.boundary_vertexes.size(); id_vertex++)
					{
						auto boundary = &(solver_grid.boundary_vertexes[id_vertex]);

						int id_type = boundary->id_type;
						Point<bool> is_take;
						Point<double> boundary_value = boundary->boundary_value(is_take, solver_grid.boundary_vertexes[id_vertex].GetIdNode(0));
						int global_id = boundary->GetDOFInLocalID(0);

						if (global_id < MassMatrix_doubleSLAE.GetMatrixSize() / 3)
						{
							auto enter_value = [&result, &F](int value_id, double value) ->void
							{
								result[value_id] = value;
								F[value_id] = value;
								return;
							};

							if (is_take.x) enter_value(global_id * 3 + 0, boundary_value.x);
							if (is_take.y) enter_value(global_id * 3 + 1, boundary_value.y);
							if (is_take.z) enter_value(global_id * 3 + 2, boundary_value.z);
						}

					}
					//симметризация
					omp_set_num_threads(math::NUM_THREADS);
#pragma omp parallel for schedule(dynamic) 
					for (int id_row = 0; id_row < MassMatrix_doubleSLAE.GetMatrixSize(); id_row++)
					{
						if (id_row % 100 == 0)
						{
							printf("\tcurrent %d/%d\r", id_row, MassMatrix_doubleSLAE.GetMatrixSize());
						}
						int iterator_in_boundary = 0;
						for (int jj = 0; jj < MassMatrix_doubleSLAE.id_column_for_A_up[id_row].size(); jj++)
						{
							int id_column = MassMatrix_doubleSLAE.id_column_for_A_up[id_row][jj];
							for (/*iterator_in_boundary = 0*/; iterator_in_boundary < solver_grid.boundary_vertexes.size(); iterator_in_boundary++)
							{
								if (solver_grid.boundary_vertexes[iterator_in_boundary].GetDOFInLocalID(0) * 3 + 0 == id_column)
								{
									Point<bool> is_take;
									Point<double> boundary_value = solver_grid.boundary_vertexes[iterator_in_boundary].boundary_value(is_take, solver_grid.boundary_vertexes[iterator_in_boundary].GetIdNode(0));

									if (is_take.x)
									{
										F[id_row] -= MassMatrix_doubleSLAE.A_up[id_row][jj] * F[id_column];
										MassMatrix_doubleSLAE.A_up[id_row][jj] = 0;
									}

									break;
								}
								if (solver_grid.boundary_vertexes[iterator_in_boundary].GetDOFInLocalID(0) * 3 + 1 == id_column)
								{
									Point<bool> is_take;
									Point<double> boundary_value = solver_grid.boundary_vertexes[iterator_in_boundary].boundary_value(is_take, solver_grid.boundary_vertexes[iterator_in_boundary].GetIdNode(0));

									if (is_take.y)
									{
										F[id_row] -= MassMatrix_doubleSLAE.A_up[id_row][jj] * F[id_column];
										MassMatrix_doubleSLAE.A_up[id_row][jj] = 0;
									}

									break;
								}
								if (solver_grid.boundary_vertexes[iterator_in_boundary].GetDOFInLocalID(0) * 3 + 2 == id_column)
								{
									Point<bool> is_take;
									Point<double> boundary_value = solver_grid.boundary_vertexes[iterator_in_boundary].boundary_value(is_take, solver_grid.boundary_vertexes[iterator_in_boundary].GetIdNode(0));

									if (is_take.z)
									{
										F[id_row] -= MassMatrix_doubleSLAE.A_up[id_row][jj] * F[id_column];
										MassMatrix_doubleSLAE.A_up[id_row][jj] = 0;
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
						for (int jj = 0; jj < MassMatrix_doubleSLAE.id_column_for_A_down[id_row].size(); jj++)
						{
							int id_column = MassMatrix_doubleSLAE.id_column_for_A_down[id_row][jj];
							for (/*iterator_in_boundary = 0*/; iterator_in_boundary < solver_grid.boundary_vertexes.size(); iterator_in_boundary++)
							{
								if (solver_grid.boundary_vertexes[iterator_in_boundary].GetDOFInLocalID(0) * 3 + 0 == id_column)
								{
									Point<bool> is_take;
									Point<double> boundary_value = solver_grid.boundary_vertexes[iterator_in_boundary].boundary_value(is_take, solver_grid.boundary_vertexes[iterator_in_boundary].GetIdNode(0));

									if (is_take.x)
									{
										F[id_row] -= MassMatrix_doubleSLAE.A_down[id_row][jj] * F[id_column];
										MassMatrix_doubleSLAE.A_down[id_row][jj] = 0;
									}

									break;
								}
								if (solver_grid.boundary_vertexes[iterator_in_boundary].GetDOFInLocalID(0) * 3 + 1 == id_column)
								{
									Point<bool> is_take;
									Point<double> boundary_value = solver_grid.boundary_vertexes[iterator_in_boundary].boundary_value(is_take, solver_grid.boundary_vertexes[iterator_in_boundary].GetIdNode(0));

									if (is_take.y)
									{
										F[id_row] -= MassMatrix_doubleSLAE.A_down[id_row][jj] * F[id_column];
										MassMatrix_doubleSLAE.A_down[id_row][jj] = 0;
									}

									break;
								}
								if (solver_grid.boundary_vertexes[iterator_in_boundary].GetDOFInLocalID(0) * 3 + 2 == id_column)
								{
									Point<bool> is_take;
									Point<double> boundary_value = solver_grid.boundary_vertexes[iterator_in_boundary].boundary_value(is_take, solver_grid.boundary_vertexes[iterator_in_boundary].GetIdNode(0));

									if (is_take.z)
									{
										F[id_row] -= MassMatrix_doubleSLAE.A_down[id_row][jj] * F[id_column];
										MassMatrix_doubleSLAE.A_down[id_row][jj] = 0;
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
				//краевые в правую часть большим числом 
				if (true) {
					//добавляем значения в правую часть
					for (int id_vertex = 0; id_vertex < solver_grid.boundary_vertexes.size(); id_vertex++)
					{
						auto boundary = &(solver_grid.boundary_vertexes[id_vertex]);

						int id_type = boundary->id_type;
						Point<bool> is_take;
						Point<double> boundary_value = boundary->boundary_value(is_take, solver_grid.boundary_vertexes[id_vertex].GetIdNode(0));
						int global_id = boundary->GetDOFInLocalID(0);

						if (global_id < MassMatrix_doubleSLAE.GetMatrixSize() / 3)
						{
							auto enter_value = [&result, &F, &MassMatrix_doubleSLAE](int value_id, double value) ->void
							{
								result[value_id] = value;
								F[value_id] = value * MassMatrix_doubleSLAE.Diag[value_id];
								return;
							};

							if (is_take.x) enter_value(global_id * 3 + 0, boundary_value.x);
							if (is_take.y) enter_value(global_id * 3 + 1, boundary_value.y);
							if (is_take.z) enter_value(global_id * 3 + 2, boundary_value.z);
						}

					}
				}

				printf("Soluting SLAY... (%d)\n", MassMatrix_doubleSLAE.GetMatrixSize());
				int MaxSize = MassMatrix_doubleSLAE.GetMatrixSize();
				MaxSize = MaxSize / 10 < 1000 ? 1000 : MaxSize / 10;
				std::vector<double> best_solution;
				math::MakeCopyVector_A_into_B(result, best_solution);
				double current_residual = 1, best_residual = 1e+25;
				int MAX_STEPS = 4;
				int ii = 1;
				
				for (int i = 0; i <= MAX_STEPS; i++)
				{
					printf_s("//---> I = %d/%d (full size %d) - needed residual %.2e\n", i, MAX_STEPS, MassMatrix_doubleSLAE.GetMatrixSize(), critical_residual);

					MassMatrix_doubleSLAE.print_logs = true;
					current_residual = abs(MassMatrix_doubleSLAE.MSG_PreconditioningSSOR(F, result, MaxSize, critical_residual, Precond));

					if (current_residual <= critical_residual)
					{
						best_residual = current_residual;
						math::MakeCopyVector_A_into_B(result, best_solution);
						break;
					}
					if (current_residual < best_residual)
					{
						best_residual = current_residual;
						math::MakeCopyVector_A_into_B(result, best_solution);
					}
				}
				math::MakeCopyVector_A_into_B(best_solution, result);

				if (best_residual >= 1)
				{
					printf_s("There are problems in the solution\n");
					int a;
					scanf_s("%d", &a);
				}
			};

			//M*u''=F-Ku0 => u''=K[0]
			{
				StiffnessMatrix.MultiplicationMatrixVector(U_prev, Ku0_prev);
				for (int i = 0; i < Ku0_prev.size(); i++)
				{
					Point<double> res = Point<double>(0, 0, 0)
						- Ku0_prev[i]
						+ VolumeForceVector.F[i]
						+ SourseVector.F[i];
					F_slae[i * 3 + 0] = res.x;
					F_slae[i * 3 + 1] = res.y;
					F_slae[i * 3 + 2] = res.z;
				}
				K[0].resize(F_slae.size());
				math::InitializationVector(K[0], 1e-10);
				SolveSLAE(F_slae, K[0], 1e-10);

				L[0].resize(F_slae.size());
				for (int i = 0; i < dU_dt_prev.size(); i++)
				{
					K[0][i * 3 + 0] *= TIME_h;
					K[0][i * 3 + 1] *= TIME_h;
					K[0][i * 3 + 2] *= TIME_h;

					L[0][i * 3 + 0] = dU_dt_prev[i].x * TIME_h;
					L[0][i * 3 + 1] = dU_dt_prev[i].y * TIME_h;
					L[0][i * 3 + 2] = dU_dt_prev[i].z * TIME_h;
				}
			}

			//M*u''=F-K(u0+L0*h/2) => u''=K1
			{
				std::vector<Point<double>> U_tmp(U_prev.size());
				for (int i = 0; i < U_prev.size(); i++)
					U_tmp[i] = U_prev[i] + Point<double>(L[0][i * 3 + 0], L[0][i * 3 + 1], L[0][i * 3 + 2]) / 2.0;
				StiffnessMatrix.MultiplicationMatrixVector(U_tmp, Ku0_prev);
				for (int i = 0; i < Ku0_prev.size(); i++)
				{
					Point<double> res = Point<double>(0, 0, 0)
						- Ku0_prev[i]
						+ VolumeForceVector.F[i]
						+ SourseVector.F[i];
					F_slae[i * 3 + 0] = res.x;
					F_slae[i * 3 + 1] = res.y;
					F_slae[i * 3 + 2] = res.z;
				}
				math::MakeCopyVector_A_into_B(K[0], K[1]);
				SolveSLAE(F_slae, K[1], 1e-10);

				L[1].resize(F_slae.size());
				for (int i = 0; i < dU_dt_prev.size(); i++)
				{
					K[1][i * 3 + 0] *= TIME_h;
					K[1][i * 3 + 1] *= TIME_h;
					K[1][i * 3 + 2] *= TIME_h;

					L[1][i * 3 + 0] = (dU_dt_prev[i].x + K[0][i * 3 + 0] / 2.0) * TIME_h;
					L[1][i * 3 + 1] = (dU_dt_prev[i].y + K[0][i * 3 + 1] / 2.0) * TIME_h;
					L[1][i * 3 + 2] = (dU_dt_prev[i].z + K[0][i * 3 + 3] / 2.0) * TIME_h;
				}
			}

			//M*u''=F-K(u0+K1*h/2) => u''=K2
			{
				std::vector<Point<double>> U_tmp(U_prev.size());
				for (int i = 0; i < U_prev.size(); i++)
					U_tmp[i] = U_prev[i] + Point<double>(L[1][i * 3 + 0], L[1][i * 3 + 1], L[1][i * 3 + 2]) / 2.0;
				StiffnessMatrix.MultiplicationMatrixVector(U_tmp, Ku0_prev);
				for (int i = 0; i < Ku0_prev.size(); i++)
				{
					Point<double> res = Point<double>(0, 0, 0)
						- Ku0_prev[i]
						+ VolumeForceVector.F[i]
						+ SourseVector.F[i];
					F_slae[i * 3 + 0] = res.x;
					F_slae[i * 3 + 1] = res.y;
					F_slae[i * 3 + 2] = res.z;
				}
				math::MakeCopyVector_A_into_B(K[1], K[2]);
				SolveSLAE(F_slae, K[2], 1e-10);

				L[2].resize(F_slae.size());
				for (int i = 0; i < dU_dt_prev.size(); i++)
				{
					K[2][i * 3 + 0] *= TIME_h;
					K[2][i * 3 + 1] *= TIME_h;
					K[2][i * 3 + 2] *= TIME_h;

					L[2][i * 3 + 0] = (dU_dt_prev[i].x + K[1][i * 3 + 0] / 2.0) * TIME_h;
					L[2][i * 3 + 1] = (dU_dt_prev[i].y + K[1][i * 3 + 1] / 2.0) * TIME_h;
					L[2][i * 3 + 2] = (dU_dt_prev[i].z + K[1][i * 3 + 3] / 2.0) * TIME_h;
				}
			}

			//M*u''=F-K(u0+K2*h) => u''=K3
			{
				std::vector<Point<double>> U_tmp(U_prev.size());
				for (int i = 0; i < U_prev.size(); i++)
					U_tmp[i] = U_prev[i] + Point<double>(L[2][i * 3 + 0], L[2][i * 3 + 1], L[2][i * 3 + 2]);
				StiffnessMatrix.MultiplicationMatrixVector(U_tmp, Ku0_prev);
				for (int i = 0; i < Ku0_prev.size(); i++)
				{
					Point<double> res = Point<double>(0, 0, 0)
						- Ku0_prev[i]
						+ VolumeForceVector.F[i]
						+ SourseVector.F[i];
					F_slae[i * 3 + 0] = res.x;
					F_slae[i * 3 + 1] = res.y;
					F_slae[i * 3 + 2] = res.z;
				}
				math::MakeCopyVector_A_into_B(K[2], K[3]);
				SolveSLAE(F_slae, K[3], 1e-10);

				L[3].resize(F_slae.size());
				for (int i = 0; i < dU_dt_prev.size(); i++)
				{
					K[3][i * 3 + 0] *= TIME_h;
					K[3][i * 3 + 1] *= TIME_h;
					K[3][i * 3 + 2] *= TIME_h;

					L[3][i * 3 + 0] = (dU_dt_prev[i].x + K[2][i * 3 + 0]) * TIME_h;
					L[3][i * 3 + 1] = (dU_dt_prev[i].y + K[2][i * 3 + 1]) * TIME_h;
					L[3][i * 3 + 2] = (dU_dt_prev[i].z + K[2][i * 3 + 3]) * TIME_h;
				}
			}

			for (int i = 0; i < U_curr.size(); i++)
			{
				dU_dt_curr[i].x = dU_dt_curr[i].x + (K[0][i * 3 + 0] + 2 * K[1][i * 3 + 0] + 2 * K[2][i * 3 + 0] + K[3][i * 3 + 0]) / 6.0;
				dU_dt_curr[i].y = dU_dt_curr[i].y + (K[0][i * 3 + 1] + 2 * K[1][i * 3 + 1] + 2 * K[2][i * 3 + 1] + K[3][i * 3 + 1]) / 6.0;
				dU_dt_curr[i].z = dU_dt_curr[i].z + (K[0][i * 3 + 2] + 2 * K[1][i * 3 + 2] + 2 * K[2][i * 3 + 2] + K[3][i * 3 + 2]) / 6.0;
				
				U_curr[i].x = U_curr[i].x + (L[0][i * 3 + 0] + 2 * L[1][i * 3 + 0] + 2 * L[2][i * 3 + 0] + L[3][i * 3 + 0]) / 6.0;
				U_curr[i].y = U_curr[i].y + (L[0][i * 3 + 1] + 2 * L[1][i * 3 + 1] + 2 * L[2][i * 3 + 1] + L[3][i * 3 + 1]) / 6.0;
				U_curr[i].z = U_curr[i].z + (L[0][i * 3 + 2] + 2 * L[1][i * 3 + 2] + 2 * L[2][i * 3 + 2] + L[3][i * 3 + 2]) / 6.0;
			}
		}

		//вывод сейсмограммы
		if (TIME_curr > TIME_L || true) {
			Point<double> U_curr_in_point = solver_grid.GetSolutionInPoint(id_elem_for_out, point_for_out, U_curr);
			Point<double> U_prev_in_point = solver_grid.GetSolutionInPoint(id_elem_for_out, point_for_out, U_prev);
			Uz_in_point.push_back((U_curr_in_point.z + U_prev_in_point.z) / 2.0);
			Time_for_Uz.push_back(TIME_curr - TIME_h / 2.0);
		}

		//output solution
		//without deformations
		if (is_print_result) {
			printf_s("Print the mech result into .dat file... ");

			FILE* fout_tech;
			char name_u_tech[5000];
			sprintf_s(name_u_tech, "%s/U_deformation_%d_t%.2e.dat", result_directory, id_STEP, TIME_curr);
			fopen_s(&fout_tech, name_u_tech, "w");
			char name_in_file[1000];
			sprintf_s(name_in_file, "step%d_time%.2e", id_STEP, TIME_curr);
			std::vector<std::vector<char>> name_value(6);
			char name_v_tmp[6][100];
			sprintf_s(name_v_tmp[0], "sigma_Mises");
			sprintf_s(name_v_tmp[1], "sigma_zz");
			sprintf_s(name_v_tmp[2], "Ux");
			sprintf_s(name_v_tmp[3], "Uy");
			sprintf_s(name_v_tmp[4], "Uz");
			sprintf_s(name_v_tmp[5], "Vz");
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
			double sigma_inv_max = 0;
			int elem_sigma_max = 0;
			for (int i = 0; i < solver_grid.GetElementsCount(); i++)
			{
				auto element = solver_grid.GetElement(i);

				Point<double> Centr = element->GetWeightCentr();

				//Point<double> U = solver_grid.GetSolutionInPoint(i, Centr, U_curr);
				Point<Point<double>> dU = solver_grid.GetDerevativeFromSolutionInPoint(i, Centr, U_curr);
				auto Eps = solver_grid.GetStrainTensorFromSolutionInPoint(i, Centr, dU);
				auto Sigma = solver_grid.GetStressTensorFromSolutionInPoint(i, Centr, Eps);
				auto MisesSigma = solver_grid.GetVonMisesStress(Sigma);

				auto Eps_inv = sqrt((Eps.val[0][0] - Eps.val[1][1]) * (Eps.val[0][0] - Eps.val[1][1])
					+ (Eps.val[1][1] - Eps.val[2][2]) * (Eps.val[1][1] - Eps.val[2][2])
					+ (Eps.val[0][0] - Eps.val[2][2]) * (Eps.val[0][0] - Eps.val[2][2])
					+ 3 * (Eps.val[0][1] * Eps.val[1][0] + Eps.val[0][2] * Eps.val[2][0] + Eps.val[1][2] * Eps.val[2][1]) / 2.0) * sqrt(2.) / 3.;

				value[0][i] = MisesSigma;
				value[1][i] = Sigma.val[2][2];
			}
			value[2].resize(solver_grid.GetDOFsCount());
			value[3].resize(solver_grid.GetDOFsCount());
			value[4].resize(solver_grid.GetDOFsCount());
			value[5].resize(solver_grid.GetDOFsCount());
			for (int i = 0; i < value[3].size(); i++)
			{
				value[2][i] = U_curr[i].x;
				value[3][i] = U_curr[i].y;
				value[4][i] = U_curr[i].z;
				value[5][i] = dU_dt_curr[i].z;
			}

			//solver_grid.MoveCoordinates(U_curr);
			solver_grid.printTecPlot3D_DiffDomains(fout_tech, value, name_value, name_in_file, TIME_curr);
			//solver_grid.ReMoveCoordinates(U_curr);
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
			sprintf_s(name_u_tech, "%s/U_plane_YZ_s%d_t%.2e.dat", result_directory, id_STEP, TIME_curr);
			fopen_s(&fout_tech, name_u_tech, "w");
			char name_in_file[1000];
			sprintf_s(name_in_file, "Time_%.4e_YZ", TIME_curr);
			std::vector<std::vector<char>> name_value(6);
			char name_v_tmp[6][100];
			sprintf_s(name_v_tmp[0], "sigma_Mises");
			sprintf_s(name_v_tmp[1], "sigma_zz");
			sprintf_s(name_v_tmp[2], "Ux");
			sprintf_s(name_v_tmp[3], "Uy");
			sprintf_s(name_v_tmp[4], "Uz");
			sprintf_s(name_v_tmp[5], "Vz");
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
					Point<double> dU_dt = solver_grid.GetSolutionInPoint(id_elem, Centr, dU_dt_curr);
					Point<Point<double>> dU = solver_grid.GetDerevativeFromSolutionInPoint(id_elem, Centr, U_curr);
					auto Eps = solver_grid.GetStrainTensorFromSolutionInPoint(id_elem, Centr, dU);
					auto Sigma = solver_grid.GetStressTensorFromSolutionInPoint(id_elem, Centr, Eps);
					auto MisesSigma = solver_grid.GetVonMisesStress(Sigma);

					auto Eps_inv = sqrt((Eps.val[0][0] - Eps.val[1][1]) * (Eps.val[0][0] - Eps.val[1][1])
						+ (Eps.val[1][1] - Eps.val[2][2]) * (Eps.val[1][1] - Eps.val[2][2])
						+ (Eps.val[0][0] - Eps.val[2][2]) * (Eps.val[0][0] - Eps.val[2][2])
						+ 3 * (Eps.val[0][1] * Eps.val[1][0] + Eps.val[0][2] * Eps.val[2][0] + Eps.val[1][2] * Eps.val[2][1]) / 2.0) * sqrt(2.) / 3.;

					value[0][i] = MisesSigma;
					value[1][i] = Sigma.val[2][2];

					value[2][i] = U.x;
					value[3][i] = U.y;
					value[4][i] = U.z;
					value[5][i] = dU_dt.z;
				}
			}

			grid_YZ_plane.printTecPlot3D(fout_tech, value, name_value, name_in_file, TIME_curr);
			fclose(fout_tech);
		}
		//Uz in point
		if (is_print_result && TIME_curr > TIME_L)
		{
			FILE* fout_Uz;
			fopen_s(&fout_Uz, name_uz, "a");
			for (int i = 0; i < Uz_in_point.size(); i++)
			{
				fprintf_s(fout_Uz, "%.10e %.10e\n", Time_for_Uz[i], Uz_in_point[i]);
			}
			Time_for_Uz.clear();
			Uz_in_point.clear();
			fclose(fout_Uz);
		}

		printf_s("\t complite\n");

		math::MakeCopyVector_A_into_B(U_curr, U_prev);
		math::MakeCopyVector_A_into_B(dU_dt_curr, dU_dt_prev);

		//обновляем сетку
		if (false) {
			solver_grid.MoveCoordinates(U_curr);
			for (int i = 0; i < solver_grid.GetElementsCount(); i++)
			{
				solver_grid.GetElement(i)->SolveAlphaMatrix();
			}
		}
	}
}
void ElastodynamicsProblem_ExplicitPredCorrSimple(char properties_file[1000])
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
	//char properties_file[1000] = { "D:/Elastodynamic/homocyl/67k/param_for_solver.txt" };
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


	{
		char _line[1000];
		std::vector<int> val;

		math::ReadNonEmptyLine_forNumbers(f_properties, _line);
		math::ParserStringToVectorInt(_line, val, " ");
		if (val[0] == 1)
			is_print_logFile = true;


		math::ReadNonEmptyLine_forNumbers(f_properties, _line);
		math::ParserStringToVectorInt(_line, val, " ");
		math::NUM_THREADS = val[0];
	}

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

	wchar_t _tmp_wc[1000];
	math::Char_To_Wchar_t(base_result_directory, _tmp_wc, 1000);
	CreateDirectory((LPCTSTR)_tmp_wc, NULL);

	bool res = false;
	double TIME_h = step_size;
	if (is_STATIONARY) {
		end_iteration = 1;
		start_iteration = 0;
	}
	std::vector<Point<double>> U_prev(geo_grid.xyz.size()), U_prevprev(geo_grid.xyz.size()), U_curr(geo_grid.xyz.size()), U_pred(geo_grid.xyz.size());
	math::InitializationVector(U_curr, 0);
	math::InitializationVector(U_prev, 0);
	math::InitializationVector(U_prevprev, 0);
	math::InitializationVector(U_pred, 0);
	std::vector<Point<double>> dU_dt_prev(geo_grid.xyz.size()), dU_dt_prevprev(geo_grid.xyz.size()), dU_dt_curr(geo_grid.xyz.size()), dU_dt_pred(geo_grid.xyz.size());
	math::InitializationVector(dU_dt_curr, 0);
	math::InitializationVector(dU_dt_prev, 0);
	math::InitializationVector(dU_dt_prevprev, 0);
	math::InitializationVector(dU_dt_pred, 0);

	CSSD_Matrix<Tensor2Rank3D, Point<double>> StiffnessMatrix;
	CSSD_Matrix<Tensor2Rank3D, Point<double>> MassMatrix;
	CSSD_Matrix<double, double> MassMatrix_doubleSLAE;
	CSSD_Matrix<double, double> Precond;
	CSSD_Matrix<Tensor2Rank3D, Point<double>> VolumeForceVector;
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
	MassMatrix.SetMatrix(StiffnessMatrix);
	VolumeForceVector.F.resize(StiffnessMatrix.F.size());
	printf_s("\t\tcomplite\n");

	//строим базовые матрицы жесткости/массы/вектор силы тяжести
	//матрица массы идет с коэффициентом rpho
	printf_s("================= Create matrix ================\n");
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
			return solver_grid.GetDomain(solver_grid.GetElement(elem)->GetIdDomain())->forThermal.rpho;
		};
		std::function<Point<double>(int, Point<double>)> VolumeForсe = [&](int elem, Point<double> X)->Point<double>
		{
			return Point<double>(0, 0, 0);
		};

		printf("Matrix assembling...\n");
		std::vector<DenseMatrix<Tensor2Rank3D, Point<double>>> local_SLAE_stiffness(solver_grid.GetElementsCount());
		std::vector<DenseMatrix<Tensor2Rank3D, Point<double>>> local_SLAE_mass(solver_grid.GetElementsCount());
		std::vector< std::vector<Point<double>>> local_force(solver_grid.GetElementsCount());
		omp_set_num_threads(math::NUM_THREADS);
#pragma omp parallel for schedule(dynamic)
		for (int id_elem = 0; id_elem < solver_grid.GetElementsCount(); id_elem++)
		{
			if (id_elem % 1000 == 0)
				printf("Solve element[%d]\r", id_elem);
			auto element = solver_grid.GetElement(id_elem);

			std::function<std::vector<std::vector<double>>(Point<double>)> D = [&](Point<double> X) {return StiffnessCoef(id_elem, X); };
			std::function<double(Point<double>)> M = [&](Point<double> X) {return MassCoef(id_elem, X); };
			std::function<Point<double>(Point<double>)> F = [&](Point<double> X) {return VolumeForсe(id_elem, X); };

			element->SolveLocalMatrix(local_SLAE_stiffness[id_elem], D);
			element->SolveMassMatrix(local_SLAE_mass[id_elem], M);
			element->SolveRightSide(local_force[id_elem], F);
		}
		for (int id_elem = 0; id_elem < solver_grid.GetElementsCount(); id_elem++)
		{
			if (id_elem % 1000 == 0)
				printf("Add the local matrix of element[%d]\r", id_elem);
			StiffnessMatrix.SummPartOfMatrix(local_SLAE_stiffness[id_elem], *solver_grid.GetElementDOFs(id_elem));
			MassMatrix.SummPartOfMatrix(local_SLAE_mass[id_elem], *solver_grid.GetElementDOFs(id_elem));
			VolumeForceVector.SummPartOfVector(local_force[id_elem], *solver_grid.GetElementDOFs(id_elem));
		}
		printf_s("                                                                                    \r");
		printf_s("\t\tcomplite\n");

		//учитываем первые краевые ТОЛЬКО в матрице массы БЕЗ симметризации
		if (false) {
			//обнуляем строки
			for (int id_vertex = 0; id_vertex < solver_grid.boundary_vertexes.size(); id_vertex++)
			{
				auto boundary = &(solver_grid.boundary_vertexes[id_vertex]);

				int id_type = boundary->id_type;
				Point<bool> is_take;
				Point<double> boundary_value = boundary->boundary_value(is_take, solver_grid.boundary_vertexes[id_vertex].GetIdNode(0));
				int global_id = boundary->GetDOFInLocalID(0);

				auto enter_boundary = [global_id, &boundary_value](int position, CSSD_Matrix<Tensor2Rank3D, Point<double>>& global_SLAE) {
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
				};
				if (is_take.x == true) enter_boundary(0, MassMatrix);
				if (is_take.y == true) enter_boundary(1, MassMatrix);
				if (is_take.z == true) enter_boundary(2, MassMatrix);
			}
		}
		if (true) {
			//большое число на диагональ
			for (int id_vertex = 0; id_vertex < solver_grid.boundary_vertexes.size(); id_vertex++)
			{
				auto boundary = &(solver_grid.boundary_vertexes[id_vertex]);

				int id_type = boundary->id_type;
				Point<bool> is_take;
				Point<double> boundary_value = boundary->boundary_value(is_take, solver_grid.boundary_vertexes[id_vertex].GetIdNode(0));
				int global_id = boundary->GetDOFInLocalID(0);

				auto enter_boundary = [global_id, &boundary_value](int position, CSSD_Matrix<Tensor2Rank3D, Point<double>>& global_SLAE) {
					global_SLAE.Diag[global_id].val[position][position] = 1e+30;
				};
				if (is_take.x == true) enter_boundary(0, MassMatrix);
				if (is_take.y == true) enter_boundary(1, MassMatrix);
				if (is_take.z == true) enter_boundary(2, MassMatrix);
			}
		}

		//переводим матрицу массы в действительную
		math::MakeCopyMatrix_A_into_B(MassMatrix, MassMatrix_doubleSLAE);
		Precond.PrecondorSSOR(0.75, MassMatrix_doubleSLAE);
	}

	///-------------------
	Point<double> point_for_out(0, 0, 0.5 - 0.01);
	std::vector<double> Uz_in_point;
	std::vector<double> Time_for_Uz;
	double TIME_L;
	int id_elem_for_out;
	char name_uz[5000];
	{
		double v = solver_grid.GetDomain(0)->forMech.GetLongitudinalWaveVelocity_Vp(solver_grid.GetDomain(0)->forThermal.rpho);
		double L = 0.5;
		TIME_L = L / v;
		double len;
		id_elem_for_out = solver_grid.GetNearestElementID(point_for_out, len);

		wchar_t _tmp_wc[1000];
		math::Char_To_Wchar_t(base_result_directory, _tmp_wc, 1000);
		CreateDirectory((LPCTSTR)_tmp_wc, NULL);

		FILE* fout_Uz;
		sprintf_s(name_uz, "%s/Uz_in_point.txt", base_result_directory);
		fopen_s(&fout_Uz, name_uz, "w");
		fprintf_s(fout_Uz, "Point: (%.2e, %.2e, %.2e)\n", point_for_out.x, point_for_out.y, point_for_out.z);
		fprintf_s(fout_Uz, "Id element: %d\n", id_elem_for_out);
		fprintf_s(fout_Uz, "T0: %.5e sec\n", TIME_L);
		fclose(fout_Uz);

		if (id_elem_for_out < 0) return;
	}
	///-------------------

	TIME_h /= 10000;
	double TIME_curr = TIME_h;
	for (int id_STEP = start_iteration; id_STEP <= end_iteration; id_STEP++)
	{
		if (id_STEP == 3) TIME_h *= 10000;
		TIME_curr += TIME_h;
		printf_s("\n================= Start solution of %d STEP (time = %.2e) ================\n", id_STEP, TIME_curr);
		bool is_print_result = false;
		char result_directory[1000];
		if (id_STEP % step_for_out == 0)
		{
			is_print_result = true;
			sprintf_s(result_directory, sizeof(result_directory), "%s", base_result_directory);
			//sprintf_s(result_directory, sizeof(result_directory), "%s/STEP_%d_t=%.2e", base_result_directory, id_STEP, TIME_curr);
			wchar_t _tmp_wc[1000];
			math::Char_To_Wchar_t(result_directory, _tmp_wc, 1000);
			CreateDirectory((LPCTSTR)_tmp_wc, NULL);
		}

		//переопределяем источник и строим вектор в правую часть (через вторые краевые)
		CSSD_Matrix<Tensor2Rank3D, Point<double>> SourseVector;
		SourseVector.F.resize(StiffnessMatrix.F.size());
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
			std::function<Point<double>(Point<double>)> new_sourse_value_sin = [TIME_curr, power, &solver_grid, TIME_h](Point<double> X) -> Point<double>
			{
				//double sourse = exp(-pow(TIME_curr - tay_0, degree) / pow(tay_plus, degree));
				//return power;
				//return Point<double>(power.x * sourse, power.y * sourse, power.z * sourse);
				double v = solver_grid.GetDomain(0)->forMech.GetLongitudinalWaveVelocity_Vp(solver_grid.GetDomain(0)->forThermal.rpho);
				double L = 0.5;
				double t0 = v / L;
				double w = 2 * M_PI / t0 * 10;
				double alpha = 5e+4;

				double sourse = 10e+6 * exp(-1 * alpha * TIME_curr) * sin(w * TIME_curr);
				if (TIME_curr > t0) sourse = 0;
				return Point<double>(0 * sourse, 0 * sourse, -1 * sourse /** TIME_h / 14.0e-7*/);
			};
			if (second_boundary.size() > 0)
			{
				second_boundary[0].value = new_sourse_value_sin;
				//second_boundary[1].value = new_sourse_value_sin;
				//second_boundary[2].value = new_sourse_value;

				for (int i = 0; i < solver_grid.boundary_faces.size(); i++)
				{
					solver_grid.boundary_faces[i].boundary_value = second_boundary[solver_grid.boundary_faces[i].id_type].value;
				}
			}

			for (int id_triangle = 0; id_triangle < solver_grid.boundary_faces.size(); id_triangle++)
			{
				std::vector<Point<double>> local_vector_SLAE;
				solver_grid.boundary_faces[id_triangle].SolveLocalBoundaryVector(local_vector_SLAE, solver_grid.boundary_faces[id_triangle].boundary_value);
				SourseVector.SummPartOfVector(local_vector_SLAE, *solver_grid.boundary_faces[id_triangle].GetElementDOFs());
			};
		}

		auto SolveSLAE = [&MassMatrix_doubleSLAE, &solver_grid, &Precond](std::vector<double>& F, std::vector<double>& result, double critical_residual) -> void
		{
			//краевые в правую часть + симметризация 
			if (false) {
				//добавляем значения в правую часть
				for (int id_vertex = 0; id_vertex < solver_grid.boundary_vertexes.size(); id_vertex++)
				{
					auto boundary = &(solver_grid.boundary_vertexes[id_vertex]);

					int id_type = boundary->id_type;
					Point<bool> is_take;
					Point<double> boundary_value = boundary->boundary_value(is_take, solver_grid.boundary_vertexes[id_vertex].GetIdNode(0));
					int global_id = boundary->GetDOFInLocalID(0);

					if (global_id < MassMatrix_doubleSLAE.GetMatrixSize() / 3)
					{
						auto enter_value = [&result, &F](int value_id, double value) ->void
						{
							result[value_id] = value;
							F[value_id] = value;
							return;
						};

						if (is_take.x) enter_value(global_id * 3 + 0, boundary_value.x);
						if (is_take.y) enter_value(global_id * 3 + 1, boundary_value.y);
						if (is_take.z) enter_value(global_id * 3 + 2, boundary_value.z);
					}

				}
				//симметризация
				omp_set_num_threads(math::NUM_THREADS);
#pragma omp parallel for schedule(dynamic) 
				for (int id_row = 0; id_row < MassMatrix_doubleSLAE.GetMatrixSize(); id_row++)
				{
					if (id_row % 100 == 0)
					{
						printf("\tcurrent %d/%d\r", id_row, MassMatrix_doubleSLAE.GetMatrixSize());
					}
					int iterator_in_boundary = 0;
					for (int jj = 0; jj < MassMatrix_doubleSLAE.id_column_for_A_up[id_row].size(); jj++)
					{
						int id_column = MassMatrix_doubleSLAE.id_column_for_A_up[id_row][jj];
						for (/*iterator_in_boundary = 0*/; iterator_in_boundary < solver_grid.boundary_vertexes.size(); iterator_in_boundary++)
						{
							if (solver_grid.boundary_vertexes[iterator_in_boundary].GetDOFInLocalID(0) * 3 + 0 == id_column)
							{
								Point<bool> is_take;
								Point<double> boundary_value = solver_grid.boundary_vertexes[iterator_in_boundary].boundary_value(is_take, solver_grid.boundary_vertexes[iterator_in_boundary].GetIdNode(0));

								if (is_take.x)
								{
									F[id_row] -= MassMatrix_doubleSLAE.A_up[id_row][jj] * F[id_column];
									MassMatrix_doubleSLAE.A_up[id_row][jj] = 0;
								}

								break;
							}
							if (solver_grid.boundary_vertexes[iterator_in_boundary].GetDOFInLocalID(0) * 3 + 1 == id_column)
							{
								Point<bool> is_take;
								Point<double> boundary_value = solver_grid.boundary_vertexes[iterator_in_boundary].boundary_value(is_take, solver_grid.boundary_vertexes[iterator_in_boundary].GetIdNode(0));

								if (is_take.y)
								{
									F[id_row] -= MassMatrix_doubleSLAE.A_up[id_row][jj] * F[id_column];
									MassMatrix_doubleSLAE.A_up[id_row][jj] = 0;
								}

								break;
							}
							if (solver_grid.boundary_vertexes[iterator_in_boundary].GetDOFInLocalID(0) * 3 + 2 == id_column)
							{
								Point<bool> is_take;
								Point<double> boundary_value = solver_grid.boundary_vertexes[iterator_in_boundary].boundary_value(is_take, solver_grid.boundary_vertexes[iterator_in_boundary].GetIdNode(0));

								if (is_take.z)
								{
									F[id_row] -= MassMatrix_doubleSLAE.A_up[id_row][jj] * F[id_column];
									MassMatrix_doubleSLAE.A_up[id_row][jj] = 0;
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
					for (int jj = 0; jj < MassMatrix_doubleSLAE.id_column_for_A_down[id_row].size(); jj++)
					{
						int id_column = MassMatrix_doubleSLAE.id_column_for_A_down[id_row][jj];
						for (/*iterator_in_boundary = 0*/; iterator_in_boundary < solver_grid.boundary_vertexes.size(); iterator_in_boundary++)
						{
							if (solver_grid.boundary_vertexes[iterator_in_boundary].GetDOFInLocalID(0) * 3 + 0 == id_column)
							{
								Point<bool> is_take;
								Point<double> boundary_value = solver_grid.boundary_vertexes[iterator_in_boundary].boundary_value(is_take, solver_grid.boundary_vertexes[iterator_in_boundary].GetIdNode(0));

								if (is_take.x)
								{
									F[id_row] -= MassMatrix_doubleSLAE.A_down[id_row][jj] * F[id_column];
									MassMatrix_doubleSLAE.A_down[id_row][jj] = 0;
								}

								break;
							}
							if (solver_grid.boundary_vertexes[iterator_in_boundary].GetDOFInLocalID(0) * 3 + 1 == id_column)
							{
								Point<bool> is_take;
								Point<double> boundary_value = solver_grid.boundary_vertexes[iterator_in_boundary].boundary_value(is_take, solver_grid.boundary_vertexes[iterator_in_boundary].GetIdNode(0));

								if (is_take.y)
								{
									F[id_row] -= MassMatrix_doubleSLAE.A_down[id_row][jj] * F[id_column];
									MassMatrix_doubleSLAE.A_down[id_row][jj] = 0;
								}

								break;
							}
							if (solver_grid.boundary_vertexes[iterator_in_boundary].GetDOFInLocalID(0) * 3 + 2 == id_column)
							{
								Point<bool> is_take;
								Point<double> boundary_value = solver_grid.boundary_vertexes[iterator_in_boundary].boundary_value(is_take, solver_grid.boundary_vertexes[iterator_in_boundary].GetIdNode(0));

								if (is_take.z)
								{
									F[id_row] -= MassMatrix_doubleSLAE.A_down[id_row][jj] * F[id_column];
									MassMatrix_doubleSLAE.A_down[id_row][jj] = 0;
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
			//краевые в правую часть большим числом 
			if (true) {
				//добавляем значения в правую часть
				for (int id_vertex = 0; id_vertex < solver_grid.boundary_vertexes.size(); id_vertex++)
				{
					auto boundary = &(solver_grid.boundary_vertexes[id_vertex]);

					int id_type = boundary->id_type;
					Point<bool> is_take;
					Point<double> boundary_value = boundary->boundary_value(is_take, solver_grid.boundary_vertexes[id_vertex].GetIdNode(0));
					int global_id = boundary->GetDOFInLocalID(0);

					if (global_id < MassMatrix_doubleSLAE.GetMatrixSize() / 3)
					{
						auto enter_value = [&result, &F, &MassMatrix_doubleSLAE](int value_id, double value) ->void
						{
							result[value_id] = value;
							F[value_id] = value * MassMatrix_doubleSLAE.Diag[value_id];
							return;
						};

						if (is_take.x) enter_value(global_id * 3 + 0, boundary_value.x);
						if (is_take.y) enter_value(global_id * 3 + 1, boundary_value.y);
						if (is_take.z) enter_value(global_id * 3 + 2, boundary_value.z);
					}

				}
			}

			printf("Soluting SLAY... (%d)\n", MassMatrix_doubleSLAE.GetMatrixSize());
			int MaxSize = MassMatrix_doubleSLAE.GetMatrixSize();
			MaxSize = MaxSize / 10 < 1000 ? 1000 : MaxSize / 10;
			std::vector<double> best_solution;
			math::MakeCopyVector_A_into_B(result, best_solution);
			double current_residual = 1, best_residual = 1e+25;
			int MAX_STEPS = 4;
			int ii = 1;

			for (int i = 0; i <= MAX_STEPS; i++)
			{
				printf_s("//---> I = %d/%d (full size %d) - needed residual %.2e\n", i, MAX_STEPS, MassMatrix_doubleSLAE.GetMatrixSize(), critical_residual);

				MassMatrix_doubleSLAE.print_logs = true;
				current_residual = abs(MassMatrix_doubleSLAE.MSG_PreconditioningSSOR(F, result, MaxSize, critical_residual, Precond));

				if (current_residual <= critical_residual)
				{
					best_residual = current_residual;
					math::MakeCopyVector_A_into_B(result, best_solution);
					break;
				}
				if (current_residual < best_residual)
				{
					best_residual = current_residual;
					math::MakeCopyVector_A_into_B(result, best_solution);
				}
			}
			math::MakeCopyVector_A_into_B(best_solution, result);

			if (best_residual >= 1)
			{
				printf_s("There are problems in the solution\n");
				int a;
				scanf_s("%d", &a);
			}
		};

		//первые итерации Рунге-Кутта, дальше прогноз-коррекция
		if (id_STEP <= 2)
		{
			std::vector<double> K[4];
			std::vector<double> L[4];
			std::vector<Point<double>> Ku0_prev;
			std::vector<double> F_slae(U_prev.size() * 3);
			
			//M*u''=F-Ku0 => u''=K[0]
			{
				StiffnessMatrix.MultiplicationMatrixVector(U_prev, Ku0_prev);
				for (int i = 0; i < Ku0_prev.size(); i++)
				{
					Point<double> res = Point<double>(0, 0, 0)
						- Ku0_prev[i]
						+ VolumeForceVector.F[i]
						+ SourseVector.F[i];
					F_slae[i * 3 + 0] = res.x;
					F_slae[i * 3 + 1] = res.y;
					F_slae[i * 3 + 2] = res.z;
				}
				K[0].resize(F_slae.size());
				math::InitializationVector(K[0], 1e-10);
				SolveSLAE(F_slae, K[0], 1e-10);

				L[0].resize(F_slae.size());
				for (int i = 0; i < dU_dt_prev.size(); i++)
				{
					K[0][i * 3 + 0] *= TIME_h;
					K[0][i * 3 + 1] *= TIME_h;
					K[0][i * 3 + 2] *= TIME_h;

					L[0][i * 3 + 0] = dU_dt_prev[i].x * TIME_h;
					L[0][i * 3 + 1] = dU_dt_prev[i].y * TIME_h;
					L[0][i * 3 + 2] = dU_dt_prev[i].z * TIME_h;
				}
			}

			//M*u''=F-K(u0+L0*h/2) => u''=K1
			{
				std::vector<Point<double>> U_tmp(U_prev.size());
				for (int i = 0; i < U_prev.size(); i++)
					U_tmp[i] = U_prev[i] + Point<double>(L[0][i * 3 + 0], L[0][i * 3 + 1], L[0][i * 3 + 2]) / 2.0;
				StiffnessMatrix.MultiplicationMatrixVector(U_tmp, Ku0_prev);
				for (int i = 0; i < Ku0_prev.size(); i++)
				{
					Point<double> res = Point<double>(0, 0, 0)
						- Ku0_prev[i]
						+ VolumeForceVector.F[i]
						+ SourseVector.F[i];
					F_slae[i * 3 + 0] = res.x;
					F_slae[i * 3 + 1] = res.y;
					F_slae[i * 3 + 2] = res.z;
				}
				math::MakeCopyVector_A_into_B(K[0], K[1]);
				SolveSLAE(F_slae, K[1], 1e-10);

				L[1].resize(F_slae.size());
				for (int i = 0; i < dU_dt_prev.size(); i++)
				{
					K[1][i * 3 + 0] *= TIME_h;
					K[1][i * 3 + 1] *= TIME_h;
					K[1][i * 3 + 2] *= TIME_h;

					L[1][i * 3 + 0] = (dU_dt_prev[i].x + K[0][i * 3 + 0] / 2.0) * TIME_h;
					L[1][i * 3 + 1] = (dU_dt_prev[i].y + K[0][i * 3 + 1] / 2.0) * TIME_h;
					L[1][i * 3 + 2] = (dU_dt_prev[i].z + K[0][i * 3 + 3] / 2.0) * TIME_h;
				}
			}

			//M*u''=F-K(u0+K1*h/2) => u''=K2
			{
				std::vector<Point<double>> U_tmp(U_prev.size());
				for (int i = 0; i < U_prev.size(); i++)
					U_tmp[i] = U_prev[i] + Point<double>(L[1][i * 3 + 0], L[1][i * 3 + 1], L[1][i * 3 + 2]) / 2.0;
				StiffnessMatrix.MultiplicationMatrixVector(U_tmp, Ku0_prev);
				for (int i = 0; i < Ku0_prev.size(); i++)
				{
					Point<double> res = Point<double>(0, 0, 0)
						- Ku0_prev[i]
						+ VolumeForceVector.F[i]
						+ SourseVector.F[i];
					F_slae[i * 3 + 0] = res.x;
					F_slae[i * 3 + 1] = res.y;
					F_slae[i * 3 + 2] = res.z;
				}
				math::MakeCopyVector_A_into_B(K[1], K[2]);
				SolveSLAE(F_slae, K[2], 1e-10);

				L[2].resize(F_slae.size());
				for (int i = 0; i < dU_dt_prev.size(); i++)
				{
					K[2][i * 3 + 0] *= TIME_h;
					K[2][i * 3 + 1] *= TIME_h;
					K[2][i * 3 + 2] *= TIME_h;

					L[2][i * 3 + 0] = (dU_dt_prev[i].x + K[1][i * 3 + 0] / 2.0) * TIME_h;
					L[2][i * 3 + 1] = (dU_dt_prev[i].y + K[1][i * 3 + 1] / 2.0) * TIME_h;
					L[2][i * 3 + 2] = (dU_dt_prev[i].z + K[1][i * 3 + 3] / 2.0) * TIME_h;
				}
			}

			//M*u''=F-K(u0+K2*h) => u''=K3
			{
				std::vector<Point<double>> U_tmp(U_prev.size());
				for (int i = 0; i < U_prev.size(); i++)
					U_tmp[i] = U_prev[i] + Point<double>(L[2][i * 3 + 0], L[2][i * 3 + 1], L[2][i * 3 + 2]);
				StiffnessMatrix.MultiplicationMatrixVector(U_tmp, Ku0_prev);
				for (int i = 0; i < Ku0_prev.size(); i++)
				{
					Point<double> res = Point<double>(0, 0, 0)
						- Ku0_prev[i]
						+ VolumeForceVector.F[i]
						+ SourseVector.F[i];
					F_slae[i * 3 + 0] = res.x;
					F_slae[i * 3 + 1] = res.y;
					F_slae[i * 3 + 2] = res.z;
				}
				math::MakeCopyVector_A_into_B(K[2], K[3]);
				SolveSLAE(F_slae, K[3], 1e-10);

				L[3].resize(F_slae.size());
				for (int i = 0; i < dU_dt_prev.size(); i++)
				{
					K[3][i * 3 + 0] *= TIME_h;
					K[3][i * 3 + 1] *= TIME_h;
					K[3][i * 3 + 2] *= TIME_h;

					L[3][i * 3 + 0] = (dU_dt_prev[i].x + K[2][i * 3 + 0]) * TIME_h;
					L[3][i * 3 + 1] = (dU_dt_prev[i].y + K[2][i * 3 + 1]) * TIME_h;
					L[3][i * 3 + 2] = (dU_dt_prev[i].z + K[2][i * 3 + 3]) * TIME_h;
				}
			}

			for (int i = 0; i < U_curr.size(); i++)
			{
				dU_dt_curr[i].x = dU_dt_curr[i].x + (K[0][i * 3 + 0] + 2 * K[1][i * 3 + 0] + 2 * K[2][i * 3 + 0] + K[3][i * 3 + 0]) / 6.0;
				dU_dt_curr[i].y = dU_dt_curr[i].y + (K[0][i * 3 + 1] + 2 * K[1][i * 3 + 1] + 2 * K[2][i * 3 + 1] + K[3][i * 3 + 1]) / 6.0;
				dU_dt_curr[i].z = dU_dt_curr[i].z + (K[0][i * 3 + 2] + 2 * K[1][i * 3 + 2] + 2 * K[2][i * 3 + 2] + K[3][i * 3 + 2]) / 6.0;

				U_curr[i].x = U_curr[i].x + (L[0][i * 3 + 0] + 2 * L[1][i * 3 + 0] + 2 * L[2][i * 3 + 0] + L[3][i * 3 + 0]) / 6.0;
				U_curr[i].y = U_curr[i].y + (L[0][i * 3 + 1] + 2 * L[1][i * 3 + 1] + 2 * L[2][i * 3 + 1] + L[3][i * 3 + 1]) / 6.0;
				U_curr[i].z = U_curr[i].z + (L[0][i * 3 + 2] + 2 * L[1][i * 3 + 2] + 2 * L[2][i * 3 + 2] + L[3][i * 3 + 2]) / 6.0;
			}
		}
		else
		{
			std::vector<Point<double>> Ku;
			std::vector<Point<double>> Phi_i(U_prev.size());
			std::vector<double> F_slae(U_prev.size() * 3);

			//M*u''=F-Ku_i => u''=phi_i
			{
				StiffnessMatrix.MultiplicationMatrixVector(U_prev, Ku);
				for (int i = 0; i < Ku.size(); i++)
				{
					Point<double> res = Point<double>(0, 0, 0)
						- Ku[i]
						+ VolumeForceVector.F[i]
						+ SourseVector.F[i];
					F_slae[i * 3 + 0] = res.x;
					F_slae[i * 3 + 1] = res.y;
					F_slae[i * 3 + 2] = res.z;
				}
				math::InitializationVector(MassMatrix_doubleSLAE.X, 1e-10);
				SolveSLAE(F_slae, MassMatrix_doubleSLAE.X, 1e-10);
			}

			//прогноз
			for (int i = 0; i < U_pred.size(); i++)
			{
				U_pred[i] = U_prevprev[i] + dU_dt_prev[i] * 2 * TIME_h;

				Phi_i[i] = Point<double>(MassMatrix_doubleSLAE.X[i * 3 + 0], MassMatrix_doubleSLAE.X[i * 3 + 1], MassMatrix_doubleSLAE.X[i * 3 + 2]);
				dU_dt_pred[i] = dU_dt_prevprev[i] + Phi_i[i] * 2 * TIME_h;
			}

			//M*u''=F-Ku_pred => u''=phi_pred
			{
				StiffnessMatrix.MultiplicationMatrixVector(U_pred, Ku);
				for (int i = 0; i < Ku.size(); i++)
				{
					Point<double> res = Point<double>(0, 0, 0)
						- Ku[i]
						+ VolumeForceVector.F[i]
						+ SourseVector.F[i];
					F_slae[i * 3 + 0] = res.x;
					F_slae[i * 3 + 1] = res.y;
					F_slae[i * 3 + 2] = res.z;
				}
				//math::InitializationVector(MassMatrix_doubleSLAE.X, 1e-10);
				SolveSLAE(F_slae, MassMatrix_doubleSLAE.X, 1e-10);
			}

			//коррекция
			for (int i = 0; i < U_pred.size(); i++)
			{
				U_curr[i] = U_prev[i] + (dU_dt_prev[i] + dU_dt_pred[i]) * TIME_h / 2.0;

				Point<double> phi_pred(MassMatrix_doubleSLAE.X[i * 3 + 0], MassMatrix_doubleSLAE.X[i * 3 + 1], MassMatrix_doubleSLAE.X[i * 3 + 2]);
				dU_dt_curr[i] = dU_dt_prev[i] + (Phi_i[i] + phi_pred) * TIME_h / 2.0;
			}
		}

		//вывод сейсмограммы
		if (TIME_curr > TIME_L || true) {
			Point<double> U_curr_in_point = solver_grid.GetSolutionInPoint(id_elem_for_out, point_for_out, U_curr);
			Point<double> U_prev_in_point = solver_grid.GetSolutionInPoint(id_elem_for_out, point_for_out, U_prev);
			Uz_in_point.push_back((U_curr_in_point.z + U_prev_in_point.z) / 2.0);
			Time_for_Uz.push_back(TIME_curr - TIME_h / 2.0);
		}

		//output solution
		//without deformations
		if (is_print_result) {
			printf_s("Print the mech result into .dat file... ");

			FILE* fout_tech;
			char name_u_tech[5000];
			sprintf_s(name_u_tech, "%s/U_deformation_%d_t%.2e.dat", result_directory, id_STEP, TIME_curr);
			fopen_s(&fout_tech, name_u_tech, "w");
			char name_in_file[1000];
			sprintf_s(name_in_file, "step%d_time%.2e", id_STEP, TIME_curr);
			std::vector<std::vector<char>> name_value(6);
			char name_v_tmp[6][100];
			sprintf_s(name_v_tmp[0], "sigma_Mises");
			sprintf_s(name_v_tmp[1], "sigma_zz");
			sprintf_s(name_v_tmp[2], "Ux");
			sprintf_s(name_v_tmp[3], "Uy");
			sprintf_s(name_v_tmp[4], "Uz");
			sprintf_s(name_v_tmp[5], "Vz");
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
			double sigma_inv_max = 0;
			int elem_sigma_max = 0;
			for (int i = 0; i < solver_grid.GetElementsCount(); i++)
			{
				auto element = solver_grid.GetElement(i);

				Point<double> Centr = element->GetWeightCentr();

				//Point<double> U = solver_grid.GetSolutionInPoint(i, Centr, U_curr);
				Point<Point<double>> dU = solver_grid.GetDerevativeFromSolutionInPoint(i, Centr, U_curr);
				auto Eps = solver_grid.GetStrainTensorFromSolutionInPoint(i, Centr, dU);
				auto Sigma = solver_grid.GetStressTensorFromSolutionInPoint(i, Centr, Eps);
				auto MisesSigma = solver_grid.GetVonMisesStress(Sigma);

				auto Eps_inv = sqrt((Eps.val[0][0] - Eps.val[1][1]) * (Eps.val[0][0] - Eps.val[1][1])
					+ (Eps.val[1][1] - Eps.val[2][2]) * (Eps.val[1][1] - Eps.val[2][2])
					+ (Eps.val[0][0] - Eps.val[2][2]) * (Eps.val[0][0] - Eps.val[2][2])
					+ 3 * (Eps.val[0][1] * Eps.val[1][0] + Eps.val[0][2] * Eps.val[2][0] + Eps.val[1][2] * Eps.val[2][1]) / 2.0) * sqrt(2.) / 3.;

				value[0][i] = MisesSigma;
				value[1][i] = Sigma.val[2][2];
			}
			value[2].resize(solver_grid.GetDOFsCount());
			value[3].resize(solver_grid.GetDOFsCount());
			value[4].resize(solver_grid.GetDOFsCount());
			value[5].resize(solver_grid.GetDOFsCount());
			for (int i = 0; i < value[3].size(); i++)
			{
				value[2][i] = U_curr[i].x;
				value[3][i] = U_curr[i].y;
				value[4][i] = U_curr[i].z;
				value[5][i] = dU_dt_curr[i].z;
			}

			//solver_grid.MoveCoordinates(U_curr);
			solver_grid.printTecPlot3D_DiffDomains(fout_tech, value, name_value, name_in_file, TIME_curr);
			//solver_grid.ReMoveCoordinates(U_curr);
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
			sprintf_s(name_u_tech, "%s/U_plane_YZ_s%d_t%.2e.dat", result_directory, id_STEP, TIME_curr);
			fopen_s(&fout_tech, name_u_tech, "w");
			char name_in_file[1000];
			sprintf_s(name_in_file, "Time_%.4e_YZ", TIME_curr);
			std::vector<std::vector<char>> name_value(6);
			char name_v_tmp[6][100];
			sprintf_s(name_v_tmp[0], "sigma_Mises");
			sprintf_s(name_v_tmp[1], "sigma_zz");
			sprintf_s(name_v_tmp[2], "Ux");
			sprintf_s(name_v_tmp[3], "Uy");
			sprintf_s(name_v_tmp[4], "Uz");
			sprintf_s(name_v_tmp[5], "Vz");
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
					Point<double> dU_dt = solver_grid.GetSolutionInPoint(id_elem, Centr, dU_dt_curr);
					Point<Point<double>> dU = solver_grid.GetDerevativeFromSolutionInPoint(id_elem, Centr, U_curr);
					auto Eps = solver_grid.GetStrainTensorFromSolutionInPoint(id_elem, Centr, dU);
					auto Sigma = solver_grid.GetStressTensorFromSolutionInPoint(id_elem, Centr, Eps);
					auto MisesSigma = solver_grid.GetVonMisesStress(Sigma);

					auto Eps_inv = sqrt((Eps.val[0][0] - Eps.val[1][1]) * (Eps.val[0][0] - Eps.val[1][1])
						+ (Eps.val[1][1] - Eps.val[2][2]) * (Eps.val[1][1] - Eps.val[2][2])
						+ (Eps.val[0][0] - Eps.val[2][2]) * (Eps.val[0][0] - Eps.val[2][2])
						+ 3 * (Eps.val[0][1] * Eps.val[1][0] + Eps.val[0][2] * Eps.val[2][0] + Eps.val[1][2] * Eps.val[2][1]) / 2.0) * sqrt(2.) / 3.;

					value[0][i] = MisesSigma;
					value[1][i] = Sigma.val[2][2];

					value[2][i] = U.x;
					value[3][i] = U.y;
					value[4][i] = U.z;
					value[5][i] = dU_dt.z;
				}
			}

			grid_YZ_plane.printTecPlot3D(fout_tech, value, name_value, name_in_file, TIME_curr);
			fclose(fout_tech);
		}
		//Uz in point
		if (is_print_result && TIME_curr > TIME_L)
		{
			FILE* fout_Uz;
			fopen_s(&fout_Uz, name_uz, "a");
			for (int i = 0; i < Uz_in_point.size(); i++)
			{
				fprintf_s(fout_Uz, "%.10e %.10e\n", Time_for_Uz[i], Uz_in_point[i]);
			}
			Time_for_Uz.clear();
			Uz_in_point.clear();
			fclose(fout_Uz);
		}

		printf_s("\t complite\n");

		math::MakeCopyVector_A_into_B(U_prev, U_prevprev);
		math::MakeCopyVector_A_into_B(U_curr, U_prev);
		math::MakeCopyVector_A_into_B(dU_dt_prev, dU_dt_prevprev);
		math::MakeCopyVector_A_into_B(dU_dt_curr, dU_dt_prev);

		//обновляем сетку
		if (false) {
			solver_grid.MoveCoordinates(U_curr);
			for (int i = 0; i < solver_grid.GetElementsCount(); i++)
			{
				solver_grid.GetElement(i)->SolveAlphaMatrix();
			}
		}
	}
}
void ElastodynamicsProblem_ImplicitNewmark_fast(char properties_file[1000])
{
	FEM::Grid_forMech solver_grid; //output
	std::vector<Point<double>> Solution; //output

	double BETTA_Newmark = 1.0 / 4.0;
	double GAMMA_Newmark = 1.0 / 2.0;

	bool is_STATIONARY = false;


	//char properties_file[1000] = { "E:/+cyl/800el/param.txt" };
	//char properties_file[1000] = { "E:/Box/100x10x200/BoxWithCracks/67k/param.txt" };
	//char properties_file[1000] = { "E:/Box/100x10x200/BoxWithCracks/638k/param.txt" };
	//char properties_file[1000] = { "E:/Box/100x50x200/param.txt" };
	//char properties_file[1000] = { "E:/Box/200x200x200/200k_fem/param_for_solver.txt" };
	//char properties_file[1000] = { "./param_for_solver.txt" };
	//char properties_file[1000] = { "D:/Elastodynamic/homocyl/67k/param_for_solver.txt" };
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


	{
		char _line[1000];
		std::vector<int> val;

		math::ReadNonEmptyLine_forNumbers(f_properties, _line);
		math::ParserStringToVectorInt(_line, val, " ");
		if (val[0] == 1)
			is_print_logFile = true;


		math::ReadNonEmptyLine_forNumbers(f_properties, _line);
		math::ParserStringToVectorInt(_line, val, " ");
		math::NUM_THREADS = val[0];
	}

	math::ReadNonEmptyLine(f_properties, mesh_directory);
	math::ReadNonEmptyLine(f_properties, base_result_directory);
	math::SimpleGrid geo_grid; //input
	std::vector<std::vector<std::vector<int>>> boundary_faces;
	geo_grid.ReadFromNVTR(mesh_directory, 4);
	geo_grid.ReadFacesBoundaryNVTR(mesh_directory, boundary_faces);

	struct Receiver
	{
		Point<double> point;
		std::vector<Point<double>> U_in_point;
		int id_node;
		char file_name[1000];
	};
	std::vector<Receiver> receivers;
	std::vector<double> receivers_times;
	/*{
		char _line[1000];
		math::ReadNonEmptyLine_forNumbers(f_properties, _line);
		std::vector<double> val;
		math::ParserStringToVectorDouble(_line, val, " ");

		int N = (int)val[0];

		receivers.resize(N);
		for (int i = 0; i < N; i++)
		{
			receivers[i].point.x = val[1 + 0 + i * 3];
			receivers[i].point.y = val[1 + 1 + i * 3];
			receivers[i].point.z = val[1 + 2 + i * 3];
		}
	}*/
	{
		char _line[1000];
		math::ReadNonEmptyLine_forNumbers(f_properties, _line);
		std::vector<std::string> string_val;
		math::ParserStringToVectorString(_line, string_val, " ");

		for (int group = 0; group < string_val.size() / 3; group++)
		{
			std::vector<double> receivers_x;
			std::vector<double> receivers_y;
			std::vector<double> receivers_z;

			math::ParserStringToDouble(string_val[group], receivers_x);
			math::ParserStringToDouble(string_val[group+1], receivers_y);
			math::ParserStringToDouble(string_val[group+2], receivers_z);

			for (int k = 0; k < receivers_z.size(); k++)
			{
				for (int j = 0; j < receivers_y.size(); j++)
				{
					for (int i = 0; i < receivers_x.size(); i++)
					{
						Receiver rec;
						rec.point.x = receivers_x[i];
						rec.point.y = receivers_y[j];
						rec.point.z = receivers_z[k];

						receivers.push_back(rec);
					}
				}
			}
		}

	}
	Point<double> crack_top;
	Point<double> crack_bottom;
	{
		char _line[1000];
		math::ReadNonEmptyLine_forNumbers(f_properties, _line);
		std::vector<double> val;
		math::ParserStringToVectorDouble(_line, val, " ");

		crack_bottom.x = val[0];
		crack_bottom.y = val[1];
		crack_bottom.z = val[2];

		crack_top.x = val[3];
		crack_top.y = val[4];
		crack_top.z = val[5];
	}

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
	struct Gauss_source {
		bool is_used = false;
		int target_boundary = 0;
		int degree = 0;
		double tay = 0;
		double t0 = 0;
		double power = 0;
		Point<double> direction;

		Point<double> value(double Time_curr)
		{
			double arg = (Time_curr - t0) / tay;
			return direction * power * exp(-1.0 * pow(arg, degree));
		}
	} gauss_source;
	struct LinearModulationSin_source
	{
		bool is_used = false;
		int target_boundary = 0;

		double w_min = 1; //Hz
		double w_max = 1000; //Hz
		double tay = 1; //signal duration
		double t0; //start signal
		double power = 1;
		Point<double> direction;

		Point<double> value(double Time_curr)
		{
			double modulation = w_min + (w_max - w_min) * (Time_curr - t0) / tay;
			if (Time_curr <= t0 + tay)
				return direction * power * sin(2 * M_PI * Time_curr * modulation);

			return direction * power * 0.0;
		}
	} linearModulationSin_source;
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
			int id_boundary = math::ParserCharToInt(line[0]);

			if (line[0] == '-')
				id_boundary = math::ParserCharToInt(line[1]) * -1;

			id_boundaries[i] = abs(id_boundary);
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

					//special boundary
					if (id_boundary < 0)
					{
						switch ((int)(val[1]))
						{
						case 1:
							gauss_source.is_used = true;
							gauss_source.power = val[2];
							gauss_source.direction.x = val[3];
							gauss_source.direction.y = val[4];
							gauss_source.direction.z = val[5];
							gauss_source.degree = (int)val[6];
							gauss_source.t0 = val[7];
							gauss_source.tay = val[8];
							gauss_source.target_boundary = i;
							break;
						case 2:
							linearModulationSin_source.is_used = true;
							linearModulationSin_source.power = val[2];
							linearModulationSin_source.direction.x = val[3];
							linearModulationSin_source.direction.y = val[4];
							linearModulationSin_source.direction.z = val[5];
							linearModulationSin_source.t0 = val[6];
							linearModulationSin_source.tay = val[7];
							linearModulationSin_source.w_min = val[8];
							linearModulationSin_source.w_max = val[9];
							linearModulationSin_source.target_boundary = i;
							break;
						default:
							break;
						}
						
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

	{
		N_domain++;
		double _E = 0;
		double _v = 0;
		double _rpho = 0;
		_materials.push_back(_material(_E, _v, _rpho));

		for (int i = 0; i < geo_grid.nvtr.size(); i++)
		{
			Point<double> centr;
			for (int j = 0; j < geo_grid.nvtr[i].size(); j++)
			{
				centr += geo_grid.xyz[geo_grid.nvtr[i][j]] / geo_grid.nvtr[i].size();
			}
			if (crack_bottom.x < centr.x && centr.x < crack_top.x &&
				crack_bottom.y < centr.y && centr.y < crack_top.y &&
				crack_bottom.z < centr.z && centr.z < crack_top.z)
			{
				geo_grid.nvkat[i] = N_domain - 1;
			}
		}
	}

	fclose(f_properties);

	wchar_t _tmp_wc[1000];
	math::Char_To_Wchar_t(base_result_directory, _tmp_wc, 1000);
	CreateDirectory((LPCTSTR)_tmp_wc, NULL);

	bool res = false;
	double TIME_h = step_size;
	if (is_STATIONARY) {
		end_iteration = 1;
		start_iteration = 0;
	}
	std::vector<Point<double>> U_prev(geo_grid.xyz.size()), U_prevprev(geo_grid.xyz.size()), U_curr(geo_grid.xyz.size());
	std::vector<double> U_prev_d(geo_grid.xyz.size()*3), U_prevprev_d(geo_grid.xyz.size()*3);
	math::InitializationVector(U_curr, 0);
	math::InitializationVector(U_prev, 0);
	math::InitializationVector(U_prevprev, 0);

	CSSD_Matrix<Tensor2Rank3D, Point<double>> StiffnessMatrix;
	CSSD_Matrix<Tensor2Rank3D, Point<double>> MassMatrix;
	CSSD_Matrix<double, double> MassMatrix_doubleSLAE;
	CSSD_Matrix<double, double> StiffnessMatrix_doubleSLAE;
	CSSD_Matrix<double, double> SLAE_n, Tmp_matrix;
	CSSD_Matrix<Tensor2Rank3D, Point<double>> VolumeForceVector;
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
	MassMatrix.SetMatrix(StiffnessMatrix);
	VolumeForceVector.F.resize(StiffnessMatrix.F.size());
	printf_s("\t\tcomplite\n");

	std::vector<int> clear_vertexes;

	//строим базовые матрицы жесткости/массы/вектор силы тяжести
	//матрица массы идет с коэффициентом rpho
	printf_s("================= Create matrix ================\n");
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
			return solver_grid.GetDomain(solver_grid.GetElement(elem)->GetIdDomain())->forThermal.rpho;
		};
		std::function<Point<double>(int, Point<double>)> VolumeForсe = [&](int elem, Point<double> X)->Point<double>
		{
			return Point<double>(0, 0, 0);
		};

		printf("Matrix assembling...\n");
		std::vector<DenseMatrix<Tensor2Rank3D, Point<double>>> local_SLAE_stiffness(solver_grid.GetElementsCount());
		std::vector<DenseMatrix<Tensor2Rank3D, Point<double>>> local_SLAE_mass(solver_grid.GetElementsCount());
		std::vector< std::vector<Point<double>>> local_force(solver_grid.GetElementsCount());
		omp_set_num_threads(math::NUM_THREADS);
#pragma omp parallel for schedule(dynamic)
		for (int id_elem = 0; id_elem < solver_grid.GetElementsCount(); id_elem++)
		{
			if (id_elem % 1000 == 0)
				printf("Solve element[%d]\r", id_elem);
			auto element = solver_grid.GetElement(id_elem);

			if (element->GetIdDomain() != solver_grid.GetDomainsCount() - 1)
			{
				std::function<std::vector<std::vector<double>>(Point<double>)> D = [&](Point<double> X) {return StiffnessCoef(id_elem, X); };
				std::function<double(Point<double>)> M = [&](Point<double> X) {return MassCoef(id_elem, X); };
				std::function<Point<double>(Point<double>)> F = [&](Point<double> X) {return VolumeForсe(id_elem, X); };

				element->SolveLocalMatrix(local_SLAE_stiffness[id_elem], D);
				element->SolveMassMatrix(local_SLAE_mass[id_elem], M);
				element->SolveRightSide(local_force[id_elem], F);
			}
			else
			{
				local_SLAE_mass[id_elem].SetSize(element->GetDOFsCount());
				local_SLAE_stiffness[id_elem].SetSize(element->GetDOFsCount());
				local_force[id_elem].resize(element->GetDOFsCount());
			}
		}
		for (int id_elem = 0; id_elem < solver_grid.GetElementsCount(); id_elem++)
		{
			if (id_elem % 1000 == 0)
				printf("Add the local matrix of element[%d]\r", id_elem);
			StiffnessMatrix.SummPartOfMatrix(local_SLAE_stiffness[id_elem], *solver_grid.GetElementDOFs(id_elem));
			MassMatrix.SummPartOfMatrix(local_SLAE_mass[id_elem], *solver_grid.GetElementDOFs(id_elem));
			VolumeForceVector.SummPartOfVector(local_force[id_elem], *solver_grid.GetElementDOFs(id_elem));
		}
		printf_s("                                                                                    \r");
		printf_s("\t\tcomplite\n");

		//учитываем первые краевые ТОЛЬКО в матрице массы БЕЗ симметризации
		if (false) {
			//обнуляем строки
			for (int id_vertex = 0; id_vertex < solver_grid.boundary_vertexes.size(); id_vertex++)
			{
				auto boundary = &(solver_grid.boundary_vertexes[id_vertex]);

				int id_type = boundary->id_type;
				Point<bool> is_take;
				Point<double> boundary_value = boundary->boundary_value(is_take, solver_grid.boundary_vertexes[id_vertex].GetIdNode(0));
				int global_id = boundary->GetDOFInLocalID(0);

				auto enter_boundary = [global_id, &boundary_value](int position, CSSD_Matrix<Tensor2Rank3D, Point<double>>& global_SLAE) {
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
				};
				if (is_take.x == true) enter_boundary(0, MassMatrix);
				if (is_take.y == true) enter_boundary(1, MassMatrix);
				if (is_take.z == true) enter_boundary(2, MassMatrix);
			}
		}
		if (false) {
			//большое число на диагональ
			for (int id_vertex = 0; id_vertex < solver_grid.boundary_vertexes.size(); id_vertex++)
			{
				auto boundary = &(solver_grid.boundary_vertexes[id_vertex]);

				int id_type = boundary->id_type;
				Point<bool> is_take;
				Point<double> boundary_value = boundary->boundary_value(is_take, solver_grid.boundary_vertexes[id_vertex].GetIdNode(0));
				int global_id = boundary->GetDOFInLocalID(0);

				auto enter_boundary = [global_id, &boundary_value](int position, CSSD_Matrix<Tensor2Rank3D, Point<double>>& global_SLAE) {
					global_SLAE.Diag[global_id].val[position][position] = 1e+30;
				};
				if (is_take.x == true) enter_boundary(0, MassMatrix);
				if (is_take.y == true) enter_boundary(1, MassMatrix);
				if (is_take.z == true) enter_boundary(2, MassMatrix);
			}
		}
		for (int i = 0; i < MassMatrix.Diag.size(); i++)
		{
			if (abs(MassMatrix.Diag[i].val[0][0]) < 1e-12 || abs(MassMatrix.Diag[i].val[1][1]) < 1e-12 || abs(MassMatrix.Diag[i].val[2][2]) < 1e-12)
			{
				MassMatrix.Diag[i].InitializationAsI();
				clear_vertexes.push_back(i);
			}
		}

		//переводим матрицу массы и жесткости в действительные
		math::MakeCopyMatrix_A_into_B(MassMatrix, MassMatrix_doubleSLAE);
		math::MakeCopyMatrix_A_into_B(StiffnessMatrix, StiffnessMatrix_doubleSLAE);
		SLAE_n.Initialization(MassMatrix_doubleSLAE.id_column_for_A_up, MassMatrix_doubleSLAE.id_column_for_A_down);
		Tmp_matrix.Initialization(MassMatrix_doubleSLAE.id_column_for_A_up, MassMatrix_doubleSLAE.id_column_for_A_down);

		double lambda_Mass = MassMatrix_doubleSLAE.SolveMaxEigenvalue();
		double lambda_Stiffness = StiffnessMatrix_doubleSLAE.SolveMaxEigenvalue();
		
		printf_s("\nlambda_M = %.2e\t lambda_K = %.2e\n\n", lambda_Mass, lambda_Stiffness);
	}

	///-------------------
	//Point<double> point_for_out(0, 0, 0.5 - 0.01);

	double TIME_L;
	{
		double v = solver_grid.GetDomain(0)->forMech.GetLongitudinalWaveVelocity_Vp(solver_grid.GetDomain(0)->forThermal.rpho);
		double L = 0.5;
		TIME_L = L / v;

		wchar_t _tmp_wc[1000];
		math::Char_To_Wchar_t(base_result_directory, _tmp_wc, 1000);
		CreateDirectory((LPCTSTR)_tmp_wc, NULL);

		for (int i = 0; i < receivers.size(); i++)
		{
			double len;
			int id_elem_for_out = solver_grid.GetNearestElementID(receivers[i].point, len);

			if (id_elem_for_out >= 0)
			{
				double min_len = 1e+50;
				int target_point = -1;
				auto elem = solver_grid.GetElement(id_elem_for_out);
				for (int j = 0; j < elem->GetNodesCount(); j++)
				{
					double _len = math::SolveLengthVector(receivers[i].point, elem->GetNode(j));
					if (_len < min_len)
					{
						min_len = _len;
						target_point = elem->GetIdNode(j);
					}
				}
				receivers[i].id_node = target_point;
			}
			else
			{
				return;
			}

			FILE* fout_Uz;
			sprintf_s(receivers[i].file_name, "%s/Uz_in_point_%d.txt", base_result_directory, i);
			fopen_s(&fout_Uz, receivers[i].file_name, "w");
			fprintf_s(fout_Uz, "Point: (%.2e, %.2e, %.2e)\n", receivers[i].point.x, receivers[i].point.y, receivers[i].point.z);
			fprintf_s(fout_Uz, "Id vertex: %d\n", receivers[i].id_node);
			fprintf_s(fout_Uz, "T0(Z): %.5e sec\n", TIME_L);
			fclose(fout_Uz);
		}
	}
	///-------------------

	double TIME_curr = TIME_h;
	for (int id_STEP = start_iteration; id_STEP <= end_iteration; id_STEP++)
	{
		TIME_curr += TIME_h;
		printf_s("\n================= Start solution of %d STEP (time = %.2e) ================\n", id_STEP, TIME_curr);
		bool is_print_result = false;
		char result_directory[1000];
		char result_directory_XZ[1000];
		char result_directory_XY[1000];
		char result_directory_YZ[1000];
		char result_directory_save[1000];
		if (id_STEP % step_for_out == 0)
		{
			is_print_result = true;
			sprintf_s(result_directory, sizeof(result_directory), "%s", base_result_directory);
			//sprintf_s(result_directory, sizeof(result_directory), "%s/STEP_%d_t=%.2e", base_result_directory, id_STEP, TIME_curr);
			wchar_t _tmp_wc[1000];
			math::Char_To_Wchar_t(result_directory, _tmp_wc, 1000);
			CreateDirectory((LPCTSTR)_tmp_wc, NULL);


			sprintf_s(result_directory_XY, sizeof(result_directory_XY), "%s/Plane_XY", base_result_directory);
			math::Char_To_Wchar_t(result_directory_XY, _tmp_wc, 1000);
			CreateDirectory((LPCTSTR)_tmp_wc, NULL);
			sprintf_s(result_directory_XZ, sizeof(result_directory_XZ), "%s/Plane_XZ", base_result_directory);
			math::Char_To_Wchar_t(result_directory_XZ, _tmp_wc, 1000);
			CreateDirectory((LPCTSTR)_tmp_wc, NULL);
			sprintf_s(result_directory_YZ, sizeof(result_directory_YZ), "%s/Plane_YZ", base_result_directory);
			math::Char_To_Wchar_t(result_directory_YZ, _tmp_wc, 1000);
			CreateDirectory((LPCTSTR)_tmp_wc, NULL);
			sprintf_s(result_directory_save, sizeof(result_directory_save), "%s/Saves", base_result_directory);
			math::Char_To_Wchar_t(result_directory_save, _tmp_wc, 1000);
			CreateDirectory((LPCTSTR)_tmp_wc, NULL);
		}

		//переопределяем источник и строим вектор в правую часть (через вторые краевые)
		CSSD_Matrix<Tensor2Rank3D, Point<double>> SourseVector;
		SourseVector.F.resize(StiffnessMatrix.F.size());
		if (true) {
			std::function<Point<double>(Point<double>)> time_source_condition1 = [TIME_curr, &gauss_source](Point<double> X) -> Point<double>
			{
				return gauss_source.value(TIME_curr);
			};
			if (gauss_source.is_used)
			{
				second_boundary[gauss_source.target_boundary].value = time_source_condition1;
				for (int i = 0; i < solver_grid.boundary_faces.size(); i++)
				{
					solver_grid.boundary_faces[i].boundary_value = second_boundary[solver_grid.boundary_faces[i].id_type].value;
				}
			}

			std::function<Point<double>(Point<double>)> time_source_condition2 = [TIME_curr, &linearModulationSin_source](Point<double> X) -> Point<double>
			{
				return linearModulationSin_source.value(TIME_curr);
			};
			if (linearModulationSin_source.is_used)
			{
				second_boundary[linearModulationSin_source.target_boundary].value = time_source_condition2;
				for (int i = 0; i < solver_grid.boundary_faces.size(); i++)
				{
					solver_grid.boundary_faces[i].boundary_value = second_boundary[solver_grid.boundary_faces[i].id_type].value;
				}
			}

			for (int id_triangle = 0; id_triangle < solver_grid.boundary_faces.size(); id_triangle++)
			{
				std::vector<Point<double>> local_vector_SLAE;
				solver_grid.boundary_faces[id_triangle].SolveLocalBoundaryVector(local_vector_SLAE, solver_grid.boundary_faces[id_triangle].boundary_value);
				SourseVector.SummPartOfVector(local_vector_SLAE, *solver_grid.boundary_faces[id_triangle].GetElementDOFs());
			};
		}

		auto SolveSLAE = [&SLAE_n, &solver_grid, &clear_vertexes](std::vector<double>& F, std::vector<double>& result, double critical_residual) -> void
		{
			//краевые в правую часть и на диагональ большим числом 
			if (true) {
				//большое число на диагональ
				for (int id_vertex = 0; id_vertex < solver_grid.boundary_vertexes.size(); id_vertex++)
				{
					auto boundary = &(solver_grid.boundary_vertexes[id_vertex]);

					int id_type = boundary->id_type;
					Point<bool> is_take;
					Point<double> boundary_value = boundary->boundary_value(is_take, solver_grid.boundary_vertexes[id_vertex].GetIdNode(0));
					int global_id = boundary->GetDOFInLocalID(0);

					auto enter_boundary = [global_id, &boundary_value](int position, CSSD_Matrix<Tensor2Rank3D, Point<double>>& global_SLAE) {
						global_SLAE.Diag[global_id].val[position][position] = 1e+30;
					};
					if (is_take.x == true) SLAE_n.Diag[global_id * 3 + 0] = 1e+30;
					if (is_take.y == true) SLAE_n.Diag[global_id * 3 + 1] = 1e+30;
					if (is_take.z == true) SLAE_n.Diag[global_id * 3 + 2] = 1e+30;
				}
					
				//добавляем значения в правую часть
				for (int id_vertex = 0; id_vertex < solver_grid.boundary_vertexes.size(); id_vertex++)
				{
					auto boundary = &(solver_grid.boundary_vertexes[id_vertex]);

					int id_type = boundary->id_type;
					Point<bool> is_take;
					Point<double> boundary_value = boundary->boundary_value(is_take, solver_grid.boundary_vertexes[id_vertex].GetIdNode(0));
					int global_id = boundary->GetDOFInLocalID(0);

					if (global_id < SLAE_n.GetMatrixSize() / 3)
					{
						auto enter_value = [&result, &F, &SLAE_n](int value_id, double value) ->void
						{
							result[value_id] = value;
							F[value_id] = value * SLAE_n.Diag[value_id];
							return;
						};

						if (is_take.x) enter_value(global_id * 3 + 0, boundary_value.x);
						if (is_take.y) enter_value(global_id * 3 + 1, boundary_value.y);
						if (is_take.z) enter_value(global_id * 3 + 2, boundary_value.z);
					}

				}
				for (int i = 0; i < clear_vertexes.size(); i++)
				{
					result[clear_vertexes[i]] = 0;
					F[clear_vertexes[i]] = 0;
				}
			}

			printf("Soluting SLAY... (%d)\n", SLAE_n.GetMatrixSize());
			int MaxSize = SLAE_n.GetMatrixSize();
			MaxSize = MaxSize / 10 < 1000 ? 1000 : MaxSize / 10;
			std::vector<double> best_solution;
			math::MakeCopyVector_A_into_B(result, best_solution);
			double current_residual = 1, best_residual = 1e+25;
			int MAX_STEPS = 4;
			int ii = 1;

			CSSD_Matrix<double, double> Precond;
			//Precond.PrecondorSSOR(0.75, SLAE_n);

			double aa = math::MakeInnerProduct(F, F);
			printf_s("//---> %.2e\n", aa);

			for (int i = 0; i <= MAX_STEPS; i++)
			{
				printf_s("//---> I = %d/%d (full size %d) - needed residual %.2e\n", i, MAX_STEPS, SLAE_n.GetMatrixSize(), critical_residual);

				SLAE_n.print_logs = true;
				//current_residual = abs(SLAE_n.MSG_PreconditioningSSOR(F, result, MaxSize, critical_residual, Precond));
				current_residual = abs(SLAE_n.MSG(MaxSize, critical_residual, result, F, true));

				if (current_residual <= critical_residual)
				{
					best_residual = current_residual;
					math::MakeCopyVector_A_into_B(result, best_solution);
					break;
				}
				if (current_residual < best_residual)
				{
					best_residual = current_residual;
					math::MakeCopyVector_A_into_B(result, best_solution);
				}
			}
			math::MakeCopyVector_A_into_B(best_solution, result);

			if (best_residual >= 1)
			{
				printf_s("There are problems in the solution\n");
				int a;
				scanf_s("%d", &a);
			}
		};

		//собираем СЛАУ текущего шага
		{
			SLAE_n.SummMatrix(MassMatrix_doubleSLAE, 1.0, StiffnessMatrix_doubleSLAE, BETTA_Newmark * TIME_h * TIME_h);

			Tmp_matrix.SummMatrix(MassMatrix_doubleSLAE, -2.0, StiffnessMatrix_doubleSLAE, (1.0 / 2.0 - 2.0 * BETTA_Newmark + GAMMA_Newmark) * TIME_h * TIME_h);
			std::vector<double> TMPu_prev;
			Tmp_matrix.MultiplicationMatrixVector(U_prev_d, TMPu_prev);

			Tmp_matrix.SummMatrix(MassMatrix_doubleSLAE, 1.0, StiffnessMatrix_doubleSLAE, (1.0 / 2.0 + BETTA_Newmark - GAMMA_Newmark)* TIME_h* TIME_h);
			std::vector<double> TMPu_prevprev;
			Tmp_matrix.MultiplicationMatrixVector(U_prevprev_d, TMPu_prevprev);

			//double aa = math::MakeInnerProduct(SourseVector.F, SourseVector.F);
			//printf_s("//---> %.2e\n", aa);

			for (int i = 0; i < SourseVector.F.size(); i++)
			{
				Point<double> F = VolumeForceVector.F[i] * TIME_h * TIME_h + SourseVector.F[i] * TIME_h * TIME_h;
				SLAE_n.F[i * 3 + 0] = F.x - TMPu_prev[i * 3 + 0] - TMPu_prevprev[i * 3 + 0];
				SLAE_n.F[i * 3 + 1] = F.y - TMPu_prev[i * 3 + 1] - TMPu_prevprev[i * 3 + 1];
				SLAE_n.F[i * 3 + 2] = F.z - TMPu_prev[i * 3 + 2] - TMPu_prevprev[i * 3 + 2];
			}

			//aa = math::MakeInnerProduct(SLAE_n.F, SLAE_n.F);
			//printf_s("//---> %.2e\n", aa);


		}
		//Решаем СЛАУ
		if (id_STEP == 0)
			math::InitializationVector(SLAE_n.X, 1e-10);
		else
			math::MakeCopyVector_A_into_B(U_prev_d, SLAE_n.X);
		SolveSLAE(SLAE_n.F, SLAE_n.X, 1e-10);
		math::MakeCopyVector_A_into_B(SLAE_n.X, U_curr);

		//вывод сейсмограммы
		if (TIME_curr > TIME_L || true)
		{
			receivers_times.push_back(TIME_curr);
			for (int r = 0; r < receivers.size(); r++)
			{
				receivers[r].U_in_point.push_back(U_curr[receivers[r].id_node]);
			}
		}

		//save reserve solution
		if (id_STEP % (step_for_out*10) == 0)
		{
			FILE* fout_save;
			char name_u_save[5000];

			sprintf_s(name_u_save, "%s/U_s%d_t%.2e.dat", result_directory_save, id_STEP, TIME_curr);
			fopen_s(&fout_save, name_u_save, "w");
			for (int i = 0; i < U_curr.size(); i++)
			{
				fprintf_s(fout_save, "%.16e %.16e %.16e\n", U_curr[i].x, U_curr[i].y, U_curr[i].z);
			}
			fclose(fout_save);

			sprintf_s(name_u_save, "%s/U_s%d_t%.2e.dat", result_directory_save, id_STEP - 1, TIME_curr - TIME_h);
			fopen_s(&fout_save, name_u_save, "w");
			for (int i = 0; i < U_prev.size(); i++)
			{
				fprintf_s(fout_save, "%.16e %.16e %.16e\n", U_prev[i].x, U_prev[i].y, U_prev[i].z);
			}
			fclose(fout_save);
		
			sprintf_s(name_u_save, "%s/U_s%d_t%.2e.dat", result_directory_save, id_STEP - 2, TIME_curr - TIME_h*2);
			fopen_s(&fout_save, name_u_save, "w");
			for (int i = 0; i < U_prevprev.size(); i++)
			{
				fprintf_s(fout_save, "%.16e %.16e %.16e\n", U_prevprev[i].x, U_prevprev[i].y, U_prevprev[i].z);
			}
			fclose(fout_save);
		}

		//output solution
		//without deformations
		if (is_print_result && id_STEP == start_iteration) {
			printf_s("Print the mech result into .dat file... ");

			FILE* fout_tech;
			char name_u_tech[5000];
			sprintf_s(name_u_tech, "%s/U_deformation_%d_t%.2e.dat", result_directory, id_STEP, TIME_curr);
			fopen_s(&fout_tech, name_u_tech, "w");
			char name_in_file[1000];
			sprintf_s(name_in_file, "step%d_time%.2e", id_STEP, TIME_curr);
			std::vector<std::vector<char>> name_value(6);
			char name_v_tmp[6][100];
			sprintf_s(name_v_tmp[0], "sigma_Mises");
			sprintf_s(name_v_tmp[1], "sigma_zz");
			sprintf_s(name_v_tmp[2], "Ux");
			sprintf_s(name_v_tmp[3], "Uy");
			sprintf_s(name_v_tmp[4], "Uz");
			sprintf_s(name_v_tmp[5], "Vz");
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
			double sigma_inv_max = 0;
			int elem_sigma_max = 0;
			for (int i = 0; i < solver_grid.GetElementsCount(); i++)
			{
				auto element = solver_grid.GetElement(i);

				Point<double> Centr = element->GetWeightCentr();

				//Point<double> U = solver_grid.GetSolutionInPoint(i, Centr, U_curr);
				Point<Point<double>> dU = solver_grid.GetDerevativeFromSolutionInPoint(i, Centr, U_curr);
				auto Eps = solver_grid.GetStrainTensorFromSolutionInPoint(i, Centr, dU);
				auto Sigma = solver_grid.GetStressTensorFromSolutionInPoint(i, Centr, Eps);
				auto MisesSigma = solver_grid.GetVonMisesStress(Sigma);

				auto Eps_inv = sqrt((Eps.val[0][0] - Eps.val[1][1]) * (Eps.val[0][0] - Eps.val[1][1])
					+ (Eps.val[1][1] - Eps.val[2][2]) * (Eps.val[1][1] - Eps.val[2][2])
					+ (Eps.val[0][0] - Eps.val[2][2]) * (Eps.val[0][0] - Eps.val[2][2])
					+ 3 * (Eps.val[0][1] * Eps.val[1][0] + Eps.val[0][2] * Eps.val[2][0] + Eps.val[1][2] * Eps.val[2][1]) / 2.0) * sqrt(2.) / 3.;

				value[0][i] = MisesSigma;
				value[1][i] = Sigma.val[2][2];
			}
			value[2].resize(solver_grid.GetDOFsCount());
			value[3].resize(solver_grid.GetDOFsCount());
			value[4].resize(solver_grid.GetDOFsCount());
			value[5].resize(solver_grid.GetDOFsCount());
			for (int i = 0; i < value[3].size(); i++)
			{
				value[2][i] = U_curr[i].x;
				value[3][i] = U_curr[i].y;
				value[4][i] = U_curr[i].z;
				value[5][i] = (U_curr[i].z - U_prev[i].z) / TIME_h;;
			}

			//solver_grid.MoveCoordinates(U_curr);
			solver_grid.printTecPlot3D_DiffDomains(fout_tech, value, name_value, name_in_file, TIME_curr);
			//solver_grid.ReMoveCoordinates(U_curr);
			fclose(fout_tech);
		}
		//via XY plane
		if (is_print_result) {
			math::SimpleGrid grid_XY_plane; //input
			char in_file[1000];
			sprintf_s(in_file, "%s/Plane_XY.dat", mesh_directory);

			FILE* fin_mesh;
			fopen_s(&fin_mesh, in_file, "r");
			if (fin_mesh != NULL)
			{
				fclose(fin_mesh);
				grid_XY_plane.ReadFromSalomeDat(in_file, 2);

				FILE* fout_tech;
				char name_u_tech[5000];
				sprintf_s(name_u_tech, "%s/U_plane_XY_s%d_t%.2e.dat", result_directory_XY, id_STEP, TIME_curr);
				fopen_s(&fout_tech, name_u_tech, "w");
				char name_in_file[1000];
				sprintf_s(name_in_file, "Time_%.4e_XY", TIME_curr);
				std::vector<std::vector<char>> name_value(6);
				char name_v_tmp[6][100];
				sprintf_s(name_v_tmp[0], "sigma_Mises");
				sprintf_s(name_v_tmp[1], "sigma_xx");
				sprintf_s(name_v_tmp[2], "Ux");
				sprintf_s(name_v_tmp[3], "Uy");
				sprintf_s(name_v_tmp[4], "Uz");
				sprintf_s(name_v_tmp[5], "Vx");
				for (int i = 0; i < name_value.size(); i++)
				{
					name_value[i].resize(100);
					for (int j = 0; j < name_value[i].size(); j++)
					{
						name_value[i][j] = name_v_tmp[i][j];
					}
				}
				std::vector<std::vector<double>> value(3 * 2);
				value[0].resize(grid_XY_plane.nvtr.size());
				value[1].resize(grid_XY_plane.nvtr.size());
				value[2].resize(grid_XY_plane.nvtr.size());
				value[3].resize(grid_XY_plane.nvtr.size());
				value[4].resize(grid_XY_plane.nvtr.size());
				value[5].resize(grid_XY_plane.nvtr.size());
				double sigma_inv_max = 0;
				int elem_sigma_max = 0;
				for (int i = 0; i < grid_XY_plane.nvtr.size(); i++)
				{
					Point<double> Centr;
					for (int j = 0; j < grid_XY_plane.nvtr[i].size(); j++)
					{
						Centr += grid_XY_plane.xyz[grid_XY_plane.nvtr[i][j]];
					}
					Centr /= grid_XY_plane.nvtr[i].size();

					double len;
					int id_elem = solver_grid.GetNearestElementID(Centr, len);
					if (id_elem >= 0)
					{
						auto element = solver_grid.GetElement(id_elem);

						Point<double> U = solver_grid.GetSolutionInPoint(id_elem, Centr, U_curr);
						Point<double> U_prev_ = solver_grid.GetSolutionInPoint(id_elem, Centr, U_prev);
						Point<Point<double>> dU = solver_grid.GetDerevativeFromSolutionInPoint(id_elem, Centr, U_curr);
						auto Eps = solver_grid.GetStrainTensorFromSolutionInPoint(id_elem, Centr, dU);
						auto Sigma = solver_grid.GetStressTensorFromSolutionInPoint(id_elem, Centr, Eps);
						auto MisesSigma = solver_grid.GetVonMisesStress(Sigma);

						auto Eps_inv = sqrt((Eps.val[0][0] - Eps.val[1][1]) * (Eps.val[0][0] - Eps.val[1][1])
							+ (Eps.val[1][1] - Eps.val[2][2]) * (Eps.val[1][1] - Eps.val[2][2])
							+ (Eps.val[0][0] - Eps.val[2][2]) * (Eps.val[0][0] - Eps.val[2][2])
							+ 3 * (Eps.val[0][1] * Eps.val[1][0] + Eps.val[0][2] * Eps.val[2][0] + Eps.val[1][2] * Eps.val[2][1]) / 2.0) * sqrt(2.) / 3.;

						value[0][i] = MisesSigma;
						value[1][i] = Sigma.val[0][0];

						value[2][i] = U.x;
						value[3][i] = U.y;
						value[4][i] = U.z;
						value[5][i] = (U_prev_.x - U.x) / TIME_h;
					}
				}

				grid_XY_plane.printTecPlot3D(fout_tech, value, name_value, name_in_file, TIME_curr);
				fclose(fout_tech);
			}
		}
		//via YZ plane
		if (is_print_result) {
			math::SimpleGrid grid_YZ_plane; //input
			char in_file[1000];
			sprintf_s(in_file, "%s/Plane_YZ.dat", mesh_directory);

			FILE* fin_mesh;
			fopen_s(&fin_mesh, in_file, "r");
			if (fin_mesh != NULL)
			{
				fclose(fin_mesh);
				grid_YZ_plane.ReadFromSalomeDat(in_file, 2);

				FILE* fout_tech;
				char name_u_tech[5000];
				sprintf_s(name_u_tech, "%s/U_plane_YZ_s%d_t%.2e.dat", result_directory_YZ, id_STEP, TIME_curr);
				fopen_s(&fout_tech, name_u_tech, "w");
				char name_in_file[1000];
				sprintf_s(name_in_file, "Time_%.4e_YZ", TIME_curr);
				std::vector<std::vector<char>> name_value(6);
				char name_v_tmp[6][100];
				sprintf_s(name_v_tmp[0], "sigma_Mises");
				sprintf_s(name_v_tmp[1], "sigma_zz");
				sprintf_s(name_v_tmp[2], "Ux");
				sprintf_s(name_v_tmp[3], "Uy");
				sprintf_s(name_v_tmp[4], "Uz");
				sprintf_s(name_v_tmp[5], "Vz");
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
						Point<double> U_prev_ = solver_grid.GetSolutionInPoint(id_elem, Centr, U_prev);
						Point<Point<double>> dU = solver_grid.GetDerevativeFromSolutionInPoint(id_elem, Centr, U_curr);
						auto Eps = solver_grid.GetStrainTensorFromSolutionInPoint(id_elem, Centr, dU);
						auto Sigma = solver_grid.GetStressTensorFromSolutionInPoint(id_elem, Centr, Eps);
						auto MisesSigma = solver_grid.GetVonMisesStress(Sigma);

						auto Eps_inv = sqrt((Eps.val[0][0] - Eps.val[1][1]) * (Eps.val[0][0] - Eps.val[1][1])
							+ (Eps.val[1][1] - Eps.val[2][2]) * (Eps.val[1][1] - Eps.val[2][2])
							+ (Eps.val[0][0] - Eps.val[2][2]) * (Eps.val[0][0] - Eps.val[2][2])
							+ 3 * (Eps.val[0][1] * Eps.val[1][0] + Eps.val[0][2] * Eps.val[2][0] + Eps.val[1][2] * Eps.val[2][1]) / 2.0) * sqrt(2.) / 3.;

						value[0][i] = MisesSigma;
						value[1][i] = Sigma.val[2][2];

						value[2][i] = U.x;
						value[3][i] = U.y;
						value[4][i] = U.z;
						value[5][i] = (U_prev_.z - U.z) / TIME_h;
					}
				}

				grid_YZ_plane.printTecPlot3D(fout_tech, value, name_value, name_in_file, TIME_curr);
				fclose(fout_tech);
			}
		}
		//via XZ plane
		if (is_print_result) {
			math::SimpleGrid grid_XZ_plane; //input
			char in_file[1000];
			sprintf_s(in_file, "%s/Plane_XZ.dat", mesh_directory);

			FILE* fin_mesh;
			fopen_s(&fin_mesh, in_file, "r");
			if (fin_mesh != NULL)
			{
				fclose(fin_mesh);
				grid_XZ_plane.ReadFromSalomeDat(in_file, 2);

				FILE* fout_tech;
				char name_u_tech[5000];
				sprintf_s(name_u_tech, "%s/U_plane_XZ_s%d_t%.2e.dat", result_directory_XZ, id_STEP, TIME_curr);
				fopen_s(&fout_tech, name_u_tech, "w");
				char name_in_file[1000];
				sprintf_s(name_in_file, "Time_%.4e_XZ", TIME_curr);
				std::vector<std::vector<char>> name_value(6);
				char name_v_tmp[6][100];
				sprintf_s(name_v_tmp[0], "sigma_Mises");
				sprintf_s(name_v_tmp[1], "sigma_zz");
				sprintf_s(name_v_tmp[2], "Ux");
				sprintf_s(name_v_tmp[3], "Uy");
				sprintf_s(name_v_tmp[4], "Uz");
				sprintf_s(name_v_tmp[5], "Vz");
				for (int i = 0; i < name_value.size(); i++)
				{
					name_value[i].resize(100);
					for (int j = 0; j < name_value[i].size(); j++)
					{
						name_value[i][j] = name_v_tmp[i][j];
					}
				}
				std::vector<std::vector<double>> value(3 * 2);
				value[0].resize(grid_XZ_plane.nvtr.size());
				value[1].resize(grid_XZ_plane.nvtr.size());
				value[2].resize(grid_XZ_plane.nvtr.size());
				value[3].resize(grid_XZ_plane.nvtr.size());
				value[4].resize(grid_XZ_plane.nvtr.size());
				value[5].resize(grid_XZ_plane.nvtr.size());
				double sigma_inv_max = 0;
				int elem_sigma_max = 0;
				for (int i = 0; i < grid_XZ_plane.nvtr.size(); i++)
				{
					Point<double> Centr;
					for (int j = 0; j < grid_XZ_plane.nvtr[i].size(); j++)
					{
						Centr += grid_XZ_plane.xyz[grid_XZ_plane.nvtr[i][j]];
					}
					Centr /= grid_XZ_plane.nvtr[i].size();

					double len;
					int id_elem = solver_grid.GetNearestElementID(Centr, len);
					if (id_elem >= 0)
					{
						auto element = solver_grid.GetElement(id_elem);

						Point<double> U = solver_grid.GetSolutionInPoint(id_elem, Centr, U_curr);
						Point<double> U_prev_ = solver_grid.GetSolutionInPoint(id_elem, Centr, U_prev);
						Point<Point<double>> dU = solver_grid.GetDerevativeFromSolutionInPoint(id_elem, Centr, U_curr);
						auto Eps = solver_grid.GetStrainTensorFromSolutionInPoint(id_elem, Centr, dU);
						auto Sigma = solver_grid.GetStressTensorFromSolutionInPoint(id_elem, Centr, Eps);
						auto MisesSigma = solver_grid.GetVonMisesStress(Sigma);

						auto Eps_inv = sqrt((Eps.val[0][0] - Eps.val[1][1]) * (Eps.val[0][0] - Eps.val[1][1])
							+ (Eps.val[1][1] - Eps.val[2][2]) * (Eps.val[1][1] - Eps.val[2][2])
							+ (Eps.val[0][0] - Eps.val[2][2]) * (Eps.val[0][0] - Eps.val[2][2])
							+ 3 * (Eps.val[0][1] * Eps.val[1][0] + Eps.val[0][2] * Eps.val[2][0] + Eps.val[1][2] * Eps.val[2][1]) / 2.0) * sqrt(2.) / 3.;

						value[0][i] = MisesSigma;
						value[1][i] = Sigma.val[2][2];

						value[2][i] = U.x;
						value[3][i] = U.y;
						value[4][i] = U.z;
						value[5][i] = (U_prev_.z - U.z) / TIME_h;
					}
				}

				grid_XZ_plane.printTecPlot3D(fout_tech, value, name_value, name_in_file, TIME_curr);
				fclose(fout_tech);
			}
		}
		//U in point
		if (is_print_result /*&& TIME_curr > TIME_L*/)
		{
			for (int r = 0; r < receivers.size(); r++)
			{
				FILE* fout_Uz;
				fopen_s(&fout_Uz, receivers[r].file_name, "a");
				for (int i = 0; i < receivers[r].U_in_point.size(); i++)
				{
					fprintf_s(fout_Uz, "%.10e %.10e %.10e %.10e\n", receivers_times[i], receivers[r].U_in_point[i].x, receivers[r].U_in_point[i].y, receivers[r].U_in_point[i].z);
				}
				receivers[r].U_in_point.clear();
				fclose(fout_Uz);
			}
			receivers_times.clear();
		}

		printf_s("\t complite\n");

		math::MakeCopyVector_A_into_B(U_prev, U_prevprev);
		math::MakeCopyVector_A_into_B(U_curr, U_prev);

		math::MakeCopyVector_A_into_B(U_prev_d, U_prevprev_d);
		math::MakeCopyVector_A_into_B(SLAE_n.X, U_prev_d);

		//обновляем сетку
		if (false) {
			solver_grid.MoveCoordinates(U_curr);
			for (int i = 0; i < solver_grid.GetElementsCount(); i++)
			{
				solver_grid.GetElement(i)->SolveAlphaMatrix();
			}
		}
	}
}

void FrostCorrection()
{
	char mesh_name[100];
	sprintf_s(mesh_name, "Mesh4_for_lab1104_time=1.msh");

	math::SimpleGrid geo_grid; //input
	std::vector<std::vector<std::vector<int>>> boundary_faces;
	geo_grid.ReadFromMSH(mesh_name, 1, 4, 4, boundary_faces);
	double minZ = geo_grid.xyz[0].z, maxZ = geo_grid.xyz[0].z;
	for (int i = 0; i < geo_grid.xyz.size(); i++)
	{
		if (geo_grid.xyz[i].z < minZ) minZ = geo_grid.xyz[i].z;
		if (geo_grid.xyz[i].z > maxZ) maxZ = geo_grid.xyz[i].z;
	}

	int target_material = 2;
	int new_material = 3;
	int level[2] = { 0.3, 0.6 }; //по оси Z считая снизу

	for (int i = 0; i < 2; i++)
	{
		double Zmax_curr = minZ + (maxZ - minZ) * level[i];

		std::vector<int> nvkat_new;
	}
}

void main()
{
	//ElasticDeformation_MsFEM_Poly();
	
	//char base_name[1000] = { "D:/Problems/Elastodynamic/252k_tst/param_for_solver" };
	char base_name[1000] = {"./param_for_solver"};
	//char base_name[1000] = { "D:/Problems/Elastodynamic/252k/param_for_solver" };
	//char base_name[1000] = {"D:/Problems/Elastodynamic/29k/param_for_solver"};
	char properties_file[1000];
	FILE* fparam;
	for(int I = 0; I < 1000; I++)
	{
		sprintf_s(properties_file, sizeof(properties_file), "%s_%d.txt", base_name, I);
		fopen_s(&fparam, properties_file, "r");
		if (fparam != NULL)
		{
			fclose(fparam);
			//EffectiveElastisity_MSH(properties_file);
			//Solve_ElasticDeformationProblem_MSH(properties_file);
			//ElasticDeformation_MsFEM_Poly(properties_file);

			//ElastodynamicsProblem_Explicit(properties_file);
			//ElastodynamicsProblem_ExplicitSimple(properties_file);
			//ElastodynamicsProblem_ExplicitRungeKutta4(properties_file);
			//ElastodynamicsProblem_ExplicitPredCorrSimple(properties_file);

			//ElastodynamicsProblem_ExplicitSimple_fast(properties_file);
			ElastodynamicsProblem_ImplicitNewmark_fast(properties_file);
		}
	}

	int a;
	scanf_s("%d", &a);
}