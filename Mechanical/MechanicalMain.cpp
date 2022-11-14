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
#pragma omp parallel for
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
#pragma omp parallel for 
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

			printf_s("Compare solution with FEM... %s", TEST_mesh_dir);
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
				fopen_s(&fin_fem, name_in_fem, "w");
				for (int i = 0; i < FEM_mesh.GetDOFsCount(); i++)
				{
					fscanf_s(fin_fem, "%lf %lf %lf", &FEM_solution[i].x, &FEM_solution[i].y, &FEM_solution[i].z);
				}
				fclose(fin_fem);


				for (int id_macro = 0; id_macro < solver_grid.GetElementsCount(); id_macro++)
				{
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
					solver_grid.GetElement(id_macro)->SetIntegrationLaw(4);
					result += solver_grid.GetElement(id_macro)->SolveIntegral(integral_value);
				}
				result.x = sqrt(result.x);
				result.y = sqrt(result.y);
				result.z = sqrt(result.z);
			}

			FILE* fout_U;
			char name_out_time[1000];
			sprintf_s(name_out_time, "%s/Time_result.txt", base_result_directory);
			fopen_s(&fout_U, name_out_time, "w");
			fprintf_s(fout_U, "TIME = %lf sec\n", (t_before - t_after) * 1.0 / CLK_TCK);
			fprintf_s(fout_U, "start = %.16g\nend = %.16g\ndiff = %.16g\n", start, end, end - start);
			fprintf_s(fout_U, "Error(X) %.5e\n", result.x);
			fprintf_s(fout_U, "Error(Y) %.5e\n", result.y);
			fprintf_s(fout_U, "Error(Z) %.5e\n", result.z);
			fprintf_s(fout_U, "Error(middle) %.5e\n", (result.x + result.y + result.z) / 3.0);
			fclose(fout_U);

			//output solution
			printf_s("Print the mech result into .dat file... ");
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


	char base_name[1000] = {"./param_for_solver"};
	//char base_name[1000] = {"D:/P-ref_MsFEM/Sample_for_convergense/PolyMesh_1Order_new/H1_h1/param_for_solver"};
	char properties_file[1000];
	FILE* fparam;
	for(int I = 0; I < 100; I++)
	{
		sprintf_s(properties_file, sizeof(properties_file), "%s_%d.txt", base_name, I);
		fopen_s(&fparam, properties_file, "r");
		if (fparam != NULL)
		{
			fclose(fparam);
			//EffectiveElastisity_MSH(properties_file);
			ElasticDeformation_MsFEM_Poly(properties_file);
		}
	}
}