#pragma once
#include "../Library/GeometryGrid.h"
#include "../Library/Point.h"
#include <functional>
#include "../Library/GeometryShape.h"
#include "MsFEM_Grid.h"
#include <stdio.h>
#include <time.h>
#include "omp.h"
#include "../Library/CSSD_Matrix.h"

namespace MsFEM {

	template <typename Dirichlet, typename Neumann>
	void MsFEM_forElasticDeformation_OrderBF2_forPoly(
		bool is_print_logFile, //input
		double MIN_RESIDUAL,
		char* fine_mesh_directory,
		math::SimpleGrid& geo_grid, //input
		std::vector<Dirichlet>& first_boundary,//input
		std::vector<Neumann>& second_boundary,//input
		char* result_directory, //output
		MsFEM::Grid_forMech& solver_grid, //output
		std::vector<Point<double>>& Solution, //output
		CSSD_Matrix<Tensor2Rank3D, Point<double>>& global_SLAE //output
	)
	{
		FILE* stream;
		char name_out_log[1000];
		sprintf_s(name_out_log, sizeof(name_out_log), "%s/log.txt", result_directory);
		/*if (is_print_logFile == false)
		{
			freopen_s(&stream, name_out_log, "w", stdout);
		}
		else {
			fopen_s(&stream, name_out_log, "w");
		}*/

		clock_t t_after = clock();
		double start = omp_get_wtime();

		printf("Initialization of grid...\n");
		solver_grid.is_print_logFile = is_print_logFile;
		solver_grid.Initialization(geo_grid, first_boundary, second_boundary, fine_mesh_directory);
		printf_s("complite                   \n\n");

		/*for (int id_elem = 0; id_elem < solver_grid.GetElementsCount(); id_elem++)
		{
			solver_grid.GetElement(id_elem)->dofs_count = solver_grid.GetElement(id_elem)->GetNodesCount();
		}
		solver_grid.SetDOFsCount(solver_grid.GetVertexCount());*/

		printf("Creation the SLAE portrait...");
		solver_grid.CreationPortrait(global_SLAE);
		printf_s("complite\n\n");

		//second condition
		for (int id_triangle = 0; id_triangle < solver_grid.boundary_faces.size(); id_triangle++)
		{
			std::vector<Point<double>> local_vector_SLAE;
			solver_grid.boundary_faces[id_triangle].SolveLocalBoundaryVector(local_vector_SLAE, solver_grid.boundary_faces[id_triangle].boundary_value);
			global_SLAE.SummPartOfVector(local_vector_SLAE, *solver_grid.boundary_faces[id_triangle].GetElementDOFs());
		};

		//SLAE assembling
		{
			printf("SLAE assembling...\n");
			std::vector<DenseMatrix<Tensor2Rank3D, Point<double>>> local_SLAE(solver_grid.GetElementsCount());
			for (int id_elem = 0; id_elem < solver_grid.GetElementsCount(); id_elem++)
			{
				if (id_elem % 1 == 0)
					printf("Solve element[%d]\r", id_elem);
				auto element = solver_grid.GetElement(id_elem);
								
				FILE* fin_matrix;
				char name_matrix[5000];
				sprintf_s(name_matrix, "%s/local_matrix.txt", element->self_direction);
				fopen_s(&fin_matrix, name_matrix, "r");
				if (fin_matrix == NULL)
				{
					std::function< std::vector<std::vector<double>>(Point<double> X) > koefD = [&](Point<double> X)->std::vector<std::vector<double>>
					{
						return (solver_grid.GetDomain(0))->forMech.GetD(3);
					};
					element->SolveLocalMatrix(local_SLAE[id_elem], koefD);

					fopen_s(&fin_matrix, name_matrix, "w");
					for (int I = 0; I < element->GetDOFsCount(); I++)
					{
						for (int ii = 0; ii < 3; ii++)
						{
							for (int J = 0; J < element->GetDOFsCount(); J++)
							{
								for (int jj = 0; jj < 3; jj++)
								{
									fprintf_s(fin_matrix, "%.8e ", local_SLAE[id_elem].A[I][J].val[ii][jj]);
								}
							}
							fprintf_s(fin_matrix, "\n");
						}
					}
					fclose(fin_matrix);
				}
				else {
					local_SLAE[id_elem].SetSize(element->GetDOFsCount());
					for (int I = 0; I < element->GetDOFsCount(); I++)
					{
						for (int ii = 0; ii < 3; ii++)
						{
							for (int J = 0; J < element->GetDOFsCount(); J++)
							{
								for (int jj = 0; jj < 3; jj++)
								{
									fscanf_s(fin_matrix, "%lf", &local_SLAE[id_elem].A[I][J].val[ii][jj]);
								}
							}
						}
					}
					fclose(fin_matrix);
				}
			}
			for (int id_elem = 0; id_elem < solver_grid.GetElementsCount(); id_elem++)
			{
				if (id_elem % 1000 == 0)
					printf("Add the local matrix of element[%d]\r", id_elem);
				global_SLAE.SummPartOfMatrix(local_SLAE[id_elem], *solver_grid.GetElementDOFs(id_elem));

				//char name[1000];
				//sprintf_s(name, sizeof(name), "%s/Matrix_%d.txt", problem_directory, id_elem);
				//PrintMatrix(name, global_SLAE);
			}
			printf_s("                                                                                    \r");
			printf_s("complite\n\n");
			/*{
				char name[1000];
				sprintf_s(name, sizeof(name), "%s/Matrix_full.txt", problem_directory);
				PrintMatrix(name, global_SLAE);
			}*/
		}

		//Boundary condition
		if(false){
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
					Point<double> boundary_value = boundary->boundary_value(is_take);
					int global_id = boundary->GetDOFInLocalID(0);
					if (global_id < solver_grid.GetDOFsCount())
					{
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
			}

		}

		CSSD_Matrix<double, double> newSLAE;
		math::MakeCopyMatrix_A_into_B(global_SLAE, newSLAE);
		//by double matrix;
		if (true)
		{
			printf("First boundary conditions...");
			for (int id_vertex = 0; id_vertex < solver_grid.boundary_vertexes.size(); id_vertex++)
			{
				auto boundary = &(solver_grid.boundary_vertexes[id_vertex]);

				Point<bool> is_take;
				Point<double> boundary_value = boundary->boundary_value(is_take);
				int global_id = boundary->GetDOFInLocalID(0);

				if (global_id < global_SLAE.GetMatrixSize())
				{
					auto enter_value = [&newSLAE](int value_id, double value) ->void
					{
						newSLAE.X[value_id] = value;
						newSLAE.F[value_id] = value;
						newSLAE.Diag[value_id] = 1;
						//Обнуляем строку
						for (int i = 0; i < newSLAE.A_down[value_id].size(); i++)
						{
							newSLAE.A_down[value_id][i] = 0;
						}
						for (int i = 0; i < newSLAE.A_up[value_id].size(); i++)
						{
							newSLAE.A_up[value_id][i] = 0;
						}
						//обнуляем столбец
						for (int i = 0; i < newSLAE.GetMatrixSize(); i++)
						{
							//верхний треугольник
							if (i < value_id)
							{
								for (int jj = 0; jj < newSLAE.A_up[i].size(); jj++)
								{
									if (newSLAE.id_column_for_A_up[i][jj] == value_id)
									{
										newSLAE.F[i] -= newSLAE.A_up[i][jj] * newSLAE.F[value_id];
										newSLAE.A_up[i][jj] = 0;
									}
								}
							}
							//нижний треугольник
							if (i > value_id)
							{
								for (int jj = 0; jj < newSLAE.A_down[i].size(); jj++)
								{
									if (newSLAE.id_column_for_A_down[i][jj] == value_id)
									{
										newSLAE.F[i] -= newSLAE.A_down[i][jj] * newSLAE.F[value_id];
										newSLAE.A_down[i][jj] = 0;
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
			for (int id = solver_grid.GetVertexCount(); false && id < solver_grid.GetDOFsCount(); id++)
			{
				Point<bool> is_take(true, true, true);
				Point<double> boundary_value(0,0,0);
				int global_id = id;

				if (global_id < global_SLAE.GetMatrixSize())
				{
					auto enter_value = [&newSLAE](int value_id, double value) ->void
					{
						newSLAE.X[value_id] = value;
						newSLAE.F[value_id] = value;
						newSLAE.Diag[value_id] = 1;
						//Обнуляем строку
						for (int i = 0; i < newSLAE.A_down[value_id].size(); i++)
						{
							newSLAE.A_down[value_id][i] = 0;
						}
						for (int i = 0; i < newSLAE.A_up[value_id].size(); i++)
						{
							newSLAE.A_up[value_id][i] = 0;
						}
						//обнуляем столбец
						for (int i = 0; i < newSLAE.GetMatrixSize(); i++)
						{
							//верхний треугольник
							if (i < value_id)
							{
								for (int jj = 0; jj < newSLAE.A_up[i].size(); jj++)
								{
									if (newSLAE.id_column_for_A_up[i][jj] == value_id)
									{
										newSLAE.F[i] -= newSLAE.A_up[i][jj] * newSLAE.F[value_id];
										newSLAE.A_up[i][jj] = 0;
									}
								}
							}
							//нижний треугольник
							if (i > value_id)
							{
								for (int jj = 0; jj < newSLAE.A_down[i].size(); jj++)
								{
									if (newSLAE.id_column_for_A_down[i][jj] == value_id)
									{
										newSLAE.F[i] -= newSLAE.A_down[i][jj] * newSLAE.F[value_id];
										newSLAE.A_down[i][jj] = 0;
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
		}

		clock_t t_before = clock();
		printf("TIME = %lf sec\n", (t_before - t_after) * 1.0 / CLK_TCK);
		double end = omp_get_wtime();
		printf("start = %.16g\nend = %.16g\ndiff = %.16g\n", start, end, end - start);

		math::MakeCopyVector_A_into_B(global_SLAE.X, Solution);
		//BiCG_Stab_double matrix
		if (true) 
		{
			printf("Soluting SLAY... (%d)\n", newSLAE.GetMatrixSize());
			int MaxSize = newSLAE.GetMatrixSize() * 100;
			std::vector<double> best_solution;
			math::MakeCopyVector_A_into_B(newSLAE.X, best_solution);
			double current_residual = 1, best_residual = 1e+25;
			int MAX_STEPS = 15;
			for (int i = 0; i <= MAX_STEPS; i++)
			{
				double needed_residual = pow(10., -1 * (i + 1));
				//current_residual /= 2.;
				printf_s("//---> I = %d/%d (full size %d) - needed residual %.2e\n", i, MAX_STEPS, newSLAE.GetMatrixSize(), needed_residual);
				//current_residual = newSLAE.GMRES(MaxSize, needed_residual);
				//current_residual = abs(newSLAE.BiCG_Stab(MaxSize, needed_residual));
				current_residual = abs(newSLAE.MSG(MaxSize, needed_residual));
				if (current_residual <= MIN_RESIDUAL)
					break;
				if (current_residual > best_residual + (1e-10))
				{
					/*math::MakeCopyVector_A_into_B(best_solution, newSLAE.X);
					printf_s("//---> BEST residual %.2e\n", best_residual);
					break;*/
				}
				if (current_residual < best_residual)
				{
					best_residual = current_residual;
					math::MakeCopyVector_A_into_B(newSLAE.X, best_solution);
				}
			}
			math::MakeCopyVector_A_into_B(best_solution, newSLAE.X);
			printf_s("//---> BEST residual %.2e\n", best_residual);
			math::MakeCopyVector_A_into_B(newSLAE.X, Solution);
		}
		//BiCG_Stab_block matrix
		if (false)
		{
			printf("Soluting SLAY... (%d)\n", global_SLAE.GetMatrixSize());
			int MaxSize = global_SLAE.GetMatrixSize() / 100;
			MaxSize = MaxSize < 10 ? 10 : MaxSize;
			std::vector<Point<double>> best_solution;
			for (int i = 0; i < global_SLAE.X.size(); i++)
			{
				global_SLAE.X[i].x = 1;
				global_SLAE.X[i].y = 1;
				global_SLAE.X[i].z = 1;
			}
			math::MakeCopyVector_A_into_B(global_SLAE.X, best_solution);
			double current_residual = 1, best_residual = 1e+25;
			int MAX_STEPS = 15;
			for (int i = 0; i <= MAX_STEPS; i++)
			{
				current_residual = pow(10., -1 * (i + 1));
				//current_residual /= 2.;
				printf_s("//---> I = %d/%d (full size %d) - needed residual %.2e\n", i, MAX_STEPS, global_SLAE.GetMatrixSize(), current_residual);
				current_residual = abs(global_SLAE.BiCG_Stab_blocks(MaxSize, current_residual));
				if (current_residual <= MIN_RESIDUAL)
					break;
				if (current_residual > best_residual + (1e-10))
				{
					/*math::MakeCopyVector_A_into_B(best_solution, global_SLAE.X);
					printf_s("//---> BEST residual %.2e\n", best_residual);
					break;*/
				}
				if (current_residual < best_residual)
				{
					best_residual = current_residual;
					math::MakeCopyVector_A_into_B(global_SLAE.X, best_solution);
				}
			}
			printf_s("//---> BEST residual %.2e\n", best_residual);
			math::MakeCopyVector_A_into_B(best_solution, global_SLAE.X);
			math::MakeCopyVector_A_into_B(global_SLAE.X, Solution);

			if (best_residual >= 1)
			{
				printf_s("There are problems in the solution\n");
				char a[23];
				scanf_s("%s", &a);
			}
		}
		//Gauss solution
		if (false)
		{
			DenseMatrix_Tensor new_SLAE;
			new_SLAE.SetSize(global_SLAE.GetMatrixSize());
			math::MakeCopyVector_A_into_B(global_SLAE.F, new_SLAE.F);
			for (int i = 0; i < global_SLAE.Diag.size(); i++)
			{
				new_SLAE.A[i][i] = global_SLAE.Diag[i];
			}
			for (int i = 0; i < global_SLAE.A_down.size(); i++)
			{
				for (int jj = 0; jj < global_SLAE.id_column_for_A_down[i].size(); jj++)
				{
					new_SLAE.A[i][global_SLAE.id_column_for_A_down[i][jj]] = global_SLAE.A_down[i][jj];
				}
				for (int jj = 0; jj < global_SLAE.id_column_for_A_up[i].size(); jj++)
				{
					new_SLAE.A[i][global_SLAE.id_column_for_A_up[i][jj]] = global_SLAE.A_up[i][jj];
				}
			}
			new_SLAE.Gauss();
			math::MakeCopyVector_A_into_B(new_SLAE.X, Solution);
		}

		printf("\tcomplit\n\n");
		t_before = clock();
		printf("TIME = %lf sec\n", (t_before - t_after) * 1.0 / CLK_TCK);
		end = omp_get_wtime();
		printf("start = %.16g\nend = %.16g\ndiff = %.16g\n", start, end, end - start);

		//update grid
		/*{
			std::vector<Point<double>> new_x(solver_grid.GetVertexCount());
			for (int id_vertex = 0; id_vertex < solver_grid.GetVertexCount(); id_vertex++)
			{
				new_x[id_vertex] = solver_grid.GetCoordinateViaID(id_vertex) + global_SLAE.X[id_vertex];
			}
			solver_grid.UpdateCoordinates(new_x);
		}*/


		printf_s("complite\n\n");
	};

}