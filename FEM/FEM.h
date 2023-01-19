#pragma once
#pragma once
#include "../Library/GeometryGrid.h"
#include "../Library/Point.h"
#include <functional>
#include "../Library/GeometryShape.h"
#include "FEM_Grid.h"
#include <stdio.h>
#include <time.h>
#include "omp.h"
#include "../Library/CSSD_Matrix.h"

namespace FEM {
	void PrintMatrix(char* name, CSSD_Matrix<Tensor2Rank3D, Point<double>> &matrix)
	{
		std::vector<std::vector<double>> dense_matrix;
		std::vector<double> X(matrix.GetMatrixSize() * 3), F(matrix.GetMatrixSize() * 3);
		math::ResizeVector(dense_matrix, matrix.GetMatrixSize() * 3, matrix.GetMatrixSize() * 3);

		for (int id_string = 0; id_string < matrix.GetMatrixSize(); id_string++)
		{
			X[id_string * 3 + 0] = matrix.X[id_string].x;
			X[id_string * 3 + 1] = matrix.X[id_string].y;
			X[id_string * 3 + 2] = matrix.X[id_string].z;

			F[id_string * 3 + 0] = matrix.F[id_string].x;
			F[id_string * 3 + 1] = matrix.F[id_string].y;
			F[id_string * 3 + 2] = matrix.F[id_string].z;

			for (int ii = 0; ii < 3; ii++)
			{
				for (int jj = 0; jj < 3; jj++)
				{
					dense_matrix[id_string * 3 + ii][id_string * 3 + jj] = matrix.Diag[id_string].val[ii][jj];
				}
			}

			for (int j = 0; j < matrix.id_column_for_A_up[id_string].size(); j++)
			{
				int id_row = matrix.id_column_for_A_up[id_string][j];
				for (int ii = 0; ii < 3; ii++)
				{
					for (int jj = 0; jj < 3; jj++)
					{
						dense_matrix[id_string * 3 + ii][id_row * 3 + jj] = matrix.A_up[id_string][j].val[ii][jj];
					}
				}
			}
			for (int j = 0; j < matrix.id_column_for_A_down[id_string].size(); j++)
			{
				int id_row = matrix.id_column_for_A_down[id_string][j];
				for (int ii = 0; ii < 3; ii++)
				{
					for (int jj = 0; jj < 3; jj++)
					{
						dense_matrix[id_string * 3 + ii][id_row * 3 + jj] = matrix.A_down[id_string][j].val[ii][jj];
					}
				}
			}
		}

		FILE *fout;
		fopen_s(&fout, name, "w");
		for (int ii = 0; ii < dense_matrix.size(); ii++)
		{
			fprintf_s(fout, "| ");
			for (int jj = 0; jj < dense_matrix.size(); jj++)
			{
				if (dense_matrix[ii][jj] >= 0) fprintf_s(fout, " ");
				fprintf_s(fout, "%.2e ", dense_matrix[ii][jj]);

				if (dense_matrix[ii][jj] != dense_matrix[jj][ii])
				{
					printf_s("Error");
				}
			}
			fprintf_s(fout, "|\t\t|%8.2e |\t\t|%8.2e |", X[ii], F[ii]);
			fprintf_s(fout, "\n");
		}
		fclose(fout);
	}
	void PrintMatrix(char* name, CSSD_Matrix<double, double> &matrix)
	{
		std::vector<std::vector<double>> dense_matrix;
		std::vector<double> X(matrix.GetMatrixSize()), F(matrix.GetMatrixSize());
		math::ResizeVector(dense_matrix, matrix.GetMatrixSize(), matrix.GetMatrixSize());

		for (int id_string = 0; id_string < matrix.GetMatrixSize(); id_string++)
		{
			X[id_string] = matrix.X[id_string];

			F[id_string] = matrix.F[id_string];

			dense_matrix[id_string][id_string] = matrix.Diag[id_string];

			for (int j = 0; j < matrix.id_column_for_A_up[id_string].size(); j++)
			{
				int id_row = matrix.id_column_for_A_up[id_string][j];

				dense_matrix[id_string][id_row] = matrix.A_up[id_string][j];
			}
			for (int j = 0; j < matrix.id_column_for_A_down[id_string].size(); j++)
			{
				int id_row = matrix.id_column_for_A_down[id_string][j];

				dense_matrix[id_string][id_row] = matrix.A_down[id_string][j];

			}
		}

		FILE *fout;
		fopen_s(&fout, name, "w");
		for (int ii = 0; ii < dense_matrix.size(); ii++)
		{
			fprintf_s(fout, "| ");
			for (int jj = 0; jj < dense_matrix.size(); jj++)
			{
				if (dense_matrix[ii][jj] >= 0) fprintf_s(fout, " ");
				fprintf_s(fout, "%.2e ", dense_matrix[ii][jj]);

				if (dense_matrix[ii][jj] != dense_matrix[jj][ii])
				{
					printf_s("Error");
				}
			}
			fprintf_s(fout, "|\t\t|%8.2e |\t\t|%8.2e |", X[ii], F[ii]);
			fprintf_s(fout, "\n");
		}
		fclose(fout);
	}
	void PrintPortrait(char* name, CSSD_Matrix<Tensor2Rank3D, Point<double>> &matrix)
	{
		std::vector<std::vector<int>> dense_matrix;
		math::ResizeVector(dense_matrix, matrix.GetMatrixSize(), matrix.GetMatrixSize());

		for (int id_string = 0; id_string < matrix.GetMatrixSize(); id_string++)
		{
			dense_matrix[id_string][id_string] = 1;// matrix.Diag[id_string].val[ii][jj];


			for (int j = 0; j < matrix.id_column_for_A_up[id_string].size(); j++)
			{
				int id_row = matrix.id_column_for_A_up[id_string][j];
				dense_matrix[id_string][id_row] = 1;
			}
			for (int j = 0; j < matrix.id_column_for_A_down[id_string].size(); j++)
			{
				int id_row = matrix.id_column_for_A_down[id_string][j];
				dense_matrix[id_string][id_row] = 1;
			}
		}

		FILE *fout;
		fopen_s(&fout, name, "w");
		for (int ii = 0; ii < dense_matrix.size(); ii++)
		{
			fprintf_s(fout, "| ");
			for (int jj = 0; jj < dense_matrix.size(); jj++)
			{
				if (dense_matrix[ii][jj] == 1) fprintf_s(fout, "*");
				else fprintf_s(fout, " ");
			}
			fprintf_s(fout, "\n");
		}
		fclose(fout);
	}
	void PrintPortrait(char* name, CSSD_Matrix<double, double> &matrix)
	{
		std::vector<std::vector<int>> dense_matrix;
		math::ResizeVector(dense_matrix, matrix.GetMatrixSize(), matrix.GetMatrixSize());

		for (int id_string = 0; id_string < matrix.GetMatrixSize(); id_string++)
		{
			dense_matrix[id_string][id_string] = 1;// matrix.Diag[id_string].val[ii][jj];


			for (int j = 0; j < matrix.id_column_for_A_up[id_string].size(); j++)
			{
				int id_row = matrix.id_column_for_A_up[id_string][j];
				dense_matrix[id_string][id_row] = 1;
			}
			for (int j = 0; j < matrix.id_column_for_A_down[id_string].size(); j++)
			{
				int id_row = matrix.id_column_for_A_down[id_string][j];
				dense_matrix[id_string][id_row] = 1;
			}
		}

		FILE *fout;
		fopen_s(&fout, name, "w");
		for (int ii = 0; ii < dense_matrix.size(); ii++)
		{
			fprintf_s(fout, "| ");
			for (int jj = 0; jj < dense_matrix.size(); jj++)
			{
				if (dense_matrix[ii][jj] == 1) fprintf_s(fout, "*");
				else fprintf_s(fout, " ");
			}
			fprintf_s(fout, "\n");
		}
		fclose(fout);
	}

		
	template <typename Dirichlet, typename Neumann>
	void FEM_forBaseEliptic(
		bool is_print_logFile, //input
		double MIN_RESIDUAL,
		math::SimpleGrid &geo_grid, //input
		std::vector<Dirichlet> &first_boundary,//input
		std::vector<Neumann> &second_boundary,//input
		std::function< double(int id_elem, Point<double> X) > &koef_forStiffnessMatrix,
		std::function< double(int id_elem, Point<double> X) > &koef_forMassMatrix,
		std::function< double(int id_elem, Point<double> X) > &koef_forRightSide,
		char* result_directory, //output
		FEM::Grid_forScal &solver_grid, //output
		std::vector<double> &Solution //output
	)
	{
		FILE *stream;
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
		solver_grid.Initialization(geo_grid, first_boundary, second_boundary);
		printf_s("complite                   \n\n");

		CSSD_Matrix<double, double> global_SLAE;

		printf("Creation the SLAE portrait...");
		solver_grid.CreationPortrait(global_SLAE);
		printf_s("complite\n\n");

		if (Solution.size() == global_SLAE.X.size())
		{
			math::MakeCopyVector_A_into_B(Solution, global_SLAE.X);
		}

		//{
		//	char name[1000];
		//	sprintf_s(name, sizeof(name), "%s/Matrix_portrait.txt", result_directory);
		//	PrintPortrait(name, global_SLAE);
		//}

		//second condition
		/*for (int id_triangle = 0; id_triangle < solver_grid.boundary_faces.size(); id_triangle++)
		{
			std::vector<Point<double>> local_vector_SLAE;
			solver_grid.boundary_faces[id_triangle].SolveLocalBoundaryVector(local_vector_SLAE, solver_grid.boundary_faces[id_triangle].boundary_value);
			global_SLAE.SummPartOfVector(local_vector_SLAE, *solver_grid.boundary_faces[id_triangle].GetElementDOFs());
		};*/

		//SLAE assembling
		{
			printf("SLAE assembling...\n");
			std::vector<DenseMatrix<double, double>> local_SLAE(solver_grid.GetElementsCount());
			for (int id_elem = 0; id_elem < solver_grid.GetElementsCount(); id_elem++)
			{
				if (id_elem % 1000 == 0)
					printf("Solve element[%d]\r", id_elem);
				auto element = solver_grid.GetElement(id_elem);
				element->SolveLocalMatrix(local_SLAE[id_elem], koef_forStiffnessMatrix, koef_forMassMatrix, koef_forRightSide);
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
		{
			//first condition
			printf("First boundary conditions...");
			for (int id_vertex = 0; id_vertex < solver_grid.boundary_vertexes.size(); id_vertex++)
			{
				auto boundary = &(solver_grid.boundary_vertexes[id_vertex]);

				int global_id = boundary->GetDOFInLocalID(0);
				Point<double> x = solver_grid.GetCoordinateViaID(global_id);
				double boundary_value = boundary->boundary_value(x);

				global_SLAE.X[global_id] = boundary_value;
				global_SLAE.F[global_id] = boundary_value;
				global_SLAE.Diag[global_id] = 1.0;
				for (int j = 0; j < global_SLAE.A_down[global_id].size(); j++)
				{
					global_SLAE.A_down[global_id][j] = 0.0;
				}
				for (int j = 0; j < global_SLAE.A_up[global_id].size(); j++)
				{
					global_SLAE.A_up[global_id][j] = 0.0;
				}
			}
			printf_s("complite\n\n");
			/*char name[1000];
			sprintf_s(name, sizeof(name), "%s/Matrix_with_boundary.txt", problem_directory);
			PrintMatrix(name, global_SLAE);*/
		}

		clock_t t_before = clock();
		printf("TIME = %lf sec\n", (t_before - t_after) * 1.0 / CLK_TCK);
		double end = omp_get_wtime();
		printf("start = %.16g\nend = %.16g\ndiff = %.16g\n", start, end, end - start);

		//SLAE solution
		//{
		//	printf("Soluting SLAY... (%d)\n", global_SLAE.GetMatrixSize());
		//	int MaxSize = global_SLAE.GetMatrixSize();
		//	int N_matrix = MaxSize / 10000 + 100;
		//	for (int i = 0; i <= N_matrix; i++)
		//	{
		//		printf_s("//---> I = %d/%d (full size %d)\n", i, N_matrix - 1, global_SLAE.GetMatrixSize());
		//		double residual = global_SLAE.BiCG_Stab(MaxSize /*/ N_matrix*/, MIN_RESIDUAL);
		//		if (residual <= MIN_RESIDUAL)
		//			break;
		//	}
		//	math::MakeCopyVector_A_into_B(global_SLAE.X, Solution);
		//}
		{
			

			printf("Soluting SLAY... (%d)\n", global_SLAE.GetMatrixSize());
			int MaxSize = global_SLAE.GetMatrixSize();
			std::vector<double> best_solution;
			math::MakeCopyVector_A_into_B(global_SLAE.X, best_solution);
			double current_residual = 1, best_residual = 1e+25;
			int MAX_STEPS = 100;
			for (int i = 0; i <= MAX_STEPS; i++)
			{
				current_residual = pow(10., -20 * i);
				//current_residual /= 2.;
				printf_s("//---> I = %d/%d (full size %d) - needed residual %.2e\n", i, MAX_STEPS, global_SLAE.GetMatrixSize(), current_residual);
				current_residual = abs(global_SLAE.BiCG_Stab(MaxSize, current_residual));
				if (current_residual <= MIN_RESIDUAL)
					break;
				if (current_residual > best_residual + (1e-10))
				{
					math::MakeCopyVector_A_into_B(best_solution, global_SLAE.X);
					printf_s("//---> BEST residual %.2e\n", best_residual);
					break;
				}
				if (current_residual < best_residual)
				{
					best_residual = current_residual;
					math::MakeCopyVector_A_into_B(global_SLAE.X, best_solution);
				}
			}
			math::MakeCopyVector_A_into_B(global_SLAE.X, Solution);
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
	template <typename Dirichlet, typename Neumann>
	void FEM_forBaseEliptic(
		bool is_print_logFile, //input
		double MIN_RESIDUAL,
		math::SimpleGrid& geo_grid, //input
		std::vector<Dirichlet>& first_boundary,//input
		std::vector<Neumann>& second_boundary,//input
		std::function< double(int id_elem, Point<double> X) >& koef_forStiffnessMatrix, //input
		char* result_directory, //output
		FEM::Grid_forScal& solver_grid, //output
		std::vector<double>& Solution //output
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
		solver_grid.Initialization(geo_grid, first_boundary, second_boundary);
		printf_s("complite                   \n\n");

		CSSD_Matrix<double, double> global_SLAE;

		FILE* fsolution;
		char name_solution[1000];
		sprintf_s(name_solution, "%s/Solution.txt", result_directory);
		fopen_s(&fsolution, name_solution, "r");
		if (fsolution != NULL)
		{
			printf_s("\n============READ current solution===========\n");
			Solution.resize(solver_grid.GetDOFsCount());
			for (int id_dof = 0; id_dof < Solution.size(); id_dof++)
			{
				fscanf_s(fsolution, "%lf", &Solution[id_dof]);
			}
			fclose(fsolution);
		}
		else {

			printf("Creation the SLAE portrait...");
			solver_grid.CreationPortrait(global_SLAE);
			printf_s("complite\n\n");

			if (Solution.size() == global_SLAE.X.size())
			{
				math::MakeCopyVector_A_into_B(Solution, global_SLAE.X);
			}

			//{
			//	char name[1000];
			//	sprintf_s(name, sizeof(name), "%s/Matrix_portrait.txt", result_directory);
			//	PrintPortrait(name, global_SLAE);
			//}

			//second condition
			/*for (int id_triangle = 0; id_triangle < solver_grid.boundary_faces.size(); id_triangle++)
			{
				std::vector<Point<double>> local_vector_SLAE;
				solver_grid.boundary_faces[id_triangle].SolveLocalBoundaryVector(local_vector_SLAE, solver_grid.boundary_faces[id_triangle].boundary_value);
				global_SLAE.SummPartOfVector(local_vector_SLAE, *solver_grid.boundary_faces[id_triangle].GetElementDOFs());
			};*/

			//SLAE assembling
			{
				printf("SLAE assembling...\n");
				std::vector<DenseMatrix<double, double>> local_SLAE(solver_grid.GetElementsCount());
				omp_set_num_threads(math::NUM_THREADS);
#pragma omp parallel for schedule(dynamic)
				for (int id_elem = 0; id_elem < solver_grid.GetElementsCount(); id_elem++)
				{
					if (id_elem % 1000 == 0)
						printf("Solve element[%d]\r", id_elem);
					auto element = solver_grid.GetElement(id_elem);
					element->SolveLocalMatrix(local_SLAE[id_elem], koef_forStiffnessMatrix);
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
			{
				//first condition
				printf("First boundary conditions...\n");
				for (int id_vertex = 0; id_vertex < solver_grid.boundary_vertexes.size(); id_vertex++)
				{
					auto boundary = &(solver_grid.boundary_vertexes[id_vertex]);

					int global_id = boundary->GetDOFInLocalID(0);
					Point<double> x = solver_grid.GetCoordinateViaID(global_id);
					double boundary_value = boundary->boundary_value(x);

					global_SLAE.X[global_id] = boundary_value;
					global_SLAE.F[global_id] = boundary_value;
					global_SLAE.Diag[global_id] = 1.0;

					//обнуляем строку
					for (int j = 0; j < global_SLAE.A_down[global_id].size(); j++)
					{
						global_SLAE.A_down[global_id][j] = 0.0;
					}
					for (int j = 0; j < global_SLAE.A_up[global_id].size(); j++)
					{
						global_SLAE.A_up[global_id][j] = 0.0;
					}
					//обнуляем столбец
					//for (int i = 0; i < global_SLAE.GetMatrixSize(); i++)
					//{
					//	//верхний треугольник
					//	if (i < global_id)
					//	{
					//		for (int jj = 0; jj < global_SLAE.A_up[i].size(); jj++)
					//		{
					//			if (global_SLAE.id_column_for_A_up[i][jj] == global_id)
					//			{
					//				global_SLAE.F[i] -= global_SLAE.A_up[i][jj] * global_SLAE.F[global_id];
					//				global_SLAE.A_up[i][jj] = 0;
					//			}
					//		}
					//	}
					//	//нижний треугольник
					//	if (i > global_id)
					//	{
					//		for (int jj = 0; jj < global_SLAE.A_down[i].size(); jj++)
					//		{
					//			if (global_SLAE.id_column_for_A_down[i][jj] == global_id)
					//			{
					//				global_SLAE.F[i] -= global_SLAE.A_down[i][jj] * global_SLAE.F[global_id];
					//				global_SLAE.A_down[i][jj] = 0;
					//			}
					//		}
					//	}
					//}
				}
				//симметризация
				omp_set_num_threads(math::NUM_THREADS);
#pragma omp parallel for 
				for (int id_row = 0; id_row < global_SLAE.GetMatrixSize(); id_row++)
				{
					if (id_row % 1000 == 0)
					{
						printf("\tsymmetrization: current row %d/%d\r", id_row, global_SLAE.GetMatrixSize());

					}
					int iterator_in_boundary = 0;
					for (int jj = 0; jj < global_SLAE.id_column_for_A_up[id_row].size(); jj++)
					{
						int id_column = global_SLAE.id_column_for_A_up[id_row][jj];
						for (/*iterator_in_boundary = 0*/; iterator_in_boundary < solver_grid.boundary_vertexes.size(); iterator_in_boundary++)
						{
							if (solver_grid.boundary_vertexes[iterator_in_boundary].GetDOFInLocalID(0) == id_column)
							{
								//double boundary_value = solver_grid.boundary_vertexes[iterator_in_boundary].boundary_value(solver_grid.boundary_vertexes[iterator_in_boundary].GetIdNode(0));

								global_SLAE.F[id_row] -= global_SLAE.A_up[id_row][jj] * global_SLAE.F[id_column];
								global_SLAE.A_up[id_row][jj] = 0;

								break;
							}
							if (solver_grid.boundary_vertexes[iterator_in_boundary].GetDOFInLocalID(0) > id_column)
							{
								//iterator_in_boundary--;
								break;
							}
						}
					}

					iterator_in_boundary = 0;
					for (int jj = 0; jj < global_SLAE.id_column_for_A_down[id_row].size(); jj++)
					{
						int id_column = global_SLAE.id_column_for_A_down[id_row][jj];
						for (/*iterator_in_boundary = 0*/; iterator_in_boundary < solver_grid.boundary_vertexes.size(); iterator_in_boundary++)
						{
							if (solver_grid.boundary_vertexes[iterator_in_boundary].GetDOFInLocalID(0) == id_column)
							{
								//Point<double> boundary_value = solver_grid.boundary_vertexes[iterator_in_boundary].boundary_value(solver_grid.boundary_vertexes[iterator_in_boundary].GetIdNode(0));

								global_SLAE.F[id_row] -= global_SLAE.A_down[id_row][jj] * global_SLAE.F[id_column];
								global_SLAE.A_down[id_row][jj] = 0;

								break;
							}
							if (solver_grid.boundary_vertexes[iterator_in_boundary].GetDOFInLocalID(0) > id_column)
							{
								//iterator_in_boundary--;
								break;
							}
						}
					}
				}
				printf_s("complite\n\n");
				/*char name[1000];
				sprintf_s(name, sizeof(name), "%s/Matrix_with_boundary.txt", problem_directory);
				PrintMatrix(name, global_SLAE);*/
			}

			clock_t t_before = clock();
			printf("TIME = %lf sec\n", (t_before - t_after) * 1.0 / CLK_TCK);
			double end = omp_get_wtime();
			printf("start = %.16g\nend = %.16g\ndiff = %.16g\n", start, end, end - start);

			//SLAE solution
			{

				global_SLAE.print_logs = true;
				printf("Soluting SLAY... (%d)\n", global_SLAE.GetMatrixSize());
				printf("Soluting SLAY... (%d)\n", global_SLAE.GetMatrixSize());
				int MaxSize = global_SLAE.GetMatrixSize();
				std::vector<double> best_solution;
				math::MakeCopyVector_A_into_B(global_SLAE.X, best_solution);
				math::MakeCopyVector_A_into_B(global_SLAE.F, global_SLAE.X);
				double current_residual = 1, best_residual = 1e+25;
				int MAX_STEPS = 5;
				int ii = 0;
				CSSD_Matrix<double, double> Predcondor;
				Predcondor.PrecondorSSOR(0.75, global_SLAE);
				for (int i = 0; i <= MAX_STEPS; i++)
				{
					double needed_residual = pow(10., -1 * (ii + 1));
					//current_residual /= 2.;
					printf_s("//---> I = %d/%d (full size %d) - needed residual %.2e\n", i, MAX_STEPS, global_SLAE.GetMatrixSize(), needed_residual);

					if (current_residual > needed_residual)
						current_residual = abs(global_SLAE.MSG_PreconditioningSSOR(MaxSize, needed_residual, Predcondor));

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
					if (current_residual <= MIN_RESIDUAL)
					{
						best_residual = current_residual;
						break;
					}
					if (current_residual < best_residual)
					{
						best_residual = current_residual;
						math::MakeCopyVector_A_into_B(global_SLAE.X, best_solution);

						FILE* fout;
						char name_out[1000];
						sprintf_s(name_out, "%s/Solution_residual_%.2e.txt", result_directory, current_residual);
						fopen_s(&fout, name_out, "w");
						for (int id_dof = 0; id_dof < best_solution.size(); id_dof++)
						{
							fprintf_s(fout, "%.15e\n", best_solution[id_dof]);
						}
						fclose(fout);
					}
				}

				FILE* fout;
				char name_out[1000];
				sprintf_s(name_out, "%s/Solution_residual_%.2e.txt", result_directory, current_residual);
				fopen_s(&fout, name_out, "w");
				for (int id_dof = 0; id_dof < best_solution.size(); id_dof++)
				{
					fprintf_s(fout, "%.15e\n", best_solution[id_dof]);
				}
				fclose(fout);

				math::MakeCopyVector_A_into_B(best_solution, global_SLAE.X);
				math::MakeCopyVector_A_into_B(global_SLAE.X, Solution);
			}
		}
		printf("\tcomplit\n\n");
		clock_t t_before = clock();
		printf("TIME = %lf sec\n", (t_before - t_after) * 1.0 / CLK_TCK);
		clock_t end = omp_get_wtime();
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
	template <typename Dirichlet, typename Neumann>
	void FEM_forBaseEliptic_tensor(
		bool is_print_logFile, //input
		double MIN_RESIDUAL,
		math::SimpleGrid& geo_grid, //input
		std::vector<Dirichlet>& first_boundary,//input
		std::vector<Neumann>& second_boundary,//input
		std::function< Tensor2Rank3D(int id_elem, Point<double> X) >& koef_forStiffnessMatrix, //input
		char* result_directory, //output
		FEM::Grid_forScal& solver_grid, //output
		std::vector<double>& Solution //output
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
		solver_grid.Initialization(geo_grid, first_boundary, second_boundary);
		printf_s("complite                   \n\n");

		CSSD_Matrix<double, double> global_SLAE;

		printf("Creation the SLAE portrait...");
		solver_grid.CreationPortrait(global_SLAE);
		printf_s("complite\n\n");

		if (Solution.size() == global_SLAE.X.size())
		{
			math::MakeCopyVector_A_into_B(Solution, global_SLAE.X);
		}

		//{
		//	char name[1000];
		//	sprintf_s(name, sizeof(name), "%s/Matrix_portrait.txt", result_directory);
		//	PrintPortrait(name, global_SLAE);
		//}

		//second condition
		/*for (int id_triangle = 0; id_triangle < solver_grid.boundary_faces.size(); id_triangle++)
		{
			std::vector<Point<double>> local_vector_SLAE;
			solver_grid.boundary_faces[id_triangle].SolveLocalBoundaryVector(local_vector_SLAE, solver_grid.boundary_faces[id_triangle].boundary_value);
			global_SLAE.SummPartOfVector(local_vector_SLAE, *solver_grid.boundary_faces[id_triangle].GetElementDOFs());
		};*/

		//SLAE assembling
		{
			printf("SLAE assembling...\n");
			std::vector<DenseMatrix<double, double>> local_SLAE(solver_grid.GetElementsCount());
			omp_set_num_threads(math::NUM_THREADS);
#pragma omp parallel for schedule(dynamic)
			for (int id_elem = 0; id_elem < solver_grid.GetElementsCount(); id_elem++)
			{
				if (id_elem % 1000 == 0)
					printf("Solve element[%d]\r", id_elem);
				auto element = solver_grid.GetElement(id_elem);
				element->SolveLocalMatrix(local_SLAE[id_elem], koef_forStiffnessMatrix);
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
		{
			//first condition
			printf("First boundary conditions...\n");
			for (int id_vertex = 0; id_vertex < solver_grid.boundary_vertexes.size(); id_vertex++)
			{
				auto boundary = &(solver_grid.boundary_vertexes[id_vertex]);

				int global_id = boundary->GetDOFInLocalID(0);
				Point<double> x = solver_grid.GetCoordinateViaID(global_id);
				double boundary_value = boundary->boundary_value(x);

				global_SLAE.X[global_id] = boundary_value;
				global_SLAE.F[global_id] = boundary_value;
				global_SLAE.Diag[global_id] = 1.0;

				//обнуляем строку
				for (int j = 0; j < global_SLAE.A_down[global_id].size(); j++)
				{
					global_SLAE.A_down[global_id][j] = 0.0;
				}
				for (int j = 0; j < global_SLAE.A_up[global_id].size(); j++)
				{
					global_SLAE.A_up[global_id][j] = 0.0;
				}
			}
			//симметризация
			omp_set_num_threads(math::NUM_THREADS);
#pragma omp parallel for schedule(dynamic) 
			for (int id_row = 0; id_row < global_SLAE.GetMatrixSize(); id_row++)
			{
				if (id_row % 1000 == 0)
				{
					printf("\tsymmetrization: current row %d/%d\r", id_row, global_SLAE.GetMatrixSize());

				}
				int iterator_in_boundary = 0;
				for (int jj = 0; jj < global_SLAE.id_column_for_A_up[id_row].size(); jj++)
				{
					int id_column = global_SLAE.id_column_for_A_up[id_row][jj];
					for (/*iterator_in_boundary = 0*/; iterator_in_boundary < solver_grid.boundary_vertexes.size(); iterator_in_boundary++)
					{
						if (solver_grid.boundary_vertexes[iterator_in_boundary].GetDOFInLocalID(0) == id_column)
						{
							//double boundary_value = solver_grid.boundary_vertexes[iterator_in_boundary].boundary_value(solver_grid.boundary_vertexes[iterator_in_boundary].GetIdNode(0));

							global_SLAE.F[id_row] -= global_SLAE.A_up[id_row][jj] * global_SLAE.F[id_column];
							global_SLAE.A_up[id_row][jj] = 0;

							break;
						}
						if (solver_grid.boundary_vertexes[iterator_in_boundary].GetDOFInLocalID(0) > id_column)
						{
							//iterator_in_boundary--;
							break;
						}
					}
				}

				iterator_in_boundary = 0;
				for (int jj = 0; jj < global_SLAE.id_column_for_A_down[id_row].size(); jj++)
				{
					int id_column = global_SLAE.id_column_for_A_down[id_row][jj];
					for (/*iterator_in_boundary = 0*/; iterator_in_boundary < solver_grid.boundary_vertexes.size(); iterator_in_boundary++)
					{
						if (solver_grid.boundary_vertexes[iterator_in_boundary].GetDOFInLocalID(0) == id_column)
						{
							//Point<double> boundary_value = solver_grid.boundary_vertexes[iterator_in_boundary].boundary_value(solver_grid.boundary_vertexes[iterator_in_boundary].GetIdNode(0));

							global_SLAE.F[id_row] -= global_SLAE.A_down[id_row][jj] * global_SLAE.F[id_column];
							global_SLAE.A_down[id_row][jj] = 0;

							break;
						}
						if (solver_grid.boundary_vertexes[iterator_in_boundary].GetDOFInLocalID(0) > id_column)
						{
							//iterator_in_boundary--;
							break;
						}
					}
				}
			}
			printf_s("complite\n\n");
			/*char name[1000];
			sprintf_s(name, sizeof(name), "%s/Matrix_with_boundary.txt", problem_directory);
			PrintMatrix(name, global_SLAE);*/
		}

		clock_t t_before = clock();
		printf("TIME = %lf sec\n", (t_before - t_after) * 1.0 / CLK_TCK);
		double end = omp_get_wtime();
		printf("start = %.16g\nend = %.16g\ndiff = %.16g\n", start, end, end - start);

		{
			printf("Soluting SLAY... (%d)\n", global_SLAE.GetMatrixSize());
			printf("Soluting SLAY... (%d)\n", global_SLAE.GetMatrixSize());
			int MaxSize = global_SLAE.GetMatrixSize();
			std::vector<double> best_solution;
			math::MakeCopyVector_A_into_B(global_SLAE.X, best_solution);
			math::MakeCopyVector_A_into_B(global_SLAE.F, global_SLAE.X);
			double current_residual = 1, best_residual = 1e+25;
			int MAX_STEPS = 5;
			int ii = 0;
			CSSD_Matrix<double, double> Predcondor;
			Predcondor.PrecondorSSOR(0.75, global_SLAE);
			for (int i = 0; i <= MAX_STEPS; i++)
			{
				double needed_residual = pow(10., -1 * (ii + 1));
				//current_residual /= 2.;
				printf_s("//---> I = %d/%d (full size %d) - needed residual %.2e\n", i, MAX_STEPS, global_SLAE.GetMatrixSize(), needed_residual);

				if (current_residual > needed_residual)
					current_residual = abs(global_SLAE.MSG_PreconditioningSSOR(MaxSize, needed_residual, Predcondor));

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
				if (current_residual <= MIN_RESIDUAL)
				{
					best_residual = current_residual;
					break;
				}
				if (current_residual < best_residual)
				{
					best_residual = current_residual;
					math::MakeCopyVector_A_into_B(global_SLAE.X, best_solution);
				}
			}
			math::MakeCopyVector_A_into_B(best_solution, global_SLAE.X);
			math::MakeCopyVector_A_into_B(global_SLAE.X, Solution);
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

	template <typename Dirichlet, typename Neumann>
	void FEM_forBaseEliptic_OrderBF2(
		bool is_print_logFile, //input
		double MIN_RESIDUAL,
		math::SimpleGrid& geo_grid, //input
		std::vector<Dirichlet>& first_boundary,//input
		std::vector<Neumann>& second_boundary,//input
		std::function< double(int id_elem, Point<double> X) >& koef_forStiffnessMatrix,
		char* result_directory, //output
		FEM::Grid_forScal_OrderBF2& solver_grid, //output
		std::vector<double>& Solution //output
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
		solver_grid.Initialization(geo_grid, first_boundary, second_boundary);
		printf_s("complite                   \n\n");

		CSSD_Matrix<double, double> global_SLAE;

		printf("Creation the SLAE portrait...");
		solver_grid.CreationPortrait(global_SLAE);
		printf_s("complite\n\n");

		if (Solution.size() == global_SLAE.X.size())
		{
			math::MakeCopyVector_A_into_B(Solution, global_SLAE.X);
		}

		//{
		//	char name[1000];
		//	sprintf_s(name, sizeof(name), "%s/Matrix_portrait.txt", result_directory);
		//	PrintPortrait(name, global_SLAE);
		//}

		//second condition
		/*for (int id_triangle = 0; id_triangle < solver_grid.boundary_faces.size(); id_triangle++)
		{
			std::vector<Point<double>> local_vector_SLAE;
			solver_grid.boundary_faces[id_triangle].SolveLocalBoundaryVector(local_vector_SLAE, solver_grid.boundary_faces[id_triangle].boundary_value);
			global_SLAE.SummPartOfVector(local_vector_SLAE, *solver_grid.boundary_faces[id_triangle].GetElementDOFs());
		};*/

		//SLAE assembling
		{
			printf("SLAE assembling...\n");
			std::vector<DenseMatrix<double, double>> local_SLAE(solver_grid.GetElementsCount());
			omp_set_num_threads(math::NUM_THREADS);
#pragma omp parallel for schedule(dynamic)
			for (int id_elem = 0; id_elem < solver_grid.GetElementsCount(); id_elem++)
			{
				if (id_elem % 1000 == 0)
					printf("Solve element[%d]\r", id_elem);
				auto element = solver_grid.GetElement(id_elem);
				element->SolveLocalMatrix(local_SLAE[id_elem], koef_forStiffnessMatrix);
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


			std::vector<DenseMatrix<double, double>> v_b1;
			std::vector<DenseMatrix<double, double>>(v_b1).swap(local_SLAE);
		}

		//Boundary condition
		{
			//first condition
			printf("First boundary conditions...\n");
			for (int id_vertex = 0; id_vertex < solver_grid.boundary_vertexes.size(); id_vertex++)
			{
				auto boundary = &(solver_grid.boundary_vertexes[id_vertex]);

				int global_id = boundary->GetDOFInLocalID(0);
				Point<double> x = boundary->GetNode(0);
				double boundary_value = boundary->boundary_value(x);

				global_SLAE.X[global_id] = boundary_value;
				global_SLAE.F[global_id] = boundary_value;
				global_SLAE.Diag[global_id] = 1.0;

				//обнуляем строку
				for (int j = 0; j < global_SLAE.A_down[global_id].size(); j++)
				{
					global_SLAE.A_down[global_id][j] = 0.0;
				}
				for (int j = 0; j < global_SLAE.A_up[global_id].size(); j++)
				{
					global_SLAE.A_up[global_id][j] = 0.0;
				}
				//обнуляем столбец
				//for (int i = 0; i < global_SLAE.GetMatrixSize(); i++)
				//{
				//	//верхний треугольник
				//	if (i < global_id)
				//	{
				//		for (int jj = 0; jj < global_SLAE.A_up[i].size(); jj++)
				//		{
				//			if (global_SLAE.id_column_for_A_up[i][jj] == global_id)
				//			{
				//				global_SLAE.F[i] -= global_SLAE.A_up[i][jj] * global_SLAE.F[global_id];
				//				global_SLAE.A_up[i][jj] = 0;
				//			}
				//		}
				//	}
				//	//нижний треугольник
				//	if (i > global_id)
				//	{
				//		for (int jj = 0; jj < global_SLAE.A_down[i].size(); jj++)
				//		{
				//			if (global_SLAE.id_column_for_A_down[i][jj] == global_id)
				//			{
				//				global_SLAE.F[i] -= global_SLAE.A_down[i][jj] * global_SLAE.F[global_id];
				//				global_SLAE.A_down[i][jj] = 0;
				//			}
				//		}
				//	}
				//}
			}
			//симметризация
			omp_set_num_threads(math::NUM_THREADS);
#pragma omp parallel for schedule(dynamic) 
			for (int id_row = 0; id_row < global_SLAE.GetMatrixSize(); id_row++)
			{
				if (id_row % 1000 == 0)
				{
					printf("\tsymmetrization: current row %d/%d\r", id_row, global_SLAE.GetMatrixSize());
				}
				int iterator_in_boundary = 0;
				for (int jj = 0; jj < global_SLAE.id_column_for_A_up[id_row].size(); jj++)
				{
					int id_column = global_SLAE.id_column_for_A_up[id_row][jj];
					for (/*iterator_in_boundary = 0*/; iterator_in_boundary < solver_grid.boundary_vertexes.size(); iterator_in_boundary++)
					{
						if (solver_grid.boundary_vertexes[iterator_in_boundary].GetDOFInLocalID(0) == id_column)
						{
							//double boundary_value = solver_grid.boundary_vertexes[iterator_in_boundary].boundary_value(solver_grid.boundary_vertexes[iterator_in_boundary].GetIdNode(0));

							global_SLAE.F[id_row] -= global_SLAE.A_up[id_row][jj] * global_SLAE.F[id_column];
							global_SLAE.A_up[id_row][jj] = 0;

							break;
						}
						if (solver_grid.boundary_vertexes[iterator_in_boundary].GetDOFInLocalID(0) > id_column)
						{
							//iterator_in_boundary--;
							break;
						}
					}
				}

				iterator_in_boundary = 0;
				for (int jj = 0; jj < global_SLAE.id_column_for_A_down[id_row].size(); jj++)
				{
					int id_column = global_SLAE.id_column_for_A_down[id_row][jj];
					for (/*iterator_in_boundary = 0*/; iterator_in_boundary < solver_grid.boundary_vertexes.size(); iterator_in_boundary++)
					{
						if (solver_grid.boundary_vertexes[iterator_in_boundary].GetDOFInLocalID(0) == id_column)
						{
							//Point<double> boundary_value = solver_grid.boundary_vertexes[iterator_in_boundary].boundary_value(solver_grid.boundary_vertexes[iterator_in_boundary].GetIdNode(0));

							global_SLAE.F[id_row] -= global_SLAE.A_down[id_row][jj] * global_SLAE.F[id_column];
							global_SLAE.A_down[id_row][jj] = 0;

							break;
						}
						if (solver_grid.boundary_vertexes[iterator_in_boundary].GetDOFInLocalID(0) > id_column)
						{
							//iterator_in_boundary--;
							break;
						}
					}
				}
			}
			printf_s("complite\n\n");
			/*char name[1000];
			sprintf_s(name, sizeof(name), "%s/Matrix_with_boundary.txt", problem_directory);
			PrintMatrix(name, global_SLAE);*/
		}

		clock_t t_before = clock();
		printf("TIME = %lf sec\n", (t_before - t_after) * 1.0 / CLK_TCK);
		double end = omp_get_wtime();
		printf("start = %.16g\nend = %.16g\ndiff = %.16g\n", start, end, end - start);

		//SLAE solution
		//{
		//	printf("Soluting SLAY... (%d)\n", global_SLAE.GetMatrixSize());
		//	int MaxSize = global_SLAE.GetMatrixSize();
		//	int N_matrix = MaxSize / 10000 + 100;
		//	for (int i = 0; i <= N_matrix; i++)
		//	{
		//		printf_s("//---> I = %d/%d (full size %d)\n", i, N_matrix - 1, global_SLAE.GetMatrixSize());
		//		double residual = global_SLAE.BiCG_Stab(MaxSize /*/ N_matrix*/, MIN_RESIDUAL);
		//		if (residual <= MIN_RESIDUAL)
		//			break;
		//	}
		//	math::MakeCopyVector_A_into_B(global_SLAE.X, Solution);
		//}
		{


			printf("Soluting SLAY... (%d)\n", global_SLAE.GetMatrixSize());
			int MaxSize = global_SLAE.GetMatrixSize();
			std::vector<double> best_solution;
			math::MakeCopyVector_A_into_B(global_SLAE.X, best_solution);
			math::MakeCopyVector_A_into_B(global_SLAE.F, global_SLAE.X);
			double current_residual = 1, best_residual = 1e+25;
			int MAX_STEPS = 5;
			int ii = 0;
			CSSD_Matrix<double, double> Predcondor; 
			//Predcondor.PrecondorSSOR_summetric(0.75, global_SLAE);
			Predcondor.PrecondorSSOR(0.75, global_SLAE);
			for (int i = 0; i <= MAX_STEPS; i++)
			{
				double needed_residual = pow(10., -1 * (ii + 1));
				//current_residual /= 2.;
				printf_s("//---> I = %d/%d (full size %d) - needed residual %.2e\n", i, MAX_STEPS, global_SLAE.GetMatrixSize(), needed_residual);
				//current_residual = abs(global_SLAE.BCG_Stab2(MaxSize, needed_residual));
				//current_residual = abs(global_SLAE.MSG(MaxSize, needed_residual));

				//for (int i = 0; i < 5 && current_residual > needed_residual; i++)
				//{
				//	current_residual = abs(global_SLAE.BCG_Stab2(MaxSize, needed_residual));

				//	/*current_residual = abs(global_SLAE.BiCG_Stab(MaxSize, needed_residual));
				//	if (current_residual < needed_residual) break;
				//	current_residual = abs(global_SLAE.BiCG_Stab(MaxSize, needed_residual));
				//	if (current_residual < needed_residual) break;
				//	current_residual = abs(global_SLAE.BiCG_Stab(MaxSize, needed_residual));*/
				//	
				//	if (current_residual < needed_residual) break;
				//	//current_residual = abs(global_SLAE.MSG(1000, needed_residual));
				//}
				if(current_residual > needed_residual)
					current_residual = abs(global_SLAE.MSG_PreconditioningSSOR(MaxSize, needed_residual, Predcondor));
				//if(current_residual > needed_residual)
				//	current_residual = abs(global_SLAE.BiCG_Stab(MaxSize, needed_residual));

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
				if (current_residual <= MIN_RESIDUAL)
				{
					best_residual = current_residual;
					break;
				}
				if (current_residual < best_residual)
				{
					best_residual = current_residual;
					math::MakeCopyVector_A_into_B(global_SLAE.X, best_solution);
				}
			}
			math::MakeCopyVector_A_into_B(best_solution, global_SLAE.X);
			math::MakeCopyVector_A_into_B(global_SLAE.X, Solution);
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

	template <typename Dirichlet, typename Neumann>
	void FEM_forElasticDeformation(
		bool is_print_logFile, //input
		double MIN_RESIDUAL,
		math::SimpleGrid &geo_grid, //input
		std::vector<Dirichlet> &first_boundary,//input
		std::vector<Neumann> &second_boundary,//input
		char* result_directory, //output
		FEM::Grid_forMech &solver_grid, //output
		std::vector<Point<double>> &Solution //output
	)
	{
		FILE *stream;
		char name_out_log[1000];
		sprintf_s(name_out_log, sizeof(name_out_log), "%s/log.txt", result_directory);
		/*if (is_print_logFile == false)
		{
			freopen_s(&stream, name_out_log, "w", stdout);
		}
		else {
			fopen_s(&stream, name_out_log, "w");
		}*/
		//freopen_s(&stream, name_out_log, "w", stdout);

		clock_t t_after = clock();
		double start = omp_get_wtime();

		printf("Initialization of grid...\n");
		solver_grid.Initialization(geo_grid, first_boundary, second_boundary);
		printf_s("complite                   \n\n");

		CSSD_Matrix<Tensor2Rank3D, Point<double>> global_SLAE;

		FILE* fsolution;
		char name_solution[1000];
		sprintf_s(name_solution, "%s/Solution.txt", result_directory);
		fopen_s(&fsolution, name_solution, "r");
		if (fsolution != NULL)
		{
			printf_s("\n============READ current solution===========\n");
			Solution.resize(solver_grid.GetDOFsCount());
			for (int id_dof = 0; id_dof < Solution.size(); id_dof++)
			{
				fscanf_s(fsolution, "%lf", &Solution[id_dof]);
			}
			fclose(fsolution);
		}
		else {
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
						/*for (int i = 0; i < element->GetDOFsCount(); i++)
						{
							local_SLAE[id_elem].A[i][i] = 1.0;
						}*/
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
					global_SLAE.SummPartOfMatrix(local_SLAE[id_elem], *solver_grid.GetElementDOFs(id_elem));

					//char name[1000];
					//sprintf_s(name, sizeof(name), "%s/Matrix_%d.txt", problem_directory, id_elem);
					//PrintMatrix(name, global_SLAE);
				}
				/*for (int i = 0; i < global_SLAE.Diag.size(); i++)
				{
					for(int ii = 0; ii < 3; ii++)
					if (math::IsEqual(global_SLAE.Diag[i].val[ii][ii], 0.0))
					{
						global_SLAE.Diag[i].val[ii][ii] = 1.0;
					}
				}*/
				printf_s("                                                                                    \r");
				printf_s("complite\n\n");
				/*{
					char name[1000];
					sprintf_s(name, sizeof(name), "%s/Matrix_full.txt", problem_directory);
					PrintMatrix(name, global_SLAE);
				}*/



			}

			//Boundary condition
			if (false)
			{
				global_SLAE.MakeDivideOnSqrtFormSumm();

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

					//if (is_take.x == true)
					//{
					//	global_SLAE.X[global_id].x = boundary_value.x;
					//	global_SLAE.F[global_id].x = boundary_value.x;
					//	global_SLAE.Diag[global_id].val[0][0] = 1.0;
					//	global_SLAE.Diag[global_id].val[0][1] = 0.0;
					//	global_SLAE.Diag[global_id].val[0][2] = 0.0;
					//	for (int j = 0; j < global_SLAE.A_down[global_id].size(); j++)
					//	{
					//		global_SLAE.A_down[global_id][j].val[0][0] = 0.0;
					//		global_SLAE.A_down[global_id][j].val[0][1] = 0.0;
					//		global_SLAE.A_down[global_id][j].val[0][2] = 0.0;
					//	}
					//	for (int j = 0; j < global_SLAE.A_up[global_id].size(); j++)
					//	{
					//		global_SLAE.A_up[global_id][j].val[0][0] = 0.0;
					//		global_SLAE.A_up[global_id][j].val[0][1] = 0.0;
					//		global_SLAE.A_up[global_id][j].val[0][2] = 0.0;
					//	}
					//	//simmetrisation
					//	int position = 0;
					//	for (int i = 0; i < global_id; i++)
					//	{
					//		for (int j = 0; j < global_SLAE.id_column_for_A_up[i].size() &&
					//			global_SLAE.id_column_for_A_up[i][j] <= global_id; j++)
					//		{
					//			if (global_id == global_SLAE.id_column_for_A_up[i][j])
					//			{
					//				global_SLAE.F[i].x -= global_SLAE.F[global_id].x *global_SLAE.A_up[i][j].val[0][position];
					//				global_SLAE.F[i].y -= global_SLAE.F[global_id].y *global_SLAE.A_up[i][j].val[1][position];
					//				global_SLAE.F[i].z -= global_SLAE.F[global_id].z *global_SLAE.A_up[i][j].val[2][position];
					//				global_SLAE.A_up[i][j].val[0][position] = 0.0;
					//				global_SLAE.A_up[i][j].val[1][position] = 0.0;
					//				global_SLAE.A_up[i][j].val[2][position] = 0.0;
					//			}
					//		}
					//	}
					//	for (int i = global_id+1; i < global_SLAE.id_column_for_A_down.size(); i++)
					//	{
					//		for (int j = 0; j < global_SLAE.id_column_for_A_down[i].size() &&
					//			global_SLAE.id_column_for_A_down[i][j] <= global_id; j++)
					//		{
					//			if (global_id == global_SLAE.id_column_for_A_down[i][j])
					//			{
					//				global_SLAE.F[i].x -= global_SLAE.F[global_id].x *global_SLAE.A_down[i][j].val[0][position];
					//				global_SLAE.F[i].y -= global_SLAE.F[global_id].y *global_SLAE.A_down[i][j].val[1][position];
					//				global_SLAE.F[i].z -= global_SLAE.F[global_id].z *global_SLAE.A_down[i][j].val[2][position];
					//				global_SLAE.A_down[i][j].val[0][position] = 0.0;
					//				global_SLAE.A_down[i][j].val[1][position] = 0.0;
					//				global_SLAE.A_down[i][j].val[2][position] = 0.0;
					//			}
					//		}
					//	}
					//}
					//if (is_take.y == true)
					//{
					//	global_SLAE.X[global_id].y = boundary_value.y;
					//	global_SLAE.F[global_id].y = boundary_value.y;
					//	global_SLAE.Diag[global_id].val[1][0] = 0.0;
					//	global_SLAE.Diag[global_id].val[1][1] = 1.0;
					//	global_SLAE.Diag[global_id].val[1][2] = 0.0;
					//	for (int j = 0; j < global_SLAE.A_down[global_id].size(); j++)
					//	{
					//		global_SLAE.A_down[global_id][j].val[1][0] = 0.0;
					//		global_SLAE.A_down[global_id][j].val[1][1] = 0.0;
					//		global_SLAE.A_down[global_id][j].val[1][2] = 0.0;
					//	}
					//	for (int j = 0; j < global_SLAE.A_up[global_id].size(); j++)
					//	{
					//		global_SLAE.A_up[global_id][j].val[1][0] = 0.0;
					//		global_SLAE.A_up[global_id][j].val[1][1] = 0.0;
					//		global_SLAE.A_up[global_id][j].val[1][2] = 0.0;
					//	}
					//}
					//if (is_take.z == true)
					//{
					//	global_SLAE.X[global_id].z = boundary_value.z;
					//	global_SLAE.F[global_id].z = boundary_value.z;
					//	global_SLAE.Diag[global_id].val[2][0] = 0.0;
					//	global_SLAE.Diag[global_id].val[2][1] = 0.0;
					//	global_SLAE.Diag[global_id].val[2][2] = 1.0;
					//	for (int j = 0; j < global_SLAE.A_down[global_id].size(); j++)
					//	{
					//		global_SLAE.A_down[global_id][j].val[2][0] = 0.0;
					//		global_SLAE.A_down[global_id][j].val[2][1] = 0.0;
					//		global_SLAE.A_down[global_id][j].val[2][2] = 0.0;
					//	}
					//	for (int j = 0; j < global_SLAE.A_up[global_id].size(); j++)
					//	{
					//		global_SLAE.A_up[global_id][j].val[2][0] = 0.0;
					//		global_SLAE.A_up[global_id][j].val[2][1] = 0.0;
					//		global_SLAE.A_up[global_id][j].val[2][2] = 0.0;
					//	}
					//}

					/*printf_s("complite\n\n");
					char name[1000];
					sprintf_s(name, sizeof(name), "%s/Matrix_with_boundary_id_%d.txt", result_directory, global_id);
					PrintMatrix(name, global_SLAE);*/
				}

			}
			CSSD_Matrix<double, double> newSLAE;
			math::MakeCopyMatrix_A_into_B(global_SLAE, newSLAE);
			//by double matrix;
			if (true)
			{
				std::vector<int> use_id;
				printf("First boundary conditions...");
				for (int id_vertex = 0; id_vertex < solver_grid.boundary_vertexes.size(); id_vertex++)
				{
					auto boundary = &(solver_grid.boundary_vertexes[id_vertex]);

					Point<bool> is_take;
					Point<double> boundary_value = boundary->boundary_value(is_take, solver_grid.boundary_vertexes[id_vertex].GetIdNode(0));

					int global_id = boundary->GetDOFInLocalID(0);

					/*boundary_value.x = 0;
					boundary_value.y = 0;
					Point<double> H = solver_grid.GetMaxCoordinate() - solver_grid.GetMinCoordinate();
					boundary_value.z = 0.01 * (solver_grid.GetCoordinateViaID(global_id).z - H.z / 2.0);
					is_take.x = true;
					is_take.y = true;
					is_take.z = true;*/

					if (global_id < global_SLAE.GetMatrixSize())
					{
						auto enter_value = [&newSLAE, &use_id](int value_id, double value) ->void
						{
							use_id.push_back(value_id);
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
							//for (int i = 0; i < newSLAE.GetMatrixSize(); i++)
							//{
							//	//верхний треугольник
							//	if (i < value_id)
							//	{
							//		for (int jj = 0; jj < newSLAE.A_up[i].size(); jj++)
							//		{
							//			if (newSLAE.id_column_for_A_up[i][jj] == value_id)
							//			{
							//				newSLAE.F[i] -= newSLAE.A_up[i][jj] * newSLAE.F[value_id];
							//				newSLAE.A_up[i][jj] = 0;
							//			}
							//		}
							//	}
							//	//нижний треугольник
							//	if (i > value_id)
							//	{
							//		for (int jj = 0; jj < newSLAE.A_down[i].size(); jj++)
							//		{
							//			if (newSLAE.id_column_for_A_down[i][jj] == value_id)
							//			{
							//				newSLAE.F[i] -= newSLAE.A_down[i][jj] * newSLAE.F[value_id];
							//				newSLAE.A_down[i][jj] = 0;
							//			}
							//		}
							//	}
							//}
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

				//
						//симметризация
				omp_set_num_threads(math::NUM_THREADS);
#pragma omp parallel for schedule(dynamic) 
				for (int id_row = 0; id_row < newSLAE.GetMatrixSize(); id_row++)
				{
					if (id_row % 100 == 0)
					{
						printf("\tcurrent %d/%d\r", id_row, newSLAE.GetMatrixSize());
					}
					int iterator_in_boundary = 0;
					for (int jj = 0; jj < newSLAE.id_column_for_A_up[id_row].size(); jj++)
					{
						int id_column = newSLAE.id_column_for_A_up[id_row][jj];
						for (/*iterator_in_boundary = 0*/; iterator_in_boundary < use_id.size(); iterator_in_boundary++)
						{
							if (use_id[iterator_in_boundary] == id_column)
							{
								//double boundary_value = solver_grid.boundary_vertexes[iterator_in_boundary].boundary_value(solver_grid.boundary_vertexes[iterator_in_boundary].GetIdNode(0));

								newSLAE.F[id_row] -= newSLAE.A_up[id_row][jj] * newSLAE.F[id_column];
								newSLAE.A_up[id_row][jj] = 0;

								break;
							}
							if (use_id[iterator_in_boundary] > id_column)
							{
								//iterator_in_boundary--;
								break;
							}
						}
					}

					iterator_in_boundary = 0;
					for (int jj = 0; jj < newSLAE.id_column_for_A_down[id_row].size(); jj++)
					{
						int id_column = newSLAE.id_column_for_A_down[id_row][jj];
						for (/*iterator_in_boundary = 0*/; iterator_in_boundary < use_id.size(); iterator_in_boundary++)
						{
							if (use_id[iterator_in_boundary] == id_column)
							{
								//Point<double> boundary_value = solver_grid.boundary_vertexes[iterator_in_boundary].boundary_value(solver_grid.boundary_vertexes[iterator_in_boundary].GetIdNode(0));

								newSLAE.F[id_row] -= newSLAE.A_down[id_row][jj] * newSLAE.F[id_column];
								newSLAE.A_down[id_row][jj] = 0;

								break;
							}
							if (use_id[iterator_in_boundary] > id_column)
							{
								//iterator_in_boundary--;
								break;
							}
						}
					}
				}
			}

			clock_t t_before = clock();
			printf("TIME = %lf sec\n", (t_before - t_after) * 1.0 / CLK_TCK);
			double end = omp_get_wtime();
			printf("start = %.16g\nend = %.16g\ndiff = %.16g\n", start, end, end - start);
			math::MakeCopyVector_A_into_B(newSLAE.X, Solution);

			if (false)
			{
				printf("Soluting SLAY... (%d)\n", global_SLAE.GetMatrixSize());
				int MaxSize = global_SLAE.GetMatrixSize() / 100;
				std::vector<Point<double>> best_solution;
				math::MakeCopyVector_A_into_B(global_SLAE.X, best_solution);
				double current_residual = 1, best_residual = 1e+25;
				int MAX_STEPS = 15;
				for (int i = 0; i <= MAX_STEPS; i++)
				{
					current_residual = pow(10., -1 * (i + 1));
					//current_residual /= 2.;
					printf_s("//---> I = %d/%d (full size %d) - needed residual %.2e\n", i, MAX_STEPS, global_SLAE.GetMatrixSize(), current_residual);
					current_residual = abs(global_SLAE.BiCG_Stab(MaxSize, current_residual));
					if (current_residual <= MIN_RESIDUAL)
						break;
					if (current_residual > best_residual + (1e-10))
					{
						math::MakeCopyVector_A_into_B(best_solution, global_SLAE.X);
						printf_s("//---> BEST residual %.2e\n", best_residual);
						break;
					}
					if (current_residual < best_residual)
					{
						best_residual = current_residual;
						math::MakeCopyVector_A_into_B(global_SLAE.X, best_solution);
					}
				}
				printf_s("//---> BEST residual %.2e\n", best_residual);
				math::MakeCopyVector_A_into_B(global_SLAE.X, Solution);
			}

			if (true)
			{
				//math::MakeCopyVector_A_into_B(newSLAE.F, newSLAE.X);
				/*for (int i = 0; i < newSLAE.GetMatrixSize(); i++)
					newSLAE.X[i] = 1e-5;*/

					//newSLAE.log_out = stream;
					//newSLAE.print_logs = is_print_logFile;

				printf("Soluting SLAY... (%d)\n", newSLAE.GetMatrixSize());
				int MaxSize = newSLAE.GetMatrixSize();
				MaxSize = MaxSize / 10 < 10 ? 10 : MaxSize / 10;
				std::vector<double> best_solution;
				math::MakeCopyVector_A_into_B(newSLAE.X, best_solution);
				double current_residual = 1, best_residual = 1e+25;
				int MAX_STEPS = 4;
				int ii = 1;
				//newSLAE.MakeDivideOnSqrtFormSumm();
				CSSD_Matrix<double, double> Precond;
				Precond.PrecondorSSOR(0.75, newSLAE);
				std::vector<double> Diag_forSolver(newSLAE.Diag.size());
				math::InitializationVector(Diag_forSolver, 1.0);
				for (int i = 0; i <= MAX_STEPS; i++)
				{
					double needed_residual = pow(10., -1 * (ii + 1));
					//current_residual /= 2.;
					printf_s("//---> I = %d/%d (full size %d) - needed residual %.2e\n", i, MAX_STEPS, newSLAE.GetMatrixSize(), needed_residual);

					newSLAE.print_logs = true;
					current_residual = abs(newSLAE.MSG_PreconditioningSSOR(newSLAE.GetMatrixSize(), needed_residual, Precond));
					//current_residual = abs(newSLAE.BCG_Stab2(newSLAE.GetMatrixSize(), needed_residual));
					//current_residual = abs(newSLAE.BiCG_Diag(newSLAE.GetMatrixSize(), needed_residual, Diag_forSolver));

					if (current_residual < needed_residual)
					{
						i = 0;
						ii++;
					}
					if (current_residual > best_residual * 100)
					{
						//math::MakeCopyVector_A_into_B(best_solution, newSLAE.X);
						//printf_s("//---> BEST residual %.2e\n", best_residual);
						//break;
					}
					if (current_residual <= MIN_RESIDUAL)
					{
						best_residual = current_residual;
						break;
					}
					if (current_residual < best_residual)
					{
						best_residual = current_residual;
						math::MakeCopyVector_A_into_B(newSLAE.X, best_solution);

						FILE* fout;
						char name_out[1000];
						sprintf_s(name_out, "%s/Solution_residual_%.2e.txt", result_directory, current_residual);
						fopen_s(&fout, name_out, "w");
						for (int id_dof = 0; id_dof < best_solution.size(); id_dof++)
						{
							fprintf_s(fout, "%.15e\n", best_solution[id_dof]);
						}
						fclose(fout);
					}
				}
				math::MakeCopyVector_A_into_B(best_solution, newSLAE.X);
				printf_s("//---> BEST residual %.2e\n", best_residual);
				math::MakeCopyVector_A_into_B(newSLAE.X, Solution);

				FILE* fout;
				char name_out[1000];
				sprintf_s(name_out, "%s/Solution_residual_%.2e.txt", result_directory, current_residual);
				fopen_s(&fout, name_out, "w");
				for (int id_dof = 0; id_dof < best_solution.size(); id_dof++)
				{
					fprintf_s(fout, "%.15e\n", best_solution[id_dof]);
				}
				fclose(fout);


				if (best_residual >= 1)
				{
					printf_s("There are problems in the solution\n");
					int a;
					scanf_s("%d", &a);
				}
			}
		}
		printf("\tcomplit\n\n");
		clock_t t_before = clock();
		printf("TIME = %lf sec\n", (t_before - t_after) * 1.0 / CLK_TCK);
		clock_t end = omp_get_wtime();
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

	template <typename Dirichlet, typename Neumann>
	void FEM_forElasticDeformation_renumeration(
		bool is_print_logFile, //input
		double MIN_RESIDUAL,
		math::SimpleGrid& geo_grid, //input
		std::vector<Dirichlet>& first_boundary,//input
		std::vector<Neumann>& second_boundary,//input
		char* result_directory, //output
		FEM::Grid_forMech& solver_grid, //output
		std::vector<Point<double>>& Solution //output
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
		//freopen_s(&stream, name_out_log, "w", stdout);

		clock_t t_after = clock();
		double start = omp_get_wtime();

		//renumeration
		math::SimpleGrid geo_grid_renum;
		std::vector<Dirichlet> first_boundary_renum;
		std::vector<Neumann> second_boundary_renum;
		{
			math::MakeCopyVector_A_into_B(geo_grid.nvkat, geo_grid_renum.nvkat);
			
			std::vector<int> map_renum(geo_grid.xyz.size());
			math::InitializationVector(map_renum, -1);
			int k = 0;

			//ренумерация по краевым
			for (int id_type = 0; id_type < first_boundary.size(); id_type++)
			{
				for (int i = 0; i < first_boundary[id_type].id_vertexes.size(); i++)
				{
					//такую вершину еще не смотрели
					if (map_renum[first_boundary[id_type].id_vertexes[i]] == -1)
					{
						map_renum[first_boundary[id_type].id_vertexes[i]] = k;
						k++;
					}
				}
			}
			/*for(int i = 0; i < map_renum.size(); i++)
				if (map_renum[i] == -1)
				{
					map_renum[i] = k;
					k++;
				}*/

			//ренумерация по материалам
			int material_nums = 0;
			for(int i =0; i < geo_grid.nvkat.size(); i++)
				if (geo_grid.nvkat[i] > material_nums) material_nums = geo_grid.nvkat[i];

			for (int id_material = 0; id_material <= material_nums; id_material++)
			{
				for(int id_elem = 0; id_elem < geo_grid.nvkat.size(); id_elem++)
					if (geo_grid.nvkat[id_elem] == id_material)
					{
						for(int j = 0; j < geo_grid.nvtr[id_elem].size(); j++)
							if (map_renum[geo_grid.nvtr[id_elem][j]] == -1)
							{
								map_renum[geo_grid.nvtr[id_elem][j]] = k;
								k++;
							}
					}
			}

			geo_grid_renum.xyz.resize(geo_grid.xyz.size());
			for (int i = 0; i < map_renum.size(); i++)
			{
				geo_grid_renum.xyz[map_renum[i]] = geo_grid.xyz[i];
			}

			geo_grid_renum.nvtr.resize(geo_grid.nvtr.size());
			for (int i = 0; i < geo_grid.nvtr.size(); i++)
			{
				geo_grid_renum.nvtr[i].resize(geo_grid.nvtr[i].size());
				for (int j = 0; j < geo_grid.nvtr[i].size(); j++)
				{
					geo_grid_renum.nvtr[i][j] = map_renum[geo_grid.nvtr[i][j]];
				}
			}

			first_boundary_renum.resize(first_boundary.size());
			for (int id_type = 0; id_type < first_boundary.size(); id_type++)
			{
				first_boundary_renum[id_type].value = first_boundary[id_type].value;
				first_boundary_renum[id_type].id_vertexes.resize(first_boundary[id_type].id_vertexes.size());
				for (int i = 0; i < first_boundary[id_type].id_vertexes.size(); i++)
				{
					first_boundary_renum[id_type].id_vertexes[i] = map_renum[first_boundary[id_type].id_vertexes[i]];
				}
			}
			second_boundary_renum.resize(second_boundary.size());
			for (int id_type = 0; id_type < second_boundary.size(); id_type++)
			{
				second_boundary_renum[id_type].value = second_boundary[id_type].value;
				math::MakeCopyVector_A_into_B(second_boundary[id_type].values, second_boundary_renum[id_type].values);
				math::MakeCopyVector_A_into_B(second_boundary[id_type].id_base_element, second_boundary_renum[id_type].id_base_element);

				second_boundary_renum[id_type].id_vertexes_as_triangle.resize(second_boundary[id_type].id_vertexes_as_triangle.size());
				for (int i = 0; i < second_boundary[id_type].id_vertexes_as_triangle.size(); i++)
				{
					second_boundary_renum[id_type].id_vertexes_as_triangle[i].resize(second_boundary[id_type].id_vertexes_as_triangle[i].size());
					for(int j = 0; j < second_boundary[id_type].id_vertexes_as_triangle[i].size(); j++)
						second_boundary_renum[id_type].id_vertexes_as_triangle[i][j] = map_renum[second_boundary[id_type].id_vertexes_as_triangle[i][j]];
				}
			}
		} //id_vertexes_as_triangle

		printf("Initialization of grid...\n");
		solver_grid.Initialization(geo_grid_renum, first_boundary_renum, second_boundary_renum);
		printf_s("complite                   \n\n");

		CSSD_Matrix<Tensor2Rank3D, Point<double>> global_SLAE;

		FILE* fsolution;
		char name_solution[1000];
		sprintf_s(name_solution, "%s/Solution.txt", result_directory);
		fopen_s(&fsolution, name_solution, "r");
		if (fsolution != NULL)
		{
			printf_s("\n============READ current solution===========\n");
			Solution.resize(solver_grid.GetDOFsCount());
			for (int id_dof = 0; id_dof < Solution.size(); id_dof++)
			{
				fscanf_s(fsolution, "%lf", &Solution[id_dof]);
			}
			fclose(fsolution);
		}
		else {
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
						/*for (int i = 0; i < element->GetDOFsCount(); i++)
						{
							local_SLAE[id_elem].A[i][i] = 1.0;
						}*/
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
					global_SLAE.SummPartOfMatrix(local_SLAE[id_elem], *solver_grid.GetElementDOFs(id_elem));

					//char name[1000];
					//sprintf_s(name, sizeof(name), "%s/Matrix_%d.txt", problem_directory, id_elem);
					//PrintMatrix(name, global_SLAE);
				}
				/*for (int i = 0; i < global_SLAE.Diag.size(); i++)
				{
					for(int ii = 0; ii < 3; ii++)
					if (math::IsEqual(global_SLAE.Diag[i].val[ii][ii], 0.0))
					{
						global_SLAE.Diag[i].val[ii][ii] = 1.0;
					}
				}*/
				printf_s("                                                                                    \r");
				printf_s("complite\n\n");
				/*{
					char name[1000];
					sprintf_s(name, sizeof(name), "%s/Matrix_full.txt", problem_directory);
					PrintMatrix(name, global_SLAE);
				}*/



			}

			//Boundary condition
			//Big_number
			if (true)
			{
				//first condition
				printf("First boundary conditions...");
				for (int id_vertex = 0; id_vertex < solver_grid.boundary_vertexes.size(); id_vertex++)
				{
					Tensor2Rank3D diag, non_diag;
					Point<double> X, F;
					double big_number = 1e+20;

					auto boundary = &(solver_grid.boundary_vertexes[id_vertex]);

					Point<bool> is_take;
					Point<double> boundary_value = boundary->boundary_value(is_take, boundary->GetIdNode(0));
					///-->>>>>
					/*boundary_value.x = 0;
					boundary_value.y = 0;
					boundary_value.z = 20*(*solver_grid.GetPtrCoordinateViaID(id_vertex)).z - 20*200*0.5;*/
					///-->>>>>
					int global_id = boundary->GetDOFInLocalID(0);
					//printf_s("\n%d\r", global_id);

					if (is_take.x == true)
					{
						global_SLAE.Diag[global_id].val[0][0] = big_number;
						global_SLAE.F[global_id].x = big_number * boundary_value.x;
						global_SLAE.X[global_id].x = boundary_value.x;
					}
					if (is_take.y == true) 
					{
						global_SLAE.Diag[global_id].val[1][1] = big_number;
						global_SLAE.F[global_id].y = big_number * boundary_value.y;
						global_SLAE.X[global_id].y = boundary_value.y;
					}
					if (is_take.z == true)
					{
						global_SLAE.Diag[global_id].val[2][2] = big_number;
						global_SLAE.F[global_id].z = big_number * boundary_value.z;
						global_SLAE.X[global_id].z = boundary_value.z;
					}

				}

			}
			CSSD_Matrix<double, double> newSLAE;
			math::MakeCopyMatrix_A_into_B(global_SLAE, newSLAE);
			//by double matrix;
			if (false)
			{
				std::vector<int> use_id;
				printf("First boundary conditions...");
				for (int id_vertex = 0; id_vertex < solver_grid.boundary_vertexes.size(); id_vertex++)
				{
					auto boundary = &(solver_grid.boundary_vertexes[id_vertex]);

					Point<bool> is_take;
					Point<double> boundary_value = boundary->boundary_value(is_take, solver_grid.boundary_vertexes[id_vertex].GetIdNode(0));

					int global_id = boundary->GetDOFInLocalID(0);

					/*boundary_value.x = 0;
					boundary_value.y = 0;
					Point<double> H = solver_grid.GetMaxCoordinate() - solver_grid.GetMinCoordinate();
					boundary_value.z = 0.01 * (solver_grid.GetCoordinateViaID(global_id).z - H.z / 2.0);*/

					if (global_id < global_SLAE.GetMatrixSize())
					{
						auto enter_value = [&newSLAE, &use_id](int value_id, double value) ->void
						{
							use_id.push_back(value_id);
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
							//for (int i = 0; i < newSLAE.GetMatrixSize(); i++)
							//{
							//	//верхний треугольник
							//	if (i < value_id)
							//	{
							//		for (int jj = 0; jj < newSLAE.A_up[i].size(); jj++)
							//		{
							//			if (newSLAE.id_column_for_A_up[i][jj] == value_id)
							//			{
							//				newSLAE.F[i] -= newSLAE.A_up[i][jj] * newSLAE.F[value_id];
							//				newSLAE.A_up[i][jj] = 0;
							//			}
							//		}
							//	}
							//	//нижний треугольник
							//	if (i > value_id)
							//	{
							//		for (int jj = 0; jj < newSLAE.A_down[i].size(); jj++)
							//		{
							//			if (newSLAE.id_column_for_A_down[i][jj] == value_id)
							//			{
							//				newSLAE.F[i] -= newSLAE.A_down[i][jj] * newSLAE.F[value_id];
							//				newSLAE.A_down[i][jj] = 0;
							//			}
							//		}
							//	}
							//}
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

				//
						//симметризация
				omp_set_num_threads(math::NUM_THREADS);
#pragma omp parallel for schedule(dynamic) 
				for (int id_row = 0; id_row < newSLAE.GetMatrixSize(); id_row++)
				{
					if (id_row % 100 == 0)
					{
						printf("\tcurrent %d/%d\r", id_row, newSLAE.GetMatrixSize());
					}
					int iterator_in_boundary = 0;
					for (int jj = 0; jj < newSLAE.id_column_for_A_up[id_row].size(); jj++)
					{
						int id_column = newSLAE.id_column_for_A_up[id_row][jj];
						for (/*iterator_in_boundary = 0*/; iterator_in_boundary < use_id.size(); iterator_in_boundary++)
						{
							if (use_id[iterator_in_boundary] == id_column)
							{
								//double boundary_value = solver_grid.boundary_vertexes[iterator_in_boundary].boundary_value(solver_grid.boundary_vertexes[iterator_in_boundary].GetIdNode(0));

								newSLAE.F[id_row] -= newSLAE.A_up[id_row][jj] * newSLAE.F[id_column];
								newSLAE.A_up[id_row][jj] = 0;

								break;
							}
							if (use_id[iterator_in_boundary] > id_column)
							{
								//iterator_in_boundary--;
								break;
							}
						}
					}

					iterator_in_boundary = 0;
					for (int jj = 0; jj < newSLAE.id_column_for_A_down[id_row].size(); jj++)
					{
						int id_column = newSLAE.id_column_for_A_down[id_row][jj];
						for (/*iterator_in_boundary = 0*/; iterator_in_boundary < use_id.size(); iterator_in_boundary++)
						{
							if (use_id[iterator_in_boundary] == id_column)
							{
								//Point<double> boundary_value = solver_grid.boundary_vertexes[iterator_in_boundary].boundary_value(solver_grid.boundary_vertexes[iterator_in_boundary].GetIdNode(0));

								newSLAE.F[id_row] -= newSLAE.A_down[id_row][jj] * newSLAE.F[id_column];
								newSLAE.A_down[id_row][jj] = 0;

								break;
							}
							if (use_id[iterator_in_boundary] > id_column)
							{
								//iterator_in_boundary--;
								break;
							}
						}
					}
				}
			}

			clock_t t_before = clock();
			printf("TIME = %lf sec\n", (t_before - t_after) * 1.0 / CLK_TCK);
			double end = omp_get_wtime();
			printf("start = %.16g\nend = %.16g\ndiff = %.16g\n", start, end, end - start);
			math::MakeCopyVector_A_into_B(newSLAE.X, Solution);

			if (false)
			{
				printf("Soluting SLAY... (%d)\n", global_SLAE.GetMatrixSize());
				int MaxSize = global_SLAE.GetMatrixSize() / 100;
				std::vector<Point<double>> best_solution;
				math::MakeCopyVector_A_into_B(global_SLAE.X, best_solution);
				double current_residual = 1, best_residual = 1e+25;
				int MAX_STEPS = 15;
				for (int i = 0; i <= MAX_STEPS; i++)
				{
					current_residual = pow(10., -1 * (i + 1));
					//current_residual /= 2.;
					printf_s("//---> I = %d/%d (full size %d) - needed residual %.2e\n", i, MAX_STEPS, global_SLAE.GetMatrixSize(), current_residual);
					current_residual = abs(global_SLAE.BiCG_Stab(MaxSize, current_residual));
					if (current_residual <= MIN_RESIDUAL)
						break;
					if (current_residual > best_residual + (1e-10))
					{
						math::MakeCopyVector_A_into_B(best_solution, global_SLAE.X);
						printf_s("//---> BEST residual %.2e\n", best_residual);
						break;
					}
					if (current_residual < best_residual)
					{
						best_residual = current_residual;
						math::MakeCopyVector_A_into_B(global_SLAE.X, best_solution);
					}
				}
				printf_s("//---> BEST residual %.2e\n", best_residual);
				math::MakeCopyVector_A_into_B(global_SLAE.X, Solution);
			}

			if (true)
			{
				//math::MakeCopyVector_A_into_B(newSLAE.F, newSLAE.X);
				for (int i = 0; i < newSLAE.GetMatrixSize(); i++)
					newSLAE.X[i] = 1e-5;

					//newSLAE.log_out = stream;
					//newSLAE.print_logs = is_print_logFile;

				printf("Soluting SLAY... (%d)\n", newSLAE.GetMatrixSize());
				int MaxSize = newSLAE.GetMatrixSize();
				MaxSize = MaxSize / 10 < 10 ? 10 : MaxSize / 10;
				std::vector<double> best_solution;
				math::MakeCopyVector_A_into_B(newSLAE.X, best_solution);
				double current_residual = 1, best_residual = 1e+25;
				int MAX_STEPS = 4;
				int ii = 1;
				//newSLAE.MakeDivideOnSqrtFormSumm();
				CSSD_Matrix<double, double> Precond;
				Precond.PrecondorSSOR(0.75, newSLAE);
				std::vector<double> _diag(newSLAE.Diag.size());
				math::InitializationVector(_diag, 1.0);
				for (int i = 0; i <= MAX_STEPS; i++)
				{
					double needed_residual = pow(10., -1 * (ii + 1));
					//current_residual /= 2.;
					printf_s("//---> I = %d/%d (full size %d) - needed residual %.2e\n", i, MAX_STEPS, newSLAE.GetMatrixSize(), needed_residual);

					newSLAE.print_logs = true;
					//current_residual = abs(newSLAE.MSG_PreconditioningSSOR(newSLAE.GetMatrixSize(), needed_residual, Precond));
					//current_residual = abs(newSLAE.BCG_Stab2(newSLAE.GetMatrixSize(), needed_residual));
					current_residual = abs(newSLAE.BiCG_Diag(newSLAE.GetMatrixSize(), needed_residual, _diag));

					if (current_residual < needed_residual)
					{
						i = 0;
						ii++;
					}
					if (current_residual > best_residual * 100)
					{
						//math::MakeCopyVector_A_into_B(best_solution, newSLAE.X);
						//printf_s("//---> BEST residual %.2e\n", best_residual);
						//break;
					}
					if (current_residual <= MIN_RESIDUAL)
					{
						best_residual = current_residual;
						break;
					}
					if (current_residual < best_residual)
					{
						best_residual = current_residual;
						math::MakeCopyVector_A_into_B(newSLAE.X, best_solution);

						FILE* fout;
						char name_out[1000];
						sprintf_s(name_out, "%s/Solution_residual_%.2e.txt", result_directory, current_residual);
						fopen_s(&fout, name_out, "w");
						for (int id_dof = 0; id_dof < best_solution.size(); id_dof++)
						{
							fprintf_s(fout, "%.15e\n", best_solution[id_dof]);
						}
						fclose(fout);
					}
				}
				math::MakeCopyVector_A_into_B(best_solution, newSLAE.X);
				printf_s("//---> BEST residual %.2e\n", best_residual);
				math::MakeCopyVector_A_into_B(newSLAE.X, Solution);

				FILE* fout;
				char name_out[1000];
				sprintf_s(name_out, "%s/Solution_residual_%.2e.txt", result_directory, current_residual);
				fopen_s(&fout, name_out, "w");
				for (int id_dof = 0; id_dof < best_solution.size(); id_dof++)
				{
					fprintf_s(fout, "%.15e\n", best_solution[id_dof]);
				}
				fclose(fout);


				if (best_residual >= 1)
				{
					printf_s("There are problems in the solution\n");
					int a;
					scanf_s("%d", &a);
				}
			}
		}
		printf("\tcomplit\n\n");
		clock_t t_before = clock();
		printf("TIME = %lf sec\n", (t_before - t_after) * 1.0 / CLK_TCK);
		clock_t end = omp_get_wtime();
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

	template <typename Dirichlet, typename Neumann>
	void FEM_forElasticDeformation(
		bool is_print_logFile, //input
		double MIN_RESIDUAL,
		math::SimpleGrid& geo_grid, //input
		std::vector<Dirichlet>& first_boundary,//input
		std::vector<Neumann>& second_boundary,//input
		std::function< std::vector<std::vector<double>>(int id_elem, Point<double> X) >& koef_forStiffnessMatrix,
		std::function< double(int id_elem, Point<double> X) >& koef_forMassMatrix,
		std::function< Point<double>(int id_elem, Point<double> X) >& koef_forRightSide,
		char* result_directory, //output
		FEM::Grid_forMech& solver_grid, //output
		std::vector<Point<double>>& Solution //output
	)
	{
		FILE* stream;
		char name_out_log[1000];
		sprintf_s(name_out_log, sizeof(name_out_log), "%s/log.txt", result_directory);
		if (is_print_logFile != false)
		{
			freopen_s(&stream, name_out_log, "w", stdout);
		}
		else {
			fopen_s(&stream, name_out_log, "w");
		}

		clock_t t_after = clock();
		double start = omp_get_wtime();

		printf("Initialization of grid...\n");
		solver_grid.Initialization(geo_grid, first_boundary, second_boundary);
		printf_s("complite                   \n\n");

		CSSD_Matrix<Tensor2Rank3D, Point<double>> global_SLAE;

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
				if (id_elem % 1000 == 0)
					printf("Solve element[%d]\r", id_elem);
				auto element = solver_grid.GetElement(id_elem);

				std::function<std::vector<std::vector<double>>(Point<double>)> D = [&](Point<double> X) {return koef_forStiffnessMatrix(id_elem, X); };
				std::function<double(Point<double>)> M = [&](Point<double> X) {return koef_forMassMatrix(id_elem, X); };
				std::function<Point<double>(Point<double>)> F = [&](Point<double> X) {return koef_forRightSide(id_elem, X); };

				element->SolveLocalMatrix(local_SLAE[id_elem], D, M, F);
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
			/*for (int i = 0; i < global_SLAE.Diag.size(); i++)
			{
				for(int ii = 0; ii < 3; ii++)
				if (math::IsEqual(global_SLAE.Diag[i].val[ii][ii], 0.0))
				{
					global_SLAE.Diag[i].val[ii][ii] = 1.0;
				}
			}*/
			printf_s("                                                                                    \r");
			printf_s("complite\n\n");
			/*{
				char name[1000];
				sprintf_s(name, sizeof(name), "%s/Matrix_full.txt", problem_directory);
				PrintMatrix(name, global_SLAE);
			}*/



		}

		//Boundary condition
		if (false)
		{
			global_SLAE.MakeDivideOnSqrtFormSumm();

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
				Point<double> boundary_value = boundary->boundary_value(is_take, solver_grid.boundary_vertexes[id_vertex].GetIdNode(0));
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
		}

		clock_t t_before = clock();
		printf("TIME = %lf sec\n", (t_before - t_after) * 1.0 / CLK_TCK);
		double end = omp_get_wtime();
		printf("start = %.16g\nend = %.16g\ndiff = %.16g\n", start, end, end - start);

		if (false)
		{
			printf("Soluting SLAY... (%d)\n", global_SLAE.GetMatrixSize());
			int MaxSize = global_SLAE.GetMatrixSize() / 100;
			std::vector<Point<double>> best_solution;
			math::MakeCopyVector_A_into_B(global_SLAE.X, best_solution);
			double current_residual = 1, best_residual = 1e+25;
			int MAX_STEPS = 15;
			for (int i = 0; i <= MAX_STEPS; i++)
			{
				current_residual = pow(10., -1 * (i + 1));
				//current_residual /= 2.;
				printf_s("//---> I = %d/%d (full size %d) - needed residual %.2e\n", i, MAX_STEPS, global_SLAE.GetMatrixSize(), current_residual);
				current_residual = abs(global_SLAE.BiCG_Stab(MaxSize, current_residual));
				if (current_residual <= MIN_RESIDUAL)
					break;
				if (current_residual > best_residual + (1e-10))
				{
					math::MakeCopyVector_A_into_B(best_solution, global_SLAE.X);
					printf_s("//---> BEST residual %.2e\n", best_residual);
					break;
				}
				if (current_residual < best_residual)
				{
					best_residual = current_residual;
					math::MakeCopyVector_A_into_B(global_SLAE.X, best_solution);
				}
			}
			printf_s("//---> BEST residual %.2e\n", best_residual);
			math::MakeCopyVector_A_into_B(global_SLAE.X, Solution);
		}
		if (true)
		{
			for (int i = 0; i < newSLAE.GetMatrixSize() / 3; i++)
			{
				newSLAE.X[i * 3 + 0] = Solution[i].x;
				newSLAE.X[i * 3 + 1] = Solution[i].y;
				newSLAE.X[i * 3 + 2] = Solution[i].z;
			}

			newSLAE.log_out = stream;
			newSLAE.print_logs = is_print_logFile;

			printf("Soluting SLAY... (%d)\n", newSLAE.GetMatrixSize());
			int MaxSize = newSLAE.GetMatrixSize();
			MaxSize = MaxSize / 10 < 10 ? 10 : MaxSize / 10;
			std::vector<double> best_solution;
			math::MakeCopyVector_A_into_B(newSLAE.X, best_solution);
			double current_residual = 1, best_residual = 1e+25;
			int MAX_STEPS = 15;
			int ii = 1;
			//newSLAE.MakeDivideOnSqrtFormSumm();
			CSSD_Matrix<double, double> Predcondor;
			Predcondor.PrecondorSSOR_summetric(0.75, newSLAE);
			for (int i = 0; i <= MAX_STEPS; i++)
			{
				double needed_residual = pow(10., -1 * (ii + 1));
				//current_residual /= 2.;
				printf_s("//---> I = %d/%d (full size %d) - needed residual %.2e\n", i, MAX_STEPS, newSLAE.GetMatrixSize(), needed_residual);
				//current_residual = newSLAE.GMRES(5, needed_residual);

				//current_residual = abs(newSLAE.BiCG_Stab(MaxSize, needed_residual));
				//current_residual = abs(newSLAE.MSG(MaxSize, needed_residual));
				current_residual = abs(newSLAE.MSG_PreconditioningSSOR(MaxSize, needed_residual, Predcondor));
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
				if (current_residual <= MIN_RESIDUAL)
				{
					best_residual = current_residual;
					break;
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

			if (best_residual >= 1)
			{
				printf_s("There are problems in the solution\n");
				int a;
				scanf_s("%d", &a);
			}
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

	template <typename Dirichlet, typename Neumann>
	void FEM_3D_forElasticDeformation(
		bool is_print_logFile, //input
		double MIN_RESIDUAL,
		math::SimpleGrid& geo_grid, //input
		std::vector<Dirichlet>& first_boundary,//input
		std::vector<Neumann>& second_boundary,//input
		std::function<Point<double>(bool&, int, Point<double>)> &sourse,
		char* result_directory, //output
		FEM::Grid_forMech& solver_grid, //output
		std::vector<Point<double>>& Solution, //output
		CSSD_Matrix<Tensor2Rank3D, Point<double>>& Stiffness_matrix //output
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

		if (solver_grid.GetElementsCount() == 0)
		{
			printf("Initialization of grid...\n");
			solver_grid.Initialization(geo_grid, first_boundary, second_boundary);
			printf_s("complite                   \n\n");
		}
		else {
			//update value functions for 1st boundary

			int id = 0;
			for (int id_type = 0; id_type < first_boundary.size(); id_type++)
			{
				for (int id_vertex = 0; id_vertex < first_boundary[id_type].id_vertexes.size(); id_vertex++)
				{
					solver_grid.boundary_vertexes[id].boundary_value = first_boundary[id_type].value;
					id++;
				}
			}
		}

		CSSD_Matrix<Tensor2Rank3D, Point<double>> global_SLAE;

		//SLAE assembling
		if(Stiffness_matrix.GetMatrixSize() == 0)
		{
			printf("Creation the SLAE portrait...");
			solver_grid.CreationPortrait(global_SLAE);
			printf_s("complite\n\n");

			printf("SLAE assembling...\n");
			std::vector<DenseMatrix<Tensor2Rank3D, Point<double>>> local_SLAE(solver_grid.GetElementsCount());
			for (int id_elem = 0; id_elem < solver_grid.GetElementsCount(); id_elem++)
			{
				if (id_elem % 1000 == 0)
					printf("Solve element[%d]\r", id_elem);
				auto element = solver_grid.GetElement(id_elem);
				std::function< std::vector<std::vector<double>>(Point<double> X) > koefD = [&](Point<double> X)->std::vector<std::vector<double>>
				{
					return (solver_grid.GetDomain(solver_grid.GetElement(id_elem)->GetIdDomain()))->forMech.GetD(3);
				};
				element->SolveLocalMatrix(local_SLAE[id_elem], koefD);

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

			Stiffness_matrix.SetMatrix(global_SLAE);
		}
		else {
			global_SLAE.SetMatrix(Stiffness_matrix);
		}

		printf("Right side...\n");
		std::vector< std::vector<Point<double>>> local_right_side(solver_grid.GetElementsCount());
		for (int id_elem = 0; id_elem < solver_grid.GetElementsCount(); id_elem++)
		{
			if (id_elem % 1000 == 0)
				printf("Solve element[%d]\r", id_elem);
			auto element = solver_grid.GetElement(id_elem);
			
			Point<double> result;
			bool is_solve;
			result = sourse(is_solve, id_elem, Point<double>(0,0,0));
			if (is_solve == true)
			{
				std::function< Point<double>(Point<double> X) > f = [&](Point<double> X)->Point<double>
				{
					Point<double> result;
					result = sourse(is_solve, id_elem, X);
					return result;
				};
				element->SolveRightSide(local_right_side[id_elem], f);
			}

		}
		for (int id_elem = 0; id_elem < solver_grid.GetElementsCount(); id_elem++)
		{
			if (id_elem % 1000 == 0)
				printf("Add the local matrix of element[%d]\r", id_elem);
			global_SLAE.SummPartOfVector(local_right_side[id_elem], *solver_grid.GetElementDOFs(id_elem));

			//char name[1000];
			//sprintf_s(name, sizeof(name), "%s/Matrix_%d.txt", problem_directory, id_elem);
			//PrintMatrix(name, global_SLAE);
		}
		printf_s("                                                                                    \r");
		printf_s("complite\n\n");

		//second condition
		for (int id_triangle = 0; id_triangle < solver_grid.boundary_faces.size(); id_triangle++)
		{
			std::vector<Point<double>> local_vector_SLAE;
			solver_grid.boundary_faces[id_triangle].SolveLocalBoundaryVector(local_vector_SLAE, solver_grid.boundary_faces[id_triangle].boundary_value);
			global_SLAE.SummPartOfVector(local_vector_SLAE, *solver_grid.boundary_faces[id_triangle].GetElementDOFs());
		};

		for (int i = 0; i < global_SLAE.X.size(); i++)
		{
			global_SLAE.X[i] = 1.0;
		}

		//Boundary condition
		if(false)
		{
			global_SLAE.MakeDivideOnSqrtFormSumm();

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
				Point<double> boundary_value = boundary->boundary_value(is_take, solver_grid.boundary_vertexes[id_vertex].GetIdNode(0));
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

				//if (is_take.x == true)
				//{
				//	global_SLAE.X[global_id].x = boundary_value.x;
				//	global_SLAE.F[global_id].x = boundary_value.x;
				//	global_SLAE.Diag[global_id].val[0][0] = 1.0;
				//	global_SLAE.Diag[global_id].val[0][1] = 0.0;
				//	global_SLAE.Diag[global_id].val[0][2] = 0.0;
				//	for (int j = 0; j < global_SLAE.A_down[global_id].size(); j++)
				//	{
				//		global_SLAE.A_down[global_id][j].val[0][0] = 0.0;
				//		global_SLAE.A_down[global_id][j].val[0][1] = 0.0;
				//		global_SLAE.A_down[global_id][j].val[0][2] = 0.0;
				//	}
				//	for (int j = 0; j < global_SLAE.A_up[global_id].size(); j++)
				//	{
				//		global_SLAE.A_up[global_id][j].val[0][0] = 0.0;
				//		global_SLAE.A_up[global_id][j].val[0][1] = 0.0;
				//		global_SLAE.A_up[global_id][j].val[0][2] = 0.0;
				//	}
				//	//simmetrisation
				//	int position = 0;
				//	for (int i = 0; i < global_id; i++)
				//	{
				//		for (int j = 0; j < global_SLAE.id_column_for_A_up[i].size() &&
				//			global_SLAE.id_column_for_A_up[i][j] <= global_id; j++)
				//		{
				//			if (global_id == global_SLAE.id_column_for_A_up[i][j])
				//			{
				//				global_SLAE.F[i].x -= global_SLAE.F[global_id].x *global_SLAE.A_up[i][j].val[0][position];
				//				global_SLAE.F[i].y -= global_SLAE.F[global_id].y *global_SLAE.A_up[i][j].val[1][position];
				//				global_SLAE.F[i].z -= global_SLAE.F[global_id].z *global_SLAE.A_up[i][j].val[2][position];
				//				global_SLAE.A_up[i][j].val[0][position] = 0.0;
				//				global_SLAE.A_up[i][j].val[1][position] = 0.0;
				//				global_SLAE.A_up[i][j].val[2][position] = 0.0;
				//			}
				//		}
				//	}
				//	for (int i = global_id+1; i < global_SLAE.id_column_for_A_down.size(); i++)
				//	{
				//		for (int j = 0; j < global_SLAE.id_column_for_A_down[i].size() &&
				//			global_SLAE.id_column_for_A_down[i][j] <= global_id; j++)
				//		{
				//			if (global_id == global_SLAE.id_column_for_A_down[i][j])
				//			{
				//				global_SLAE.F[i].x -= global_SLAE.F[global_id].x *global_SLAE.A_down[i][j].val[0][position];
				//				global_SLAE.F[i].y -= global_SLAE.F[global_id].y *global_SLAE.A_down[i][j].val[1][position];
				//				global_SLAE.F[i].z -= global_SLAE.F[global_id].z *global_SLAE.A_down[i][j].val[2][position];
				//				global_SLAE.A_down[i][j].val[0][position] = 0.0;
				//				global_SLAE.A_down[i][j].val[1][position] = 0.0;
				//				global_SLAE.A_down[i][j].val[2][position] = 0.0;
				//			}
				//		}
				//	}
				//}
				//if (is_take.y == true)
				//{
				//	global_SLAE.X[global_id].y = boundary_value.y;
				//	global_SLAE.F[global_id].y = boundary_value.y;
				//	global_SLAE.Diag[global_id].val[1][0] = 0.0;
				//	global_SLAE.Diag[global_id].val[1][1] = 1.0;
				//	global_SLAE.Diag[global_id].val[1][2] = 0.0;
				//	for (int j = 0; j < global_SLAE.A_down[global_id].size(); j++)
				//	{
				//		global_SLAE.A_down[global_id][j].val[1][0] = 0.0;
				//		global_SLAE.A_down[global_id][j].val[1][1] = 0.0;
				//		global_SLAE.A_down[global_id][j].val[1][2] = 0.0;
				//	}
				//	for (int j = 0; j < global_SLAE.A_up[global_id].size(); j++)
				//	{
				//		global_SLAE.A_up[global_id][j].val[1][0] = 0.0;
				//		global_SLAE.A_up[global_id][j].val[1][1] = 0.0;
				//		global_SLAE.A_up[global_id][j].val[1][2] = 0.0;
				//	}
				//}
				//if (is_take.z == true)
				//{
				//	global_SLAE.X[global_id].z = boundary_value.z;
				//	global_SLAE.F[global_id].z = boundary_value.z;
				//	global_SLAE.Diag[global_id].val[2][0] = 0.0;
				//	global_SLAE.Diag[global_id].val[2][1] = 0.0;
				//	global_SLAE.Diag[global_id].val[2][2] = 1.0;
				//	for (int j = 0; j < global_SLAE.A_down[global_id].size(); j++)
				//	{
				//		global_SLAE.A_down[global_id][j].val[2][0] = 0.0;
				//		global_SLAE.A_down[global_id][j].val[2][1] = 0.0;
				//		global_SLAE.A_down[global_id][j].val[2][2] = 0.0;
				//	}
				//	for (int j = 0; j < global_SLAE.A_up[global_id].size(); j++)
				//	{
				//		global_SLAE.A_up[global_id][j].val[2][0] = 0.0;
				//		global_SLAE.A_up[global_id][j].val[2][1] = 0.0;
				//		global_SLAE.A_up[global_id][j].val[2][2] = 0.0;
				//	}
				//}

				/*printf_s("complite\n\n");
				char name[1000];
				sprintf_s(name, sizeof(name), "%s/Matrix_with_boundary_id_%d.txt", result_directory, global_id);
				PrintMatrix(name, global_SLAE);*/
			}

		}
		CSSD_Matrix<double, double> newSLAE;
		math::MakeCopyMatrix_A_into_B(global_SLAE, newSLAE);
		//by double matrix;
		if (true)
		{
			printf("First boundary conditions...");
			std::vector<int> use_id;
			//работа со строками
			for (int id_vertex = 0; id_vertex < solver_grid.boundary_vertexes.size(); id_vertex++)
			{
				int global_id = solver_grid.boundary_vertexes[id_vertex].GetDOFInLocalID(0);
				Point<double> vertex = solver_grid.GetCoordinateViaID(solver_grid.boundary_vertexes[id_vertex].GetIdNode(0));
				Point<bool> is_take;
				Point<double> boundary_value = solver_grid.boundary_vertexes[id_vertex].boundary_value(is_take, global_id);

				if (global_id < global_SLAE.GetMatrixSize())
				{
					auto enter_value = [&newSLAE, &use_id](int value_id, double value) ->void
					{
						use_id.push_back(value_id);
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
			for (int id_row = 0; id_row < newSLAE.GetMatrixSize(); id_row++)
			{
				if (id_row % 1000 == 0)
				{
					printf("\tcurrent %d/%d\r", id_row, newSLAE.GetMatrixSize());
				}
				int iterator_in_boundary = 0;
				for (int jj = 0; jj < newSLAE.id_column_for_A_up[id_row].size(); jj++)
				{
					int id_column = newSLAE.id_column_for_A_up[id_row][jj];
					for (; iterator_in_boundary < use_id.size(); iterator_in_boundary++)
					{
						if (use_id[iterator_in_boundary] == id_column)
						{
							newSLAE.F[id_row] -= newSLAE.A_up[id_row][jj] * newSLAE.F[id_column];
							newSLAE.A_up[id_row][jj] = 0;

							break;
						}
						if (use_id[iterator_in_boundary] > id_column)
						{
							break;
						}
					}
				}

				iterator_in_boundary = 0;
				for (int jj = 0; jj < newSLAE.id_column_for_A_down[id_row].size(); jj++)
				{
					int id_column = newSLAE.id_column_for_A_down[id_row][jj];
					for (; iterator_in_boundary < use_id.size(); iterator_in_boundary++)
					{
						if (use_id[iterator_in_boundary] == id_column)
						{
							newSLAE.F[id_row] -= newSLAE.A_down[id_row][jj] * newSLAE.F[id_column];
							newSLAE.A_down[id_row][jj] = 0;

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

		//перестройка стпеней свободы
		CSSD_Matrix<double, double> SLAE_in_old_DOF;
		if (true)
		{
			SLAE_in_old_DOF.Diag.resize(solver_grid.GetDOFsCount() * 3);
			SLAE_in_old_DOF.X.resize(solver_grid.GetDOFsCount() * 3);
			SLAE_in_old_DOF.F.resize(solver_grid.GetDOFsCount() * 3);

			std::vector<std::vector<int>> tmp_down_columns(solver_grid.GetDOFsCount() * 3), tmp_up_columns(solver_grid.GetDOFsCount() * 3);
			std::vector<std::vector<int>> down_columns(solver_grid.GetDOFsCount() * 3), up_columns(solver_grid.GetDOFsCount() * 3);
			for (int id_elem = 0; id_elem < solver_grid.GetElementsCount(); id_elem++)
			{
				if (id_elem % 1000 == 0)
					printf_s("Work with element[%d]\r", id_elem);
				auto element = solver_grid.GetElement(id_elem);
				for (int i = 0; i < element->GetDOFsCount(); i++)
				{
					int global_dof_i = element->GetDOFInLocalID(i);
					for (int ii = 0; ii < 3; ii++)
					{
						int new_global_dof_i = global_dof_i + ii * solver_grid.GetDOFsCount();
						for (int j = i; j < element->GetDOFsCount(); j++)
						{
							int global_dof_j = element->GetDOFInLocalID(j);
							for (int jj = 0; jj < 3; jj++)
							{
								int new_global_dof_j = global_dof_j + jj * solver_grid.GetDOFsCount();
								if (new_global_dof_i != new_global_dof_j)
								{
									if (new_global_dof_i > new_global_dof_j) //down elements of matrix
									{
										tmp_down_columns[new_global_dof_i].push_back(new_global_dof_j);
										tmp_up_columns[new_global_dof_j].push_back(new_global_dof_i);
									}
									else
									{
										tmp_down_columns[new_global_dof_j].push_back(new_global_dof_i);
										tmp_up_columns[new_global_dof_i].push_back(new_global_dof_j);
									}
								}
							}
						}
					}
				}
			}
			printf_s("                                  \r");
			for (int id_string = 0; id_string < tmp_down_columns.size(); id_string++)
			{
				if (id_string % 1000 == 0)
					printf_s("Work with row[%d]\r", id_string);

				math::MakeQuickSort(tmp_down_columns[id_string], 0, (int)tmp_down_columns[id_string].size() - 1);
				math::MakeQuickSort(tmp_up_columns[id_string], 0, (int)tmp_up_columns[id_string].size() - 1);

				math::MakeRemovalOfDuplication(tmp_down_columns[id_string], down_columns[id_string]);
				math::MakeRemovalOfDuplication(tmp_up_columns[id_string], up_columns[id_string]);
			}

			SLAE_in_old_DOF.Initialization(up_columns, down_columns);

			//переписываем матрицы
			for (int id_str = 0; id_str < solver_grid.GetDOFsCount(); id_str++)
			{
				for (int ii = 0; ii < 3; ii++)
				{
					int new_id_str = id_str + ii * solver_grid.GetDOFsCount();

					SLAE_in_old_DOF.X[new_id_str] = newSLAE.X[id_str * 3 + ii];
					SLAE_in_old_DOF.F[new_id_str] = newSLAE.F[id_str * 3 + ii];
					//SLAE_in_old_DOF.Diag[new_id_str] = newSLAE.Diag[id_str * 3 + ii];
				}


				//диагональ
				for (int ii = 0; ii < 3; ii++) 
				{
					int new_id_str = id_str + ii * solver_grid.GetDOFsCount();
					int xyz_id_str = id_str * 3 + ii;
					for (int jj = 0; jj < 3; jj++)
					{
						int new_id_col = id_str + jj * solver_grid.GetDOFsCount();
						int xyz_id_col = id_str * 3 + jj;
						double val;
						if (!newSLAE.GetValue(xyz_id_str, xyz_id_col, val))
						{
							printf_s("ERROR in translation portrait!\n");
							int aa;
							scanf_s("%d", &aa);
						}
						if (!SLAE_in_old_DOF.SetValue(new_id_str, new_id_col,val))
						{
							printf_s("ERROR in translation portrait!\n");
							int aa;
							scanf_s("%d", &aa);
						}
					}
				}

				//верхний треугольник
				for (int j = 0; j < global_SLAE.id_column_for_A_up[id_str].size(); j++)
				{
					int id_col = global_SLAE.id_column_for_A_up[id_str][j];

					for (int ii = 0; ii < 3; ii++)
					{
						int new_id_str = id_str + ii * solver_grid.GetDOFsCount();
						int xyz_id_str = id_str * 3 + ii;

						for (int jj = 0; jj < 3; jj++)
						{
							int new_id_col = id_col + jj * solver_grid.GetDOFsCount();
							int xyz_id_col = id_col * 3 + jj;
							double val;
							if (!newSLAE.GetValue(xyz_id_str, xyz_id_col, val))
							{
								printf_s("ERROR in translation portrait!\n");
								int aa;
								scanf_s("%d", &aa);
							}
							if (!SLAE_in_old_DOF.SetValue(new_id_str, new_id_col, val))
							{
								printf_s("ERROR in translation portrait!\n");
								int aa;
								scanf_s("%d", &aa);
							}
						}
					}
				}

				//нижний треугольник
				for (int j = 0; j < global_SLAE.id_column_for_A_down[id_str].size(); j++)
				{
					int id_col = global_SLAE.id_column_for_A_down[id_str][j];

					for (int ii = 0; ii < 3; ii++)
					{
						int new_id_str = id_str + ii * solver_grid.GetDOFsCount();
						int xyz_id_str = id_str * 3 + ii;

						for (int jj = 0; jj < 3; jj++)
						{
							int new_id_col = id_col + jj * solver_grid.GetDOFsCount();
							int xyz_id_col = id_col * 3 + jj;
							double val;
							if (!newSLAE.GetValue(xyz_id_str, xyz_id_col, val))
							{
								printf_s("ERROR in translation portrait!\n");
								int aa;
								scanf_s("%d", &aa);
							}
							if (!SLAE_in_old_DOF.SetValue(new_id_str, new_id_col, val))
							{
								printf_s("ERROR in translation portrait!\n");
								int aa;
								scanf_s("%d", &aa);
							}
						}
					}
				}
			}

			//печатаем портреты
			//char name_Tensor[1000];
			//sprintf_s(name_Tensor, "%s/global_SLAE.txt", result_directory);
			//PrintPortrait(name_Tensor, global_SLAE);
//
			//char name_Old[1000];
			//sprintf_s(name_Old, "%s/SLAE_in_old_DOF.txt", result_directory);
			//PrintPortrait(name_Old, SLAE_in_old_DOF);
		}

		clock_t t_before = clock();
		printf("TIME = %lf sec\n", (t_before - t_after) * 1.0 / CLK_TCK);
		double end = omp_get_wtime();
		printf("start = %.16g\nend = %.16g\ndiff = %.16g\n", start, end, end - start);

		if (false)
		{
			printf("Soluting SLAY... (%d)\n", newSLAE.GetMatrixSize());
			int MaxSize = newSLAE.GetMatrixSize();
			MaxSize = MaxSize / 10 < 1000 ? 1000 : MaxSize / 10;
			std::vector<double> best_solution;
			math::MakeCopyVector_A_into_B(newSLAE.X, best_solution);
			double current_residual = 1, best_residual = 1e+25;
			int MAX_STEPS = 4;
			int ii = 1;
			CSSD_Matrix<double, double> Precond;
			Precond.PrecondorSSOR(0.75, newSLAE);
			for (int i = 0; i <= MAX_STEPS; i++)
			{
				double needed_residual = pow(10., -1 * (ii + 1));
				needed_residual = needed_residual < MIN_RESIDUAL ? MIN_RESIDUAL : needed_residual;
				//current_residual /= 2.;
				printf_s("//---> I = %d/%d (full size %d) - needed residual %.2e\n", i, MAX_STEPS, newSLAE.GetMatrixSize(), needed_residual);

				current_residual = abs(newSLAE.MSG_PreconditioningSSOR(MaxSize, 
					needed_residual,
					Precond));

				if (current_residual < needed_residual)
				{
					i = 0;
					ii++;
				}
				if (current_residual > best_residual * 100)
				{
					break;
				}
				if (current_residual <= MIN_RESIDUAL)
				{
					best_residual = current_residual;
					break;
				}
				if (current_residual < best_residual)
				{
					best_residual = current_residual;
					math::MakeCopyVector_A_into_B(newSLAE.X, best_solution);
				}
			}
			math::MakeCopyVector_A_into_B(best_solution, newSLAE.X);
			printf_s("//---> BEST residual %.2e\n", best_residual);
			//переписываем лучшее решение в векторный вид
			math::MakeCopyVector_A_into_B(newSLAE.X, Solution);

			if (best_residual >= 1)
			{
				printf_s("There are problems in the solution\n");
				int a;
				scanf_s("%d", &a);
			}
		}
		if (true)
		{
			//SLAE_in_old_DOF.print_logs = true;
			printf("Soluting SLAY... (%d)\n", SLAE_in_old_DOF.GetMatrixSize());
			int MaxSize = SLAE_in_old_DOF.GetMatrixSize();
			MaxSize = MaxSize / 10 < 1000 ? 1000 : MaxSize / 10;
			std::vector<double> best_solution;
			math::MakeCopyVector_A_into_B(SLAE_in_old_DOF.X, best_solution);
			double current_residual = 1, best_residual = 1e+25;
			int MAX_STEPS = 4;
			int ii = 1;
			CSSD_Matrix<double, double> Precond;
			Precond.PrecondorSSOR(0.75, SLAE_in_old_DOF);
			SLAE_in_old_DOF.print_logs = is_print_logFile;

			for (int i = 0; i <= MAX_STEPS; i++)
			{
				double needed_residual = pow(10., -1 * (ii + 1));
				needed_residual = needed_residual < MIN_RESIDUAL ? MIN_RESIDUAL : needed_residual;
				//current_residual /= 2.;
				if (is_print_logFile) printf_s("//---> I = %d/%d (full size %d) - needed residual %.2e\n", i, MAX_STEPS, newSLAE.GetMatrixSize(), needed_residual);

				current_residual = abs(SLAE_in_old_DOF.MSG_PreconditioningSSOR(MaxSize,
					needed_residual,
					Precond));

				if (current_residual < needed_residual)
				{
					i = 0;
					ii++;
				}
				else {
					current_residual = abs(SLAE_in_old_DOF.BCG_Stab2(MaxSize, needed_residual));
				}
				if (current_residual > best_residual * 100)
				{
					break;
				}
				if (current_residual <= MIN_RESIDUAL)
				{
					best_residual = current_residual;
					break;
				}
				if (current_residual < best_residual)
				{
					best_residual = current_residual;
					math::MakeCopyVector_A_into_B(SLAE_in_old_DOF.X, best_solution);
				}
			}
			math::MakeCopyVector_A_into_B(best_solution, SLAE_in_old_DOF.X);
			printf_s("//---> BEST residual %.2e\n", best_residual);
			//переписываем лучшее решение в векторный вид
			//math::MakeCopyVector_A_into_B(SLAE_in_old_DOF.X, Solution);
			Solution.resize(solver_grid.GetDOFsCount());
			for (int i = 0; i < Solution.size(); i++)
			{
				Solution[i].x = SLAE_in_old_DOF.X[i + 0 * Solution.size()];
				Solution[i].y = SLAE_in_old_DOF.X[i + 1 * Solution.size()];
				Solution[i].z = SLAE_in_old_DOF.X[i + 2 * Solution.size()];
			}

			if (best_residual >= 1)
			{
				printf_s("There are problems in the solution\n");
				int a;
				scanf_s("%d", &a);
			}
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

	template <typename Dirichlet, typename Neumann>
	void FEM_1D_forElasticDeformation(
		bool is_print_logFile, //input
		double MIN_RESIDUAL,
		math::SimpleGrid& geo_grid, //input
		std::vector<Dirichlet>& first_boundary,//input
		std::vector<Neumann>& second_boundary,//input
		std::function<Point<double>(bool &is_solve, int id_element, Point<double> X)> &sourse,
		char* result_directory, //output
		FEM::Grid_1D_forMech& solver_grid, //output
		std::vector<Point<double>>& Solution, //output
		CSSD_Matrix<Tensor2Rank3D, Point<double>>& Stiffness_matrix //output
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

		if (solver_grid.GetElementsCount() == 0)
		{
			if (is_print_logFile) printf("Initialization of grid...\n");
			solver_grid.Initialization(geo_grid, first_boundary, second_boundary);
			if (is_print_logFile) printf_s("complite                   \n\n");
		}
		else {
			//update value functions for 1st boundary

			int id = 0;
			for (int id_type = 0; id_type < first_boundary.size(); id_type++)
			{
				for (int id_vertex = 0; id_vertex < first_boundary[id_type].id_vertexes.size(); id_vertex++)
				{
					solver_grid.boundary_vertexes[id].boundary_value = first_boundary[id_type].value;
					id++;
				}
			}
		}

		CSSD_Matrix<Tensor2Rank3D, Point<double>> global_SLAE;

		if (Stiffness_matrix.GetMatrixSize() == 0)
		{
			if (is_print_logFile) printf("Creation the SLAE portrait...");
			solver_grid.CreationPortrait(global_SLAE);
			if (is_print_logFile) printf_s("complite\n\n");

			//SLAE assembling
			{
				if (is_print_logFile) printf("SLAE assembling...\n");
				std::vector<DenseMatrix<Tensor2Rank3D, Point<double>>> local_SLAE(solver_grid.GetElementsCount());
				for (int id_elem = 0; id_elem < solver_grid.GetElementsCount(); id_elem++)
				{
					if (id_elem % 1000 == 0)
						printf("Solve element[%d]\r", id_elem);
					auto element = solver_grid.GetElement(id_elem);
					std::function< std::vector<std::vector<double>>(Point<double> X) > koefD = [&](Point<double> X)->std::vector<std::vector<double>>
					{
						return solver_grid.GetDomain(solver_grid.GetElement(id_elem)->GetIdDomain())->forMech.GetD(3);
					};
					element->SolveLocalMatrix(local_SLAE[id_elem], koefD);

				}
				for (int id_elem = 0; id_elem < solver_grid.GetElementsCount(); id_elem++)
				{
					if (id_elem % 1000 == 0)
						printf("Add the local matrix of element[%d]\r", id_elem);
					global_SLAE.SummPartOfMatrix(local_SLAE[id_elem], *solver_grid.GetElementDOFs(id_elem));
				}
				if (is_print_logFile) printf_s("                                                                                    \r");
				if (is_print_logFile) printf_s("complite\n\n");
			}
			//copy base matrix
			Stiffness_matrix.SetMatrix(global_SLAE);
		}
		else {
			global_SLAE.SetMatrix(Stiffness_matrix);
		}

		//add right side
		{
			if (is_print_logFile) printf("SLAE assembling (right side...\n");
			std::vector<std::vector<Point<double>>> local_right_side(solver_grid.GetElementsCount());
			for (int id_elem = 0; id_elem < solver_grid.GetElementsCount(); id_elem++)
			{
				if (id_elem % 1000 == 0)
					printf("Solve element[%d]\r", id_elem);
				auto element = solver_grid.GetElement(id_elem);
				bool is_solve = false;
				sourse(is_solve, id_elem, element->GetWeightCentr());
				if (is_solve)
				{
					std::function< Point<double>(Point<double> X) > sourse_elem = [&is_solve, id_elem, &sourse](Point<double> X)->Point<double>
					{
						return sourse(is_solve, id_elem, X);
					};
					element->SolveRightSide(local_right_side[id_elem], sourse_elem);
				}

			}
			for (int id_elem = 0; id_elem < solver_grid.GetElementsCount(); id_elem++)
			{
				if (id_elem % 1000 == 0)
					printf("Add the local matrix of element[%d]\r", id_elem);
				bool is_solve = false;
				sourse(is_solve, id_elem, solver_grid.GetElement(id_elem)->GetWeightCentr());
				if (is_solve)
				{
					global_SLAE.SummPartOfVector(local_right_side[id_elem], *solver_grid.GetElementDOFs(id_elem));
				}
			}
			if (is_print_logFile) printf_s("                                                                                    \r");
			if (is_print_logFile) printf_s("complite\n\n");
		}

		//Boundary condition
		if(false){
			global_SLAE.MakeDivideOnSqrtFormSumm();

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

				//if (is_take.x == true)
				//{
				//	global_SLAE.X[global_id].x = boundary_value.x;
				//	global_SLAE.F[global_id].x = boundary_value.x;
				//	global_SLAE.Diag[global_id].val[0][0] = 1.0;
				//	global_SLAE.Diag[global_id].val[0][1] = 0.0;
				//	global_SLAE.Diag[global_id].val[0][2] = 0.0;
				//	for (int j = 0; j < global_SLAE.A_down[global_id].size(); j++)
				//	{
				//		global_SLAE.A_down[global_id][j].val[0][0] = 0.0;
				//		global_SLAE.A_down[global_id][j].val[0][1] = 0.0;
				//		global_SLAE.A_down[global_id][j].val[0][2] = 0.0;
				//	}
				//	for (int j = 0; j < global_SLAE.A_up[global_id].size(); j++)
				//	{
				//		global_SLAE.A_up[global_id][j].val[0][0] = 0.0;
				//		global_SLAE.A_up[global_id][j].val[0][1] = 0.0;
				//		global_SLAE.A_up[global_id][j].val[0][2] = 0.0;
				//	}
				//	//simmetrisation
				//	int position = 0;
				//	for (int i = 0; i < global_id; i++)
				//	{
				//		for (int j = 0; j < global_SLAE.id_column_for_A_up[i].size() &&
				//			global_SLAE.id_column_for_A_up[i][j] <= global_id; j++)
				//		{
				//			if (global_id == global_SLAE.id_column_for_A_up[i][j])
				//			{
				//				global_SLAE.F[i].x -= global_SLAE.F[global_id].x *global_SLAE.A_up[i][j].val[0][position];
				//				global_SLAE.F[i].y -= global_SLAE.F[global_id].y *global_SLAE.A_up[i][j].val[1][position];
				//				global_SLAE.F[i].z -= global_SLAE.F[global_id].z *global_SLAE.A_up[i][j].val[2][position];
				//				global_SLAE.A_up[i][j].val[0][position] = 0.0;
				//				global_SLAE.A_up[i][j].val[1][position] = 0.0;
				//				global_SLAE.A_up[i][j].val[2][position] = 0.0;
				//			}
				//		}
				//	}
				//	for (int i = global_id+1; i < global_SLAE.id_column_for_A_down.size(); i++)
				//	{
				//		for (int j = 0; j < global_SLAE.id_column_for_A_down[i].size() &&
				//			global_SLAE.id_column_for_A_down[i][j] <= global_id; j++)
				//		{
				//			if (global_id == global_SLAE.id_column_for_A_down[i][j])
				//			{
				//				global_SLAE.F[i].x -= global_SLAE.F[global_id].x *global_SLAE.A_down[i][j].val[0][position];
				//				global_SLAE.F[i].y -= global_SLAE.F[global_id].y *global_SLAE.A_down[i][j].val[1][position];
				//				global_SLAE.F[i].z -= global_SLAE.F[global_id].z *global_SLAE.A_down[i][j].val[2][position];
				//				global_SLAE.A_down[i][j].val[0][position] = 0.0;
				//				global_SLAE.A_down[i][j].val[1][position] = 0.0;
				//				global_SLAE.A_down[i][j].val[2][position] = 0.0;
				//			}
				//		}
				//	}
				//}
				//if (is_take.y == true)
				//{
				//	global_SLAE.X[global_id].y = boundary_value.y;
				//	global_SLAE.F[global_id].y = boundary_value.y;
				//	global_SLAE.Diag[global_id].val[1][0] = 0.0;
				//	global_SLAE.Diag[global_id].val[1][1] = 1.0;
				//	global_SLAE.Diag[global_id].val[1][2] = 0.0;
				//	for (int j = 0; j < global_SLAE.A_down[global_id].size(); j++)
				//	{
				//		global_SLAE.A_down[global_id][j].val[1][0] = 0.0;
				//		global_SLAE.A_down[global_id][j].val[1][1] = 0.0;
				//		global_SLAE.A_down[global_id][j].val[1][2] = 0.0;
				//	}
				//	for (int j = 0; j < global_SLAE.A_up[global_id].size(); j++)
				//	{
				//		global_SLAE.A_up[global_id][j].val[1][0] = 0.0;
				//		global_SLAE.A_up[global_id][j].val[1][1] = 0.0;
				//		global_SLAE.A_up[global_id][j].val[1][2] = 0.0;
				//	}
				//}
				//if (is_take.z == true)
				//{
				//	global_SLAE.X[global_id].z = boundary_value.z;
				//	global_SLAE.F[global_id].z = boundary_value.z;
				//	global_SLAE.Diag[global_id].val[2][0] = 0.0;
				//	global_SLAE.Diag[global_id].val[2][1] = 0.0;
				//	global_SLAE.Diag[global_id].val[2][2] = 1.0;
				//	for (int j = 0; j < global_SLAE.A_down[global_id].size(); j++)
				//	{
				//		global_SLAE.A_down[global_id][j].val[2][0] = 0.0;
				//		global_SLAE.A_down[global_id][j].val[2][1] = 0.0;
				//		global_SLAE.A_down[global_id][j].val[2][2] = 0.0;
				//	}
				//	for (int j = 0; j < global_SLAE.A_up[global_id].size(); j++)
				//	{
				//		global_SLAE.A_up[global_id][j].val[2][0] = 0.0;
				//		global_SLAE.A_up[global_id][j].val[2][1] = 0.0;
				//		global_SLAE.A_up[global_id][j].val[2][2] = 0.0;
				//	}
				//}

				/*printf_s("complite\n\n");
				char name[1000];
				sprintf_s(name, sizeof(name), "%s/Matrix_with_boundary_id_%d.txt", result_directory, global_id);
				PrintMatrix(name, global_SLAE);*/
			}

		}
		CSSD_Matrix<double, double> newSLAE;
		math::MakeCopyMatrix_A_into_B(global_SLAE, newSLAE);
		//by double matrix;
		if (true)
		{
			if (is_print_logFile) printf("First boundary conditions...");
			std::vector<int> use_id;
			//работа со строками
			for (int id_vertex = 0; id_vertex < solver_grid.boundary_vertexes.size(); id_vertex++)
			{
				int global_id = solver_grid.boundary_vertexes[id_vertex].GetDOFInLocalID(0);
				Point<double> vertex = solver_grid.GetCoordinateViaID(solver_grid.boundary_vertexes[id_vertex].GetIdNode(0));
				Point<bool> is_take;
				Point<double> boundary_value = solver_grid.boundary_vertexes[id_vertex].boundary_value(is_take);

				if (global_id < global_SLAE.GetMatrixSize())
				{
					auto enter_value = [&newSLAE, &use_id](int value_id, double value) ->void
					{
						use_id.push_back(value_id);
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
			for (int id_row = 0; id_row < newSLAE.GetMatrixSize(); id_row++)
			{
				if (id_row % 1000 == 0)
				{
					if (is_print_logFile) printf("\tcurrent %d/%d\r", id_row, newSLAE.GetMatrixSize());
				}
				int iterator_in_boundary = 0;
				for (int jj = 0; jj < newSLAE.id_column_for_A_up[id_row].size(); jj++)
				{
					int id_column = newSLAE.id_column_for_A_up[id_row][jj];
					for (; iterator_in_boundary < use_id.size(); iterator_in_boundary++)
					{
						if (use_id[iterator_in_boundary] == id_column)
						{
							newSLAE.F[id_row] -= newSLAE.A_up[id_row][jj] * newSLAE.F[id_column];
							newSLAE.A_up[id_row][jj] = 0;

							break;
						}
						if (use_id[iterator_in_boundary] > id_column)
						{
							break;
						}
					}
				}

				iterator_in_boundary = 0;
				for (int jj = 0; jj < newSLAE.id_column_for_A_down[id_row].size(); jj++)
				{
					int id_column = newSLAE.id_column_for_A_down[id_row][jj];
					for (; iterator_in_boundary < use_id.size(); iterator_in_boundary++)
					{
						if (use_id[iterator_in_boundary] == id_column)
						{
							newSLAE.F[id_row] -= newSLAE.A_down[id_row][jj] * newSLAE.F[id_column];
							newSLAE.A_down[id_row][jj] = 0;

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

		clock_t t_before = clock();
		if (is_print_logFile) printf("TIME = %lf sec\n", (t_before - t_after) * 1.0 / CLK_TCK);
		double end = omp_get_wtime();
		if (is_print_logFile) printf("start = %.16g\nend = %.16g\ndiff = %.16g\n", start, end, end - start);

		if (true)
		{
			if (is_print_logFile) printf("Soluting SLAY... (%d)\n", newSLAE.GetMatrixSize());
			int MaxSize = newSLAE.GetMatrixSize();
			MaxSize = MaxSize / 10 < 1000 ? 1000 : MaxSize / 10;
			std::vector<double> best_solution;
			math::MakeCopyVector_A_into_B(newSLAE.X, best_solution);
			double current_residual = 1, best_residual = 1e+25;
			int MAX_STEPS = 4;
			int ii = 1;
			CSSD_Matrix<double, double> Precond;
			Precond.PrecondorSSOR(0.75, newSLAE);
			newSLAE.print_logs = is_print_logFile;

			for (int i = 0; i <= MAX_STEPS; i++)
			{
				double needed_residual = pow(10., -1 * (ii + 1));
				needed_residual = needed_residual < MIN_RESIDUAL ? MIN_RESIDUAL : needed_residual;
				//current_residual /= 2.;
				if (is_print_logFile) printf_s("//---> I = %d/%d (full size %d) - needed residual %.2e\n", i, MAX_STEPS, newSLAE.GetMatrixSize(), needed_residual);

				current_residual = abs(newSLAE.MSG_PreconditioningSSOR(MaxSize,
					needed_residual,
					Precond));

				if (current_residual < needed_residual)
				{
					i = 0;
					ii++;
				}
				if (current_residual > best_residual * 100)
				{
					break;
				}
				if (current_residual <= MIN_RESIDUAL)
				{
					best_residual = current_residual;
					break;
				}
				if (current_residual < best_residual)
				{
					best_residual = current_residual;
					math::MakeCopyVector_A_into_B(newSLAE.X, best_solution);
				}
			}
			if(current_residual > 1e-15)
				printf_s("//---> BEST residual %.2e\n", best_residual);

			math::MakeCopyVector_A_into_B(best_solution, newSLAE.X);
			if (is_print_logFile) printf_s("//---> BEST residual %.2e\n", best_residual);
			//переписываем лучшее решение в векторный вид
			math::MakeCopyVector_A_into_B(newSLAE.X, Solution);

			if (best_residual >= 1)
			{
				printf_s("There are problems in the solution\n");
				int a;
				scanf_s("%d", &a);
			}
		}

		if (is_print_logFile) printf("\tcomplit\n\n");
		t_before = clock();
		if (is_print_logFile) printf("TIME = %lf sec\n", (t_before - t_after) * 1.0 / CLK_TCK);
		end = omp_get_wtime();
		if (is_print_logFile) printf("start = %.16g\nend = %.16g\ndiff = %.16g\n", start, end, end - start);



		//update grid
		/*{
			std::vector<Point<double>> new_x(solver_grid.GetVertexCount());
			for (int id_vertex = 0; id_vertex < solver_grid.GetVertexCount(); id_vertex++)
			{
				new_x[id_vertex] = solver_grid.GetCoordinateViaID(id_vertex) + global_SLAE.X[id_vertex];
			}
			solver_grid.UpdateCoordinates(new_x);
		}*/


		if (is_print_logFile) printf_s("complite\n\n");
	};

	template <typename Dirichlet, typename Neumann>
	void FEM_2D_forElasticDeformation(
		bool is_print_logFile, //input
		double MIN_RESIDUAL,
		math::SimpleGrid& geo_grid, //input
		std::vector<Dirichlet>& first_boundary,//input
		std::vector<Neumann>& second_boundary,//input
		std::function<Point<double>(bool& is_solve, int id_element, Point<double> X)>& sourse,
		char* result_directory, //output
		FEM::Grid_2D_forMech& solver_grid, //output
		std::vector<Point<double>>& Solution, //output
		CSSD_Matrix<Tensor2Rank3D, Point<double>>& Stiffness_matrix //output
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

		if (solver_grid.GetElementsCount() == 0)
		{
			printf("Initialization of grid...\n");
			solver_grid.Initialization(geo_grid, first_boundary, second_boundary);
			printf_s("complite                   \n\n");
		}
		else {
			//update value functions for 1st boundary

			int id = 0;
			for (int id_type = 0; id_type < first_boundary.size(); id_type++)
			{
				for (int id_vertex = 0; id_vertex < first_boundary[id_type].id_vertexes.size(); id_vertex++)
				{
					solver_grid.boundary_vertexes[id].boundary_value = first_boundary[id_type].value;
					id++;
				}
			}
		}

		CSSD_Matrix<Tensor2Rank3D, Point<double>> global_SLAE;

		if (Stiffness_matrix.GetMatrixSize() == 0)
		{
			printf("Creation the SLAE portrait...");
			solver_grid.CreationPortrait(global_SLAE);
			printf_s("complite\n\n");

			//SLAE assembling
			{
				printf("SLAE assembling...\n");
				std::vector<DenseMatrix<Tensor2Rank3D, Point<double>>> local_SLAE(solver_grid.GetElementsCount());
				for (int id_elem = 0; id_elem < solver_grid.GetElementsCount(); id_elem++)
				{
					if (id_elem % 1000 == 0)
						printf("Solve element[%d]\r", id_elem);
					auto element = solver_grid.GetElement(id_elem);
					std::function< std::vector<std::vector<double>>(Point<double> X) > koefD = [&](Point<double> X)->std::vector<std::vector<double>>
					{
						return solver_grid.GetDomain(solver_grid.GetElement(id_elem)->GetIdDomain())->forMech.GetD(3);
					};
					element->SolveLocalMatrix(local_SLAE[id_elem], koefD);

				}
				for (int id_elem = 0; id_elem < solver_grid.GetElementsCount(); id_elem++)
				{
					if (id_elem % 1000 == 0)
						printf("Add the local matrix of element[%d]\r", id_elem);
					global_SLAE.SummPartOfMatrix(local_SLAE[id_elem], *solver_grid.GetElementDOFs(id_elem));
				}
				printf_s("                                                                                    \r");
				printf_s("complite\n\n");
			}
			//copy base matrix
			Stiffness_matrix.SetMatrix(global_SLAE);
		}
		else {
			global_SLAE.SetMatrix(Stiffness_matrix);
		}

		for (int i = 0; i < global_SLAE.X.size(); i++)
		{
			global_SLAE.X[i] = 0.0;
		}

		//add right side
		{
			if (is_print_logFile)printf("SLAE assembling (right side...\n");
			std::vector<std::vector<Point<double>>> local_RightSide(solver_grid.GetElementsCount());
			for (int id_elem = 0; id_elem < solver_grid.GetElementsCount(); id_elem++)
			{
				if (id_elem % 1000 == 0)
					printf("Solve element[%d]\r", id_elem);
				auto element = solver_grid.GetElement(id_elem);
				bool is_solve = false;
				sourse(is_solve, id_elem, element->GetWeightCentr());
				if (is_solve)
				{
					std::function< Point<double>(Point<double> X) > sourse_elem = [&is_solve, id_elem, &sourse](Point<double> X)->Point<double>
					{
						return sourse(is_solve, id_elem, X);
					};
					element->SolveRightSide(local_RightSide[id_elem], sourse_elem);
				}

			}
			if (is_print_logFile)for (int id_elem = 0; id_elem < solver_grid.GetElementsCount(); id_elem++)
			{
				if (id_elem % 1000 == 0)
					printf("Add the local matrix of element[%d]\r", id_elem);
				bool is_solve = false;
				sourse(is_solve, id_elem, solver_grid.GetElement(id_elem)->GetWeightCentr());
				if (is_solve)
				{
					global_SLAE.SummPartOfVector(local_RightSide[id_elem], *solver_grid.GetElementDOFs(id_elem));
				}
			}
			if (is_print_logFile) printf_s("                                                                                    \r");
			if (is_print_logFile) printf_s("complite\n\n");
		}

		//Boundary condition
		if(false){
			global_SLAE.MakeDivideOnSqrtFormSumm();

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
				Point<double> boundary_value = boundary->boundary_value(is_take, solver_grid.boundary_vertexes[id_vertex].GetIdNode(0));
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

				//if (is_take.x == true)
				//{
				//	global_SLAE.X[global_id].x = boundary_value.x;
				//	global_SLAE.F[global_id].x = boundary_value.x;
				//	global_SLAE.Diag[global_id].val[0][0] = 1.0;
				//	global_SLAE.Diag[global_id].val[0][1] = 0.0;
				//	global_SLAE.Diag[global_id].val[0][2] = 0.0;
				//	for (int j = 0; j < global_SLAE.A_down[global_id].size(); j++)
				//	{
				//		global_SLAE.A_down[global_id][j].val[0][0] = 0.0;
				//		global_SLAE.A_down[global_id][j].val[0][1] = 0.0;
				//		global_SLAE.A_down[global_id][j].val[0][2] = 0.0;
				//	}
				//	for (int j = 0; j < global_SLAE.A_up[global_id].size(); j++)
				//	{
				//		global_SLAE.A_up[global_id][j].val[0][0] = 0.0;
				//		global_SLAE.A_up[global_id][j].val[0][1] = 0.0;
				//		global_SLAE.A_up[global_id][j].val[0][2] = 0.0;
				//	}
				//	//simmetrisation
				//	int position = 0;
				//	for (int i = 0; i < global_id; i++)
				//	{
				//		for (int j = 0; j < global_SLAE.id_column_for_A_up[i].size() &&
				//			global_SLAE.id_column_for_A_up[i][j] <= global_id; j++)
				//		{
				//			if (global_id == global_SLAE.id_column_for_A_up[i][j])
				//			{
				//				global_SLAE.F[i].x -= global_SLAE.F[global_id].x *global_SLAE.A_up[i][j].val[0][position];
				//				global_SLAE.F[i].y -= global_SLAE.F[global_id].y *global_SLAE.A_up[i][j].val[1][position];
				//				global_SLAE.F[i].z -= global_SLAE.F[global_id].z *global_SLAE.A_up[i][j].val[2][position];
				//				global_SLAE.A_up[i][j].val[0][position] = 0.0;
				//				global_SLAE.A_up[i][j].val[1][position] = 0.0;
				//				global_SLAE.A_up[i][j].val[2][position] = 0.0;
				//			}
				//		}
				//	}
				//	for (int i = global_id+1; i < global_SLAE.id_column_for_A_down.size(); i++)
				//	{
				//		for (int j = 0; j < global_SLAE.id_column_for_A_down[i].size() &&
				//			global_SLAE.id_column_for_A_down[i][j] <= global_id; j++)
				//		{
				//			if (global_id == global_SLAE.id_column_for_A_down[i][j])
				//			{
				//				global_SLAE.F[i].x -= global_SLAE.F[global_id].x *global_SLAE.A_down[i][j].val[0][position];
				//				global_SLAE.F[i].y -= global_SLAE.F[global_id].y *global_SLAE.A_down[i][j].val[1][position];
				//				global_SLAE.F[i].z -= global_SLAE.F[global_id].z *global_SLAE.A_down[i][j].val[2][position];
				//				global_SLAE.A_down[i][j].val[0][position] = 0.0;
				//				global_SLAE.A_down[i][j].val[1][position] = 0.0;
				//				global_SLAE.A_down[i][j].val[2][position] = 0.0;
				//			}
				//		}
				//	}
				//}
				//if (is_take.y == true)
				//{
				//	global_SLAE.X[global_id].y = boundary_value.y;
				//	global_SLAE.F[global_id].y = boundary_value.y;
				//	global_SLAE.Diag[global_id].val[1][0] = 0.0;
				//	global_SLAE.Diag[global_id].val[1][1] = 1.0;
				//	global_SLAE.Diag[global_id].val[1][2] = 0.0;
				//	for (int j = 0; j < global_SLAE.A_down[global_id].size(); j++)
				//	{
				//		global_SLAE.A_down[global_id][j].val[1][0] = 0.0;
				//		global_SLAE.A_down[global_id][j].val[1][1] = 0.0;
				//		global_SLAE.A_down[global_id][j].val[1][2] = 0.0;
				//	}
				//	for (int j = 0; j < global_SLAE.A_up[global_id].size(); j++)
				//	{
				//		global_SLAE.A_up[global_id][j].val[1][0] = 0.0;
				//		global_SLAE.A_up[global_id][j].val[1][1] = 0.0;
				//		global_SLAE.A_up[global_id][j].val[1][2] = 0.0;
				//	}
				//}
				//if (is_take.z == true)
				//{
				//	global_SLAE.X[global_id].z = boundary_value.z;
				//	global_SLAE.F[global_id].z = boundary_value.z;
				//	global_SLAE.Diag[global_id].val[2][0] = 0.0;
				//	global_SLAE.Diag[global_id].val[2][1] = 0.0;
				//	global_SLAE.Diag[global_id].val[2][2] = 1.0;
				//	for (int j = 0; j < global_SLAE.A_down[global_id].size(); j++)
				//	{
				//		global_SLAE.A_down[global_id][j].val[2][0] = 0.0;
				//		global_SLAE.A_down[global_id][j].val[2][1] = 0.0;
				//		global_SLAE.A_down[global_id][j].val[2][2] = 0.0;
				//	}
				//	for (int j = 0; j < global_SLAE.A_up[global_id].size(); j++)
				//	{
				//		global_SLAE.A_up[global_id][j].val[2][0] = 0.0;
				//		global_SLAE.A_up[global_id][j].val[2][1] = 0.0;
				//		global_SLAE.A_up[global_id][j].val[2][2] = 0.0;
				//	}
				//}

				/*printf_s("complite\n\n");
				char name[1000];
				sprintf_s(name, sizeof(name), "%s/Matrix_with_boundary_id_%d.txt", result_directory, global_id);
				PrintMatrix(name, global_SLAE);*/
			}

		}
		CSSD_Matrix<double, double> newSLAE;
		math::MakeCopyMatrix_A_into_B(global_SLAE, newSLAE);
		//by double matrix;
		if (true)
		{
			if (is_print_logFile) printf("First boundary conditions...");
			std::vector<int> use_id;
			//работа со строками
			for (int id_vertex = 0; id_vertex < solver_grid.boundary_vertexes.size(); id_vertex++)
			{
				int global_id = solver_grid.boundary_vertexes[id_vertex].GetDOFInLocalID(0);
				Point<double> vertex = solver_grid.GetCoordinateViaID(solver_grid.boundary_vertexes[id_vertex].GetIdNode(0));
				Point<bool> is_take;
				Point<double> boundary_value = solver_grid.boundary_vertexes[id_vertex].boundary_value(is_take, global_id);

				if (global_id < global_SLAE.GetMatrixSize())
				{
					auto enter_value = [&newSLAE, &use_id](int value_id, double value) ->void
					{
						use_id.push_back(value_id);
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
			for (int id_row = 0; id_row < newSLAE.GetMatrixSize(); id_row++)
			{
				if (id_row % 1000 == 0)
				{
					if (is_print_logFile) printf("\tcurrent %d/%d\r", id_row, newSLAE.GetMatrixSize());
				}
				int iterator_in_boundary = 0;
				for (int jj = 0; jj < newSLAE.id_column_for_A_up[id_row].size(); jj++)
				{
					int id_column = newSLAE.id_column_for_A_up[id_row][jj];
					for (; iterator_in_boundary < use_id.size(); iterator_in_boundary++)
					{
						if (use_id[iterator_in_boundary] == id_column)
						{
							newSLAE.F[id_row] -= newSLAE.A_up[id_row][jj] * newSLAE.F[id_column];
							newSLAE.A_up[id_row][jj] = 0;

							break;
						}
						if (use_id[iterator_in_boundary] > id_column)
						{
							break;
						}
					}
				}

				iterator_in_boundary = 0;
				for (int jj = 0; jj < newSLAE.id_column_for_A_down[id_row].size(); jj++)
				{
					int id_column = newSLAE.id_column_for_A_down[id_row][jj];
					for (; iterator_in_boundary < use_id.size(); iterator_in_boundary++)
					{
						if (use_id[iterator_in_boundary] == id_column)
						{
							newSLAE.F[id_row] -= newSLAE.A_down[id_row][jj] * newSLAE.F[id_column];
							newSLAE.A_down[id_row][jj] = 0;

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

		//перестройка стпеней свободы
		CSSD_Matrix<double, double> SLAE_in_old_DOF;
		if (true)
		{
			SLAE_in_old_DOF.Diag.resize(solver_grid.GetDOFsCount() * 3);
			SLAE_in_old_DOF.X.resize(solver_grid.GetDOFsCount() * 3);
			SLAE_in_old_DOF.F.resize(solver_grid.GetDOFsCount() * 3);

			std::vector<std::vector<int>> tmp_down_columns(solver_grid.GetDOFsCount() * 3), tmp_up_columns(solver_grid.GetDOFsCount() * 3);
			std::vector<std::vector<int>> down_columns(solver_grid.GetDOFsCount() * 3), up_columns(solver_grid.GetDOFsCount() * 3);
			for (int id_elem = 0; id_elem < solver_grid.GetElementsCount(); id_elem++)
			{
				if (id_elem % 1000 == 0)
					printf_s("Work with element[%d]\r", id_elem);
				auto element = solver_grid.GetElement(id_elem);
				for (int i = 0; i < element->GetDOFsCount(); i++)
				{
					int global_dof_i = element->GetDOFInLocalID(i);
					for (int ii = 0; ii < 3; ii++)
					{
						int new_global_dof_i = global_dof_i + ii * solver_grid.GetDOFsCount();
						for (int j = i; j < element->GetDOFsCount(); j++)
						{
							int global_dof_j = element->GetDOFInLocalID(j);
							for (int jj = 0; jj < 3; jj++)
							{
								int new_global_dof_j = global_dof_j + jj * solver_grid.GetDOFsCount();
								if (new_global_dof_i != new_global_dof_j)
								{
									if (new_global_dof_i > new_global_dof_j) //down elements of matrix
									{
										tmp_down_columns[new_global_dof_i].push_back(new_global_dof_j);
										tmp_up_columns[new_global_dof_j].push_back(new_global_dof_i);
									}
									else
									{
										tmp_down_columns[new_global_dof_j].push_back(new_global_dof_i);
										tmp_up_columns[new_global_dof_i].push_back(new_global_dof_j);
									}
								}
							}
						}
					}
				}
			}
			printf_s("                                  \r");
			for (int id_string = 0; id_string < tmp_down_columns.size(); id_string++)
			{
				if (id_string % 1000 == 0)
					printf_s("Work with row[%d]\r", id_string);

				math::MakeQuickSort(tmp_down_columns[id_string], 0, (int)tmp_down_columns[id_string].size() - 1);
				math::MakeQuickSort(tmp_up_columns[id_string], 0, (int)tmp_up_columns[id_string].size() - 1);

				math::MakeRemovalOfDuplication(tmp_down_columns[id_string], down_columns[id_string]);
				math::MakeRemovalOfDuplication(tmp_up_columns[id_string], up_columns[id_string]);
			}

			SLAE_in_old_DOF.Initialization(up_columns, down_columns);

			//переписываем матрицы
			for (int id_str = 0; id_str < solver_grid.GetDOFsCount(); id_str++)
			{
				for (int ii = 0; ii < 3; ii++)
				{
					int new_id_str = id_str + ii * solver_grid.GetDOFsCount();

					SLAE_in_old_DOF.X[new_id_str] = newSLAE.X[id_str * 3 + ii];
					SLAE_in_old_DOF.F[new_id_str] = newSLAE.F[id_str * 3 + ii];
					//SLAE_in_old_DOF.Diag[new_id_str] = newSLAE.Diag[id_str * 3 + ii];
				}


				//диагональ
				for (int ii = 0; ii < 3; ii++)
				{
					int new_id_str = id_str + ii * solver_grid.GetDOFsCount();
					int xyz_id_str = id_str * 3 + ii;
					for (int jj = 0; jj < 3; jj++)
					{
						int new_id_col = id_str + jj * solver_grid.GetDOFsCount();
						int xyz_id_col = id_str * 3 + jj;
						double val;
						if (!newSLAE.GetValue(xyz_id_str, xyz_id_col, val))
						{
							printf_s("ERROR in translation portrait!\n");
							int aa;
							scanf_s("%d", &aa);
						}
						if (!SLAE_in_old_DOF.SetValue(new_id_str, new_id_col, val))
						{
							printf_s("ERROR in translation portrait!\n");
							int aa;
							scanf_s("%d", &aa);
						}
					}
				}

				//верхний треугольник
				for (int j = 0; j < global_SLAE.id_column_for_A_up[id_str].size(); j++)
				{
					int id_col = global_SLAE.id_column_for_A_up[id_str][j];

					for (int ii = 0; ii < 3; ii++)
					{
						int new_id_str = id_str + ii * solver_grid.GetDOFsCount();
						int xyz_id_str = id_str * 3 + ii;

						for (int jj = 0; jj < 3; jj++)
						{
							int new_id_col = id_col + jj * solver_grid.GetDOFsCount();
							int xyz_id_col = id_col * 3 + jj;
							double val;
							if (!newSLAE.GetValue(xyz_id_str, xyz_id_col, val))
							{
								printf_s("ERROR in translation portrait!\n");
								int aa;
								scanf_s("%d", &aa);
							}
							if (!SLAE_in_old_DOF.SetValue(new_id_str, new_id_col, val))
							{
								printf_s("ERROR in translation portrait!\n");
								int aa;
								scanf_s("%d", &aa);
							}
						}
					}
				}

				//нижний треугольник
				for (int j = 0; j < global_SLAE.id_column_for_A_down[id_str].size(); j++)
				{
					int id_col = global_SLAE.id_column_for_A_down[id_str][j];

					for (int ii = 0; ii < 3; ii++)
					{
						int new_id_str = id_str + ii * solver_grid.GetDOFsCount();
						int xyz_id_str = id_str * 3 + ii;

						for (int jj = 0; jj < 3; jj++)
						{
							int new_id_col = id_col + jj * solver_grid.GetDOFsCount();
							int xyz_id_col = id_col * 3 + jj;
							double val;
							if (!newSLAE.GetValue(xyz_id_str, xyz_id_col, val))
							{
								printf_s("ERROR in translation portrait!\n");
								int aa;
								scanf_s("%d", &aa);
							}
							if (!SLAE_in_old_DOF.SetValue(new_id_str, new_id_col, val))
							{
								printf_s("ERROR in translation portrait!\n");
								int aa;
								scanf_s("%d", &aa);
							}
						}
					}
				}
			}

			//печатаем портреты
			//char name_Tensor[1000];
			//sprintf_s(name_Tensor, "%s/global_SLAE.txt", result_directory);
			//PrintPortrait(name_Tensor, global_SLAE);
//
			//char name_Old[1000];
			//sprintf_s(name_Old, "%s/SLAE_in_old_DOF.txt", result_directory);
			//PrintPortrait(name_Old, SLAE_in_old_DOF);
		}

		clock_t t_before = clock();
		if (is_print_logFile) printf("TIME = %lf sec\n", (t_before - t_after) * 1.0 / CLK_TCK);
		double end = omp_get_wtime();
		if (is_print_logFile) printf("start = %.16g\nend = %.16g\ndiff = %.16g\n", start, end, end - start);

		if (false)
		{
			if (is_print_logFile) printf("Soluting SLAY... (%d)\n", newSLAE.GetMatrixSize());
			int MaxSize = newSLAE.GetMatrixSize();
			MaxSize = MaxSize / 10 < 1000 ? 1000 : MaxSize / 10;
			std::vector<double> best_solution;
			math::MakeCopyVector_A_into_B(newSLAE.X, best_solution);
			double current_residual = 1, best_residual = 1e+25;
			int MAX_STEPS = 4;
			int ii = 1;
			CSSD_Matrix<double, double> Precond;
			Precond.PrecondorSSOR(0.75, newSLAE);
			for (int i = 0; i <= MAX_STEPS; i++)
			{
				double needed_residual = pow(10., -1 * (ii + 1));
				needed_residual = needed_residual < MIN_RESIDUAL ? MIN_RESIDUAL : needed_residual;
				//current_residual /= 2.;
				if(is_print_logFile) printf_s("//---> I = %d/%d (full size %d) - needed residual %.2e\n", i, MAX_STEPS, newSLAE.GetMatrixSize(), needed_residual);

				current_residual = abs(newSLAE.MSG_PreconditioningSSOR(MaxSize,
					needed_residual,
					Precond));

				if (current_residual < needed_residual)
				{
					i = 0;
					ii++;
				}
				if (current_residual > best_residual * 100)
				{
					break;
				}
				if (current_residual <= MIN_RESIDUAL)
				{
					best_residual = current_residual;
					break;
				}
				if (current_residual < best_residual)
				{
					best_residual = current_residual;
					math::MakeCopyVector_A_into_B(newSLAE.X, best_solution);
				}
			}
			math::MakeCopyVector_A_into_B(best_solution, newSLAE.X);
			if (is_print_logFile) printf_s("//---> BEST residual %.2e\n", best_residual);
			//переписываем лучшее решение в векторный вид
			math::MakeCopyVector_A_into_B(newSLAE.X, Solution);

			if (best_residual >= 1)
			{
				printf_s("There are problems in the solution\n");
				int a;
				scanf_s("%d", &a);
			}
		}
		if (true)
		{
			//SLAE_in_old_DOF.print_logs = true;
			printf("Soluting SLAY... (%d)\n", SLAE_in_old_DOF.GetMatrixSize());
			int MaxSize = SLAE_in_old_DOF.GetMatrixSize();
			MaxSize = MaxSize / 10 < 1000 ? 1000 : MaxSize / 10;
			std::vector<double> best_solution;
			math::MakeCopyVector_A_into_B(SLAE_in_old_DOF.X, best_solution);
			double current_residual = 1, best_residual = 1e+25;
			int MAX_STEPS = 4;
			int ii = 1;
			CSSD_Matrix<double, double> Precond;
			Precond.PrecondorSSOR(0.75, SLAE_in_old_DOF);
			SLAE_in_old_DOF.print_logs = is_print_logFile;
			for (int i = 0; i <= MAX_STEPS; i++)
			{
				double needed_residual = pow(10., -1 * (ii + 1));
				needed_residual = needed_residual < MIN_RESIDUAL ? MIN_RESIDUAL : needed_residual;
				//current_residual /= 2.;
				if (is_print_logFile) printf_s("//---> I = %d/%d (full size %d) - needed residual %.2e\n", i, MAX_STEPS, newSLAE.GetMatrixSize(), needed_residual);

				current_residual = abs(SLAE_in_old_DOF.MSG_PreconditioningSSOR(MaxSize,
					needed_residual,
					Precond));

				if (current_residual < needed_residual)
				{
					i = 0;
					ii++;
				}
				else {
					current_residual = abs(SLAE_in_old_DOF.BCG_Stab2(MaxSize, 1e-18));
					current_residual = abs(SLAE_in_old_DOF.MSG(MaxSize, 1e-18));
				}
				if (current_residual > best_residual * 100)
				{
					break;
				}
				if (current_residual <= MIN_RESIDUAL)
				{
					best_residual = current_residual;
					break;
				}
				if (current_residual < best_residual)
				{
					best_residual = current_residual;
					math::MakeCopyVector_A_into_B(SLAE_in_old_DOF.X, best_solution);
				}
			}
			if (current_residual > 1e-15)
				printf_s("//---> BEST residual %.2e\n", best_residual);

			math::MakeCopyVector_A_into_B(best_solution, SLAE_in_old_DOF.X);
			printf_s("//---> BEST residual %.2e\n", best_residual);
			//переписываем лучшее решение в векторный вид
			//math::MakeCopyVector_A_into_B(SLAE_in_old_DOF.X, Solution);
			Solution.resize(solver_grid.GetDOFsCount());
			for (int i = 0; i < Solution.size(); i++)
			{
				Solution[i].x = SLAE_in_old_DOF.X[i + 0 * Solution.size()];
				Solution[i].y = SLAE_in_old_DOF.X[i + 1 * Solution.size()];
				Solution[i].z = SLAE_in_old_DOF.X[i + 2 * Solution.size()];
			}

			if (best_residual >= 1)
			{
				printf_s("There are problems in the solution\n");
				int a;
				scanf_s("%d", &a);
			}
		}

		if (is_print_logFile) printf("\tcomplit\n\n");
		t_before = clock();
		if (is_print_logFile) printf("TIME = %lf sec\n", (t_before - t_after) * 1.0 / CLK_TCK);
		end = omp_get_wtime();
		if (is_print_logFile) printf("start = %.16g\nend = %.16g\ndiff = %.16g\n", start, end, end - start);



		//update grid
		/*{
			std::vector<Point<double>> new_x(solver_grid.GetVertexCount());
			for (int id_vertex = 0; id_vertex < solver_grid.GetVertexCount(); id_vertex++)
			{
				new_x[id_vertex] = solver_grid.GetCoordinateViaID(id_vertex) + global_SLAE.X[id_vertex];
			}
			solver_grid.UpdateCoordinates(new_x);
		}*/


		if (is_print_logFile) printf_s("complite\n\n");
	};

	template <typename Dirichlet, typename Neumann>
	void FEM_forElasticDeformation(
		bool is_print_logFile, //input
		double MIN_RESIDUAL,
		math::SimpleGrid &geo_grid, //input
		std::vector<Dirichlet> &first_boundary,//input
		std::vector<Neumann> &second_boundary,//input
		char* result_directory, //output
		FEM::Grid_forMech &solver_grid, //output
		std::vector<Point<double>> &Solution, //output
		CSSD_Matrix<Tensor2Rank3D, Point<double>> &global_SLAE //output
	)
	{
		FILE *stream;
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
		solver_grid.Initialization(geo_grid, first_boundary, second_boundary);
		printf_s("complite                   \n\n");

		

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
				if (id_elem % 1000 == 0)
					printf("Solve element[%d]\r", id_elem);
				auto element = solver_grid.GetElement(id_elem);
				std::function< std::vector<std::vector<double>>(Point<double> X) > koefD = [&](Point<double> X)->std::vector<std::vector<double>>
				{
					return (solver_grid.GetDomain(0))->forMech.GetD(3);
				};
				element->SolveLocalMatrix(local_SLAE[id_elem], koefD);

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
		{
			//global_SLAE.MakeDivideOnSqrtFormSumm();

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

				//if (is_take.x == true)
				//{
				//	global_SLAE.X[global_id].x = boundary_value.x;
				//	global_SLAE.F[global_id].x = boundary_value.x;
				//	global_SLAE.Diag[global_id].val[0][0] = 1.0;
				//	global_SLAE.Diag[global_id].val[0][1] = 0.0;
				//	global_SLAE.Diag[global_id].val[0][2] = 0.0;
				//	for (int j = 0; j < global_SLAE.A_down[global_id].size(); j++)
				//	{
				//		global_SLAE.A_down[global_id][j].val[0][0] = 0.0;
				//		global_SLAE.A_down[global_id][j].val[0][1] = 0.0;
				//		global_SLAE.A_down[global_id][j].val[0][2] = 0.0;
				//	}
				//	for (int j = 0; j < global_SLAE.A_up[global_id].size(); j++)
				//	{
				//		global_SLAE.A_up[global_id][j].val[0][0] = 0.0;
				//		global_SLAE.A_up[global_id][j].val[0][1] = 0.0;
				//		global_SLAE.A_up[global_id][j].val[0][2] = 0.0;
				//	}
				//	//simmetrisation
				//	int position = 0;
				//	for (int i = 0; i < global_id; i++)
				//	{
				//		for (int j = 0; j < global_SLAE.id_column_for_A_up[i].size() &&
				//			global_SLAE.id_column_for_A_up[i][j] <= global_id; j++)
				//		{
				//			if (global_id == global_SLAE.id_column_for_A_up[i][j])
				//			{
				//				global_SLAE.F[i].x -= global_SLAE.F[global_id].x *global_SLAE.A_up[i][j].val[0][position];
				//				global_SLAE.F[i].y -= global_SLAE.F[global_id].y *global_SLAE.A_up[i][j].val[1][position];
				//				global_SLAE.F[i].z -= global_SLAE.F[global_id].z *global_SLAE.A_up[i][j].val[2][position];
				//				global_SLAE.A_up[i][j].val[0][position] = 0.0;
				//				global_SLAE.A_up[i][j].val[1][position] = 0.0;
				//				global_SLAE.A_up[i][j].val[2][position] = 0.0;
				//			}
				//		}
				//	}
				//	for (int i = global_id+1; i < global_SLAE.id_column_for_A_down.size(); i++)
				//	{
				//		for (int j = 0; j < global_SLAE.id_column_for_A_down[i].size() &&
				//			global_SLAE.id_column_for_A_down[i][j] <= global_id; j++)
				//		{
				//			if (global_id == global_SLAE.id_column_for_A_down[i][j])
				//			{
				//				global_SLAE.F[i].x -= global_SLAE.F[global_id].x *global_SLAE.A_down[i][j].val[0][position];
				//				global_SLAE.F[i].y -= global_SLAE.F[global_id].y *global_SLAE.A_down[i][j].val[1][position];
				//				global_SLAE.F[i].z -= global_SLAE.F[global_id].z *global_SLAE.A_down[i][j].val[2][position];
				//				global_SLAE.A_down[i][j].val[0][position] = 0.0;
				//				global_SLAE.A_down[i][j].val[1][position] = 0.0;
				//				global_SLAE.A_down[i][j].val[2][position] = 0.0;
				//			}
				//		}
				//	}
				//}
				//if (is_take.y == true)
				//{
				//	global_SLAE.X[global_id].y = boundary_value.y;
				//	global_SLAE.F[global_id].y = boundary_value.y;
				//	global_SLAE.Diag[global_id].val[1][0] = 0.0;
				//	global_SLAE.Diag[global_id].val[1][1] = 1.0;
				//	global_SLAE.Diag[global_id].val[1][2] = 0.0;
				//	for (int j = 0; j < global_SLAE.A_down[global_id].size(); j++)
				//	{
				//		global_SLAE.A_down[global_id][j].val[1][0] = 0.0;
				//		global_SLAE.A_down[global_id][j].val[1][1] = 0.0;
				//		global_SLAE.A_down[global_id][j].val[1][2] = 0.0;
				//	}
				//	for (int j = 0; j < global_SLAE.A_up[global_id].size(); j++)
				//	{
				//		global_SLAE.A_up[global_id][j].val[1][0] = 0.0;
				//		global_SLAE.A_up[global_id][j].val[1][1] = 0.0;
				//		global_SLAE.A_up[global_id][j].val[1][2] = 0.0;
				//	}
				//}
				//if (is_take.z == true)
				//{
				//	global_SLAE.X[global_id].z = boundary_value.z;
				//	global_SLAE.F[global_id].z = boundary_value.z;
				//	global_SLAE.Diag[global_id].val[2][0] = 0.0;
				//	global_SLAE.Diag[global_id].val[2][1] = 0.0;
				//	global_SLAE.Diag[global_id].val[2][2] = 1.0;
				//	for (int j = 0; j < global_SLAE.A_down[global_id].size(); j++)
				//	{
				//		global_SLAE.A_down[global_id][j].val[2][0] = 0.0;
				//		global_SLAE.A_down[global_id][j].val[2][1] = 0.0;
				//		global_SLAE.A_down[global_id][j].val[2][2] = 0.0;
				//	}
				//	for (int j = 0; j < global_SLAE.A_up[global_id].size(); j++)
				//	{
				//		global_SLAE.A_up[global_id][j].val[2][0] = 0.0;
				//		global_SLAE.A_up[global_id][j].val[2][1] = 0.0;
				//		global_SLAE.A_up[global_id][j].val[2][2] = 0.0;
				//	}
				//}

				/*printf_s("complite\n\n");
				char name[1000];
				sprintf_s(name, sizeof(name), "%s/Matrix_with_boundary_id_%d.txt", result_directory, global_id);
				PrintMatrix(name, global_SLAE);*/
			}

		}

		clock_t t_before = clock();
		printf("TIME = %lf sec\n", (t_before - t_after) * 1.0 / CLK_TCK);
		double end = omp_get_wtime();
		printf("start = %.16g\nend = %.16g\ndiff = %.16g\n", start, end, end - start);

		//SLAE solution
		//{
		//	printf("Soluting SLAY... (%d)\n", global_SLAE.GetMatrixSize());
		//	int MaxSize = global_SLAE.GetMatrixSize();
		//	int N_matrix = MaxSize / 10000 + 100;
		//	for (int i = 0; i <= N_matrix; i++)
		//	{
		//		printf_s("//---> I = %d/%d (full size %d)\n", i, N_matrix - 1, global_SLAE.GetMatrixSize());
		//		double residual = global_SLAE.BiCG_Stab(MaxSize /*/ N_matrix*/, MIN_RESIDUAL);
		//		if (residual <= MIN_RESIDUAL)
		//			break;
		//	}
		//	math::MakeCopyVector_A_into_B(global_SLAE.X, Solution);
		//}
		if (false)
		{
			printf("Soluting SLAY... (%d)\n", global_SLAE.GetMatrixSize());
			int MaxSize = global_SLAE.GetMatrixSize() / 100;
			std::vector<Point<double>> best_solution;
			math::MakeCopyVector_A_into_B(global_SLAE.X, best_solution);
			double current_residual = 1, best_residual = 1e+25;
			int MAX_STEPS = 12;
			for (int i = 0; i <= MAX_STEPS; i++)
			{
				current_residual = pow(10., -1 * (i + 1));
				//current_residual /= 2.;
				printf_s("//---> I = %d/%d (full size %d) - needed residual %.2e\n", i, MAX_STEPS, global_SLAE.GetMatrixSize(), current_residual);
				current_residual = abs(global_SLAE.BiCG_Stab(MaxSize, current_residual));
				if (current_residual <= MIN_RESIDUAL)
					break;
				if (current_residual > best_residual + (1e-10))
				{
					math::MakeCopyVector_A_into_B(best_solution, global_SLAE.X);
					printf_s("//---> BEST residual %.2e\n", best_residual);
					break;
				}
				if (current_residual < best_residual)
				{
					best_residual = current_residual;
					math::MakeCopyVector_A_into_B(global_SLAE.X, best_solution);
				}
			}
			math::MakeCopyVector_A_into_B(global_SLAE.X, Solution);
		}
		if (true) {
			CSSD_Matrix<double, double> newSLAE;
			char file_matrix[1000];
			sprintf_s(file_matrix, sizeof(file_matrix), "%s/matrix_portrait.txt", result_directory);
			//FEM::PrintMatrix(file_matrix, global_SLAE);
			math::MakeCopyMatrix_A_into_B(global_SLAE, newSLAE);
			sprintf_s(file_matrix, sizeof(file_matrix), "%s/matrix_portrait_NEW.txt", result_directory);

			newSLAE.MakeDivideOnSqrtFormSumm();
			//FEM::PrintMatrix(file_matrix, newSLAE);

			printf("Soluting SLAY... (%d)\n", newSLAE.GetMatrixSize());
			int MaxSize = newSLAE.GetMatrixSize();
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
				current_residual = abs(newSLAE.BiCG_Stab(MaxSize, needed_residual));
				//current_residual = abs(newSLAE.MSG(MaxSize, needed_residual));
				if (current_residual <= MIN_RESIDUAL)
					break;
				if (current_residual > best_residual + (1e-10))
				{
					//math::MakeCopyVector_A_into_B(best_solution, newSLAE.X);
					//printf_s("//---> BEST residual %.2e\n", best_residual);
					//break;
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
	template <typename Dirichlet, typename Neumann>
	void FEM_forElasticDeformation_OrderBF2(
		bool is_print_logFile, //input
		double MIN_RESIDUAL,
		math::SimpleGrid &geo_grid, //input
		std::vector<Dirichlet> &first_boundary,//input
		std::vector<Neumann> &second_boundary,//input
		char* result_directory, //output
		FEM::Grid_forMech &solver_grid, //output
		std::vector<Point<double>> &Solution, //output
		CSSD_Matrix<Tensor2Rank3D, Point<double>> &global_SLAE //output
	)
	{
		FILE *stream;
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
		//creation basis functions
		/*std::vector< std::vector<std::function< Point<double>(Point<double> X)>>> bf;
		std::vector < std::vector<std::function< Point<Point<double>>(Point<double> X)>>> derivative_bf;

		math::ResizeVector(bf, geo_grid.nvtr.size(), 4+6);
		auto element = this->GetElement(id_element);
		element->ResizeDOF(element->GetNodesCount());
		element->SolveAlphaMatrix();

		for (int _id_vertex = 0; _id_vertex < element->GetNodesCount(); _id_vertex++)
		{
			std::function< Point<double>(Point<double> X)> bf = [element, _id_vertex](Point<double> X) -> Point<double>
			{
				Point<double> result;
				result.x = element->alpha[_id_vertex][0] + element->alpha[_id_vertex][1] * X.x + element->alpha[_id_vertex][2] * X.y + element->alpha[_id_vertex][3] * X.z;
				result.y = result.x;
				result.z = result.x;
				return result;
			};
			std::function< Point<Point<double>>(Point<double> X)> derivative_bf = [element, _id_vertex](Point<double> X) -> Point<Point<double>>
			{
				Point<Point<double>> result;

				result.x.x = element->alpha[_id_vertex][1];
				result.x.y = element->alpha[_id_vertex][2];
				result.x.z = element->alpha[_id_vertex][3];
				result.y.x = result.x.x;
				result.y.y = result.x.y;
				result.y.z = result.x.z;
				result.z.x = result.x.x;
				result.z.y = result.x.y;
				result.z.z = result.x.z;

				return result;
			};
			element->SetDOF(_id_vertex, element->GetIdNode(_id_vertex), bf, derivative_bf);
		}
		solver_grid.Initialization(geo_grid, first_boundary, second_boundary);*/
		printf_s("complite                   \n\n");



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
				if (id_elem % 1000 == 0)
					printf("Solve element[%d]\r", id_elem);
				auto element = solver_grid.GetElement(id_elem);
				std::function< std::vector<std::vector<double>>(Point<double> X) > koefD = [&](Point<double> X)->std::vector<std::vector<double>>
				{
					return (solver_grid.GetDomain(0))->forMech.GetD(3);
				};
				element->SolveLocalMatrix(local_SLAE[id_elem], koefD);

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
		{
			//global_SLAE.MakeDivideOnSqrtFormSumm();

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

				//if (is_take.x == true)
				//{
				//	global_SLAE.X[global_id].x = boundary_value.x;
				//	global_SLAE.F[global_id].x = boundary_value.x;
				//	global_SLAE.Diag[global_id].val[0][0] = 1.0;
				//	global_SLAE.Diag[global_id].val[0][1] = 0.0;
				//	global_SLAE.Diag[global_id].val[0][2] = 0.0;
				//	for (int j = 0; j < global_SLAE.A_down[global_id].size(); j++)
				//	{
				//		global_SLAE.A_down[global_id][j].val[0][0] = 0.0;
				//		global_SLAE.A_down[global_id][j].val[0][1] = 0.0;
				//		global_SLAE.A_down[global_id][j].val[0][2] = 0.0;
				//	}
				//	for (int j = 0; j < global_SLAE.A_up[global_id].size(); j++)
				//	{
				//		global_SLAE.A_up[global_id][j].val[0][0] = 0.0;
				//		global_SLAE.A_up[global_id][j].val[0][1] = 0.0;
				//		global_SLAE.A_up[global_id][j].val[0][2] = 0.0;
				//	}
				//	//simmetrisation
				//	int position = 0;
				//	for (int i = 0; i < global_id; i++)
				//	{
				//		for (int j = 0; j < global_SLAE.id_column_for_A_up[i].size() &&
				//			global_SLAE.id_column_for_A_up[i][j] <= global_id; j++)
				//		{
				//			if (global_id == global_SLAE.id_column_for_A_up[i][j])
				//			{
				//				global_SLAE.F[i].x -= global_SLAE.F[global_id].x *global_SLAE.A_up[i][j].val[0][position];
				//				global_SLAE.F[i].y -= global_SLAE.F[global_id].y *global_SLAE.A_up[i][j].val[1][position];
				//				global_SLAE.F[i].z -= global_SLAE.F[global_id].z *global_SLAE.A_up[i][j].val[2][position];
				//				global_SLAE.A_up[i][j].val[0][position] = 0.0;
				//				global_SLAE.A_up[i][j].val[1][position] = 0.0;
				//				global_SLAE.A_up[i][j].val[2][position] = 0.0;
				//			}
				//		}
				//	}
				//	for (int i = global_id+1; i < global_SLAE.id_column_for_A_down.size(); i++)
				//	{
				//		for (int j = 0; j < global_SLAE.id_column_for_A_down[i].size() &&
				//			global_SLAE.id_column_for_A_down[i][j] <= global_id; j++)
				//		{
				//			if (global_id == global_SLAE.id_column_for_A_down[i][j])
				//			{
				//				global_SLAE.F[i].x -= global_SLAE.F[global_id].x *global_SLAE.A_down[i][j].val[0][position];
				//				global_SLAE.F[i].y -= global_SLAE.F[global_id].y *global_SLAE.A_down[i][j].val[1][position];
				//				global_SLAE.F[i].z -= global_SLAE.F[global_id].z *global_SLAE.A_down[i][j].val[2][position];
				//				global_SLAE.A_down[i][j].val[0][position] = 0.0;
				//				global_SLAE.A_down[i][j].val[1][position] = 0.0;
				//				global_SLAE.A_down[i][j].val[2][position] = 0.0;
				//			}
				//		}
				//	}
				//}
				//if (is_take.y == true)
				//{
				//	global_SLAE.X[global_id].y = boundary_value.y;
				//	global_SLAE.F[global_id].y = boundary_value.y;
				//	global_SLAE.Diag[global_id].val[1][0] = 0.0;
				//	global_SLAE.Diag[global_id].val[1][1] = 1.0;
				//	global_SLAE.Diag[global_id].val[1][2] = 0.0;
				//	for (int j = 0; j < global_SLAE.A_down[global_id].size(); j++)
				//	{
				//		global_SLAE.A_down[global_id][j].val[1][0] = 0.0;
				//		global_SLAE.A_down[global_id][j].val[1][1] = 0.0;
				//		global_SLAE.A_down[global_id][j].val[1][2] = 0.0;
				//	}
				//	for (int j = 0; j < global_SLAE.A_up[global_id].size(); j++)
				//	{
				//		global_SLAE.A_up[global_id][j].val[1][0] = 0.0;
				//		global_SLAE.A_up[global_id][j].val[1][1] = 0.0;
				//		global_SLAE.A_up[global_id][j].val[1][2] = 0.0;
				//	}
				//}
				//if (is_take.z == true)
				//{
				//	global_SLAE.X[global_id].z = boundary_value.z;
				//	global_SLAE.F[global_id].z = boundary_value.z;
				//	global_SLAE.Diag[global_id].val[2][0] = 0.0;
				//	global_SLAE.Diag[global_id].val[2][1] = 0.0;
				//	global_SLAE.Diag[global_id].val[2][2] = 1.0;
				//	for (int j = 0; j < global_SLAE.A_down[global_id].size(); j++)
				//	{
				//		global_SLAE.A_down[global_id][j].val[2][0] = 0.0;
				//		global_SLAE.A_down[global_id][j].val[2][1] = 0.0;
				//		global_SLAE.A_down[global_id][j].val[2][2] = 0.0;
				//	}
				//	for (int j = 0; j < global_SLAE.A_up[global_id].size(); j++)
				//	{
				//		global_SLAE.A_up[global_id][j].val[2][0] = 0.0;
				//		global_SLAE.A_up[global_id][j].val[2][1] = 0.0;
				//		global_SLAE.A_up[global_id][j].val[2][2] = 0.0;
				//	}
				//}

				/*printf_s("complite\n\n");
				char name[1000];
				sprintf_s(name, sizeof(name), "%s/Matrix_with_boundary_id_%d.txt", result_directory, global_id);
				PrintMatrix(name, global_SLAE);*/
			}

		}

		clock_t t_before = clock();
		printf("TIME = %lf sec\n", (t_before - t_after) * 1.0 / CLK_TCK);
		double end = omp_get_wtime();
		printf("start = %.16g\nend = %.16g\ndiff = %.16g\n", start, end, end - start);

		//SLAE solution
		//{
		//	printf("Soluting SLAY... (%d)\n", global_SLAE.GetMatrixSize());
		//	int MaxSize = global_SLAE.GetMatrixSize();
		//	int N_matrix = MaxSize / 10000 + 100;
		//	for (int i = 0; i <= N_matrix; i++)
		//	{
		//		printf_s("//---> I = %d/%d (full size %d)\n", i, N_matrix - 1, global_SLAE.GetMatrixSize());
		//		double residual = global_SLAE.BiCG_Stab(MaxSize /*/ N_matrix*/, MIN_RESIDUAL);
		//		if (residual <= MIN_RESIDUAL)
		//			break;
		//	}
		//	math::MakeCopyVector_A_into_B(global_SLAE.X, Solution);
		//}
		if (false)
		{
			printf("Soluting SLAY... (%d)\n", global_SLAE.GetMatrixSize());
			int MaxSize = global_SLAE.GetMatrixSize() / 100;
			std::vector<Point<double>> best_solution;
			math::MakeCopyVector_A_into_B(global_SLAE.X, best_solution);
			double current_residual = 1, best_residual = 1e+25;
			int MAX_STEPS = 12;
			for (int i = 0; i <= MAX_STEPS; i++)
			{
				current_residual = pow(10., -1 * (i + 1));
				//current_residual /= 2.;
				printf_s("//---> I = %d/%d (full size %d) - needed residual %.2e\n", i, MAX_STEPS, global_SLAE.GetMatrixSize(), current_residual);
				current_residual = abs(global_SLAE.BiCG_Stab(MaxSize, current_residual));
				if (current_residual <= MIN_RESIDUAL)
					break;
				if (current_residual > best_residual + (1e-10))
				{
					math::MakeCopyVector_A_into_B(best_solution, global_SLAE.X);
					printf_s("//---> BEST residual %.2e\n", best_residual);
					break;
				}
				if (current_residual < best_residual)
				{
					best_residual = current_residual;
					math::MakeCopyVector_A_into_B(global_SLAE.X, best_solution);
				}
			}
			math::MakeCopyVector_A_into_B(global_SLAE.X, Solution);
		}
		if (true) {
			CSSD_Matrix<double, double> newSLAE;
			char file_matrix[1000];
			sprintf_s(file_matrix, sizeof(file_matrix), "%s/matrix_portrait.txt", result_directory);
			//FEM::PrintMatrix(file_matrix, global_SLAE);
			math::MakeCopyMatrix_A_into_B(global_SLAE, newSLAE);
			sprintf_s(file_matrix, sizeof(file_matrix), "%s/matrix_portrait_NEW.txt", result_directory);

			newSLAE.MakeDivideOnSqrtFormSumm();
			//FEM::PrintMatrix(file_matrix, newSLAE);

			printf("Soluting SLAY... (%d)\n", newSLAE.GetMatrixSize());
			int MaxSize = newSLAE.GetMatrixSize();
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
				current_residual = abs(newSLAE.BiCG_Stab(MaxSize, needed_residual));
				//current_residual = abs(newSLAE.MSG(MaxSize, needed_residual));
				if (current_residual <= MIN_RESIDUAL)
					break;
				if (current_residual > best_residual + (1e-10))
				{
					//math::MakeCopyVector_A_into_B(best_solution, newSLAE.X);
					//printf_s("//---> BEST residual %.2e\n", best_residual);
					//break;
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