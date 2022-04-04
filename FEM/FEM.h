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
		if (is_print_logFile == false)
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

		if(false)
		{
			printf("Soluting SLAY... (%d)\n", global_SLAE.GetMatrixSize());
			int MaxSize = global_SLAE.GetMatrixSize()/100;
			std::vector<Point<double>> best_solution;
			math::MakeCopyVector_A_into_B(global_SLAE.X, best_solution);
			double current_residual = 1, best_residual = 1e+25;
			int MAX_STEPS = 15;
			for (int i = 0; i <= MAX_STEPS; i++)
			{
				current_residual = pow(10., -1 * (i+1));
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

			newSLAE.log_out = stream;
			newSLAE.print_logs = is_print_logFile;

			printf("Soluting SLAY... (%d)\n", newSLAE.GetMatrixSize());
			int MaxSize = newSLAE.GetMatrixSize();
			MaxSize = MaxSize/10 < 10 ? 10 : MaxSize/10;
			std::vector<double> best_solution;
			math::MakeCopyVector_A_into_B(newSLAE.X, best_solution);
			double current_residual = 1, best_residual = 1e+25;
			int MAX_STEPS = 15;
			int ii = 1;
			//newSLAE.MakeDivideOnSqrtFormSumm();
			for (int i = 0; i <= MAX_STEPS; i++)
			{
				double needed_residual = pow(10., -1 * (ii + 1));
				//current_residual /= 2.;
				printf_s("//---> I = %d/%d (full size %d) - needed residual %.2e\n", i, MAX_STEPS, newSLAE.GetMatrixSize(), needed_residual);
				//current_residual = newSLAE.GMRES(5, needed_residual);
				
				//current_residual = abs(newSLAE.BiCG_Stab(MaxSize, needed_residual));
				current_residual = abs(newSLAE.MSG(MaxSize, needed_residual));

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

			//if (best_residual >= 1e-8)
			//{
			//	MaxSize = 100;
			//	MaxSize = MaxSize < 10 ? 10 : MaxSize;
			//	current_residual = 1;
			//	int MAX_STEPS = 11;
			//	double needed_residual = best_residual / 10;
			//	for (int i = 0; i <= MAX_STEPS; i++)
			//	{
			//		needed_residual /= 10;
			//		printf_s("//---> I = %d/%d (full size %d) - needed residual %.2e\n", i, MAX_STEPS, newSLAE.GetMatrixSize(), needed_residual);
			//		current_residual = abs(newSLAE.MSG(MaxSize, needed_residual));
			//		if (current_residual <= MIN_RESIDUAL)
			//			break;
			//		if (current_residual > best_residual + (1e-10))
			//		{
			//			//math::MakeCopyVector_A_into_B(best_solution, newSLAE.X);
			//			//printf_s("//---> BEST residual %.2e\n", best_residual);
			//			break;
			//		}
			//		if (current_residual < best_residual)
			//		{
			//			best_residual = current_residual;
			//			math::MakeCopyVector_A_into_B(newSLAE.X, best_solution);
			//		}
			//	}
			//	math::MakeCopyVector_A_into_B(best_solution, newSLAE.X);
			//	printf_s("//---> BEST residual %.2e\n", best_residual);
			//	math::MakeCopyVector_A_into_B(newSLAE.X, Solution);
			//}

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
				current_residual = abs(newSLAE.MSG_Preconditioning(MaxSize, needed_residual, Predcondor));
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
			MaxSize = MaxSize < 10 ? 10 : MaxSize;
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

			if (best_residual >= 1)
			{
				printf_s("There are problems in the solution\n");
				int a;
				scanf_s("%d", &a);
			}
		}
		if (true) 
		{
			math::MakeCopyVector_A_into_B(newSLAE.F, newSLAE.X);

			printf("Soluting SLAY... (%d)\n", newSLAE.GetMatrixSize());
			int MaxSize = newSLAE.GetMatrixSize();
			MaxSize = MaxSize < 10 ? 10 : MaxSize;
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

			if (best_residual >= 1e-8)
			{
				MaxSize *= 100;
				MaxSize = MaxSize < 10 ? 10 : MaxSize;
				current_residual = 1;
				int MAX_STEPS = 15;
				double needed_residual = best_residual / 10;
				for (int i = 0; i <= MAX_STEPS; i++)
				{
					needed_residual /= 10;
					printf_s("//---> I = %d/%d (full size %d) - needed residual %.2e\n", i, MAX_STEPS, newSLAE.GetMatrixSize(), needed_residual);
					current_residual = abs(newSLAE.MSG(MaxSize, needed_residual));
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

		//add right side
		{
			printf("SLAE assembling (right side...\n");
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
			printf_s("                                                                                    \r");
			printf_s("complite\n\n");
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
		}

		clock_t t_before = clock();
		printf("TIME = %lf sec\n", (t_before - t_after) * 1.0 / CLK_TCK);
		double end = omp_get_wtime();
		printf("start = %.16g\nend = %.16g\ndiff = %.16g\n", start, end, end - start);

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

			if (best_residual >= 1)
			{
				printf_s("There are problems in the solution\n");
				char a[23];
				scanf_s("%s", &a);
			}
		}
		if (true) {
			math::MakeCopyVector_A_into_B(newSLAE.F, newSLAE.X);

			printf("Soluting SLAY... (%d)\n", newSLAE.GetMatrixSize());
			int MaxSize = newSLAE.GetMatrixSize();
			MaxSize = MaxSize < 10 ? 10 : MaxSize;
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

			if (best_residual >= 1e-8)
			{
				MaxSize *= 100;
				MaxSize = MaxSize < 10 ? 10 : MaxSize;
				current_residual = 1;
				int MAX_STEPS = 15;
				double needed_residual = best_residual / 10;
				for (int i = 0; i <= MAX_STEPS; i++)
				{
					needed_residual /= 10;
					printf_s("//---> I = %d/%d (full size %d) - needed residual %.2e\n", i, MAX_STEPS, newSLAE.GetMatrixSize(), needed_residual);
					current_residual = abs(newSLAE.MSG(MaxSize, needed_residual));
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
			printf("SLAE assembling (right side...\n");
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
			for (int id_elem = 0; id_elem < solver_grid.GetElementsCount(); id_elem++)
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
			printf_s("                                                                                    \r");
			printf_s("complite\n\n");
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
			MaxSize = MaxSize < 10 ? 10 : MaxSize;
			std::vector<Point<double>> best_solution;
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

			if (best_residual >= 1)
			{
				printf_s("There are problems in the solution\n");
				char a[23];
				scanf_s("%s", &a);
			}
		}
		if (true) 
		{
			math::MakeCopyVector_A_into_B(newSLAE.F, newSLAE.X);

			printf("Soluting SLAY... (%d)\n", newSLAE.GetMatrixSize());
			int MaxSize = newSLAE.GetMatrixSize();
			MaxSize = MaxSize < 10 ? 10 : MaxSize;
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

			if (best_residual >= 1e-8)
			{
				MaxSize *= 100;
				MaxSize = MaxSize < 10 ? 10 : MaxSize;
				current_residual = 1;
				int MAX_STEPS = 15;
				double needed_residual = best_residual / 10;
				for (int i = 0; i <= MAX_STEPS; i++)
				{
					needed_residual /= 10;
					printf_s("//---> I = %d/%d (full size %d) - needed residual %.2e\n", i, MAX_STEPS, newSLAE.GetMatrixSize(), needed_residual);
					current_residual = abs(newSLAE.MSG(MaxSize, needed_residual));
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