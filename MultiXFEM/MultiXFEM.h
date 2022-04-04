#pragma once
#include "../Library/GeometryGrid.h"
#include "../Library/Point.h"
#include <functional>
#include "../Library/GeometryShape.h"
#include "MultiXFEM_Grid.h"
#include <stdio.h>
#include <time.h>
#include "omp.h"
#include "../Library/CSSD_Matrix.h"
#include<complex>

namespace MultiXFEM {
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
					dense_matrix[id_string * 3 + ii][id_string * 3 + jj] = 1;// matrix.Diag[id_string].val[ii][jj];
				}
			}

			for (int j = 0; j < matrix.id_column_for_A_up[id_string].size(); j++)
			{
				int id_row = matrix.id_column_for_A_up[id_string][j];
				for (int ii = 0; ii < 3; ii++)
				{
					for (int jj = 0; jj < 3; jj++)
					{
						dense_matrix[id_string * 3 + ii][id_row * 3 + jj] = 6;// matrix.A_up[id_string][j].val[ii][jj];
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
						dense_matrix[id_string * 3 + ii][id_row * 3 + jj] = 6;// matrix.A_down[id_string][j].val[ii][jj];
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
	
	void CrackPropagation_3D(MultiXFEM::Grid &grid, std::vector<Point<double>> &U, char *dir, int curr_STEP)
	{
		//to find the propagation vector for each a front point
		for (int Cr = 0; Cr < grid.cracks.size(); Cr++)
		{
			double R_for_integration = grid.cracks[Cr].R;
			for (int Fr = 0; Fr < grid.cracks[Cr].id_points_in_fronts.size(); Fr++)
			{
				grid.cracks[Cr].propogation.resize(grid.cracks[Cr].id_points_in_fronts.size());
				grid.cracks[Cr].propogation[Fr].resize(grid.cracks[Cr].id_points_in_fronts[Fr].size() - 2);
				for (int segm = 1; segm < grid.cracks[Cr].id_points_in_fronts[Fr].size() - 2; segm++)
				{
					Point<double> P0, P1, P2, P3;
					P0 = grid.cracks[Cr].xyz[grid.cracks[Cr].id_points_in_fronts[Fr][segm - 1]];
					P1 = grid.cracks[Cr].xyz[grid.cracks[Cr].id_points_in_fronts[Fr][segm]];
					P2 = grid.cracks[Cr].xyz[grid.cracks[Cr].id_points_in_fronts[Fr][segm + 1]];
					P3 = grid.cracks[Cr].xyz[grid.cracks[Cr].id_points_in_fronts[Fr][segm + 2]];

					if (segm - 1 == 0 && grid.cracks[Cr].id_points_in_fronts[Fr][segm - 1] == grid.cracks[Cr].id_points_in_fronts[Fr][segm])
					{
						P0 = P1 + (P1 - P2);
					}
					if ((segm + 2) == (grid.cracks[Cr].id_points_in_fronts[Fr].size() - 1) &&
						grid.cracks[Cr].id_points_in_fronts[Fr][segm + 1] == grid.cracks[Cr].id_points_in_fronts[Fr][segm + 2])
					{
						P3 = P2 + (P2 - P1);
					}

					Point<double> r3 = P1 * 3 - P2 * 3 + P3 + P0;
					Point<double> r2 = P1 * (-5) + P2 * 4 - P3;
					Point<double> r1 = P0 * (-1) + P2;
					Point<double> r0 = P1 * 2;
					r3 = P0 * (-1) + P1 * 3 - P2 * 3 + P3;
					r2 = P0 * 2 + P1 * (-5) + P2 * 4 - P3;
					r1 = P0 * (-1) + P2;
					r0 = P1 * 2;
					auto F = [&](double t) -> Point<double>
					{
						return (r3*pow(t, 3) + r2 * pow(t, 2) + r1 * t + r0) / 2.;
					};

					Point<double> _r2 = P1 * 9 - P2 * 9 + P3 * 3 + P0 * 3;
					Point<double> _r1 = P1 * (-10) + P2 * 8 - P3 * 2;
					Point<double> _r0 = P0 * (-1) + P2;
					_r2 = r3 * 3.;
					_r1 = r2 * 2.;
					_r0 = r1;
					auto dF_dt = [&](double t) -> Point<double>
					{
						return (_r2*pow(t, 2) + _r1 * t + _r0) / 2.;
					};

					struct newCoordSystem {
						std::vector<Point<double>> _X; //new basis in OLD coord
						Point<double> O; //New coord centr in OLD coord
						Point<double> _O; //New coord centr in NEW coord
						Tensor2Rank3D A; //transfer from Old into New
						Tensor2Rank3D _A; //transfer from New into Old
						Point<double> Offset; //offset from Old into New
						Point<double> _Offset; //offset from New into Old 

						Point<double> transfer_from_OLD_into_NEW(Point<double> &T)
						{
							return math::MakeTransferPoint(T, A, Offset);
						}
						Point<double> transfer_from_NEW_into_OLD(Point<double> &_T)
						{
							return math::MakeTransferPoint(_T, _A, _Offset);
						}
						Tensor2Rank3D transfer_from_OLD_into_NEW(Tensor2Rank3D &T)
						{
							return math::MakeTransferTensor(T, A, Offset);
						}
					};
					auto NewSystem = [Cr, Fr, &grid, &dF_dt, &F](double t, newCoordSystem &System) -> void
					{
						System._X.resize(3);
						int curr_segm;
						System._X[2] = dF_dt(t);
						curr_segm = 0;

						Point<double> Line[2];
						Line[0] = F(t);
						Line[1] = F(t) + System._X[2];
						Point<double> test;
						for (int tt = 0; tt < 3; tt++)
						{
							if (grid.cracks[Cr].fronts[Fr][curr_segm].id_left !=
								grid.cracks[Cr].triangles[grid.cracks[Cr].fronts[Fr][curr_segm].id_base_triangles].GetIdNode(tt)
								&& grid.cracks[Cr].fronts[Fr][curr_segm].id_right !=
								grid.cracks[Cr].triangles[grid.cracks[Cr].fronts[Fr][curr_segm].id_base_triangles].GetIdNode(tt))
							{
								test = grid.cracks[Cr].xyz[grid.cracks[Cr].triangles[grid.cracks[Cr].fronts[Fr][curr_segm].id_base_triangles].GetIdNode(tt)];
								break;
							}
						}
						DenseMatrix<double, double> test_matrix;
						System._X[0] = math::MakeNormalForLine(Line, test, test_matrix);
						System._X[1] = grid.cracks[Cr].triangles[grid.cracks[Cr].fronts[Fr][curr_segm].id_base_triangles].GetNormal();

						System._X[0] /= math::SolveLengthVector(System._X[0]);
						System._X[1] /= math::SolveLengthVector(System._X[1]);
						System._X[2] /= math::SolveLengthVector(System._X[2]);

						System.A.val[0][0] = System._X[0].x; System.A.val[0][1] = System._X[0].y; System.A.val[0][2] = System._X[0].z;
						System.A.val[1][0] = System._X[1].x; System.A.val[1][1] = System._X[1].y; System.A.val[1][2] = System._X[1].z;
						System.A.val[2][0] = System._X[2].x; System.A.val[2][1] = System._X[2].y; System.A.val[2][2] = System._X[2].z;
						System._A = math::SolveInverseMatrix3x3(System.A);

						Point<double> _p;
						System.Offset = F(t)*(-1.0);
						System._Offset = math::MakeTransferPoint(_p, System.A, System.Offset)*(-1);

						Point<double> _f;
						_f = F(t);
						System.O = F(t);
						System._O = System.transfer_from_OLD_into_NEW(_f);
						System.O = System.transfer_from_NEW_into_OLD(System._O);
					};

					//integration in IxJxK via NEW СOORDINATES
					std::vector<double> inJ(3);
					int k = 0;
					int sizeI = 4;
					int sizeJ = 4;
					int sizeK = 4;
					int stepsOfIntegr_K = 2;
					double step_t_viaK = 1.0 / (sizeK);
					auto f = [](Point<double> X) ->double {return 1; };
					for (int K = 0; K < sizeK; K++)
					{
						std::vector<std::vector<std::vector<double>>> data_store_for_integrals_in_K(3);
						for (int k = 0; k < 3; k++)
						{
							data_store_for_integrals_in_K[k].resize(sizeJ);
							for (int j = 0; j < data_store_for_integrals_in_K[k].size(); j++)
							{
								data_store_for_integrals_in_K[k][j].resize(sizeI);
							}
						}

						std::vector < double > xj_integr1D, q_integr1D;
						geometry::Rectangle tmp_elem;
						tmp_elem.SetIntegrationLaw(stepsOfIntegr_K, xj_integr1D, q_integr1D);

						double current_t = 0;
						double start_t = K * step_t_viaK;
						double end_t = (K + 1)*step_t_viaK;
						for (int T = 0; T < stepsOfIntegr_K; T++)
						{
							double _currentZ = 0.0;
							double current_t = (start_t + end_t + xj_integr1D[T] * (end_t - start_t)) / 2.0;
							newCoordSystem System_for_t;
							NewSystem(current_t, System_for_t);

							Point<double> _Up(R_for_integration, R_for_integration, 0.0), _Down(-R_for_integration, -R_for_integration, 0.0);
							Point<double> _H((_Up.x - _Down.x) / sizeI, (_Up.y - _Down.y) / sizeJ, 0.0);
							for (int J = 0; J < sizeJ; J++)
							{
								for (int I = 0; I < sizeI; I++)
								{
									for (int k = 0; k < 3; k++)
									{
										auto Func_q = [&grid, &U, &_Up, &_Down](Point<double> &X, newCoordSystem &System_for_t, int k) ->double
										{
											double res = 0;
											Point<double> _X = System_for_t.transfer_from_OLD_into_NEW(X);

											int curr_el = grid.GetElementID(X);
											if (curr_el == -1) return 0.0;
											auto current_element = grid.GetElement(curr_el);


											/*if (current_element->temperature_loc.size() != 0)
											{
												std::vector<std::vector<double>> D_loc;
												grid.domain[current_element->get_id_domain()].forMech.get_D(current_element->get_temperature(X), 3, D_loc);
											}*/

											/*std::vector<std::vector<double>> SIG, EPS, _SIG, _EPS, dU_dX;
											dU_dX = current_element->get_value_dU_dX_mech(X, value_in_DOF);
											EPS = current_element->get_value_EPS_mech(X, value_in_DOF);
											SIG = current_element->get_value_SIGMA_mech(EPS, grid.domain[current_element->get_id_domain()]);*/

											Point<Point<double>> dU = grid.GetDerevativeFromSolutionInPoint(curr_el, X, U);
											Tensor2Rank3D _dU;
											_dU.val[0][0] = dU.x.x; _dU.val[0][1] = dU.x.y; _dU.val[0][2] = dU.x.z;
											_dU.val[1][0] = dU.y.x; _dU.val[1][1] = dU.y.y; _dU.val[1][2] = dU.y.z;
											_dU.val[2][0] = dU.z.x; _dU.val[2][1] = dU.z.y; _dU.val[2][2] = dU.z.z;

											Tensor2Rank3D EPS = grid.GetStressTensorFromSolutionInPoint(curr_el, dU);
											Tensor2Rank3D SIG = grid.GetStrainTensorFromSolutionInPoint(curr_el, X, EPS);

											Tensor2Rank3D _SIG = System_for_t.transfer_from_OLD_into_NEW(SIG);
											Tensor2Rank3D _EPS = System_for_t.transfer_from_OLD_into_NEW(SIG);

											auto W = [](double ksi, double ksi_max, double ksi_min, int atr) -> double
											{
												double res;
												switch (atr)
												{
												case 0:
													res = (ksi_max - ksi) / (ksi_max - ksi_min);
													break;
												case 1:
													res = (ksi - ksi_min) / (ksi_max - ksi_min);
													break;
												default:
													res = 0;
													break;
												}
												return res;
											};
											auto dW = [](double ksi, double ksi_max, double ksi_min, int atr) -> double
											{
												double res;
												switch (atr)
												{
												case 0:
													res = -1 / (ksi_max - ksi_min);
													break;
												case 1:
													res = 1 / (ksi_max - ksi_min);
													break;
												default:
													res = 0;
													break;
												}
												return res;
											};
											auto dq = [&](Point<double> &_X, Point<double> &_o, Point<double> &_up, Point<double> &_down) -> std::vector<double>
											{
												std::vector<double> dq_dxi(3);

												//quadr
												int base_func[2];
												Point<double> _up_loc, _down_loc;
												bool find = false;
												if (!find && _X.x <= _o.x && _X.y <= _o.y)//III
												{
													base_func[0] = 1;
													base_func[1] = 1;
													_down_loc = _down;
													_up_loc = _o;
													find = true;
												}
												if (!find && _X.x >= _o.x && _X.y <= _o.y)//IV
												{
													base_func[0] = 0;
													base_func[1] = 1;
													_down_loc = Point<double>(_o.x, _down.y, _o.z);
													_up_loc = Point<double>(_up.x, _o.y, _o.z);
													find = true;
												}
												if (!find && _X.x <= _o.x && _X.y >= _o.y)//II
												{
													base_func[0] = 1;
													base_func[1] = 0;
													_down_loc = Point<double>(_down.x, _o.y, _o.z);
													_up_loc = Point<double>(_o.x, _up.y, _o.z);
													find = true;
												}
												if (!find && _X.x >= _o.x && _X.y >= _o.y) //I
												{
													base_func[0] = 0;
													base_func[1] = 0;
													_down_loc = _o;
													_up_loc = _up;
													find = true;
												}

												dq_dxi[0] = dW(_X.x, _up_loc.x, _down_loc.x, base_func[0]) * W(_X.y, _up_loc.y, _down_loc.y, base_func[1]);
												dq_dxi[1] = W(_X.x, _up_loc.x, _down_loc.x, base_func[0]) * dW(_X.y, _up_loc.y, _down_loc.y, base_func[1]);
												dq_dxi[2] = 0;

												return dq_dxi;
											};

											//left part
											std::vector<double> left_res(3);
											for (int m = 0; m < 3; m++)
											{
												//for (int n = 0; n < 3; n++)
												//{
												//left_res[k] += (_SIG[m][n] * _EPS[m][n]) / 2.;
												//}
												left_res[k] += (_SIG.val[m][k] * _EPS.val[m][k]) / 2.;
											}

											//right_part
											std::vector<double> right_res(3);
											for (int i = 0; i < 3; i++)
											{
												for (int j = 0; j < 3; j++)
												{
													double duj_dxk = 0;
													for (int m = 0; m < 3; m++)
													{
														duj_dxk += _dU.val[j][m] * System_for_t._A.val[m][k];
													}
													right_res[i] += _SIG.val[i][j] * duj_dxk;
												}
											}

											//q
											std::vector<double> dq_dxi = dq(_X, System_for_t._O, _Up, _Down);

											for (int i = 0; i < 3; i++)
											{
												res += -1 * (left_res[i] - right_res[i])*dq_dxi[i];
											}

											//printf("elem[%d/%d/%d](%d): left_res=%.5lf right_res=%.5lf\n", I, J, K, k, left_res, right_res);
											return res;
										};

										Point<double> _up(_Down.x + _H.x*(I + 1), _Down.y + _H.y*(J + 1), 0.0);
										Point<double> _down(_Down.x + _H.x*(I), _Down.y + _H.y*(J), 0.0);
										double temp_res = 0;
										for (int i = 0; i < stepsOfIntegr_K; ++i)
										{
											double _currentX = (_down.x + _up.x + xj_integr1D[i] * (_up.x - _down.x)) / 2.0;
											for (int j = 0; j < stepsOfIntegr_K; ++j)
											{
												double _currentY = (_down.y + _up.y + xj_integr1D[i] * (_up.y - _down.y)) / 2.0;

												Point<double> _p(_currentX, _currentY, _currentZ);
												Point<double> X = System_for_t.transfer_from_NEW_into_OLD(_p);
												//temp_res +=	q_integr1D[i] * q_integr1D[j] * q_integr1D[T] *
												//	f(X)*det3(System_for_t._A);
												temp_res += q_integr1D[i] * q_integr1D[j] * q_integr1D[T] *
													Func_q(X, System_for_t, k)*math::GetDeterminantForMatrix3x3(System_for_t._A);
											}
										}

										double len_x = _up.x - _down.x;
										double len_y = _up.y - _down.y;
										double len_z = math::SolveLengthVector(F(end_t) - F(start_t));
										data_store_for_integrals_in_K[k][J][I] += temp_res * len_x * len_y * len_z / 8.0;
									}
								}
							}
						}

						for (int k = 0; k < 3; k++)
						{
							for (int j = 0; j < data_store_for_integrals_in_K[k].size(); j++)
							{
								for (int i = 0; i < data_store_for_integrals_in_K[k][j].size(); i++)
								{
									inJ[k] += data_store_for_integrals_in_K[k][j][i] /*/ SIGMA_PA_F*/;
								}
							}
						}
					}

					double Qrad = atan(inJ[1] / inJ[0]); //in rad
					//double Qrad = atan(inJ[0] / inJ[1]); //in rad
					double Qgrad = 180 * Qrad / M_PI; //in grad

					//t=0
					newCoordSystem s0;
					NewSystem(0, s0);
					Point<double> _R, R;
					_R.z = s0._O.z;
					_R.x = cos(Qrad) * R_for_integration;
					_R.y = sqrt(R_for_integration*R_for_integration - _R.x*_R.x) *math::GetSignum(Qrad);
					R = s0.transfer_from_NEW_into_OLD(_R);
					R -= s0.O;
					//printf("FR[%d], Point<double>[%d], J1[%.2e], J2[%.2e], J3[%.2e], Qgrad[%.2lf], R(%.2e, %.2e, %.2e)=%.2e\n", Fr, segm - 1, inJ[0], inJ[1], inJ[2], Qgrad, R.x, R.y, R.z, math::SolveLengthVector(R));
					//printf("FR[%d], Point<double>[%d], J1[%.2e], J2[%.2e], J3[%.2e], Qgrad[%.2lf], R=%.2e, _R=%.2e\n", Fr, segm - 1, inJ[0], inJ[1], inJ[2], Qgrad, math::SolveLengthVector(R), math::SolveLengthVector(_R));
					grid.cracks[Cr].propogation[Fr][segm - 1].vect_of_propogation.push_back(R);
					grid.cracks[Cr].propogation[Fr][segm - 1].id_old_point = grid.cracks[Cr].id_points_in_fronts[Fr][segm];
					grid.cracks[Cr].propogation[Fr][segm - 1].base_coord_system.push_back(newCoordSystem(s0));
					grid.cracks[Cr].propogation[Fr][segm - 1].Qgrad.push_back(Qgrad);
					
					//t=1
					newCoordSystem s1;
					NewSystem(1, s1);
					_R.z = s1._O.z;
					_R.x = cos(Qrad) * R_for_integration;
					_R.y = sqrt(R_for_integration*R_for_integration - _R.x*_R.x) *math::GetSignum(Qrad);
					R = s1.transfer_from_NEW_into_OLD(_R) - s1.O;
					//printf("FR[%d], Point<double>[%d], J1[%.2e], J2[%.2e], J3[%.2e], Qgrad[%.2lf], R(%.2e, %.2e, %.2e)=%.2e\n", Fr, segm, inJ[0], inJ[1], inJ[2], Qgrad, R.x, R.y, R.z, math::SolveLengthVector(R));
					grid.cracks[Cr].propogation[Fr][segm].vect_of_propogation.push_back(R);
					grid.cracks[Cr].propogation[Fr][segm].id_old_point = grid.cracks[Cr].id_points_in_fronts[Fr][segm + 1];
					grid.cracks[Cr].propogation[Fr][segm].base_coord_system.push_back(newCoordSystem(s1));
					grid.cracks[Cr].propogation[Fr][segm].Qgrad.push_back(Qgrad);
				}

				//create new cracks
				//new nodes
				int old_size_xyz = grid.cracks[Cr].xyz.size();
				if (false) {
					for (int p = 0; p < grid.cracks[Cr].propogation[Fr].size(); p++)
					{
						Point<double> middle_vect;
						for (int v = 0; v < grid.cracks[Cr].propogation[Fr][p].vect_of_propogation.size(); v++)
						{
							middle_vect += grid.cracks[Cr].propogation[Fr][p].vect_of_propogation[v] / grid.cracks[Cr].propogation[Fr][p].vect_of_propogation.size();
						}
						middle_vect = middle_vect / math::SolveLengthVector(middle_vect) * grid.cracks[Cr].R;
						//printf("Len(FR[%d], Point[%d]) = %.2e\n", Fr, p, math::SolveLengthVector(middle_vect));
						grid.cracks[Cr].xyz.push_back(grid.cracks[Cr].xyz[grid.cracks[Cr].propogation[Fr][p].id_old_point] + middle_vect);
						grid.cracks[Cr].propogation[Fr][p].id_new_point = grid.cracks[Cr].xyz.size() - 1;
					}
				}

				//FULL MIDDLE MOVE!!!
				if (false)
				{
					Point<double> FULL_MIDDLE_VECTOR;
					for (int p = 0; p < false && grid.cracks[Cr].propogation[Fr].size(); p++)
					{
						Point<double> middle_vect;
						for (int v = 0; v < grid.cracks[Cr].propogation[Fr][p].vect_of_propogation.size(); v++)
						{
							middle_vect += grid.cracks[Cr].propogation[Fr][p].vect_of_propogation[v] / grid.cracks[Cr].propogation[Fr][p].vect_of_propogation.size();
						}
						middle_vect = middle_vect / math::SolveLengthVector(middle_vect) * grid.cracks[Cr].R;
						FULL_MIDDLE_VECTOR += middle_vect;
					}
					FULL_MIDDLE_VECTOR /= grid.cracks[Cr].propogation[Fr].size();
					for (int p = 0; p < grid.cracks[Cr].propogation[Fr].size(); p++)
					{
						grid.cracks[Cr].xyz.push_back(grid.cracks[Cr].xyz[grid.cracks[Cr].propogation[Fr][p].id_old_point] + FULL_MIDDLE_VECTOR);
						grid.cracks[Cr].propogation[Fr][p].id_new_point = grid.cracks[Cr].xyz.size() - 1;
					}
				}

				//smoothing via Ns nodes
				int Ns = grid.cracks[Cr].smoothing_coefficient;
				for (int id_base_vertex = 0; id_base_vertex < grid.cracks[Cr].propogation[Fr].size(); id_base_vertex++)
				{
					int start_id = id_base_vertex;
					if (id_base_vertex + Ns >= grid.cracks[Cr].propogation[Fr].size())
					{
						start_id = grid.cracks[Cr].propogation[Fr].size() - Ns;
						if (start_id < 0) start_id = 0;
					}

					double middle_Qgrad = 0;
					int n_i = 0;

					for (int i = 0; i < Ns && (start_id + i) < grid.cracks[Cr].propogation[Fr].size(); i++)
					{
						for (int q = 0; q < grid.cracks[Cr].propogation[Fr][start_id+i].Qgrad.size(); q++)
						{
							middle_Qgrad += grid.cracks[Cr].propogation[Fr][start_id+i].Qgrad[q];
							n_i++;
						}
					}
					middle_Qgrad /= n_i;

					Point<double> propogation_vector;
					Point<double> simple_vector;
					double simpleQ = 0;
					for (int q = 0; q < grid.cracks[Cr].propogation[Fr][id_base_vertex].base_coord_system.size(); q++)
					{
						propogation_vector += grid.cracks[Cr].propogation[Fr][id_base_vertex].base_coord_system[q].MakePropogation(middle_Qgrad, grid.cracks[Cr].R)/
							grid.cracks[Cr].propogation[Fr][id_base_vertex].base_coord_system.size();
						simple_vector += grid.cracks[Cr].propogation[Fr][id_base_vertex].vect_of_propogation[q] /
							grid.cracks[Cr].propogation[Fr][id_base_vertex].base_coord_system.size();
						simpleQ += grid.cracks[Cr].propogation[Fr][id_base_vertex].Qgrad[q]/
							grid.cracks[Cr].propogation[Fr][id_base_vertex].base_coord_system.size();
					}
					propogation_vector = propogation_vector / math::SolveLengthVector(propogation_vector) * grid.cracks[Cr].R;
					simple_vector = simple_vector / math::SolveLengthVector(simple_vector) * grid.cracks[Cr].R;

					/*printf_s("Front[%d], old_id[%d], Q=%.2lf; OldPoint(%.2lf, %.2lf, %.2lf); Vector=(%.2lf, %.2lf, %.2lf); SimpleQ=%.2lf; SimpleVector=(%.2lf, %.2lf, %.2lf)\n", 
						Fr, 
						grid.cracks[Cr].propogation[Fr][id_base_vertex].id_old_point, 
						middle_Qgrad, 
						grid.cracks[Cr].xyz[grid.cracks[Cr].propogation[Fr][id_base_vertex].id_old_point].x,
						grid.cracks[Cr].xyz[grid.cracks[Cr].propogation[Fr][id_base_vertex].id_old_point].y,
						grid.cracks[Cr].xyz[grid.cracks[Cr].propogation[Fr][id_base_vertex].id_old_point].z,
						propogation_vector.x,
						propogation_vector.y,
						propogation_vector.z,
						simpleQ,
						simple_vector.x,
						simple_vector.y,
						simple_vector.z);*/

					grid.cracks[Cr].xyz.push_back(grid.cracks[Cr].xyz[grid.cracks[Cr].propogation[Fr][id_base_vertex].id_old_point] + propogation_vector);
					grid.cracks[Cr].propogation[Fr][id_base_vertex].id_new_point = grid.cracks[Cr].xyz.size() - 1;
				}


				//new triangles and fronts
				Point<double> true_normal = grid.cracks[Cr].triangles[0].GetNormal();
				for (int segm = 0; segm < grid.cracks[Cr].fronts[Fr].size(); segm++)
				{
					Point<double> new_normal;
					int cur_segm[3] = { grid.cracks[Cr].fronts[Fr][segm].id_base_triangles, grid.cracks[Cr].fronts[Fr][segm].id_left, grid.cracks[Cr].fronts[Fr][segm].id_right };
					int new_segm[3];
					for (int ii = 1; ii < 3; ii++)
					{
						for (int p = 0; p < grid.cracks[Cr].propogation[Fr].size(); p++)
						{
							if (cur_segm[ii] == grid.cracks[Cr].propogation[Fr][p].id_old_point)
							{
								new_segm[ii] = grid.cracks[Cr].propogation[Fr][p].id_new_point;
								break;
							}
						}
					}

					for (int ii = 0; ii < grid.cracks[Cr].triangles.size(); ii++)
					{
						for (int jj = 0; jj < grid.cracks[Cr].triangles[ii].GetNodesCount(); jj++)
						{
							grid.cracks[Cr].triangles[ii].SetNode(jj, &grid.cracks[Cr].xyz[grid.cracks[Cr].triangles[ii].GetIdNode(jj)]);
						}
					}
					
					grid.cracks[Cr].triangles.push_back(geometry::Triangle());
					grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].SetIdNode(0, cur_segm[1]);
					grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].SetIdNode(1, cur_segm[2]);
					grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].SetIdNode(2, new_segm[2]);
					grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].SetNode(0, &grid.cracks[Cr].xyz[grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].GetIdNode(0)]);
					grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].SetNode(1, &grid.cracks[Cr].xyz[grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].GetIdNode(1)]);
					grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].SetNode(2, &grid.cracks[Cr].xyz[grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].GetIdNode(2)]);
					new_normal = grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].GetNormal();
					if (true_normal*new_normal < 0)
					{
						grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].SetIdNode(0, cur_segm[2]);
						grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].SetIdNode(1, cur_segm[1]);
						grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].SetIdNode(2, new_segm[2]);
						grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].SetNode(0, &grid.cracks[Cr].xyz[grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].GetIdNode(0)]);
						grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].SetNode(1, &grid.cracks[Cr].xyz[grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].GetIdNode(1)]);
						grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].SetNode(2, &grid.cracks[Cr].xyz[grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].GetIdNode(2)]);
					}

					grid.cracks[Cr].triangles.push_back(geometry::Triangle());
					grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].SetIdNode(0, cur_segm[1]);
					grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].SetIdNode(1, new_segm[1]);
					grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].SetIdNode(2, new_segm[2]);
					grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].SetNode(0, &grid.cracks[Cr].xyz[grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].GetIdNode(0)]);
					grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].SetNode(1, &grid.cracks[Cr].xyz[grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].GetIdNode(1)]);
					grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].SetNode(2, &grid.cracks[Cr].xyz[grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].GetIdNode(2)]);
					new_normal = grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].GetNormal();
					if (true_normal*new_normal < 0)
					{
						grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].SetIdNode(0, cur_segm[1]);
						grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].SetIdNode(1, new_segm[2]);
						grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].SetIdNode(2, new_segm[1]);
						grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].SetNode(0, &grid.cracks[Cr].xyz[grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].GetIdNode(0)]);
						grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].SetNode(1, &grid.cracks[Cr].xyz[grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].GetIdNode(1)]);
						grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].SetNode(2, &grid.cracks[Cr].xyz[grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].GetIdNode(2)]);
					}

					new_segm[0] = grid.cracks[Cr].triangles.size() - 1;
					grid.cracks[Cr].fronts[Fr][segm].id_base_triangles = new_segm[0];
					grid.cracks[Cr].fronts[Fr][segm].id_left = new_segm[1];
					grid.cracks[Cr].fronts[Fr][segm].id_right = new_segm[2];
				}
			}
		}

		FILE* fout;
		char new_crack[1000];
		for (int Cr = 0; Cr < grid.cracks.size(); Cr++)
		{
			//sprintf_s(new_crack, "%s/Crack_%d_step%d.txt", dir, Cr, curr_STEP + 1);
			sprintf_s(new_crack, "%s/Crack_%d.txt", dir, Cr, curr_STEP + 1);
			fopen_s(&fout, new_crack, "w");
			fprintf_s(fout, "%d\n", (int)(grid.cracks[Cr].xyz.size()));
			for (int i = 0; i < grid.cracks[Cr].xyz.size(); i++)
			{
				fprintf_s(fout, "%.6e\t%.6e\t%.6e\n",
					grid.cracks[Cr].xyz[i].x,
					grid.cracks[Cr].xyz[i].y,
					grid.cracks[Cr].xyz[i].z);
			}

			fprintf_s(fout, "\n%d\n", (int)grid.cracks[Cr].triangles.size());
			for (int i = 0; i < grid.cracks[Cr].triangles.size(); i++)
			{
				fprintf_s(fout, "%d\t%d\t%d\n",
					grid.cracks[Cr].triangles[i].GetIdNode(0),
					grid.cracks[Cr].triangles[i].GetIdNode(1),
					grid.cracks[Cr].triangles[i].GetIdNode(2));
			}

			fprintf_s(fout, "\n%d\n", (int)grid.cracks[Cr].fronts.size());
			for (int i = 0; i < grid.cracks[Cr].fronts.size(); i++)
			{
				fprintf_s(fout, "%d\n", (int)grid.cracks[Cr].fronts[i].size());
				for (int j = 0; j < grid.cracks[Cr].fronts[i].size(); j++)
				{
					fprintf_s(fout, "%d\t%d\t%d\n",
						grid.cracks[Cr].fronts[i][j].id_base_triangles,
						grid.cracks[Cr].fronts[i][j].id_left,
						grid.cracks[Cr].fronts[i][j].id_right);
				}
			}
			fclose(fout);

			{ //tech_plot - as middle plane
				sprintf_s(new_crack, "%s/Crack_%d_step%d.dat", dir, Cr, curr_STEP + 1);
				fopen_s(&fout, new_crack, "w");

				fprintf_s(fout, "TITLE     = \"numerical\"\n");
				fprintf_s(fout, "VARIABLES = \"x\"\n \"y\"\n \"z\"\n");
				fprintf_s(fout, "\"crack\"\n");
				fprintf_s(fout, "ZONE T=\"Crack_%d_step%d\"\n", Cr, curr_STEP + 1);
				fprintf_s(fout, " N=%d,  E=%d, F=FEBLOCK ET=Triangle \n", (int)grid.cracks[Cr].xyz.size(), (int)grid.cracks[Cr].triangles.size());
				fprintf_s(fout, " VARLOCATION=(NODAL NODAL NODAL CELLCENTERED)\n");

				for (int i = 0; i < grid.cracks[Cr].xyz.size(); i++)
					fprintf_s(fout, "%.10e\n", grid.cracks[Cr].xyz[i].x);
				fprintf_s(fout, "\n");
				for (int i = 0; i < grid.cracks[Cr].xyz.size(); i++)
					fprintf_s(fout, "%.10e\n", grid.cracks[Cr].xyz[i].y);
				fprintf_s(fout, "\n");
				for (int i = 0; i < grid.cracks[Cr].xyz.size(); i++)
					fprintf_s(fout, "%.10e\n", grid.cracks[Cr].xyz[i].z);
				fprintf_s(fout, "\n");

				for (int i = 0; i < grid.cracks[Cr].triangles.size(); i++)
					fprintf_s(fout, "%d\n", curr_STEP + 1);
				fprintf_s(fout, "\n");

				for (int i = 0; i < grid.cracks[Cr].triangles.size(); i++)
				{
					for (int j = 0; j < 3; j++)
						fprintf_s(fout, "%d ", grid.cracks[Cr].triangles[i].GetIdNode(j) + 1);
					fprintf_s(fout, "\n");
				}

				fclose(fout);
			}

			{ //tech_plot - as face planes
				std::vector<Point<double>> xyz_up(grid.cracks[Cr].xyz.size()), xyz_down(grid.cracks[Cr].xyz.size());
				std::vector<int> elem_into_nodes(grid.cracks[Cr].xyz.size());
				double MaxSize = 0, MinSize = 1e+15;

				for (int id = 0; id < grid.cracks[Cr].triangles.size(); id++)
				{
					Point<double> normal = grid.cracks[Cr].triangles[id].GetNormal();
					normal /= math::SolveLengthVector(normal);
					Point<double> O = grid.cracks[Cr].triangles[id].GetWeightCentr();
					double move_koef = 0.5;
					
					Point<double> up_X = O + normal * move_koef;
					int up_elem = grid.GetElementID(up_X);
					Point<double> move_vector_UP = grid.GetSolutionInPoint(up_elem, up_X, U);
					
					Point<double> down_X = O - normal * move_koef;
					int down_elem = grid.GetElementID(down_X);
					Point<double> move_vector_DOWN = grid.GetSolutionInPoint(down_elem, down_X, U);

					double move_size = math::SolveLengthVector(move_vector_DOWN) + math::SolveLengthVector(move_vector_UP);
					if (move_size > MaxSize) MaxSize = move_size;
					if (move_size < MinSize) MinSize = move_size;

					for (int _i = 0; _i < 3; _i++)
					{
						elem_into_nodes[grid.cracks[Cr].triangles[id].GetIdNode(_i)] ++;
						xyz_up[grid.cracks[Cr].triangles[id].GetIdNode(_i)] += move_vector_UP;
						xyz_down[grid.cracks[Cr].triangles[id].GetIdNode(_i)] += move_vector_DOWN;
					}
				}
				for (int i = 0; i < grid.cracks[Cr].fronts.size(); i++)
				{
					for (int j = 0; j < grid.cracks[Cr].fronts[i].size(); j++)
					{
						xyz_up[grid.cracks[Cr].fronts[i][j].id_left] = Point<double>(0, 0, 0);
						xyz_up[grid.cracks[Cr].fronts[i][j].id_right] = Point<double>(0, 0, 0);

						xyz_down[grid.cracks[Cr].fronts[i][j].id_left] = Point<double>(0, 0, 0);
						xyz_down[grid.cracks[Cr].fronts[i][j].id_right] = Point<double>(0, 0, 0);
					}
				}


				for (int id = 0; id < xyz_up.size(); id++)
				{
					xyz_up[id] /= elem_into_nodes[id];
					xyz_down[id] /= elem_into_nodes[id];
					xyz_up[id] += grid.cracks[Cr].xyz[id];
					xyz_down[id] += grid.cracks[Cr].xyz[id];
				}



				sprintf_s(new_crack, "%s/Crack_%d_step%d_in3D_MaxSize_%.2e_MinSize_%.2e.dat", dir, Cr, curr_STEP + 1, MaxSize, MinSize);
				fopen_s(&fout, new_crack, "w");

				fprintf_s(fout, "TITLE     = \"numerical\"\n");
				fprintf_s(fout, "VARIABLES = \"x\"\n \"y\"\n \"z\"\n");
				fprintf_s(fout, "\"crack\"\n");
				fprintf_s(fout, "ZONE T=\"Crack_%d_step%d\"\n", Cr, curr_STEP + 1);
				fprintf_s(fout, " N=%d,  E=%d, F=FEBLOCK ET=Triangle \n", (int)xyz_up.size()*2, (int)grid.cracks[Cr].triangles.size()*2);
				fprintf_s(fout, " VARLOCATION=(NODAL NODAL NODAL CELLCENTERED)\n");

				for (int i = 0; i < xyz_up.size(); i++)
					fprintf_s(fout, "%.10e\n", xyz_up[i].x);
				for (int i = 0; i < xyz_down.size(); i++)
					fprintf_s(fout, "%.10e\n", xyz_down[i].x);
				fprintf_s(fout, "\n");
				for (int i = 0; i < xyz_up.size(); i++)
					fprintf_s(fout, "%.10e\n", xyz_up[i].y);
				for (int i = 0; i < xyz_down.size(); i++)
					fprintf_s(fout, "%.10e\n", xyz_down[i].y);
				fprintf_s(fout, "\n");
				for (int i = 0; i < xyz_up.size(); i++)
					fprintf_s(fout, "%.10e\n", xyz_up[i].z);
				for (int i = 0; i < xyz_down.size(); i++)
					fprintf_s(fout, "%.10e\n", xyz_down[i].z);
				fprintf_s(fout, "\n");
								

				for (int i = 0; i < grid.cracks[Cr].triangles.size()*2; i++)
					fprintf_s(fout, "%d\n", curr_STEP + 1);
				fprintf_s(fout, "\n");

				for (int i = 0; i < grid.cracks[Cr].triangles.size(); i++)
				{
					for (int j = 0; j < 3; j++)
						fprintf_s(fout, "%d ", grid.cracks[Cr].triangles[i].GetIdNode(j) + 1);
					fprintf_s(fout, "\n");
				}
				for (int i = 0; i < grid.cracks[Cr].triangles.size(); i++)
				{
					for (int j = 0; j < 3; j++)
						fprintf_s(fout, "%d ", (int)grid.cracks[Cr].xyz.size() + grid.cracks[Cr].triangles[i].GetIdNode(j) + 1);
					fprintf_s(fout, "\n");
				}

				fclose(fout);
			}
		}
	}
	void CrackPropagation_3D_Cherepanov(MultiXFEM::Grid &grid, std::vector<Point<double>> &U, char *dir, int curr_STEP)
	{
		//to find the propagation vector for each a front point
		for (int Cr = 0; Cr < grid.cracks.size(); Cr++)
		{
			double R_for_integration = grid.cracks[Cr].R;
			for (int Fr = 0; Fr < grid.cracks[Cr].id_points_in_fronts.size(); Fr++)
			{
				grid.cracks[Cr].propogation.resize(grid.cracks[Cr].id_points_in_fronts.size());
				grid.cracks[Cr].propogation[Fr].resize(grid.cracks[Cr].id_points_in_fronts[Fr].size() - 2);
				for (int segm = 1; segm < grid.cracks[Cr].id_points_in_fronts[Fr].size() - 2; segm++)
				{
					Point<double> P0, P1, P2, P3;
					P0 = grid.cracks[Cr].xyz[grid.cracks[Cr].id_points_in_fronts[Fr][segm - 1]];
					P1 = grid.cracks[Cr].xyz[grid.cracks[Cr].id_points_in_fronts[Fr][segm]];
					P2 = grid.cracks[Cr].xyz[grid.cracks[Cr].id_points_in_fronts[Fr][segm + 1]];
					P3 = grid.cracks[Cr].xyz[grid.cracks[Cr].id_points_in_fronts[Fr][segm + 2]];

					if (segm - 1 == 0 && grid.cracks[Cr].id_points_in_fronts[Fr][segm - 1] == grid.cracks[Cr].id_points_in_fronts[Fr][segm])
					{
						P0 = P1 + (P1 - P2);
					}
					if ((segm + 2) == (grid.cracks[Cr].id_points_in_fronts[Fr].size() - 1) &&
						grid.cracks[Cr].id_points_in_fronts[Fr][segm + 1] == grid.cracks[Cr].id_points_in_fronts[Fr][segm + 2])
					{
						P3 = P2 + (P2 - P1);
					}

					Point<double> r3 = P1 * 3 - P2 * 3 + P3 + P0;
					Point<double> r2 = P1 * (-5) + P2 * 4 - P3;
					Point<double> r1 = P0 * (-1) + P2;
					Point<double> r0 = P1 * 2;
					r3 = P0 * (-1) + P1 * 3 - P2 * 3 + P3;
					r2 = P0 * 2 + P1 * (-5) + P2 * 4 - P3;
					r1 = P0 * (-1) + P2;
					r0 = P1 * 2;
					auto F = [&](double t) -> Point<double>
					{
						return (r3*pow(t, 3) + r2 * pow(t, 2) + r1 * t + r0) / 2.;
					};

					Point<double> _r2 = P1 * 9 - P2 * 9 + P3 * 3 + P0 * 3;
					Point<double> _r1 = P1 * (-10) + P2 * 8 - P3 * 2;
					Point<double> _r0 = P0 * (-1) + P2;
					_r2 = r3 * 3.;
					_r1 = r2 * 2.;
					_r0 = r1;
					auto dF_dt = [&](double t) -> Point<double>
					{
						return (_r2*pow(t, 2) + _r1 * t + _r0) / 2.;
					};

					struct newCoordSystem {
						std::vector<Point<double>> _X; //new basis in OLD coord
						Point<double> O; //New coord centr in OLD coord
						Point<double> _O; //New coord centr in NEW coord
						Tensor2Rank3D A; //transfer from Old into New
						Tensor2Rank3D _A; //transfer from New into Old
						Point<double> Offset; //offset from Old into New
						Point<double> _Offset; //offset from New into Old 

						Point<double> transfer_from_OLD_into_NEW(Point<double> &T)
						{
							return math::MakeTransferPoint(T, A, Offset);
						}
						Point<double> transfer_from_NEW_into_OLD(Point<double> &_T)
						{
							return math::MakeTransferPoint(_T, _A, _Offset);
						}
						Tensor2Rank3D transfer_from_OLD_into_NEW(Tensor2Rank3D &T)
						{
							return math::MakeTransferTensor(T, A, Offset);
						}
					};
					auto NewSystem = [Cr, Fr, &grid, &dF_dt, &F](double t, newCoordSystem &System) -> void
					{
						System._X.resize(3);
						int curr_segm;
						System._X[2] = dF_dt(t);
						curr_segm = 0;

						Point<double> Line[2];
						Line[0] = F(t);
						Line[1] = F(t) + System._X[2];
						Point<double> test;
						for (int tt = 0; tt < 3; tt++)
						{
							if (grid.cracks[Cr].fronts[Fr][curr_segm].id_left !=
								grid.cracks[Cr].triangles[grid.cracks[Cr].fronts[Fr][curr_segm].id_base_triangles].GetIdNode(tt)
								&& grid.cracks[Cr].fronts[Fr][curr_segm].id_right !=
								grid.cracks[Cr].triangles[grid.cracks[Cr].fronts[Fr][curr_segm].id_base_triangles].GetIdNode(tt))
							{
								test = grid.cracks[Cr].xyz[grid.cracks[Cr].triangles[grid.cracks[Cr].fronts[Fr][curr_segm].id_base_triangles].GetIdNode(tt)];
								break;
							}
						}
						DenseMatrix<double, double> test_matrix;
						System._X[0] = math::MakeNormalForLine(Line, test, test_matrix);
						System._X[1] = grid.cracks[Cr].triangles[grid.cracks[Cr].fronts[Fr][curr_segm].id_base_triangles].GetNormal();

						System._X[0] /= math::SolveLengthVector(System._X[0]);
						System._X[1] /= math::SolveLengthVector(System._X[1]);
						System._X[2] /= math::SolveLengthVector(System._X[2]);

						System.A.val[0][0] = System._X[0].x; System.A.val[0][1] = System._X[0].y; System.A.val[0][2] = System._X[0].z;
						System.A.val[1][0] = System._X[1].x; System.A.val[1][1] = System._X[1].y; System.A.val[1][2] = System._X[1].z;
						System.A.val[2][0] = System._X[2].x; System.A.val[2][1] = System._X[2].y; System.A.val[2][2] = System._X[2].z;
						System._A = math::SolveInverseMatrix3x3(System.A);

						Point<double> _p;
						System.Offset = F(t)*(-1.0);
						System._Offset = math::MakeTransferPoint(_p, System.A, System.Offset)*(-1);

						Point<double> _f;
						_f = F(t);
						System.O = F(t);
						System._O = System.transfer_from_OLD_into_NEW(_f);
						System.O = System.transfer_from_NEW_into_OLD(System._O);
					};

					//integration in IxJxK via NEW СOORDINATES
					std::vector<double> inJ(3);
					int k = 0;
					int sizeI = 2;
					int sizeJ = 2;
					int sizeK = 2;
					int stepsOfIntegr_K = 2;
					double step_t_viaK = 1.0 / (sizeK);
					auto f = [](Point<double> X) ->double {return 1; };
					for (int K = 0; K < sizeK; K++)
					{
						std::vector<std::vector<std::vector<double>>> data_store_for_integrals_in_K(3);
						for (int k = 0; k < 3; k++)
						{
							data_store_for_integrals_in_K[k].resize(sizeJ);
							for (int j = 0; j < data_store_for_integrals_in_K[k].size(); j++)
							{
								data_store_for_integrals_in_K[k][j].resize(sizeI);
							}
						}

						std::vector < double > xj_integr1D, q_integr1D;
						geometry::Rectangle tmp_elem;
						tmp_elem.SetIntegrationLaw(stepsOfIntegr_K, xj_integr1D, q_integr1D);

						double current_t = 0;
						double start_t = K * step_t_viaK;
						double end_t = (K + 1)*step_t_viaK;
						for (int T = 0; T < stepsOfIntegr_K; T++)
						{
							double _currentZ = 0.0;
							double current_t = (start_t + end_t + xj_integr1D[T] * (end_t - start_t)) / 2.0;
							newCoordSystem System_for_t;
							NewSystem(current_t, System_for_t);

							Point<double> _Up(R_for_integration, R_for_integration, 0.0), _Down(-R_for_integration, -R_for_integration, 0.0);
							Point<double> _H((_Up.x - _Down.x) / sizeI, (_Up.y - _Down.y) / sizeJ, 0.0);
							for (int J = 0; J < sizeJ; J++)
							{
								for (int I = 0; I < sizeI; I++)
								{
									for (int k = 0; k < 3; k++)
									{
										auto Func_q = [&grid, &U, &_Up, &_Down](Point<double> &X, newCoordSystem &System_for_t, int k) ->double
										{
											double res = 0;
											Point<double> _X = System_for_t.transfer_from_OLD_into_NEW(X);

											int curr_el = grid.GetElementID(X);
											if (curr_el == -1) return 0.0;
											auto current_element = grid.GetElement(curr_el);


											/*if (current_element->temperature_loc.size() != 0)
											{
												std::vector<std::vector<double>> D_loc;
												grid.domain[current_element->get_id_domain()].forMech.get_D(current_element->get_temperature(X), 3, D_loc);
											}*/

											/*std::vector<std::vector<double>> SIG, EPS, _SIG, _EPS, dU_dX;
											dU_dX = current_element->get_value_dU_dX_mech(X, value_in_DOF);
											EPS = current_element->get_value_EPS_mech(X, value_in_DOF);
											SIG = current_element->get_value_SIGMA_mech(EPS, grid.domain[current_element->get_id_domain()]);*/

											Point<Point<double>> dU = grid.GetDerevativeFromSolutionInPoint(curr_el, X, U);
											Tensor2Rank3D _dU_old, _dU_new;
											Point<Point<double>> dU_new;
											_dU_old.val[0][0] = dU.x.x; _dU_old.val[0][1] = dU.x.y; _dU_old.val[0][2] = dU.x.z;
											_dU_old.val[1][0] = dU.y.x; _dU_old.val[1][1] = dU.y.y; _dU_old.val[1][2] = dU.y.z;
											_dU_old.val[2][0] = dU.z.x; _dU_old.val[2][1] = dU.z.y; _dU_old.val[2][2] = dU.z.z;
											for (int i = 0; i < 3; i++)
											{
												for (int j = 0; j < 3; j++)
												{
													for (int m = 0; m < 3; m++)
													{
														_dU_new.val[i][j] += _dU_old.val[i][m] * System_for_t._A.val[m][j];
													}
												}
											}
											dU_new.x.x = _dU_new.val[0][0]; dU_new.x.y = _dU_new.val[0][1]; dU_new.x.z = _dU_new.val[0][2];
											dU_new.y.x = _dU_new.val[1][0]; dU_new.y.y = _dU_new.val[1][1]; dU_new.y.z = _dU_new.val[1][2];
											dU_new.z.x = _dU_new.val[2][0]; dU_new.z.y = _dU_new.val[2][1]; dU_new.z.z = _dU_new.val[2][2];

											Tensor2Rank3D _EPS = grid.GetStressTensorFromSolutionInPoint(curr_el, dU_new);
											Tensor2Rank3D _SIG = grid.GetStrainTensorFromSolutionInPoint(curr_el, X, _EPS);

											/*Tensor2Rank3D EPS = grid.GetStressTensorFromSolutionInPoint(curr_el, dU);
											Tensor2Rank3D SIG = grid.GetStrainTensorFromSolutionInPoint(curr_el, X, EPS);
											Tensor2Rank3D __SIG = System_for_t.transfer_from_OLD_into_NEW(SIG);
											Tensor2Rank3D __EPS = System_for_t.transfer_from_OLD_into_NEW(EPS);*/

											auto W = [](double ksi, double ksi_max, double ksi_min, int atr) -> double
											{
												double res;
												switch (atr)
												{
												case 0:
													res = (ksi_max - ksi) / (ksi_max - ksi_min);
													break;
												case 1:
													res = (ksi - ksi_min) / (ksi_max - ksi_min);
													break;
												default:
													res = 0;
													break;
												}
												return res;
											};
											auto dW = [](double ksi, double ksi_max, double ksi_min, int atr) -> double
											{
												double res;
												switch (atr)
												{
												case 0:
													res = -1 / (ksi_max - ksi_min);
													break;
												case 1:
													res = 1 / (ksi_max - ksi_min);
													break;
												default:
													res = 0;
													break;
												}
												return res;
											};
											auto dq = [&](Point<double> &_X, Point<double> &_o, Point<double> &_up, Point<double> &_down) -> std::vector<double>
											{
												std::vector<double> dq_dxi(3);

												//quadr
												int base_func[2];
												Point<double> _up_loc, _down_loc;
												bool find = false;
												if (!find && _X.x <= _o.x && _X.y <= _o.y)//III
												{
													base_func[0] = 1;
													base_func[1] = 1;
													_down_loc = _down;
													_up_loc = _o;
													find = true;
												}
												if (!find && _X.x >= _o.x && _X.y <= _o.y)//IV
												{
													base_func[0] = 0;
													base_func[1] = 1;
													_down_loc = Point<double>(_o.x, _down.y, _o.z);
													_up_loc = Point<double>(_up.x, _o.y, _o.z);
													find = true;
												}
												if (!find && _X.x <= _o.x && _X.y >= _o.y)//II
												{
													base_func[0] = 1;
													base_func[1] = 0;
													_down_loc = Point<double>(_down.x, _o.y, _o.z);
													_up_loc = Point<double>(_o.x, _up.y, _o.z);
													find = true;
												}
												if (!find && _X.x >= _o.x && _X.y >= _o.y) //I
												{
													base_func[0] = 0;
													base_func[1] = 0;
													_down_loc = _o;
													_up_loc = _up;
													find = true;
												}

												dq_dxi[0] = dW(_X.x, _up_loc.x, _down_loc.x, base_func[0]) * W(_X.y, _up_loc.y, _down_loc.y, base_func[1]);
												dq_dxi[1] = W(_X.x, _up_loc.x, _down_loc.x, base_func[0]) * dW(_X.y, _up_loc.y, _down_loc.y, base_func[1]);
												dq_dxi[2] = 0; // W(_X.x, _up_loc.x, _down_loc.x, base_func[0])*W(_X.y, _up_loc.y, _down_loc.y, base_func[1]);

												return dq_dxi;
											};

											//left part
											std::vector<double> _left_res(3);
											double left_res = 0;
											for (int m = 0; m < 3; m++)
											{
												_left_res[m] = (_SIG.val[m][k] * _EPS.val[m][k]) / 2.;
												left_res += (_SIG.val[m][k] * _EPS.val[m][k]) / 2.;
											}

											//right_part
											std::vector<double> _right_res(3);
											double right_res =0;
											for (int j = 0; j < 3; j++)
											{
												for (int i = 0; i < 3; i++)
												{
													double dui_dxk = 0;
													for (int m = 0; m < 3; m++)
													{
														dui_dxk += _dU_old.val[j][m] * System_for_t._A.val[m][k];
													}
													_right_res[j] += _SIG.val[i][j] * dui_dxk;
												}
												right_res += _right_res[j];
											}

											//q
											std::vector<double> dq_dxj = dq(_X, System_for_t._O, _Up, _Down);
											//if(k == 0) printf("%.5lf %.5lf %.5lf %.5lf\n", _X.x, _X.y, _X.z, dq_dxi[2]);

											for (int j = 0; j < 3; j++)
											{
												res += -1 * (_left_res[j] - _right_res[j])*dq_dxj[j]; //??
												//res += -1 * (left_res - _right_res[j])*dq_dxj[j]; //??
											}
											double res_new = 0;
											res_new = -1 * (left_res - right_res)*dq_dxj[k];//??

											//printf("elem[%d/%d/%d](%d): left_res=%.5lf right_res=%.5lf\n", I, J, K, k, left_res, right_res);
											return res_new;
										};

										Point<double> _up(_Down.x + _H.x*(I + 1), _Down.y + _H.y*(J + 1), 0.0);
										Point<double> _down(_Down.x + _H.x*(I), _Down.y + _H.y*(J), 0.0);
										double temp_res = 0;
										for (int i = 0; i < stepsOfIntegr_K; ++i)
										{
											double _currentX = (_down.x + _up.x + xj_integr1D[i] * (_up.x - _down.x)) / 2.0;
											for (int j = 0; j < stepsOfIntegr_K; ++j)
											{
												double _currentY = (_down.y + _up.y + xj_integr1D[j] * (_up.y - _down.y)) / 2.0;

												Point<double> _p(_currentX, _currentY, _currentZ);
												Point<double> X = System_for_t.transfer_from_NEW_into_OLD(_p);
												//temp_res +=	q_integr1D[i] * q_integr1D[j] * q_integr1D[T] *
												//	f(X)*det3(System_for_t._A);
												temp_res += q_integr1D[i] * q_integr1D[j] * q_integr1D[T] *
													Func_q(X, System_for_t, k)*math::GetDeterminantForMatrix3x3(System_for_t._A);
											}
										}

										double len_x = _up.x - _down.x;
										double len_y = _up.y - _down.y;
										double len_z = math::SolveLengthVector(F(end_t) - F(start_t));
										data_store_for_integrals_in_K[k][J][I] += temp_res * len_x * len_y * len_z / 8.0;
									}
								}
							}
						}

						for (int k = 0; k < 3; k++)
						{
							for (int j = 0; j < data_store_for_integrals_in_K[k].size(); j++)
							{
								for (int i = 0; i < data_store_for_integrals_in_K[k][j].size(); i++)
								{
									inJ[k] += data_store_for_integrals_in_K[k][j][i] /*/ SIGMA_PA_F*/;
								}
							}
						}
					}

					//solve the SIF: K_I, K_II
					double K_I = 0, K_II = 0;
					{
						//solving a point with progection into both side of crack faces
						Point<double> _P_up, P_up; // _P - in LOCAL coordinates, P - in GLOBAL
						Point<double> _P_down, P_down; // _P - in LOCAL coordinates, P - in GLOBAL
						newCoordSystem s_middle;
						NewSystem(0.5, s_middle); //t=0.5
						_P_up = s_middle._O - Point<double>(R_for_integration/2., 0, 0.01); //The point is opposite to positive front orientation
						_P_down = s_middle._O - Point<double>(R_for_integration/2., 0, -0.01); //The point is opposite to positive front orientation
						P_up = s_middle.transfer_from_NEW_into_OLD(_P_up);
						P_down = s_middle.transfer_from_NEW_into_OLD(_P_down);
						double _len = 0;
						int k_up = -1, k_down = -1;
						k_up = grid.GetNearestElementID(P_up, _len);
						k_down = grid.GetNearestElementID(P_down, _len);
						Point<double> U_up, U_down;
						U_up = grid.GetSolutionInPoint(k_up, P_up, U) + P_up;
						U_down = grid.GetSolutionInPoint(k_up, P_up, U) + P_down;
						

						//solving delta in LOCAL system
						Point<double> _U_up = s_middle.transfer_from_OLD_into_NEW(U_up);
						Point<double> _U_down = s_middle.transfer_from_OLD_into_NEW(U_down);
						double _delta_I = _U_up.y - _U_down.y;
						double _delta_II = _U_up.x - _U_down.x;

						//solving SIF (p.698: J.W. Eischen "An improved method for computing the J2 integral")
						std::complex<double> _compl_K_I, _compl_K_II, _compl_in;
						double _in = 1 - (inJ[1] / inJ[0])*(inJ[1] / inJ[0]);
						_compl_in = _in;
						_compl_K_I = sqrt(_compl_in);
						_compl_K_II = _compl_K_I;
						K_I = _compl_K_I.real();
						K_II = _compl_K_II.real();
						if (abs(_delta_I) < abs(_delta_II))
						{
							K_I = 1 - K_I;//signum??
							K_II = 1 + K_II; //signum??
						}
						else {
							K_I = 1 + K_I;//signum??
							K_II = 1 - K_II;//signum??
						}
						double E = grid.GetDomain(0)->forMech.GetE(0);
						K_I = math::GetSignum(_delta_I)*sqrt(E*inJ[1] / 2 * K_I);
						K_II = math::GetSignum(_delta_II)*sqrt(E*inJ[1] / 2 * K_I);
					}
					//Cherepanov's criterion (с.108: В.И. Гришин, В.Ю. Донченко "Метод расчета коэффициентов интенсивности напряжений в элементах конструкций с криволинейными трещинами")
					double Qrad = atan(-2 * K_I*K_II / (K_I*K_I + K_II * K_II)); //in rad

					double Qrad_simple = atan(inJ[1] / inJ[0]); //in rad
					printf("FR[%d], Point<double>[%d], J1[%.2e], J2[%.2e], J3[%.2e], Qgrad_new[%.2lf], Qgrad_old[%.2lf]\n", Fr, segm - 1, inJ[0], inJ[1], inJ[2], 180 * Qrad / M_PI, 180 * Qrad_simple / M_PI);
					Qrad = Qrad_simple;
					
					//solving critical value
					double G = abs(inJ[0] * cos(Qrad) + inJ[1] * sin(Qrad));
					double K_Ic = (grid.GetDomain(0))->forMech.GetK_Ic(); // песчаник
					double Step_crack_front = R_for_integration;
					//if (G < K_Ic) Step_crack_front /= 10;
					printf("\tG[%.2e], K_Ic[%.2e], Step_crack_front[%.2lf]\n", G, K_Ic, Step_crack_front);

					double Qgrad = 180 * Qrad / M_PI; //in grad

					//t=0
					newCoordSystem s0;
					NewSystem(0, s0);
					Point<double> _R, R;
					_R.z = s0._O.z;
					_R.x = cos(Qrad) * R_for_integration;
					_R.y = sqrt(R_for_integration*R_for_integration - _R.x*_R.x) *math::GetSignum(Qrad);
					R = s0.transfer_from_NEW_into_OLD(_R);
					R -= s0.O; //it is a VECTOR therefore we have the minus O
					//printf("FR[%d], Point<double>[%d], J1[%.2e], J2[%.2e], J3[%.2e], Qgrad[%.2lf], R(%.2e, %.2e, %.2e)=%.2e\n", Fr, segm - 1, inJ[0], inJ[1], inJ[2], Qgrad, R.x, R.y, R.z, math::SolveLengthVector(R));
					//printf("FR[%d], Point<double>[%d], J1[%.2e], J2[%.2e], J3[%.2e], Qgrad[%.2lf], R=%.2e, _R=%.2e\n", Fr, segm - 1, inJ[0], inJ[1], inJ[2], Qgrad, math::SolveLengthVector(R), math::SolveLengthVector(_R));
					grid.cracks[Cr].propogation[Fr][segm - 1].vect_of_propogation.push_back(R);
					grid.cracks[Cr].propogation[Fr][segm - 1].length_of_propogation.push_back(Step_crack_front);
					grid.cracks[Cr].propogation[Fr][segm - 1].id_old_point = grid.cracks[Cr].id_points_in_fronts[Fr][segm];
					grid.cracks[Cr].propogation[Fr][segm - 1].base_coord_system.push_back(newCoordSystem(s0));
					grid.cracks[Cr].propogation[Fr][segm - 1].Qgrad.push_back(Qgrad);
					grid.cracks[Cr].propogation[Fr][segm - 1].G.push_back(G);

					//t=1
					newCoordSystem s1;
					NewSystem(1, s1);
					_R.z = s1._O.z;
					_R.x = cos(Qrad) * R_for_integration;
					_R.y = sqrt(R_for_integration*R_for_integration - _R.x*_R.x) *math::GetSignum(Qrad);
					R = s1.transfer_from_NEW_into_OLD(_R) - s1.O;
					//printf("FR[%d], Point<double>[%d], J1[%.2e], J2[%.2e], J3[%.2e], Qgrad[%.2lf], R(%.2e, %.2e, %.2e)=%.2e\n", Fr, segm, inJ[0], inJ[1], inJ[2], Qgrad, R.x, R.y, R.z, math::SolveLengthVector(R));
					grid.cracks[Cr].propogation[Fr][segm].vect_of_propogation.push_back(R);
					grid.cracks[Cr].propogation[Fr][segm].length_of_propogation.push_back(Step_crack_front);
					grid.cracks[Cr].propogation[Fr][segm].id_old_point = grid.cracks[Cr].id_points_in_fronts[Fr][segm + 1];
					grid.cracks[Cr].propogation[Fr][segm].base_coord_system.push_back(newCoordSystem(s1));
					grid.cracks[Cr].propogation[Fr][segm].Qgrad.push_back(Qgrad);
					grid.cracks[Cr].propogation[Fr][segm].G.push_back(G);
				}

				//create new cracks
				//new nodes
				int old_size_xyz = grid.cracks[Cr].xyz.size();
				if (false) {
					for (int p = 0; p < grid.cracks[Cr].propogation[Fr].size(); p++)
					{
						Point<double> middle_vect;
						for (int v = 0; v < grid.cracks[Cr].propogation[Fr][p].vect_of_propogation.size(); v++)
						{
							middle_vect += grid.cracks[Cr].propogation[Fr][p].vect_of_propogation[v] / grid.cracks[Cr].propogation[Fr][p].vect_of_propogation.size();
						}
						middle_vect = middle_vect / math::SolveLengthVector(middle_vect) * grid.cracks[Cr].R;
						//printf("Len(FR[%d], Point[%d]) = %.2e\n", Fr, p, math::SolveLengthVector(middle_vect));
						grid.cracks[Cr].xyz.push_back(grid.cracks[Cr].xyz[grid.cracks[Cr].propogation[Fr][p].id_old_point] + middle_vect);
						grid.cracks[Cr].propogation[Fr][p].id_new_point = grid.cracks[Cr].xyz.size() - 1;
					}
				}

				//FULL MIDDLE MOVE!!!
				if (false)
				{
					Point<double> FULL_MIDDLE_VECTOR;
					for (int p = 0; p < false && grid.cracks[Cr].propogation[Fr].size(); p++)
					{
						Point<double> middle_vect;
						for (int v = 0; v < grid.cracks[Cr].propogation[Fr][p].vect_of_propogation.size(); v++)
						{
							middle_vect += grid.cracks[Cr].propogation[Fr][p].vect_of_propogation[v] / grid.cracks[Cr].propogation[Fr][p].vect_of_propogation.size();
						}
						middle_vect = middle_vect / math::SolveLengthVector(middle_vect) * grid.cracks[Cr].R;
						FULL_MIDDLE_VECTOR += middle_vect;
					}
					FULL_MIDDLE_VECTOR /= grid.cracks[Cr].propogation[Fr].size();
					for (int p = 0; p < grid.cracks[Cr].propogation[Fr].size(); p++)
					{
						grid.cracks[Cr].xyz.push_back(grid.cracks[Cr].xyz[grid.cracks[Cr].propogation[Fr][p].id_old_point] + FULL_MIDDLE_VECTOR);
						grid.cracks[Cr].propogation[Fr][p].id_new_point = grid.cracks[Cr].xyz.size() - 1;
					}
				}

				//smoothing via Ns nodes
				int Ns = grid.cracks[Cr].smoothing_coefficient;
				int end_of_update = grid.cracks[Cr].propogation[Fr].size();
				if(grid.cracks[Cr].propogation[Fr][0].id_old_point == grid.cracks[Cr].propogation[Fr][grid.cracks[Cr].propogation[Fr].size()-1].id_old_point)
				{
					end_of_update--;
				}
				for (int id_base_vertex = 0; id_base_vertex < end_of_update; id_base_vertex++)
				{
					int start_id = id_base_vertex;
					if (id_base_vertex + Ns >= grid.cracks[Cr].propogation[Fr].size())
					{
						start_id = grid.cracks[Cr].propogation[Fr].size() - Ns;
						if (start_id < 0) start_id = 0;
					}

					double middle_Qgrad = 0;
					double middle_G = 0;
					int n_i = 0;

					for (int i = 0; i < Ns && (start_id + i) < grid.cracks[Cr].propogation[Fr].size(); i++)
					{
						for (int q = 0; q < grid.cracks[Cr].propogation[Fr][start_id + i].Qgrad.size(); q++)
						{
							middle_Qgrad += grid.cracks[Cr].propogation[Fr][start_id + i].Qgrad[q];
							middle_G += grid.cracks[Cr].propogation[Fr][start_id + i].G[q];
							n_i++;
						}
					}
					middle_Qgrad /= n_i;
					middle_G /= n_i;

					Point<double> propogation_vector;
					Point<double> simple_vector;
					double simpleQ = 0;
					double simpleG = 0;
					double length = 0;
					for (int q = 0; q < grid.cracks[Cr].propogation[Fr][id_base_vertex].base_coord_system.size(); q++)
					{
						propogation_vector += grid.cracks[Cr].propogation[Fr][id_base_vertex].base_coord_system[q].MakePropogation(middle_Qgrad, grid.cracks[Cr].R) /
							grid.cracks[Cr].propogation[Fr][id_base_vertex].base_coord_system.size();
						simple_vector += grid.cracks[Cr].propogation[Fr][id_base_vertex].vect_of_propogation[q] /
							grid.cracks[Cr].propogation[Fr][id_base_vertex].base_coord_system.size();
						simpleQ += grid.cracks[Cr].propogation[Fr][id_base_vertex].Qgrad[q] /
							grid.cracks[Cr].propogation[Fr][id_base_vertex].base_coord_system.size();
						simpleG += grid.cracks[Cr].propogation[Fr][id_base_vertex].G[q] /
							grid.cracks[Cr].propogation[Fr][id_base_vertex].base_coord_system.size();
						length += grid.cracks[Cr].propogation[Fr][id_base_vertex].length_of_propogation[q] /
							grid.cracks[Cr].propogation[Fr][id_base_vertex].base_coord_system.size();
					}
					if (middle_G < (grid.GetDomain(0))->forMech.GetK_Ic())
						length /= 10;

					propogation_vector = propogation_vector / math::SolveLengthVector(propogation_vector) * length /*grid.cracks[Cr].R*/;
					simple_vector = simple_vector / math::SolveLengthVector(simple_vector) * grid.cracks[Cr].R;

					/*printf_s("Front[%d], old_id[%d], Q=%.2lf; OldPoint(%.2lf, %.2lf, %.2lf); Vector=(%.2lf, %.2lf, %.2lf); SimpleQ=%.2lf; SimpleVector=(%.2lf, %.2lf, %.2lf)\n",
						Fr,
						grid.cracks[Cr].propogation[Fr][id_base_vertex].id_old_point,
						middle_Qgrad,
						grid.cracks[Cr].xyz[grid.cracks[Cr].propogation[Fr][id_base_vertex].id_old_point].x,
						grid.cracks[Cr].xyz[grid.cracks[Cr].propogation[Fr][id_base_vertex].id_old_point].y,
						grid.cracks[Cr].xyz[grid.cracks[Cr].propogation[Fr][id_base_vertex].id_old_point].z,
						propogation_vector.x,
						propogation_vector.y,
						propogation_vector.z,
						simpleQ,
						simple_vector.x,
						simple_vector.y,
						simple_vector.z);*/

					grid.cracks[Cr].xyz.push_back(grid.cracks[Cr].xyz[grid.cracks[Cr].propogation[Fr][id_base_vertex].id_old_point] + propogation_vector);
					grid.cracks[Cr].propogation[Fr][id_base_vertex].id_new_point = grid.cracks[Cr].xyz.size() - 1;
					grid.cracks[Cr].propogation[Fr][id_base_vertex]._vect_of_propogation = propogation_vector;
					grid.cracks[Cr].propogation[Fr][id_base_vertex]._length_of_propogation = length;
					grid.cracks[Cr].propogation[Fr][id_base_vertex]._Qgrad = middle_Qgrad;
					grid.cracks[Cr].propogation[Fr][id_base_vertex]._G = middle_G;
				}
				if (grid.cracks[Cr].propogation[Fr][0].id_old_point == grid.cracks[Cr].propogation[Fr][grid.cracks[Cr].propogation[Fr].size() - 1].id_old_point)
				{
					grid.cracks[Cr].propogation[Fr][grid.cracks[Cr].propogation[Fr].size() - 1].id_new_point = grid.cracks[Cr].propogation[Fr][0].id_new_point;
				}


				//new triangles and fronts
				Point<double> true_normal = grid.cracks[Cr].triangles[0].GetNormal();
				for (int segm = 0; segm < grid.cracks[Cr].fronts[Fr].size(); segm++)
				{
					Point<double> new_normal;
					int cur_segm[3] = { grid.cracks[Cr].fronts[Fr][segm].id_base_triangles, grid.cracks[Cr].fronts[Fr][segm].id_left, grid.cracks[Cr].fronts[Fr][segm].id_right };
					int new_segm[3];
					for (int ii = 1; ii < 3; ii++)
					{
						for (int p = 0; p < grid.cracks[Cr].propogation[Fr].size(); p++)
						{
							if (cur_segm[ii] == grid.cracks[Cr].propogation[Fr][p].id_old_point)
							{
								new_segm[ii] = grid.cracks[Cr].propogation[Fr][p].id_new_point;
								break;
							}
						}
					}

					for (int ii = 0; ii < grid.cracks[Cr].triangles.size(); ii++)
					{
						for (int jj = 0; jj < grid.cracks[Cr].triangles[ii].GetNodesCount(); jj++)
						{
							grid.cracks[Cr].triangles[ii].SetNode(jj, &grid.cracks[Cr].xyz[grid.cracks[Cr].triangles[ii].GetIdNode(jj)]);
						}
					}

					grid.cracks[Cr].triangles.push_back(geometry::Triangle());
					grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].SetIdNode(0, cur_segm[1]);
					grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].SetIdNode(1, cur_segm[2]);
					grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].SetIdNode(2, new_segm[2]);
					grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].SetNode(0, &grid.cracks[Cr].xyz[grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].GetIdNode(0)]);
					grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].SetNode(1, &grid.cracks[Cr].xyz[grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].GetIdNode(1)]);
					grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].SetNode(2, &grid.cracks[Cr].xyz[grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].GetIdNode(2)]);
					new_normal = grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].GetNormal();
					if (true_normal*new_normal < 0)
					{
						grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].SetIdNode(0, cur_segm[2]);
						grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].SetIdNode(1, cur_segm[1]);
						grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].SetIdNode(2, new_segm[2]);
						grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].SetNode(0, &grid.cracks[Cr].xyz[grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].GetIdNode(0)]);
						grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].SetNode(1, &grid.cracks[Cr].xyz[grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].GetIdNode(1)]);
						grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].SetNode(2, &grid.cracks[Cr].xyz[grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].GetIdNode(2)]);
					}

					grid.cracks[Cr].triangles.push_back(geometry::Triangle());
					grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].SetIdNode(0, cur_segm[1]);
					grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].SetIdNode(1, new_segm[1]);
					grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].SetIdNode(2, new_segm[2]);
					grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].SetNode(0, &grid.cracks[Cr].xyz[grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].GetIdNode(0)]);
					grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].SetNode(1, &grid.cracks[Cr].xyz[grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].GetIdNode(1)]);
					grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].SetNode(2, &grid.cracks[Cr].xyz[grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].GetIdNode(2)]);
					new_normal = grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].GetNormal();
					if (true_normal*new_normal < 0)
					{
						grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].SetIdNode(0, cur_segm[1]);
						grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].SetIdNode(1, new_segm[2]);
						grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].SetIdNode(2, new_segm[1]);
						grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].SetNode(0, &grid.cracks[Cr].xyz[grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].GetIdNode(0)]);
						grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].SetNode(1, &grid.cracks[Cr].xyz[grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].GetIdNode(1)]);
						grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].SetNode(2, &grid.cracks[Cr].xyz[grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].GetIdNode(2)]);
					}

					new_segm[0] = grid.cracks[Cr].triangles.size() - 1;
					grid.cracks[Cr].fronts[Fr][segm].id_base_triangles = new_segm[0];
					grid.cracks[Cr].fronts[Fr][segm].id_left = new_segm[1];
					grid.cracks[Cr].fronts[Fr][segm].id_right = new_segm[2];
				}
			}
		}

		FILE* fout;
		char new_crack[1000];
		for (int Cr = 0; Cr < grid.cracks.size(); Cr++)
		{
			//sprintf_s(new_crack, "%s/Crack_%d_step%d.txt", dir, Cr, curr_STEP + 1);

			sprintf_s(new_crack, "%s/Crack_%d_properties.txt", dir, Cr, curr_STEP + 1);
			fopen_s(&fout, new_crack, "w");
			for (int Fr = 0; Fr < grid.cracks[Cr].propogation.size(); Fr++)
			{
				fprintf_s(fout, "Front %d\n", Fr);
				for (int i = 0; i < grid.cracks[Cr].propogation[Fr].size(); i++)
				{
					fprintf_s(fout, "point[%d]: id_old_point[%d], id_new_point[%d], vect(%.2lf, %.2lf, %.2lf), len[%.2lf], Qgrad[%.2lf], G[%.2e]\n", 
						i, 
						grid.cracks[Cr].propogation[Fr][i].id_old_point,
						grid.cracks[Cr].propogation[Fr][i].id_new_point, 
						grid.cracks[Cr].propogation[Fr][i]._vect_of_propogation.x, grid.cracks[Cr].propogation[Fr][i]._vect_of_propogation.y, grid.cracks[Cr].propogation[Fr][i]._vect_of_propogation.z,
						grid.cracks[Cr].propogation[Fr][i]._length_of_propogation,
						grid.cracks[Cr].propogation[Fr][i]._Qgrad,
						grid.cracks[Cr].propogation[Fr][i]._G);
				}
				fprintf_s(fout, "\n");
			}
			fclose(fout);

			sprintf_s(new_crack, "%s/Crack_%d_%d.txt", dir, Cr, curr_STEP + 1);
			fopen_s(&fout, new_crack, "w");
			fprintf_s(fout, "%d\n", (int)grid.cracks[Cr].xyz.size());
			for (int i = 0; i < grid.cracks[Cr].xyz.size(); i++)
			{
				fprintf_s(fout, "%.6e\t%.6e\t%.6e\n",
					grid.cracks[Cr].xyz[i].x,
					grid.cracks[Cr].xyz[i].y,
					grid.cracks[Cr].xyz[i].z);
			}

			fprintf_s(fout, "\n%d\n", (int)grid.cracks[Cr].triangles.size());
			for (int i = 0; i < grid.cracks[Cr].triangles.size(); i++)
			{
				fprintf_s(fout, "%d\t%d\t%d\n",
					grid.cracks[Cr].triangles[i].GetIdNode(0),
					grid.cracks[Cr].triangles[i].GetIdNode(1),
					grid.cracks[Cr].triangles[i].GetIdNode(2));
			}

			fprintf_s(fout, "\n%d\n", (int)grid.cracks[Cr].fronts.size());
			for (int i = 0; i < grid.cracks[Cr].fronts.size(); i++)
			{
				fprintf_s(fout, "%d\n", (int)grid.cracks[Cr].fronts[i].size());
				for (int j = 0; j < grid.cracks[Cr].fronts[i].size(); j++)
				{
					fprintf_s(fout, "%d\t%d\t%d\n",
						grid.cracks[Cr].fronts[i][j].id_base_triangles,
						grid.cracks[Cr].fronts[i][j].id_left,
						grid.cracks[Cr].fronts[i][j].id_right);
				}
			}
			fclose(fout);

			{ //tech_plot - as middle plane
				sprintf_s(new_crack, "%s/Crack_%d_step%d.dat", dir, Cr, curr_STEP + 1);
				fopen_s(&fout, new_crack, "w");

				fprintf_s(fout, "TITLE     = \"numerical\"\n");
				fprintf_s(fout, "VARIABLES = \"x\"\n \"y\"\n \"z\"\n");
				fprintf_s(fout, "\"crack\"\n");
				fprintf_s(fout, "ZONE T=\"Crack_%d_step%d\"\n", Cr, curr_STEP + 1);
				fprintf_s(fout, " N=%d,  E=%d, F=FEBLOCK ET=Triangle \n", (int)grid.cracks[Cr].xyz.size(), (int)grid.cracks[Cr].triangles.size());
				fprintf_s(fout, " VARLOCATION=(NODAL NODAL NODAL CELLCENTERED)\n");

				for (int i = 0; i < grid.cracks[Cr].xyz.size(); i++)
					fprintf_s(fout, "%.10e\n", grid.cracks[Cr].xyz[i].x);
				fprintf_s(fout, "\n");
				for (int i = 0; i < grid.cracks[Cr].xyz.size(); i++)
					fprintf_s(fout, "%.10e\n", grid.cracks[Cr].xyz[i].y);
				fprintf_s(fout, "\n");
				for (int i = 0; i < grid.cracks[Cr].xyz.size(); i++)
					fprintf_s(fout, "%.10e\n", grid.cracks[Cr].xyz[i].z);
				fprintf_s(fout, "\n");

				for (int i = 0; i < grid.cracks[Cr].triangles.size(); i++)
					fprintf_s(fout, "%d\n", curr_STEP + 1);
				fprintf_s(fout, "\n");

				for (int i = 0; i < grid.cracks[Cr].triangles.size(); i++)
				{
					for (int j = 0; j < 3; j++)
						fprintf_s(fout, "%d ", grid.cracks[Cr].triangles[i].GetIdNode(j) + 1);
					fprintf_s(fout, "\n");
				}

				fclose(fout);
			}

			{ //tech_plot - as face planes
				std::vector<Point<double>> xyz_up(grid.cracks[Cr].xyz.size()), xyz_down(grid.cracks[Cr].xyz.size());
				std::vector<int> elem_into_nodes(grid.cracks[Cr].xyz.size());
				double MaxSize = 0, MinSize = 1e+15;

				for (int id = 0; id < grid.cracks[Cr].triangles.size(); id++)
				{
					Point<double> normal = grid.cracks[Cr].triangles[id].GetNormal();
					normal /= math::SolveLengthVector(normal);
					Point<double> O = grid.cracks[Cr].triangles[id].GetWeightCentr();
					double move_koef = 0.5;

					Point<double> up_X = O + normal * move_koef;
					double len;
					int up_elem = grid.GetNearestElementID(up_X, len);
					Point<double> move_vector_UP = grid.GetSolutionInPoint(up_elem, up_X, U);

					Point<double> down_X = O - normal * move_koef;
					int down_elem = grid.GetNearestElementID(down_X, len);
					Point<double> move_vector_DOWN = grid.GetSolutionInPoint(down_elem, down_X, U);

					double move_size = math::SolveLengthVector(move_vector_DOWN) + math::SolveLengthVector(move_vector_UP);
					if (move_size > MaxSize) MaxSize = move_size;
					if (move_size < MinSize) MinSize = move_size;

					for (int _i = 0; _i < 3; _i++)
					{
						elem_into_nodes[grid.cracks[Cr].triangles[id].GetIdNode(_i)] ++;
						xyz_up[grid.cracks[Cr].triangles[id].GetIdNode(_i)] += move_vector_UP;
						xyz_down[grid.cracks[Cr].triangles[id].GetIdNode(_i)] += move_vector_DOWN;
					}
				}
				for (int i = 0; i < grid.cracks[Cr].fronts.size(); i++)
				{
					for (int j = 0; j < grid.cracks[Cr].fronts[i].size(); j++)
					{
						xyz_up[grid.cracks[Cr].fronts[i][j].id_left] = Point<double>(0, 0, 0);
						xyz_up[grid.cracks[Cr].fronts[i][j].id_right] = Point<double>(0, 0, 0);

						xyz_down[grid.cracks[Cr].fronts[i][j].id_left] = Point<double>(0, 0, 0);
						xyz_down[grid.cracks[Cr].fronts[i][j].id_right] = Point<double>(0, 0, 0);
					}
				}


				for (int id = 0; id < xyz_up.size(); id++)
				{
					xyz_up[id] /= elem_into_nodes[id];
					xyz_down[id] /= elem_into_nodes[id];
					xyz_up[id] += grid.cracks[Cr].xyz[id];
					xyz_down[id] += grid.cracks[Cr].xyz[id];
				}



				sprintf_s(new_crack, "%s/Crack_%d_step%d_in3D_MaxSize_%.2e_MinSize_%.2e.dat", dir, Cr, curr_STEP + 1, MaxSize, MinSize);
				fopen_s(&fout, new_crack, "w");

				fprintf_s(fout, "TITLE     = \"numerical\"\n");
				fprintf_s(fout, "VARIABLES = \"x\"\n \"y\"\n \"z\"\n");
				fprintf_s(fout, "\"crack\"\n");
				fprintf_s(fout, "ZONE T=\"Crack_%d_step%d\"\n", Cr, curr_STEP + 1);
				fprintf_s(fout, " N=%d,  E=%d, F=FEBLOCK ET=Triangle \n", xyz_up.size() * 2, grid.cracks[Cr].triangles.size() * 2);
				fprintf_s(fout, " VARLOCATION=(NODAL NODAL NODAL CELLCENTERED)\n");

				for (int i = 0; i < xyz_up.size(); i++)
					fprintf_s(fout, "%.10e\n", xyz_up[i].x);
				for (int i = 0; i < xyz_down.size(); i++)
					fprintf_s(fout, "%.10e\n", xyz_down[i].x);
				fprintf_s(fout, "\n");
				for (int i = 0; i < xyz_up.size(); i++)
					fprintf_s(fout, "%.10e\n", xyz_up[i].y);
				for (int i = 0; i < xyz_down.size(); i++)
					fprintf_s(fout, "%.10e\n", xyz_down[i].y);
				fprintf_s(fout, "\n");
				for (int i = 0; i < xyz_up.size(); i++)
					fprintf_s(fout, "%.10e\n", xyz_up[i].z);
				for (int i = 0; i < xyz_down.size(); i++)
					fprintf_s(fout, "%.10e\n", xyz_down[i].z);
				fprintf_s(fout, "\n");


				for (int i = 0; i < grid.cracks[Cr].triangles.size() * 2; i++)
					fprintf_s(fout, "%d\n", curr_STEP + 1);
				fprintf_s(fout, "\n");

				for (int i = 0; i < grid.cracks[Cr].triangles.size(); i++)
				{
					for (int j = 0; j < 3; j++)
						fprintf_s(fout, "%d ", grid.cracks[Cr].triangles[i].GetIdNode(j) + 1);
					fprintf_s(fout, "\n");
				}
				for (int i = 0; i < grid.cracks[Cr].triangles.size(); i++)
				{
					for (int j = 0; j < 3; j++)
						fprintf_s(fout, "%d ", grid.cracks[Cr].xyz.size() + grid.cracks[Cr].triangles[i].GetIdNode(j) + 1);
					fprintf_s(fout, "\n");
				}

				fclose(fout);
			}
		}
	}
	void CrackPropagation_3D_Cherepanov_v2(MultiXFEM::Grid &grid, std::vector<Point<double>> &U, char *dir, int curr_STEP)
	{
		//to find the propagation vector for each a front point
		for (int Cr = 0; Cr < grid.cracks.size(); Cr++)
		{
			double R_for_integration = grid.cracks[Cr].R;
			for (int Fr = 0; Fr < grid.cracks[Cr].id_points_in_fronts.size(); Fr++)
			{
				grid.cracks[Cr].propogation.resize(grid.cracks[Cr].id_points_in_fronts.size());
				grid.cracks[Cr].propogation[Fr].resize(grid.cracks[Cr].id_points_in_fronts[Fr].size() - 2);
				for (int segm = 1; segm < grid.cracks[Cr].id_points_in_fronts[Fr].size() - 2; segm++)
				{
					Point<double> P0, P1, P2, P3;
					P0 = grid.cracks[Cr].xyz[grid.cracks[Cr].id_points_in_fronts[Fr][segm - 1]];
					P1 = grid.cracks[Cr].xyz[grid.cracks[Cr].id_points_in_fronts[Fr][segm]];
					P2 = grid.cracks[Cr].xyz[grid.cracks[Cr].id_points_in_fronts[Fr][segm + 1]];
					P3 = grid.cracks[Cr].xyz[grid.cracks[Cr].id_points_in_fronts[Fr][segm + 2]];

					if (segm - 1 == 0 && grid.cracks[Cr].id_points_in_fronts[Fr][segm - 1] == grid.cracks[Cr].id_points_in_fronts[Fr][segm])
					{
						P0 = P1 + (P1 - P2);
					}
					if ((segm + 2) == (grid.cracks[Cr].id_points_in_fronts[Fr].size() - 1) &&
						grid.cracks[Cr].id_points_in_fronts[Fr][segm + 1] == grid.cracks[Cr].id_points_in_fronts[Fr][segm + 2])
					{
						P3 = P2 + (P2 - P1);
					}

					Point<double> r3 = P1 * 3 - P2 * 3 + P3 + P0;
					Point<double> r2 = P1 * (-5) + P2 * 4 - P3;
					Point<double> r1 = P0 * (-1) + P2;
					Point<double> r0 = P1 * 2;
					r3 = P0 * (-1) + P1 * 3 - P2 * 3 + P3;
					r2 = P0 * 2 + P1 * (-5) + P2 * 4 - P3;
					r1 = P0 * (-1) + P2;
					r0 = P1 * 2;
					auto F = [&](double t) -> Point<double>
					{
						return (r3*pow(t, 3) + r2 * pow(t, 2) + r1 * t + r0) / 2.;
					};

					Point<double> _r2 = P1 * 9 - P2 * 9 + P3 * 3 + P0 * 3;
					Point<double> _r1 = P1 * (-10) + P2 * 8 - P3 * 2;
					Point<double> _r0 = P0 * (-1) + P2;
					_r2 = r3 * 3.;
					_r1 = r2 * 2.;
					_r0 = r1;
					auto dF_dt = [&](double t) -> Point<double>
					{
						return (_r2*pow(t, 2) + _r1 * t + _r0) / 2.;
					};

					struct newCoordSystem {
						std::vector<Point<double>> _X; //new basis in OLD coord
						Point<double> O; //New coord centr in OLD coord
						Point<double> _O; //New coord centr in NEW coord
						Tensor2Rank3D A; //transfer from Old into New
						Tensor2Rank3D _A; //transfer from New into Old
						Point<double> Offset; //offset from Old into New
						Point<double> _Offset; //offset from New into Old 

						Point<double> transfer_from_OLD_into_NEW(Point<double> &T)
						{
							return math::MakeTransferPoint(T, A, Offset);
						}
						Point<double> transfer_from_NEW_into_OLD(Point<double> &_T)
						{
							return math::MakeTransferPoint(_T, _A, _Offset);
						}
						Tensor2Rank3D transfer_from_OLD_into_NEW(Tensor2Rank3D &T)
						{
							return math::MakeTransferTensor(T, A, Offset);
						}
					};
					auto NewSystem = [Cr, Fr, &grid, &dF_dt, &F](double t, newCoordSystem &System) -> void
					{
						System._X.resize(3);
						int curr_segm;
						System._X[2] = dF_dt(t);
						curr_segm = 0;

						Point<double> Line[2];
						Line[0] = F(t);
						Line[1] = F(t) + System._X[2];
						Point<double> test;
						for (int tt = 0; tt < 3; tt++)
						{
							if (grid.cracks[Cr].fronts[Fr][curr_segm].id_left !=
								grid.cracks[Cr].triangles[grid.cracks[Cr].fronts[Fr][curr_segm].id_base_triangles].GetIdNode(tt)
								&& grid.cracks[Cr].fronts[Fr][curr_segm].id_right !=
								grid.cracks[Cr].triangles[grid.cracks[Cr].fronts[Fr][curr_segm].id_base_triangles].GetIdNode(tt))
							{
								test = grid.cracks[Cr].xyz[grid.cracks[Cr].triangles[grid.cracks[Cr].fronts[Fr][curr_segm].id_base_triangles].GetIdNode(tt)];
								break;
							}
						}
						DenseMatrix<double, double> test_matrix;
						System._X[0] = math::MakeNormalForLine(Line, test, test_matrix);
						System._X[1] = grid.cracks[Cr].triangles[grid.cracks[Cr].fronts[Fr][curr_segm].id_base_triangles].GetNormal();

						System._X[0] /= math::SolveLengthVector(System._X[0]);
						System._X[1] /= math::SolveLengthVector(System._X[1]);
						System._X[2] /= math::SolveLengthVector(System._X[2]);

						System.A.val[0][0] = System._X[0].x; System.A.val[0][1] = System._X[0].y; System.A.val[0][2] = System._X[0].z;
						System.A.val[1][0] = System._X[1].x; System.A.val[1][1] = System._X[1].y; System.A.val[1][2] = System._X[1].z;
						System.A.val[2][0] = System._X[2].x; System.A.val[2][1] = System._X[2].y; System.A.val[2][2] = System._X[2].z;
						System._A = math::SolveInverseMatrix3x3(System.A);

						Point<double> _p;
						System.Offset = F(t)*(-1.0);
						System._Offset = math::MakeTransferPoint(_p, System.A, System.Offset)*(-1);

						Point<double> _f;
						_f = F(t);
						System.O = F(t);
						System._O = System.transfer_from_OLD_into_NEW(_f);
						System.O = System.transfer_from_NEW_into_OLD(System._O);
					};

					//integration in IxJxK via NEW СOORDINATES
					std::vector<double> inJ(3);
					int k = 0;
					int sizeI = 2;
					int sizeJ = 2;
					int sizeK = 2;
					int stepsOfIntegr_K = 2;
					double step_t_viaK = 1.0 / (sizeK);
					auto f = [](Point<double> X) ->double {return 1; };
					for (int K = 0; K < sizeK; K++)
					{
						std::vector<std::vector<std::vector<double>>> data_store_for_integrals_in_K(3);
						for (int k = 0; k < 3; k++)
						{
							data_store_for_integrals_in_K[k].resize(sizeJ);
							for (int j = 0; j < data_store_for_integrals_in_K[k].size(); j++)
							{
								data_store_for_integrals_in_K[k][j].resize(sizeI);
							}
						}

						std::vector < double > xj_integr1D, q_integr1D;
						geometry::Rectangle tmp_elem;
						tmp_elem.SetIntegrationLaw(stepsOfIntegr_K, xj_integr1D, q_integr1D);

						double current_t = 0;
						double start_t = K * step_t_viaK;
						double end_t = (K + 1)*step_t_viaK;
						for (int T = 0; T < stepsOfIntegr_K; T++)
						{
							double _currentZ = 0.0;
							double current_t = (start_t + end_t + xj_integr1D[T] * (end_t - start_t)) / 2.0;
							newCoordSystem System_for_t;
							NewSystem(current_t, System_for_t);

							Point<double> _Up(R_for_integration, R_for_integration, 0.0), _Down(-R_for_integration, -R_for_integration, 0.0);
							Point<double> _H((_Up.x - _Down.x) / sizeI, (_Up.y - _Down.y) / sizeJ, 0.0);
							for (int J = 0; J < sizeJ; J++)
							{
								for (int I = 0; I < sizeI; I++)
								{
									for (int k = 0; k < 3; k++)
									{
										auto Func_q = [&grid, &U, &_Up, &_Down](Point<double> &X, newCoordSystem &System_for_t, int k) ->double
										{
											double res = 0;
											Point<double> _X = System_for_t.transfer_from_OLD_into_NEW(X);

											int curr_el = grid.GetElementID(X);
											if (curr_el == -1) return 0.0;
											auto current_element = grid.GetElement(curr_el);


											/*if (current_element->temperature_loc.size() != 0)
											{
												std::vector<std::vector<double>> D_loc;
												grid.domain[current_element->get_id_domain()].forMech.get_D(current_element->get_temperature(X), 3, D_loc);
											}*/

											/*std::vector<std::vector<double>> SIG, EPS, _SIG, _EPS, dU_dX;
											dU_dX = current_element->get_value_dU_dX_mech(X, value_in_DOF);
											EPS = current_element->get_value_EPS_mech(X, value_in_DOF);
											SIG = current_element->get_value_SIGMA_mech(EPS, grid.domain[current_element->get_id_domain()]);*/

											Point<Point<double>> dU = grid.GetDerevativeFromSolutionInPoint(curr_el, X, U);
											Tensor2Rank3D _dU_old, _dU_new;
											Point<Point<double>> dU_new;
											_dU_old.val[0][0] = dU.x.x; _dU_old.val[0][1] = dU.x.y; _dU_old.val[0][2] = dU.x.z;
											_dU_old.val[1][0] = dU.y.x; _dU_old.val[1][1] = dU.y.y; _dU_old.val[1][2] = dU.y.z;
											_dU_old.val[2][0] = dU.z.x; _dU_old.val[2][1] = dU.z.y; _dU_old.val[2][2] = dU.z.z;
											for (int i = 0; i < 3; i++)
											{
												for (int j = 0; j < 3; j++)
												{
													for (int m = 0; m < 3; m++)
													{
														_dU_new.val[i][j] += _dU_old.val[i][m] * System_for_t._A.val[m][j];
													}
												}
											}
											dU_new.x.x = _dU_new.val[0][0]; dU_new.x.y = _dU_new.val[0][1]; dU_new.x.z = _dU_new.val[0][2];
											dU_new.y.x = _dU_new.val[1][0]; dU_new.y.y = _dU_new.val[1][1]; dU_new.y.z = _dU_new.val[1][2];
											dU_new.z.x = _dU_new.val[2][0]; dU_new.z.y = _dU_new.val[2][1]; dU_new.z.z = _dU_new.val[2][2];

											//Tensor2Rank3D _EPS = grid.GetStressTensorFromSolutionInPoint(curr_el, dU_new);
											//Tensor2Rank3D _SIG = grid.GetStrainTensorFromSolutionInPoint(curr_el, X, _EPS);

											Tensor2Rank3D EPS = grid.GetStressTensorFromSolutionInPoint(curr_el, dU);
											Tensor2Rank3D SIG = grid.GetStrainTensorFromSolutionInPoint(curr_el, X, EPS);

											auto W = [](double ksi, double ksi_max, double ksi_min, int atr) -> double
											{
												double res;
												switch (atr)
												{
												case 0:
													res = (ksi_max - ksi) / (ksi_max - ksi_min);
													break;
												case 1:
													res = (ksi - ksi_min) / (ksi_max - ksi_min);
													break;
												default:
													res = 0;
													break;
												}
												return res;
											};
											auto dW = [](double ksi, double ksi_max, double ksi_min, int atr) -> double
											{
												double res;
												switch (atr)
												{
												case 0:
													res = -1 / (ksi_max - ksi_min);
													break;
												case 1:
													res = 1 / (ksi_max - ksi_min);
													break;
												default:
													res = 0;
													break;
												}
												return res;
											};
											auto dq = [&](Point<double> &_X, Point<double> &_o, Point<double> &_up, Point<double> &_down) -> std::vector<double>
											{
												std::vector<double> dq_dxi(3);

												//quadr
												int base_func[2];
												Point<double> _up_loc, _down_loc;
												bool find = false;
												if (!find && _X.x <= _o.x && _X.y <= _o.y)//III
												{
													base_func[0] = 1;
													base_func[1] = 1;
													_down_loc = _down;
													_up_loc = _o;
													find = true;
												}
												if (!find && _X.x >= _o.x && _X.y <= _o.y)//IV
												{
													base_func[0] = 0;
													base_func[1] = 1;
													_down_loc = Point<double>(_o.x, _down.y, _o.z);
													_up_loc = Point<double>(_up.x, _o.y, _o.z);
													find = true;
												}
												if (!find && _X.x <= _o.x && _X.y >= _o.y)//II
												{
													base_func[0] = 1;
													base_func[1] = 0;
													_down_loc = Point<double>(_down.x, _o.y, _o.z);
													_up_loc = Point<double>(_o.x, _up.y, _o.z);
													find = true;
												}
												if (!find && _X.x >= _o.x && _X.y >= _o.y) //I
												{
													base_func[0] = 0;
													base_func[1] = 0;
													_down_loc = _o;
													_up_loc = _up;
													find = true;
												}

												dq_dxi[0] = dW(_X.x, _up_loc.x, _down_loc.x, base_func[0]) * W(_X.y, _up_loc.y, _down_loc.y, base_func[1]);
												dq_dxi[1] = W(_X.x, _up_loc.x, _down_loc.x, base_func[0]) * dW(_X.y, _up_loc.y, _down_loc.y, base_func[1]);
												dq_dxi[2] = 0; // W(_X.x, _up_loc.x, _down_loc.x, base_func[0])*W(_X.y, _up_loc.y, _down_loc.y, base_func[1]);

												return dq_dxi;
											};

											//left part
											std::vector<double> _left_res(3);
											double _w = 0;
											double _w_kj = 0;
											for (int i = 0; i < 3; i++)
											{
												for (int j = 0; j < 3; j++)
												{
													_w += (SIG.val[i][j] * EPS.val[i][j]) / 2.;
												}
												_w_kj = (SIG.val[i][k] * EPS.val[i][k]) / 2.;
											}
											_left_res[k] = _w;
											_left_res[k] = _w_kj;

											//right_part
											std::vector<double> _right_res(3);
											for (int j = 0; j < 3; j++)
											{
												for (int i = 0; i < 3; i++)
												{
													double dui_dxk = 0;
													for (int m = 0; m < 3; m++)
													{
														dui_dxk += _dU_old.val[i][m] * System_for_t._A.val[m][k];
													}
													_right_res[j] += SIG.val[i][j] * dui_dxk;
												}
											}

											//q
											std::vector<double> dq_dxj = dq(_X, System_for_t._O, _Up, _Down);
											//if(k == 0) printf("%.5lf %.5lf %.5lf %.5lf\n", _X.x, _X.y, _X.z, dq_dxi[2]);

											res = 0;
											for (int j = 0; j < 3; j++)
											{
												res += -1 * (_left_res[j] - _right_res[j])*dq_dxj[j]; //??
												//res += -1 * (left_res - _right_res[j])*dq_dxj[j]; //??
											}

											//printf("elem[%d/%d/%d](%d): left_res=%.5lf right_res=%.5lf\n", I, J, K, k, left_res, right_res);
											return res;
										};

										Point<double> _up(_Down.x + _H.x*(I + 1), _Down.y + _H.y*(J + 1), 0.0);
										Point<double> _down(_Down.x + _H.x*(I), _Down.y + _H.y*(J), 0.0);
										double temp_res = 0;
										for (int i = 0; i < stepsOfIntegr_K; ++i)
										{
											double _currentX = (_down.x + _up.x + xj_integr1D[i] * (_up.x - _down.x)) / 2.0;
											for (int j = 0; j < stepsOfIntegr_K; ++j)
											{
												double _currentY = (_down.y + _up.y + xj_integr1D[j] * (_up.y - _down.y)) / 2.0;

												Point<double> _p(_currentX, _currentY, _currentZ);
												Point<double> X = System_for_t.transfer_from_NEW_into_OLD(_p);
												//temp_res +=	q_integr1D[i] * q_integr1D[j] * q_integr1D[T] *
												//	f(X)*det3(System_for_t._A);
												temp_res += q_integr1D[i] * q_integr1D[j] * q_integr1D[T] *
													Func_q(X, System_for_t, k)*math::GetDeterminantForMatrix3x3(System_for_t._A);
											}
										}

										double len_x = _up.x - _down.x;
										double len_y = _up.y - _down.y;
										double len_z = math::SolveLengthVector(F(end_t) - F(start_t));
										data_store_for_integrals_in_K[k][J][I] += temp_res * len_x * len_y * len_z / 8.0;
									}
								}
							}
						}

						for (int k = 0; k < 3; k++)
						{
							for (int j = 0; j < data_store_for_integrals_in_K[k].size(); j++)
							{
								for (int i = 0; i < data_store_for_integrals_in_K[k][j].size(); i++)
								{
									inJ[k] += data_store_for_integrals_in_K[k][j][i] /*/ SIGMA_PA_F*/;
								}
							}
						}
					}

					//solve the SIF: K_I, K_II
					double K_I = 0, K_II = 0;
					{
						//solving a point with progection into both side of crack faces
						Point<double> _P_up, P_up; // _P - in LOCAL coordinates, P - in GLOBAL
						Point<double> _P_down, P_down; // _P - in LOCAL coordinates, P - in GLOBAL
						newCoordSystem s_middle;
						NewSystem(0.5, s_middle); //t=0.5
						_P_up = s_middle._O - Point<double>(R_for_integration / 2., 0, 0.01); //The point is opposite to positive front orientation
						_P_down = s_middle._O - Point<double>(R_for_integration / 2., 0, -0.01); //The point is opposite to positive front orientation
						P_up = s_middle.transfer_from_NEW_into_OLD(_P_up);
						P_down = s_middle.transfer_from_NEW_into_OLD(_P_down);
						double _len = 0;
						int k_up = -1, k_down = -1;
						k_up = grid.GetNearestElementID(P_up, _len);
						k_down = grid.GetNearestElementID(P_down, _len);
						Point<double> U_up, U_down;
						U_up = grid.GetSolutionInPoint(k_up, P_up, U) + P_up;
						U_down = grid.GetSolutionInPoint(k_up, P_up, U) + P_down;


						//solving delta in LOCAL system
						Point<double> _U_up = s_middle.transfer_from_OLD_into_NEW(U_up);
						Point<double> _U_down = s_middle.transfer_from_OLD_into_NEW(U_down);
						double _delta_I = _U_up.y - _U_down.y;
						double _delta_II = _U_up.x - _U_down.x;

						//solving SIF (p.698: J.W. Eischen "An improved method for computing the J2 integral")
						std::complex<double> _compl_K_I, _compl_K_II, _compl_in;
						double _in = 1 - (inJ[1] / inJ[0])*(inJ[1] / inJ[0]);
						_compl_in = _in;
						_compl_K_I = sqrt(_compl_in);
						_compl_K_II = _compl_K_I;
						K_I = _compl_K_I.real();
						K_II = _compl_K_II.real();
						if (abs(_delta_I) < abs(_delta_II))
						{
							K_I = 1 - K_I;//signum??
							K_II = 1 + K_II; //signum??
						}
						else {
							K_I = 1 + K_I;//signum??
							K_II = 1 - K_II;//signum??
						}
						double E = grid.GetDomain(0)->forMech.GetE(0);
						K_I = math::GetSignum(_delta_I)*sqrt(E*inJ[1] / 2 * K_I);
						K_II = math::GetSignum(_delta_II)*sqrt(E*inJ[1] / 2 * K_I);
					}
					//Cherepanov's criterion (с.108: В.И. Гришин, В.Ю. Донченко "Метод расчета коэффициентов интенсивности напряжений в элементах конструкций с криволинейными трещинами")
					double Qrad = atan(-2 * K_I*K_II / (K_I*K_I + K_II * K_II)); //in rad

					double Qrad_simple = atan(inJ[1] / inJ[0]); //in rad
					printf("FR[%d], Point<double>[%d], J1[%.2e], J2[%.2e], J3[%.2e], Qgrad_new[%.2lf], Qgrad_old[%.2lf]\n", Fr, segm - 1, inJ[0], inJ[1], inJ[2], 180 * Qrad / M_PI, 180 * Qrad_simple / M_PI);
					Qrad = Qrad_simple;

					//solving critical value
					double G = abs(inJ[0] * cos(Qrad) + inJ[1] * sin(Qrad));
					double K_Ic = (grid.GetDomain(0))->forMech.GetK_Ic(); // песчаник
					double Step_crack_front = R_for_integration;
					//if (G < K_Ic) Step_crack_front /= 10;
					printf("\tG[%.2e], K_Ic[%.2e], Step_crack_front[%.2lf]\n", G, K_Ic, Step_crack_front);

					double Qgrad = 180 * Qrad / M_PI; //in grad

					//t=0
					newCoordSystem s0;
					NewSystem(0, s0);
					Point<double> _R, R;
					_R.z = s0._O.z;
					_R.x = cos(Qrad) * R_for_integration;
					_R.y = sqrt(R_for_integration*R_for_integration - _R.x*_R.x) *math::GetSignum(Qrad);
					R = s0.transfer_from_NEW_into_OLD(_R);
					R -= s0.O; //it is a VECTOR therefore we have the minus O
					//printf("FR[%d], Point<double>[%d], J1[%.2e], J2[%.2e], J3[%.2e], Qgrad[%.2lf], R(%.2e, %.2e, %.2e)=%.2e\n", Fr, segm - 1, inJ[0], inJ[1], inJ[2], Qgrad, R.x, R.y, R.z, math::SolveLengthVector(R));
					//printf("FR[%d], Point<double>[%d], J1[%.2e], J2[%.2e], J3[%.2e], Qgrad[%.2lf], R=%.2e, _R=%.2e\n", Fr, segm - 1, inJ[0], inJ[1], inJ[2], Qgrad, math::SolveLengthVector(R), math::SolveLengthVector(_R));
					grid.cracks[Cr].propogation[Fr][segm - 1].vect_of_propogation.push_back(R);
					grid.cracks[Cr].propogation[Fr][segm - 1].length_of_propogation.push_back(Step_crack_front);
					grid.cracks[Cr].propogation[Fr][segm - 1].id_old_point = grid.cracks[Cr].id_points_in_fronts[Fr][segm];
					grid.cracks[Cr].propogation[Fr][segm - 1].base_coord_system.push_back(newCoordSystem(s0));
					grid.cracks[Cr].propogation[Fr][segm - 1].Qgrad.push_back(Qgrad);
					grid.cracks[Cr].propogation[Fr][segm - 1].G.push_back(G);
					grid.cracks[Cr].propogation[Fr][segm - 1].J1.push_back(inJ[0]);
					grid.cracks[Cr].propogation[Fr][segm - 1].J2.push_back(inJ[1]);

					//t=1
					newCoordSystem s1;
					NewSystem(1, s1);
					_R.z = s1._O.z;
					_R.x = cos(Qrad) * R_for_integration;
					_R.y = sqrt(R_for_integration*R_for_integration - _R.x*_R.x) *math::GetSignum(Qrad);
					R = s1.transfer_from_NEW_into_OLD(_R) - s1.O;
					//printf("FR[%d], Point<double>[%d], J1[%.2e], J2[%.2e], J3[%.2e], Qgrad[%.2lf], R(%.2e, %.2e, %.2e)=%.2e\n", Fr, segm, inJ[0], inJ[1], inJ[2], Qgrad, R.x, R.y, R.z, math::SolveLengthVector(R));
					grid.cracks[Cr].propogation[Fr][segm].vect_of_propogation.push_back(R);
					grid.cracks[Cr].propogation[Fr][segm].length_of_propogation.push_back(Step_crack_front);
					grid.cracks[Cr].propogation[Fr][segm].id_old_point = grid.cracks[Cr].id_points_in_fronts[Fr][segm + 1];
					grid.cracks[Cr].propogation[Fr][segm].base_coord_system.push_back(newCoordSystem(s1));
					grid.cracks[Cr].propogation[Fr][segm].Qgrad.push_back(Qgrad);
					grid.cracks[Cr].propogation[Fr][segm].G.push_back(G);
					grid.cracks[Cr].propogation[Fr][segm].J1.push_back(inJ[0]);
					grid.cracks[Cr].propogation[Fr][segm].J2.push_back(inJ[1]);
				}

				//create new cracks
				//new nodes
				int old_size_xyz = grid.cracks[Cr].xyz.size();
				if (false) {
					for (int p = 0; p < grid.cracks[Cr].propogation[Fr].size(); p++)
					{
						Point<double> middle_vect;
						for (int v = 0; v < grid.cracks[Cr].propogation[Fr][p].vect_of_propogation.size(); v++)
						{
							middle_vect += grid.cracks[Cr].propogation[Fr][p].vect_of_propogation[v] / grid.cracks[Cr].propogation[Fr][p].vect_of_propogation.size();
						}
						middle_vect = middle_vect / math::SolveLengthVector(middle_vect) * grid.cracks[Cr].R;
						//printf("Len(FR[%d], Point[%d]) = %.2e\n", Fr, p, math::SolveLengthVector(middle_vect));
						grid.cracks[Cr].xyz.push_back(grid.cracks[Cr].xyz[grid.cracks[Cr].propogation[Fr][p].id_old_point] + middle_vect);
						grid.cracks[Cr].propogation[Fr][p].id_new_point = grid.cracks[Cr].xyz.size() - 1;
					}
				}

				//FULL MIDDLE MOVE!!!
				if (false)
				{
					Point<double> FULL_MIDDLE_VECTOR;
					for (int p = 0; p < false && grid.cracks[Cr].propogation[Fr].size(); p++)
					{
						Point<double> middle_vect;
						for (int v = 0; v < grid.cracks[Cr].propogation[Fr][p].vect_of_propogation.size(); v++)
						{
							middle_vect += grid.cracks[Cr].propogation[Fr][p].vect_of_propogation[v] / grid.cracks[Cr].propogation[Fr][p].vect_of_propogation.size();
						}
						middle_vect = middle_vect / math::SolveLengthVector(middle_vect) * grid.cracks[Cr].R;
						FULL_MIDDLE_VECTOR += middle_vect;
					}
					FULL_MIDDLE_VECTOR /= grid.cracks[Cr].propogation[Fr].size();
					for (int p = 0; p < grid.cracks[Cr].propogation[Fr].size(); p++)
					{
						grid.cracks[Cr].xyz.push_back(grid.cracks[Cr].xyz[grid.cracks[Cr].propogation[Fr][p].id_old_point] + FULL_MIDDLE_VECTOR);
						grid.cracks[Cr].propogation[Fr][p].id_new_point = grid.cracks[Cr].xyz.size() - 1;
					}
				}

				//smoothing via Ns nodes
				int Ns = grid.cracks[Cr].smoothing_coefficient;
				int end_of_update = grid.cracks[Cr].propogation[Fr].size();
				if (grid.cracks[Cr].propogation[Fr][0].id_old_point == grid.cracks[Cr].propogation[Fr][grid.cracks[Cr].propogation[Fr].size() - 1].id_old_point)
				{
					end_of_update--;
				}
				for (int id_base_vertex = 0; id_base_vertex < end_of_update; id_base_vertex++)
				{
					int start_id = id_base_vertex;
					if (id_base_vertex + Ns >= grid.cracks[Cr].propogation[Fr].size())
					{
						start_id = grid.cracks[Cr].propogation[Fr].size() - Ns;
						if (start_id < 0) start_id = 0;
					}

					double middle_Qgrad = 0;
					double middle_G = 0;
					double middle_J1 = 0;
					double middle_J2 = 0;
					int n_i = 0;

					for (int i = 0; i < Ns && (start_id + i) < grid.cracks[Cr].propogation[Fr].size(); i++)
					{
						for (int q = 0; q < grid.cracks[Cr].propogation[Fr][start_id + i].Qgrad.size(); q++)
						{
							middle_Qgrad += grid.cracks[Cr].propogation[Fr][start_id + i].Qgrad[q];
							middle_G += grid.cracks[Cr].propogation[Fr][start_id + i].G[q];
							middle_J1 += grid.cracks[Cr].propogation[Fr][start_id + i].J1[q];
							middle_J2 += grid.cracks[Cr].propogation[Fr][start_id + i].J2[q];
							n_i++;
						}
					}
					middle_Qgrad /= n_i;
					middle_G /= n_i;
					middle_J1 /= n_i;
					middle_J2 /= n_i;

					Point<double> propogation_vector;
					Point<double> simple_vector;
					double simpleQ = 0;
					double simpleG = 0;
					double length = 0;
					for (int q = 0; q < grid.cracks[Cr].propogation[Fr][id_base_vertex].base_coord_system.size(); q++)
					{
						propogation_vector += grid.cracks[Cr].propogation[Fr][id_base_vertex].base_coord_system[q].MakePropogation(middle_Qgrad, grid.cracks[Cr].R) /
							grid.cracks[Cr].propogation[Fr][id_base_vertex].base_coord_system.size();
						simple_vector += grid.cracks[Cr].propogation[Fr][id_base_vertex].vect_of_propogation[q] /
							grid.cracks[Cr].propogation[Fr][id_base_vertex].base_coord_system.size();
						simpleQ += grid.cracks[Cr].propogation[Fr][id_base_vertex].Qgrad[q] /
							grid.cracks[Cr].propogation[Fr][id_base_vertex].base_coord_system.size();
						simpleG += grid.cracks[Cr].propogation[Fr][id_base_vertex].G[q] /
							grid.cracks[Cr].propogation[Fr][id_base_vertex].base_coord_system.size();
						length += grid.cracks[Cr].propogation[Fr][id_base_vertex].length_of_propogation[q] /
							grid.cracks[Cr].propogation[Fr][id_base_vertex].base_coord_system.size();
					}
					if (middle_G < (grid.GetDomain(0))->forMech.GetK_Ic())
						length /= 10;

					propogation_vector = propogation_vector / math::SolveLengthVector(propogation_vector) * length /*grid.cracks[Cr].R*/;
					simple_vector = simple_vector / math::SolveLengthVector(simple_vector) * grid.cracks[Cr].R;

					/*printf_s("Front[%d], old_id[%d], Q=%.2lf; OldPoint(%.2lf, %.2lf, %.2lf); Vector=(%.2lf, %.2lf, %.2lf); SimpleQ=%.2lf; SimpleVector=(%.2lf, %.2lf, %.2lf)\n",
						Fr,
						grid.cracks[Cr].propogation[Fr][id_base_vertex].id_old_point,
						middle_Qgrad,
						grid.cracks[Cr].xyz[grid.cracks[Cr].propogation[Fr][id_base_vertex].id_old_point].x,
						grid.cracks[Cr].xyz[grid.cracks[Cr].propogation[Fr][id_base_vertex].id_old_point].y,
						grid.cracks[Cr].xyz[grid.cracks[Cr].propogation[Fr][id_base_vertex].id_old_point].z,
						propogation_vector.x,
						propogation_vector.y,
						propogation_vector.z,
						simpleQ,
						simple_vector.x,
						simple_vector.y,
						simple_vector.z);*/

					grid.cracks[Cr].xyz.push_back(grid.cracks[Cr].xyz[grid.cracks[Cr].propogation[Fr][id_base_vertex].id_old_point] + propogation_vector);
					grid.cracks[Cr].propogation[Fr][id_base_vertex].id_new_point = grid.cracks[Cr].xyz.size() - 1;
					grid.cracks[Cr].propogation[Fr][id_base_vertex]._vect_of_propogation = propogation_vector;
					grid.cracks[Cr].propogation[Fr][id_base_vertex]._length_of_propogation = length;
					grid.cracks[Cr].propogation[Fr][id_base_vertex]._Qgrad = middle_Qgrad;
					grid.cracks[Cr].propogation[Fr][id_base_vertex]._G = middle_G;
					grid.cracks[Cr].propogation[Fr][id_base_vertex]._J1 = middle_J1;
					grid.cracks[Cr].propogation[Fr][id_base_vertex]._J2 = middle_J2;
				}
				if (grid.cracks[Cr].propogation[Fr][0].id_old_point == grid.cracks[Cr].propogation[Fr][grid.cracks[Cr].propogation[Fr].size() - 1].id_old_point)
				{
					grid.cracks[Cr].propogation[Fr][grid.cracks[Cr].propogation[Fr].size() - 1].id_new_point = grid.cracks[Cr].propogation[Fr][0].id_new_point;
				}


				//new triangles and fronts
				Point<double> true_normal = grid.cracks[Cr].triangles[0].GetNormal();
				for (int segm = 0; segm < grid.cracks[Cr].fronts[Fr].size(); segm++)
				{
					Point<double> new_normal;
					int cur_segm[3] = { grid.cracks[Cr].fronts[Fr][segm].id_base_triangles, grid.cracks[Cr].fronts[Fr][segm].id_left, grid.cracks[Cr].fronts[Fr][segm].id_right };
					int new_segm[3];
					for (int ii = 1; ii < 3; ii++)
					{
						for (int p = 0; p < grid.cracks[Cr].propogation[Fr].size(); p++)
						{
							if (cur_segm[ii] == grid.cracks[Cr].propogation[Fr][p].id_old_point)
							{
								new_segm[ii] = grid.cracks[Cr].propogation[Fr][p].id_new_point;
								break;
							}
						}
					}

					for (int ii = 0; ii < grid.cracks[Cr].triangles.size(); ii++)
					{
						for (int jj = 0; jj < grid.cracks[Cr].triangles[ii].GetNodesCount(); jj++)
						{
							grid.cracks[Cr].triangles[ii].SetNode(jj, &grid.cracks[Cr].xyz[grid.cracks[Cr].triangles[ii].GetIdNode(jj)]);
						}
					}

					grid.cracks[Cr].triangles.push_back(geometry::Triangle());
					grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].SetIdNode(0, cur_segm[1]);
					grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].SetIdNode(1, cur_segm[2]);
					grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].SetIdNode(2, new_segm[2]);
					grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].SetNode(0, &grid.cracks[Cr].xyz[grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].GetIdNode(0)]);
					grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].SetNode(1, &grid.cracks[Cr].xyz[grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].GetIdNode(1)]);
					grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].SetNode(2, &grid.cracks[Cr].xyz[grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].GetIdNode(2)]);
					new_normal = grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].GetNormal();
					if (true_normal*new_normal < 0)
					{
						grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].SetIdNode(0, cur_segm[2]);
						grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].SetIdNode(1, cur_segm[1]);
						grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].SetIdNode(2, new_segm[2]);
						grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].SetNode(0, &grid.cracks[Cr].xyz[grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].GetIdNode(0)]);
						grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].SetNode(1, &grid.cracks[Cr].xyz[grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].GetIdNode(1)]);
						grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].SetNode(2, &grid.cracks[Cr].xyz[grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].GetIdNode(2)]);
					}

					grid.cracks[Cr].triangles.push_back(geometry::Triangle());
					grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].SetIdNode(0, cur_segm[1]);
					grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].SetIdNode(1, new_segm[1]);
					grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].SetIdNode(2, new_segm[2]);
					grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].SetNode(0, &grid.cracks[Cr].xyz[grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].GetIdNode(0)]);
					grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].SetNode(1, &grid.cracks[Cr].xyz[grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].GetIdNode(1)]);
					grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].SetNode(2, &grid.cracks[Cr].xyz[grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].GetIdNode(2)]);
					new_normal = grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].GetNormal();
					if (true_normal*new_normal < 0)
					{
						grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].SetIdNode(0, cur_segm[1]);
						grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].SetIdNode(1, new_segm[2]);
						grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].SetIdNode(2, new_segm[1]);
						grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].SetNode(0, &grid.cracks[Cr].xyz[grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].GetIdNode(0)]);
						grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].SetNode(1, &grid.cracks[Cr].xyz[grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].GetIdNode(1)]);
						grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].SetNode(2, &grid.cracks[Cr].xyz[grid.cracks[Cr].triangles[grid.cracks[Cr].triangles.size() - 1].GetIdNode(2)]);
					}

					new_segm[0] = grid.cracks[Cr].triangles.size() - 1;
					grid.cracks[Cr].fronts[Fr][segm].id_base_triangles = new_segm[0];
					grid.cracks[Cr].fronts[Fr][segm].id_left = new_segm[1];
					grid.cracks[Cr].fronts[Fr][segm].id_right = new_segm[2];
				}
			}
		}

		FILE* fout;
		char new_crack[1000];
		for (int Cr = 0; Cr < grid.cracks.size(); Cr++)
		{
			//sprintf_s(new_crack, "%s/Crack_%d_step%d.txt", dir, Cr, curr_STEP + 1);

			sprintf_s(new_crack, "%s/Crack_%d_properties.txt", dir, Cr, curr_STEP + 1);
			fopen_s(&fout, new_crack, "w");
			for (int Fr = 0; Fr < grid.cracks[Cr].propogation.size(); Fr++)
			{
				fprintf_s(fout, "Front %d\n", Fr);
				for (int i = 0; i < grid.cracks[Cr].propogation[Fr].size(); i++)
				{
					fprintf_s(fout, "point[%d]: id_old_point[%d], id_new_point[%d], vect(%.2lf, %.2lf, %.2lf), len[%.2lf], Qgrad[%.2lf], G[%.2e], J1[%.2e], J2[%.2e]\n",
						i,
						grid.cracks[Cr].propogation[Fr][i].id_old_point,
						grid.cracks[Cr].propogation[Fr][i].id_new_point,
						grid.cracks[Cr].propogation[Fr][i]._vect_of_propogation.x, grid.cracks[Cr].propogation[Fr][i]._vect_of_propogation.y, grid.cracks[Cr].propogation[Fr][i]._vect_of_propogation.z,
						grid.cracks[Cr].propogation[Fr][i]._length_of_propogation,
						grid.cracks[Cr].propogation[Fr][i]._Qgrad,
						grid.cracks[Cr].propogation[Fr][i]._G,
						grid.cracks[Cr].propogation[Fr][i]._J1,
						grid.cracks[Cr].propogation[Fr][i]._J2);
				}
				fprintf_s(fout, "\n");
			}
			fclose(fout);

			sprintf_s(new_crack, "%s/Crack_%d.txt", dir, Cr, curr_STEP + 1);
			fopen_s(&fout, new_crack, "w");
			fprintf_s(fout, "%d\n", grid.cracks[Cr].xyz.size());
			for (int i = 0; i < grid.cracks[Cr].xyz.size(); i++)
			{
				fprintf_s(fout, "%.6e\t%.6e\t%.6e\n",
					grid.cracks[Cr].xyz[i].x,
					grid.cracks[Cr].xyz[i].y,
					grid.cracks[Cr].xyz[i].z);
			}

			fprintf_s(fout, "\n%d\n", grid.cracks[Cr].triangles.size());
			for (int i = 0; i < grid.cracks[Cr].triangles.size(); i++)
			{
				fprintf_s(fout, "%d\t%d\t%d\n",
					grid.cracks[Cr].triangles[i].GetIdNode(0),
					grid.cracks[Cr].triangles[i].GetIdNode(1),
					grid.cracks[Cr].triangles[i].GetIdNode(2));
			}

			fprintf_s(fout, "\n%d\n", grid.cracks[Cr].fronts.size());
			for (int i = 0; i < grid.cracks[Cr].fronts.size(); i++)
			{
				fprintf_s(fout, "%d\n", grid.cracks[Cr].fronts[i].size());
				for (int j = 0; j < grid.cracks[Cr].fronts[i].size(); j++)
				{
					fprintf_s(fout, "%d\t%d\t%d\n",
						grid.cracks[Cr].fronts[i][j].id_base_triangles,
						grid.cracks[Cr].fronts[i][j].id_left,
						grid.cracks[Cr].fronts[i][j].id_right);
				}
			}
			fclose(fout);

			{ //tech_plot - as middle plane
				sprintf_s(new_crack, "%s/Crack_%d_step%d.dat", dir, Cr, curr_STEP + 1);
				fopen_s(&fout, new_crack, "w");

				fprintf_s(fout, "TITLE     = \"numerical\"\n");
				fprintf_s(fout, "VARIABLES = \"x\"\n \"y\"\n \"z\"\n");
				fprintf_s(fout, "\"crack\"\n");
				fprintf_s(fout, "ZONE T=\"Crack_%d_step%d\"\n", Cr, curr_STEP + 1);
				fprintf_s(fout, " N=%d,  E=%d, F=FEBLOCK ET=Triangle \n", grid.cracks[Cr].xyz.size(), grid.cracks[Cr].triangles.size());
				fprintf_s(fout, " VARLOCATION=(NODAL NODAL NODAL CELLCENTERED)\n");

				for (int i = 0; i < grid.cracks[Cr].xyz.size(); i++)
					fprintf_s(fout, "%.10e\n", grid.cracks[Cr].xyz[i].x);
				fprintf_s(fout, "\n");
				for (int i = 0; i < grid.cracks[Cr].xyz.size(); i++)
					fprintf_s(fout, "%.10e\n", grid.cracks[Cr].xyz[i].y);
				fprintf_s(fout, "\n");
				for (int i = 0; i < grid.cracks[Cr].xyz.size(); i++)
					fprintf_s(fout, "%.10e\n", grid.cracks[Cr].xyz[i].z);
				fprintf_s(fout, "\n");

				for (int i = 0; i < grid.cracks[Cr].triangles.size(); i++)
					fprintf_s(fout, "%d\n", curr_STEP + 1);
				fprintf_s(fout, "\n");

				for (int i = 0; i < grid.cracks[Cr].triangles.size(); i++)
				{
					for (int j = 0; j < 3; j++)
						fprintf_s(fout, "%d ", grid.cracks[Cr].triangles[i].GetIdNode(j) + 1);
					fprintf_s(fout, "\n");
				}

				fclose(fout);
			}

			{ //tech_plot - as face planes
				std::vector<Point<double>> xyz_up(grid.cracks[Cr].xyz.size()), xyz_down(grid.cracks[Cr].xyz.size());
				std::vector<int> elem_into_nodes(grid.cracks[Cr].xyz.size());
				double MaxSize = 0, MinSize = 1e+15;

				for (int id = 0; id < grid.cracks[Cr].triangles.size(); id++)
				{
					Point<double> normal = grid.cracks[Cr].triangles[id].GetNormal();
					normal /= math::SolveLengthVector(normal);
					Point<double> O = grid.cracks[Cr].triangles[id].GetWeightCentr();
					double move_koef = 0.5;

					Point<double> up_X = O + normal * move_koef;
					double len;
					int up_elem = grid.GetNearestElementID(up_X, len);
					Point<double> move_vector_UP = grid.GetSolutionInPoint(up_elem, up_X, U);

					Point<double> down_X = O - normal * move_koef;
					int down_elem = grid.GetNearestElementID(down_X, len);
					Point<double> move_vector_DOWN = grid.GetSolutionInPoint(down_elem, down_X, U);

					double move_size = math::SolveLengthVector(move_vector_DOWN) + math::SolveLengthVector(move_vector_UP);
					if (move_size > MaxSize) MaxSize = move_size;
					if (move_size < MinSize) MinSize = move_size;

					for (int _i = 0; _i < 3; _i++)
					{
						elem_into_nodes[grid.cracks[Cr].triangles[id].GetIdNode(_i)] ++;
						xyz_up[grid.cracks[Cr].triangles[id].GetIdNode(_i)] += move_vector_UP;
						xyz_down[grid.cracks[Cr].triangles[id].GetIdNode(_i)] += move_vector_DOWN;
					}
				}
				for (int i = 0; i < grid.cracks[Cr].fronts.size(); i++)
				{
					for (int j = 0; j < grid.cracks[Cr].fronts[i].size(); j++)
					{
						xyz_up[grid.cracks[Cr].fronts[i][j].id_left] = Point<double>(0, 0, 0);
						xyz_up[grid.cracks[Cr].fronts[i][j].id_right] = Point<double>(0, 0, 0);

						xyz_down[grid.cracks[Cr].fronts[i][j].id_left] = Point<double>(0, 0, 0);
						xyz_down[grid.cracks[Cr].fronts[i][j].id_right] = Point<double>(0, 0, 0);
					}
				}


				for (int id = 0; id < xyz_up.size(); id++)
				{
					xyz_up[id] /= elem_into_nodes[id];
					xyz_down[id] /= elem_into_nodes[id];
					xyz_up[id] += grid.cracks[Cr].xyz[id];
					xyz_down[id] += grid.cracks[Cr].xyz[id];
				}



				sprintf_s(new_crack, "%s/Crack_%d_step%d_in3D_MaxSize_%.2e_MinSize_%.2e.dat", dir, Cr, curr_STEP + 1, MaxSize, MinSize);
				fopen_s(&fout, new_crack, "w");

				fprintf_s(fout, "TITLE     = \"numerical\"\n");
				fprintf_s(fout, "VARIABLES = \"x\"\n \"y\"\n \"z\"\n");
				fprintf_s(fout, "\"crack\"\n");
				fprintf_s(fout, "ZONE T=\"Crack_%d_step%d\"\n", Cr, curr_STEP + 1);
				fprintf_s(fout, " N=%d,  E=%d, F=FEBLOCK ET=Triangle \n", xyz_up.size() * 2, grid.cracks[Cr].triangles.size() * 2);
				fprintf_s(fout, " VARLOCATION=(NODAL NODAL NODAL CELLCENTERED)\n");

				for (int i = 0; i < xyz_up.size(); i++)
					fprintf_s(fout, "%.10e\n", xyz_up[i].x);
				for (int i = 0; i < xyz_down.size(); i++)
					fprintf_s(fout, "%.10e\n", xyz_down[i].x);
				fprintf_s(fout, "\n");
				for (int i = 0; i < xyz_up.size(); i++)
					fprintf_s(fout, "%.10e\n", xyz_up[i].y);
				for (int i = 0; i < xyz_down.size(); i++)
					fprintf_s(fout, "%.10e\n", xyz_down[i].y);
				fprintf_s(fout, "\n");
				for (int i = 0; i < xyz_up.size(); i++)
					fprintf_s(fout, "%.10e\n", xyz_up[i].z);
				for (int i = 0; i < xyz_down.size(); i++)
					fprintf_s(fout, "%.10e\n", xyz_down[i].z);
				fprintf_s(fout, "\n");


				for (int i = 0; i < grid.cracks[Cr].triangles.size() * 2; i++)
					fprintf_s(fout, "%d\n", curr_STEP + 1);
				fprintf_s(fout, "\n");

				for (int i = 0; i < grid.cracks[Cr].triangles.size(); i++)
				{
					for (int j = 0; j < 3; j++)
						fprintf_s(fout, "%d ", grid.cracks[Cr].triangles[i].GetIdNode(j) + 1);
					fprintf_s(fout, "\n");
				}
				for (int i = 0; i < grid.cracks[Cr].triangles.size(); i++)
				{
					for (int j = 0; j < 3; j++)
						fprintf_s(fout, "%d ", grid.cracks[Cr].xyz.size() + grid.cracks[Cr].triangles[i].GetIdNode(j) + 1);
					fprintf_s(fout, "\n");
				}

				fclose(fout);
			}
		}
	}

	template <typename Dirichlet, typename Neumann>
	void MultiXFEM_forElasticDeformation(
		bool is_print_logFile, //input
		double MIN_RESIDUAL,
		char* cracks_directory, //input
		double step_size_for_crack,
		int crack_smoothing_coefficient,
		math::SimpleGrid &geo_grid, //input
		std::vector<Dirichlet> &first_boundary,//input
		std::vector<Neumann> &second_boundary,//input
		char* result_directory, //output
		MultiXFEM::Grid &solver_grid, //output
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

		clock_t t_after = clock();
		double start = omp_get_wtime();

		bool flag = true;
		for (int i = 0; flag; i++)
		{
			char name_crack[1000];
			sprintf_s(name_crack, sizeof(name_crack), "%s/Crack_%d.txt", cracks_directory, i);
			FILE *fin;
			fopen_s(&fin, name_crack, "r");
			if (fin != NULL)
			{
				fclose(fin);
				solver_grid.cracks.push_back(geometry::Crack(name_crack));
				solver_grid.cracks[solver_grid.cracks.size() - 1].R = step_size_for_crack;
				solver_grid.cracks[solver_grid.cracks.size() - 1].smoothing_coefficient = crack_smoothing_coefficient;
			}
			else {
				flag = false;
				break;
			}
		}

		printf("Initialization of grid...\n");
		solver_grid.Initialization(geo_grid, first_boundary, second_boundary);
		printf_s("complite                   \n\n");

		CSSD_Matrix<Tensor2Rank3D, Point<double>> global_SLAE;

		printf("Creation the SLAE portrait...");
		solver_grid.CreationPortrait(global_SLAE);
		printf_s("complite\n\n");

		//{
		//	char name[1000];
		//	sprintf_s(name, sizeof(name), "%s/Matrix_portrait.txt", result_directory);
		//	PrintPortrait(name, global_SLAE);
		//}

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
				std::function< double(Point<double> X) > koefCohesive = [&](Point<double> X)->double
				{
					return (solver_grid.GetDomain(0))->forMech.GetPressureOfFluid();
				};
				element->SolveLocalMatrix(local_SLAE[id_elem], koefD);
				element->SolveLocalCohesiveMatrix(local_SLAE[id_elem], koefCohesive);

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
					if (true)
					{
						double target_coef_in_firgt_side = 0;
						switch (position)
						{
						case 0:
							target_coef_in_firgt_side = global_SLAE.F[global_id].x;
							global_SLAE.F[global_id].y -= target_coef_in_firgt_side *global_SLAE.Diag[global_id].val[1][position];
							global_SLAE.F[global_id].z -= target_coef_in_firgt_side *global_SLAE.Diag[global_id].val[2][position];
							global_SLAE.Diag[global_id].val[0][position] = 0.0;
							global_SLAE.Diag[global_id].val[1][position] = 0.0;
							break;
						case 1:
							target_coef_in_firgt_side = global_SLAE.F[global_id].y;
							global_SLAE.F[global_id].x -= target_coef_in_firgt_side *global_SLAE.Diag[global_id].val[0][position];
							global_SLAE.F[global_id].z -= target_coef_in_firgt_side *global_SLAE.Diag[global_id].val[2][position];
							global_SLAE.Diag[global_id].val[0][position] = 0.0;
							global_SLAE.Diag[global_id].val[2][position] = 0.0;
							break;
						case 2:
							target_coef_in_firgt_side = global_SLAE.F[global_id].z;
							global_SLAE.F[global_id].x -= target_coef_in_firgt_side *global_SLAE.Diag[global_id].val[0][position];
							global_SLAE.F[global_id].y -= target_coef_in_firgt_side *global_SLAE.Diag[global_id].val[1][position];
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
									global_SLAE.F[i].x -= target_coef_in_firgt_side *global_SLAE.A_down[i][j].val[0][position];
									global_SLAE.F[i].y -= target_coef_in_firgt_side *global_SLAE.A_down[i][j].val[1][position];
									global_SLAE.F[i].z -= target_coef_in_firgt_side *global_SLAE.A_down[i][j].val[2][position];
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
		if(true){ //via block matrix
			printf("Soluting SLAY... (%d)\n", global_SLAE.GetMatrixSize());
			int MaxSize = global_SLAE.GetMatrixSize()*3;
			std::vector<Point<double>> best_solution;
			math::MakeCopyVector_A_into_B(global_SLAE.X, best_solution);
			double current_residual = 1, best_residual = 1e+25;
			int MAX_STEPS = 14;
			for (int i = 0; i <= MAX_STEPS; i++)
			{
				current_residual = pow(10., -1 * (i + 1));
				//current_residual /= 2.;
				printf_s("//---> I = %d/%d (full size %d) - needed residual %.2e\n", i, MAX_STEPS, global_SLAE.GetMatrixSize(), current_residual);
				current_residual = abs(global_SLAE.BiCG_Stab(MaxSize, current_residual));
				//current_residual = abs(global_SLAE.MSG(MaxSize, current_residual));
				if (current_residual <= MIN_RESIDUAL)
					break;
				if (current_residual > best_residual + (1e-10))
				{
					//math::MakeCopyVector_A_into_B(best_solution, global_SLAE.X);
					//printf_s("//---> BEST residual %.2e\n", best_residual);
					//break;
				}
				if (current_residual < best_residual)
				{
					best_residual = current_residual;
					math::MakeCopyVector_A_into_B(global_SLAE.X, best_solution);
					printf_s("//---> BEST residual %.2e\n", best_residual);
				}
			}
			printf_s("//---> BEST residual %.2e\n", best_residual);
			math::MakeCopyVector_A_into_B(best_solution, global_SLAE.X);
			math::MakeCopyVector_A_into_B(global_SLAE.X, Solution);
		}
		if(false){ //classic matrix
			printf("Soluting SLAY... (%d)\n", global_SLAE.GetMatrixSize());
			int MaxSize = global_SLAE.GetMatrixSize();
			

			CSSD_Matrix<double, double> tmp_matrix;
			math::MakeCopyMatrix_A_into_B(global_SLAE, tmp_matrix);
			std::vector<double> best_solution;
			math::MakeCopyVector_A_into_B(tmp_matrix.X, best_solution);

			double current_residual = 1, best_residual = 1e+25;
			int MAX_STEPS = 12;
			for (int i = 0; i <= MAX_STEPS; i++)
			{
				current_residual = pow(10., -1 * (i + 1));
				//current_residual /= 2.;
				printf_s("//---> I = %d/%d (full size %d) - needed residual %.2e\n", i, MAX_STEPS, tmp_matrix.GetMatrixSize(), current_residual);
				//current_residual = abs(tmp_matrix.BiCG_Stab(MaxSize, current_residual));
				current_residual = abs(tmp_matrix.MSG(MaxSize, current_residual));
				if (current_residual <= MIN_RESIDUAL)
					break;
				if (current_residual > best_residual + (1e-10))
				{
					math::MakeCopyVector_A_into_B(best_solution, tmp_matrix.X);
					printf_s("//---> BEST residual %.2e\n", best_residual);
					break;
				}
				if (current_residual < best_residual)
				{
					best_residual = current_residual;
					math::MakeCopyVector_A_into_B(tmp_matrix.X, best_solution);
				}
			}
			printf_s("//---> BEST residual %.2e\n", best_residual);
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
}