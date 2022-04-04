#pragma once
#include "stdio.h"
#include <vector>
#include "Math.h"

template <typename TForMatrix, typename TForVector> class DenseMatrix
{
private:
	std::vector<std::vector<TForMatrix>> LU;

	void LUdecomposition()
	{
		LU.resize( A.size() );
		for( int i = 0; i < A.size(); i++ )
			LU[i].resize( A.size() );
		TForMatrix summ;

		for( int i = 0; i < LU.size(); i++ )
			for( int j = 0; j < LU.size(); j++ )
			{
				if( i >= j)
				{
					summ = 0;
					for( int k = 0; k < j; k++ )
						summ += LU[i][k] * LU[k][j];
					LU[i][j] = A[i][j] - summ;
				} else {
					summ = 0;
					for( int k = 0; k < i; k++ )
						summ += LU[i][k] * LU[k][j];
					LU[i][j] = (A[i][j] - summ) / LU[i][i];
				}
			}
	};

	void zero_in_diag()
	{
		for( int i = 0; i < this->A.size(); i++ )
		{
			if( math::IsEqual(A[i][i],0.) && this->A.size() != 1)
			{
				if( i == this->A.size()-1 )
				{
					for( int k = this->A.size()-1; k > 0; k-- )
						if( !math::IsEqual(A[k][i],0.) )
						{
							for( int j = 0; j < this->A.size(); j++ )
								A[i][j] += A[k][j];
							F[i] += F[k];
							break;
						}
				} else {
					for( int k = i+1; k < this->A.size(); k++ )
						if( !math::IsEqual(A[k][i],0.) )
						{
							for( int j = 0; j < this->A.size(); j++ )
								A[i][j] += A[k][j];
							F[i] += F[k];
							break;
						}
				}
			}
		}
	};

public:

	std::vector<std::vector<TForMatrix>> A;
	std::vector<TForVector> F;
	std::vector<TForVector> X;
	std::vector<TForVector> Y;

	//N - размерность матрицы
	DenseMatrix( int N )
	{
		A.resize(N);
		for( int i = 0; i < N; i++ )
			A[i].resize(N);

		F.resize(N);

		X.resize(N);
	};
	//N-строки, M-столбцы
	DenseMatrix(int N, int M)
	{
		A.resize(N);
		for (int i = 0; i < N; i++)
			A[i].resize(M);

		F.resize(M);

		X.resize(M);
	};
	DenseMatrix(){};

	int GetSize(){return (int)X.size();};
	void SetSize(int N)
	{
		A.resize(N);
		for (int i = 0; i < N; i++)
			A[i].resize(N);

		F.resize(N);

		X.resize(N);
	}

	void SetElementOfA ( int i, int j, TForMatrix value )
	{
		A[i][j] = value;
	};
	TForMatrix GetElementOfA ( int i, int j )
	{
		return A[i][j];
	};
	std::vector<std::vector<TForMatrix>> GetA ()
	{
		return A;
	};

	void SetF ( int i, TForVector value )
	{
		F[i] = value;
	};
	void SetX ( int i, TForVector value)
	{
		X[i] = value;
	};
	TForVector GetElementOfF ( int i )
	{
		return F[i];
	};
	std::vector<TForVector> GetF()
	{
		return F;
	};

	TForVector GetElementOfX ( int i )
	{
		return X[i];
	};
	std::vector<TForVector> GetX ()
	{
		return X;
	};

	std::vector<TForVector> Gauss_LU()
	{
		zero_in_diag();

		for( int j = 0; j < X.size(); j++ )
		{
			X[j] = 0;
		}

		LUdecomposition();
		//
		for( int i = 0; i < LU.size(); i++ )
		{
			TForVector summ = 0;
			for( int j = 0; j < i; j++ )
				summ += LU[i][j] * X[j];
			X[i] = (F[i] - summ) / LU[i][i]; 
		}
		//
		for( int j = LU.size()-1; j >= 0; j-- )
		{
			for( int i = 0; i < j; i++ )
				X[i] -= LU[i][j] * X[j];
		}
		//

		Y.resize(F.size());
		for( int i = 0; i < F.size(); i++ )
		{
			Y[i] = 0;
			for( int j = 0; j < F.size(); j++ )
			{
				Y[i]+=A[i][j]*X[j];
			}
			if(!math::IsEqual(Y[i],F[i]))
			{
				printf("error in SLAE (Gauss_LU)!!!(%.5e, %.5e) => %.5e\n", Y[i], F[i], Y[i] - F[i]);
				//Sleep(100);
			}
		}

		return X;
	};

	std::vector<TForVector> Gauss()
	{
		for( int j = 0; j < X.size(); j++ )
		{
			X[j] = 0;
		}

		std::vector<std::vector<TForMatrix>> A_new;
		std::vector<TForVector> b_new(this->F.size());
		A_new.resize(A.size());
		std::vector<int> line(A.size());
		for( int i = 0; i < X.size(); i++ )
		{
			A_new[i].resize(A[i].size());
			for( int j = 0; j < X.size(); j++)
			{
				/*if( math::IsEqual( A[i][j], 0.0 ) )
					A[i][j] = 0.0;*/
				A_new[i][j] = A[i][j];
			}
			line[i] = i;
			/*if( math::IsEqual( this->F[i], 0.0 ) )
				this->F[i] = 0.0;*/
			b_new[i] = this->F[i];
		}

		//{
		//	for( int j = 0; j < X.size(); j++)
		//	{
		//		/*if( math::IsEqual( A[i][j], 0.0 ) )
		//			A[i][j] = 0.0;*/
		//		A_new[1][j] = A[2][j];
		//	}
		//	b_new[1] = this->F[2];
		//	for( int j = 0; j < X.size(); j++)
		//	{
		//		/*if( math::IsEqual( A[i][j], 0.0 ) )
		//			A[i][j] = 0.0;*/
		//		A_new[2][j] = A[1][j];
		//	}
		//	b_new[2] = this->F[1];
		//}
		int current_row = 0;
		TForVector zero = 0;
		for( int i = 0; i < line.size(); i++ )
		{
			//if (i % 100 == 0)
			//	printf(">>current i = %d\r", i);
			if( math::IsEqual( A_new[line[i]][current_row], zero ) )
			{
				for( int ii = i+1; ii < line.size(); ii++ )
				{
					if( !math::IsEqual( A_new[line[ii]][current_row], zero ) )
					{
						int _tmp = line[ii]; line[ii] = line[i]; line[i] = _tmp;
						/*printf("----------A_new------\n");
						for (int iii = 0; iii < F.size(); iii++)
						{
							for (int jjj = 0; jjj < F.size(); jjj++)
							{
								printf("%.5e ", A_new[line[iii]][jjj]);
							}
							printf("\n");
						}
						printf("----------x------\n");*/
						break;
					}
				}
			}
			if( !math::IsEqual( A_new[line[i]][current_row], zero ) )
			{
				TForMatrix koef = A_new[line[i]][current_row];
				for( int jj = current_row; jj < line.size(); jj++ )
					A_new[line[i]][jj] /= koef;
				b_new[line[i]] /= koef;

				for( int ii = i+1; ii < line.size(); ii++ )
				{
					TForMatrix mult = A_new[line[ii]][current_row];
					for( int jj = current_row; jj < line.size(); jj++ )
					{
						A_new[line[ii]][jj] -= mult * A_new[line[i]][jj];
					}
					b_new[line[ii]] -= mult*b_new[line[i]];
				}
			} else {i--;}

			current_row++;
			if( current_row >= line.size() )
				break;
		}

		for( int i = (int)A_new.size()-1; i >= 0; i-- )
		{
			int j = 0;
			for(; j < A_new.size(); j++)
			{
				if(!math::IsEqual(A_new[line[i]][j],zero))
				{
					break;
				}
			}
			if( j < A_new.size())
			{
				TForVector _tmp = 0.0;
				for( int jj = j+1; jj < A_new.size(); jj++ )
				{
					_tmp += X[jj]*A_new[line[i]][jj];
				}
				X[j] = b_new[line[i]] - _tmp;
			}
		}
		//
		/*for( int i = 0; i < LU.size(); i++ )
		{
			double summ = 0;
			for( int j = 0; j < i; j++ )
				summ += LU[i][j] * X[j];
			X[i] = (F[i] - summ) / LU[i][i]; 
		}*/
		//
		/*for( int j = LU.size()-1; j >= 0; j-- )
		{
			for( int i = 0; i < j; i++ )
				X[i] -= A_new[line[i]][j] * X[j];
		}*/
		//

		Y.resize(F.size());
		for( int i = 0; i < F.size(); i++ )
		{
			Y[i] = 0;
			for( int j = 0; j < F.size(); j++ )
			{
				Y[i]+=A[i][j]*X[j];
			}
			if(abs(Y[i] - F[i]) > 1E-5)
			{
				printf("error in SLAE (Gauss)!!!(%.5e, %.5e) => %.5e\n", Y[i], F[i], abs(Y[i] - F[i]));
				printf("----------A------\n");
				for( int ii = 0; ii < F.size(); ii++ )
				{
					for( int jj = 0; jj < F.size(); jj++ )
					{
						printf("%.5e ", A[ii][jj]);
					}
					printf("\n");
				}
				printf("----------x------\n");
				for( int ii = 0; ii < F.size(); ii++ )
				{
					printf("%.5e ", X[ii]);
				}
				printf("\n");
				printf("----------b------\n");
				for( int ii = 0; ii < F.size(); ii++ )
				{
					printf("%.5e ", F[ii]);
				}
				printf("\n");
				printf("----------A_new------\n");
				for( int ii = 0; ii < F.size(); ii++ )
				{
					for( int jj = 0; jj < F.size(); jj++ )
					{
						printf("%.5e ", A_new[ii][jj]);
					}
					printf("\n");
				}
				printf("----------b_new------\n");
				for( int ii = 0; ii < F.size(); ii++ )
				{
					printf("%.5e ", b_new[ii]);
				}
				printf("\n");

				//Sleep(100);
			}
		}

		return X;
	};
	std::vector<TForVector> Kramer()
	{
		if(this->A.size()!=3)
		{
			return X;
		}

		for (int j = 0; j < X.size(); j++)
		{
			X[j] = 0;
		}

		double D = this->det();
		std::vector<std::vector<TForMatrix>> A_new;
		A_new.resize(A.size());
		for (int i = 0; i < A_new.size(); i++)
		{
			A_new[i].resize(A.size());
			A_new[i][0] = this->F[i];
			A_new[i][1] = this->A[i][1];
			A_new[i][2] = this->A[i][2];
		}
		this->X[0] = det3(A_new) / D;

		for (int i = 0; i < A_new.size(); i++)
		{
			A_new[i][0] = this->A[i][0];
			A_new[i][1] = this->F[i];
			A_new[i][2] = this->A[i][2];
		}
		this->X[1] = det3(A_new) / D;

		for (int i = 0; i < A_new.size(); i++)
		{
			A_new[i][0] = this->A[i][0];
			A_new[i][1] = this->A[i][1];
			A_new[i][2] = this->F[i];
		}
		this->X[2] = det3(A_new) / D;

		Y.resize(F.size());
		for (int i = 0; i < F.size(); i++)
		{
			Y[i] = 0;
			for (int j = 0; j < F.size(); j++)
			{
				Y[i] += A[i][j] * X[j];
			}
			if (abs(Y[i] - F[i]) > 1E-5)
			{
				printf("error in SLAE (Kramer)!!!(%.5e, %.5e) => %.5e\n", Y[i], F[i], abs(Y[i] - F[i]));
				/*printf("----------A------\n");
				for( int ii = 0; ii < F.size(); ii++ )
				{
				for( int jj = 0; jj < F.size(); jj++ )
				{
				printf("%.5e ", A[ii][jj]);
				}
				printf("\n");
				}
				printf("----------x------\n");
				for( int ii = 0; ii < F.size(); ii++ )
				{
				printf("%.5e ", X[ii]);
				}
				printf("\n");
				printf("----------b------\n");
				for( int ii = 0; ii < F.size(); ii++ )
				{
				printf("%.5e ", F[ii]);
				}
				printf("\n");
				printf("----------A_new------\n");
				for( int ii = 0; ii < F.size(); ii++ )
				{
				for( int jj = 0; jj < F.size(); jj++ )
				{
				printf("%.5e ", A_new[ii][jj]);
				}
				printf("\n");
				}
				printf("----------b_new------\n");
				for( int ii = 0; ii < F.size(); ii++ )
				{
				printf("%.5e ", b_new[ii]);
				}
				printf("\n");*/

				//Sleep(100);
			}
		}

		return X;
	};

	DenseMatrix& operator = ( DenseMatrix& Matr )
	{
		this->set_size(Matr.size());
		for( int i = 0; i < this->A.size(); i++  )
			for( int j = 0; j < this->A.size(); j++  )
				this->A[i][j] = Matr.get_A( i,j );

		return *this;
	}
	DenseMatrix& operator + (DenseMatrix& right) 
	{
		for( int i = 0; i < this->A.size(); i++  )
			for( int j = 0; j < this->A.size(); j++  )
				this->A[i][j] += right.get_A( i,j );

		return *this;
	}

	DenseMatrix& operator / (const TForMatrix& A) 
	{
		for( int i = 0; i < this->A.size(); i++  )
			for( int j = 0; j < this->A.size(); j++  )
				this->A[i][j] /= A;

		return *this;
	}

	double det()
	{
		double deet;
		if(X.size()==3)
		{
			deet  = det3( this->A );
		}
		return deet;
	}
};

class DenseMatrix_Tensor
{
private:
	std::vector<std::vector<double>> A_simple;
	std::vector<double> F_simple, X_simple, Y_simple;
	std::vector<std::vector<double>> LU;

	void CreateSimpleMatrix()
	{
		math::ResizeVector(A_simple, A.size() * 3, A.size() * 3);
		F_simple.resize(A_simple.size());
		X_simple.resize(A_simple.size());

		for (int i = 0; i < A.size(); i++)
		{
			for (int j = 0; j < A.size(); j++)
			{
				for (int ii = 0; ii < 3; ii++)
				{
					for (int jj = 0; jj < 3; jj++)
					{
						A_simple[i * 3 + ii][j * 3 + jj] = A[i][j].val[ii][jj];
					}
				}
			}
			F_simple[i * 3 + 0] = F[i].x;
			F_simple[i * 3 + 1] = F[i].y;
			F_simple[i * 3 + 2] = F[i].z;

			X_simple[i * 3 + 0] = X[i].x;
			X_simple[i * 3 + 1] = X[i].y;
			X_simple[i * 3 + 2] = X[i].z;
		}
	}
	void LUdecomposition()
	{
		LU.resize(A_simple.size());
		for (int i = 0; i < A_simple.size(); i++)
			LU[i].resize(A_simple.size());
		double summ;

		for (int i = 0; i < LU.size(); i++)
			for (int j = 0; j < LU.size(); j++)
			{
				if (i >= j)
				{
					summ = 0;
					for (int k = 0; k < j; k++)
						summ += LU[i][k] * LU[k][j];
					LU[i][j] = A_simple[i][j] - summ;
				}
				else {
					summ = 0;
					for (int k = 0; k < i; k++)
						summ += LU[i][k] * LU[k][j];
					LU[i][j] = (A_simple[i][j] - summ) / LU[i][i];
				}
			}
	};
	void zero_in_diag()
	{
		for (int i = 0; i < this->A.size(); i++)
		{
			if (math::IsEqual(A_simple[i][i], 0.) && this->A_simple.size() != 1)
			{
				if (i == this->A_simple.size() - 1)
				{
					for (int k = this->A_simple.size() - 1; k > 0; k--)
						if (!math::IsEqual(A_simple[k][i], 0.))
						{
							for (int j = 0; j < this->A_simple.size(); j++)
								A_simple[i][j] += A_simple[k][j];
							F_simple[i] += F_simple[k];
							break;
						}
				}
				else {
					for (int k = i + 1; k < this->A_simple.size(); k++)
						if (!math::IsEqual(A_simple[k][i], 0.))
						{
							for (int j = 0; j < this->A_simple.size(); j++)
								A_simple[i][j] += A_simple[k][j];
							F_simple[i] += F_simple[k];
							break;
						}
				}
			}
		}
	};

public:

	std::vector<std::vector<Tensor2Rank3D>> A;
	std::vector<Point<double>> F;
	std::vector<Point<double>> X;
	
	//N - размерность матрицы
	DenseMatrix_Tensor(int N)
	{
		A.resize(N);
		for (int i = 0; i < N; i++)
			A[i].resize(N);

		F.resize(N);

		X.resize(N);
	};
	DenseMatrix_Tensor() {};

	int GetSize() { return (int)X.size(); };
	void SetSize(int N)
	{
		A.resize(N);
		for (int i = 0; i < N; i++)
			A[i].resize(N);

		F.resize(N);

		X.resize(N);
	}

	void SetElementOfA(int i, int j, Tensor2Rank3D &value)
	{
		A[i][j] = value;
	};
	Tensor2Rank3D GetElementOfA(int i, int j)
	{
		return A[i][j];
	};
	std::vector<std::vector<Tensor2Rank3D>> GetA()
	{
		return A;
	};

	void SetF(int i, Point<double> value)
	{
		F[i] = value;
	};
	void SetX(int i, Point<double> value)
	{
		X[i] = value;
	};
	Point<double> GetElementOfF(int i)
	{
		return F[i];
	};
	std::vector<Point<double>> GetF()
	{
		return F;
	};
	Point<double> GetElementOfX(int i)
	{
		return X[i];
	};
	std::vector<Point<double>> GetX()
	{
		return X;
	};

	void Gauss_LU()
	{
		CreateSimpleMatrix();

		zero_in_diag();

		for (int j = 0; j < X_simple.size(); j++)
		{
			X_simple[j] = 0;
		}

		LUdecomposition();
		//
		for (int i = 0; i < LU.size(); i++)
		{
			double summ = 0;
			for (int j = 0; j < i; j++)
				summ += LU[i][j] * X_simple[j];
			X_simple[i] = (F_simple[i] - summ) / LU[i][i];
		}
		//
		for (int j = LU.size() - 1; j >= 0; j--)
		{
			for (int i = 0; i < j; i++)
				X_simple[i] -= LU[i][j] * X_simple[j];
		}
		//

		for (int i = 0; i < X.size(); i++)
		{
			X[i].x = X_simple[i * 3 + 0];
			X[i].y = X_simple[i * 3 + 1];
			X[i].z = X_simple[i * 3 + 2];
		}

		Y_simple.resize(F_simple.size());
		for (int i = 0; i < F_simple.size(); i++)
		{
			Y_simple[i] = 0;
			for (int j = 0; j < F_simple.size(); j++)
			{
				Y_simple[i] += A_simple[i][j] * X_simple[j];
			}
			if (!math::IsEqual(Y_simple[i], F_simple[i]))
			{
				printf("error in SLA_simpleE (Gauss_LU)!!!(%.5e, %.5e) => %.5e\n", Y_simple[i], F_simple[i], Y_simple[i] - F_simple[i]);
				//Sleep(100);
			}
		}

	};
	void Gauss()
	{
		CreateSimpleMatrix();

		for (int j = 0; j < X.size(); j++)
		{
			X_simple[j] = 0;
		}

		std::vector<std::vector<double>> A_new;
		std::vector<double> b_new(this->F_simple.size());
		A_new.resize(A_simple.size());
		std::vector<int> line(A_simple.size());
		for (int i = 0; i < X_simple.size(); i++)
		{
			A_new[i].resize(A_simple[i].size());
			for (int j = 0; j < X_simple.size(); j++)
			{
				A_new[i][j] = A_simple[i][j];
			}
			line[i] = i;
			b_new[i] = this->F_simple[i];
		}

		int current_row = 0;
		double zero = 0;
		for (int i = 0; i < line.size(); i++)
		{
			if (math::IsEqual(A_new[line[i]][current_row], zero))
			{
				for (int ii = i + 1; ii < line.size(); ii++)
				{
					if (!math::IsEqual(A_new[line[ii]][current_row], zero))
					{
						int _tmp = line[ii]; line[ii] = line[i]; line[i] = _tmp;
						/*printf("----------A_new------\n");
						for (int iii = 0; iii < F_simple.size(); iii++)
						{
							for (int jjj = 0; jjj < F_simple.size(); jjj++)
							{
								printf("%.5e ", A_new[line[iii]][jjj]);
							}
							printf("\n");
						}
						printf("----------x------\n");*/
						break;
					}
				}
			}
			if (!math::IsEqual(A_new[line[i]][current_row], zero))
			{
				double koef = A_new[line[i]][current_row];
				for (int jj = current_row; jj < line.size(); jj++)
					A_new[line[i]][jj] /= koef;
				b_new[line[i]] /= koef;

				for (int ii = i + 1; ii < line.size(); ii++)
				{
					double mult = A_new[line[ii]][current_row];
					for (int jj = current_row; jj < line.size(); jj++)
					{
						A_new[line[ii]][jj] -= mult * A_new[line[i]][jj];
					}
					b_new[line[ii]] -= mult * b_new[line[i]];
				}
			}
			else { i--; }

			current_row++;
			if (current_row >= line.size())
				break;
		}

		for (int i = (int)A_new.size() - 1; i >= 0; i--)
		{
			int j = 0;
			for (; j < A_new.size(); j++)
			{
				if (!math::IsEqual(A_new[line[i]][j], zero))
				{
					break;
				}
			}
			if (j < A_new.size())
			{
				double _tmp = 0.0;
				for (int jj = j + 1; jj < A_new.size(); jj++)
				{
					_tmp += X_simple[jj] * A_new[line[i]][jj];
				}
				X_simple[j] = b_new[line[i]] - _tmp;
			}
		}

		Y_simple.resize(F_simple.size());
		for (int i = 0; i < F_simple.size(); i++)
		{
			Y_simple[i] = 0;
			for (int j = 0; j < F_simple.size(); j++)
			{
				Y_simple[i] += A_simple[i][j] * X_simple[j];
			}
			if (abs(Y_simple[i] - F_simple[i]) > 1E-5)
			{
				printf("error in SLAE (Gauss)!!!(%.5e, %.5e) => %.5e\n", Y_simple[i], F_simple[i], abs(Y_simple[i] - F_simple[i]));
				printf("----------A------\n");
				for (int ii = 0; ii < F_simple.size(); ii++)
				{
					for (int jj = 0; jj < F_simple.size(); jj++)
					{
						printf("%.5e ", A_simple[ii][jj]);
					}
					printf("\n");
				}
				printf("----------x------\n");
				for (int ii = 0; ii < F_simple.size(); ii++)
				{
					printf("%.5e ", X_simple[ii]);
				}
				printf("\n");
				printf("----------b------\n");
				for (int ii = 0; ii < F_simple.size(); ii++)
				{
					printf("%.5e ", F_simple[ii]);
				}
				printf("\n");
				printf("----------A_new------\n");
				for (int ii = 0; ii < F_simple.size(); ii++)
				{
					for (int jj = 0; jj < F_simple.size(); jj++)
					{
						printf("%.5e ", A_new[ii][jj]);
					}
					printf("\n");
				}
				printf("----------b_new------\n");
				for (int ii = 0; ii < F_simple.size(); ii++)
				{
					printf("%.5e ", b_new[ii]);
				}
				printf("\n");

				//Sleep(100);
			}
		}

		for (int i = 0; i < X.size(); i++)
		{
			X[i].x = X_simple[i * 3 + 0];
			X[i].y = X_simple[i * 3 + 1];
			X[i].z = X_simple[i * 3 + 2];
		}
	};
};