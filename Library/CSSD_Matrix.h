
#pragma once
#include <vector>

//SLAE solution AX=F 
template <typename TForMatrix, typename TForVector> 
class CSSD_Matrix{
private:
public:
	std::vector<TForMatrix> Diag;
	std::vector<std::vector<TForMatrix>> A_up, A_down;
	std::vector<std::vector<int>> id_column_for_A_up, id_column_for_A_down;

	std::vector<TForVector> X, F;
	bool print_logs;
	FILE *log_out;

	CSSD_Matrix() {};
	~CSSD_Matrix() {
		std::vector<TForMatrix> v_Diag;
		std::vector<TForMatrix>(v_Diag).swap(this->Diag);
		std::vector<std::vector<TForMatrix>> v_A_up;
		std::vector<std::vector<TForMatrix>>(v_A_up).swap(this->A_up);
		std::vector<std::vector<TForMatrix>>(v_A_up).swap(this->A_down);
		std::vector<std::vector<int>> v_id_column_for_A_up;
		std::vector<std::vector<int>>(v_id_column_for_A_up).swap(this->id_column_for_A_up);
		std::vector<std::vector<int>>(v_id_column_for_A_up).swap(this->id_column_for_A_down);
		std::vector<TForVector> v_X;
		std::vector<TForVector>(v_X).swap(this->X);
		std::vector<TForVector>(v_X).swap(this->F);
	};

	int GetMatrixSize()
	{
		return (int)Diag.size();
	}
	void SetMatrix(CSSD_Matrix<TForMatrix, TForVector>& A)
	{
		math::MakeCopyVector_A_into_B(A.Diag, this->Diag);
		math::MakeCopyVector_A_into_B(A.X, this->X);
		math::MakeCopyVector_A_into_B(A.F, this->F);

		this->print_logs = A.print_logs;
		this->log_out = A.log_out;
		
		this->A_up.resize(A.A_up.size());
		for (int i = 0; i < A.A_up.size(); i++)
		{
			math::MakeCopyVector_A_into_B(A.A_up[i], this->A_up[i]);
		}
		this->A_down.resize(A.A_down.size());
		for (int i = 0; i < A.A_down.size(); i++)
		{
			math::MakeCopyVector_A_into_B(A.A_down[i], this->A_down[i]);
		}

		this->id_column_for_A_up.resize(A.id_column_for_A_up.size());
		for (int i = 0; i < A.id_column_for_A_up.size(); i++)
		{
			math::MakeCopyVector_A_into_B(A.id_column_for_A_up[i], this->id_column_for_A_up[i]);
		}
		this->id_column_for_A_down.resize(A.id_column_for_A_down.size());
		for (int i = 0; i < A.id_column_for_A_down.size(); i++)
		{
			math::MakeCopyVector_A_into_B(A.id_column_for_A_down[i], this->id_column_for_A_down[i]);
		}
	}

	std::vector<int>* GetIDsColumnsInDownString(int id_string)
	{
		return &(this->id_column_for_A_down[id_string]);
	}
	std::vector<int>* GetIDsColumnsInUpString(int id_string)
	{
		return &(this->id_column_for_A_up[id_string]);
	}

	void Initialization(std::vector<std::vector<int>> &id_column_for_A_up, std::vector<std::vector<int>> &id_column_for_A_down)
	{
		int string_count = (int)id_column_for_A_up.size();

		this->id_column_for_A_up.resize(string_count);
		this->id_column_for_A_down.resize(string_count);
		this->Diag.resize(string_count);
		this->A_up.resize(string_count);
		this->A_down.resize(string_count);
		this->X.resize(string_count);
		this->F.resize(string_count);

		for (int i = 0; i < string_count; i++)
		{
			this->A_up[i].resize(id_column_for_A_up[i].size());
			this->A_down[i].resize(id_column_for_A_down[i].size());
			
			math::MakeCopyVector_A_into_B(id_column_for_A_up[i], this->id_column_for_A_up[i]);
			math::MakeCopyVector_A_into_B(id_column_for_A_down[i], this->id_column_for_A_down[i]);
		}
	}
	void UpdateEquation(int id_equation, TForMatrix &diag, TForMatrix &non_diag, TForVector &X, TForVector &F)
	{
		this->Diag[id_equation] = diag;
		this->X[id_equation] = X;
		this->F[id_equation] = F;
		for (int j = 0; j < this->A_up[id_equation].size(); j++)
		{
			this->A_up[id_equation][j] = non_diag;
		}
		for (int j = 0; j < this->A_down[id_equation].size(); j++)
		{
			this->A_down[id_equation][j] = non_diag;
		}
	}
	void Symmetrization()
	{
		this->id_column_for_A_down.resize(this->Diag.size());
		this->A_down.resize(this->Diag.size());
		for (int i_up = 0; i_up < this->id_column_for_A_up.size(); i_up++)
		{
			for (int jj_up = 0; jj_up < this->id_column_for_A_up[i_up].size(); jj_up++)
			{
				int j_up = this->id_column_for_A_up[i_up][jj_up];
				this->id_column_for_A_down[j_up].push_back(i_up);
				this->A_down[j_up].push_back(this->A_up[i_up][jj_up]);
			}
		}
	}
	void RandomMatrix(int N, double X_val)
	{
		this->Diag.resize(N);
		this->X.resize(N);
		this->F.resize(N);
		this->id_column_for_A_up.resize(N);
		this->A_up.resize(N);

		for (int i = 0; i < N; i++)
		{
			this->Diag[i] = rand() * 1000.0 / RAND_MAX;
			this->X[i] = X_val;

			int num_elem_in_string = (rand() % (N - i)) % (N/10);
			if (num_elem_in_string != 0)
			{
				this->id_column_for_A_up[i].resize(num_elem_in_string);
				this->A_up[i].resize(num_elem_in_string);
				for (int jj = 0; jj < num_elem_in_string; jj++)
				{
					this->id_column_for_A_up[i][jj] = rand() % (N - i - 1) + (i + 1);
					this->A_up[i][jj] = rand() * 1.0 / RAND_MAX;
				}

				bool is_duplicate = false;
				do {
					is_duplicate = false;
					math::MakeQuickSort(this->id_column_for_A_up[i]);
					for (int ii = 0; ii < this->id_column_for_A_up[i].size() - 1; ii++)
					{
						if (this->id_column_for_A_up[i][ii] == this->id_column_for_A_up[i][ii + 1])
						{
							is_duplicate = true;
							this->id_column_for_A_up[i][ii + 1] = rand() % (N - i - 1) + (i + 1);
							break;
						}
					}
				} while (is_duplicate == true);
			}
		}

		this->Symmetrization();
		this->MultiplicationMatrixVector(this->X, this->F);
	}
	void ClearVariables()
	{
		math::InitializationVector(this->Diag, 0);
		math::InitializationVector(this->A_up, 0);
		math::InitializationVector(this->A_down, 0);
		math::InitializationVector(this->X, 0);
		math::InitializationVector(this->F, 0);
	}
	bool SetValue(int id_string, int id_column, TForMatrix& val)
	{
		if (id_string >= this->GetMatrixSize() || id_column >= this->GetMatrixSize())
			return false;

		if (id_column == id_string)
		{
			this->Diag[id_string] = val;
			return true;
		}
		if (id_column > id_string)
		{
			int jj = math::GetPositionInSortVector(this->id_column_for_A_up[id_string], id_column);
			if (jj >= 0)
			{
				this->A_up[id_string][jj] = val;
				return true;
			}
		}
		if (id_column < id_string)
		{
			int jj = math::GetPositionInSortVector(this->id_column_for_A_down[id_string], id_column);
			if (jj >= 0)
			{
				this->A_down[id_string][jj] = val;
				return true;
			}
		}

		return false;
	}
	bool GetValue(int id_string, int id_column, TForMatrix& val)
	{
		if (id_string >= this->GetMatrixSize() || id_column >= this->GetMatrixSize())
			return false;

		if (id_column == id_string)
		{
			val = this->Diag[id_string];
			return true;
		}
		if (id_column > id_string)
		{
			int jj = math::GetPositionInSortVector(this->id_column_for_A_up[id_string], id_column);
			if (jj >= 0)
			{
				val = this->A_up[id_string][jj];
				return true;
			}
		}
		if (id_column < id_string)
		{
			int jj = math::GetPositionInSortVector(this->id_column_for_A_down[id_string], id_column);
			if (jj >= 0)
			{
				val = this->A_down[id_string][jj];
				return true;
			}
		}

		return false;
	}

	void SummPartOfMatrix(DenseMatrix<TForMatrix, TForVector> &matrix, std::vector<int> &id_row)
	{
		try
		{
			for (int i = 0; i < matrix.GetSize(); i++)
			{
				this->F[id_row[i]] += matrix.GetElementOfF(i);
				this->X[id_row[i]] += matrix.GetElementOfX(i);
				for (int j = 0; j < matrix.GetSize(); j++)
				{
					int _id_string = id_row[i];
					int _id_row = id_row[j];

					bool find = false;
					if (_id_string > _id_row)
					{
						for (int jj = 0; jj < this->id_column_for_A_down[_id_string].size() /*&& j >= this->id_column_for_A_down[_id_string][jj]*/; jj++)
						{
							if (this->id_column_for_A_down[_id_string][jj] == _id_row)
							{
								A_down[_id_string][jj] += matrix.GetElementOfA(i, j);
								find = true;
							}
						}

					}
					else {
						if (_id_string < _id_row)
						{
							for (int jj = 0; jj < this->id_column_for_A_up[_id_string].size() /*&& j >= this->id_column_for_A_up[_id_string][jj]*/; jj++)
							{
								if (this->id_column_for_A_up[_id_string][jj] == _id_row)
								{
									A_up[_id_string][jj] += matrix.GetElementOfA(i, j);
									find = true;
								}
							}
						}
						else {
							this->Diag[_id_string] += matrix.GetElementOfA(i, i);
							find = true;
						}
					}

					if (find == false)
					{
						printf_s("We don’t find the element [%d][%d].\n", _id_string, _id_row);
					}
				}
			}
			return;
		}
		catch (const std::exception&)
		{
			printf_s("Error: Library/CSSD_Matrix.h/void SummPartOfMatrix(DenseMatrix<Tensor2Rank3D, Point<double>> &matrix, std::vector<int> *id_row)\n");
		}
	}
	void SummPartOfVector(std::vector<TForVector> &vector, std::vector<int> &id_row)
	{
		try
		{
			for (int i = 0; i < vector.size(); i++)
			{
				this->F[id_row[i]] += vector[i];
				//this->X[id_row[i]] += vector[i];
			}
			return;
		}
		catch (const std::exception&)
		{
			printf_s("Error: Library/CSSD_Matrix.h/void SummPartOfMatrix(DenseMatrix<Tensor2Rank3D, Point<double>> &matrix, std::vector<int> *id_row)\n");
		}
	}


	void MultiplicationMatrixVector(std::vector<TForVector> &B, std::vector<TForVector> &result)
	{
		try
		{
			if (result.size() != B.size()) result.resize(B.size());
			omp_set_num_threads(math::NUM_THREADS);
#pragma omp parallel for schedule(dynamic) 
			for (int i = 0; i < result.size(); i++)
			{
				TForVector tmp;
				tmp = 0.0;
				for (int j = 0; j < A_up[i].size(); j++)
				{
					tmp += A_up[i][j] * B[id_column_for_A_up[i][j]];
				}
				for (int j = 0; j < A_down[i].size(); j++)
				{
					tmp += A_down[i][j] * B[id_column_for_A_down[i][j]];
				}
				tmp += Diag[i] * B[i];
				result[i] = tmp;
			}
		}
		catch (const std::exception&)
		{
			printf_s("Error: Library/Math.h/void MultiplicationMatrixVector(std::vector<std::vector<T1>> &A, std::vector<T2> &B, std::vector<T2> &result)\n");

		}
	}
	void MakeDivideOnSqrtFormSumm()
	{
		for (int i = 0; true && i < this->GetMatrixSize(); i++)
		{
			TForVector koef;
			koef = 0.0;
			koef += Diag[i] * Diag[i];

			for (int j = 0; j < this->A_up[i].size(); j++)
			{
				koef += A_up[i][j] * A_up[i][j];
			}
			for (int j = 0; j < this->A_down[i].size(); j++)
			{
				koef += A_down[i][j] * A_down[i][j];
			}

			koef = math::SolveSqrt(koef);

			X[i] /= koef;
			F[i] /= koef;

			Diag[i] /= koef;

			for (int j = 0; j < this->A_up[i].size(); j++)
			{
				A_up[i][j] /= koef;
			}
			for (int j = 0; j < this->A_down[i].size(); j++)
			{
				A_down[i][j] /= koef;
			}
		}
	}
	void PrecondorSSOR(double w, CSSD_Matrix<double, double>& A)
	{
		this->Diag.resize(A.Diag.size());
		this->A_up.resize(A.A_up.size());
		for (int i = 0; i < this->A_up.size(); i++)
		{
			this->A_up[i].resize(A.A_up[i].size());
		}
		math::MakeCopyVector_A_into_B(A.id_column_for_A_up, this->id_column_for_A_up);
		math::MakeCopyVector_A_into_B(A.A_up, this->A_up);
		if (A.A_down.size() == 0)
		{
			this->A_down.resize(A.A_up.size());
			this->id_column_for_A_down.resize(A.A_up.size());
			for (int i_down = 1; i_down < A_down.size(); i_down++)
			{
				for (int i_up = 0; i_up < i_down; i_up++)
				{
					for (int j_up = 0; j_up < A.id_column_for_A_up[i_up].size(); j_up++)
					{
						if (i_down == A.id_column_for_A_up[i_up][j_up])
						{
							this->id_column_for_A_down[i_down].push_back(i_up);
							this->A_down[i_down].push_back(A.A_up[i_up][j_up]);
						}
						if (i_down < A.id_column_for_A_up[i_up][j_up]) break;
					}
				}
			}
		}
		else
		{
			this->A_down.resize(A.A_down.size());
			for (int i = 0; i < this->A_down.size(); i++)
			{
				this->A_down[i].resize(A.A_down[i].size());
			}
			math::MakeCopyVector_A_into_B(A.id_column_for_A_down, this->id_column_for_A_down);
			math::MakeCopyVector_A_into_B(A.A_down, this->A_down);
		}

		for (int I = 0; I < A.GetMatrixSize(); I++)
		{
			this->Diag[I] = 1.0 / A.Diag[I] / w / (2 - w);
			//this->Diag[I] = A.Diag[I];

			for (int j_up = 0; j_up < this->id_column_for_A_up[I].size(); j_up++)
			{
				this->A_up[I][j_up] = -1 * w * this->A_up[I][j_up];
			}
			for (int j_down = 0; j_down < this->id_column_for_A_down[I].size(); j_down++)
			{
				this->A_down[I][j_down] = -1 * w * this->A_down[I][j_down];
			}
		}
	}
	void SolvePrecondorSSOR_SLAE(std::vector<double>& free_diag, std::vector<double>& X, std::vector<double>& Z)
	{
		X.resize(Z.size());

		//STEP 1: Ex=f
		for (int I = this->GetMatrixSize()-1; I >=0 ; I--)
		{
			X[I] = Z[I];
			for (int j_up = 0; j_up < this->id_column_for_A_up[I].size(); j_up++)
			{
				int J = this->id_column_for_A_up[I][j_up];
				X[I] -= X[J] * this->A_up[I][j_up];
			}
			X[I] /= free_diag[I];
		}

		//STEP 2: Dx_2=x_1
		for (int I = 0; I < this->GetMatrixSize(); I++)
		{
			X[I] = X[I] / this->Diag[I];
		}

		//STEP 3: Fx_3=x_2
		for (int I = 0; I < this->GetMatrixSize(); I++)
		{
			X[I] = X[I];
			for (int j_down = 0; j_down < this->id_column_for_A_down[I].size(); j_down++)
			{
				int J = this->id_column_for_A_down[I][j_down];
				X[I] -= X[J] * this->A_down[I][j_down];
			}
			X[I] /= free_diag[I];
		}
	}
	void PrecondorSSOR_summetric(double w, CSSD_Matrix<double, double> &A)
	{
		this->Diag.resize(A.Diag.size());
		this->A_up.resize(A.A_up.size());
		for (int i = 0; i < this->A_up.size(); i++)
		{
			this->A_up[i].resize(A.A_up[i].size());
		}
		math::MakeCopyVector_A_into_B(A.id_column_for_A_up, this->id_column_for_A_up);
		printf_s("\nPrecondorSSOR_summetric\n");

		/*for (int i = 0; i < this->Diag.size(); i++)
		{
			this->Diag[i] = A.Diag[i];
			for (int jj = 0; jj < A.id_column_for_A_up[i].size(); jj++)
			{
				this->A_up[i][jj] = w * A.A_up[i][jj];
				
				int j = A.id_column_for_A_up[j][jj];
				if (j > i)
				{
					this->Diag[i] += w * w / A.Diag[j] * A.A_up[i][jj] * A.A_up[i][jj];

					for (int ii = 0; ii < A.id_column_for_A_up[j].size(); ii++)
						if (A.id_column_for_A_up[j][ii] == j)
						{
							this->A_up[i][jj] += w * w / A.Diag[j] * A.A_up[i][jj] A.A_up[j][ii];
							break;
						}
				}
			}
		}*/
//#pragma omp parallel for schedule(dynamic) 
		for (int i = 0; i < this->Diag.size(); i++)
		{
			if(omp_get_thread_num() == 0 && i %1000 == 0)
				printf_s("\t\t%d / %d\r", i, this->Diag.size());

			this->Diag[i] = A.Diag[i];
			for (int p = i + 1, pp = 0; p < this->Diag.size(); p++)
			{
				for (; pp < this->id_column_for_A_up[i].size(); pp++)
				{
					if (this->id_column_for_A_up[i][pp] == p)
					{
						this->Diag[i] += (w * w / A.Diag[p]) * A.A_up[i][pp] * A.A_up[i][pp];
						break;
					}
				}
			}

			for (int j = i+1, jj = 0; j < this->Diag.size(); j++)
			{
				for (; jj < this->id_column_for_A_up[i].size(); jj++)
				{
					if (this->id_column_for_A_up[i][jj] == j)
					{
						this->A_up[i][jj] = w * A.A_up[i][jj];

						for (int p = j + 1; p < this->Diag.size(); p++)
						{
							for (int pp_i = 0; pp_i < this->id_column_for_A_up[i].size(); pp_i++)
							{
								if (this->id_column_for_A_up[i][pp_i] == p)
								{
									for (int pp_j = 0; pp_j < this->id_column_for_A_up[j].size(); pp_j++)
									{
										if (this->id_column_for_A_up[j][pp_j] == p)
										{
											this->A_up[i][jj] += (w * w / A.Diag[p]) * A.A_up[i][pp_i] * A.A_up[j][pp_j];
											break;
										}
									}
									break;
								}
							}
						}
						break;
					}
				}
			}
		}
		this->Symmetrization();
	}


	double BCG_Stab2(int maxiter, double E)
	{
		int N = this->GetMatrixSize();

		double enorma, fnorma;
		double eps3 = 1.E-16;
		std::vector <double> r;
		std::vector <double> rt;
		std::vector <double> p;
		std::vector <double> s;
		std::vector <double> q;
		std::vector <double> h;
		std::vector <double> x0;
		std::vector <double> x;

		r.resize(N);
		rt.resize(N);
		p.resize(N);
		s.resize(N);
		q.resize(N);
		h.resize(N);
		x0.resize(N);

		//!нулевая итерация
		this->MultiplicationMatrixVector(this->X, q);
		for (auto i = 0; i < N; i++)
		{
			r[i] = this->F[i] - q[i];
			rt[i] = r[i]; // !1d0;
			p[i] = r[i];
		}

		double f_norma = math::MakeInnerProduct(this->F, this->F);
		f_norma = sqrt(f_norma);

		double a_norma = 1;

		//!основной цикл
		int i = 0;
		for (i = 0; i < maxiter; i++)
		{
			double a1 = math::MakeInnerProduct(r, rt);
			this->MultiplicationMatrixVector(p, q);//!A * p = > q
			double a2 = math::MakeInnerProduct(q, rt);
			double alfa = a1 / a2;
			for (auto j = 0; j < N; j++)
				s[j] = r[j] - alfa * q[j];
			
			this->MultiplicationMatrixVector(s, h); //!A * s = > h
			a2 = math::MakeInnerProduct(s, h);
			double a3 = math::MakeInnerProduct(h, h);
			double w = a2 / a3;

			for (auto j = 0; j < N; j++)
			{
				x0[j] = this->X[j];
				this->X[j] += alfa * p[j] + w * s[j];
				r[j] = s[j] - w * h[j];
			}

			a2 = math::MakeInnerProduct(r, rt);
			double beta = a2 / a1 * alfa / w;

			for (auto j = 0; j < N; j++)
				p[j] = r[j] + beta * (p[j] - w * q[j]);

			double a_norma_new = math::MakeInnerProduct(r, r);
			a_norma_new = sqrt(a_norma_new) / f_norma;
			/*if (a_norma_new > a_norma*10)
			{
				a_norma = a_norma_new;
				break;
			}*/
			a_norma = a_norma_new;

			if (a_norma < E) break;
			if (abs(beta) <= eps3)
			{
				printf("\tBSGstab\t>\tNO SOLUTION\t>\t%d\t-\t%.5e\n", i, a_norma);
				break;
			}
			if (i % 1000 == 0 && this->print_logs == true)
			{
				printf_s("\tBSGstab\t\t>\t%d\t-\t%.5e\n", i, a_norma);
			}
		}

		this->MultiplicationMatrixVector(this->X, r);
		double temp = 0;
		for (int i = 0; i < this->GetMatrixSize(); i++)
			temp += (F[i] - r[i]) * (F[i] - r[i]);
		double result = sqrt(temp) / f_norma;
		if (this->print_logs) printf("\tBSGstab\t\t>\t%d\t-\t%.5e (%.5e resolve residual)\n", i, a_norma, result);

		return result;
	}
	double BiCG_Stab(int maxiter, double E)
	{
		//------>>>>>!!!!!
		//print_logs = false;

		if (print_logs == true)
		{
			if (log_out == NULL) print_logs = false;
			else fprintf_s(log_out, "Matrix size %d\n", this->GetMatrixSize());
		}

		std::vector<TForVector> residual, mult_result, mult_result2;
		double norma_residual_curr = 0, norma_residual_prev = 0, norma_F = 0;
		residual.resize(this->GetMatrixSize());
		mult_result.resize(this->GetMatrixSize());
		mult_result2.resize(this->GetMatrixSize());

		///----->>>
		/*for (int i = 0; i < this->GetMatrixSize(); i++)
		{
			math::InitializationVector(A_up[i], i*1.0);
			math::InitializationVector(A_down[i], i*1.0);
		}
		math::InitializationVector(Diag, 1000.0);
		math::InitializationVector(X, 1.0);
		this->MultiplicationMatrixVector(this->X, F);*/
		//math::InitializationVector(X, 0.0);
		///----->>>

		///->>>>>
		//MakeDivideOnSqrtFormSumm();
		///->>>>

		//initialisation
		this->MultiplicationMatrixVector(this->X, mult_result);
		math::MakeSummVectors(this->F, 1.0, mult_result, -1.0, residual);
		norma_residual_curr = math::MakeInnerProduct(residual, residual);
		//norma_residual_prev = 1;
		norma_F = math::MakeInnerProduct(this->F, this->F);

		double betta, alpha = 1, omega = 1;
		std::vector<TForVector> P(this->GetMatrixSize()), 
			V(this->GetMatrixSize()),
			S(this->GetMatrixSize()),
			T(this->GetMatrixSize()),
			residual_add(this->GetMatrixSize());
		
		math::MakeCopyVector_A_into_B(residual, residual_add);
		math::MakeCopyVector_A_into_B(residual, P);
		//math::InitializationVector(P, 1.0);
		//math::InitializationVector(V, 0.0);
		//math::InitializationVector(S, 0.0);
		//math::InitializationVector(T, 0.0);

		int k;
		double real_residual_prev = norma_residual_curr;
		double real_residual = norma_residual_curr;
		for (k = 1; k < maxiter && sqrt(norma_residual_curr / norma_F) > E && norma_residual_curr == norma_residual_curr; k++)
		{
			///solve Vk = A*Pk
			this->MultiplicationMatrixVector(P, V);
			//solve alpha_k = (r_k, r_add) / (A*Pk, r_add)
			norma_residual_prev = math::MakeInnerProduct(residual, residual_add);
			alpha = norma_residual_prev / math::MakeInnerProduct(residual_add, V);
			//solve Sk = residual - alpha*(A*Pk)
			math::MakeSummVectors(residual, 1.0, V, -alpha, S);
			///solve Tk = A*Sk
			this->MultiplicationMatrixVector(S, T);
			//solve omega = (A*Sk, Sk) / (A*Sk, A*Sk)
			omega = math::MakeInnerProduct(T, S) / math::MakeInnerProduct(T, T);
			//solve Xk = X(k-1) + omega*Sk + alpha*Pk
			math::MakeSummVectors(S, omega, P, alpha, mult_result);
			math::MakeSummVectors(this->X, 1.0, mult_result, 1.0, this->X);
			//residual=F-A*Xk
			this->MultiplicationMatrixVector(this->X, mult_result);
			math::MakeSummVectors(this->F, 1.0, mult_result, -1.0, residual);
			//betta=(residual_k, residual_add)/(residuall_k-1, residual_add)
			norma_residual_curr = math::MakeInnerProduct(residual, residual_add);
			//betta = (norma_residual_curr) / (norma_residual_prev);
			betta = ((norma_residual_curr) / (norma_residual_prev)) * (alpha / omega);
			//solve Pk = residual + betta * (P(k-1) - omega*V(k-1))
			math::MakeSummVectors(P, betta, V, -betta*omega, mult_result);
			math::MakeSummVectors(residual, 1, mult_result, 1, P);

			//solve residual = Sk - omega*Tk
			//math::MakeSummVectors(S, 1.0, T, -omega, residual);

			//norma_residual_prev = norma_residual_curr;
			//norma_residual_curr = math::MakeInnerProduct(residual_add, residual);

			real_residual = math::MakeInnerProduct(residual, residual);
			if (real_residual > real_residual_prev)
			{
				break;
			}
			else {
				real_residual_prev = real_residual;
			}

			if (print_logs && (k % 1000 == 0 || k == 1))
			{
				//double scal = math::MakeInnerProduct(residual, residual);
				printf_s(">%d -  %.5e \t(real: %.5e)\n", k, sqrt(norma_residual_curr / norma_F), sqrt(real_residual / norma_F));
			}
			if (print_logs)
			{
				fprintf_s(log_out, "%d\t%.15e\t(real: %.5e)\n", k, sqrt(norma_residual_curr / norma_F), sqrt(real_residual / norma_F));
			}
		}
		
		this->MultiplicationMatrixVector(this->X, mult_result);
		double temp = 0;
		for (int i = 0; i < this->GetMatrixSize(); i++)
			temp += (F[i] - mult_result[i])*(F[i] - mult_result[i]);
		double result = sqrt(temp / norma_F);
		printf(">REAL RESULT> %d - %.15e\n", k, result);

		if (print_logs) fprintf_s(log_out, "REAL RESULT %.15e\n\n", result);
		if (k == 1)
		{
			result *= -1;
		}
		return result;
	}
	double BiCG_Stab_blocks(int maxiter, double E)
	{
		if (print_logs == true)
		{
			if (log_out == NULL) print_logs = false;
			else fprintf_s(log_out, "Matrix size %d\n", this->GetMatrixSize());
		}

		std::vector<TForVector> residual, mult_result, mult_result2;
		double norma_residual_curr = 0, norma_residual_prev = 0, norma_F = 0;
		residual.resize(this->GetMatrixSize());
		mult_result.resize(this->GetMatrixSize());
		mult_result2.resize(this->GetMatrixSize());

		///----->>>
		/*for (int i = 0; i < this->GetMatrixSize(); i++)
		{
			math::InitializationVector(A_up[i], i*1.0);
			math::InitializationVector(A_down[i], i*1.0);
		}
		math::InitializationVector(Diag, 1000.0);
		math::InitializationVector(X, 1.0);
		this->MultiplicationMatrixVector(this->X, F);*/
		//math::InitializationVector(X, 0.0);
		///----->>>

		///->>>>>
		//MakeDivideOnSqrtFormSumm();
		for (int i = 0; true && i < this->GetMatrixSize(); i++)
		{
			Point<double> koef;
			koef.x += Diag[i].val[0][0] * Diag[i].val[0][0];
			koef.x += Diag[i].val[0][1] * Diag[i].val[0][1];
			koef.x += Diag[i].val[0][2] * Diag[i].val[0][2];

			koef.y += Diag[i].val[1][0] * Diag[i].val[1][0];
			koef.y += Diag[i].val[1][1] * Diag[i].val[1][1];
			koef.y += Diag[i].val[1][2] * Diag[i].val[1][2];

			koef.z += Diag[i].val[2][0] * Diag[i].val[2][0];
			koef.z += Diag[i].val[2][1] * Diag[i].val[2][1];
			koef.z += Diag[i].val[2][2] * Diag[i].val[2][2];

			for (int j = 0; j < this->A_up[i].size(); j++)
			{
				koef.x += A_up[i][j].val[0][0] * A_up[i][j].val[0][0];
				koef.x += A_up[i][j].val[0][1] * A_up[i][j].val[0][1];
				koef.x += A_up[i][j].val[0][2] * A_up[i][j].val[0][2];

				koef.y += A_up[i][j].val[1][0] * A_up[i][j].val[1][0];
				koef.y += A_up[i][j].val[1][1] * A_up[i][j].val[1][1];
				koef.y += A_up[i][j].val[1][2] * A_up[i][j].val[1][2];

				koef.z += A_up[i][j].val[2][0] * A_up[i][j].val[2][0];
				koef.z += A_up[i][j].val[2][1] * A_up[i][j].val[2][1];
				koef.z += A_up[i][j].val[2][2] * A_up[i][j].val[2][2];
			}
			for (int j = 0; j < this->A_down[i].size(); j++)
			{
				koef.x += A_down[i][j].val[0][0] * A_down[i][j].val[0][0];
				koef.x += A_down[i][j].val[0][1] * A_down[i][j].val[0][1];
				koef.x += A_down[i][j].val[0][2] * A_down[i][j].val[0][2];

				koef.y += A_down[i][j].val[1][0] * A_down[i][j].val[1][0];
				koef.y += A_down[i][j].val[1][1] * A_down[i][j].val[1][1];
				koef.y += A_down[i][j].val[1][2] * A_down[i][j].val[1][2];

				koef.z += A_down[i][j].val[2][0] * A_down[i][j].val[2][0];
				koef.z += A_down[i][j].val[2][1] * A_down[i][j].val[2][1];
				koef.z += A_down[i][j].val[2][2] * A_down[i][j].val[2][2];
			}

			koef.x = sqrt(koef.x);
			koef.y = sqrt(koef.y);
			koef.z = sqrt(koef.z);

			X[i].x /= koef.x;
			X[i].y /= koef.y;
			X[i].z /= koef.z;

			F[i].x /= koef.x;
			F[i].y /= koef.y;
			F[i].z /= koef.z;

			Diag[i].val[0][0] /= koef.x;
			Diag[i].val[0][1] /= koef.x;
			Diag[i].val[0][2] /= koef.x;

			Diag[i].val[1][0] /= koef.y;
			Diag[i].val[1][1] /= koef.y;
			Diag[i].val[1][2] /= koef.y;

			Diag[i].val[2][0] /= koef.z;
			Diag[i].val[2][1] /= koef.z;
			Diag[i].val[2][2] /= koef.z;

			for (int j = 0; j < this->A_up[i].size(); j++)
			{
				A_up[i][j].val[0][0] /= koef.x;
				A_up[i][j].val[0][1] /= koef.x;
				A_up[i][j].val[0][2] /= koef.x;

				A_up[i][j].val[1][0] /= koef.y;
				A_up[i][j].val[1][1] /= koef.y;
				A_up[i][j].val[1][2] /= koef.y;

				A_up[i][j].val[2][0] /= koef.z;
				A_up[i][j].val[2][1] /= koef.z;
				A_up[i][j].val[2][2] /= koef.z;
			}
			for (int j = 0; j < this->A_down[i].size(); j++)
			{
				A_down[i][j].val[0][0] /= koef.x;
				A_down[i][j].val[0][1] /= koef.x;
				A_down[i][j].val[0][2] /= koef.x;

				A_down[i][j].val[1][0] /= koef.y;
				A_down[i][j].val[1][1] /= koef.y;
				A_down[i][j].val[1][2] /= koef.y;

				A_down[i][j].val[2][0] /= koef.z;
				A_down[i][j].val[2][1] /= koef.z;
				A_down[i][j].val[2][2] /= koef.z;
			}
		}
		///->>>>

		//initialisation
		this->MultiplicationMatrixVector(this->X, mult_result);
		math::MakeSummVectors(this->F, 1.0, mult_result, -1.0, residual);
		norma_residual_curr = math::MakeInnerProduct(residual, residual);
		//norma_residual_prev = 1;
		norma_F = math::MakeInnerProduct(this->F, this->F);

		double betta, alpha = 1, omega = 1;
		std::vector<TForVector> P(this->GetMatrixSize()),
			V(this->GetMatrixSize()),
			S(this->GetMatrixSize()),
			T(this->GetMatrixSize()),
			residual_add(this->GetMatrixSize());

		math::MakeCopyVector_A_into_B(residual, residual_add);
		math::MakeCopyVector_A_into_B(residual, P);
		//math::InitializationVector(P, 1.0);
		//math::InitializationVector(V, 0.0);
		//math::InitializationVector(S, 0.0);
		//math::InitializationVector(T, 0.0);

		int k;
		for (k = 1; k < maxiter && sqrt(norma_residual_curr / norma_F) > E; k++)
		{
			///solve Vk = A*Pk
			this->MultiplicationMatrixVector(P, V);
			//solve alpha_k = (r_k, r_add) / (A*Pk, r_add)
			norma_residual_prev = math::MakeInnerProduct(residual, residual_add);
			alpha = norma_residual_prev / math::MakeInnerProduct(residual_add, V);
			//solve Sk = residual - alpha*(A*Pk)
			math::MakeSummVectors(residual, 1.0, V, -alpha, S);
			///solve Tk = A*Sk
			this->MultiplicationMatrixVector(S, T);
			//solve omega = (A*Sk, Sk) / (A*Sk, A*Sk)
			omega = math::MakeInnerProduct(T, S) / math::MakeInnerProduct(T, T);
			//solve Xk = X(k-1) + omega*Sk + alpha*Pk
			math::MakeSummVectors(S, omega, P, alpha, mult_result);
			math::MakeSummVectors(this->X, 1.0, mult_result, 1.0, this->X);
			//residual=F-A*Xk
			this->MultiplicationMatrixVector(this->X, mult_result);
			math::MakeSummVectors(this->F, 1.0, mult_result, -1.0, residual);
			//betta=(residual_k, residual_add)/(residuall_k-1, residual_add)
			norma_residual_curr = math::MakeInnerProduct(residual, residual_add);
			//betta = (norma_residual_curr) / (norma_residual_prev);
			betta = ((norma_residual_curr) / (norma_residual_prev)) * (alpha / omega);
			//solve Pk = residual + betta * (P(k-1) - omega*V(k-1))
			math::MakeSummVectors(P, betta, V, -betta * omega, mult_result);
			math::MakeSummVectors(residual, 1, mult_result, 1, P);

			//solve residual = Sk - omega*Tk
			//math::MakeSummVectors(S, 1.0, T, -omega, residual);

			//norma_residual_prev = norma_residual_curr;
			//norma_residual_curr = math::MakeInnerProduct(residual_add, residual);

			if (k % 1000 == 0 || k == 1)
			{
				double scal = math::MakeInnerProduct(residual, residual);
				printf_s(">%d -  %.5e \t(real: %.5e)\n", k, sqrt(norma_residual_curr / norma_F), sqrt(scal / norma_F));
			}
			if (print_logs) fprintf_s(log_out, "%d\t%.15e\n", k, sqrt(norma_residual_curr / norma_F));
		}

		this->MultiplicationMatrixVector(this->X, mult_result);
		double temp = 0;
		for (int i = 0; i < this->GetMatrixSize(); i++)
			temp += (F[i] - mult_result[i]) * (F[i] - mult_result[i]);
		double result = sqrt(temp / norma_F);
		printf(">REAL RESULT> %d - %.15e\n", k, result);

		if (print_logs) fprintf_s(log_out, "REAL RESULT %.15e\n", result);
		if (k == 1)
		{
			result *= -1;
		}
		return result;
	}
	double BiCG_Stab_Preconditioning(int maxiter, double E, CSSD_Matrix<double, double>& M)
	{
		//------>>>>>!!!!!
		//print_logs = false;

		if (print_logs == true)
		{
			if (log_out == NULL) print_logs = false;
			else fprintf_s(log_out, "Matrix size %d\n", this->GetMatrixSize());
		}

		double alpha = 1, betta = 1, omega = 1, norma_F, norma_residual, norma_residual_old;
		std::vector<TForVector> residual, residual_0, mult_result, mult_result2;
		residual.resize(this->GetMatrixSize());
		residual_0.resize(this->GetMatrixSize());
		mult_result.resize(this->GetMatrixSize());
		mult_result2.resize(this->GetMatrixSize());

		//initialisation
		this->MultiplicationMatrixVector(this->X, mult_result);
		math::MakeSummVectors(this->F, 1.0, mult_result, -1.0, residual);
		
		std::vector<TForVector> P(this->GetMatrixSize()), P_add(this->GetMatrixSize()),
			S(this->GetMatrixSize()), S_add(this->GetMatrixSize());

		math::MakeCopyVector_A_into_B(residual, residual_0);
		math::MakeCopyVector_A_into_B(residual, P);
		
		norma_F = math::MakeInnerProduct(this->F, this->F);
		norma_residual = math::MakeInnerProduct(residual, residual_0);
		norma_residual_old = norma_residual;

		int k;
		for (k = 1; k < maxiter && norma_residual == norma_residual && sqrt(norma_residual/norma_F) > 1e-16; k++)
		{
			//solve M*p_add=p
			if (!(k % 10) || k == 1)
				M.MSG(maxiter, 1e-10, P_add, P, true);
			else
				M.MSG(maxiter, 1e-10, P_add, P, false);
			//solve mult_result = A*P_add;
			this->MultiplicationMatrixVector(P_add, mult_result);
			//solve alpha_k = (r_k, r_0) / (A*P_add, r_0)
			alpha = math::MakeInnerProduct(residual, residual_0) / math::MakeInnerProduct(residual_0, mult_result);
			//solve Sk = residual - alpha*mult_result
			math::MakeSummVectors(residual, 1.0, mult_result, -alpha, S);
			//solve M*S_add=S
			if (!(k % 10) || k == 1)
				M.MSG(maxiter, 1e-10, S_add, S, true);
			else
				M.MSG(maxiter, 1e-10, S_add, S, false);
			//solve mult_result = A*S_add;
			this->MultiplicationMatrixVector(S_add, mult_result2);
			//solve omega = (A*S_add, S) / (A*S_add, A*S_add)
			omega = math::MakeInnerProduct(mult_result2, S) / math::MakeInnerProduct(mult_result2, mult_result2);
			//solve Xk = X(k-1) + alpha*P_add + omega*S_add
			math::MakeSummVectors(this->X, 1.0, P_add, alpha, this->X);
			math::MakeSummVectors(this->X, 1.0, S_add, omega, this->X);
			//solve r_new = s - omega * A*S_add
			math::MakeSummVectors(S, 1.0, mult_result2, -1 * omega, residual);
			norma_residual_old = norma_residual;
			norma_residual = math::MakeInnerProduct(residual, residual_0);
			//solve betta = (r_new,r_0) / (r_old, r_0) * alpha/omega
			betta = (norma_residual / norma_residual_old) * (alpha / omega);
			//solve P = r + betta*p - betta*omega*A*P_add
			math::MakeSummVectors(residual, 1, P, betta, P);
			math::MakeSummVectors(P, 1, mult_result, -1 * betta * omega, P);

			if (/*print_logs &&*/ (k % 10 == 0 || k == 1))
			{
				this->MultiplicationMatrixVector(this->X, mult_result);
				double temp = 0;
				for (int i = 0; i < this->GetMatrixSize(); i++)
					temp += (F[i] - mult_result[i]) * (F[i] - mult_result[i]);
				double real_residual = sqrt(temp / norma_F);

				if (real_residual <= E) break;

				printf_s(">%d -  %.5e \t(real: %.5e)\n", k, sqrt(norma_residual / norma_F), real_residual);
			}
			if (print_logs)
			{
				this->MultiplicationMatrixVector(this->X, mult_result);
				double temp = 0;
				for (int i = 0; i < this->GetMatrixSize(); i++)
					temp += (F[i] - mult_result[i]) * (F[i] - mult_result[i]);
				double real_residual = sqrt(temp / norma_F);

				fprintf_s(log_out, "%d\t%.15e\t(real: %.5e)\n", k, sqrt(norma_residual / norma_F), real_residual);
			}
		}

		this->MultiplicationMatrixVector(this->X, mult_result);
		double temp = 0;
		for (int i = 0; i < this->GetMatrixSize(); i++)
			temp += (F[i] - mult_result[i]) * (F[i] - mult_result[i]);
		double result = sqrt(temp / norma_F);
		printf(">REAL RESULT> %d - %.15e\n", k, result);

		if (print_logs) fprintf_s(log_out, "REAL RESULT %.15e\n\n", result);
		if (k == 1)
		{
			result *= -1;
		}
		return result;
	}

	double MSG_PreconditioningSSOR(int maxiter, double E, CSSD_Matrix<double, double> &Precondor)
	{
		int k, n_matr = this->GetMatrixSize();
		std::vector<double> r(n_matr); //невязка
		std::vector<double> z(n_matr); //невязка предобусловленной СЛАУ
		std::vector<double> p(n_matr);
		std::vector<double> temp_mult(n_matr);

		//подготовка к итерациям
		k = 0;
		this->MultiplicationMatrixVector(this->X, temp_mult);
		for (int i = 0; i < n_matr; i++)
		{
			r[i] = this->F[i] - temp_mult[i];
			if (r[i] != r[i])
			{
				printf("r[%d] = INF\n", i);
				//Sleep(10000);
			}
			if (this->F[i] != this->F[i] && this->print_logs)
			{
				printf("this->F[%d] = INF\n", i);

				///Sleep(10000);
			}
		}
	//	CSSD_Matrix<double, double> Precondor;
		//Precondor.PrecondorSSOR(0.75, this);
		Precondor.SolvePrecondorSSOR_SLAE(this->Diag, z, r);
		math::MakeCopyVector_A_into_B(z, p);

		//основные итерации
		double alpha, betta, scal_r_z;
		scal_r_z = math::MakeInnerProduct(r, z);
		for (k = 1; k < maxiter && sqrt(math::MakeInnerProduct(r, r) / math::MakeInnerProduct(this->F, this->F)) > E; k++)
		{
			//temp_mult = A*p
			this->MultiplicationMatrixVector(p, temp_mult);
			//alpha = (r,z)/(A*p,p)
			alpha = scal_r_z / math::MakeInnerProduct(temp_mult, p);
			//x(+1) = x + alpha*p
			//r(+1) = r - alpha*p
			for (int i = 0; i < n_matr; i++)
			{
				this->X[i] += alpha * p[i];
				r[i] -= alpha * temp_mult[i];
			}
			//z(+1) = M^(-1)*r(+1)
			if (!(k % 10) || k == 1)
				Precondor.SolvePrecondorSSOR_SLAE(this->Diag, z, r);
			else
				Precondor.SolvePrecondorSSOR_SLAE(this->Diag, z, r);
			//betta = (r(+1),z(+1))/(r,z)
			betta = math::MakeInnerProduct(r, z) / scal_r_z;
			for (int i = 0; i < n_matr; i++)
			{
				p[i] = z[i] + betta * p[i];
			}

			//готовим на следующий шаг
			scal_r_z = math::MakeInnerProduct(r, z);
			if (this->print_logs && (!(k % 100) || k == 1))
			{
				printf("\tMSG_SSOR\t>\t%d\t-\t%.5e\n", k, sqrt(math::MakeInnerProduct(r, r) / math::MakeInnerProduct(this->F, this->F)));
			}

			//if (sqrt(math::MakeInnerProduct(r, r) / math::MakeInnerProduct(this->F, this->F)) <= E)
			//{
			//	//считаем настоящую невязку
			//	this->MultiplicationMatrixVector(this->X, temp_mult);
			//	double tmp1 = 0, tmp2 = 0;
			//	for (int i = 0; i < this->X.size(); i++)
			//	{
			//		tmp1 += (temp_mult[i] - this->F[i]) * (temp_mult[i] - this->F[i]);
			//		tmp2 += (this->F[i]) * (this->F[i]);
			//	}
			//	if (sqrt(tmp1 / tmp2) <= E) break;
			//	
			//}
		}

		//считаем настоящую невязку
		this->MultiplicationMatrixVector(this->X, temp_mult);
		double tmp1 = 0, tmp2 = 0;
		for (int i = 0; i < this->X.size(); i++)
		{
			tmp1 += (temp_mult[i] - this->F[i]) * (temp_mult[i] - this->F[i]);
			tmp2 += (this->F[i]) * (this->F[i]);
		}
		if(this->print_logs)
		printf("\tMSG_SSOR\t>\t%d\t-\t%.5e (%.5e resolve residual)\n",
			k,
			sqrt(math::MakeInnerProduct(r, r) / math::MakeInnerProduct(this->F, this->F)),
			sqrt(tmp1 / tmp2));

		return sqrt(tmp1 / tmp2);
	}
	double MSG(int maxiter, double E)
	{
		//maxiter *= 10;
		int i, k, n_matr = this->GetMatrixSize();
		double r_scal, a, b_msg;
		std::vector<double>r, z, temp_mult;
		r.resize(n_matr);
		z.resize(n_matr);
		temp_mult.resize(n_matr);

		this->MultiplicationMatrixVector(this->X, temp_mult);
		for (i = 0; i < n_matr; i++)
		{
			r[i] = this->F[i] - temp_mult[i];
			z[i] = r[i];
			if (r[i] != r[i])
			{
				printf("r[%d] = INF\n", i);
				//Sleep(10000);
			}

			if (this->F[i] != this->F[i])
			{
				printf("this->F[%d] = INF\n", i);

				///Sleep(10000);
			}
		}
		for (k = 1; k <= maxiter && sqrt(math::MakeInnerProduct(r, r) / math::MakeInnerProduct(this->F, this->F)) > E; k++)
		{
			r_scal = math::MakeInnerProduct(r, r);
			this->MultiplicationMatrixVector(z, temp_mult);
			a = r_scal / math::MakeInnerProduct(temp_mult, z);
			for (i = 0; i < n_matr; i++)
			{
				this->X[i] += a * z[i];
				r[i] -= a * temp_mult[i];
			}
			b_msg = math::MakeInnerProduct(r, r) / r_scal;
			for (i = 0; i < n_matr; i++)
			{
				z[i] = r[i] + b_msg * z[i];
			}
			/*mult(U);
			summ_v(r, 1, temp_mult, -1, this->F);*/

			//if( k == maxiter ) printf("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA!!!!!!!\n");
			if (!(k % 1000) || k == 1)
			{
				if(this->print_logs) printf("\tMSG\t>\t%d\t-\t%.5e\n", k-1, sqrt(math::MakeInnerProduct(r, r) / math::MakeInnerProduct(this->F, this->F)));
				/*mult(U);
				double tmp1 = 0, tmp2 = 0;
				for (int i = 0; i < U.size(); i++)
				{
					tmp1 += (temp_mult[i] - this->F[i])*(temp_mult[i] - this->F[i]);
					tmp2 += (this->F[i])*(this->F[i]);
				}
				printf("-%d - %.15lf\n", k, sqrt(tmp1 / tmp2));*/
			}
			/*if (k % 100 == 0)
			{
				mult(U, temp_mult);
				summ_v(r, -1, temp_mult, 1, this->F);
			}*/
		}
		//printf(">>%d - %.15e\n", k, sqrt(math::MakeInnerProduct(r,r)/math::MakeInnerProduct(this->F,this->F)));

		this->MultiplicationMatrixVector(this->X, temp_mult);
		double tmp1 = 0, tmp2 = 0;
		for (int i = 0; i < X.size(); i++)
		{
			tmp1 += (temp_mult[i] - this->F[i])*(temp_mult[i] - this->F[i]);
			tmp2 += (this->F[i])*(this->F[i]);
		}
		if (this->print_logs)printf("\tMSG\t>\t%d\t-\t%.5e (resolve residual)\n", k, sqrt(tmp1 / tmp2));
		double current_r = sqrt(tmp1 / tmp2);

		temp_mult.clear();

		return current_r;
	}
	double MSG(int maxiter, double E, std::vector<TForVector>& X_new, std::vector<TForVector>& F_new, bool print_information)
	{
		//maxiter *= 10;
		int i, k, n_matr = this->GetMatrixSize();
		double r_scal, a, b_msg;
		std::vector<double>r, z, temp_mult;
		r.resize(n_matr);
		z.resize(n_matr);
		temp_mult.resize(n_matr);

		this->MultiplicationMatrixVector(X_new, temp_mult);
		for (i = 0; i < n_matr; i++)
		{
			r[i] = F_new[i] - temp_mult[i];
			z[i] = r[i];
			if (r[i] != r[i])
			{
				printf("r[%d] = INF\n", i);
				//Sleep(10000);
			}

			if (F_new[i] != F_new[i])
			{
				printf("F_new[%d] = INF\n", i);

				///Sleep(10000);
			}
		}
		for (k = 1; k <= maxiter && sqrt(math::MakeInnerProduct(r, r) / math::MakeInnerProduct(F_new, F_new)) > E; k++)
		{
			r_scal = math::MakeInnerProduct(r, r);
			this->MultiplicationMatrixVector(z, temp_mult);
			a = r_scal / math::MakeInnerProduct(temp_mult, z);
			for (i = 0; i < n_matr; i++)
			{
				X_new[i] += a * z[i];
				r[i] -= a * temp_mult[i];
			}
			b_msg = math::MakeInnerProduct(r, r) / r_scal;
			for (i = 0; i < n_matr; i++)
			{
				z[i] = r[i] + b_msg * z[i];
			}
			/*mult(U);
			summ_v(r, 1, temp_mult, -1, F_new);*/

			//if( k == maxiter ) printf("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA!!!!!!!\n");
			if ((!(k % 1000) || k == 1) && print_information == true)
			{
			printf(">%d - %.15lf\n", k, sqrt(math::MakeInnerProduct(r, r) / math::MakeInnerProduct(F_new, F_new)));
			/*mult(U);
			double tmp1 = 0, tmp2 = 0;
			for (int i = 0; i < U.size(); i++)
			{
				tmp1 += (temp_mult[i] - F_new[i])*(temp_mult[i] - F_new[i]);
				tmp2 += (F_new[i])*(F_new[i]);
			}
			printf("-%d - %.15lf\n", k, sqrt(tmp1 / tmp2));*/
			}
			/*if (k % 100 == 0)
			{
				mult(U, temp_mult);
				summ_v(r, -1, temp_mult, 1, F_new);
			}*/
		}
		//printf(">>%d - %.15e\n", k, sqrt(math::MakeInnerProduct(r,r)/math::MakeInnerProduct(F_new,F_new)));



		this->MultiplicationMatrixVector(X_new, temp_mult);
		double tmp1 = 0, tmp2 = 0;
		for (int i = 0; i < X_new.size(); i++)
		{
			tmp1 += (temp_mult[i] - F_new[i]) * (temp_mult[i] - F_new[i]);
			tmp2 += (F_new[i]) * (F_new[i]);
		}
		if (print_information == true)
		{
			printf(">>%d - %.15e\n", k, sqrt(tmp1 / tmp2));
		}
		double current_r = sqrt(tmp1 / tmp2);

		temp_mult.clear();


		return current_r;
	}
	//double MSG_LLt(int maxiter, double e)
	//{
	//	int i, k;
	//	double r_scal, a, b_msg;
	//	std::vector<double>r, z, temp_mult2;
	//	al_LU.resize(al.size());
	//	diag_LL.resize(n_matr);
	//	r.resize(n_matr);
	//	z.resize(n_matr);
	//	temp_mult.resize(n_matr);
	//	tempT_mult.resize(n_matr);
	//	temp_mult2.resize(n_matr);
	//	LL();
	//	//for(i=0; i < n_matr;i++)
	//		//U[i] = 1;
	//	mult(U);
	//	for (i = 0; i < n_matr; i++)
	//	{
	//		r[i] = b[i] - temp_mult[i];
	//		temp_mult2[i] = 0;
	//		z[i] = 0;
	//	}
	//	down(r, temp_mult);
	//	middle(temp_mult, temp_mult2);
	//	up(temp_mult2, z);
	//	for (k = 0; k<maxiter && sqrt(scal(r, r) / scal(b, b))>e; k++)
	//	{
	//		down(r, temp_mult2);
	//		middle(temp_mult2, temp_mult);
	//		up(temp_mult, temp_mult2);
	//		mult(z);
	//		r_scal = scal(temp_mult2, r);
	//		a = r_scal / scal(temp_mult, z);
	//		for (i = 0; i < n_matr; i++)
	//		{
	//			U[i] += a * z[i];
	//			r[i] -= a * temp_mult[i];
	//		}
	//		down(r, temp_mult2);
	//		middle(temp_mult2, temp_mult);
	//		up(temp_mult, temp_mult2);
	//		b_msg = scal(temp_mult2, r) / r_scal;
	//		for (i = 0; i < n_matr; i++)
	//			z[i] = temp_mult2[i] + b_msg * z[i];
	//		if (!(k % 1000) || k == 0)
	//			printf(">%d - %.15lf\n", k, sqrt(scal(r, r) / scal(b, b)));
	//	}
	//	mult(U);
	//	//up(X,X,au_LU,diag_U);
	//	//printf(">%d - %.15lf\n", k, sqrt(scal(r,r)/scal(b,b)));
	//	//al_LL.clear();
	//	//diag_LL.clear();
	//	//temp_mult.clear();
	//	//tempT_mult.clear();
	//}

	/// <summary>
		/// процедура минимизации
		/// возвращает норму невязки
		/// </summary>
		/// <param name="Norma_r - текущая норма невязки"></param>
		/// <param name="f - результат минимизации"></param>
	//m-глубина метода
	double Minimization_Problem(int m, double Norma_r, std::vector<double> &f, std::vector<std::vector<double>> &H)
	{
		//вектор правой части для решения СЛАУ с треуг.матрицей f = ||r|| * e1
		f.resize(m + 1);
		for (int i = 1; i < m + 1; i++) 
			f[i] = 0.0;
		f[0] = Norma_r;

		//вспомогательные переменные
		double help1, help2;

		//косинус, синус, тангенс
		double c, s, t;

		//преобразования Гивенса: приведение матрицы Хессенберга H к верхней треугольной форме
		for (int i = 0; i < m; i++)
		{
			t = H[i + 1][i] / H[i][i];
			c = 1.0 / sqrt(1 + t * t);
			s = t * c;

			//H_new = Gt * H (минус у матрицы вращения внизу)
			for (int k = i; k < H[0].size(); k++)
			{
				help1 = c * H[i][k] + s * H[i + 1][k];
				help2 = c * H[i + 1][k] - s * H[i][k];
				H[i][k] = help1;
				H[i + 1][k] = help2;
			}

			//перемножаем слева вектор правой части на трансп.матрицу преобразования
			help1 = c * f[i] + s * f[i + 1];
			help2 = c * f[i + 1] - s * f[i];
			f[i] = help1;
			f[i + 1] = help2;
		}

		//решается система m X m с верхней треугольной матрицей
		for (int i = m - 1; i >= 0; i--)
		{
			f[i] /= H[i][i];
			for (int j = i - 1; j >= 0; j--)
				f[j] -= H[j][i] * f[i];
		}

		return abs(f[m]);
	}
	double GMRES(int maxiter, double E)
	{
		int M = 5; //глубина метода
		int n = this->GetMatrixSize();
		//вспомогательные матрицы метода
		std::vector<std::vector<double>> H, H_, V, W;
		math::ResizeVector(H, M + 1, M);
		math::ResizeVector(H_, M + 1, M);
		math::ResizeVector(V, M + 1, this->GetMatrixSize());
		math::ResizeVector(W, M, this->GetMatrixSize());

		//вспомогательные векторы
		std::vector<double> r(this->GetMatrixSize());
		std::vector<double> Help(this->GetMatrixSize());
		std::vector<double> Minimization_Result(M+1);

		//параметры метода: флаг останова и глубина GMRES
		int flag_m = 0, old_m = M;

		//норма невязкаи
		double Norma_r, Norma_F;

		//деление на корень из квадратов
		for (int i = 0; true && i < this->GetMatrixSize(); i++)
		{
			TForVector koef;
			koef = 0.0;
			koef += Diag[i] * Diag[i];

			for (int j = 0; j < this->A_up[i].size(); j++)
			{
				koef += A_up[i][j] * A_up[i][j];
			}
			for (int j = 0; j < this->A_down[i].size(); j++)
			{
				koef += A_down[i][j] * A_down[i][j];
			}

			koef = math::SolveSqrt(koef);

			X[i] /= koef;
			F[i] /= koef;

			Diag[i] /= koef;

			for (int j = 0; j < this->A_up[i].size(); j++)
			{
				A_up[i][j] /= koef;
			}
			for (int j = 0; j < this->A_down[i].size(); j++)
			{
				A_down[i][j] /= koef;
			}
		}

		//вычисляем вектор невязки
		this->MultiplicationMatrixVector(this->X, Help);
		math::MakeSummVectors(this->F, 1.0, Help, -1.0, r);

		///процедура предобусловливания r = M^(-1) * (F - A * x)
		///Preconditioner.Start_Preconditioner(Help, r);

		//норма вектора невязки
		Norma_r = math::MakeNorma_inH1(r);
		Norma_F = math::MakeNorma_inH1(this->F);

		//итерации метода GMRES
		int Flag = 0, Iter= 0;
		while (Flag == 0 && Iter < maxiter)
		{
			//восстановление параметров
			M = old_m;
			flag_m = 0;

			//ЗДЕСЬ БЫЛО ЗАНУЛЕНИЕ МАТРИЦ H V W
			math::InitializationVector(H, 0.0);
			math::InitializationVector(V, 0.0);
			math::InitializationVector(W, 0.0);
			math::InitializationVector(H_, 0.0);


			//строим первый вектор подпространства Крылова v1 = r / ||r||
			for (int i = 0; i < n; i++)
				V[0][i] = r[i] / Norma_r;

			//строим оставшиеся базисные векторы
			for (int j = 0; j < M && flag_m == 0; j++)
			{
				//Help = A * Vj
				this->MultiplicationMatrixVector(V[j], Help);

				///w = Help = M^(-1) * Help
				///Preconditioner.Start_Preconditioner(Help, Help);

				//процедура ортогонализации Арнольди
				for (int i = 0; i <= j; i++)
				{
					//скалярное произведение Hij = Wj * Vi, где Wj = Help
					H[i][j] = 0.0;
					for (int l = 0; l < n; l++)
					{
						W[j][l] = Help[l];
						H[i][j] += W[j][l] * V[i][l];
					}

					//в H_ хранится не модифицированная матрица коэффициентов ортогонализации
					H_[i][j] = H[i][j];

					//Wj = Wj - Hij * Vi
					for (int l = 0; l < n; l++)
					{
						W[j][l] -= H[i][j] * V[i][l];
					}
				}

				//расширение матрицы H и H_
				H[j + 1][j] = math::MakeNorma_inH1(W[j]);
				H_[j + 1][j] = H[j + 1][j];

				if (abs(H[j + 1][j]) < E)
				{
					M = j;
					flag_m = 1;
				}
				else
				{
					for (int l = 0; l < n; l++) 
						V[j + 1][l] = W[j][l] / H[j + 1][j];
				}
			}

			/*
			//проверка на ортогональность
			for (int I = 0; I < V.M; I++)
			{
				for (int J = I; J < V.M; J++)
				{
					Console.WriteLine("V{0} * V{1} = {2}", I+1, J+1, Operation.ScalMult(V.Elem[I], V.Elem[J]));
					Console.ReadLine();
				}
			}
			*/

			//решение задачи минимизации (МНК) и норма невязки
			double Norma_Residual = Minimization_Problem(M, Norma_r, Minimization_Result, H);
			if (Iter % 1 == 0)
			{
				std::vector<double> tmp(Help.size());
				this->MultiplicationMatrixVector(this->X, Help);
				math::MakeSummVectors(this->F, 1.0, Help, -1.0, tmp);
				double real_r = math::MakeNorma_inH1(tmp);
				printf_s(">%d -  %.5e \t(real: %.5e)\n", Iter, Norma_Residual, real_r / Norma_F);
			}

			if (Norma_Residual < E) Flag = 1;
			else
			{
				//вычисляем произведение Vy = V(t) * y
				math::MultiplicationTMatrixVector(V, Minimization_Result, Help);

				//вычисление нового результата
				for (int l = 0; l < this->GetMatrixSize(); l++) 
				{ 
					X[l] += Help[l]; 
				}

				/*
				//вычисление вектора невязки r = V * (||r||e1 - (H_)m * y)
				for (int i = 0; i < m; i++)
				{
					r.Elem[i] = 0.0;
					for (int j = 0; j < m; j++) r.Elem[i] -= H_.Elem[i][j] * Minimization_Result.Elem[j];
				}
				r.Elem[0] += Norma_r;

				for (int i = 0; i < n; i++)
				{
					Help.Elem[i] = 0.0;
					for (int j = 0; j < m; j++)
					{
						Help.Elem[i] += V.Elem[j][i] * r.Elem[j];
					}
				}
				*/

				this->MultiplicationMatrixVector(X, Help);
				//math::MakeSummVectors(this->F, 1.0, Help, -1.0, Help); ///!!
				math::MakeSummVectors(this->F, 1.0, Help, -1.0, r); ///!!


				//процедура предобусловливания r = M^(-1) * Help, где Help = f - A * x
				//Preconditioner.Start_Preconditioner(Help, r);

				Norma_r = math::MakeNorma_inH1(r);

				Iter++;
			}
		}

		this->MultiplicationMatrixVector(this->X, Help);
		double temp = 0;
		for (int i = 0; i < this->GetMatrixSize(); i++)
			temp += (F[i] - Help[i])*(F[i] - Help[i]);
		double result = sqrt(temp) / Norma_F;
		printf(">REAL RESULT> %d - %.15e\n", Iter, result);

		H.shrink_to_fit();
		H_.shrink_to_fit();
		V.shrink_to_fit();
		W.shrink_to_fit();
		r.shrink_to_fit();
		Help.shrink_to_fit();
		Minimization_Result.shrink_to_fit();

		return result;
	}
};