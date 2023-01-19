#pragma once
#include "stdio.h"
#include "Point.h"
#include <vector>
#include "Math.h"

#define mu0 4*M_PI*1E-7
#define eps0 8.85*1E-12

class Domain
{
private:
	

public:
	struct forMechanics
	{
	private:
		double K_Ic; //critical stress intensity factor for I mode
		double pressure_of_fluid; //for hidrofracture
		double epsilon; //Young's modulus
		std::vector<std::vector<double>> Epsilon;
		
		std::vector<std::vector<double>> D; //elastisity tensor
	public:
		double v; //Poisson's ratio
		std::vector<double> alpha_termal;

		void resize_D()
				{
					D.resize(6);
					for (int i = 0; i < 6; i++)
					{
						D[i].resize(6);
					}
				}
		void resize_D_2D()
		{
			D.resize(3);
			for (int i = 0; i < 3; i++)
			{
				D[i].resize(3);
			}
		}
		void resize_D_1D()
		{
			D.resize(1);
			for (int i = 0; i < 1; i++)
			{
				D[i].resize(1);
			}
		}
		void SetIsotropic_D()
		{
			double k = epsilon*(1. - v) / (1. + v) / (1. - 2. * v);
			double b = (1. - 2. * v) / (2. * (1. - v));
			double a = v / (1. - v);

			resize_D();
			for (int i = 0; i < 3; i++)
			{
				for (int j = 0; j < 3; j++)
				{
					if (i == j)
					{
						D[i][j] = 1.;
					}
					else {
						D[i][j] = a;
					}
					D[i][j] *= k;
				}
			}
			for (int i = 3; i < 6; i++)
			{
				D[i][i] = b*k;
			}
		}
		void SetIsotropic_D_2D()
		{
			double k = epsilon / (1.0 - v*v);
			double b = (1 - v) / (2.0);
			double a = v;

			resize_D_2D();
			for (int i = 0; i < 2; i++)
			{
				for (int j = 0; j < 2; j++)
				{
					if (i == j)
					{
						D[i][j] = 1;
					}
					else {
						D[i][j] = a;
					}
					D[i][j] *= k;
				}
			}
			for (int i = 2; i < 3; i++)
			{
				D[i][i] = b*k;
			}
		}
		void SetIsotropic_D_1D()
		{
			double k = epsilon*(1 - v) / (1 - v*v);
			double b = (1 - 2 * v) / (2 * (1 - v));
			double a = v / (1 - v);

			resize_D_1D();
			D[0][0] = k;			
		}
		void SetIsotropic_D(int Dim)
		{
			switch (Dim)
			{
			case 1: SetIsotropic_D_1D(); break;
			case 2: SetIsotropic_D_2D(); break;
			case 3: SetIsotropic_D(); break;
			default:
				break;
			}
			return;
		}
		void SetNonlinear_E(std::vector<std::vector<double>> &E)
		{
			Epsilon.resize(2);
			Epsilon[0].resize(E[0].size());
			Epsilon[1].resize(E[1].size());
			for (int i = 0; i < E[0].size(); i++)
			{
				Epsilon[0][i] = E[0][i];
				Epsilon[1][i] = E[1][i];
			}
		}
		void GetD(double temperature, int dimension, std::vector<std::vector<double>> &res_D)
		{
			if (Epsilon.size() != 0)
			{
				bool find = false;
				if (temperature <= Epsilon[0][0])
				{
					this->epsilon = Epsilon[1][0];
				}
				for (int i = 0; !find && i < Epsilon[0].size()-1; i++)
				{
					if (Epsilon[0][i] <= temperature && temperature <= Epsilon[0][i + 1])
					{
						double coef = (temperature - Epsilon[0][i]) / (Epsilon[0][i + 1] - Epsilon[0][i]);
						this->epsilon = Epsilon[1][i] + (Epsilon[1][i + 1] - Epsilon[1][i]) * coef;
						find = true;
						break;
					}
				}
				if (!find && temperature >= Epsilon[0][Epsilon[0].size()-1])
				{
					this->epsilon = Epsilon[1][Epsilon[0].size() - 1];
					find = true;
				}
				SetIsotropic_D(dimension);
			}
			else {
				if (D.size() == 0)
					SetIsotropic_D(dimension);
			}

			res_D = this->D;
		};
		std::vector<std::vector<double>> GetD(int dimension)
		{
			if (D.size() == 0)
				SetIsotropic_D(dimension);
			return this->D;
		}
		double GetE(double t)
		{
			/*if (Epsilon.size() != 0 )
			{
				bool find = false;
				if (t <= Epsilon[0][0])
				{
					this->epsilon = Epsilon[1][0];
				}
				for (int i = 0; !find && i < Epsilon[0].size() - 1; i++)
				{
					if (Epsilon[0][i] <= t && t <= Epsilon[0][i + 1])
					{
						double coef = (t - Epsilon[0][i]) / (Epsilon[0][i + 1] - Epsilon[0][i]);
						this->epsilon = Epsilon[1][i] + (Epsilon[1][i + 1] - Epsilon[1][i]) * coef;
						find = true;
						break;
					}
				}
				if (!find && t >= Epsilon[0][Epsilon[0].size() - 1])
				{
					this->epsilon = Epsilon[1][Epsilon[0].size() - 1];
					find = true;
				}
			}*/
			return this->epsilon;
		}
		void SetE(double eps) { this->epsilon = eps; };
		void SetV(double v) { this->v = v; };
		void SetK_Ic(double K_Ic) { this->K_Ic = K_Ic; };
		void SetPressureOfFluid(double pressure_of_fluid) { this->pressure_of_fluid = pressure_of_fluid; };
		void SetThermal_Alpha(double al) 
		{
			this->alpha_termal.resize(6);
			alpha_termal[0] = al;
			alpha_termal[1] = al;
			alpha_termal[2] = al;
		};
		double GetLongitudinalWaveVelocity_Vp(double rpho)
		{
			double lambda = this->epsilon * this->v / ((1 + this->v) * (1 - 2 * this->v));
			double mu = this->epsilon / (2 * (1 + this->v));
			double velocity = sqrt((lambda + 2 * mu) / rpho);
			return velocity;
		}

		void GetLameCoefficients(double& lambda, double& mu)
		{
			lambda = epsilon * v / ((1 + v) * (1 - 2 * v));
			mu = epsilon / (2 * (1 + v));
		}

		double GetK_Ic() { return this->K_Ic; };
		double GetPressureOfFluid() { return this->pressure_of_fluid; };

		std::vector<double> transfer_tensor_into_vector_EPS(Tensor2Rank3D &T)
		{
			std::vector<double> vT(6);
			vT[0] = T.val[0][0];
			vT[1] = T.val[1][1];
			vT[2] = T.val[2][2];
			vT[3] = T.val[1][2]*2;
			vT[4] = T.val[0][2]*2;
			vT[5] = T.val[0][1]*2;
			return vT;
		}
		Tensor2Rank3D transfer_vector_into_tensor_EPS(std::vector<double> &vT)
		{
			Tensor2Rank3D T;
			T.val[0][0] = vT[0];
			T.val[1][1] = vT[1];
			T.val[2][2] = vT[2];
			T.val[1][2] = vT[3]/2.;
			T.val[2][1] = vT[3]/2.;
			T.val[0][2] = vT[4]/2.;
			T.val[2][0] = vT[4]/2.;
			T.val[0][1] = vT[5]/2.;
			T.val[1][0] = vT[5]/2.;
			return T;
		}
		std::vector<double> transfer_tensor_into_vector_SIG(Tensor2Rank3D &T)
		{
			std::vector<double> vT(6);
			vT[0] = T.val[0][0];
			vT[1] = T.val[1][1];
			vT[2] = T.val[2][2];
			vT[3] = T.val[1][2];
			vT[4] = T.val[0][2];
			vT[5] = T.val[0][1];
			return vT;
		}
		Tensor2Rank3D transfer_vector_into_tensor_SIG(std::vector<double> &vT)
		{
			Tensor2Rank3D T;
			T.val[0][0] = vT[0];
			T.val[1][1] = vT[1];
			T.val[2][2] = vT[2];
			T.val[1][2] = vT[3];
			T.val[2][1] = vT[3];
			T.val[0][2] = vT[4];
			T.val[2][0] = vT[4];
			T.val[0][1] = vT[5];
			T.val[1][0] = vT[5];
			return T;
		}
		Tensor2Rank3D solve_SIGMA(Tensor2Rank3D &EPS)
		{
			std::vector<double> vSIG;
			auto vEPS = transfer_tensor_into_vector_EPS(EPS);
			if (this->D.size() == 0)
				this->SetIsotropic_D();
			math::MultiplicationMatrixVector(this->D, vEPS, vSIG);

			return transfer_vector_into_tensor_SIG(vSIG);
		}
		std::vector<double> solve_vSIGMA(Tensor2Rank3D &EPS)
		{
			std::vector<double> vSIG;
			auto vEPS = transfer_tensor_into_vector_EPS(EPS);
			math::MultiplicationMatrixVector(this->D, vEPS, vSIG);

			return vSIG;
		}
		double mult_SIGMA_EPS(Tensor2Rank3D &EPS, Tensor2Rank3D &SIG)
		{
			auto vEPS = transfer_tensor_into_vector_EPS(EPS);
			auto vSIG = transfer_tensor_into_vector_SIG(SIG);
			return math::MakeInnerProduct(vEPS, vSIG);
		}
		void print_D(FILE *f)
		{
			for (int i = 0; i < D.size(); i++)
			{
				for (int j = 0; j < D[0].size(); j++)
				{
					if (D[i][j] > 0)
					{
						fprintf_s(f, " ");
					}
					fprintf_s(f, " %.2e ", D[i][j]);
				}
				fprintf_s(f, "\n");
			}
		}
		void scanf_D(FILE *f)
		{
			resize_D();
			for (int i = 0; i < D.size(); i++)
			{
				for (int j = 0; j < D[0].size(); j++)
				{
					fscanf_s(f, "%lf", &D[i][j]);
				}
			}
		}
	} forMech;
	struct forAcoustics
	{
		double k; //wave number
	} forAcoust;
	struct forThermals
	{
		double lambda;
		double c;
		double rpho;
	} forThermal;
	struct forElectrical {
		double sigma;
		double mu;
		double eps;
		double w;
		Tensor2Rank3D sigma_tensor;

		void SetSigmaTensor()
		{
			sigma_tensor.InitializationAs0();
			sigma_tensor.val[0][0] = sigma;
			sigma_tensor.val[1][1] = sigma;
			sigma_tensor.val[2][2] = sigma;
		}
		void SetSigmaTensor(double s_xx, double s_yy, double s_zz)
		{
			sigma_tensor.InitializationAs0();
			sigma_tensor.val[0][0] = s_xx;
			sigma_tensor.val[1][1] = s_yy;
			sigma_tensor.val[2][2] = s_zz;
		}
		void SetSigmaTensor(std::vector<std::vector<double>> sigma)
		{
			sigma_tensor.InitializationAs0();
			for (int i = 0; i < 3; i++)
			{
				for (int j = 0; j < 3; j++)
				{
					sigma_tensor.val[i][j] = sigma[i][j];
				}
			}
		}
	} forElectrical;
	
	Domain()
	{};
	~Domain()
	{

	};
};