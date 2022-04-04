#pragma once
#include <vector>
#include "Point.h"
#include <cstdlib>

class Tensor2Rank3D{
private:
	bool IsEqual(double A, double B)
	{
		if (abs(A - B) < 1e-10) return true;
		return false;
	}
	double GetDeterminant(const Tensor2Rank3D& A)
	{
		return A.val[0][0] * A.val[1][1] * A.val[2][2] + A.val[0][1] * A.val[1][2] * A.val[2][0] + A.val[1][0] * A.val[0][2] * A.val[2][1]
			- A.val[0][2] * A.val[1][1] * A.val[2][0] - A.val[0][0] * A.val[2][1] * A.val[1][2] - A.val[0][1] * A.val[1][0] * A.val[2][2];
	}
	void MultMatrixMatrix(std::vector<std::vector<double>>& A, std::vector<std::vector<double>>& B, std::vector<std::vector<double>>& result)
	{
		try
		{
			result.resize(A.size());
			for (int i = 0; i < result.size(); i++)
			{
				result[i].resize(B[0].size());
				for (int j = 0; j < result[i].size(); j++)
				{
					double tmp = 0;
					for (int k = 0; k < A[0].size(); k++)
					{
						tmp += A[i][k] * B[k][j];
					}
					result[i][j] = tmp;
				}
			}
		}
		catch (const std::exception&)
		{
			printf_s("Error: Library/Math.h/void MultMatrixMatrix(std::vector<std::vector<T>> &A, std::vector<std::vector<T>> &B, std::vector<std::vector<T>> &result)\n");
			if (A[0].size() != B.size())
			{
				printf("\t\tError in mult matrix (%d,%d)x(%d,%d)\n", (int)A.size(), (int)A[0].size(), (int)B.size(), (int)B[0].size());
			}
		}
	}

public:
	std::vector<std::vector<double>> val;
	
	Tensor2Rank3D()
	{
		this->val.resize(3);
		for (int i = 0; i < 3; i++)
		{
			val[i].resize(3);
		}
	}
	Tensor2Rank3D& operator = ( const Tensor2Rank3D& A )
	{
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				val[i][j] = A.val[i][j];
			}
		}
		return *this;
	}
	Tensor2Rank3D& operator = (const double& A)
	{
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				val[i][j] = A;
			}
		}
		return *this;
	}
	Tensor2Rank3D& operator = (const std::vector<std::vector<double>>& A)
	{
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				val[i][j] = A[i][j];
			}
		}
		return *this;
	}
	Tensor2Rank3D& operator /= (const double& A)
	{
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				val[i][j] = val[i][j] / A;
			}
		}
		return *this;
	}
	Tensor2Rank3D& operator /= ( const Point<double>& A )
	{
		for (int j = 0; j < 3; j++)
		{
			val[0][j] = val[0][j] / A.x;
			val[1][j] = val[1][j] / A.y;
			val[2][j] = val[2][j] / A.z;
		}
		return *this;
	}
	Tensor2Rank3D& operator *= ( const double& A )
	{
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				val[i][j] = val[i][j] * A;
			}
		}
		return *this;
	}
	Tensor2Rank3D& operator += ( const Tensor2Rank3D& A )
	{
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				val[i][j] = val[i][j] + A.val[i][j];
			}
		}
		return *this;
	}
	Tensor2Rank3D& operator -= (const Tensor2Rank3D& A)
	{
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				val[i][j] = val[i][j] - A.val[i][j];
			}
		}
		return *this;
	}
	Tensor2Rank3D operator + (const Tensor2Rank3D& right)
	{
		Tensor2Rank3D result;
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				result.val[i][j] = this->val[i][j] + right.val[i][j];
			}
		}
		return result;
	}
	Tensor2Rank3D operator - (const Tensor2Rank3D& right)
	{
		Tensor2Rank3D result;
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				result.val[i][j] = this->val[i][j] - right.val[i][j];
			}
		}
		return result;
	}
	Tensor2Rank3D operator / (const double& right)
	{
		Tensor2Rank3D result;
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				result.val[i][j] = this->val[i][j] / right;
			}
		}
		return result;
	}
	//A/B = A*B^(-1)
	Tensor2Rank3D operator / (const Tensor2Rank3D& right)
	{
		Tensor2Rank3D result, B;
		double det_B = GetDeterminant(right);
		B.val[0][0] = right.val[1][1] * right.val[2][2] - right.val[1][2] * right.val[2][1];
		B.val[0][0] = - right.val[0][1] * right.val[2][2] + right.val[0][2] * right.val[2][1];
		B.val[0][1] = right.val[0][1] * right.val[1][2] - right.val[0][2] * right.val[1][1];
		B.val[0][0] = -right.val[1][0] * right.val[2][2] + right.val[1][2] * right.val[2][0];
		B.val[0][0] = right.val[0][0] * right.val[2][2] - right.val[0][2] * right.val[2][0];
		B.val[0][1] = -right.val[0][0] * right.val[1][2] + right.val[0][2] * right.val[1][0];
		B.val[1][0] = right.val[1][0] * right.val[2][1] - right.val[1][1] * right.val[2][0];
		B.val[1][0] = -right.val[0][0] * right.val[2][1] + right.val[0][1] * right.val[2][0];
		B.val[1][1] = right.val[0][0] * right.val[1][1] - right.val[0][1] * right.val[1][0];
		
		MultMatrixMatrix(this->val, B.val, result.val);
		return result;
	}
	Tensor2Rank3D operator * (const double& right)
	{
		Tensor2Rank3D result;
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				result.val[i][j] = this->val[i][j] * right;
			}
		}
		return result;
	}
	Point<double> operator * (const Point<double>& right)
	{
		Point<double> result;
		result.x += this->val[0][0] * right.x;
		result.x += this->val[0][1] * right.y;
		result.x += this->val[0][2] * right.z;

		result.y += this->val[1][0] * right.x;
		result.y += this->val[1][1] * right.y;
		result.y += this->val[1][2] * right.z;

		result.z += this->val[2][0] * right.x;
		result.z += this->val[2][1] * right.y;
		result.z += this->val[2][2] * right.z;
		return result;
	}
	//свертка
	Point<double> operator * (const Tensor2Rank3D& right)
	{
		Point<double> result;

		result.x += val[0][0] * val[0][0];
		result.x += val[0][1] * val[0][1];
		result.x += val[0][2] * val[0][2];
		
		result.y += val[1][0] * val[1][0];
		result.y += val[1][1] * val[1][1];
		result.y += val[1][2] * val[1][2];
		
		result.z += val[2][0] * val[2][0];
		result.z += val[2][1] * val[2][1];
		result.z += val[2][2] * val[2][2];

		return result;
	}

	bool operator == ( const Tensor2Rank3D &A )
	{
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				if (!IsEqual(this->val[i][j], A.val[i][j]))
					return false;
			}
		}
		return true;
	}
	bool operator != (const Tensor2Rank3D &A)
	{
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				if (!IsEqual(this->val[i][j], A.val[i][j]))
					return true;
			}
		}
		return false;
	}
	void print()
	{
		for (int i = 0; i < 3; i++)
		{
			printf_s("| ");
			for (int j = 0; j < 3; j++)
			{
				if(this->val[i][j] > 0)
					printf_s(" %.2e ", this->val[i][j]);
				else
					printf_s("%.2e ", this->val[i][j]);
			}
			printf_s("|");
		}
		
	}
	void InitializationAs0()
	{
		for (int i = 0; i < val.size(); i++)
		{
			for (int j = 0; j < val[i].size(); j++)
			{
				val[i][j] = 0.0;
			}
		}
	}
	void InitializationAsI()
	{
		InitializationAs0();
		for (int i = 0; i < val.size(); i++)
		{
			val[i][i] = 1.0;
		}
	}
	void SetValue(double v)
	{
		for (int i = 0; i < val.size(); i++)
		{
			for (int j = 0; j < val[i].size(); j++)
			{
				val[i][j] = v;
			}
		}
	}
	
	Tensor2Rank3D T()
	{
		Tensor2Rank3D res;
		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				res.val[i][j] = this->val[j][i];
			}
		}
		return res;
	}
};