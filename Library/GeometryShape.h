#pragma once
#include <vector>
#include "Point.h"
#include "Math.h"


namespace geometry
{
	class Plane {
	public:
		Point<double> A, B, C;
		Plane()
		{}
		Plane(Point<double> A, Point<double> B, Point<double> C)
		{
			this->A = A;
			this->B = B;
			this->C = C;
		}
		bool IsContain(Point<double> X)
		{
			Point<double> tmp;
			tmp.x = math::GetRound(X.x, 6);
			tmp.y = math::GetRound(X.y, 6);
			tmp.z = math::GetRound(X.z, 6);
			double matr[3][3] = { { tmp.x - A.x, tmp.y - A.y, tmp.z - A.z },
								  { B.x - A.x, B.y - A.y, B.z - A.z },
								  { C.x - A.x, C.y - A.y, C.z - A.z } };
			double d = abs(math::GetDeterminantForMatrix3x3(matr));
			if (d <= 1E-7)
				return true;
			return false;
		}
		double SolveDistanceToPoint(Point<double> X)
		{
			double a = (B.y - A.y)*(C.z - A.z) - (C.y - A.y)*(B.z - A.z);
			double b = (B.z - A.z)*(C.x - A.x) - (C.z - A.z)*(B.x - A.x);
			double c = (B.x - A.x)*(C.y - A.y) - (C.x - A.x)*(B.y - A.y);
			double d = -A.x*(B.y - A.y)*(C.z - A.z) - A.y*(C.x - A.x)*(B.z - A.z) - A.z*(B.x - A.x)*(C.y - A.y)
				+ A.x*(C.y - A.y)*(B.z - A.z) + A.y*(B.x - A.x)*(C.z - A.z) + A.z*(C.x - A.x)*(B.y - A.y);

			return abs(a*X.x + b * X.y + c * X.z + d) / sqrt(a*a + b * b + c * c);
		}
		void GetCoefficients(double &a, double &b, double &c, double &d)
		{
			a = (B.y - A.y)*(C.z - A.z) - (C.y - A.y)*(B.z - A.z);
			b = (B.z - A.z)*(C.x - A.x) - (C.z - A.z)*(B.x - A.x);
			c = (B.x - A.x)*(C.y - A.y) - (C.x - A.x)*(B.y - A.y);
			d = -A.x*(B.y - A.y)*(C.z - A.z) - A.y*(C.x - A.x)*(B.z - A.z) - A.z*(B.x - A.x)*(C.y - A.y)
				+ A.x*(C.y - A.y)*(B.z - A.z) + A.y*(B.x - A.x)*(C.z - A.z) + A.z*(C.x - A.x)*(B.y - A.y);
		}
	};

	class Shape
	{
	public:
		Shape()
		{
			return;
		}

		Shape(std::vector<int> &id_nodes)
		{
			try
			{
				this->nodes.resize(id_nodes.size());
				math::MakeCopyVector_A_into_B(id_nodes, this->id_nodes);
			}
			catch (const std::exception&)
			{
				printf_s("Error: GeometryShape.cpp/geometry::Shape(std::vector<int> &id_nodes) \n");
			}
		};
		Shape(std::vector<int> &id_nodes, std::vector<Point<double>*> &xyz_nodes)
		{
			try
			{
				math::MakeCopyVector_A_into_B(id_nodes, this->id_nodes);
				math::MakeCopyVector_A_into_B(xyz_nodes, this->nodes);
			}
			catch (const std::exception&)
			{
				printf_s("Error: GeometryShape.cpp/geometry::Shape(std::vector<int> &id_nodes, std::vector<Point<double>*> &xyz_nodes)\n");
			}
		};
		~Shape()
		{
			std::vector<Point<double>*> v_n;
			std::vector<Point<double>*>(v_n).swap(this->nodes);
			std::vector<int> v_idn;
			std::vector<int>(v_idn).swap(this->id_nodes);
		};

		Point<double> GetNode(int _id)
		{
			try
			{
				return *this->nodes[_id];
			}
			catch (const std::exception&)
			{
				printf_s("Error: GeometryShape.cpp/geometry::Shape::GetNode(int _id)\n");
			}
		}
		Point<double>* GetPtrNode(int _id)
		{
			return this->nodes[_id];
		}
		void SetNode(int local_id, Point<double> *X)
		{
			try
			{
				this->nodes[local_id] = X;
			}
			catch (const std::exception&)
			{
				printf_s("Error: GeometryShape.cpp/void geometry::Shape::SetNode(int local_id, Point<double> *X)");
			}
		}
		void UpdateNodes(std::vector<Point<double>*>& P)
		{
			try
			{
				this->nodes.resize(P.size());
				for (int i = 0; i < P.size(); i++)
				{
					this->nodes[i] = P[i];
				}
			}
			catch (const std::exception&)
			{
				printf_s("Error: GeometryShape.cpp/Tetrahedron(std::vector<Point<double> *> &P)");
			}
		}

		Point<double> GetWeightCentr()
		{
			try
			{
				Point<double> result;
				for(int i  =0 ; i < this->GetNodesCount(); i++)
				{
					result += *this->nodes[i];
				}
				result /= 1.0* (int)this->nodes.size();
				return result;
			}
			catch (const std::exception&)
			{
				printf_s("Error: GeometryShape.cpp/geometry::Shape::GetWeightCentr()\n");
			}
		}
		int GetNodesCount()
		{
			return (int)this->nodes.size();
		}
		int GetIdNode(int local_id)
		{
			try
			{
				return this->id_nodes[local_id];
			}
			catch (const std::exception&)
			{
				printf_s("Error: GeometryShape.cpp/geometry::Shape::SetIdNode(int local_id, int global_id)\n");
			}
		};
		std::vector<Point<double>> GetNodes()
		{
			try
			{
				std::vector<Point<double>> res(this->GetNodesCount());
				for (int i = 0; i < res.size(); i++)
				{
					res[i] = *this->nodes[i];
				}
				return res;
			}
			catch (const std::exception&)
			{
				printf_s("Error: GeometryShape.cpp/geometry::Shape::SetIdNode(int local_id, int global_id)\n");
			}
		};
		std::vector<int>* GetIdNodes()
		{
			return &(this->id_nodes);
		}
		int GetLocalID_forVertexID(int global_id)
		{
			int result = -1;
			for (int id = 0; id < this->GetNodesCount(); id++)
			{
				if (this->GetIdNode(id) == global_id)
				{
					result = id;
					break;
				}
			}
			return result;
		}
		void SetIdNodes(std::vector<int> &id_nodes)
		{
			try
			{
				this->id_nodes.resize(id_nodes.size());
				for (int i = 0; i < id_nodes.size(); i++)
				{
					this->id_nodes[i] = id_nodes[i];
				}

				if (this->nodes.size() == 0)
					this->nodes.resize(id_nodes.size());
			}
			catch (const std::exception&)
			{
				printf_s("Error: GeometryShape.cpp/geometry::Shape::SetIdNodes(std::vector<int> &id_nodes)\n");
			}
		}
		void SetIdNode(int local_id, int global_id)
		{
			try
			{
				for (int i = (int)this->id_nodes.size(); i <= local_id; i++)
				{
					this->id_nodes.push_back(-1);
				}
				this->id_nodes[local_id] = global_id;
			}
			catch (const std::exception&)
			{
				printf_s("Error: GeometryShape.cpp/geometry::Shape::SetIdNode(int local_id, int global_id)\n");
			}
		}
		
		void SetIdDomain(int id)
		{
			this->id_domain = id;
		}
		int GetIdDomain()
		{
			return this->id_domain;
		}

		virtual double GetVolume() = 0;
		virtual bool IsContainThePoint(Point<double> A) = 0;
		virtual bool IsContainThePoint(Point<double> A, double &length) = 0;
		bool IsElemInBox(Point<double> down, Point <double> up)
		{
			for (int i = 0; i < nodes.size(); i++)
			{
				if (!(down < *nodes[i] && *nodes[i] < up))
					return false;
			}
			return true;
		}
		bool IsPartOfElemInBox(Point<double> bottom, Point <double> top)
		{
			Point<double> down, up;
			down = *nodes[0];
			up = *nodes[0];
			for (int i = 1; i < nodes.size(); i++)
			{
				if ((*nodes[i]).x < down.x) down.x = (*nodes[i]).x;
				if ((*nodes[i]).y < down.y) down.y = (*nodes[i]).y;
				if ((*nodes[i]).z < down.z) down.z = (*nodes[i]).z;

				if ((*nodes[i]).x > up.x) up.x = (*nodes[i]).x;
				if ((*nodes[i]).y > up.y) up.y = (*nodes[i]).y;
				if ((*nodes[i]).z > up.z) up.z = (*nodes[i]).z;
			}


			auto qwe = [](Point<double> &x, Point<double> y)->bool {
				if (x.x >= y.x && x.y >= y.y && x.z >= y.z)
					return true;
				return false;
			};

			if ((qwe(down, bottom) || qwe(up, bottom)) && (qwe(top, down) || qwe(top, up)))
				return true;
			return false;
		}

		virtual void SetIntegrationLaw(int points_count) = 0;
		double GetIntegrationPointsCount() {
			return this->integration_parameters.GetPointNumber();
		}
		double GetIntegrationWeight(int id)
		{
			return this->integration_parameters.w[id];
		}
		Point<double> GetIntegrationPoint(int id)
		{
			return this->integration_parameters.x[id];
		}

		template <typename T>
		T SolveIntegral(std::function< T(Point<double>) > &function)
		{
			T result;
			result = 0.0;
			if (this->integration_parameters.GetPointNumber() != 0)
			{
				for (int i = 0; i < this->integration_parameters.GetPointNumber(); i++)
				{
					result += function(this->integration_parameters.x[i]) * this->integration_parameters.w[i];
				}
			}
			else {
				this->SetIntegrationLaw(-1);
			}
			return result;
		}

	protected:
		struct Integral {
			std::vector<double> w;
			std::vector<Point<double>> x;
			int GetPointNumber() { return int(x.size()); }
		} integration_parameters;
		std::vector<Point<double>*> nodes;
		std::vector<int> id_nodes; //in global grid
		int id_domain;
	};
	class Vertex : public Shape
	{
	public:
		Vertex()
		{
			this->id_nodes.resize(1);
			this->nodes.resize(1);
		};
		Vertex(std::vector<Point<double> *> &P)
		{
			try
			{
				this->UpdateNodes(P);
			}
			catch (const std::exception&)
			{
				printf_s("Error: GeometryShape.cpp/Triangle(std::vector<Point<double> *> &P)");
			}
		}
		double GetVolume()
		{
			try
			{
				return 0;
			}
			catch (const std::exception&)
			{
				printf_s("Error: GeometryShape.cpp/GetVolume()");
			}
		}

		void SetGeometry(std::vector<int>& id_nodes, std::vector<Point<double>*>& P)
		{
			try
			{
				this->SetIdNodes(id_nodes);
				this->UpdateNodes(P);
			}
			catch (const std::exception&)
			{
				printf_s("Error: GeometryShape.h/SetGeometry(std::vector<std::vector<int>>& id_nodes, std::vector<Point<double>*>& P)");
			}
		};


		void SetIntegrationLaw(int points_count)
		{
		}

		bool IsContainThePoint(Point<double> A)
		{
			try
			{
				return false;
			}
			catch (const std::exception&)
			{
				printf_s("Error: GeometryShape.cpp/IsContainThePoint(Point<double> A)");
			}
		}
		bool IsContainThePoint(Point<double> A, double &length)
		{
			try
			{
				length = math::SolveLengthVector(A, *this->nodes[0]);
				return false;
			}
			catch (const std::exception&)
			{
				printf_s("Error: GeometryShape.cpp/IsContainThePoint(Point<double> A)");
			}
		}
	};
	class Segment : public Shape
	{
		std::vector<std::vector<double>> self_basis;
		std::vector<std::vector<double>> self_basis_revers;
		std::vector<Point<double>> self_nodes;

		void CreateSelfBasis()
		{
			double self_basis[3][3];
			for (int i = 0; i < 3; i++)
				for (int j = 0; j < 3; j++)
					self_basis[i][j] = 0.0;

			Point<double> X, Y, Z; //��� ���� ����� ���� (�� ��������� ������ ��������� � (0,0,0))

			//������ ���� X
			X = *this->nodes[1] - *this->nodes[0];
			X /= math::SolveLengthVector(X, Point<double>(0, 0, 0)); //����������

			double t = math::SolveLengthVector(X, Point<double>(0, 0, 0));

			//������ ���� Y �� ����������� ��������������� � �
			if (X.z != 0)
			{
				Y.x = 1.0;
				Y.y = 1.0;
				Y.z = (-X.x - X.y) / X.z;
			}
			else {
				if (X.y != 0)
				{
					Y.x = 1.0;
					Y.y = (-X.x - X.z) / X.y;
					Y.z = 1.0;
				}
				else {
					Y.x = (-X.y - X.z) / X.x;
					Y.y = 1.0;
					Y.z = 1.0;
				}
			}
			Y /= math::SolveLengthVector(Y, Point<double>(0, 0, 0)); //����������

			t = math::SolveLengthVector(Y, Point<double>(0, 0, 0));

			//������ ���� Z �� ����������� ��������������� � � � Y
			if (X.z != 0 && (Y.y - X.y * Y.z / X.z) != 0)
			{
				Z.x = 1.0;
				Z.y = (X.x * Y.z / X.z - Y.x) * Z.x / (Y.y - X.y * Y.z / X.z);
				Z.z = (-X.x * Z.x - X.y * Z.y) / X.z;
			}
			else {
				if (X.y != 0 && (Y.z - X.z * Y.y / X.y) != 0)
				{
					Z.x = 1.0;
					Z.z = (X.x * Y.y / X.y - Y.x) * Z.x / (Y.z - X.z * Y.y / X.y);
					Z.y = (-X.x * Z.x - X.z * Z.z) / X.y;
				}
				else {
					Z.z = 1.0;
					Z.y = (X.z * Y.x / X.x - Y.z) / (Y.y - X.y * Y.x / X.x);
					Z.x = (-X.z - X.y * Z.y) / X.x;
				}
			}
			Z /= math::SolveLengthVector(Z, Point<double>(0, 0, 0)); //����������

			t = math::SolveLengthVector(Z, Point<double>(0, 0, 0));

			t = X * Y;
			t = X * Z;
			t = Z * Y;

			//��������
			if (!(X * Y <= 1E-6 && X * Z <= 1e-8 && Y * Z <= 1e-8))
			{
				printf_s("ERROR in self basis 1D!!");
				//Sleep(10000);
			}

			self_basis[0][0] = X.x; self_basis[0][1] = Y.x; self_basis[0][2] = Z.x;
			self_basis[1][0] = X.y; self_basis[1][1] = Y.y; self_basis[1][2] = Z.y;
			self_basis[2][0] = X.z; self_basis[2][1] = Y.z; self_basis[2][2] = Z.z;

			this->SetSelfBasis(self_basis);
		};
		bool CreateSelfBasis(Point<double> base_vect)
		{
			double self_basis[3][3];
			for (int i = 0; i < 3; i++)
				for (int j = 0; j < 3; j++)
					self_basis[i][j] = 0.0;

			Point<double> X, Y, Z; //��� ���� ����� ���� (�� ��������� ������ ��������� � (0,0,0))

			//������ ���� X
			X = base_vect;
			X /= math::SolveLengthVector(X, Point<double>(0, 0, 0)); //����������

			double t = math::SolveLengthVector(X, Point<double>(0, 0, 0));

			//������ ���� Y �� ����������� ��������������� � �
			if (X.z != 0)
			{
				Y.x = 1.0;
				Y.y = 1.0;
				Y.z = (-X.x - X.y) / X.z;
			}
			else {
				if (X.y != 0)
				{
					Y.x = 1.0;
					Y.y = (-X.x - X.z) / X.y;
					Y.z = 1.0;
				}
				else {
					Y.x = (-X.y - X.z) / X.x;
					Y.y = 1.0;
					Y.z = 1.0;
				}
			}
			Y /= math::SolveLengthVector(Y, Point<double>(0, 0, 0)); //����������

			t = math::SolveLengthVector(Y, Point<double>(0, 0, 0));

			//������ ���� Z �� ����������� ��������������� � � � Y
			if (X.z != 0 && (Y.y - X.y * Y.z / X.z) != 0)
			{
				Z.x = 1.0;
				Z.y = (X.x * Y.z / X.z - Y.x) * Z.x / (Y.y - X.y * Y.z / X.z);
				Z.z = (-X.x * Z.x - X.y * Z.y) / X.z;
			}
			else {
				if (X.y != 0 && (Y.z - X.z * Y.y / X.y) != 0)
				{
					Z.x = 1.0;
					Z.z = (X.x * Y.y / X.y - Y.x) * Z.x / (Y.z - X.z * Y.y / X.y);
					Z.y = (-X.x * Z.x - X.z * Z.z) / X.y;
				}
				else {
					Z.z = 1.0;
					Z.y = (X.z * Y.x / X.x - Y.z) / (Y.y - X.y * Y.x / X.x);
					Z.x = (-X.z - X.y * Z.y) / X.x;
				}
			}
			Z /= math::SolveLengthVector(Z, Point<double>(0, 0, 0)); //����������

			t = math::SolveLengthVector(Z, Point<double>(0, 0, 0));

			t = X * Y;
			t = X * Z;
			t = Z * Y;

			//��������
			if (!(X * Y <= 1E-6 && X * Z <= 1e-8 && Y * Z <= 1e-8))
			{
				printf_s("ERROR in self basis 1D!!");
				//Sleep(10000);
			}

			self_basis[0][0] = X.x; self_basis[0][1] = Y.x; self_basis[0][2] = Z.x;
			self_basis[1][0] = X.y; self_basis[1][1] = Y.y; self_basis[1][2] = Z.y;
			self_basis[2][0] = X.z; self_basis[2][1] = Y.z; self_basis[2][2] = Z.z;

			this->SetSelfBasis(self_basis);
			return true;
		};
	public:
		Segment()
		{
			this->id_nodes.resize(2);
			this->nodes.resize(2);
		};
		Segment(std::vector<Point<double>*>& P)
		{
			try
			{
				this->UpdateNodes(P);
			}
			catch (const std::exception&)
			{
				printf_s("Error: GeometryShape.cpp/Triangle(std::vector<Point<double> *> &P)");
			}
		}
		Segment(Point<double> P1, Point<double> P2)
		{
			try
			{
				std::vector<Point<double>*> P(2);
				P[0] = &P1;
				P[1] = &P2;
				this->UpdateNodes(P);
			}
			catch (const std::exception&)
			{
				printf_s("Error: GeometryShape.cpp/Triangle(std::vector<Point<double> *> &P)");
			}
		}
		Segment(std::vector<int>& id_nodes, std::vector<Point<double>*>& xyz_nodes)
		{
			try
			{
				math::MakeCopyVector_A_into_B(id_nodes, this->id_nodes);
				math::MakeCopyVector_A_into_B(xyz_nodes, this->nodes);
			}
			catch (const std::exception&)
			{
				printf_s("Error: GeometryShape.cpp/geometry::Shape(std::vector<int> &id_nodes, std::vector<Point<double>*> &xyz_nodes)\n");
			}
		};
		double GetVolume()
		{
			try
			{
				return math::SolveLengthVector(this->GetNode(0), this->GetNode(1));
			}
			catch (const std::exception&)
			{
				printf_s("Error: GeometryShape.cpp/GetVolume()");
				return 0.0;
			}
		}

		void SetGeometry(std::vector<int>& id_nodes, std::vector<Point<double>*>& P)
		{
			try
			{

				this->SetIdNodes(id_nodes);
				this->UpdateNodes(P);
			}
			catch (const std::exception&)
			{
				printf_s("Error: GeometryShape.h/SetGeometry(std::vector<std::vector<int>>& id_nodes, std::vector<Point<double>*>& P)");
			}
		};

		void CreateBasis()
		{
			if (this->self_basis.size() != 0)
				return;

			this->self_basis.resize(3);
			this->self_basis[0].resize(3);
			this->self_basis[1].resize(3);
			this->self_basis[2].resize(3);

			CreateSelfBasis();
		}
		void CreateReversBasis()
		{
			if (self_basis.size() == 0)
				this->CreateSelfBasis();

			this->self_basis_revers.resize(3);
			this->self_basis_revers[0].resize(3);
			this->self_basis_revers[1].resize(3);
			this->self_basis_revers[2].resize(3);

			this->self_basis_revers[0][0] = (this->self_basis[1][1] * this->self_basis[2][2] - this->self_basis[1][2] * this->self_basis[2][1]) / (this->self_basis[2][0] * this->self_basis[0][1] * this->self_basis[1][2] - this->self_basis[2][0] * this->self_basis[0][2] * this->self_basis[1][1] - this->self_basis[1][0] * this->self_basis[0][1] * this->self_basis[2][2] + this->self_basis[1][0] * this->self_basis[0][2] * this->self_basis[2][1] + this->self_basis[0][0] * this->self_basis[1][1] * this->self_basis[2][2] - this->self_basis[0][0] * this->self_basis[1][2] * this->self_basis[2][1]);
			this->self_basis_revers[0][1] = -(this->self_basis[0][1] * this->self_basis[2][2] - this->self_basis[0][2] * this->self_basis[2][1]) / (this->self_basis[2][0] * this->self_basis[0][1] * this->self_basis[1][2] - this->self_basis[2][0] * this->self_basis[0][2] * this->self_basis[1][1] - this->self_basis[1][0] * this->self_basis[0][1] * this->self_basis[2][2] + this->self_basis[1][0] * this->self_basis[0][2] * this->self_basis[2][1] + this->self_basis[0][0] * this->self_basis[1][1] * this->self_basis[2][2] - this->self_basis[0][0] * this->self_basis[1][2] * this->self_basis[2][1]);
			this->self_basis_revers[0][2] = (this->self_basis[0][1] * this->self_basis[1][2] - this->self_basis[0][2] * this->self_basis[1][1]) / (this->self_basis[2][0] * this->self_basis[0][1] * this->self_basis[1][2] - this->self_basis[2][0] * this->self_basis[0][2] * this->self_basis[1][1] - this->self_basis[1][0] * this->self_basis[0][1] * this->self_basis[2][2] + this->self_basis[1][0] * this->self_basis[0][2] * this->self_basis[2][1] + this->self_basis[0][0] * this->self_basis[1][1] * this->self_basis[2][2] - this->self_basis[0][0] * this->self_basis[1][2] * this->self_basis[2][1]);

			this->self_basis_revers[1][0] = -(-this->self_basis[2][0] * this->self_basis[1][2] + this->self_basis[1][0] * this->self_basis[2][2]) / (this->self_basis[2][0] * this->self_basis[0][1] * this->self_basis[1][2] - this->self_basis[2][0] * this->self_basis[0][2] * this->self_basis[1][1] - this->self_basis[1][0] * this->self_basis[0][1] * this->self_basis[2][2] + this->self_basis[1][0] * this->self_basis[0][2] * this->self_basis[2][1] + this->self_basis[0][0] * this->self_basis[1][1] * this->self_basis[2][2] - this->self_basis[0][0] * this->self_basis[1][2] * this->self_basis[2][1]);
			this->self_basis_revers[1][1] = (-this->self_basis[2][0] * this->self_basis[0][2] + this->self_basis[0][0] * this->self_basis[2][2]) / (this->self_basis[2][0] * this->self_basis[0][1] * this->self_basis[1][2] - this->self_basis[2][0] * this->self_basis[0][2] * this->self_basis[1][1] - this->self_basis[1][0] * this->self_basis[0][1] * this->self_basis[2][2] + this->self_basis[1][0] * this->self_basis[0][2] * this->self_basis[2][1] + this->self_basis[0][0] * this->self_basis[1][1] * this->self_basis[2][2] - this->self_basis[0][0] * this->self_basis[1][2] * this->self_basis[2][1]);
			this->self_basis_revers[1][2] = -(-this->self_basis[1][0] * this->self_basis[0][2] + this->self_basis[0][0] * this->self_basis[1][2]) / (this->self_basis[2][0] * this->self_basis[0][1] * this->self_basis[1][2] - this->self_basis[2][0] * this->self_basis[0][2] * this->self_basis[1][1] - this->self_basis[1][0] * this->self_basis[0][1] * this->self_basis[2][2] + this->self_basis[1][0] * this->self_basis[0][2] * this->self_basis[2][1] + this->self_basis[0][0] * this->self_basis[1][1] * this->self_basis[2][2] - this->self_basis[0][0] * this->self_basis[1][2] * this->self_basis[2][1]);

			this->self_basis_revers[2][0] = (-this->self_basis[2][0] * this->self_basis[1][1] + this->self_basis[1][0] * this->self_basis[2][1]) / (this->self_basis[2][0] * this->self_basis[0][1] * this->self_basis[1][2] - this->self_basis[2][0] * this->self_basis[0][2] * this->self_basis[1][1] - this->self_basis[1][0] * this->self_basis[0][1] * this->self_basis[2][2] + this->self_basis[1][0] * this->self_basis[0][2] * this->self_basis[2][1] + this->self_basis[0][0] * this->self_basis[1][1] * this->self_basis[2][2] - this->self_basis[0][0] * this->self_basis[1][2] * this->self_basis[2][1]);
			this->self_basis_revers[2][1] = -(-this->self_basis[2][0] * this->self_basis[0][1] + this->self_basis[0][0] * this->self_basis[2][1]) / (this->self_basis[2][0] * this->self_basis[0][1] * this->self_basis[1][2] - this->self_basis[2][0] * this->self_basis[0][2] * this->self_basis[1][1] - this->self_basis[1][0] * this->self_basis[0][1] * this->self_basis[2][2] + this->self_basis[1][0] * this->self_basis[0][2] * this->self_basis[2][1] + this->self_basis[0][0] * this->self_basis[1][1] * this->self_basis[2][2] - this->self_basis[0][0] * this->self_basis[1][2] * this->self_basis[2][1]);
			this->self_basis_revers[2][2] = (-this->self_basis[1][0] * this->self_basis[0][1] + this->self_basis[0][0] * this->self_basis[1][1]) / (this->self_basis[2][0] * this->self_basis[0][1] * this->self_basis[1][2] - this->self_basis[2][0] * this->self_basis[0][2] * this->self_basis[1][1] - this->self_basis[1][0] * this->self_basis[0][1] * this->self_basis[2][2] + this->self_basis[1][0] * this->self_basis[0][2] * this->self_basis[2][1] + this->self_basis[0][0] * this->self_basis[1][1] * this->self_basis[2][2] - this->self_basis[0][0] * this->self_basis[1][2] * this->self_basis[2][1]);
		};
		void SetSelfBasis(double basis[3][3])
		{
			this->self_basis.resize(3);
			this->self_basis[0].resize(3);
			this->self_basis[1].resize(3);
			this->self_basis[2].resize(3);
			for (int i = 0; i < 3; i++)
				for (int j = 0; j < 3; j++)
					this->self_basis[i][j] = basis[i][j];

			CreateReversBasis();

			this->self_nodes.resize(this->GetNodesCount());
			for (int i = 0; i < this->nodes.size(); i++)
			{
				self_nodes[i].x = this->self_basis_revers[0][0] * this->nodes[i]->x + this->self_basis_revers[0][1] * this->nodes[i]->y + this->self_basis_revers[0][2] * this->nodes[i]->z;
				self_nodes[i].y = this->self_basis_revers[1][0] * this->nodes[i]->x + this->self_basis_revers[1][1] * this->nodes[i]->y + this->self_basis_revers[1][2] * this->nodes[i]->z;
				self_nodes[i].z = this->self_basis_revers[2][0] * this->nodes[i]->x + this->self_basis_revers[2][1] * this->nodes[i]->y + this->self_basis_revers[2][2] * this->nodes[i]->z;
			}
		};
		std::vector<std::vector<double>>* GetSelfReverseBasis()
		{
			if (this->self_basis_revers.size() == 0)
			{
				this->CreateReversBasis();
			}
			return &(this->self_basis_revers);
		}
		std::vector<std::vector<double>>* GetSelfBasis()
		{
			if (this->self_basis.size() == 0)
			{
				this->CreateBasis();
			}
			return &(this->self_basis);
		}
		Point<double> MakeTransferSelfIntoXyz(Point<double> X)
		{
			Point<double> tmp;
			if (this->self_basis.size() == 0)
			{
				CreateBasis();
			}

			tmp.x = this->self_basis[0][0] * X.x + this->self_basis[0][1] * X.y + this->self_basis[0][2] * X.z;
			tmp.y = this->self_basis[1][0] * X.x + this->self_basis[1][1] * X.y + this->self_basis[1][2] * X.z;
			tmp.z = this->self_basis[2][0] * X.x + this->self_basis[2][1] * X.y + this->self_basis[2][2] * X.z;
			return tmp;
		};
		Point<double> MakeTransferXyzIntoSelf(Point<double> X)
		{
			if (this->self_basis.size() == 0)
			{
				CreateBasis();
			}
			if (this->self_basis_revers.size() == 0)
				CreateReversBasis();

			Point<double> tmp;
			tmp.x = this->self_basis_revers[0][0] * X.x + this->self_basis_revers[0][1] * X.y + this->self_basis_revers[0][2] * X.z;
			tmp.y = this->self_basis_revers[1][0] * X.x + this->self_basis_revers[1][1] * X.y + this->self_basis_revers[1][2] * X.z;
			tmp.z = this->self_basis_revers[2][0] * X.x + this->self_basis_revers[2][1] * X.y + this->self_basis_revers[2][2] * X.z;
			return tmp;
		};
		Point<double> GetSelfNode(int id_local)
		{
			if (self_nodes.size() == 0)
			{
				this->self_nodes.resize(this->GetNodesCount());
				for (int i = 0; i < this->self_nodes.size(); i++)
				{
					this->self_nodes[i] = this->MakeTransferXyzIntoSelf(this->GetNode(i));
				}
			}
			return self_nodes[id_local];
		}

		void SetIntegrationLaw(int points_count)
		{
			this->integration_parameters.w.resize(points_count);
			this->integration_parameters.x.resize(points_count);

			std::vector<double> xj(points_count);
			std::vector<double> q(points_count);

			/*    Points     */
			switch (points_count)
			{
			case 10:
				xj[0] = 0.9739065318;
				xj[1] = 0.8650633629;
				xj[2] = 0.6794095691;
				xj[3] = 0.4333953941;
				xj[4] = 0.1488743390;
				xj[5] = -xj[0];
				xj[6] = -xj[1];
				xj[7] = -xj[2];
				xj[8] = -xj[3];
				xj[9] = -xj[3];
				break;
			case 9:
				xj[0] = 0.0;
				xj[1] = 0.3242534234;
				xj[2] = 0.6133714327;
				xj[3] = 0.8360311073;
				xj[4] = 0.9681602395;
				xj[5] = -xj[1];
				xj[6] = -xj[2];
				xj[7] = -xj[3];
				xj[8] = -xj[4];
				break;
			case 8:
				xj[0] = 0.9602898565;
				xj[1] = 0.7966664774;
				xj[2] = 0.5255324099;
				xj[3] = 0.1834346425;
				xj[4] = -xj[0];
				xj[5] = -xj[1];
				xj[6] = -xj[2];
				xj[7] = -xj[3];
				break;
			case 7:
				xj[0] = 0.0;
				xj[1] = 0.4058451514;
				xj[2] = 0.7415311856;
				xj[3] = 0.9491079123;
				xj[4] = -xj[1];
				xj[5] = -xj[2];
				xj[6] = -xj[3];
				break;
			case 6:
				xj[0] = 0.2386191861;
				xj[1] = 0.6612093865;
				xj[2] = 0.9324695142;
				xj[3] = -xj[0];
				xj[4] = -xj[1];
				xj[5] = -xj[2];
				break;
			case 5:
				xj[0] = 0.;
				xj[1] = 1 / 21. * sqrt(245 - 14 * sqrt(70));
				xj[2] = 1 / 21. * sqrt(245 + 14 * sqrt(70));
				xj[3] = -xj[1];
				xj[4] = -xj[2];
				break;
			case 4:
				xj[0] = sqrt((3.0 - 2.0 * sqrt(6.0 / 5.0)) / 7.0);
				xj[1] = sqrt((3.0 + 2.0 * sqrt(6.0 / 5.0)) / 7.0);
				xj[2] = -xj[1];
				xj[3] = -xj[0];
				break;
			case 3:
				xj[0] = 0;
				xj[1] = 0.774596669241483;
				xj[2] = -xj[1];
				break;
			case 2:
				xj[0] = 1 / 3. * sqrt(3.);
				xj[1] = -xj[0];
				break;
			case 1:
				xj[0] = 0;
				break;
			}
			/*    Weights    */
			switch (points_count)
			{
			case 10:
				q[0] = 0.06667134202;
				q[1] = 0.1494513565;
				q[2] = 0.2190863541;
				q[3] = 0.2692667242;
				q[4] = 0.2955242232;
				q[5] = q[0];
				q[6] = q[1];
				q[7] = q[2];
				q[8] = q[3];
				q[9] = q[3];
				break;
			case 9:
				q[0] = 0.3356031191;
				q[1] = 0.3048968104;
				q[2] = 0.2710276447;
				q[3] = 0.1638265280;
				q[4] = 0.09199001997;
				q[5] = q[1];
				q[6] = q[2];
				q[7] = q[3];
				q[8] = q[4];
				break;
			case 8:
				q[0] = 0.1012285363;
				q[1] = 0.2223810344;
				q[2] = 0.3137066459;
				q[3] = 0.3626837833;
				q[4] = q[0];
				q[5] = q[1];
				q[6] = q[2];
				q[7] = q[3];
				break;
			case 7:
				q[0] = 0.4179591838;
				q[1] = 0.3818300507;
				q[2] = 0.2797053914;
				q[3] = 0.1294849662;
				q[4] = q[1];
				q[5] = q[2];
				q[6] = q[3];
				break;
			case 6:
				q[0] = 0.4679139346;
				q[1] = 0.3607615729;
				q[2] = 0.1713244924;
				q[3] = q[0];
				q[4] = q[1];
				q[5] = q[2];
				break;
			case 5:
				q[0] = 0.5688888888888888888;
				q[1] = 0.4786286705;
				q[2] = 0.2369268851;
				q[3] = q[1];
				q[4] = q[2];
				break;
			case 4:
				q[0] = (18 + sqrt(30.0)) / 36.0;
				q[1] = (18 - sqrt(30.0)) / 36.0;
				q[2] = q[1];
				q[3] = q[0];
				break;
			case 3:
				q[0] = 0.888888888888889;
				q[1] = 0.555555555555556;
				q[2] = q[1];
				break;
			case 2:
				q[0] = 1.0;
				q[1] = q[0];
				break;
			case 1:
				q[0] = 2.0;
				break;
			}

			this->CreateBasis();

			double a = self_nodes[0].x, b = self_nodes[1].x;
			if (self_nodes[0].x > self_nodes[1].x)
			{
				a = self_nodes[1].x; b = self_nodes[0].x;
			}
			double c = self_nodes[0].y, d = self_nodes[1].y;
			for (int i = 0; i < points_count; ++i)
			{
				double xi = (a + b + xj[i] * (b - a)) / 2.0;
				this->integration_parameters.w[i] = q[i] * (b - a) / 2.0;
				this->integration_parameters.x[i] = this->MakeTransferSelfIntoXyz(Point<double>(xi, self_nodes[0].y, self_nodes[0].z));
			}
		}

		bool IsContainThePoint(Point<double> A)
		{
			try
			{
				Point<double> Xtmp = this->MakeTransferXyzIntoSelf(A);
				if (self_nodes[0].x <= Xtmp.x && self_nodes[1].x >= Xtmp.x || self_nodes[1].x <= Xtmp.x && self_nodes[0].x >= Xtmp.x)
					return true;
				return false;
			}
			catch (const std::exception&)
			{
				printf_s("Error: GeometryShape.cpp/IsContainThePoint(Point<double> A)");
				return false;
			}
		}
		bool IsContainThePoint(Point<double> A, double& length)
		{
			try
			{
				Point<double> Xtmp = this->MakeTransferXyzIntoSelf(A);
				if (self_nodes[0].x <= Xtmp.x && self_nodes[1].x >= Xtmp.x || self_nodes[1].x <= Xtmp.x && self_nodes[0].x >= Xtmp.x)
					return true;

				Point<double> Xtmp_progection = math::MakeProgectionOfPointIntoLine(self_nodes[0], self_nodes[1], Xtmp);
				if (this->IsContainThePoint(Xtmp_progection))
				{
					length = math::SolveLengthVector(Xtmp, Xtmp_progection);
				}
				else
				{
					double len0 = math::SolveLengthVector(Xtmp, self_nodes[0]);
					double len1 = math::SolveLengthVector(Xtmp, self_nodes[1]);
					length = len0 < len1 ? len0 : len1;
				}
				return false;
			}
			catch (const std::exception&)
			{
				printf_s("Error: GeometryShape.cpp/IsContainThePoint(Point<double> A)");
				return false;
			}
		}
	};

	class Triangle : public Shape
	{
		std::vector<std::vector<double>> self_basis;
		std::vector<std::vector<double>> self_basis_revers;
		std::vector<Point<double>> self_nodes;

		void CreateSelfBasis()
		{
			double self_basis[3][3];
			for (int i = 0; i < 3; i++)
				for (int j = 0; j < 3; j++)
					self_basis[i][j] = 0.0;

			Point<double> X, Y, Z; //��� ���� ����� ���� (�� ��������� ������ ��������� � (0,0,0))

			//������ ���� X
			X = *this->nodes[1] - *this->nodes[0];
			X /= math::SolveLengthVector(X, Point<double>(0, 0, 0)); //����������

			double t = math::SolveLengthVector(X, Point<double>(0, 0, 0));

			//������ ���� Y �� ����������� ��������������� � � � �������������� � ��������� ������
			Plane tr(*this->nodes[0], *this->nodes[1], *this->nodes[2]);
			double _kx, _ky, _kz, _kd, kx, ky, kz, kd;
			tr.GetCoefficients(_kx, _ky, _kz, _kd);
			kx = _kx / sqrt(_kx*_kx + _ky * _ky + _kz * _kz);
			ky = _ky / sqrt(_kx*_kx + _ky * _ky + _kz * _kz);
			kz = _kz / sqrt(_kx*_kx + _ky * _ky + _kz * _kz);
			kd = _kd / sqrt(_kx*_kx + _ky * _ky + _kz * _kz);

			t = tr.SolveDistanceToPoint(Point<double>(0, 0, 0));

			if (!math::IsEqual(X.z, 0.0))
			{
				if (!math::IsEqual((ky - X.y*kz / X.z), 0.0))
				{
					Y.x = 1.0;
					Y.y = Y.x*(-kx + X.x*kz / X.z) / (ky - X.y*kz / X.z);
					Y.z = (-X.y*Y.y - X.x*Y.x) / X.z;
				}
				else {
					if (!math::IsEqual((kx - X.x*kz / X.z), 0.0))
					{
						Y.y = 1.0;
						Y.x = Y.y*(-ky + X.y*kz / X.z) / (kx - X.x*kz / X.z);
						Y.z = (-X.y*Y.y - X.x*Y.x) / X.z;
					}
					else goto next_if1;
				}
			}
			else {
			next_if1:
				if (!math::IsEqual(X.y, 0.0))
				{
					if (!math::IsEqual((kz - X.z*ky / X.y), 0.0))
					{
						Y.x = 1.0;
						Y.z = Y.x*(-kx + X.x*ky / X.y) / (kz - X.z*ky / X.y);
						Y.y = (-X.z*Y.z - X.x*Y.x) / X.y;
					}
					else {
						if (!math::IsEqual((kx - X.x*ky / X.y), 0.0))
						{
							Y.z = 1.0;
							Y.x = Y.z*(-kz + X.z*ky / X.y) / (kx - X.x*ky / X.y);
							Y.y = (-X.z*Y.z - X.x*Y.x) / X.y;
						}
						else goto next_if2;
					}
				}
				else {
				next_if2:
					if (!math::IsEqual((ky - X.y*kx / X.x), 0.0))
					{
						Y.z = 1.0;
						Y.y = Y.z*(-kz + X.z*kx / X.x) / (ky - X.y*kx / X.x);
						Y.x = (-X.z*Y.z - X.y*Y.y) / X.x;
					}
					else {
						if (!math::IsEqual((kz - X.z*kx / X.x), 0.0))
						{
							Y.y = 1.0;
							Y.z = Y.y*(-ky + X.y*kx / X.x) / (kz - X.z*kx / X.x);
							Y.x = (-X.z*Y.z - X.y*Y.y) / X.x;
						}
						else goto next_if2;
					}
				}
			}
			Y /= math::SolveLengthVector(Y, Point<double>(0, 0, 0)); //����������

			t = math::SolveLengthVector(Y, Point<double>(0, 0, 0));

			//������ ���� Z �� ����������� ��������������� � � � Y
			if (X.z != 0)
			{
				Z.x = 1.0;
				Z.y = (X.x*Y.z / X.z - Y.x)*Z.x / (Y.y - X.y*Y.z / X.z);
				Z.z = (-X.x*Z.x - X.y*Z.y) / X.z;

				if (Z.z != Z.z || Z.x != Z.x || Z.y != Z.y || math::IsEqual((Y.y - X.y*Y.z / X.z), 0.0))
				{
					Z.y = 1.0;
					Z.x = (X.y*Y.z / X.z - Y.y)*Z.y / (Y.x - X.x*Y.z / X.z);
					Z.z = (-X.x*Z.x - X.y*Z.y) / X.z;
				}
			}
			else {
				if (X.y != 0)
				{
					Z.x = 1.0;
					Z.z = (X.x*Y.y / X.y - Y.x)*Z.x / (Y.z - X.z*Y.y / X.y);
					Z.y = (-X.x*Z.x - X.z*Z.z) / X.y;

					if (Z.z != Z.z || Z.x != Z.x || Z.y != Z.y || math::IsEqual((Y.z - X.z*Y.y / X.y), 0.0))
					{
						Z.z = 1.0;
						Z.x = (X.z*Y.y / X.y - Y.z)*Z.z / (Y.x - X.x*Y.y / X.y);
						Z.y = (-X.x*Z.x - X.z*Z.z) / X.y;
					}
				}
				else {
					Z.z = 1.0;
					Z.y = (X.z*Y.x / X.x - Y.z) / (Y.y - X.y*Y.x / X.x);
					Z.x = (-X.z - X.y*Z.y) / X.x;

					if (Z.z != Z.z || Z.x != Z.x || Z.y != Z.y || math::IsEqual((Y.y - X.y*Y.x / X.x), 0.0))
					{
						Z.y = 1.0;
						Z.z = (X.y*Y.x / X.x - Y.y) / (Y.z - X.z*Y.x / X.x);
						Z.x = (-X.z - X.y*Z.y) / X.x;
					}
				}
			}
			Z /= math::SolveLengthVector(Z, Point<double>(0, 0, 0)); //����������


			t = X * Y;
			t = X * Z;
			t = Z * Y;

			//��������
			if (!(abs(X*Y) <= 1E-5 && abs(X*Z) <= 1e-5 && abs(Y*Z) <= 1e-5))
			{
				
				//printf_s("X(%.8lf, %.8lf, %.8lf)\n", X.x, X.y, X.z);
				//printf_s("Y(%.8lf, %.8lf, %.8lf)\n", Y.x, Y.y, Y.z);
				//printf_s("Z(%.8lf, %.8lf, %.8lf)\n", Z.x, Z.y, Z.z);
				//printf_s("%.8lf, %.8lf, %.8lf, %.8lf\n", _kx, _ky, _kz, _kd);
				//printf_s("%.8lf, %.8lf, %.8lf, %.8lf\n", kx, ky, kz, kd);
				bool a = CreateSelfBasis(*this->nodes[2] - *this->nodes[0]);
				bool b = false;
				if (a == false)
					b = CreateSelfBasis(*this->nodes[2] - *this->nodes[1]);
				if (!(a == true || b == true))
					printf_s("ERROR in self basis 2D!! XY=%.2e, XZ=%.2e, YZ=%.2e\n", X*Y, X*Z, Y*Z);
				//Sleep(10000);
			}

			self_basis[0][0] = X.x; self_basis[0][1] = Y.x; self_basis[0][2] = Z.x;
			self_basis[1][0] = X.y; self_basis[1][1] = Y.y; self_basis[1][2] = Z.y;
			self_basis[2][0] = X.z; self_basis[2][1] = Y.z; self_basis[2][2] = Z.z;

			this->SetSelfBasis(self_basis);
		};
		bool CreateSelfBasis(Point<double> base_vect)
		{
			double self_basis[3][3];
			for (int i = 0; i < 3; i++)
				for (int j = 0; j < 3; j++)
					self_basis[i][j] = 0.0;

			Point<double> X, Y, Z; //��� ���� ����� ���� (�� ��������� ������ ��������� � (0,0,0))

			//������ ���� X
			X = base_vect;// *this->nodes[1] - *this->nodes[0];
			X /= math::SolveLengthVector(X, Point<double>(0, 0, 0)); //����������

			double t = math::SolveLengthVector(X, Point<double>(0, 0, 0));

			//������ ���� Y �� ����������� ��������������� � � � �������������� � ��������� ������
			Plane tr(*this->nodes[0], *this->nodes[1], *this->nodes[2]);
			double _kx, _ky, _kz, _kd, kx, ky, kz, kd;
			tr.GetCoefficients(_kx, _ky, _kz, _kd);
			kx = _kx / sqrt(_kx*_kx + _ky * _ky + _kz * _kz);
			ky = _ky / sqrt(_kx*_kx + _ky * _ky + _kz * _kz);
			kz = _kz / sqrt(_kx*_kx + _ky * _ky + _kz * _kz);
			kd = _kd / sqrt(_kx*_kx + _ky * _ky + _kz * _kz);

			t = tr.SolveDistanceToPoint(Point<double>(0, 0, 0));

			if (!math::IsEqual(X.z, 0.0))
			{
				if (!math::IsEqual((ky - X.y*kz / X.z), 0.0))
				{
					Y.x = 1.0;
					Y.y = Y.x*(-kx + X.x*kz / X.z) / (ky - X.y*kz / X.z);
					Y.z = (-X.y*Y.y - X.x*Y.x) / X.z;
				}
				else {
					if (!math::IsEqual((kx - X.x*kz / X.z), 0.0))
					{
						Y.y = 1.0;
						Y.x = Y.y*(-ky + X.y*kz / X.z) / (kx - X.x*kz / X.z);
						Y.z = (-X.y*Y.y - X.x*Y.x) / X.z;
					}
					else goto next_if1;
				}
			}
			else {
			next_if1:
				if (!math::IsEqual(X.y, 0.0))
				{
					if (!math::IsEqual((kz - X.z*ky / X.y), 0.0))
					{
						Y.x = 1.0;
						Y.z = Y.x*(-kx + X.x*ky / X.y) / (kz - X.z*ky / X.y);
						Y.y = (-X.z*Y.z - X.x*Y.x) / X.y;
					}
					else {
						if (!math::IsEqual((kx - X.x*ky / X.y), 0.0))
						{
							Y.z = 1.0;
							Y.x = Y.z*(-kz + X.z*ky / X.y) / (kx - X.x*ky / X.y);
							Y.y = (-X.z*Y.z - X.x*Y.x) / X.y;
						}
						else goto next_if2;
					}
				}
				else {
				next_if2:
					if (!math::IsEqual((ky - X.y*kx / X.x), 0.0))
					{
						Y.z = 1.0;
						Y.y = Y.z*(-kz + X.z*kx / X.x) / (ky - X.y*kx / X.x);
						Y.x = (-X.z*Y.z - X.y*Y.y) / X.x;
					}
					else {
						if (!math::IsEqual((kz - X.z*kx / X.x), 0.0))
						{
							Y.y = 1.0;
							Y.z = Y.y*(-ky + X.y*kx / X.x) / (kz - X.z*kx / X.x);
							Y.x = (-X.z*Y.z - X.y*Y.y) / X.x;
						}
						else goto next_if2;
					}
				}
			}
			Y /= math::SolveLengthVector(Y, Point<double>(0, 0, 0)); //����������

			t = math::SolveLengthVector(Y, Point<double>(0, 0, 0));

			//������ ���� Z �� ����������� ��������������� � � � Y
			if (X.z != 0)
			{
				Z.x = 1.0;
				Z.y = (X.x*Y.z / X.z - Y.x)*Z.x / (Y.y - X.y*Y.z / X.z);
				Z.z = (-X.x*Z.x - X.y*Z.y) / X.z;

				if (Z.z != Z.z || Z.x != Z.x || Z.y != Z.y || math::IsEqual((Y.y - X.y*Y.z / X.z), 0.0))
				{
					Z.y = 1.0;
					Z.x = (X.y*Y.z / X.z - Y.y)*Z.y / (Y.x - X.x*Y.z / X.z);
					Z.z = (-X.x*Z.x - X.y*Z.y) / X.z;
				}
			}
			else {
				if (X.y != 0)
				{
					Z.x = 1.0;
					Z.z = (X.x*Y.y / X.y - Y.x)*Z.x / (Y.z - X.z*Y.y / X.y);
					Z.y = (-X.x*Z.x - X.z*Z.z) / X.y;

					if (Z.z != Z.z || Z.x != Z.x || Z.y != Z.y || math::IsEqual((Y.z - X.z*Y.y / X.y), 0.0))
					{
						Z.z = 1.0;
						Z.x = (X.z*Y.y / X.y - Y.z)*Z.z / (Y.x - X.x*Y.y / X.y);
						Z.y = (-X.x*Z.x - X.z*Z.z) / X.y;
					}
				}
				else {
					Z.z = 1.0;
					Z.y = (X.z*Y.x / X.x - Y.z) / (Y.y - X.y*Y.x / X.x);
					Z.x = (-X.z - X.y*Z.y) / X.x;

					if (Z.z != Z.z || Z.x != Z.x || Z.y != Z.y || math::IsEqual((Y.y - X.y*Y.x / X.x), 0.0))
					{
						Z.y = 1.0;
						Z.z = (X.y*Y.x / X.x - Y.y) / (Y.z - X.z*Y.x / X.x);
						Z.x = (-X.z - X.y*Z.y) / X.x;
					}
				}
			}
			Z /= math::SolveLengthVector(Z, Point<double>(0, 0, 0)); //����������


			t = X * Y;
			t = X * Z;
			t = Z * Y;

			//��������
			if (!(abs(X*Y) <= 1E-6 && abs(X*Z) <= 1e-6 && abs(Y*Z) <= 1e-6))
			{
				printf_s("\t(add) ERROR in self basis 2D!! XY=%.2e, XZ=%.2e, YZ=%.2e\n", X*Y, X*Z, Y*Z);
				return false;
				//Sleep(10000);
			}

			self_basis[0][0] = X.x; self_basis[0][1] = Y.x; self_basis[0][2] = Z.x;
			self_basis[1][0] = X.y; self_basis[1][1] = Y.y; self_basis[1][2] = Z.y;
			self_basis[2][0] = X.z; self_basis[2][1] = Y.z; self_basis[2][2] = Z.z;

			this->SetSelfBasis(self_basis);
			//printf_s("\t(add) XY=%.2e, XZ=%.2e, YZ=%.2e\n", X*Y, X*Z, Y*Z);
			return true;
		};

	public:
		std::vector<std::vector<double>> alpha;

		Triangle()
		{
			this->id_nodes.resize(3);
			this->nodes.resize(3);
		};
		Triangle(std::vector<Point<double> *> &P)
		{
			try
			{
				this->UpdateNodes(P);
			}
			catch (const std::exception&)
			{
				printf_s("Error: GeometryShape.cpp/Triangle(std::vector<Point<double> *> &P)");
			}
		}
		Triangle(Point<double> P1, Point<double> P2, Point<double> P3)
		{
			try
			{
				std::vector<Point<double> *> P(3);
				P[0] = &P1;
				P[1] = &P2;
				P[2] = &P3;
				this->UpdateNodes(P);
			}
			catch (const std::exception&)
			{
				printf_s("Error: GeometryShape.cpp/Triangle(std::vector<Point<double> *> &P)");
			}
		}
		Triangle(std::vector<int>& id_nodes, std::vector<Point<double>*>& xyz_nodes)
		{
			try
			{
				math::MakeCopyVector_A_into_B(id_nodes, this->id_nodes);
				math::MakeCopyVector_A_into_B(xyz_nodes, this->nodes);
			}
			catch (const std::exception&)
			{
				printf_s("Error: GeometryShape.cpp/geometry::Shape(std::vector<int> &id_nodes, std::vector<Point<double>*> &xyz_nodes)\n");
			}
		};
		double GetVolume()
		{
			try
			{
				Point<double> _nodes[3];
				for (int i = 0; i < 3; i++)
				{
					_nodes[i] = this->GetNode(i);
				}
				double Mx[3][3] = { { 1, _nodes[0].y, _nodes[0].z },
				{ 1, _nodes[1].y, _nodes[1].z },
				{ 1, _nodes[2].y, _nodes[2].z } };
				double My[3][3] = { { _nodes[0].x, 1, _nodes[0].z },
				{ _nodes[1].x, 1, _nodes[1].z },
				{ _nodes[2].x, 1, _nodes[2].z } };
				double Mz[3][3] = { { _nodes[0].x, _nodes[0].y, 1 },
				{ _nodes[1].x, _nodes[1].y, 1 },
				{ _nodes[2].x, _nodes[2].y, 1 } };
				double Sx = 0.5 * math::GetDeterminantForMatrix3x3(Mx);
				double Sy = 0.5 * math::GetDeterminantForMatrix3x3(My);
				double Sz = 0.5 * math::GetDeterminantForMatrix3x3(Mz);

				return sqrt(Sx*Sx + Sy * Sy + Sz * Sz);
			}
			catch (const std::exception&)
			{
				printf_s("Error: GeometryShape.cpp/GetVolume()");
				return 0.0;
			}
		}
		Point<double> GetNormal()
		{
			Point<double> _n;
			double _t;
			math::GetPlaneEquation(*this->nodes[0], *this->nodes[1], *this->nodes[2], _n.x, _n.y, _n.z, _t);
			return _n;
		}

		void CreateBasis()
		{
			if (this->self_basis.size() != 0)
				return;

			this->self_basis.resize(3);
			this->self_basis[0].resize(3);
			this->self_basis[1].resize(3);
			this->self_basis[2].resize(3);

			CreateSelfBasis();
		}
		void CreateReversBasis()
		{
			this->self_basis_revers.resize(3);
			this->self_basis_revers[0].resize(3);
			this->self_basis_revers[1].resize(3);
			this->self_basis_revers[2].resize(3);

			this->self_basis_revers[0][0] = (this->self_basis[1][1] * this->self_basis[2][2] - this->self_basis[1][2] * this->self_basis[2][1]) / (this->self_basis[2][0] * this->self_basis[0][1] * this->self_basis[1][2] - this->self_basis[2][0] * this->self_basis[0][2] * this->self_basis[1][1] - this->self_basis[1][0] * this->self_basis[0][1] * this->self_basis[2][2] + this->self_basis[1][0] * this->self_basis[0][2] * this->self_basis[2][1] + this->self_basis[0][0] * this->self_basis[1][1] * this->self_basis[2][2] - this->self_basis[0][0] * this->self_basis[1][2] * this->self_basis[2][1]);
			this->self_basis_revers[0][1] = -(this->self_basis[0][1] * this->self_basis[2][2] - this->self_basis[0][2] * this->self_basis[2][1]) / (this->self_basis[2][0] * this->self_basis[0][1] * this->self_basis[1][2] - this->self_basis[2][0] * this->self_basis[0][2] * this->self_basis[1][1] - this->self_basis[1][0] * this->self_basis[0][1] * this->self_basis[2][2] + this->self_basis[1][0] * this->self_basis[0][2] * this->self_basis[2][1] + this->self_basis[0][0] * this->self_basis[1][1] * this->self_basis[2][2] - this->self_basis[0][0] * this->self_basis[1][2] * this->self_basis[2][1]);
			this->self_basis_revers[0][2] = (this->self_basis[0][1] * this->self_basis[1][2] - this->self_basis[0][2] * this->self_basis[1][1]) / (this->self_basis[2][0] * this->self_basis[0][1] * this->self_basis[1][2] - this->self_basis[2][0] * this->self_basis[0][2] * this->self_basis[1][1] - this->self_basis[1][0] * this->self_basis[0][1] * this->self_basis[2][2] + this->self_basis[1][0] * this->self_basis[0][2] * this->self_basis[2][1] + this->self_basis[0][0] * this->self_basis[1][1] * this->self_basis[2][2] - this->self_basis[0][0] * this->self_basis[1][2] * this->self_basis[2][1]);

			this->self_basis_revers[1][0] = -(-this->self_basis[2][0] * this->self_basis[1][2] + this->self_basis[1][0] * this->self_basis[2][2]) / (this->self_basis[2][0] * this->self_basis[0][1] * this->self_basis[1][2] - this->self_basis[2][0] * this->self_basis[0][2] * this->self_basis[1][1] - this->self_basis[1][0] * this->self_basis[0][1] * this->self_basis[2][2] + this->self_basis[1][0] * this->self_basis[0][2] * this->self_basis[2][1] + this->self_basis[0][0] * this->self_basis[1][1] * this->self_basis[2][2] - this->self_basis[0][0] * this->self_basis[1][2] * this->self_basis[2][1]);
			this->self_basis_revers[1][1] = (-this->self_basis[2][0] * this->self_basis[0][2] + this->self_basis[0][0] * this->self_basis[2][2]) / (this->self_basis[2][0] * this->self_basis[0][1] * this->self_basis[1][2] - this->self_basis[2][0] * this->self_basis[0][2] * this->self_basis[1][1] - this->self_basis[1][0] * this->self_basis[0][1] * this->self_basis[2][2] + this->self_basis[1][0] * this->self_basis[0][2] * this->self_basis[2][1] + this->self_basis[0][0] * this->self_basis[1][1] * this->self_basis[2][2] - this->self_basis[0][0] * this->self_basis[1][2] * this->self_basis[2][1]);
			this->self_basis_revers[1][2] = -(-this->self_basis[1][0] * this->self_basis[0][2] + this->self_basis[0][0] * this->self_basis[1][2]) / (this->self_basis[2][0] * this->self_basis[0][1] * this->self_basis[1][2] - this->self_basis[2][0] * this->self_basis[0][2] * this->self_basis[1][1] - this->self_basis[1][0] * this->self_basis[0][1] * this->self_basis[2][2] + this->self_basis[1][0] * this->self_basis[0][2] * this->self_basis[2][1] + this->self_basis[0][0] * this->self_basis[1][1] * this->self_basis[2][2] - this->self_basis[0][0] * this->self_basis[1][2] * this->self_basis[2][1]);

			this->self_basis_revers[2][0] = (-this->self_basis[2][0] * this->self_basis[1][1] + this->self_basis[1][0] * this->self_basis[2][1]) / (this->self_basis[2][0] * this->self_basis[0][1] * this->self_basis[1][2] - this->self_basis[2][0] * this->self_basis[0][2] * this->self_basis[1][1] - this->self_basis[1][0] * this->self_basis[0][1] * this->self_basis[2][2] + this->self_basis[1][0] * this->self_basis[0][2] * this->self_basis[2][1] + this->self_basis[0][0] * this->self_basis[1][1] * this->self_basis[2][2] - this->self_basis[0][0] * this->self_basis[1][2] * this->self_basis[2][1]);
			this->self_basis_revers[2][1] = -(-this->self_basis[2][0] * this->self_basis[0][1] + this->self_basis[0][0] * this->self_basis[2][1]) / (this->self_basis[2][0] * this->self_basis[0][1] * this->self_basis[1][2] - this->self_basis[2][0] * this->self_basis[0][2] * this->self_basis[1][1] - this->self_basis[1][0] * this->self_basis[0][1] * this->self_basis[2][2] + this->self_basis[1][0] * this->self_basis[0][2] * this->self_basis[2][1] + this->self_basis[0][0] * this->self_basis[1][1] * this->self_basis[2][2] - this->self_basis[0][0] * this->self_basis[1][2] * this->self_basis[2][1]);
			this->self_basis_revers[2][2] = (-this->self_basis[1][0] * this->self_basis[0][1] + this->self_basis[0][0] * this->self_basis[1][1]) / (this->self_basis[2][0] * this->self_basis[0][1] * this->self_basis[1][2] - this->self_basis[2][0] * this->self_basis[0][2] * this->self_basis[1][1] - this->self_basis[1][0] * this->self_basis[0][1] * this->self_basis[2][2] + this->self_basis[1][0] * this->self_basis[0][2] * this->self_basis[2][1] + this->self_basis[0][0] * this->self_basis[1][1] * this->self_basis[2][2] - this->self_basis[0][0] * this->self_basis[1][2] * this->self_basis[2][1]);
		};
		void SetSelfBasis(double basis[3][3])
		{
			this->self_basis.resize(3);
			this->self_basis[0].resize(3);
			this->self_basis[1].resize(3);
			this->self_basis[2].resize(3);
			for (int i = 0; i < 3; i++)
				for (int j = 0; j < 3; j++)
					this->self_basis[i][j] = basis[i][j];

			CreateReversBasis();

			this->self_nodes.resize(this->GetNodesCount());
			for (int i = 0; i < this->nodes.size(); i++)
			{
				self_nodes[i].x = this->self_basis_revers[0][0] * this->nodes[i]->x + this->self_basis_revers[0][1] * this->nodes[i]->y + this->self_basis_revers[0][2] * this->nodes[i]->z;
				self_nodes[i].y = this->self_basis_revers[1][0] * this->nodes[i]->x + this->self_basis_revers[1][1] * this->nodes[i]->y + this->self_basis_revers[1][2] * this->nodes[i]->z;
				self_nodes[i].z = this->self_basis_revers[2][0] * this->nodes[i]->x + this->self_basis_revers[2][1] * this->nodes[i]->y + this->self_basis_revers[2][2] * this->nodes[i]->z;
			}
		};
		std::vector<std::vector<double>>* GetSelfReverseBasis()
		{
			if (this->self_basis_revers.size() == 0)
			{
				this->CreateReversBasis();
			}
			return &(this->self_basis_revers);
		}
		std::vector<std::vector<double>>* GetSelfBasis()
		{
			if (this->self_basis.size() == 0)
			{
				this->CreateBasis();
			}
			return &(this->self_basis);
		}
		Point<double> MakeTransferSelfIntoXyz(Point<double> X)
		{
			Point<double> tmp;
			if (this->self_basis.size() == 0)
			{
				CreateBasis();
			}

			tmp.x = this->self_basis[0][0] * X.x + this->self_basis[0][1] * X.y + this->self_basis[0][2] * X.z;
			tmp.y = this->self_basis[1][0] * X.x + this->self_basis[1][1] * X.y + this->self_basis[1][2] * X.z;
			tmp.z = this->self_basis[2][0] * X.x + this->self_basis[2][1] * X.y + this->self_basis[2][2] * X.z;
			return tmp;
		};
		Point<double> MakeTransferXyzIntoSelf(Point<double> X)
		{
			if (this->self_basis.size() == 0)
			{
				CreateBasis();
			}
			if (this->self_basis_revers.size() == 0)
				CreateReversBasis();

			Point<double> tmp;
			tmp.x = this->self_basis_revers[0][0] * X.x + this->self_basis_revers[0][1] * X.y + this->self_basis_revers[0][2] * X.z;
			tmp.y = this->self_basis_revers[1][0] * X.x + this->self_basis_revers[1][1] * X.y + this->self_basis_revers[1][2] * X.z;
			tmp.z = this->self_basis_revers[2][0] * X.x + this->self_basis_revers[2][1] * X.y + this->self_basis_revers[2][2] * X.z;
			return tmp;
		};
		Point<double> GetSelfNode(int id_local)
		{
			if (self_nodes.size() == 0)
			{
				this->self_nodes.resize(this->GetNodesCount());
				for (int i = 0; i < this->self_nodes.size(); i++)
				{
					this->self_nodes[i] = this->MakeTransferXyzIntoSelf(this->GetNode(i));
				}
			}
			return self_nodes[id_local];
		}
		
		void SetIntegrationLaw(int points_count)
		{
			this->integration_parameters.w.resize(points_count);
			this->integration_parameters.x.resize(points_count);

			this->CreateBasis();

			switch (points_count)
			{
			case 1:
			{
				this->integration_parameters.w[0] = 1. / 2.; this->integration_parameters.x[0] = Point<double>(1. / 3., 1. / 3., 0);
			}
			break;
			case 3:
			{
				this->integration_parameters.w[0] = 1. / 6.; this->integration_parameters.x[0] = Point<double>(1. / 6., 1. / 6., 0);
				this->integration_parameters.w[1] = 1. / 6.; this->integration_parameters.x[1] = Point<double>(2. / 3., 1. / 6., 0);
				this->integration_parameters.w[2] = 1. / 6.; this->integration_parameters.x[2] = Point<double>(2. / 3., 1. / 6., 0);
			}
			break;
			case 4:
			{
				this->integration_parameters.w[0] = -9. / 32.; this->integration_parameters.x[0] = Point<double>(1. / 3., 1. / 3., 0);
				this->integration_parameters.w[1] = 25. / 96.; this->integration_parameters.x[1] = Point<double>(3. / 5., 1. / 5., 0);
				this->integration_parameters.w[2] = 25. / 96.; this->integration_parameters.x[2] = Point<double>(1. / 5., 3. / 5., 0);
				this->integration_parameters.w[3] = 25. / 96.; this->integration_parameters.x[3] = Point<double>(1. / 5., 1. / 5., 0);
			}
			break;
			default:
			{
				this->integration_parameters.w[0] = -9. / 32.; this->integration_parameters.x[0] = Point<double>(1. / 3., 1. / 3., 0);
				this->integration_parameters.w[1] = 25. / 96.; this->integration_parameters.x[1] = Point<double>(3. / 5., 1. / 5., 0);
				this->integration_parameters.w[2] = 25. / 96.; this->integration_parameters.x[2] = Point<double>(1. / 5., 3. / 5., 0);
				this->integration_parameters.w[3] = 25. / 96.; this->integration_parameters.x[3] = Point<double>(1. / 5., 1. / 5., 0);
			}
			break;
			}

			double detD = abs(this->self_nodes[0].x*this->self_nodes[1].y + this->self_nodes[1].x*this->self_nodes[2].y + this->self_nodes[2].x*this->self_nodes[0].y -
				this->self_nodes[2].x*this->self_nodes[1].y - this->self_nodes[1].x*this->self_nodes[0].y - this->self_nodes[0].x*this->self_nodes[2].y);
			//this->get_alpha();

			for (int i = 0; i < this->integration_parameters.w.size(); i++)
			{
				this->integration_parameters.w[i] *= detD;

				Point<double> Xtr;
				Xtr.x = (this->self_nodes[1].x - this->self_nodes[0].x)*this->integration_parameters.x[i].x +
					(this->self_nodes[2].x - this->self_nodes[0].x)*this->integration_parameters.x[i].y + this->self_nodes[0].x;
				Xtr.y = (this->self_nodes[1].y - this->self_nodes[0].y)*this->integration_parameters.x[i].x +
					(this->self_nodes[2].y - this->self_nodes[0].y)*this->integration_parameters.x[i].y + this->self_nodes[0].y;
				Xtr.z = this->self_nodes[0].z;

				this->integration_parameters.x[i] = Xtr;
				this->integration_parameters.x[i] = this->MakeTransferSelfIntoXyz(Xtr);
			}
		}

		bool IsContainThePoint(Point<double> A)
		{
			try
			{
				Point<double> P[3];
				P[0] = this->GetNode(0);
				P[1] = this->GetNode(1);
				P[2] = this->GetNode(2);
				double V = GetVolume(), V1, V2, V3;

				P[0] = this->GetNode(0);
				P[1] = this->GetNode(1);
				P[2] = A;
				V1 = geometry::Triangle(P[0], P[1], P[2]).GetVolume();
				P[0] = this->GetNode(0);
				P[1] = A;
				P[2] = this->GetNode(2);
				V2 = geometry::Triangle(P[0], P[1], P[2]).GetVolume();
				P[0] = A;
				P[1] = this->GetNode(1);
				P[2] = this->GetNode(2);
				V3 = geometry::Triangle(P[0], P[1], P[2]).GetVolume();

				double summ = V1 + V2 + V3;
				double res = math::GetRound(abs(summ - V) / V, 5);
				if (res <= 1E-4)
					return true;
				return false;
			}
			catch (const std::exception&)
			{
				printf_s("Error: GeometryShape.cpp/IsContainThePoint(Point<double> A)");
				return false;
			}
		}
		bool IsContainThePoint(Point<double> A, double &length)
		{
			try
			{
				Point<double> P[3];
				P[0] = this->GetNode(0);
				P[1] = this->GetNode(1);
				P[2] = this->GetNode(2);
				double V = GetVolume(), V1, V2, V3;

				P[0] = this->GetNode(0);
				P[1] = this->GetNode(1);
				P[2] = A;
				V1 = geometry::Triangle(P[0], P[1], P[2]).GetVolume();
				P[0] = this->GetNode(0);
				P[1] = A;
				P[2] = this->GetNode(2);
				V2 = geometry::Triangle(P[0], P[1], P[2]).GetVolume();
				P[0] = A;
				P[1] = this->GetNode(1);
				P[2] = this->GetNode(2);
				V3 = geometry::Triangle(P[0], P[1], P[2]).GetVolume();

				double summ = V1 + V2 + V3;
				double res = math::GetRound(abs(summ - V) / V, 5);
				length = res; //!! it is stange!!
				if (res <= 1E-4)
					return true;
				return false;
			}
			catch (const std::exception&)
			{
				printf_s("Error: GeometryShape.cpp/IsContainThePoint(Point<double> A)");
				return false;
			}
		}

		void SolveAlphaMatrix()
		{
			if (alpha.size() == 0)
			{
				alpha.resize(3);
				alpha[0].resize(3);
				alpha[1].resize(3);
				alpha[2].resize(3);
			}
			std::vector<Point<double>> node_self(3);
			node_self[0] = this->GetSelfNode(0);
			node_self[1] = this->GetSelfNode(1);
			node_self[2] = this->GetSelfNode(2);

			double al[3][3] = { {(node_self[1].x * node_self[2].y - node_self[2].x * node_self[1].y) / (node_self[0].y * node_self[2].x - node_self[0].y * node_self[1].x - node_self[0].x * node_self[2].y + node_self[0].x * node_self[1].y + node_self[1].x * node_self[2].y - node_self[2].x * node_self[1].y), (-node_self[2].y + node_self[1].y) / (node_self[0].y * node_self[2].x - node_self[0].y * node_self[1].x - node_self[0].x * node_self[2].y + node_self[0].x * node_self[1].y + node_self[1].x * node_self[2].y - node_self[2].x * node_self[1].y), -(-node_self[2].x + node_self[1].x) / (node_self[0].y * node_self[2].x - node_self[0].y * node_self[1].x - node_self[0].x * node_self[2].y + node_self[0].x * node_self[1].y + node_self[1].x * node_self[2].y - node_self[2].x * node_self[1].y)},
						   {-(node_self[0].x * node_self[2].y - node_self[0].y * node_self[2].x) / (node_self[0].y * node_self[2].x - node_self[0].y * node_self[1].x - node_self[0].x * node_self[2].y + node_self[0].x * node_self[1].y + node_self[1].x * node_self[2].y - node_self[2].x * node_self[1].y), -(node_self[0].y - node_self[2].y) / (node_self[0].y * node_self[2].x - node_self[0].y * node_self[1].x - node_self[0].x * node_self[2].y + node_self[0].x * node_self[1].y + node_self[1].x * node_self[2].y - node_self[2].x * node_self[1].y), (node_self[0].x - node_self[2].x) / (node_self[0].y * node_self[2].x - node_self[0].y * node_self[1].x - node_self[0].x * node_self[2].y + node_self[0].x * node_self[1].y + node_self[1].x * node_self[2].y - node_self[2].x * node_self[1].y)},
						   {(node_self[0].x * node_self[1].y - node_self[0].y * node_self[1].x) / (node_self[0].y * node_self[2].x - node_self[0].y * node_self[1].x - node_self[0].x * node_self[2].y + node_self[0].x * node_self[1].y + node_self[1].x * node_self[2].y - node_self[2].x * node_self[1].y), (node_self[0].y - node_self[1].y) / (node_self[0].y * node_self[2].x - node_self[0].y * node_self[1].x - node_self[0].x * node_self[2].y + node_self[0].x * node_self[1].y + node_self[1].x * node_self[2].y - node_self[2].x * node_self[1].y), -(node_self[0].x - node_self[1].x) / (node_self[0].y * node_self[2].x - node_self[0].y * node_self[1].x - node_self[0].x * node_self[2].y + node_self[0].x * node_self[1].y + node_self[1].x * node_self[2].y - node_self[2].x * node_self[1].y)} };
			
			for (int i = 0; i < 3; i++)
				for (int j = 0; j < 3; j++)
					alpha[i][j] = al[i][j];
		}
	};
	class Tetrahedron: public Shape
	{
	public:
		Tetrahedron()
		{
			this->id_nodes.resize(4);
			this->nodes.resize(4);
		};
		std::vector<std::vector<double>> alpha;

		int GetEdgesCount()
		{
			return 6;
		}

		Tetrahedron(std::vector<Point<double> *> &P)
		{
			try
			{
				this->UpdateNodes(P);
			}
			catch (const std::exception&)
			{
				printf_s("Error: GeometryShape.cpp/Tetrahedron(std::vector<Point<double> *> &P)");
			}
		}
		Tetrahedron(std::vector<int>& id_nodes, std::vector<Point<double>*>& xyz_nodes)
		{
			try
			{
				math::MakeCopyVector_A_into_B(id_nodes, this->id_nodes);
				math::MakeCopyVector_A_into_B(xyz_nodes, this->nodes);
			}
			catch (const std::exception&)
			{
				printf_s("Error: GeometryShape.cpp/geometry::Shape(std::vector<int> &id_nodes, std::vector<Point<double>*> &xyz_nodes)\n");
			}
		};
		double GetVolume()
		{
			try
			{
				double detD_1;
				Point<double> _nodes[4];
				for (int i = 0; i < 4; i++)
				{
					_nodes[i] = this->GetNode(i);
				}
				detD_1 = _nodes[1].z*_nodes[2].y*_nodes[0].x + _nodes[3].z*_nodes[0].x*_nodes[1].y + _nodes[0].x*_nodes[2].z*_nodes[3].y - _nodes[2].z*_nodes[0].x*_nodes[1].y - _nodes[0].x*_nodes[3].z*_nodes[2].y - _nodes[1].z*_nodes[3].y*_nodes[0].x
					- _nodes[3].z*_nodes[1].x*_nodes[0].y + _nodes[0].z*_nodes[3].y*_nodes[1].x - _nodes[0].z*_nodes[2].y*_nodes[1].x + _nodes[2].z*_nodes[1].x*_nodes[0].y + _nodes[3].z*_nodes[1].x*_nodes[2].y - _nodes[3].z*_nodes[1].y*_nodes[2].x + _nodes[3].z*_nodes[0].y*_nodes[2].x
					- _nodes[1].z*_nodes[3].x*_nodes[2].y - _nodes[1].z*_nodes[0].y*_nodes[2].x + _nodes[0].z*_nodes[3].x*_nodes[2].y - _nodes[0].z*_nodes[3].x*_nodes[1].y - _nodes[2].z*_nodes[1].x*_nodes[3].y + _nodes[2].z*_nodes[1].y*_nodes[3].x - _nodes[2].z*_nodes[0].y*_nodes[3].x
					+ _nodes[1].z*_nodes[2].x*_nodes[3].y + _nodes[1].z*_nodes[0].y*_nodes[3].x - _nodes[0].z*_nodes[2].x*_nodes[3].y + _nodes[0].z*_nodes[2].x*_nodes[1].y;
				detD_1 = abs(detD_1) / 6.;
				/*double detD_2;
				Point<double> P21 = *_nodes[1] - *_nodes[0];
				Point<double> P31 = *_nodes[2] - *_nodes[0];
				Point<double> P41 = *_nodes[3] - *_nodes[0];
				detD_2 = P21.x*(P31.y*P41.z - P41.y*P31.z)
				+ P21.y*(P31.z*P41.x - P41.z*P31.x)
				+ P21.z*(P31.x*P41.y - P41.x*P31.y);
				detD_2 = abs(detD_2) / 6.;*/

				return detD_1;
			}
			catch (const std::exception&)
			{
				printf_s("Error: GeometryShape.cpp/GetVolume()");
				return 0.0;
			}
		}
		double GetVolume(Point<double>* _P)
		{
			try
			{
				double detD_1 = _P[1].z*_P[2].y*_P[0].x + _P[3].z*_P[0].x*_P[1].y + _P[0].x*_P[2].z*_P[3].y - _P[2].z*_P[0].x*_P[1].y - _P[0].x*_P[3].z*_P[2].y - _P[1].z*_P[3].y*_P[0].x
					- _P[3].z*_P[1].x*_P[0].y + _P[0].z*_P[3].y*_P[1].x - _P[0].z*_P[2].y*_P[1].x + _P[2].z*_P[1].x*_P[0].y + _P[3].z*_P[1].x*_P[2].y - _P[3].z*_P[1].y*_P[2].x + _P[3].z*_P[0].y*_P[2].x
					- _P[1].z*_P[3].x*_P[2].y - _P[1].z*_P[0].y*_P[2].x + _P[0].z*_P[3].x*_P[2].y - _P[0].z*_P[3].x*_P[1].y - _P[2].z*_P[1].x*_P[3].y + _P[2].z*_P[1].y*_P[3].x - _P[2].z*_P[0].y*_P[3].x
					+ _P[1].z*_P[2].x*_P[3].y + _P[1].z*_P[0].y*_P[3].x - _P[0].z*_P[2].x*_P[3].y + _P[0].z*_P[2].x*_P[1].y;
				detD_1 = abs(detD_1) / 6.;

				return detD_1;
			}
			catch (const std::exception&)
			{
				printf_s("Error: GeometryShape.cpp/GetVolume(Point[4] _P)");
			}
		}
		static double GetVolume(Point<double> &_P0, Point<double>& _P1, Point<double>& _P2, Point<double>& _P3)
		{
			try
			{
				double detD_1 = _P1.z * _P2.y * _P0.x + _P3.z * _P0.x * _P1.y + _P0.x * _P2.z * _P3.y - _P2.z * _P0.x * _P1.y - _P0.x * _P3.z * _P2.y - _P1.z * _P3.y * _P0.x
					- _P3.z * _P1.x * _P0.y + _P0.z * _P3.y * _P1.x - _P0.z * _P2.y * _P1.x + _P2.z * _P1.x * _P0.y + _P3.z * _P1.x * _P2.y - _P3.z * _P1.y * _P2.x + _P3.z * _P0.y * _P2.x
					- _P1.z * _P3.x * _P2.y - _P1.z * _P0.y * _P2.x + _P0.z * _P3.x * _P2.y - _P0.z * _P3.x * _P1.y - _P2.z * _P1.x * _P3.y + _P2.z * _P1.y * _P3.x - _P2.z * _P0.y * _P3.x
					+ _P1.z * _P2.x * _P3.y + _P1.z * _P0.y * _P3.x - _P0.z * _P2.x * _P3.y + _P0.z * _P2.x * _P1.y;
				detD_1 = abs(detD_1) / 6.;

				return detD_1;
			}
			catch (const std::exception&)
			{
				printf_s("Error: GeometryShape.cpp/GetVolume(Point[4] _P)");
			}
		}

		void SetGeometry(std::vector<int>& id_global_nodes, std::vector<Point<double>*>& P)
		{
			try
			{
				this->SetIdNodes(id_global_nodes);
				this->UpdateNodes(P);
			}
			catch (const std::exception&)
			{
				printf_s("Error: GeometryShape.h/Tetrahedron::SetGeometry(std::vector<std::vector<int>>& id_nodes, std::vector<Point<double>*>& P)");
			}
		};

		void SetIntegrationLaw(int points_count)
		{
			this->integration_parameters.w.resize(points_count);
			this->integration_parameters.x.resize(points_count);

			switch (points_count)
			{
			case 1:
			{
				this->integration_parameters.w[0] = 1. / 6.; this->integration_parameters.x[0] = Point<double>(1. / 4., 1. / 4., 1. / 4.);
			}
			break;
			case 4:
			{
				double u = (5. - sqrt(5.)) / 20., v = 1 - 3 * u;
				this->integration_parameters.w[0] = 1. / 24.; this->integration_parameters.x[0] = Point<double>(u, u, u);
				this->integration_parameters.w[1] = 1. / 24.; this->integration_parameters.x[1] = Point<double>(v, u, u);
				this->integration_parameters.w[2] = 1. / 24.; this->integration_parameters.x[2] = Point<double>(u, v, u);
				this->integration_parameters.w[3] = 1. / 24.; this->integration_parameters.x[3] = Point<double>(u, u, v);
			}
			break;
			case 5:
			{
				double u = 1. / 6., v = 1 - 3 * u;
				this->integration_parameters.w[0] = 36. / 480.; this->integration_parameters.x[0] = Point<double>(u, u, u);
				this->integration_parameters.w[1] = 36. / 480.; this->integration_parameters.x[1] = Point<double>(v, u, u);
				this->integration_parameters.w[2] = 36. / 480.; this->integration_parameters.x[2] = Point<double>(u, v, u);
				this->integration_parameters.w[3] = 36. / 480.; this->integration_parameters.x[3] = Point<double>(u, u, v);
				this->integration_parameters.w[4] = -64. / 480.; this->integration_parameters.x[4] = Point<double>(1. / 4., 1. / 4., 1. / 4.);
			}
			break;
			case 11:
			{
				double u1 = 1. / 14., v1 = 1 - 3 * u1;
				double u2 = (14. + sqrt(70.)) / 56., v2 = 1 - 3 * u2;
				this->integration_parameters.w[0] = 686. / 90000.; this->integration_parameters.x[0] = Point<double>(u1, u1, u1);
				this->integration_parameters.w[1] = 686. / 90000.; this->integration_parameters.x[1] = Point<double>(v1, u1, u1);
				this->integration_parameters.w[2] = 686. / 90000.; this->integration_parameters.x[2] = Point<double>(u1, v1, u1);
				this->integration_parameters.w[3] = 686. / 90000.; this->integration_parameters.x[3] = Point<double>(u1, u1, v1);
				this->integration_parameters.w[4] = -1184. / 90000.; this->integration_parameters.x[4] = Point<double>(1. / 4., 1. / 4., 1. / 4.);
				this->integration_parameters.w[5] = 448. / 18000.; this->integration_parameters.x[5] = Point<double>((u2 + v2) / 2., u2, u2);
				this->integration_parameters.w[6] = 448. / 18000.; this->integration_parameters.x[6] = Point<double>(u2, (u2 + v2) / 2., u2);
				this->integration_parameters.w[7] = 448. / 18000.; this->integration_parameters.x[7] = Point<double>(u2, u2, (u2 + v2) / 2.);
				this->integration_parameters.w[8] = 448. / 18000.; this->integration_parameters.x[8] = Point<double>((u2 + v2) / 2., (u2 + v2) / 2., u2);
				this->integration_parameters.w[9] = 448. / 18000.; this->integration_parameters.x[9] = Point<double>((u2 + v2) / 2., u2, (u2 + v2) / 2.);
				this->integration_parameters.w[10] = 448. / 18000.; this->integration_parameters.x[10] = Point<double>(u2, (u2 + v2) / 2., (u2 + v2) / 2.);
			}
			break;
			case 14:
			{
				double u1 = 0.31088592, v1 = 1 - 3 * u1;
				double u2 = 0.09273525, v2 = 1 - 3 * u2;
				double u3 = 0.45449630, v3 = 1 - 3 * u3;
				this->integration_parameters.w[0] = 0.01878132; this->integration_parameters.x[0] = Point<double>(u1, u1, u1);
				this->integration_parameters.w[1] = 0.01878132; this->integration_parameters.x[1] = Point<double>(v1, u1, u1);
				this->integration_parameters.w[2] = 0.01878132; this->integration_parameters.x[2] = Point<double>(u1, v1, u1);
				this->integration_parameters.w[3] = 0.01878132; this->integration_parameters.x[3] = Point<double>(u1, u1, v1);
				this->integration_parameters.w[4] = 0.01224884; this->integration_parameters.x[4] = Point<double>(u2, u2, u2);
				this->integration_parameters.w[5] = 0.01224884; this->integration_parameters.x[5] = Point<double>(v2, u2, u2);
				this->integration_parameters.w[6] = 0.01224884; this->integration_parameters.x[6] = Point<double>(u2, v2, u2);
				this->integration_parameters.w[7] = 0.01224884; this->integration_parameters.x[7] = Point<double>(u2, u2, v2);
				this->integration_parameters.w[8] = 0.00709100; this->integration_parameters.x[8] = Point<double>((u3 + v3) / 2., u3, u3);
				this->integration_parameters.w[9] = 0.00709100; this->integration_parameters.x[9] = Point<double>(u3, (u3 + v3) / 2., u3);
				this->integration_parameters.w[10] = 0.00709100; this->integration_parameters.x[10] = Point<double>(u3, u3, (u3 + v3) / 2.);
				this->integration_parameters.w[11] = 0.00709100; this->integration_parameters.x[11] = Point<double>((u3 + v3) / 2., (u3 + v3) / 2., u3);
				this->integration_parameters.w[12] = 0.00709100; this->integration_parameters.x[12] = Point<double>((u3 + v3) / 2., u3, (u3 + v3) / 2.);
				this->integration_parameters.w[13] = 0.00709100; this->integration_parameters.x[13] = Point<double>(u3, (u3 + v3) / 2., (u3 + v3) / 2.);
			}
			break;
			default: //4
			{
				double u = (5. - sqrt(5.)) / 20., v = 1 - 3 * u;
				this->integration_parameters.w[0] = 1. / 24.; this->integration_parameters.x[0] = Point<double>(u, u, u);
				this->integration_parameters.w[1] = 1. / 24.; this->integration_parameters.x[1] = Point<double>(v, u, u);
				this->integration_parameters.w[2] = 1. / 24.; this->integration_parameters.x[2] = Point<double>(u, v, u);
				this->integration_parameters.w[3] = 1. / 24.; this->integration_parameters.x[3] = Point<double>(u, u, v);
			}
			break;
			}

			for (int i = 0; i < this->integration_parameters.w.size(); i++)
			{
				this->integration_parameters.w[i] *= 6 * this->GetVolume();

				Point<double> xtr;
				xtr.x = (this->nodes[1]->x - this->nodes[0]->x)*this->integration_parameters.x[i].x +
					(this->nodes[2]->x - this->nodes[0]->x)*this->integration_parameters.x[i].y +
					(this->nodes[3]->x - this->nodes[0]->x)*this->integration_parameters.x[i].z + this->nodes[0]->x;
				xtr.y = (this->nodes[1]->y - this->nodes[0]->y)*this->integration_parameters.x[i].x +
					(this->nodes[2]->y - this->nodes[0]->y)*this->integration_parameters.x[i].y +
					(this->nodes[3]->y - this->nodes[0]->y)*this->integration_parameters.x[i].z + this->nodes[0]->y;
				xtr.z = (this->nodes[1]->z - this->nodes[0]->z)*this->integration_parameters.x[i].x +
					(this->nodes[2]->z - this->nodes[0]->z)*this->integration_parameters.x[i].y +
					(this->nodes[3]->z - this->nodes[0]->z)*this->integration_parameters.x[i].z + this->nodes[0]->z;

				this->integration_parameters.x[i] = xtr;
			}
		}
		void GetIntegrationProperties(std::vector<double> &W, std::vector<Point<double>> &X)
		{
			W.resize(this->GetIntegrationPointsCount());
			X.resize(this->GetIntegrationPointsCount());
			for (int i = 0; i < W.size(); i++)
			{
				W[i] = this->integration_parameters.w[i];
				X[i] = this->integration_parameters.x[i];
			}
		}

		bool IsContainThePoint(Point<double> A)
		{
			try
			{
				Point<double> P[4];
				double V = GetVolume(), V1, V2, V3, V4;

				P[0] = this->GetNode(0);
				P[1] = this->GetNode(1);
				P[2] = this->GetNode(2);
				P[3] = A;
				V1 = GetVolume(P);
				P[0] = this->GetNode(0);
				P[1] = this->GetNode(1);
				P[2] = A;
				P[3] = this->GetNode(3);
				V2 = GetVolume(P);
				P[0] = this->GetNode(0);
				P[1] = A;
				P[2] = this->GetNode(2);
				P[3] = this->GetNode(3);
				V3 = GetVolume(P);
				P[0] = A;
				P[1] = this->GetNode(1);
				P[2] = this->GetNode(2);
				P[3] = this->GetNode(3);
				V4 = GetVolume(P);

				double summ = V1 + V2 + V3 + V4;
				double res = math::GetRound(abs(summ - V) / V, 5);
				if (res <= 1E-15)
					return true;
				return false;
			}
			catch (const std::exception&)
			{
				printf_s("Error: GeometryShape.cpp/IsContainThePoint(Point<double> A)");
				return false;
			}
		}
		bool IsContainThePoint(Point<double> A, double & discrepancy)
		{
			try
			{
				/*Point<double> P[4];
				double V = GetVolume(), V1, V2, V3, V4;

				P[0] = this->GetNode(0);
				P[1] = this->GetNode(1);
				P[2] = this->GetNode(2);
				P[3] = A;
				V1 = GetVolume(P);
				P[0] = this->GetNode(0);
				P[1] = this->GetNode(1);
				P[2] = A;
				P[3] = this->GetNode(3);
				V2 = GetVolume(P);
				P[0] = this->GetNode(0);
				P[1] = A;
				P[2] = this->GetNode(2);
				P[3] = this->GetNode(3);
				V3 = GetVolume(P);
				P[0] = A;
				P[1] = this->GetNode(1);
				P[2] = this->GetNode(2);
				P[3] = this->GetNode(3);
				V4 = GetVolume(P);

				double summ = V1 + V2 + V3 + V4;
				double res = math::GetRound(abs(summ - V) / V, 5);
				length = res;
				if (res <= 1E-15)
					return true;
				return false;*/

				double L1 = alpha[0][0] + alpha[0][1] * A.x + alpha[0][2] * A.y + alpha[0][3] * A.z;
				double L2 = alpha[1][0] + alpha[1][1] * A.x + alpha[1][2] * A.y + alpha[1][3] * A.z;
				double L3 = alpha[2][0] + alpha[2][1] * A.x + alpha[2][2] * A.y + alpha[2][3] * A.z;
				double L4 = alpha[3][0] + alpha[3][1] * A.x + alpha[3][2] * A.y + alpha[3][3] * A.z;

				if (abs(L1) < 1e-15 && abs(L1 - 1) < 1e-15 &&
					abs(L2) < 1e-15 && abs(L2 - 1) < 1e-15 &&
					abs(L3) < 1e-15 && abs(L3 - 1) < 1e-15 &&
					abs(L4) < 1e-15 && abs(L4 - 1) < 1e-15)
				{
					discrepancy = 0.0;
					return true;
				}
				discrepancy = fmin(abs(L1), abs(L1 - 1)) + fmin(abs(L2), abs(L2 - 1)) + fmin(abs(L3), abs(L3 - 1)) + fmin(abs(L4), abs(L4 - 1));
			}
			catch (const std::exception&)
			{
				printf_s("Error: GeometryShape.cpp/IsContainThePoint(Point<double> A)");
				return false;
			}
		}

		void SolveAlphaMatrix()
		{
			if (alpha.size() != 4) {
				alpha.resize(4);
				alpha[0].resize(4);
				alpha[1].resize(4);
				alpha[2].resize(4);
				alpha[3].resize(4);
			}

			std::vector<Point<double>> node_xyz(4);
			node_xyz[0] = this->GetNode(0);
			node_xyz[1] = this->GetNode(1);
			node_xyz[2] = this->GetNode(2);
			node_xyz[3] = this->GetNode(3);

			alpha[0][0] = -(node_xyz[3].z*node_xyz[1].x*node_xyz[2].y - node_xyz[3].z*node_xyz[1].y*node_xyz[2].x - node_xyz[1].z*node_xyz[3].x*node_xyz[2].y - node_xyz[2].z*node_xyz[1].x*node_xyz[3].y
				+ node_xyz[2].z*node_xyz[1].y*node_xyz[3].x + node_xyz[1].z*node_xyz[2].x*node_xyz[3].y) / (-node_xyz[3].z*node_xyz[0].x*node_xyz[1].y + node_xyz[0].x*node_xyz[3].z*node_xyz[2].y
					- node_xyz[0].x*node_xyz[2].z*node_xyz[3].y + node_xyz[2].z*node_xyz[0].x*node_xyz[1].y - node_xyz[1].z*node_xyz[2].y*node_xyz[0].x + node_xyz[1].z*node_xyz[3].y*node_xyz[0].x
					+ node_xyz[3].z*node_xyz[1].x*node_xyz[0].y - node_xyz[0].z*node_xyz[3].y*node_xyz[1].x + node_xyz[0].z*node_xyz[2].y*node_xyz[1].x - node_xyz[2].z*node_xyz[1].x*node_xyz[0].y
					- node_xyz[3].z*node_xyz[1].x*node_xyz[2].y + node_xyz[3].z*node_xyz[1].y*node_xyz[2].x - node_xyz[3].z*node_xyz[0].y*node_xyz[2].x + node_xyz[1].z*node_xyz[3].x*node_xyz[2].y
					+ node_xyz[1].z*node_xyz[0].y*node_xyz[2].x - node_xyz[0].z*node_xyz[3].x*node_xyz[2].y + node_xyz[0].z*node_xyz[3].x*node_xyz[1].y + node_xyz[2].z*node_xyz[1].x*node_xyz[3].y
					- node_xyz[2].z*node_xyz[1].y*node_xyz[3].x + node_xyz[2].z*node_xyz[0].y*node_xyz[3].x - node_xyz[1].z*node_xyz[2].x*node_xyz[3].y - node_xyz[1].z*node_xyz[0].y*node_xyz[3].x
					+ node_xyz[0].z*node_xyz[2].x*node_xyz[3].y - node_xyz[0].z*node_xyz[2].x*node_xyz[1].y);
			alpha[0][1] = (-node_xyz[3].z*node_xyz[1].y + node_xyz[3].z*node_xyz[2].y - node_xyz[2].z*node_xyz[3].y + node_xyz[2].z*node_xyz[1].y - node_xyz[1].z*node_xyz[2].y + node_xyz[1].z*node_xyz[3].y) /
				(-node_xyz[3].z*node_xyz[0].x*node_xyz[1].y + node_xyz[0].x*node_xyz[3].z*node_xyz[2].y - node_xyz[0].x*node_xyz[2].z*node_xyz[3].y + node_xyz[2].z*node_xyz[0].x*node_xyz[1].y
					- node_xyz[1].z*node_xyz[2].y*node_xyz[0].x + node_xyz[1].z*node_xyz[3].y*node_xyz[0].x + node_xyz[3].z*node_xyz[1].x*node_xyz[0].y - node_xyz[0].z*node_xyz[3].y*node_xyz[1].x
					+ node_xyz[0].z*node_xyz[2].y*node_xyz[1].x - node_xyz[2].z*node_xyz[1].x*node_xyz[0].y - node_xyz[3].z*node_xyz[1].x*node_xyz[2].y + node_xyz[3].z*node_xyz[1].y*node_xyz[2].x
					- node_xyz[3].z*node_xyz[0].y*node_xyz[2].x + node_xyz[1].z*node_xyz[3].x*node_xyz[2].y + node_xyz[1].z*node_xyz[0].y*node_xyz[2].x - node_xyz[0].z*node_xyz[3].x*node_xyz[2].y
					+ node_xyz[0].z*node_xyz[3].x*node_xyz[1].y + node_xyz[2].z*node_xyz[1].x*node_xyz[3].y - node_xyz[2].z*node_xyz[1].y*node_xyz[3].x + node_xyz[2].z*node_xyz[0].y*node_xyz[3].x
					- node_xyz[1].z*node_xyz[2].x*node_xyz[3].y - node_xyz[1].z*node_xyz[0].y*node_xyz[3].x + node_xyz[0].z*node_xyz[2].x*node_xyz[3].y - node_xyz[0].z*node_xyz[2].x*node_xyz[1].y);
			alpha[0][2] = -(-node_xyz[3].z*node_xyz[1].x + node_xyz[2].z*node_xyz[1].x + node_xyz[3].z*node_xyz[2].x - node_xyz[1].z*node_xyz[2].x - node_xyz[2].z*node_xyz[3].x + node_xyz[1].z*node_xyz[3].x) /
				(-node_xyz[3].z*node_xyz[0].x*node_xyz[1].y + node_xyz[0].x*node_xyz[3].z*node_xyz[2].y - node_xyz[0].x*node_xyz[2].z*node_xyz[3].y + node_xyz[2].z*node_xyz[0].x*node_xyz[1].y
					- node_xyz[1].z*node_xyz[2].y*node_xyz[0].x + node_xyz[1].z*node_xyz[3].y*node_xyz[0].x + node_xyz[3].z*node_xyz[1].x*node_xyz[0].y - node_xyz[0].z*node_xyz[3].y*node_xyz[1].x
					+ node_xyz[0].z*node_xyz[2].y*node_xyz[1].x - node_xyz[2].z*node_xyz[1].x*node_xyz[0].y - node_xyz[3].z*node_xyz[1].x*node_xyz[2].y + node_xyz[3].z*node_xyz[1].y*node_xyz[2].x
					- node_xyz[3].z*node_xyz[0].y*node_xyz[2].x + node_xyz[1].z*node_xyz[3].x*node_xyz[2].y + node_xyz[1].z*node_xyz[0].y*node_xyz[2].x - node_xyz[0].z*node_xyz[3].x*node_xyz[2].y
					+ node_xyz[0].z*node_xyz[3].x*node_xyz[1].y + node_xyz[2].z*node_xyz[1].x*node_xyz[3].y - node_xyz[2].z*node_xyz[1].y*node_xyz[3].x + node_xyz[2].z*node_xyz[0].y*node_xyz[3].x
					- node_xyz[1].z*node_xyz[2].x*node_xyz[3].y - node_xyz[1].z*node_xyz[0].y*node_xyz[3].x + node_xyz[0].z*node_xyz[2].x*node_xyz[3].y - node_xyz[0].z*node_xyz[2].x*node_xyz[1].y);
			alpha[0][3] = (-node_xyz[3].y*node_xyz[1].x + node_xyz[2].y*node_xyz[1].x - node_xyz[3].x*node_xyz[2].y + node_xyz[1].y*node_xyz[3].x + node_xyz[2].x*node_xyz[3].y - node_xyz[1].y*node_xyz[2].x) /
				(-node_xyz[3].z*node_xyz[0].x*node_xyz[1].y + node_xyz[0].x*node_xyz[3].z*node_xyz[2].y - node_xyz[0].x*node_xyz[2].z*node_xyz[3].y + node_xyz[2].z*node_xyz[0].x*node_xyz[1].y
					- node_xyz[1].z*node_xyz[2].y*node_xyz[0].x + node_xyz[1].z*node_xyz[3].y*node_xyz[0].x + node_xyz[3].z*node_xyz[1].x*node_xyz[0].y - node_xyz[0].z*node_xyz[3].y*node_xyz[1].x
					+ node_xyz[0].z*node_xyz[2].y*node_xyz[1].x - node_xyz[2].z*node_xyz[1].x*node_xyz[0].y - node_xyz[3].z*node_xyz[1].x*node_xyz[2].y + node_xyz[3].z*node_xyz[1].y*node_xyz[2].x
					- node_xyz[3].z*node_xyz[0].y*node_xyz[2].x + node_xyz[1].z*node_xyz[3].x*node_xyz[2].y + node_xyz[1].z*node_xyz[0].y*node_xyz[2].x - node_xyz[0].z*node_xyz[3].x*node_xyz[2].y
					+ node_xyz[0].z*node_xyz[3].x*node_xyz[1].y + node_xyz[2].z*node_xyz[1].x*node_xyz[3].y - node_xyz[2].z*node_xyz[1].y*node_xyz[3].x + node_xyz[2].z*node_xyz[0].y*node_xyz[3].x
					- node_xyz[1].z*node_xyz[2].x*node_xyz[3].y - node_xyz[1].z*node_xyz[0].y*node_xyz[3].x + node_xyz[0].z*node_xyz[2].x*node_xyz[3].y - node_xyz[0].z*node_xyz[2].x*node_xyz[1].y);
			alpha[1][0] = (node_xyz[0].x*node_xyz[3].z*node_xyz[2].y - node_xyz[0].x*node_xyz[2].z*node_xyz[3].y - node_xyz[3].z*node_xyz[0].y*node_xyz[2].x + node_xyz[2].z*node_xyz[0].y*node_xyz[3].x
				- node_xyz[0].z*node_xyz[3].x*node_xyz[2].y + node_xyz[0].z*node_xyz[2].x*node_xyz[3].y) / (-node_xyz[3].z*node_xyz[0].x*node_xyz[1].y + node_xyz[0].x*node_xyz[3].z*node_xyz[2].y
					- node_xyz[0].x*node_xyz[2].z*node_xyz[3].y + node_xyz[2].z*node_xyz[0].x*node_xyz[1].y - node_xyz[1].z*node_xyz[2].y*node_xyz[0].x + node_xyz[1].z*node_xyz[3].y*node_xyz[0].x
					+ node_xyz[3].z*node_xyz[1].x*node_xyz[0].y - node_xyz[0].z*node_xyz[3].y*node_xyz[1].x + node_xyz[0].z*node_xyz[2].y*node_xyz[1].x - node_xyz[2].z*node_xyz[1].x*node_xyz[0].y
					- node_xyz[3].z*node_xyz[1].x*node_xyz[2].y + node_xyz[3].z*node_xyz[1].y*node_xyz[2].x - node_xyz[3].z*node_xyz[0].y*node_xyz[2].x + node_xyz[1].z*node_xyz[3].x*node_xyz[2].y
					+ node_xyz[1].z*node_xyz[0].y*node_xyz[2].x - node_xyz[0].z*node_xyz[3].x*node_xyz[2].y + node_xyz[0].z*node_xyz[3].x*node_xyz[1].y + node_xyz[2].z*node_xyz[1].x*node_xyz[3].y
					- node_xyz[2].z*node_xyz[1].y*node_xyz[3].x + node_xyz[2].z*node_xyz[0].y*node_xyz[3].x - node_xyz[1].z*node_xyz[2].x*node_xyz[3].y - node_xyz[1].z*node_xyz[0].y*node_xyz[3].x
					+ node_xyz[0].z*node_xyz[2].x*node_xyz[3].y - node_xyz[0].z*node_xyz[2].x*node_xyz[1].y);
			alpha[1][1] = -(-node_xyz[3].z*node_xyz[0].y + node_xyz[0].z*node_xyz[3].y - node_xyz[0].z* node_xyz[2].y + node_xyz[2].z* node_xyz[0].y + node_xyz[3].z* node_xyz[2].y - node_xyz[2].z* node_xyz[3].y) /
				(-node_xyz[3].z*node_xyz[0].x*node_xyz[1].y + node_xyz[0].x* node_xyz[3].z* node_xyz[2].y - node_xyz[0].x* node_xyz[2].z* node_xyz[3].y + node_xyz[2].z* node_xyz[0].x* node_xyz[1].y
					- node_xyz[1].z*node_xyz[2].y*node_xyz[0].x + node_xyz[1].z* node_xyz[3].y* node_xyz[0].x + node_xyz[3].z* node_xyz[1].x* node_xyz[0].y - node_xyz[0].z* node_xyz[3].y* node_xyz[1].x
					+ node_xyz[0].z*node_xyz[2].y*node_xyz[1].x - node_xyz[2].z* node_xyz[1].x* node_xyz[0].y - node_xyz[3].z* node_xyz[1].x* node_xyz[2].y + node_xyz[3].z* node_xyz[1].y* node_xyz[2].x
					- node_xyz[3].z*node_xyz[0].y*node_xyz[2].x + node_xyz[1].z* node_xyz[3].x* node_xyz[2].y + node_xyz[1].z* node_xyz[0].y* node_xyz[2].x - node_xyz[0].z* node_xyz[3].x* node_xyz[2].y
					+ node_xyz[0].z*node_xyz[3].x*node_xyz[1].y + node_xyz[2].z* node_xyz[1].x* node_xyz[3].y - node_xyz[2].z* node_xyz[1].y* node_xyz[3].x + node_xyz[2].z* node_xyz[0].y* node_xyz[3].x
					- node_xyz[1].z*node_xyz[2].x*node_xyz[3].y - node_xyz[1].z* node_xyz[0].y* node_xyz[3].x + node_xyz[0].z* node_xyz[2].x* node_xyz[3].y - node_xyz[0].z* node_xyz[2].x* node_xyz[1].y);
			alpha[1][2] = (-node_xyz[3].z*node_xyz[0].x + node_xyz[2].z *node_xyz[0].x + node_xyz[3].z* node_xyz[2].x + node_xyz[0].z* node_xyz[3].x - node_xyz[2].z* node_xyz[3].x - node_xyz[0].z *node_xyz[2].x) /
				(-node_xyz[3].z*node_xyz[0].x* node_xyz[1].y + node_xyz[0].x* node_xyz[3].z *node_xyz[2].y - node_xyz[0].x* node_xyz[2].z *node_xyz[3].y + node_xyz[2].z *node_xyz[0].x *node_xyz[1].y
					- node_xyz[1].z*node_xyz[2].y* node_xyz[0].x + node_xyz[1].z* node_xyz[3].y *node_xyz[0].x + node_xyz[3].z* node_xyz[1].x *node_xyz[0].y - node_xyz[0].z *node_xyz[3].y *node_xyz[1].x
					+ node_xyz[0].z*node_xyz[2].y* node_xyz[1].x - node_xyz[2].z* node_xyz[1].x *node_xyz[0].y - node_xyz[3].z* node_xyz[1].x *node_xyz[2].y + node_xyz[3].z *node_xyz[1].y *node_xyz[2].x
					- node_xyz[3].z*node_xyz[0].y* node_xyz[2].x + node_xyz[1].z* node_xyz[3].x *node_xyz[2].y + node_xyz[1].z* node_xyz[0].y *node_xyz[2].x - node_xyz[0].z *node_xyz[3].x *node_xyz[2].y
					+ node_xyz[0].z*node_xyz[3].x* node_xyz[1].y + node_xyz[2].z* node_xyz[1].x *node_xyz[3].y - node_xyz[2].z* node_xyz[1].y *node_xyz[3].x + node_xyz[2].z *node_xyz[0].y *node_xyz[3].x
					- node_xyz[1].z*node_xyz[2].x* node_xyz[3].y - node_xyz[1].z* node_xyz[0].y *node_xyz[3].x + node_xyz[0].z* node_xyz[2].x *node_xyz[3].y - node_xyz[0].z *node_xyz[2].x *node_xyz[1].y);
			alpha[1][3] = -(node_xyz[2].y* node_xyz[0].x - node_xyz[3].y* node_xyz[0].x + node_xyz[2].x *node_xyz[3].y + node_xyz[0].y *node_xyz[3].x - node_xyz[3].x *node_xyz[2].y - node_xyz[0].y *node_xyz[2].x) /
				(-node_xyz[3].z* node_xyz[0].x* node_xyz[1].y + node_xyz[0].x *node_xyz[3].z *node_xyz[2].y - node_xyz[0].x* node_xyz[2].z* node_xyz[3].y + node_xyz[2].z* node_xyz[0].x* node_xyz[1].y
					- node_xyz[1].z* node_xyz[2].y* node_xyz[0].x + node_xyz[1].z *node_xyz[3].y *node_xyz[0].x + node_xyz[3].z* node_xyz[1].x* node_xyz[0].y - node_xyz[0].z* node_xyz[3].y* node_xyz[1].x
					+ node_xyz[0].z* node_xyz[2].y* node_xyz[1].x - node_xyz[2].z *node_xyz[1].x *node_xyz[0].y - node_xyz[3].z* node_xyz[1].x* node_xyz[2].y + node_xyz[3].z* node_xyz[1].y* node_xyz[2].x
					- node_xyz[3].z* node_xyz[0].y* node_xyz[2].x + node_xyz[1].z *node_xyz[3].x *node_xyz[2].y + node_xyz[1].z* node_xyz[0].y* node_xyz[2].x - node_xyz[0].z* node_xyz[3].x* node_xyz[2].y
					+ node_xyz[0].z* node_xyz[3].x* node_xyz[1].y + node_xyz[2].z *node_xyz[1].x *node_xyz[3].y - node_xyz[2].z* node_xyz[1].y* node_xyz[3].x + node_xyz[2].z* node_xyz[0].y* node_xyz[3].x
					- node_xyz[1].z* node_xyz[2].x* node_xyz[3].y - node_xyz[1].z *node_xyz[0].y *node_xyz[3].x + node_xyz[0].z* node_xyz[2].x* node_xyz[3].y - node_xyz[0].z* node_xyz[2].x* node_xyz[1].y);
			alpha[2][0] = -(node_xyz[3].z*node_xyz[0].x*node_xyz[1].y - node_xyz[1].z*node_xyz[3].y*node_xyz[0].x - node_xyz[0].z*node_xyz[3].x*node_xyz[1].y - node_xyz[3].z*node_xyz[1].x*node_xyz[0].y + node_xyz[1].z*node_xyz[0].y*node_xyz[3].x + node_xyz[0].z*node_xyz[3].y*node_xyz[1].x) /
				(-node_xyz[3].z*node_xyz[0].x*node_xyz[1].y + node_xyz[0].x*node_xyz[3].z*node_xyz[2].y - node_xyz[0].x*node_xyz[2].z*node_xyz[3].y + node_xyz[2].z*node_xyz[0].x*node_xyz[1].y - node_xyz[1].z*node_xyz[2].y*node_xyz[0].x + node_xyz[1].z*node_xyz[3].y*node_xyz[0].x
					+ node_xyz[3].z*node_xyz[1].x*node_xyz[0].y - node_xyz[0].z*node_xyz[3].y*node_xyz[1].x + node_xyz[0].z*node_xyz[2].y*node_xyz[1].x - node_xyz[2].z*node_xyz[1].x*node_xyz[0].y - node_xyz[3].z*node_xyz[1].x*node_xyz[2].y + node_xyz[3].z*node_xyz[1].y*node_xyz[2].x
					- node_xyz[3].z*node_xyz[0].y*node_xyz[2].x + node_xyz[1].z*node_xyz[3].x*node_xyz[2].y + node_xyz[1].z*node_xyz[0].y*node_xyz[2].x - node_xyz[0].z*node_xyz[3].x*node_xyz[2].y + node_xyz[0].z*node_xyz[3].x*node_xyz[1].y + node_xyz[2].z*node_xyz[1].x*node_xyz[3].y
					- node_xyz[2].z*node_xyz[1].y*node_xyz[3].x + node_xyz[2].z*node_xyz[0].y*node_xyz[3].x - node_xyz[1].z*node_xyz[2].x*node_xyz[3].y - node_xyz[1].z*node_xyz[0].y*node_xyz[3].x + node_xyz[0].z*node_xyz[2].x*node_xyz[3].y - node_xyz[0].z*node_xyz[2].x*node_xyz[1].y);
			alpha[2][1] = (node_xyz[3].z* node_xyz[1].y - node_xyz[3].z* node_xyz[0].y + node_xyz[1].z* node_xyz[0].y - node_xyz[1].z* node_xyz[3].y + node_xyz[0].z *node_xyz[3].y - node_xyz[0].z* node_xyz[1].y) /
				(-node_xyz[3].z* node_xyz[0].x* node_xyz[1].y + node_xyz[0].x* node_xyz[3].z* node_xyz[2].y - node_xyz[0].x* node_xyz[2].z* node_xyz[3].y + node_xyz[2].z* node_xyz[0].x* node_xyz[1].y
-node_xyz[1].z* node_xyz[2].y* node_xyz[0].x + node_xyz[1].z* node_xyz[3].y* node_xyz[0].x + node_xyz[3].z* node_xyz[1].x* node_xyz[0].y - node_xyz[0].z* node_xyz[3].y* node_xyz[1].x
+ node_xyz[0].z* node_xyz[2].y* node_xyz[1].x - node_xyz[2].z* node_xyz[1].x* node_xyz[0].y - node_xyz[3].z* node_xyz[1].x* node_xyz[2].y + node_xyz[3].z* node_xyz[1].y* node_xyz[2].x
- node_xyz[3].z* node_xyz[0].y* node_xyz[2].x + node_xyz[1].z* node_xyz[3].x* node_xyz[2].y + node_xyz[1].z* node_xyz[0].y* node_xyz[2].x - node_xyz[0].z* node_xyz[3].x* node_xyz[2].y
+ node_xyz[0].z* node_xyz[3].x* node_xyz[1].y + node_xyz[2].z* node_xyz[1].x* node_xyz[3].y - node_xyz[2].z* node_xyz[1].y* node_xyz[3].x + node_xyz[2].z* node_xyz[0].y* node_xyz[3].x
- node_xyz[1].z* node_xyz[2].x* node_xyz[3].y - node_xyz[1].z* node_xyz[0].y* node_xyz[3].x + node_xyz[0].z* node_xyz[2].x* node_xyz[3].y - node_xyz[0].z* node_xyz[2].x* node_xyz[1].y);
alpha[2][2] = -(node_xyz[3].z *node_xyz[1].x - node_xyz[3].z* node_xyz[0].x - node_xyz[0].z* node_xyz[1].x - node_xyz[1].z* node_xyz[3].x + node_xyz[1].z* node_xyz[0].x + node_xyz[0].z *node_xyz[3].x) /
(-node_xyz[3].z* node_xyz[0].x *node_xyz[1].y + node_xyz[0].x *node_xyz[3].z *node_xyz[2].y - node_xyz[0].x* node_xyz[2].z* node_xyz[3].y + node_xyz[2].z* node_xyz[0].x* node_xyz[1].y
	- node_xyz[1].z* node_xyz[2].y *node_xyz[0].x + node_xyz[1].z *node_xyz[3].y *node_xyz[0].x + node_xyz[3].z* node_xyz[1].x* node_xyz[0].y - node_xyz[0].z* node_xyz[3].y* node_xyz[1].x
	+ node_xyz[0].z* node_xyz[2].y *node_xyz[1].x - node_xyz[2].z *node_xyz[1].x *node_xyz[0].y - node_xyz[3].z* node_xyz[1].x* node_xyz[2].y + node_xyz[3].z* node_xyz[1].y* node_xyz[2].x
	- node_xyz[3].z* node_xyz[0].y *node_xyz[2].x + node_xyz[1].z *node_xyz[3].x *node_xyz[2].y + node_xyz[1].z* node_xyz[0].y* node_xyz[2].x - node_xyz[0].z* node_xyz[3].x* node_xyz[2].y
	+ node_xyz[0].z* node_xyz[3].x *node_xyz[1].y + node_xyz[2].z *node_xyz[1].x *node_xyz[3].y - node_xyz[2].z* node_xyz[1].y* node_xyz[3].x + node_xyz[2].z* node_xyz[0].y* node_xyz[3].x
	- node_xyz[1].z* node_xyz[2].x *node_xyz[3].y - node_xyz[1].z *node_xyz[0].y *node_xyz[3].x + node_xyz[0].z* node_xyz[2].x* node_xyz[3].y - node_xyz[0].z* node_xyz[2].x* node_xyz[1].y);
alpha[2][3] = (node_xyz[3].y*node_xyz[1].x - node_xyz[3].y*node_xyz[0].x - node_xyz[0].y*node_xyz[1].x - node_xyz[1].y*node_xyz[3].x + node_xyz[1].y*node_xyz[0].x + node_xyz[0].y*node_xyz[3].x) /
(-node_xyz[3].z*node_xyz[0].x*node_xyz[1].y + node_xyz[0].x*node_xyz[3].z*node_xyz[2].y - node_xyz[0].x*node_xyz[2].z*node_xyz[3].y + node_xyz[2].z*node_xyz[0].x*node_xyz[1].y
	- node_xyz[1].z*node_xyz[2].y*node_xyz[0].x + node_xyz[1].z*node_xyz[3].y*node_xyz[0].x + node_xyz[3].z*node_xyz[1].x*node_xyz[0].y - node_xyz[0].z*node_xyz[3].y*node_xyz[1].x
	+ node_xyz[0].z*node_xyz[2].y*node_xyz[1].x - node_xyz[2].z*node_xyz[1].x*node_xyz[0].y - node_xyz[3].z*node_xyz[1].x*node_xyz[2].y + node_xyz[3].z*node_xyz[1].y*node_xyz[2].x
	- node_xyz[3].z*node_xyz[0].y*node_xyz[2].x + node_xyz[1].z*node_xyz[3].x*node_xyz[2].y + node_xyz[1].z*node_xyz[0].y*node_xyz[2].x - node_xyz[0].z*node_xyz[3].x*node_xyz[2].y
	+ node_xyz[0].z*node_xyz[3].x*node_xyz[1].y + node_xyz[2].z*node_xyz[1].x*node_xyz[3].y - node_xyz[2].z*node_xyz[1].y*node_xyz[3].x + node_xyz[2].z*node_xyz[0].y*node_xyz[3].x
	- node_xyz[1].z*node_xyz[2].x*node_xyz[3].y - node_xyz[1].z*node_xyz[0].y*node_xyz[3].x + node_xyz[0].z*node_xyz[2].x*node_xyz[3].y - node_xyz[0].z*node_xyz[2].x*node_xyz[1].y);
alpha[3][0] = (-node_xyz[1].z*node_xyz[2].y*node_xyz[0].x + node_xyz[2].z*node_xyz[0].x*node_xyz[1].y + node_xyz[1].z*node_xyz[0].y*node_xyz[2].x + node_xyz[0].z*node_xyz[2].y*node_xyz[1].x - node_xyz[2].z*node_xyz[1].x*node_xyz[0].y - node_xyz[0].z*node_xyz[2].x*node_xyz[1].y) /
(-node_xyz[3].z*node_xyz[0].x*node_xyz[1].y + node_xyz[0].x*node_xyz[3].z*node_xyz[2].y - node_xyz[0].x*node_xyz[2].z*node_xyz[3].y + node_xyz[2].z*node_xyz[0].x*node_xyz[1].y - node_xyz[1].z*node_xyz[2].y*node_xyz[0].x + node_xyz[1].z*node_xyz[3].y*node_xyz[0].x
	+ node_xyz[3].z*node_xyz[1].x*node_xyz[0].y - node_xyz[0].z*node_xyz[3].y*node_xyz[1].x + node_xyz[0].z*node_xyz[2].y*node_xyz[1].x - node_xyz[2].z*node_xyz[1].x*node_xyz[0].y - node_xyz[3].z*node_xyz[1].x*node_xyz[2].y + node_xyz[3].z*node_xyz[1].y*node_xyz[2].x
	- node_xyz[3].z*node_xyz[0].y*node_xyz[2].x + node_xyz[1].z*node_xyz[3].x*node_xyz[2].y + node_xyz[1].z*node_xyz[0].y*node_xyz[2].x - node_xyz[0].z*node_xyz[3].x*node_xyz[2].y + node_xyz[0].z*node_xyz[3].x*node_xyz[1].y + node_xyz[2].z*node_xyz[1].x*node_xyz[3].y
	- node_xyz[2].z*node_xyz[1].y*node_xyz[3].x + node_xyz[2].z*node_xyz[0].y*node_xyz[3].x - node_xyz[1].z*node_xyz[2].x*node_xyz[3].y - node_xyz[1].z*node_xyz[0].y*node_xyz[3].x + node_xyz[0].z*node_xyz[2].x*node_xyz[3].y - node_xyz[0].z*node_xyz[2].x*node_xyz[1].y);
alpha[3][1] = -(-node_xyz[1].z*node_xyz[2].y + node_xyz[0].z*node_xyz[2].y - node_xyz[0].z*node_xyz[1].y + node_xyz[2].z*node_xyz[1].y - node_xyz[2].z*node_xyz[0].y + node_xyz[1].z*node_xyz[0].y) /
(-node_xyz[3].z*node_xyz[0].x*node_xyz[1].y + node_xyz[0].x*node_xyz[3].z*node_xyz[2].y - node_xyz[0].x*node_xyz[2].z*node_xyz[3].y + node_xyz[2].z*node_xyz[0].x*node_xyz[1].y
	- node_xyz[1].z*node_xyz[2].y*node_xyz[0].x + node_xyz[1].z*node_xyz[3].y*node_xyz[0].x + node_xyz[3].z*node_xyz[1].x*node_xyz[0].y - node_xyz[0].z*node_xyz[3].y*node_xyz[1].x
	+ node_xyz[0].z*node_xyz[2].y*node_xyz[1].x - node_xyz[2].z*node_xyz[1].x*node_xyz[0].y - node_xyz[3].z*node_xyz[1].x*node_xyz[2].y + node_xyz[3].z*node_xyz[1].y*node_xyz[2].x
	- node_xyz[3].z*node_xyz[0].y*node_xyz[2].x + node_xyz[1].z*node_xyz[3].x*node_xyz[2].y + node_xyz[1].z*node_xyz[0].y*node_xyz[2].x - node_xyz[0].z*node_xyz[3].x*node_xyz[2].y
	+ node_xyz[0].z*node_xyz[3].x*node_xyz[1].y + node_xyz[2].z*node_xyz[1].x*node_xyz[3].y - node_xyz[2].z*node_xyz[1].y*node_xyz[3].x + node_xyz[2].z*node_xyz[0].y*node_xyz[3].x
	- node_xyz[1].z*node_xyz[2].x*node_xyz[3].y - node_xyz[1].z*node_xyz[0].y*node_xyz[3].x + node_xyz[0].z*node_xyz[2].x*node_xyz[3].y - node_xyz[0].z*node_xyz[2].x*node_xyz[1].y);
alpha[3][2] = (node_xyz[2].z*node_xyz[1].x - node_xyz[2].z*node_xyz[0].x - node_xyz[0].z*node_xyz[1].x - node_xyz[1].z*node_xyz[2].x + node_xyz[1].z*node_xyz[0].x + node_xyz[0].z*node_xyz[2].x) /
(-node_xyz[1].z*node_xyz[2].y*node_xyz[0].x - node_xyz[3].z*node_xyz[0].x*node_xyz[1].y - node_xyz[0].x*node_xyz[2].z*node_xyz[3].y + node_xyz[2].z*node_xyz[0].x*node_xyz[1].y
	+ node_xyz[0].x*node_xyz[3].z*node_xyz[2].y + node_xyz[1].z*node_xyz[3].y*node_xyz[0].x + node_xyz[3].z*node_xyz[1].x*node_xyz[0].y - node_xyz[0].z*node_xyz[3].y*node_xyz[1].x
	+ node_xyz[0].z*node_xyz[2].y*node_xyz[1].x - node_xyz[2].z*node_xyz[1].x*node_xyz[0].y - node_xyz[3].z*node_xyz[1].x*node_xyz[2].y + node_xyz[3].z*node_xyz[1].y*node_xyz[2].x
	- node_xyz[3].z*node_xyz[0].y*node_xyz[2].x + node_xyz[1].z*node_xyz[3].x*node_xyz[2].y + node_xyz[1].z*node_xyz[0].y*node_xyz[2].x - node_xyz[0].z*node_xyz[3].x*node_xyz[2].y
	+ node_xyz[0].z*node_xyz[3].x*node_xyz[1].y + node_xyz[2].z*node_xyz[1].x*node_xyz[3].y - node_xyz[2].z*node_xyz[1].y*node_xyz[3].x + node_xyz[2].z*node_xyz[0].y*node_xyz[3].x
	- node_xyz[1].z*node_xyz[2].x*node_xyz[3].y - node_xyz[1].z*node_xyz[0].y*node_xyz[3].x + node_xyz[0].z*node_xyz[2].x*node_xyz[3].y - node_xyz[0].z*node_xyz[2].x*node_xyz[1].y);
alpha[3][3] = -(node_xyz[2].y*node_xyz[1].x - node_xyz[2].y*node_xyz[0].x - node_xyz[0].y*node_xyz[1].x - node_xyz[1].y*node_xyz[2].x + node_xyz[1].y*node_xyz[0].x + node_xyz[0].y*node_xyz[2].x) /
(-node_xyz[3].z*node_xyz[0].x*node_xyz[1].y + node_xyz[0].x*node_xyz[3].z*node_xyz[2].y - node_xyz[0].x*node_xyz[2].z*node_xyz[3].y + node_xyz[2].z*node_xyz[0].x*node_xyz[1].y
	- node_xyz[1].z*node_xyz[2].y*node_xyz[0].x + node_xyz[1].z*node_xyz[3].y*node_xyz[0].x + node_xyz[3].z*node_xyz[1].x*node_xyz[0].y - node_xyz[0].z*node_xyz[3].y*node_xyz[1].x
	+ node_xyz[0].z*node_xyz[2].y*node_xyz[1].x - node_xyz[2].z*node_xyz[1].x*node_xyz[0].y - node_xyz[3].z*node_xyz[1].x*node_xyz[2].y + node_xyz[3].z*node_xyz[1].y*node_xyz[2].x
	- node_xyz[3].z*node_xyz[0].y*node_xyz[2].x + node_xyz[1].z*node_xyz[3].x*node_xyz[2].y + node_xyz[1].z*node_xyz[0].y*node_xyz[2].x - node_xyz[0].z*node_xyz[3].x*node_xyz[2].y
	+ node_xyz[0].z*node_xyz[3].x*node_xyz[1].y + node_xyz[2].z*node_xyz[1].x*node_xyz[3].y - node_xyz[2].z*node_xyz[1].y*node_xyz[3].x + node_xyz[2].z*node_xyz[0].y*node_xyz[3].x
	- node_xyz[1].z*node_xyz[2].x*node_xyz[3].y - node_xyz[1].z*node_xyz[0].y*node_xyz[3].x + node_xyz[0].z*node_xyz[2].x*node_xyz[3].y - node_xyz[0].z*node_xyz[2].x*node_xyz[1].y);
		}
	};

	class Rectangle : public Shape
	{

	public:
		Rectangle()
		{
			this->id_nodes.resize(4);
			this->nodes.resize(4);
		};
		Rectangle(std::vector<Point<double> *> &P)
		{
			try
			{
				this->UpdateNodes(P);
			}
			catch (const std::exception&)
			{
				printf_s("Error: GeometryShape.cpp/Rectangle(std::vector<Point<double> *> &P)");
			}
		}
		Rectangle(Point<double> P1, Point<double> P2, Point<double> P3, Point<double> P4)
		{
			try
			{
				std::vector<Point<double> *> P(4);
				P[0] = &P1;
				P[1] = &P2;
				P[2] = &P3;
				P[3] = &P4;
				this->UpdateNodes(P);
			}
			catch (const std::exception&)
			{
				printf_s("Error: GeometryShape.cpp/Rectangle(std::vector<Point<double> *> &P)");
			}
		}
		double GetVolume()
		{
			try
			{
				Point<double> _nodes[3];
				for (int i = 0; i < 3; i++)
				{
					_nodes[i] = this->GetNode(i);
				}
				double Mx[3][3] = { { 1, _nodes[0].y, _nodes[0].z },
				{ 1, _nodes[1].y, _nodes[1].z },
				{ 1, _nodes[2].y, _nodes[2].z } };
				double My[3][3] = { { _nodes[0].x, 1, _nodes[0].z },
				{ _nodes[1].x, 1, _nodes[1].z },
				{ _nodes[2].x, 1, _nodes[2].z } };
				double Mz[3][3] = { { _nodes[0].x, _nodes[0].y, 1 },
				{ _nodes[1].x, _nodes[1].y, 1 },
				{ _nodes[2].x, _nodes[2].y, 1 } };
				double Sx = 0.5 * math::GetDeterminantForMatrix3x3(Mx);
				double Sy = 0.5 * math::GetDeterminantForMatrix3x3(My);
				double Sz = 0.5 * math::GetDeterminantForMatrix3x3(Mz);

				return sqrt(Sx*Sx + Sy * Sy + Sz * Sz);
			}
			catch (const std::exception&)
			{
				printf_s("Error: GeometryShape.cpp/GetVolume()");
				return 0.0;
			}
		}

		void SetIntegrationLaw(int points_count)
		{
			std::vector<double> xj(points_count);
			std::vector<double> q(points_count);

			/*    Points     */
			switch (points_count)
			{
			case 10:
				xj[0] = 0.9739065318;
				xj[1] = 0.8650633629;
				xj[2] = 0.6794095691;
				xj[3] = 0.4333953941;
				xj[4] = 0.1488743390;
				xj[5] = -xj[0];
				xj[6] = -xj[1];
				xj[7] = -xj[2];
				xj[8] = -xj[3];
				xj[9] = -xj[3];
				break;
			case 9:
				xj[0] = 0.0;
				xj[1] = 0.3242534234;
				xj[2] = 0.6133714327;
				xj[3] = 0.8360311073;
				xj[4] = 0.9681602395;
				xj[5] = -xj[1];
				xj[6] = -xj[2];
				xj[7] = -xj[3];
				xj[8] = -xj[4];
				break;
			case 8:
				xj[0] = 0.9602898565;
				xj[1] = 0.7966664774;
				xj[2] = 0.5255324099;
				xj[3] = 0.1834346425;
				xj[4] = -xj[0];
				xj[5] = -xj[1];
				xj[6] = -xj[2];
				xj[7] = -xj[3];
				break;
			case 7:
				xj[0] = 0.0;
				xj[1] = 0.4058451514;
				xj[2] = 0.7415311856;
				xj[3] = 0.9491079123;
				xj[4] = -xj[1];
				xj[5] = -xj[2];
				xj[6] = -xj[3];
				break;
			case 6:
				xj[0] = 0.2386191861;
				xj[1] = 0.6612093865;
				xj[2] = 0.9324695142;
				xj[3] = -xj[0];
				xj[4] = -xj[1];
				xj[5] = -xj[2];
				break;
			case 5:
				xj[0] = 0.;
				xj[1] = 1 / 21. * sqrt(245 - 14 * sqrt(70));
				xj[2] = 1 / 21. * sqrt(245 + 14 * sqrt(70));
				xj[3] = -xj[1];
				xj[4] = -xj[2];
				break;
			case 4:
				xj[0] = sqrt((3.0 - 2.0 * sqrt(6.0 / 5.0)) / 7.0);
				xj[1] = sqrt((3.0 + 2.0 * sqrt(6.0 / 5.0)) / 7.0);
				xj[2] = -xj[1];
				xj[3] = -xj[0];
				break;
			case 3:
				xj[0] = 0;
				xj[1] = 0.774596669241483;
				xj[2] = -xj[1];
				break;
			case 2:
				xj[0] = 1 / 3. * sqrt(3.);
				xj[1] = -xj[0];
				break;
			case 1:
				xj[0] = 0;
				break;
			}
			/*    Weights    */
			switch (points_count)
			{
			case 10:
				q[0] = 0.06667134202;
				q[1] = 0.1494513565;
				q[2] = 0.2190863541;
				q[3] = 0.2692667242;
				q[4] = 0.2955242232;
				q[5] = q[0];
				q[6] = q[1];
				q[7] = q[2];
				q[8] = q[3];
				q[9] = q[3];
				break;
			case 9:
				q[0] = 0.3356031191;
				q[1] = 0.3048968104;
				q[2] = 0.2710276447;
				q[3] = 0.1638265280;
				q[4] = 0.09199001997;
				q[5] = q[1];
				q[6] = q[2];
				q[7] = q[3];
				q[8] = q[4];
				break;
			case 8:
				q[0] = 0.1012285363;
				q[1] = 0.2223810344;
				q[2] = 0.3137066459;
				q[3] = 0.3626837833;
				q[4] = q[0];
				q[5] = q[1];
				q[6] = q[2];
				q[7] = q[3];
				break;
			case 7:
				q[0] = 0.4179591838;
				q[1] = 0.3818300507;
				q[2] = 0.2797053914;
				q[3] = 0.1294849662;
				q[4] = q[1];
				q[5] = q[2];
				q[6] = q[3];
				break;
			case 6:
				q[0] = 0.4679139346;
				q[1] = 0.3607615729;
				q[2] = 0.1713244924;
				q[3] = q[0];
				q[4] = q[1];
				q[5] = q[2];
				break;
			case 5:
				q[0] = 0.5688888888888888888;
				q[1] = 0.4786286705;
				q[2] = 0.2369268851;
				q[3] = q[1];
				q[4] = q[2];
				break;
			case 4:
				q[0] = (18 + sqrt(30.0)) / 36.0;
				q[1] = (18 - sqrt(30.0)) / 36.0;
				q[2] = q[1];
				q[3] = q[0];
				break;
			case 3:
				q[0] = 0.888888888888889;
				q[1] = 0.555555555555556;
				q[2] = q[1];
				break;
			case 2:
				q[0] = 1.0;
				q[1] = q[0];
				break;
			case 1:
				q[0] = 2.0;
				break;
			}

			///If you need update this block, you must re verify "CrackPropagation_3D"
			/*for (int i = 0; i < this->integration_parameters.w.size(); i++)
			{
				this->integration_parameters.w[i] *= detD;

				Point<double> Xtr;
				Xtr.x = (this->self_nodes[1].x - this->self_nodes[0].x)*this->integration_parameters.x[i].x +
					(this->self_nodes[2].x - this->self_nodes[0].x)*this->integration_parameters.x[i].y + this->self_nodes[0].x;
				Xtr.y = (this->self_nodes[1].y - this->self_nodes[0].y)*this->integration_parameters.x[i].x +
					(this->self_nodes[2].y - this->self_nodes[0].y)*this->integration_parameters.x[i].y + this->self_nodes[0].y;
				Xtr.z = this->self_nodes[0].z;

				this->integration_parameters.x[i] = Xtr;
			}*/
		}
		void SetIntegrationLaw(int points_count, std::vector<double> &xj, std::vector<double> &q)
		{
			xj.resize(points_count);
			q.resize(points_count);

			/*    Points     */
			switch (points_count)
			{
			case 10:
				xj[0] = 0.9739065318;
				xj[1] = 0.8650633629;
				xj[2] = 0.6794095691;
				xj[3] = 0.4333953941;
				xj[4] = 0.1488743390;
				xj[5] = -xj[0];
				xj[6] = -xj[1];
				xj[7] = -xj[2];
				xj[8] = -xj[3];
				xj[9] = -xj[3];
				break;
			case 9:
				xj[0] = 0.0;
				xj[1] = 0.3242534234;
				xj[2] = 0.6133714327;
				xj[3] = 0.8360311073;
				xj[4] = 0.9681602395;
				xj[5] = -xj[1];
				xj[6] = -xj[2];
				xj[7] = -xj[3];
				xj[8] = -xj[4];
				break;
			case 8:
				xj[0] = 0.9602898565;
				xj[1] = 0.7966664774;
				xj[2] = 0.5255324099;
				xj[3] = 0.1834346425;
				xj[4] = -xj[0];
				xj[5] = -xj[1];
				xj[6] = -xj[2];
				xj[7] = -xj[3];
				break;
			case 7:
				xj[0] = 0.0;
				xj[1] = 0.4058451514;
				xj[2] = 0.7415311856;
				xj[3] = 0.9491079123;
				xj[4] = -xj[1];
				xj[5] = -xj[2];
				xj[6] = -xj[3];
				break;
			case 6:
				xj[0] = 0.2386191861;
				xj[1] = 0.6612093865;
				xj[2] = 0.9324695142;
				xj[3] = -xj[0];
				xj[4] = -xj[1];
				xj[5] = -xj[2];
				break;
			case 5:
				xj[0] = 0.;
				xj[1] = 1 / 21. * sqrt(245 - 14 * sqrt(70));
				xj[2] = 1 / 21. * sqrt(245 + 14 * sqrt(70));
				xj[3] = -xj[1];
				xj[4] = -xj[2];
				break;
			case 4:
				xj[0] = sqrt((3.0 - 2.0 * sqrt(6.0 / 5.0)) / 7.0);
				xj[1] = sqrt((3.0 + 2.0 * sqrt(6.0 / 5.0)) / 7.0);
				xj[2] = -xj[1];
				xj[3] = -xj[0];
				break;
			case 3:
				xj[0] = 0;
				xj[1] = 0.774596669241483;
				xj[2] = -xj[1];
				break;
			case 2:
				xj[0] = 1 / 3. * sqrt(3.);
				xj[1] = -xj[0];
				break;
			case 1:
				xj[0] = 0;
				break;
			}
			/*    Weights    */
			switch (points_count)
			{
			case 10:
				q[0] = 0.06667134202;
				q[1] = 0.1494513565;
				q[2] = 0.2190863541;
				q[3] = 0.2692667242;
				q[4] = 0.2955242232;
				q[5] = q[0];
				q[6] = q[1];
				q[7] = q[2];
				q[8] = q[3];
				q[9] = q[3];
				break;
			case 9:
				q[0] = 0.3356031191;
				q[1] = 0.3048968104;
				q[2] = 0.2710276447;
				q[3] = 0.1638265280;
				q[4] = 0.09199001997;
				q[5] = q[1];
				q[6] = q[2];
				q[7] = q[3];
				q[8] = q[4];
				break;
			case 8:
				q[0] = 0.1012285363;
				q[1] = 0.2223810344;
				q[2] = 0.3137066459;
				q[3] = 0.3626837833;
				q[4] = q[0];
				q[5] = q[1];
				q[6] = q[2];
				q[7] = q[3];
				break;
			case 7:
				q[0] = 0.4179591838;
				q[1] = 0.3818300507;
				q[2] = 0.2797053914;
				q[3] = 0.1294849662;
				q[4] = q[1];
				q[5] = q[2];
				q[6] = q[3];
				break;
			case 6:
				q[0] = 0.4679139346;
				q[1] = 0.3607615729;
				q[2] = 0.1713244924;
				q[3] = q[0];
				q[4] = q[1];
				q[5] = q[2];
				break;
			case 5:
				q[0] = 0.5688888888888888888;
				q[1] = 0.4786286705;
				q[2] = 0.2369268851;
				q[3] = q[1];
				q[4] = q[2];
				break;
			case 4:
				q[0] = (18 + sqrt(30.0)) / 36.0;
				q[1] = (18 - sqrt(30.0)) / 36.0;
				q[2] = q[1];
				q[3] = q[0];
				break;
			case 3:
				q[0] = 0.888888888888889;
				q[1] = 0.555555555555556;
				q[2] = q[1];
				break;
			case 2:
				q[0] = 1.0;
				q[1] = q[0];
				break;
			case 1:
				q[0] = 2.0;
				break;
			}

			///If you need update this block, you must re verify "CrackPropagation_3D"
			/*for (int i = 0; i < this->integration_parameters.w.size(); i++)
			{
				this->integration_parameters.w[i] *= detD;

				Point<double> Xtr;
				Xtr.x = (this->self_nodes[1].x - this->self_nodes[0].x)*this->integration_parameters.x[i].x +
					(this->self_nodes[2].x - this->self_nodes[0].x)*this->integration_parameters.x[i].y + this->self_nodes[0].x;
				Xtr.y = (this->self_nodes[1].y - this->self_nodes[0].y)*this->integration_parameters.x[i].x +
					(this->self_nodes[2].y - this->self_nodes[0].y)*this->integration_parameters.x[i].y + this->self_nodes[0].y;
				Xtr.z = this->self_nodes[0].z;

				this->integration_parameters.x[i] = Xtr;
			}*/
		}

		bool IsContainThePoint(Point<double> A)
		{
			try
			{
				
			}
			catch (const std::exception&)
			{
				printf_s("Error: GeometryShape.cpp/IsContainThePoint(Point<double> A)");
				return false;
			}
		}
		bool IsContainThePoint(Point<double> A, double &length)
		{
			try
			{

			}
			catch (const std::exception&)
			{
				printf_s("Error: GeometryShape.cpp/IsContainThePoint(Point<double> A)");
				return false;
			}
		}
	};

	class Polygon : public Shape
	{
	private:
		std::vector<Triangle> sub_elements_inself;
		Point<double> PolygonCentr;

	public:
		Polygon()
		{
			this->id_nodes.resize(4);
			this->nodes.resize(4);
		};
		Polygon(std::vector<Point<double>*>& P)
		{
			try
			{
				this->UpdateNodes(P);
			}
			catch (const std::exception&)
			{
				printf_s("Error: GeometryShape.h/Polygon(std::vector<Point<double> *> &P)");
			}
		}
		//id_nodes[id_face][id vertex] in plyhedron numeration (local!)
		void CreateSubMeshInSelf()
		{
			try
			{
				this->PolygonCentr = Point<double>(0.0, 0.0, 0.0);
				for (int i = 0; i < this->GetNodesCount(); i++)
				{
					PolygonCentr += this->GetNode(i);
				}
				PolygonCentr /= this->GetNodesCount();

				int id_centr_point = this->GetNodesCount();

				for (int i = 0; i < id_nodes.size(); i++)
				{
					std::vector<int> id_nodes_in_triag(3);
					std::vector<Point<double>*> xyz_in_triag(3);
					id_nodes_in_triag[0] = id_centr_point;
					id_nodes_in_triag[1] = i;
					id_nodes_in_triag[2] = (i + 1) % this->GetNodesCount();

					xyz_in_triag[0] = &this->PolygonCentr;
					xyz_in_triag[1] = this->nodes[id_nodes_in_triag[1]];
					xyz_in_triag[2] = this->nodes[id_nodes_in_triag[2]];

					this->sub_elements_inself.push_back(Triangle(id_nodes_in_triag, xyz_in_triag));
				}
			}
			catch (const std::exception&)
			{
				printf_s("Error: GeometryShape.h/SetGeometry(std::vector<std::vector<int>>& id_nodes");
			}
		};
		void SetGeometry(std::vector<int>& id_nodes, std::vector<Point<double>*>& P)
		{
			try
			{
				this->SetIdNodes(id_nodes);
				this->UpdateNodes(P);
				this->CreateSubMeshInSelf();
			}
			catch (const std::exception&)
			{
				printf_s("Error: GeometryShape.h/SetGeometry(std::vector<std::vector<int>>& id_nodes, std::vector<Point<double>*>& P)");
			}
		};
		double GetVolume()
		{
			try
			{
				double V = 0;

				for (int i = 0; i < this->sub_elements_inself.size(); i++)
				{
					V += sub_elements_inself[i].GetVolume();
				}
				return V;
			}
			catch (const std::exception&)
			{
				printf_s("Error: GeometryShape.cpp/GetVolume()");
				return 0.0;
			}
		}

		//via sub tetrahedral grid
		void SetIntegrationLaw(int points_count)
		{
			this->integration_parameters.w.resize(this->sub_elements_inself.size() * points_count);
			this->integration_parameters.x.resize(this->sub_elements_inself.size() * points_count);
			for (int i = 0; i < this->sub_elements_inself.size(); i++)
			{
				this->sub_elements_inself[i].SetIntegrationLaw(points_count);
				for (int j = 0; j < points_count; j++)
				{
					this->integration_parameters.w[i * points_count + j] = this->sub_elements_inself[i].GetIntegrationWeight(j);
					this->integration_parameters.x[i * points_count + j] = this->sub_elements_inself[i].GetIntegrationPoint(j);
				}
			}
		}

		bool IsContainThePoint(Point<double> A)
		{
			try
			{
				bool result = false;
				for (int i = 0; result == false && i < this->sub_elements_inself.size(); i++)
				{
					result = this->sub_elements_inself[i].IsContainThePoint(A);
				}
				return result;
			}
			catch (const std::exception&)
			{
				printf_s("Error: GeometryShape.cpp/IsContainThePoint(Point<double> A)");
				return false;
			}
		}
		bool IsContainThePoint(Point<double> A, double& length)
		{
			try
			{
				bool result = false;
				length = 1e+302;
				for (int i = 0; result == false && i < this->sub_elements_inself.size(); i++)
				{
					double len;
					result = this->sub_elements_inself[i].IsContainThePoint(A, len);
					if (len < length) length = len;
				}
				return result;
			}
			catch (const std::exception&)
			{
				printf_s("Error: GeometryShape.cpp/IsContainThePoint(Point<double> A)");
				return false;
			}
		}
	};
	class Polyhedron : public Shape
	{
	private: 
		std::vector<Tetrahedron> sub_elements_inself;
		Point<double> PolyhedronCentr;
		std::vector<Point<double>> PolygonCentr;
		std::vector<Polygon> faces;

	public:
		Polyhedron()
		{
			this->id_nodes.resize(4);
			this->nodes.resize(4);
		};
		Polyhedron(std::vector<Point<double>*>& P)
		{
			try
			{
				this->UpdateNodes(P);
			}
			catch (const std::exception&)
			{
				printf_s("Error: GeometryShape.h/Polyhedron(std::vector<Point<double> *> &P)");
			}
		}
		//id_nodes[id_face][id vertex] in plyhedron numeration (local!)
		void CreateSubMeshInSelf(std::vector<std::vector<int>>& id_nodes)
		{
			try
			{
				this->PolyhedronCentr = Point<double>(0.0, 0.0, 0.0);
				for (int i = 0; i < this->GetNodesCount(); i++)
				{
					PolyhedronCentr += this->GetNode(i);
				}
				PolyhedronCentr /= this->GetNodesCount();

				int id_centr_point = this->GetNodesCount();
				
				//by faces
				this->faces.resize(id_nodes.size());
				this->PolygonCentr.resize(id_nodes.size());
				for (int f = 0; f < id_nodes.size(); f++)
				{
					int id_face_point = id_centr_point + 1 + f;

					std::vector<Point<double>*> P_faces(id_nodes[f].size());
					for (int nf = 0; nf < P_faces.size(); nf++)
					{
						P_faces[nf] = this->GetPtrNode(id_nodes[f][nf]);
					}
					faces[f].SetGeometry(id_nodes[f], P_faces);

					for (int i = 0; i < id_nodes[f].size(); i++)
					{
						this->PolygonCentr[f] += this->GetNode(id_nodes[f][i]);
					}
					this->PolygonCentr[f] /= id_nodes[f].size();

					for (int i = 0; i < id_nodes[f].size(); i++)
					{
						std::vector<int> id_nodes_in_tert(4);
						std::vector<Point<double>*> xyz_in_tert(4);
						id_nodes_in_tert[0] = id_centr_point;
						id_nodes_in_tert[1] = id_face_point;
						id_nodes_in_tert[2] = id_nodes[f][i];
						id_nodes_in_tert[3] = id_nodes[f][(i+1) % id_nodes[f].size()];

						xyz_in_tert[0] = &this->PolyhedronCentr;
						xyz_in_tert[1] = &this->PolygonCentr[f];
						xyz_in_tert[2] = this->nodes[id_nodes_in_tert[2]];
						xyz_in_tert[3] = this->nodes[id_nodes_in_tert[3]];

						this->sub_elements_inself.push_back(Tetrahedron(id_nodes_in_tert, xyz_in_tert));
					}
				}
			}
			catch (const std::exception&)
			{
				printf_s("Error: GeometryShape.h/SetGeometry(std::vector<std::vector<int>>& id_nodes");
			}
		};
		void SetGeometry(std::vector<int> &id_global_nodes, std::vector<std::vector<int>>& id_nodes_in_faces, std::vector<Point<double>*>& P)
		{
			try
			{
				this->SetIdNodes(id_global_nodes);
				this->UpdateNodes(P);
				this->CreateSubMeshInSelf(id_nodes_in_faces);
			}
			catch (const std::exception&)
			{
				printf_s("Error: GeometryShape.h/SetGeometry(std::vector<std::vector<int>>& id_nodes, std::vector<Point<double>*>& P)");
			}
		};
		int GetCountFaces()
		{
			return this->faces.size();
		}
		void GetGlobalIdFace(int face_num, std::vector<int>& ids)
		{
			ids.resize(this->faces[face_num].GetNodesCount());
			for (int i = 0; i < ids.size(); i++)
			{
				ids[i] = this->GetIdNode(this->faces[face_num].GetIdNode(i));
			}
		}
		int GetCountEdgeInFace(int face_num)
		{
			return this->faces[face_num].GetNodesCount();
		}
		void GetGlobalIdEdge(int face_num, int edge_num, std::vector<int> &ids) 
		{
			ids.resize(2);
			ids[0] = this->GetIdNode(this->faces[face_num].GetIdNode(edge_num));
			ids[1] = this->GetIdNode(this->faces[face_num].GetIdNode((edge_num+1)%this->faces[face_num].GetNodesCount()));
		}

		double GetVolume()
		{
			try
			{
				double V = 0;

				for (int i = 0; i < this->sub_elements_inself.size(); i++)
				{
					V += sub_elements_inself[i].GetVolume();
				}
				return V;
			}
			catch (const std::exception&)
			{
				printf_s("Error: GeometryShape.cpp/GetVolume()");
				return 0.0;
			}
		}

		//via sub tetrahedral grid
		void SetIntegrationLaw(int points_count)
		{
			this->integration_parameters.w.resize(this->sub_elements_inself.size() * points_count);
			this->integration_parameters.x.resize(this->sub_elements_inself.size() * points_count);
			for (int i = 0; i < this->sub_elements_inself.size(); i++)
			{
				this->sub_elements_inself[i].SetIntegrationLaw(points_count);
				for (int j = 0; j < points_count; j++)
				{
					this->integration_parameters.w[i * points_count + j] = this->sub_elements_inself[i].GetIntegrationWeight(j);
					this->integration_parameters.x[i * points_count + j] = this->sub_elements_inself[i].GetIntegrationPoint(j);
				}
			}
		}

		bool IsContainThePoint(Point<double> A)
		{
			try
			{
				bool result = false;
				for (int i = 0; result == false && i < this->sub_elements_inself.size(); i++)
				{
					result = this->sub_elements_inself[i].IsContainThePoint(A);
				}
				return result;
			}
			catch (const std::exception&)
			{
				printf_s("Error: GeometryShape.cpp/IsContainThePoint(Point<double> A)");
				return false;
			}
		}
		bool IsContainThePoint(Point<double> A, double& length)
		{
			try
			{
				bool result = false;
				length = 1e+302;
				for (int i = 0; result == false && i < this->sub_elements_inself.size(); i++)
				{
					double len;
					result = this->sub_elements_inself[i].IsContainThePoint(A, len);
					if (len < length) length = len;
				}
				return result;
			}
			catch (const std::exception&)
			{
				printf_s("Error: GeometryShape.cpp/IsContainThePoint(Point<double> A)");
				return false;
			}
		}
	};



	class Crack
	{
	private:
		struct propogation_points {
			int id_old_point;
			int id_new_point;
			std::vector<Point<double>> vect_of_propogation;
			std::vector<double> length_of_propogation;
			std::vector<double> Qgrad;
			std::vector<double> G;
			std::vector<double> J1;
			std::vector<double> J2;
			Point<double> _vect_of_propogation;
			double _length_of_propogation;
			double _Qgrad;
			double _G;
			double _J1;
			double _J2;
			
			struct newCoordSystem {
				std::vector<Point<double>> _X; //new basis in OLD coord
				Point<double> O; //New coord centr in OLD coord
				Point<double> _O; //New coord centr in NEW coord
				Tensor2Rank3D A; //transfer from Old into New
				Tensor2Rank3D _A; //transfer from New into Old
				Point<double> Offset; //offset from Old into New
				Point<double> _Offset; //offset from New into Old 

				template <typename T>
				newCoordSystem(T A)
				{
					math::MakeCopyVector_A_into_B(A._X, this->_X);
					this->O = A.O;
					this->_O = A._O;
					this->A = A.A;
					this->_A = A._A;
					this->Offset = A.Offset;
					this->_Offset = A._Offset;
				}

				Point<double> MakePropogation(double Qgrad, double lenght)
				{
					double Qrad = Qgrad * M_PI / 180;
					Point<double> _R, R;
					_R.z = this->_O.z;
					_R.x = cos(Qrad) * lenght;
					_R.y = sqrt(lenght*lenght - _R.x*_R.x) *math::GetSignum(Qrad);
					R = this->transfer_from_NEW_into_OLD(_R);
					R -= this->O;
					return R;
				}

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

			std::vector <newCoordSystem> base_coord_system;
		};

	public:
		std::vector<Triangle> triangles; //хранится 3мя узлами
		std::vector<bool> triangles_in_front;

		Point<double> min_point, max_point;

		struct segment {
			int id_base_triangles;
			int id_left, id_right;
		};
		std::vector<std::vector<segment>> fronts; //fronts[номер фронта][номер сегмента] 
		std::vector<Point<double>> xyz;
		std::vector<std::vector<int>> id_points_in_fronts;
		std::vector<std::vector<propogation_points>> propogation; //[front][points in front]

		double R; //радиус ввода сингулярных функций
		int smoothing_coefficient; 

		Crack() { return; };
		Crack(char in_file_name[1000])
		{
			this->Input(in_file_name);
		}
		bool Input(char in_file_name[1000])
		{
			try {
				FILE *f;
				fopen_s(&f, in_file_name, "r");

				if (f == NULL) return false;

				int num_nodes;
				fscanf_s(f, "%d", &num_nodes);
				this->xyz.resize(num_nodes);
				for (int i = 0; i < num_nodes; i++)
				{
					fscanf_s(f, "%lf %lf %lf", &this->xyz[i].x, &this->xyz[i].y, &this->xyz[i].z);
				}

				this->min_point = this->xyz[0];
				this->max_point = this->xyz[0];
				for (int i = 1; i < this->xyz.size(); i++)
				{
					if (min_point.x > this->xyz[i].x) min_point.x = this->xyz[i].x;
					if (min_point.y > this->xyz[i].y) min_point.y = this->xyz[i].y;
					if (min_point.z > this->xyz[i].z) min_point.z = this->xyz[i].z;

					if (max_point.x < this->xyz[i].x) max_point.x = this->xyz[i].x;
					if (max_point.y < this->xyz[i].y) max_point.y = this->xyz[i].y;
					if (max_point.z < this->xyz[i].z) max_point.z = this->xyz[i].z;
				}

				int num_Planes;
				fscanf_s(f, "%d", &num_Planes);
				this->triangles.resize(num_Planes);
				this->triangles_in_front.resize(num_Planes);
				for (int i = 0; i < num_Planes; i++)
				{
					this->triangles_in_front[i] = false;
					int id_nodes[3];
					fscanf_s(f, "%d %d %d",
						&id_nodes[0],
						&id_nodes[1],
						&id_nodes[2]);
					this->triangles[i].SetIdNode(0, id_nodes[0]);
					this->triangles[i].SetIdNode(1, id_nodes[1]);
					this->triangles[i].SetIdNode(2, id_nodes[2]);
					this->triangles[i].SetNode(0, &this->xyz[id_nodes[0]]);
					this->triangles[i].SetNode(1, &this->xyz[id_nodes[1]]);
					this->triangles[i].SetNode(2, &this->xyz[id_nodes[2]]);
				}

				int num_fronts;
				fscanf_s(f, "%d", &num_fronts);
				this->fronts.resize(num_fronts);
				this->id_points_in_fronts.resize(num_fronts);
				for (int i = 0; i < num_fronts; i++)
				{
					int num_elem_in_front;
					fscanf_s(f, "%d", &num_elem_in_front);
					this->fronts[i].resize(num_elem_in_front);
					for (int j = 0; j < num_elem_in_front; j++)
					{
						fscanf_s(f, "%d %d %d", &this->fronts[i][j].id_base_triangles,
							&this->fronts[i][j].id_left, &this->fronts[i][j].id_right);

						this->triangles_in_front[this->fronts[i][j].id_base_triangles] = true;
					}

					//closed front
					if (this->fronts[i][0].id_left == this->fronts[i][num_elem_in_front - 1].id_right)
					{
						this->id_points_in_fronts[i].push_back(this->fronts[i][num_elem_in_front - 1].id_left);
						for (int j = 0; j < num_elem_in_front; j++)
						{
							this->id_points_in_fronts[i].push_back(this->fronts[i][j].id_left);
						}
						this->id_points_in_fronts[i].push_back(this->fronts[i][num_elem_in_front - 1].id_right);
						this->id_points_in_fronts[i].push_back(this->fronts[i][0].id_right);
					}
					else {
						this->id_points_in_fronts[i].push_back(this->fronts[i][0].id_left);
						for (int j = 0; j < num_elem_in_front; j++)
						{
							this->id_points_in_fronts[i].push_back(this->fronts[i][j].id_left);
						}
						this->id_points_in_fronts[i].push_back(this->fronts[i][num_elem_in_front - 1].id_right);
						this->id_points_in_fronts[i].push_back(this->fronts[i][num_elem_in_front - 1].id_right);
					}
				}

				fclose(f);
				return true;
			}
			catch (const std::exception&)
			{
				printf_s("Error: GeometryShape.cpp/bool geometry::Crack::Input(char in_file_name[1000])");
			}
		}
	};

	
}