#pragma once
#include "../Library/GeometryGrid.h"
#include "../Library/GeometryShape.h"
#include "../Library/TopologyShape.h"
#include "../Library/FunctionalShape.h"
#include "../Library/DenseMatrix.h"
#include "../Library/Tensor.h"
#include "../Library/Math.h"
#include <stdio.h>
#include <vector>

namespace MultiXFEM {
	class FiniteElement_forMech :
		public geometry::Tetrahedron,
		//public topology::Tetrahedron<topology::lower::Triangle, topology::upper::EmptyElement>,
		public topology::Tetrahedron<topology::lower::Vertex, topology::upper::EmptyElement>,
		public functional::Shape<int, Point<double>>
	{
	public:
		math::SimpleGrid SubGrid_for_integration;
		math::SimpleGrid SubGrid_for_integr_byTriangles;

		FiniteElement_forMech() { return; };
		~FiniteElement_forMech() { return; };

		void SolveLocalMatrix(DenseMatrix<Tensor2Rank3D, Point<double>> &local_matix, std::function<std::vector<std::vector<double>>(Point<double>)> &koefD)
		{
			try {
				local_matix.SetSize(this->GetDOFsCount());

				this->SetIntegrationLaw(4);

				/*if (this->GetIdDomain() == 2) 
					printf_s("\nThis elem[%d] is singular\n", this->GetSelfGlobalId());*/

				/*int bf = this->GetDOFsCount();
				if (bf > 4)
				{
					bf *= 1;
					printf_s("In element[%d] there are additional DOFs!\n", this->GetSelfGlobalId());
				}*/

				for (int _bf_local_id_I = 0; _bf_local_id_I < this->GetDOFsCount(); _bf_local_id_I++)
				{
					for (int _bf_local_id_J = _bf_local_id_I; _bf_local_id_J < this->GetDOFsCount(); _bf_local_id_J++)
					{

						auto D_basis_function_I = this->GetDerivativeOfBasisFunctionInLocalID(_bf_local_id_I);
						auto D_basis_function_J = this->GetDerivativeOfBasisFunctionInLocalID(_bf_local_id_J);
						auto basis_function_I = this->GetBasisFunctionInLocalID(_bf_local_id_I);
						auto basis_function_J = this->GetBasisFunctionInLocalID(_bf_local_id_J);

						Point<double> o = this->GetWeightCentr();
						auto val = (*this->GetDerivativeOfBasisFunctionInLocalID(0))(o);

						std::function<Tensor2Rank3D(Point<double>)> StiffnessMatrix = [&_bf_local_id_I, &_bf_local_id_J, &koefD, &D_basis_function_I, &D_basis_function_J]
						(Point<double> X) -> Tensor2Rank3D
						{
							Tensor2Rank3D result;

							std::vector<std::vector<double>> grad_right, gradT_left;
							math::ResizeVector(grad_right, 6, 3);
							math::ResizeVector(gradT_left, 3, 6);
							auto D = koefD(X);
							Point<Point<double>> derevativeBF_I_inX = (*D_basis_function_I)(X);
							Point<Point<double>> derevativeBF_J_inX = (*D_basis_function_J)(X);

							grad_right[0][0] = derevativeBF_I_inX.x.x;
							grad_right[1][1] = derevativeBF_I_inX.y.y;
							grad_right[2][2] = derevativeBF_I_inX.z.z;
							/*grad_right[3][1] = derevativeBF_I_inX.z.y;
							grad_right[3][2] = derevativeBF_I_inX.y.z;
							grad_right[4][0] = derevativeBF_I_inX.z.x;
							grad_right[4][2] = derevativeBF_I_inX.x.z;
							grad_right[5][0] = derevativeBF_I_inX.y.x;
							grad_right[5][1] = derevativeBF_I_inX.x.y;*/
							grad_right[4][2] = derevativeBF_I_inX.z.y;
							grad_right[4][1] = derevativeBF_I_inX.y.z;
							grad_right[5][2] = derevativeBF_I_inX.z.x;
							grad_right[5][0] = derevativeBF_I_inX.x.z;
							grad_right[3][1] = derevativeBF_I_inX.y.x;
							grad_right[3][0] = derevativeBF_I_inX.x.y;

							gradT_left[0][0] = derevativeBF_J_inX.x.x;
							gradT_left[1][1] = derevativeBF_J_inX.y.y;
							gradT_left[2][2] = derevativeBF_J_inX.z.z;
							/*gradT_left[1][3] = derevativeBF_J_inX.z.y;
							gradT_left[2][3] = derevativeBF_J_inX.y.z;
							gradT_left[0][4] = derevativeBF_J_inX.z.x;
							gradT_left[2][4] = derevativeBF_J_inX.x.z;
							gradT_left[0][5] = derevativeBF_J_inX.y.x;
							gradT_left[1][5] = derevativeBF_J_inX.x.y;*/
							gradT_left[2][4] = derevativeBF_J_inX.z.y;
							gradT_left[1][4] = derevativeBF_J_inX.y.z;
							gradT_left[2][5] = derevativeBF_J_inX.z.x;
							gradT_left[0][5] = derevativeBF_J_inX.x.z;
							gradT_left[1][3] = derevativeBF_J_inX.y.x;
							gradT_left[0][3] = derevativeBF_J_inX.x.y;

							/*std::vector<std::vector<double>> mult_D_grad;
							math::MultMatrixMatrix(D, grad_right, mult_D_grad);
							std::vector<std::vector<double>> mult_gradT_D_grad;
							math::MultMatrixMatrix(gradT_left, mult_D_grad, mult_gradT_D_grad);*/
							std::vector<std::vector<double>> mult_gradT_D;
							math::MultMatrixMatrix(gradT_left, D, mult_gradT_D);
							std::vector<std::vector<double>> mult_gradT_D_grad;
							math::MultMatrixMatrix(mult_gradT_D, grad_right, mult_gradT_D_grad);

							result = mult_gradT_D_grad;
							//result = 1.0;
							return result;
						};
						std::function<Tensor2Rank3D(Point<double>)> MassMatrix = [&_bf_local_id_I, &_bf_local_id_J, &basis_function_I, &basis_function_J]
						(Point<double> X) -> Tensor2Rank3D
						{
							Tensor2Rank3D result;

							std::vector<std::vector<double>> grad_right, gradT_left;
							math::ResizeVector(grad_right, 6, 3);
							math::ResizeVector(gradT_left, 3, 6);
							Point<double> BF_I_inX = (*basis_function_I)(X);
							Point<double> BF_J_inX = (*basis_function_J)(X);

							/*std::vector<std::vector<double>> mult_D_grad;
							math::MultMatrixMatrix(D, grad_right, mult_D_grad);
							std::vector<std::vector<double>> mult_gradT_D_grad;
							math::MultMatrixMatrix(gradT_left, mult_D_grad, mult_gradT_D_grad);*/
							std::vector<std::vector<double>> mult_res;
							double cohesive_koef=-7.4e+6;
							math::ResizeVector(mult_res, 3, 3);
							mult_res[0][0] = BF_I_inX.x*BF_J_inX.x*cohesive_koef;
							mult_res[0][1] = BF_I_inX.x*BF_J_inX.y*cohesive_koef;
							mult_res[0][2] = BF_I_inX.x*BF_J_inX.z*cohesive_koef;

							mult_res[1][0] = BF_I_inX.y*BF_J_inX.x*cohesive_koef;
							mult_res[1][1] = BF_I_inX.y*BF_J_inX.y*cohesive_koef;
							mult_res[1][2] = BF_I_inX.y*BF_J_inX.z*cohesive_koef;

							mult_res[2][0] = BF_I_inX.z*BF_J_inX.x*cohesive_koef;
							mult_res[2][1] = BF_I_inX.z*BF_J_inX.y*cohesive_koef;
							mult_res[2][2] = BF_I_inX.z*BF_J_inX.z*cohesive_koef;

							result = mult_res;
							//result = 1.0;
							return result;
						};

						double V = this->GetVolume();
						local_matix.A[_bf_local_id_I][_bf_local_id_J] = this->SolveIntegral(StiffnessMatrix);
						local_matix.A[_bf_local_id_J][_bf_local_id_I] = local_matix.A[_bf_local_id_I][_bf_local_id_J].T();

						//костыль дл€ гидроразрыва!!!
						/*if (_bf_local_id_I >= 4 && _bf_local_id_J >= 4)
						{
							for (int tr = 0; tr < this->SubGrid_for_integr_byTriangles.nvtr.size(); tr++)
							{
								geometry::Triangle triangle(
									this->SubGrid_for_integr_byTriangles.xyz[this->SubGrid_for_integr_byTriangles.nvtr[tr][0]],
									this->SubGrid_for_integr_byTriangles.xyz[this->SubGrid_for_integr_byTriangles.nvtr[tr][1]],
									this->SubGrid_for_integr_byTriangles.xyz[this->SubGrid_for_integr_byTriangles.nvtr[tr][2]]);

								triangle.SetIntegrationLaw(4);
								Tensor2Rank3D res = triangle.SolveIntegral(MassMatrix);
								local_matix.A[_bf_local_id_I][_bf_local_id_J] += res;
							}
						}*/
					}
				}
				/*if (this->GetIdDomain() == 2)
					printf_s("\nThis elem[%d] is singular\n", this->GetSelfGlobalId());*/
			}
			catch (const std::exception&)
			{
				printf_s("Error: XFEM/MultiXFEM_Grid.h/void SolveLocalMatrix(DenseMatrix<Tensor2Rank3D> &local_matix, std::function<std::vector<std::vector<double>>(Point)> &koefD)\n");
			}
		}
		void SolveLocalCohesiveMatrix(DenseMatrix<Tensor2Rank3D, Point<double>> &local_matix, std::function<double(Point<double>)> &koefCohesive)
		{
			try {
				local_matix.SetSize(this->GetDOFsCount());

				this->SetIntegrationLaw(4);

				/*if (this->GetIdDomain() == 2)
					printf_s("\nThis elem[%d] is singular\n", this->GetSelfGlobalId());*/

				for (int _bf_local_id_I = 4; _bf_local_id_I < this->GetDOFsCount(); _bf_local_id_I++)
				{
					for (int _bf_local_id_J = 4; _bf_local_id_J < this->GetDOFsCount(); _bf_local_id_J++)
					{

						//auto D_basis_function_I = this->GetDerivativeOfBasisFunctionInLocalID(_bf_local_id_I);
						//auto D_basis_function_J = this->GetDerivativeOfBasisFunctionInLocalID(_bf_local_id_J);
						auto basis_function_I = this->GetBasisFunctionInLocalID(_bf_local_id_I);
						auto basis_function_J = this->GetBasisFunctionInLocalID(_bf_local_id_J);

						Point<double> o = this->GetWeightCentr();
						auto val = (*this->GetDerivativeOfBasisFunctionInLocalID(0))(o);

						std::function<Tensor2Rank3D(Point<double>)> MassMatrix = [&koefCohesive, &_bf_local_id_I, &_bf_local_id_J, &basis_function_I, &basis_function_J]
						(Point<double> X) -> Tensor2Rank3D
						{
							Tensor2Rank3D result;

							std::vector<std::vector<double>> grad_right, gradT_left;
							math::ResizeVector(grad_right, 6, 3);
							math::ResizeVector(gradT_left, 3, 6);
							Point<double> BF_I_inX = (*basis_function_I)(X);
							Point<double> BF_J_inX = (*basis_function_J)(X);

							/*std::vector<std::vector<double>> mult_D_grad;
							math::MultMatrixMatrix(D, grad_right, mult_D_grad);
							std::vector<std::vector<double>> mult_gradT_D_grad;
							math::MultMatrixMatrix(gradT_left, mult_D_grad, mult_gradT_D_grad);*/
							std::vector<std::vector<double>> mult_res;
							math::ResizeVector(mult_res, 3, 3);
							mult_res[0][0] = BF_I_inX.x*BF_J_inX.x*koefCohesive(X);
							mult_res[0][1] = BF_I_inX.x*BF_J_inX.y*koefCohesive(X);
							mult_res[0][2] = BF_I_inX.x*BF_J_inX.z*koefCohesive(X);

							mult_res[1][0] = BF_I_inX.y*BF_J_inX.x*koefCohesive(X);
							mult_res[1][1] = BF_I_inX.y*BF_J_inX.y*koefCohesive(X);
							mult_res[1][2] = BF_I_inX.y*BF_J_inX.z*koefCohesive(X);

							mult_res[2][0] = BF_I_inX.z*BF_J_inX.x*koefCohesive(X);
							mult_res[2][1] = BF_I_inX.z*BF_J_inX.y*koefCohesive(X);
							mult_res[2][2] = BF_I_inX.z*BF_J_inX.z*koefCohesive(X);

							result = mult_res;
							//result = 1.0;
							return result;
						};
						
						for (int tr = 0; tr < this->SubGrid_for_integr_byTriangles.nvtr.size(); tr++)
						{
							geometry::Triangle triangle(
								this->SubGrid_for_integr_byTriangles.xyz[this->SubGrid_for_integr_byTriangles.nvtr[tr][0]],
								this->SubGrid_for_integr_byTriangles.xyz[this->SubGrid_for_integr_byTriangles.nvtr[tr][1]],
								this->SubGrid_for_integr_byTriangles.xyz[this->SubGrid_for_integr_byTriangles.nvtr[tr][2]]);

							triangle.SetIntegrationLaw(4);
							Tensor2Rank3D res = triangle.SolveIntegral(MassMatrix);
							local_matix.A[_bf_local_id_I][_bf_local_id_J] += res;
							//local_matix.A[_bf_local_id_J][_bf_local_id_I] += res;
						}
						/*if (this->SubGrid_for_integr_byTriangles.nvtr.size() == 0)
						{
							geometry::Triangle triangle(
								*(this->nodes[0]),
								*(this->nodes[1]),
								*(this->nodes[2]));

							triangle.SetIntegrationLaw(4);
							Tensor2Rank3D res = triangle.SolveIntegral(MassMatrix);
							local_matix.A[_bf_local_id_I][_bf_local_id_J] += res;
						}*/
					}
				}
				/*if (this->GetIdDomain() == 2)
					printf_s("\nThis elem[%d] is singular\n", this->GetSelfGlobalId());*/
			}
			catch (const std::exception&)
			{
				printf_s("Error: XFEM/MultiXFEM_Grid.h/void SolveLocalMatrix(DenseMatrix<Tensor2Rank3D> &local_matix, std::function<std::vector<std::vector<double>>(Point)> &koefD)\n");
			}
		}
	};
	class BoundaryVertex_forMech :
		public geometry::Vertex,
		//public topology::Vertex<topology::lower::EmptyElement, topology::upper::Segment>,
		public topology::Vertex<topology::lower::EmptyElement, topology::upper::Tetrahedron>,
		public functional::Shape<int, Point<double>>
	{
	public:
		std::function<Point<double>(Point<bool>&)> boundary_value;

		BoundaryVertex_forMech() { return; };
		~BoundaryVertex_forMech() { return; };

		void SetTypeBoundary(int type)
		{
			this->type_boundary = type;
		}
		int GetTypeBoundary()
		{
			return this->type_boundary;
		}

	private:
		int type_boundary;

	};
	class BoundaryFace_forMech :
		public geometry::Triangle,
		public topology::Triangle<topology::lower::Vertex, topology::upper::Tetrahedron>,
		public functional::Shape<int, Point<double>>
	{
	public:
		std::function<Point<double>(Point<double>)> boundary_value;

		BoundaryFace_forMech() { return; };
		~BoundaryFace_forMech() { return; };

		void SetTypeBoundary(int type)
		{
			this->type_boundary = type;
		}
		int GetTypeBoundary()
		{
			return this->type_boundary;
		}

		void SolveLocalBoundaryVector(std::vector<Point<double>> &local_vector, std::function<Point<double>(Point<double>)> &boundary_value)
		{
			try {
				local_vector.resize(this->GetDOFsCount());
				this->SetIntegrationLaw(4);

				for (int _bf_local_id_I = 0; _bf_local_id_I < this->GetDOFsCount(); _bf_local_id_I++)
				{
					auto basis_function_I = this->GetBasisFunctionInLocalID(_bf_local_id_I);
					
					std::function<Point<double>(Point<double>)> ValueInPoint = [&basis_function_I, &boundary_value]
					(Point<double> X) -> Point<double>
					{
						Point<double> result;
						auto boundary = boundary_value(X);
						auto func = (*basis_function_I)(X);
						result.x = func.x * boundary.x;
						result.y = func.y * boundary.y;
						result.z = func.z * boundary.z;

						/*result.x = 0.;
						result.y = 0.;
						result.z = 1.;*/
						return result;
					};

					double V = this->GetVolume();
					local_vector[_bf_local_id_I] = this->SolveIntegral(ValueInPoint);

				}
			}
			catch (const std::exception&)
			{
				printf_s("Error: XFEM/MultiXFEM_Grid.h/void SolveLocalMatrix(DenseMatrix<Tensor2Rank3D> &local_matix, std::function<std::vector<std::vector<double>>(Point)> &koefD)\n");
			}
		}

	private:
		int type_boundary;
	};

	class Grid : public geometry::Grid<FiniteElement_forMech>
	{
		int DOFs_count;
		//std::vector<topology::Triangle<topology::lower::Segment, topology::upper::Tetrahedron>> faces;
		//std::vector<topology::Segment<topology::lower::Vertex, topology::upper::Triangle>> edges; // Segments
		//std::vector<topology::Vertex<topology::lower::EmptyElement, topology::upper::Segment>> vertexes; //Vertexes
		std::vector<topology::Vertex<topology::lower::EmptyElement, topology::upper::Tetrahedron>> vertexes; //Vertexes

	public:
		std::vector<geometry::Crack> cracks;
		std::vector<BoundaryVertex_forMech> boundary_vertexes;
		std::vector<BoundaryFace_forMech> boundary_faces;

		std::vector<int> accordance_DOF_and_vertex;

		Grid()
		{
			DOFs_count = 0;
		};

		std::vector<int>* GetElementDOFs(int id_element)
		{
			try {
				return this->GetElement(id_element)->GetElementDOFs();
			}
			catch (const std::exception&)
			{
				printf_s("Error: MultiXFEM/MultiXFEM_Grid.h/std::vector<int>* GetElementDOFs(int id_element)\n");
			}

		}
		int GetDOFsCount()
		{
			return this->DOFs_count;
		}
		void SetDOFsCount(int count)
		{
			this->DOFs_count = count;
			this->accordance_DOF_and_vertex.resize(count);
			for (int i = 0; i < this->accordance_DOF_and_vertex.size(); i++)
			{
				this->accordance_DOF_and_vertex[i] = i;
			}
		}
		void AddDOF_ToEnd(int id_base_vertex)
		{
			this->DOFs_count++;
			this->accordance_DOF_and_vertex.push_back(id_base_vertex);
		}


		/*void CreateTopology()
		{
			try {
				struct face_temp
				{
					int num[3];
					int elem;
					int local_num;

				public:
					face_temp()
					{
						elem = 0;
						local_num = 0;
						num[0] = 0;
						num[1] = 0;
						num[2] = 0;
					}
					face_temp(int el, int local, int n1, int n2, int n3)
					{
						elem = el;
						local_num = local;
						num[0] = n1;
						num[1] = n2;
						num[2] = n3;

						std::vector<int> n(3);
						n[0] = n1; n[1] = n2; n[2] = n3;
						math::MakeQuickSort(n, 0, int(n.size()) - 1);
						num[0] = n[0];
						num[1] = n[1];
						num[2] = n[2];
					}

					void operator= (face_temp A)
					{
						elem = A.elem;
						local_num = A.local_num;
						num[0] = A.num[0];
						num[1] = A.num[1];
						num[2] = A.num[2];
					}
					bool operator< (face_temp A)
					{
						for (int i = 0; i < 3; i++)
						{
							if (num[i] < A.num[i]) return true;
							if (num[i] > A.num[i]) return false;
						}
						return false;
					}

					bool operator>(face_temp A)
					{
						for (int i = 0; i < 3; i++)
						{
							if (num[i] > A.num[i]) return true;
							if (num[i] < A.num[i]) return false;
						}
						return false;
					}
					bool operator!= (face_temp A)
					{
						for (int i = 0; i < 3; i++)
						{
							if (num[i] != A.num[i]) return true;
						}
						return false;
					}
				};
				struct face
				{
					int num[4];
					std::vector <int> elem;
					std::vector <int> local_num;

				public:
					face(face_temp A)
					{
						elem.push_back(A.elem);
						local_num.push_back(A.local_num);
						num[0] = A.num[0];
						num[1] = A.num[1];
						num[2] = A.num[2];
						num[3] = A.num[3];
					}
					face()
					{
						num[0] = 0;
						num[1] = 0;
						num[2] = 0;
						num[3] = 0;
					}
					void set_elem(face_temp A)
					{
						elem.push_back(A.elem);
						local_num.push_back(A.local_num);
					}
				};
				struct edge_temp
				{
					int num1;
					int num2;
					int element;
					int face;
					int local_in_element;
					int local_in_face;
					std::vector<double> test_ABCD;

				public:
					edge_temp()
					{
						/elem = 0;
						//local_num = 0;
						num1 = 0;
						num2 = 0;
					}
					edge_temp(int element, int local_in_element,
						int face, int local_in_face,
						int n1, int n2)
					{
						this->element = element;
						this->local_in_element = local_in_element;
						this->face = face;
						this->local_in_face = local_in_face;
						if (n1 < n2)
						{
							num1 = n1;
							num2 = n2;
						}
						else {
							num1 = n2;
							num2 = n1;
						}
					}

					void operator= (edge_temp A)
					{
						this->element = A.element;
						this->local_in_element = A.local_in_element;
						this->face = A.face;
						this->local_in_face = A.local_in_face;

						num1 = A.num1;
						num2 = A.num2;
					}
					bool operator< (edge_temp A)
					{
						if (num1 < A.num1) return true;
						if (num1 == A.num1)
						{
							if (num2 < A.num2) return true;
							else return false;
						}

						return false;
					}

					bool operator>(edge_temp A)
					{
						if (num1 > A.num1) return true;
						if (num1 == A.num1)
						{
							if (num2 > A.num2) return true;
							else return false;
						}

						return false;
					}
					bool operator!= (edge_temp A)
					{
						if (num1 != A.num1) return true;
						if (num1 == A.num1)
							if (num2 != A.num2) return true;

						return false;
					}

					bool operator== (edge_temp A)
					{
						if (num1 == A.num1 && num2 == A.num2) return true;

						return false;
					}
				};
				struct edge
				{
					int num1;
					int num2;
					std::vector <int> element, local_in_element;
					std::vector <int> face, local_in_face;

				public:
					edge(edge_temp A)
					{
						//element.push_back(A.element);
						//local_in_element.push_back(A.local_in_element);
						face.push_back(A.face);
						local_in_face.push_back(A.local_in_face);
						num1 = A.num1;
						num2 = A.num2;
					}
					edge()
					{
						num1 = 0;
						num2 = 0;
					}
					void set_elem(edge_temp A)
					{
						//element.push_back(A.element);
						//local_in_element.push_back(A.local_in_element);
						face.push_back(A.face);
						local_in_face.push_back(A.local_in_face);
					}
				};
				struct node_temp
				{
					int num;
					int element, local_in_element;
					int face, local_in_face;
					int edge, local_in_edge;
				public:
					node_temp()
					{
						element = 0; local_in_element=0;
						face = 0; local_in_face = 0;
						edge = 0; local_in_edge=0;

						num = 0;
					}
					node_temp(int element, int local_in_element,
						int face, int local_in_face,
						int edge, int local_in_edge, int n1)
					{
						this->element = element; this->local_in_element = local_in_element;
						this->face = face; this->local_in_face = local_in_face;
						this->edge = edge; this->local_in_edge = local_in_edge;

						num = n1;
					}

					void operator= (node_temp A)
					{
						this->element = A.element; this->local_in_element = A.local_in_element;
						this->face = A.face; this->local_in_face = A.local_in_face;
						this->edge = A.edge; this->local_in_edge = A.local_in_edge;

						num = A.num;
					}
					bool operator< (node_temp A)
					{
						if (num < A.num) return true;
						if (num > A.num) return false;
						return false;
					}

					bool operator>(node_temp A)
					{
						if (num > A.num) return true;
						if (num < A.num) return false;
						return false;
					}
					bool operator!= (node_temp A)
					{
						if (num != A.num) return true;
						return false;
					}

					bool operator == (node_temp A)
					{
						if (num == A.num) return true;

						return false;
					}
				};
				struct node
				{
					int num;

					std::vector<int> element, local_in_element;
					std::vector<int> face, local_in_face;
					std::vector<int> edge, local_in_edge;
				public:
					node(node_temp A)
					{
						//element.push_back(A.element);
						//local_in_element.push_back(A.local_in_element);
						//face.push_back(A.face);
						//local_in_face.push_back(A.local_in_face);
						edge.push_back(A.edge);
						local_in_edge.push_back(A.local_in_edge);

						num = A.num;
					}
					node()
					{
						num = 0;
					}
					void set_elem(node_temp A)
					{
						//element.push_back(A.element);
						//local_in_element.push_back(A.local_in_element);
						//face.push_back(A.face);
						//local_in_face.push_back(A.local_in_face);
						edge.push_back(A.edge);
						local_in_edge.push_back(A.local_in_edge);

					}
				};

				std::vector <edge_temp> e_temp;
				std::vector <face_temp> f_temp;
				std::vector <node_temp> n_temp;
				std::vector <edge> edges_vector;
				std::vector <face> faces_vector;
				std::vector <node> nodes_vector;
				topology::lower::Tetrahedron tmp_3D;
				topology::lower::Triangle tmp_2D;
				topology::lower::Segment tmp_1D;
				topology::lower::Vertex tmp_0D;
				std::vector<std::vector<int>> face_pattern, edge_pattern, node_pattern;
				tmp_3D.GetLowerElementPatternInLocal(face_pattern);
				tmp_2D.GetLowerElementPatternInLocal(edge_pattern);
				tmp_1D.GetLowerElementPatternInLocal(node_pattern);
				f_temp.reserve(face_pattern.size()*this->GetElementsCount());
				e_temp.reserve(edge_pattern.size()* face_pattern.size()*this->GetElementsCount());
				n_temp.reserve(node_pattern.size() * edge_pattern.size()* face_pattern.size()*this->GetElementsCount());

				auto localRemovalOfDuplication_e = [](std::vector <edge_temp> &A, std::vector <edge> &B)
				{
					int k = 0, j = 0;
					B.push_back(edge(A[0]));
					for (int i = 1; i < A.size(); i++)
					{
						if (A[j] != A[i]) //не повтор
						{
							j = i;
							B.push_back(edge(A[i]));
							k++;
						}
						else {
							B[k].set_elem(A[i]);
						}
					}
				};
				auto localRemovalOfDuplication_f = [](std::vector <face_temp> &A, std::vector <face> &B)
				{
					int k = 0, j = 0;
					B.push_back(face(A[0]));
					for (int i = 1; i < A.size(); i++)
					{
						if (A[j] != A[i]) //не повтор
						{
							j = i;
							B.push_back(face(A[i]));
							k++;
						}
						else {
							B[k].set_elem(A[i]);
						}
					}
				};
				auto localRemovalOfDuplication_n = [](std::vector <node_temp> &A, std::vector <node> &B)
				{
					int k = 0, j = 0;
					B.push_back(node(A[0]));
					for (int i = 1; i < A.size(); i++)
					{
						if (A[j] != A[i]) //не повтор
						{
							j = i;
							B.push_back(node(A[i]));
							k++;
						}
						else {
							B[k].set_elem(A[i]);
						}
					}
				};


				for (int id_element = 0; id_element < this->GetElementsCount(); id_element++)
				{
					auto element3D = this->GetElement(id_element);

					for (int id_face = 0; id_face < tmp_3D.GetLowerElementCount(); id_face++)
					{
						int _id_elem = id_element;
						int _id_face_in_elem = id_face;
						f_temp.push_back(face_temp(_id_elem, _id_face_in_elem,
							element3D->GetIdNode(face_pattern[id_face][0]),
							element3D->GetIdNode(face_pattern[id_face][1]),
							element3D->GetIdNode(face_pattern[id_face][2])));
					}
				}
				math::MakeQuickSort(f_temp);
				localRemovalOfDuplication_f(f_temp, faces_vector);

				for (int id_face = 0; id_face < faces_vector.size(); id_face++)
				{
					auto _face = &(faces_vector[id_face]);
					int _id_elem = -1;
					int _id_face_in_elem = -1;

					for (int id_edge = 0; id_edge < tmp_2D.GetLowerElementCount(); id_edge++)
					{
						int _id_face = id_face;
						int _id_edge_in_face = id_edge;
						e_temp.push_back(edge_temp(_id_elem, _id_face_in_elem,
							_id_face, _id_edge_in_face,
							_face->num[edge_pattern[id_edge][0]],
							_face->num[edge_pattern[id_edge][1]]));
					}
				}
				math::MakeQuickSort(e_temp);
				localRemovalOfDuplication_e(e_temp, edges_vector);

				for (int id_edge = 0; id_edge < edges_vector.size(); id_edge++)
				{
					auto _edge = &(edges_vector[id_edge]);
					int _id_elem = -1;
					int _id_face_in_elem = -1;
					int _id_face = -1;
					int _id_edge_in_face = -1;

					{
						int _id_edge = id_edge;
						int _id_vertex_in_edge =1;
						n_temp.push_back(node_temp(_id_elem, _id_face_in_elem,
							_id_face, _id_edge_in_face,
							_id_edge, _id_vertex_in_edge,
							_edge->num1));

						_id_vertex_in_edge = 1;
						n_temp.push_back(node_temp(_id_elem, _id_face_in_elem,
							_id_face, _id_edge_in_face,
							_id_edge, _id_vertex_in_edge,
							_edge->num2));

					}
				}
				math::MakeQuickSort(n_temp);
				localRemovalOfDuplication_n(n_temp, nodes_vector);

				this->faces.resize(faces_vector.size());
				this->edges.resize(edges_vector.size());
				this->vertexes.resize(nodes_vector.size());
				for (int id_face = 0; id_face < faces_vector.size(); id_face++)
				{
					this->faces[id_face].SetUpperElementCount((int)faces_vector[id_face].elem.size());
				}
				for (int id_edge = 0; id_edge < edges_vector.size(); id_edge++)
				{
					this->edges[id_edge].SetUpperElementCount((int)edges_vector[id_edge].face.size());
				}
				for (int id_node = 0; id_node < nodes_vector.size(); id_node++)
				{
					this->vertexes[id_node].SetUpperElementCount((int)nodes_vector[id_node].edge.size());
				}

				//add Upper and Lower Elements
				for (int id_element = 0; id_element < this->GetElementsCount(); id_element++)
				{
					this->GetElement(id_element)->SetSelfGlobalId(id_element);
				}
				for (int id_face = 0; id_face < faces_vector.size(); id_face++)
				{
					this->faces[id_face].SetSelfGlobalId(id_face);
				}
				for (int id_edge = 0; id_edge < edges_vector.size(); id_edge++)
				{
					this->edges[id_edge].SetSelfGlobalId(id_edge);
				}
				for (int id_node = 0; id_node < nodes_vector.size(); id_node++)
				{
					this->vertexes[id_node].SetSelfGlobalId(id_node);
				}

				for (int id_face = 0; id_face < faces_vector.size(); id_face++)
				{
					for (int id_volume = 0; id_volume < faces_vector[id_face].elem.size(); id_volume++)
					{
						this->faces[id_face].SetUpperElement(id_volume, this->GetElement(faces_vector[id_face].elem[id_volume]));
						this->GetElement(faces_vector[id_face].elem[id_volume])->SetLowerElement(faces_vector[id_face].local_num[id_volume], &(this->faces[id_face]));
					}
				}
				for (int id_edge = 0; id_edge < edges_vector.size(); id_edge++)
				{
					for (int id_face = 0; id_face < edges_vector[id_edge].face.size(); id_face++)
					{
						this->edges[id_edge].SetUpperElement(id_face, &(this->faces[edges_vector[id_edge].face[id_face]]));
						this->faces[edges_vector[id_edge].face[id_face]].SetLowerElement(edges_vector[id_edge].local_in_face[id_face], &(this->edges[id_edge]));
					}
				}
				for (int id_node = 0; id_node < nodes_vector.size(); id_node++)
				{
					for (int id_edge = 0; id_edge < nodes_vector[id_node].edge.size(); id_edge++)
					{
						this->vertexes[id_node].SetUpperElement(id_edge, &(this->edges[nodes_vector[id_node].edge[id_edge]]));
						this->edges[nodes_vector[id_node].edge[id_edge]].SetLowerElement(nodes_vector[id_node].local_in_edge[id_edge], &(this->vertexes[id_node]));
					}
				}
			}
			catch (const std::exception&)
			{
				printf_s("Error: MultiXFEM/MultiXFEM_Grid.h/void CreateTopology()\n");
			}
		}*/
		void CreateTopology()
		{
			try {
				struct node_temp
				{
					int num;
					int element, local_in_element;
					int face, local_in_face;
					int edge, local_in_edge;
				public:
					node_temp()
					{
						element = 0; local_in_element = 0;
						face = 0; local_in_face = 0;
						edge = 0; local_in_edge = 0;

						num = 0;
					}
					node_temp(int element, int local_in_element,
						int face, int local_in_face,
						int edge, int local_in_edge, int n1)
					{
						this->element = element; this->local_in_element = local_in_element;
						this->face = face; this->local_in_face = local_in_face;
						this->edge = edge; this->local_in_edge = local_in_edge;

						num = n1;
					}

					void operator= (node_temp A)
					{
						this->element = A.element; this->local_in_element = A.local_in_element;
						this->face = A.face; this->local_in_face = A.local_in_face;
						this->edge = A.edge; this->local_in_edge = A.local_in_edge;

						num = A.num;
					}
					bool operator< (node_temp A)
					{
						if (num < A.num) return true;
						if (num > A.num) return false;
						return false;
					}

					bool operator>(node_temp A)
					{
						if (num > A.num) return true;
						if (num < A.num) return false;
						return false;
					}
					bool operator!= (node_temp A)
					{
						if (num != A.num) return true;
						return false;
					}

					bool operator == (node_temp A)
					{
						if (num == A.num) return true;

						return false;
					}
				};
				struct node
				{
					int num;

					std::vector<int> element, local_in_element;
					std::vector<int> face, local_in_face;
					std::vector<int> edge, local_in_edge;
				public:
					node(node_temp A)
					{
						//element.push_back(A.element);
						//local_in_element.push_back(A.local_in_element);
						//face.push_back(A.face);
						//local_in_face.push_back(A.local_in_face);
						edge.push_back(A.edge);
						local_in_edge.push_back(A.local_in_edge);

						num = A.num;
					}
					node()
					{
						num = 0;
					}
					void set_elem(node_temp A)
					{
						//element.push_back(A.element);
						//local_in_element.push_back(A.local_in_element);
						//face.push_back(A.face);
						//local_in_face.push_back(A.local_in_face);
						edge.push_back(A.edge);
						local_in_edge.push_back(A.local_in_edge);

					}
				};

				std::vector <node_temp> n_temp;
				std::vector <node> nodes_vector;
				topology::lower::TetrahedronToVertex tmp_3D;
				topology::lower::Vertex tmp_0D;
				std::vector<std::vector<int>> node_pattern;
				tmp_3D.GetLowerElementPatternInLocal(node_pattern);
				n_temp.reserve(node_pattern.size() * this->GetElementsCount());

				auto localRemovalOfDuplication_n = [](std::vector <node_temp> &A, std::vector <node> &B)
				{
					int k = 0, j = 0;
					B.push_back(node(A[0]));
					for (int i = 1; i < A.size(); i++)
					{
						if (A[j] != A[i]) //не повтор
						{
							j = i;
							B.push_back(node(A[i]));
							k++;
						}
						else {
							B[k].set_elem(A[i]);
						}
					}
				};


				for (int id_element = 0; id_element < this->GetElementsCount(); id_element++)
				{
					auto element3D = this->GetElement(id_element);

					for (int id_vertex = 0; id_vertex < tmp_3D.GetLowerElementCount(); id_vertex++)
					{
						int _id_elem = -1;
						int _id_face_in_elem = -1;
						int _id_face = -1;
						int _id_edge_in_face = -1;
						int _id_edge = id_element;
						int _id_vertex_in_edge = id_vertex;
						n_temp.push_back(node_temp(_id_elem, _id_face_in_elem,
							_id_face, _id_edge_in_face,
							_id_edge, _id_vertex_in_edge,
							element3D->GetIdNode(node_pattern[id_vertex][0])));
					}
				}
				math::MakeQuickSort(n_temp);
				localRemovalOfDuplication_n(n_temp, nodes_vector);

				this->vertexes.resize(nodes_vector.size());
				for (int id_node = 0; id_node < nodes_vector.size(); id_node++)
				{
					this->vertexes[id_node].SetUpperElementCount((int)nodes_vector[id_node].edge.size());
				}

				//add Upper and Lower Elements
				for (int id_element = 0; id_element < this->GetElementsCount(); id_element++)
				{
					this->GetElement(id_element)->SetSelfGlobalId(id_element);
				}
				for (int id_node = 0; id_node < nodes_vector.size(); id_node++)
				{
					this->vertexes[id_node].SetSelfGlobalId(id_node);
				}

				for (int id_node = 0; id_node < nodes_vector.size(); id_node++)
				{
					for (int id_volume = 0; id_volume < nodes_vector[id_node].edge.size(); id_volume++)
					{
						this->vertexes[id_node].SetUpperElement(id_volume, this->GetElement(nodes_vector[id_node].edge[id_volume]));
						this->GetElement(nodes_vector[id_node].edge[id_volume])->SetLowerElement(nodes_vector[id_node].local_in_edge[id_volume], &(this->vertexes[id_node]));
					}
				}
			}
			catch (const std::exception&)
			{
				printf_s("Error: MultiXFEM/MultiXFEM_Grid.h/void CreateTopology()\n");
			}
		}

		void CreateExtendedFunctionalSpace_ForCrack()
		{
			//return;
			for (int id_crack = 0; id_crack < this->cracks.size(); id_crack++)
			{
				auto current_crack = &(this->cracks[id_crack]);

				//добавка "обогащающих" функций/степеней свободы
				//разрывные
				std::vector <int> id_extended_disc_vertexes;
				if(true){
					double minH = 200;
					for (int i = 0; i < this->GetElementsCount(); i++)
					{
						bool transform = false;
						if (i % 1000 == 0)
						{
							printf("current elem is %d\r", i);
						}

						auto current_element = this->GetElement(i);

						//определение пересечени€ тетраэдра плоскостью трещины current_crack
						{
							bool find_include = false;
							//это что блин за условие?? а если трещина ѕЋќ— јя??
							/*for (int id_vertex = 0; id_vertex < current_element->GetNodesCount() && !find_include; id_vertex++)
							{
								Point<double> x = current_element->GetNode(id_vertex);
								
								if (current_crack->min_point <= x && x <= current_crack->max_point)
									find_include = true;
								
							}*/
							if (find_include || true)
							{
								Point<double> edge[6];
								int pattern_for_edge[6][2] = { { 0,1 },{ 0,2 },{ 0,3 },{ 1,2 },{ 1,3 },{ 2,3 } };
								int cross_section[6], SUMM = 0, num_of_tr[6];
								cross_section[0] = 0;
								cross_section[1] = 0;
								cross_section[2] = 0;
								cross_section[3] = 0;
								cross_section[4] = 0;
								cross_section[5] = 0;
								for (int tr = 0; tr < current_crack->triangles.size() && SUMM < 4; tr++)
								{
									double A, B, C, D;
									math::GetPlaneEquation(current_crack->triangles[tr].GetNode(0),
										current_crack->triangles[tr].GetNode(1),
										current_crack->triangles[tr].GetNode(2),
										A, B, C, D);

									for (int e = 0; e < 6; e++)
									{
										if (cross_section[e] == 0)
										{
											int test_for0 = (int)math::GetSignum(A*current_element->GetNode(pattern_for_edge[e][0]).x + B * current_element->GetNode(pattern_for_edge[e][0]).y + C * current_element->GetNode(pattern_for_edge[e][0]).z + D);
											int test_for1 = math::GetSignum(A*current_element->GetNode(pattern_for_edge[e][1]).x + B * current_element->GetNode(pattern_for_edge[e][1]).y + C * current_element->GetNode(pattern_for_edge[e][1]).z + D);
											if (test_for0 != 0 && test_for1 != 0 && test_for0 != test_for1)
											{
												edge[e] = math::GetIntersectionPointOfLineAndPlane(A, B, C, D, current_element->GetNode(pattern_for_edge[e][0]), current_element->GetNode(pattern_for_edge[e][1]));
												bool _1 = current_crack->triangles[tr].IsContainThePoint(edge[e]);
												bool _2 = current_element->IsContainThePoint(edge[e]);
												auto _3 = current_element->GetNode(pattern_for_edge[e][0]);
												auto _4 = current_element->GetNode(pattern_for_edge[e][1]);

												if (current_crack->triangles[tr].IsContainThePoint(edge[e]) &&
													current_element->IsContainThePoint(edge[e]) &&
													edge[e] != current_element->GetNode(pattern_for_edge[e][0]) &&
													edge[e] != current_element->GetNode(pattern_for_edge[e][1]))
												{
													cross_section[e] = 1;
													num_of_tr[e] = tr;
													SUMM++;
												}
											}
										}
									}
								}

								int cross_section_forNodes[4], SUMM_forNodes = 0, num_of_tr_forNodes[6];
								cross_section_forNodes[0] = 0;
								cross_section_forNodes[1] = 0;
								cross_section_forNodes[2] = 0;
								cross_section_forNodes[3] = 0;
								double _min = 1.0e+10;
								for (int tr = 0; tr < current_crack->triangles.size() && SUMM_forNodes < 4; tr++)
								{
									double A, B, C, D;
									math::GetPlaneEquation(current_crack->triangles[tr].GetNode(0),
										current_crack->triangles[tr].GetNode(1),
										current_crack->triangles[tr].GetNode(2),
										A, B, C, D);

									for (int n = 0; n < 4; n++)
									{
										Point<double> curr_node = current_element->GetNode(n);
										if (cross_section_forNodes[n] == 0)
										{
											double _tt = A * curr_node.x + B * curr_node.y + C * curr_node.z + D;
											int test = math::GetSignum(_tt);
											if (abs(_tt) < abs(_min))
											{
												_min = _tt;
												//printf("(%d:%d) %d - %.2e\n", i, tr, test, A*curr_node.x + B*curr_node.y + C*curr_node.z + D);
											}
											if (current_crack->triangles[tr].IsContainThePoint(curr_node) &&
												test == 0)
											{
												cross_section_forNodes[n] = 1;
												num_of_tr_forNodes[n] = tr;
												SUMM_forNodes++;
											}
										}
									}
								}

								//2 ребра и 1 узел
								auto _2edge_1node = [&](int node1, int node2, int node3, int node4, std::vector<int> test_tr, Point<double> e23, Point<double> e34)
								{
									double A, B, C, D;
									int signum[4];
									Point<double> test_vector;
									bool reverse = false;
									for (int t = 0; t < test_tr.size(); t++)
									{
										double _A, _B, _C, _D;
										math::GetPlaneEquation(current_crack->triangles[test_tr[t]].GetNodes(), _A, _B, _C, _D);
										test_vector += Point<double>(_A, _B, _C) / (double)test_tr.size();
									}
									math::GetPlaneEquation(
										e23,
										e34,
										current_element->GetNode(node1), A, B, C, D);
									if (Point<double>(A, B, C)*test_vector < 0)
									{
										math::GetPlaneEquation(
											e23,
											current_element->GetNode(node1),
											e34, A, B, C, D);
										reverse = true;
									}

									for (int node = 0; node < current_element->GetNodesCount(); node++)
									{
										signum[node] = math::GetSignum(A*current_element->GetNode(node).x + B * current_element->GetNode(node).y + C * current_element->GetNode(node).z + D);
									}
									if (signum[node1] == 0 && signum[node2] < 0 && signum[node3] > 0 && signum[node4] > 0 ||
										signum[node1] == 0 && signum[node2] > 0 && signum[node3] < 0 && signum[node4] < 0)
									{
										current_element->SubGrid_for_integration.xyz.resize(7);
										int pattern[3][4];
										int local_id_nodes[4] = { node1, node2, node3, node4 };
										int local_id_vertex[3][2] = { { node2, node3 },{ node3, node4 },{ node1, node1 } };

										current_element->SubGrid_for_integration.xyz[0] = this->GetCoordinateViaID(current_element->GetIdNode(node1));
										current_element->SubGrid_for_integration.xyz[1] = this->GetCoordinateViaID(current_element->GetIdNode(node2));
										current_element->SubGrid_for_integration.xyz[2] = this->GetCoordinateViaID(current_element->GetIdNode(node3));
										current_element->SubGrid_for_integration.xyz[3] = this->GetCoordinateViaID(current_element->GetIdNode(node4));

										current_element->SubGrid_for_integration.xyz[4] = e23;
										current_element->SubGrid_for_integration.xyz[5] = e34;
										current_element->SubGrid_for_integration.xyz[6] = this->GetCoordinateViaID(current_element->GetIdNode(node1));

										//edge 13
										if (current_element->GetIdNode(local_id_nodes[2]) < current_element->GetIdNode(local_id_nodes[3]))
										{
											pattern[0][0] = 2; pattern[0][1] = 3; pattern[0][2] = 4; pattern[0][3] = 6;
											pattern[1][0] = 3; pattern[1][1] = 4; pattern[1][2] = 5; pattern[1][3] = 6;
										}
										else {
											pattern[0][0] = 2; pattern[0][1] = 3; pattern[0][2] = 5; pattern[0][3] = 6;
											pattern[1][0] = 2; pattern[1][1] = 4; pattern[1][2] = 5; pattern[1][3] = 6;
										}
										pattern[2][0] = 1; pattern[2][1] = 4; pattern[2][2] = 5; pattern[2][3] = 0;

										//testing
										if (reverse)
										{
											local_id_vertex[0][0] = node2; local_id_vertex[0][1] = node3;
											local_id_vertex[1][0] = node1; local_id_vertex[1][1] = node1;
											local_id_vertex[2][0] = node3; local_id_vertex[2][1] = node4;

											current_element->SubGrid_for_integration.xyz[4] = e23;
											current_element->SubGrid_for_integration.xyz[5] = this->GetCoordinateViaID(current_element->GetIdNode(node1));
											current_element->SubGrid_for_integration.xyz[6] = e34;

											//edge 13
											if (current_element->GetIdNode(local_id_nodes[1]) < current_element->GetIdNode(local_id_nodes[3]))
											{
												pattern[0][0] = 2; pattern[0][1] = 3; pattern[0][2] = 4; pattern[0][3] = 5;
												pattern[1][0] = 3; pattern[1][1] = 4; pattern[1][2] = 6; pattern[1][3] = 5;
											}
											else {
												pattern[0][0] = 2; pattern[0][1] = 3; pattern[0][2] = 6; pattern[0][3] = 5;
												pattern[1][0] = 2; pattern[1][1] = 4; pattern[1][2] = 6; pattern[1][3] = 5;
											}
											pattern[2][0] = 1; pattern[2][1] = 4; pattern[2][2] = 6; pattern[2][3] = 0;
										}

										for (int p = 0; p < 3; p++)
										{
											std::vector<int> elem(4);
											for (int pp = 0; pp < 4; pp++)
											{
												elem[pp] = pattern[p][pp];
											}
											current_element->SubGrid_for_integration.nvtr.push_back(elem);
											current_element->SubGrid_for_integration.nvkat.push_back(0);
										}

										current_element->SubGrid_for_integr_byTriangles.xyz.resize(3);
										current_element->SubGrid_for_integr_byTriangles.xyz[0] = current_element->SubGrid_for_integration.xyz[4];
										current_element->SubGrid_for_integr_byTriangles.xyz[1] = current_element->SubGrid_for_integration.xyz[5];
										current_element->SubGrid_for_integr_byTriangles.xyz[2] = current_element->SubGrid_for_integration.xyz[6];
										std::vector<int> elem(3);
										elem[0] = 0; elem[1] = 1; elem[2] = 2;
										current_element->SubGrid_for_integr_byTriangles.nvtr.push_back(elem);
										current_element->SubGrid_for_integr_byTriangles.nvkat.push_back(id_crack);

										return true;
									}
									return false;
								};
								//1 ребро и 2 узла
								auto _1edge_2node = [&](int node1, int node2, int node3, int node4, std::vector<int> test_tr, Point<double> e34)
								{
									double A, B, C, D;
									int signum[4];
									Point<double> test_vector;
									bool reverse = false;
									for (int t = 0; t < test_tr.size(); t++)
									{
										double _A, _B, _C, _D;
										math::GetPlaneEquation(current_crack->triangles[test_tr[t]].GetNodes(), _A, _B, _C, _D);
										test_vector += Point<double>(_A, _B, _C) / (double)test_tr.size();
									}
									math::GetPlaneEquation(
										e34,
										this->GetCoordinateViaID(current_element->GetIdNode(node2)),
										this->GetCoordinateViaID(current_element->GetIdNode(node1)), A, B, C, D);
									if (Point<double>(A, B, C)*test_vector < 0)
									{
										math::GetPlaneEquation(
											e34,
											this->GetCoordinateViaID(current_element->GetIdNode(node1)),
											this->GetCoordinateViaID(current_element->GetIdNode(node2)), A, B, C, D);
										reverse = true;
									}
									for (int node = 0; node < current_element->GetNodesCount(); node++)
									{
										signum[node] = math::GetSignum(A*this->GetCoordinateViaID(current_element->GetIdNode(node)).x + B * this->GetCoordinateViaID(current_element->GetIdNode(node)).y + C * this->GetCoordinateViaID(current_element->GetIdNode(node)).z + D);
									}
									if (signum[node1] == 0 && signum[node2] == 0 && signum[node3] > 0 && signum[node4] < 0 ||
										signum[node1] == 0 && signum[node2] == 0 && signum[node3] < 0 && signum[node4] > 0)
									{
										current_element->SubGrid_for_integration.xyz.resize(7);
										int pattern[2][4];
										int local_id_nodes[4] = { node1, node2, node3, node4 };
										int local_id_vertex[3][2] = { { node3, node4 },{ node2, node2 },{ node1, node1 } };

										current_element->SubGrid_for_integration.xyz[0] = this->GetCoordinateViaID(current_element->GetIdNode(node1));
										current_element->SubGrid_for_integration.xyz[1] = this->GetCoordinateViaID(current_element->GetIdNode(node2));
										current_element->SubGrid_for_integration.xyz[2] = this->GetCoordinateViaID(current_element->GetIdNode(node3));
										current_element->SubGrid_for_integration.xyz[3] = this->GetCoordinateViaID(current_element->GetIdNode(node4));

										current_element->SubGrid_for_integration.xyz[4] = e34;
										current_element->SubGrid_for_integration.xyz[5] = this->GetCoordinateViaID(current_element->GetIdNode(node2));
										current_element->SubGrid_for_integration.xyz[6] = this->GetCoordinateViaID(current_element->GetIdNode(node1));
										pattern[0][0] = 3; pattern[0][1] = 4; pattern[0][2] = 1; pattern[0][3] = 0;
										pattern[1][0] = 2; pattern[1][1] = 4; pattern[1][2] = 5; pattern[1][3] = 6;

										//testing
										if (reverse == true)
										{
											local_id_vertex[0][0] = node3; local_id_vertex[0][1] = node4;
											local_id_vertex[1][0] = node1; local_id_vertex[1][1] = node1;
											local_id_vertex[2][0] = node2; local_id_vertex[2][1] = node2;

											current_element->SubGrid_for_integration.xyz[4] = e34;
											current_element->SubGrid_for_integration.xyz[5] = this->GetCoordinateViaID(current_element->GetIdNode(node1));
											current_element->SubGrid_for_integration.xyz[6] = this->GetCoordinateViaID(current_element->GetIdNode(node2));
											pattern[0][0] = 3; pattern[0][1] = 4; pattern[0][2] = 1; pattern[0][3] = 0;
											pattern[1][0] = 2; pattern[1][1] = 4; pattern[1][2] = 5; pattern[1][3] = 6;
										}


										for (int p = 0; p < 2; p++)
										{
											std::vector<int> elem(4);
											for (int pp = 0; pp < 4; pp++)
											{
												elem[pp] = pattern[p][pp];
											}
											current_element->SubGrid_for_integration.nvtr.push_back(elem);
											current_element->SubGrid_for_integration.nvkat.push_back(current_element->GetSelfGlobalId());
										}

										current_element->SubGrid_for_integr_byTriangles.xyz.resize(3);
										current_element->SubGrid_for_integr_byTriangles.xyz[0] = current_element->SubGrid_for_integration.xyz[4];
										current_element->SubGrid_for_integr_byTriangles.xyz[1] = current_element->SubGrid_for_integration.xyz[5];
										current_element->SubGrid_for_integr_byTriangles.xyz[2] = current_element->SubGrid_for_integration.xyz[6];
										std::vector<int> elem(3);
										elem[0] = 0; elem[1] = 1; elem[2] = 2;
										current_element->SubGrid_for_integr_byTriangles.nvtr.push_back(elem);
										current_element->SubGrid_for_integr_byTriangles.nvkat.push_back(id_crack);

										return true;
									}
									return false;
								};
								//0 ребро и 3 узла
								auto _0edge_3node = [&](int node1, int node2, int node3, int node4, std::vector<int> test_tr)
								{
									double A, B, C, D;
									int signum[4];
									Point<double> test_vector;
									bool reverse = false;
									for (int t = 0; t < test_tr.size(); t++)
									{
										double _A, _B, _C, _D;
										math::GetPlaneEquation(current_crack->triangles[test_tr[t]].GetNodes(), _A, _B, _C, _D);
										test_vector += Point<double>(_A, _B, _C) / (double)test_tr.size();
									}
									math::GetPlaneEquation(
										this->GetCoordinateViaID(current_element->GetIdNode(node1)),
										this->GetCoordinateViaID(current_element->GetIdNode(node2)),
										this->GetCoordinateViaID(current_element->GetIdNode(node3)), A, B, C, D);
									if (Point<double>(A, B, C)*test_vector < 0)
									{
										math::GetPlaneEquation(
											this->GetCoordinateViaID(current_element->GetIdNode(node1)),
											this->GetCoordinateViaID(current_element->GetIdNode(node3)),
											this->GetCoordinateViaID(current_element->GetIdNode(node2)), A, B, C, D);
										reverse = true;
									}

									for (int node = 0; node < current_element->GetNodesCount(); node++)
									{
										signum[node] = math::GetSignum(A*this->GetCoordinateViaID(current_element->GetIdNode(node)).x + B * this->GetCoordinateViaID(current_element->GetIdNode(node)).y + C * this->GetCoordinateViaID(current_element->GetIdNode(node)).z + D);
									}
									if (signum[node1] == 0 && signum[node2] == 0 && signum[node3] == 0 && signum[node4] < 0 ||
										signum[node1] == 0 && signum[node2] == 0 && signum[node3] == 0 && signum[node4] > 0)
									{
										current_element->SubGrid_for_integration.xyz.resize(7);
										int pattern[1][4];
										int local_id_nodes[4] = { node1, node2, node3, node4 };
										int local_id_vertex[3][2] = { { node1, node1 },{ node2, node2 },{ node3, node3 } };

										current_element->SubGrid_for_integration.xyz[0] = this->GetCoordinateViaID(current_element->GetIdNode(node1));
										current_element->SubGrid_for_integration.xyz[1] = this->GetCoordinateViaID(current_element->GetIdNode(node2));
										current_element->SubGrid_for_integration.xyz[2] = this->GetCoordinateViaID(current_element->GetIdNode(node3));
										current_element->SubGrid_for_integration.xyz[3] = this->GetCoordinateViaID(current_element->GetIdNode(node4));

										if (!reverse)
										{
											current_element->SubGrid_for_integration.xyz[4] = this->GetCoordinateViaID(current_element->GetIdNode(node1));
											current_element->SubGrid_for_integration.xyz[5] = this->GetCoordinateViaID(current_element->GetIdNode(node2));
											current_element->SubGrid_for_integration.xyz[6] = this->GetCoordinateViaID(current_element->GetIdNode(node3));
											pattern[0][0] = 3; pattern[0][1] = 4; pattern[0][2] = 5; pattern[0][3] = 6;
										}
										else {
											local_id_vertex[0][0] = node1; local_id_vertex[0][1] = node1;
											local_id_vertex[1][0] = node3; local_id_vertex[1][1] = node3;
											local_id_vertex[2][0] = node2; local_id_vertex[2][1] = node2;

											current_element->SubGrid_for_integration.xyz[4] = this->GetCoordinateViaID(current_element->GetIdNode(node1));
											current_element->SubGrid_for_integration.xyz[5] = this->GetCoordinateViaID(current_element->GetIdNode(node3));
											current_element->SubGrid_for_integration.xyz[6] = this->GetCoordinateViaID(current_element->GetIdNode(node2));
											pattern[0][0] = 3; pattern[0][1] = 4; pattern[0][2] = 5; pattern[0][3] = 6;
										}

										for (int p = 0; p < 1; p++)
										{
											std::vector<int> elem(4);
											for (int pp = 0; pp < 4; pp++)
											{
												elem[pp] = pattern[p][pp];
											}
											current_element->SubGrid_for_integration.nvtr.push_back(elem);
											current_element->SubGrid_for_integration.nvkat.push_back(0);
										}

										current_element->SubGrid_for_integr_byTriangles.xyz.resize(3);
										current_element->SubGrid_for_integr_byTriangles.xyz[0] = current_element->SubGrid_for_integration.xyz[4];
										current_element->SubGrid_for_integr_byTriangles.xyz[1] = current_element->SubGrid_for_integration.xyz[5];
										current_element->SubGrid_for_integr_byTriangles.xyz[2] = current_element->SubGrid_for_integration.xyz[6];
										std::vector<int> elem(3);
										elem[0] = 0; elem[1] = 1; elem[2] = 2;
										current_element->SubGrid_for_integr_byTriangles.nvtr.push_back(elem);
										current_element->SubGrid_for_integr_byTriangles.nvkat.push_back(id_crack);

										return true;
									}
									return false;
								};
								if (SUMM == 2 && SUMM_forNodes == 1)
								{
									int pattern_edge[12][2] = {
										{ 0,1 },{ 0,2 },{ 0,3 },
										{ 0,4 },{ 1,2 },{ 1,3 },
										{ 1,5 },{ 2,4 },{ 2,5 },
										{ 3,4 },{ 3,5 },{ 4,5 } };
									int pattern_node[12][4] = {
										{ 3, 1, 2, 0 },{ 2, 1, 0, 3 },{ 3, 1, 2, 0 },
										{ 2, 1, 0, 3 },{ 1, 0, 2, 3 },{ 3, 1, 2, 0 },
										{ 1, 0, 2, 3 },{ 2, 1, 0, 3 },{ 1, 0, 2, 3 },
										{ 0, 1, 2, 3 },{ 0, 1, 2, 3 },{ 0, 1, 2, 3 } };
									for (int mm = 0; mm < 12; mm++)
									{
										if (cross_section[pattern_edge[mm][0]] == 1 &&
											cross_section[pattern_edge[mm][1]] == 1 &&
											cross_section_forNodes[pattern_node[mm][0]] == 1)
										{
											std::vector<int> tr(3);
											tr[0] = num_of_tr_forNodes[pattern_node[mm][0]];
											tr[1] = num_of_tr[pattern_edge[mm][0]];
											tr[2] = num_of_tr[pattern_edge[mm][1]];
											bool update = _2edge_1node(
												pattern_node[mm][0],
												pattern_node[mm][1],
												pattern_node[mm][2],
												pattern_node[mm][3],
												tr,
												edge[pattern_edge[mm][0]],
												edge[pattern_edge[mm][1]]);

											if (update == true)
											{
												transform = true;
												//current_element->SetIdDomain(1);
												id_extended_disc_vertexes.push_back(current_element->GetIdNode(0));
												id_extended_disc_vertexes.push_back(current_element->GetIdNode(1));
												id_extended_disc_vertexes.push_back(current_element->GetIdNode(2));
												id_extended_disc_vertexes.push_back(current_element->GetIdNode(3));
												//n_disc_temp.push_back(node_temp(i, 0, c, 0, current_element->GetIdNode(0)));
												//n_disc_temp.push_back(node_temp(i, 0, c, 1, current_element->GetIdNode(1)));
												//n_disc_temp.push_back(node_temp(i, 0, c, 2, current_element->GetIdNode(2)));
												//n_disc_temp.push_back(node_temp(i, 0, c, 3, current_element->GetIdNode(3)));
											}
										}
									}
								}
								if (SUMM == 1 && SUMM_forNodes == 2)
								{
									int pattern_node[6][4] = {
										{ 2, 3, 0, 1 },{ 1, 3, 0, 2 },{ 1, 2, 0, 3 },
										{ 0, 3, 1, 2 },{ 0, 2, 1, 3 },{ 0, 1, 2, 3 } };
									for (int mm = 0; mm < 6; mm++)
									{
										if (cross_section[mm] == 1)
										{
											std::vector<int> tr(3);
											tr[0] = num_of_tr_forNodes[pattern_node[mm][0]];
											tr[1] = num_of_tr_forNodes[pattern_node[mm][1]];
											tr[2] = num_of_tr[mm];
											bool update = _1edge_2node(
												pattern_node[mm][0],
												pattern_node[mm][1],
												pattern_node[mm][2],
												pattern_node[mm][3],
												tr,
												edge[mm]);

											if (update == true)
											{
												transform = true;
												//n_disc_temp.push_back(node_temp(i, 0, c, 0, current_element->GetIdNode(0)));
												//n_disc_temp.push_back(node_temp(i, 0, c, 1, current_element->GetIdNode(1)));
												//n_disc_temp.push_back(node_temp(i, 0, c, 2, current_element->GetIdNode(2)));
												//n_disc_temp.push_back(node_temp(i, 0, c, 3, current_element->GetIdNode(3)));

												//current_element->SetIdDomain(2);
												id_extended_disc_vertexes.push_back(current_element->GetIdNode(0));
												id_extended_disc_vertexes.push_back(current_element->GetIdNode(1));
												id_extended_disc_vertexes.push_back(current_element->GetIdNode(2));
												id_extended_disc_vertexes.push_back(current_element->GetIdNode(3));
											}
										}
									}
								}
								if (SUMM == 0 && SUMM_forNodes == 3)
								{
									int pattern_node[4][4] = { { 0,1,2,3 },{ 0,1,3,2 },{ 0,2,3,1 },{ 1,2,3,0 } };
									for (int mm = 0; mm < 4; mm++)
									{
										if (cross_section_forNodes[pattern_node[mm][0]] == 1 &&
											cross_section_forNodes[pattern_node[mm][1]] == 1 &&
											cross_section_forNodes[pattern_node[mm][2]] == 1)
										{
											std::vector<int> tr(4);
											tr[0] = num_of_tr_forNodes[pattern_node[mm][0]];
											tr[1] = num_of_tr_forNodes[pattern_node[mm][1]];
											tr[2] = num_of_tr_forNodes[pattern_node[mm][2]];
											_0edge_3node(
												pattern_node[mm][0],
												pattern_node[mm][1],
												pattern_node[mm][2],
												pattern_node[mm][3],
												tr);

											//add_elem = true;
											transform = true;
											//current_element->set_id_domain(3);
											//n_disc_temp.push_back(node_temp(i, 0, c, pattern_node[mm][0], current_element->GetIdNode(pattern_node[mm][0])));
											//n_disc_temp.push_back(node_temp(i, 0, c, pattern_node[mm][1], current_element->GetIdNode(pattern_node[mm][1])));
											//n_disc_temp.push_back(node_temp(i, 0, c, pattern_node[mm][2], current_element->GetIdNode(pattern_node[mm][2])));

											//current_element->SetIdDomain(3);
											id_extended_disc_vertexes.push_back(current_element->GetIdNode(pattern_node[mm][0]));
											id_extended_disc_vertexes.push_back(current_element->GetIdNode(pattern_node[mm][1]));
											id_extended_disc_vertexes.push_back(current_element->GetIdNode(pattern_node[mm][2]));
										}
									}
								}

								//пересекает весь тетраэдр через 3 ребра
								auto _3edge = [&](int node1, int node2, int node3, int node4, std::vector<int> test_tr, Point<double> e12, Point<double> e13, Point<double> e14)
								{
									double A, B, C, D;
									int signum[4];
									Point<double> test_vector;
									bool reverse = false;
									for (int t = 0; t < test_tr.size(); t++)
									{
										double _A, _B, _C, _D;
										math::GetPlaneEquation(current_crack->triangles[test_tr[t]].GetNodes(), _A, _B, _C, _D);
										test_vector += Point<double>(_A, _B, _C) / (double)test_tr.size();
									}
									math::GetPlaneEquation(e12, e13, e14, A, B, C, D);
									if (Point<double>(A, B, C)*test_vector < 0)
									{
										math::GetPlaneEquation(e14, e13, e12, A, B, C, D);
										reverse = true;
									}
									/*printf("_A=%.2e _B=%.2e _C=%.2e _D=%.2e -> \nA=%.2e B=%.2e C=%.2e D=%.2e\n", _A, _B, _C, _D, A, B, C, D);
									printf("Node0(%.8e, %.8e, %.8e) \ne12(%.8e, %.8e, %.8e) e13(%.8e,%.8e,%.8e) e14(%.8e,%.8e,%.8e)\n", this->GetCoordinateViaID(current_element->GetIdNode(node1)).x, this->GetCoordinateViaID(current_element->GetIdNode(node1)).y, this->GetCoordinateViaID(current_element->GetIdNode(node1)).z,
									e12.x, e12.y, e12.z, e13.x, e13.y, e13.z, e14.x, e14.y, e14.z);
									printf("sign = %.2e\n", A*this->GetCoordinateViaID(current_element->GetIdNode(node1)).x + B*this->GetCoordinateViaID(current_element->GetIdNode(node1)).y + C*this->GetCoordinateViaID(current_element->GetIdNode(node1)).z + D);
									*/
									for (int node = 0; node < current_element->GetNodesCount(); node++)
									{
										signum[node] = math::GetSignum(A*this->GetCoordinateViaID(current_element->GetIdNode(node)).x + B * this->GetCoordinateViaID(current_element->GetIdNode(node)).y + C * this->GetCoordinateViaID(current_element->GetIdNode(node)).z + D);
									}
									if (signum[node1] > 0 && signum[node2] < 0 && signum[node3] < 0 && signum[node4] < 0 ||
										signum[node1] < 0 && signum[node2] > 0 && signum[node3] > 0 && signum[node4] > 0)
									{
										current_element->SubGrid_for_integration.xyz.resize(8);
										int pattern[9][4];
										int local_id_nodes[4] = { node1, node2, node3, node4 };
										int local_id_vertex[3][2] = { { node1, node2 },{ node1, node3 },{ node1, node4 } };
										//int local_id_transform_faces[3][4] = { { node1,node2,node4,node5},{ node1,node3,node4,node6},{ node2,node3,node5,node6} };

										current_element->SubGrid_for_integration.xyz[0] = this->GetCoordinateViaID(current_element->GetIdNode(node1));
										current_element->SubGrid_for_integration.xyz[1] = this->GetCoordinateViaID(current_element->GetIdNode(node2));
										current_element->SubGrid_for_integration.xyz[2] = this->GetCoordinateViaID(current_element->GetIdNode(node3));
										current_element->SubGrid_for_integration.xyz[3] = this->GetCoordinateViaID(current_element->GetIdNode(node4));

										pattern[6][0] = 1; pattern[6][1] = 2; pattern[6][2] = 3; pattern[6][3] = 7;
										pattern[7][0] = 0; pattern[7][1] = 4; pattern[7][2] = 5; pattern[7][3] = 6;
										pattern[8][0] = 7; pattern[8][1] = 4; pattern[8][2] = 5; pattern[8][3] = 6;
										//for_transform
										current_element->SubGrid_for_integration.xyz[7] = (e12 + e13 + e14 + current_element->SubGrid_for_integration.xyz[1] + current_element->SubGrid_for_integration.xyz[2] + current_element->SubGrid_for_integration.xyz[3]) / 6.;

										if (!reverse)
										{
											current_element->SubGrid_for_integration.xyz[4] = e12;
											current_element->SubGrid_for_integration.xyz[5] = e13;
											current_element->SubGrid_for_integration.xyz[6] = e14;

											//face 012
											if (current_element->GetIdNode(local_id_nodes[1]) < current_element->GetIdNode(local_id_nodes[2]))
											{
												pattern[0][0] = 2; pattern[0][1] = 4; pattern[0][2] = 5; pattern[0][3] = 7;
												pattern[1][0] = 1; pattern[1][1] = 2; pattern[1][2] = 4; pattern[1][3] = 7;
											}
											else {
												pattern[0][0] = 1; pattern[0][1] = 4; pattern[0][2] = 5; pattern[0][3] = 7;
												pattern[1][0] = 1; pattern[1][1] = 2; pattern[1][2] = 5; pattern[1][3] = 7;
											}
											//face 013
											if (current_element->GetIdNode(local_id_nodes[1]) < current_element->GetIdNode(local_id_nodes[3]))
											{
												pattern[2][0] = 3; pattern[2][1] = 4; pattern[2][2] = 6; pattern[2][3] = 7;
												pattern[3][0] = 1; pattern[3][1] = 3; pattern[3][2] = 4; pattern[3][3] = 7;
											}
											else {
												pattern[2][0] = 1; pattern[2][1] = 3; pattern[2][2] = 6; pattern[2][3] = 7;
												pattern[3][0] = 1; pattern[3][1] = 4; pattern[3][2] = 6; pattern[3][3] = 7;
											}
											//face 023
											if (current_element->GetIdNode(local_id_nodes[2]) < current_element->GetIdNode(local_id_nodes[3]))
											{
												pattern[4][0] = 3; pattern[4][1] = 5; pattern[4][2] = 6; pattern[4][3] = 7;
												pattern[5][0] = 2; pattern[5][1] = 3; pattern[5][2] = 5; pattern[5][3] = 7;
											}
											else {
												pattern[4][0] = 2; pattern[4][1] = 5; pattern[4][2] = 6; pattern[4][3] = 7;
												pattern[5][0] = 2; pattern[5][1] = 3; pattern[5][2] = 6; pattern[5][3] = 7;
											}

										}
										else {
											local_id_vertex[0][0] = node1; local_id_vertex[0][1] = node4;
											local_id_vertex[1][0] = node1; local_id_vertex[1][1] = node3;
											local_id_vertex[2][0] = node1; local_id_vertex[2][1] = node2;
											current_element->SubGrid_for_integration.xyz[4] = e14;
											current_element->SubGrid_for_integration.xyz[5] = e13;
											current_element->SubGrid_for_integration.xyz[6] = e12;

											//face 012
											if (current_element->GetIdNode(local_id_nodes[1]) < current_element->GetIdNode(local_id_nodes[2]))
											{
												pattern[0][0] = 2; pattern[0][1] = 6; pattern[0][2] = 5; pattern[0][3] = 7;
												pattern[1][0] = 1; pattern[1][1] = 2; pattern[1][2] = 6; pattern[1][3] = 7;
											}
											else {
												pattern[0][0] = 1; pattern[0][1] = 6; pattern[0][2] = 5; pattern[0][3] = 7;
												pattern[1][0] = 1; pattern[1][1] = 2; pattern[1][2] = 5; pattern[1][3] = 7;
											}
											//face 013
											if (current_element->GetIdNode(local_id_nodes[1]) < current_element->GetIdNode(local_id_nodes[3]))
											{
												pattern[2][0] = 3; pattern[2][1] = 6; pattern[2][2] = 4; pattern[2][3] = 7;
												pattern[3][0] = 1; pattern[3][1] = 3; pattern[3][2] = 6; pattern[3][3] = 7;
											}
											else {
												pattern[2][0] = 1; pattern[2][1] = 3; pattern[2][2] = 4; pattern[2][3] = 7;
												pattern[3][0] = 1; pattern[3][1] = 6; pattern[3][2] = 4; pattern[3][3] = 7;
											}
											//face 023
											if (current_element->GetIdNode(local_id_nodes[2]) < current_element->GetIdNode(local_id_nodes[3]))
											{
												pattern[4][0] = 3; pattern[4][1] = 5; pattern[4][2] = 4; pattern[4][3] = 7;
												pattern[5][0] = 2; pattern[5][1] = 3; pattern[5][2] = 5; pattern[5][3] = 7;
											}
											else {
												pattern[4][0] = 2; pattern[4][1] = 5; pattern[4][2] = 4; pattern[4][3] = 7;
												pattern[5][0] = 2; pattern[5][1] = 3; pattern[5][2] = 4; pattern[5][3] = 7;
											}
										}

										for (int p = 0; p < 9; p++)
										{
											std::vector<int> elem(4);
											for (int pp = 0; pp < 4; pp++)
											{
												elem[pp] = pattern[p][pp];
											}
											current_element->SubGrid_for_integration.nvtr.push_back(elem);
											current_element->SubGrid_for_integration.nvkat.push_back(0);
										}

										current_element->SubGrid_for_integr_byTriangles.xyz.resize(3);
										current_element->SubGrid_for_integr_byTriangles.xyz[0] = current_element->SubGrid_for_integration.xyz[4];
										current_element->SubGrid_for_integr_byTriangles.xyz[1] = current_element->SubGrid_for_integration.xyz[5];
										current_element->SubGrid_for_integr_byTriangles.xyz[2] = current_element->SubGrid_for_integration.xyz[6];
										std::vector<int> elem(3);
										elem[0] = 0; elem[1] = 1; elem[2] = 2;
										current_element->SubGrid_for_integr_byTriangles.nvtr.push_back(elem);
										current_element->SubGrid_for_integr_byTriangles.nvkat.push_back(id_crack);

										return true;
									}
									return false;
								};
								if (SUMM == 3 && SUMM_forNodes == 0)
								{
									int pattern_edge[4][3] = { { 0,1,2 },{ 0,3,4 },{ 1,3,5 },{ 2,4,5 } };
									int pattern_node[4][4] = { { 0,1,2,3 },{ 1,0,2,3 },{ 2,0,1,3 },{ 3,0,1,2 } };
									for (int mm = 0; mm < 4; mm++)
									{
										if (cross_section[pattern_edge[mm][0]] == 1 &&
											cross_section[pattern_edge[mm][1]] == 1 &&
											cross_section[pattern_edge[mm][2]] == 1)
										{
											std::vector<int> tr(3);
											tr[0] = num_of_tr[pattern_edge[mm][0]];
											tr[1] = num_of_tr[pattern_edge[mm][1]];
											tr[2] = num_of_tr[pattern_edge[mm][2]];
											_3edge(
												pattern_node[mm][0],
												pattern_node[mm][1],
												pattern_node[mm][2],
												pattern_node[mm][3],
												tr,
												edge[pattern_edge[mm][0]],
												edge[pattern_edge[mm][1]],
												edge[pattern_edge[mm][2]]);

											//current_element->set_id_domain(-3);
											//current_element->SetIdDomain(1);
											id_extended_disc_vertexes.push_back(current_element->GetIdNode(0));
											id_extended_disc_vertexes.push_back(current_element->GetIdNode(1));
											id_extended_disc_vertexes.push_back(current_element->GetIdNode(2));
											id_extended_disc_vertexes.push_back(current_element->GetIdNode(3));
										}
									}
								}
								//пересекает весь тетраэдр через 4 ребра
								auto _4edge = [&](int node1, int node2, int node3, int node4,
									std::vector<int> test_tr, Point<double> e13, Point<double> e14, Point<double> e23, Point<double> e24)
								{
									double A, B, C, D;
									int signum[4];
									Point<double> test_vector;
									bool reverse = false;
									for (int t = 0; t < test_tr.size(); t++)
									{
										double _A, _B, _C, _D;
										math::GetPlaneEquation(current_crack->triangles[test_tr[t]].GetNodes(), _A, _B, _C, _D);
										test_vector += Point<double>(_A, _B, _C) / (double)test_tr.size();
									}
									math::GetPlaneEquation(e24, e23, e13, A, B, C, D);
									if (Point<double>(A, B, C)*test_vector < 0)
									{
										math::GetPlaneEquation(e14, e13, e23, A, B, C, D);
										reverse = true;
									}

									for (int node = 0; node < current_element->GetNodesCount(); node++)
									{
										signum[node] = math::GetSignum(A*this->GetCoordinateViaID(current_element->GetIdNode(node)).x + B * this->GetCoordinateViaID(current_element->GetIdNode(node)).y + C * this->GetCoordinateViaID(current_element->GetIdNode(node)).z + D);
									}
									if (signum[node1] > 0 && signum[node2] > 0 && signum[node3] < 0 && signum[node4] < 0 ||
										signum[node1] < 0 && signum[node2] < 0 && signum[node3] > 0 && signum[node4] > 0)
									{
										int local_id_nodes[4] = { node1, node2, node3, node4 };
										int local_id_vertex[4][2] = { { node2, node4 },{ node2, node3 },{ node1, node3 },{ node1, node4 } };

										current_element->SubGrid_for_integration.xyz.resize(8);
										current_element->SubGrid_for_integration.xyz[0] = this->GetCoordinateViaID(current_element->GetIdNode(node1));
										current_element->SubGrid_for_integration.xyz[1] = this->GetCoordinateViaID(current_element->GetIdNode(node2));
										current_element->SubGrid_for_integration.xyz[2] = this->GetCoordinateViaID(current_element->GetIdNode(node3));
										current_element->SubGrid_for_integration.xyz[3] = this->GetCoordinateViaID(current_element->GetIdNode(node4));

										int pattern[6][4];
										if (!reverse)
										{
											current_element->SubGrid_for_integration.xyz[4] = e24;
											current_element->SubGrid_for_integration.xyz[5] = e23;
											current_element->SubGrid_for_integration.xyz[6] = e13;
											current_element->SubGrid_for_integration.xyz[7] = e14;
											//edge 01
											if (current_element->GetIdNode(local_id_nodes[0]) < current_element->GetIdNode(local_id_nodes[1]))
											{
												pattern[0][0] = 0; pattern[0][1] = 1; pattern[0][2] = 6; pattern[0][3] = 7;
												pattern[1][0] = 1; pattern[1][1] = 4; pattern[1][2] = 6; pattern[1][3] = 7;
												pattern[2][0] = 1; pattern[2][1] = 4; pattern[2][2] = 5; pattern[2][3] = 6;
											}
											else {
												pattern[0][0] = 0; pattern[0][1] = 1; pattern[0][2] = 4; pattern[0][3] = 5;
												pattern[1][0] = 0; pattern[1][1] = 4; pattern[1][2] = 5; pattern[1][3] = 6;
												pattern[2][0] = 0; pattern[2][1] = 4; pattern[2][2] = 6; pattern[2][3] = 7;
											}
											//edge 23
											if (current_element->GetIdNode(local_id_nodes[2]) < current_element->GetIdNode(local_id_nodes[3]))
											{
												pattern[3][0] = 2; pattern[3][1] = 3; pattern[3][2] = 5; pattern[3][3] = 6;
												pattern[4][0] = 3; pattern[4][1] = 4; pattern[4][2] = 5; pattern[4][3] = 6;
												pattern[5][0] = 3; pattern[5][1] = 4; pattern[5][2] = 6; pattern[5][3] = 7;
											}
											else {
												pattern[3][0] = 2; pattern[3][1] = 3; pattern[3][2] = 4; pattern[3][3] = 7;
												pattern[4][0] = 2; pattern[4][1] = 4; pattern[4][2] = 6; pattern[4][3] = 7;
												pattern[5][0] = 2; pattern[5][1] = 4; pattern[5][2] = 5; pattern[5][3] = 6;
											}
										}
										else {
											local_id_vertex[0][0] = node1; local_id_vertex[0][1] = node4;
											local_id_vertex[1][0] = node1; local_id_vertex[1][1] = node3;
											local_id_vertex[2][0] = node2; local_id_vertex[2][1] = node3;
											local_id_vertex[3][0] = node2; local_id_vertex[3][1] = node4;

											current_element->SubGrid_for_integration.xyz[4] = e14;
											current_element->SubGrid_for_integration.xyz[5] = e13;
											current_element->SubGrid_for_integration.xyz[6] = e23;
											current_element->SubGrid_for_integration.xyz[7] = e24;

											//edge 01
											if (current_element->GetIdNode(local_id_nodes[0]) < current_element->GetIdNode(local_id_nodes[1]))
											{
												pattern[0][0] = 0; pattern[0][1] = 1; pattern[0][2] = 5; pattern[0][3] = 4;
												pattern[1][0] = 1; pattern[1][1] = 7; pattern[1][2] = 5; pattern[1][3] = 4;
												pattern[2][0] = 1; pattern[2][1] = 7; pattern[2][2] = 6; pattern[2][3] = 5;
											}
											else {
												pattern[0][0] = 0; pattern[0][1] = 1; pattern[0][2] = 7; pattern[0][3] = 6;
												pattern[1][0] = 0; pattern[1][1] = 7; pattern[1][2] = 6; pattern[1][3] = 5;
												pattern[2][0] = 0; pattern[2][1] = 7; pattern[2][2] = 5; pattern[2][3] = 4;
											}
											//edge 23
											if (current_element->GetIdNode(local_id_nodes[2]) < current_element->GetIdNode(local_id_nodes[3]))
											{
												pattern[3][0] = 2; pattern[3][1] = 3; pattern[3][2] = 6; pattern[3][3] = 5;
												pattern[4][0] = 3; pattern[4][1] = 7; pattern[4][2] = 6; pattern[4][3] = 5;
												pattern[5][0] = 3; pattern[5][1] = 7; pattern[5][2] = 5; pattern[5][3] = 4;
											}
											else {
												pattern[3][0] = 2; pattern[3][1] = 3; pattern[3][2] = 7; pattern[3][3] = 4;
												pattern[4][0] = 2; pattern[4][1] = 7; pattern[4][2] = 5; pattern[4][3] = 4;
												pattern[5][0] = 2; pattern[5][1] = 7; pattern[5][2] = 6; pattern[5][3] = 5;
											}
										}

										for (int p = 0; p < 6; p++)
										{
											std::vector<int> elem(4);
											for (int pp = 0; pp < 4; pp++)
											{
												elem[pp] = pattern[p][pp];
											}
											current_element->SubGrid_for_integration.nvtr.push_back(elem);
											current_element->SubGrid_for_integration.nvkat.push_back(0);
										}

										current_element->SubGrid_for_integr_byTriangles.xyz.resize(4);
										current_element->SubGrid_for_integr_byTriangles.xyz[0] = current_element->SubGrid_for_integration.xyz[4];
										current_element->SubGrid_for_integr_byTriangles.xyz[1] = current_element->SubGrid_for_integration.xyz[5];
										current_element->SubGrid_for_integr_byTriangles.xyz[2] = current_element->SubGrid_for_integration.xyz[6];
										current_element->SubGrid_for_integr_byTriangles.xyz[3] = current_element->SubGrid_for_integration.xyz[7];
										std::vector<int> elem(3);
										elem[0] = 0; elem[1] = 1; elem[2] = 2;
										current_element->SubGrid_for_integr_byTriangles.nvtr.push_back(elem);
										current_element->SubGrid_for_integr_byTriangles.nvkat.push_back(id_crack);
										elem[0] = 0; elem[1] = 2; elem[2] = 3;
										current_element->SubGrid_for_integr_byTriangles.nvtr.push_back(elem);
										current_element->SubGrid_for_integr_byTriangles.nvkat.push_back(id_crack);

										return true;
									}
									
									return false;
								};
								if (SUMM == 4 && SUMM_forNodes == 0)//(_4edge(0, 1, 2, 3) || _4edge(0, 2, 1, 3) || _4edge(0, 3, 1, 2))
								{
									int pattern_edge[3][4] = { { 1,2,3,4 },{ 0,2,3,5 },{ 0,1,4,5 } };
									int pattern_node[3][4] = { { 0,1,2,3 },{ 0,2,1,3 },{ 0,3,1,2 } };
									for (int mm = 0; mm < 3; mm++)
									{
										if (cross_section[pattern_edge[mm][0]] == 1 &&
											cross_section[pattern_edge[mm][1]] == 1 &&
											cross_section[pattern_edge[mm][2]] == 1 &&
											cross_section[pattern_edge[mm][3]] == 1)
										{
											std::vector<int> tr(4);
											tr[0] = num_of_tr[pattern_edge[mm][0]];
											tr[1] = num_of_tr[pattern_edge[mm][1]];
											tr[2] = num_of_tr[pattern_edge[mm][2]];
											tr[3] = num_of_tr[pattern_edge[mm][3]];
											_4edge(
												pattern_node[mm][0],
												pattern_node[mm][1],
												pattern_node[mm][2],
												pattern_node[mm][3],
												tr,
												edge[pattern_edge[mm][0]],
												edge[pattern_edge[mm][1]],
												edge[pattern_edge[mm][2]],
												edge[pattern_edge[mm][3]]);

											//current_element->set_id_domain(-4);
											//current_element->SetIdDomain(1);
											id_extended_disc_vertexes.push_back(current_element->GetIdNode(0));
											id_extended_disc_vertexes.push_back(current_element->GetIdNode(1));
											id_extended_disc_vertexes.push_back(current_element->GetIdNode(2));
											id_extended_disc_vertexes.push_back(current_element->GetIdNode(3));
										}
									}
								}

								/*if (tranform && !add_elem)
								{
									printf("Error in updating %d elem\n", i);
									printf("P[%d]: %.5lf %.5lf %.5lf\n", current_element->get_id_DOF(0), current_element->GetNode(0).x, current_element->GetNode(0).y, current_element->GetNode(0).z);
									printf("P[%d]: %.5lf %.5lf %.5lf\n", current_element->get_id_DOF(1), current_element->GetNode(1).x, current_element->GetNode(1).y, current_element->GetNode(1).z);
									printf("P[%d]: %.5lf %.5lf %.5lf\n", current_element->get_id_DOF(2), current_element->GetNode(2).x, current_element->GetNode(2).y, current_element->GetNode(2).z);
									printf("P[%d]: %.5lf %.5lf %.5lf\n", current_element->get_id_DOF(3), current_element->GetNode(3).x, current_element->GetNode(3).y, current_element->GetNode(3).z);
								}*/

							}
						}

					}
				non_disc:
					if (id_extended_disc_vertexes.size() != 0)
					{
						math::MakeQuickSort(id_extended_disc_vertexes, 0, (int)id_extended_disc_vertexes.size() - 1);
						std::vector<int> _tmp;
						math::MakeRemovalOfDuplication(id_extended_disc_vertexes, _tmp);
						math::MakeCopyVector_A_into_B(_tmp, id_extended_disc_vertexes);
						for (int i = 0; i < id_extended_disc_vertexes.size(); i++)
						{
							int id_vertex = id_extended_disc_vertexes[i];
							this->AddDOF_ToEnd(id_vertex);
							/*std::vector<int> update_tetrahedrons;
							for (int id_edge = 0; id_edge < this->vertexes[id_vertex].GetUpperElementCount(); id_edge++)
							{
								for (int id_face = 0; id_face < this->vertexes[id_vertex].GetUpperElement(id_edge)->GetUpperElementCount(); id_face++)
								{
									for (int id_element = 0; id_element < this->vertexes[id_vertex].GetUpperElement(id_edge)
										->GetUpperElement(id_face)
										->GetUpperElementCount(); id_element++)
									{
										update_tetrahedrons.push_back(this->vertexes[id_vertex].GetUpperElement(id_edge)
											->GetUpperElement(id_face)
											->GetUpperElement(id_element)->GetSelfGlobalId());
									}
								}
							}
							math::MakeQuickSort(update_tetrahedrons);
							for (int id_element = 0; id_element < update_tetrahedrons.size(); id_element++)
							{
								if (id_element == 0 || update_tetrahedrons[id_element] != update_tetrahedrons[id_element - 1])
								{
									auto current_element = this->GetElement(update_tetrahedrons[id_element]);
									int id_base_DOF = current_element->GetLocalID_forDOF(id_vertex);
									if (id_base_DOF == -1) {
										printf_s("We have problem with element[%d]. We can't find the vertex %d.\n", update_tetrahedrons[id_element], id_vertex);
									}

									std::function<double(Point<double>)> HevisadeFunc = [current_crack, current_element](Point<double> X)->double
									{
										double H;

										int Hx = 0, Hxi = 0;

										int signum_crack;
										int summ = 0;

										Point<double> P[3];
										if (current_element->SubGrid_for_integration.xyz.size() != 0)
										{
											P[0] = current_element->SubGrid_for_integration.xyz[4];
											P[1] = current_element->SubGrid_for_integration.xyz[5];
											P[2] = current_element->SubGrid_for_integration.xyz[6];
										}
										else {
											P[0] = current_crack->triangles[0].GetNode(0);
											P[1] = current_crack->triangles[0].GetNode(1);
											P[2] = current_crack->triangles[0].GetNode(2);
										}
										double A, B, C, D;
										math::GetPlaneEquation(P, A, B, C, D);
										signum_crack = math::GetSignum(A*X.x + B * X.y + C * X.z + D);

										switch (signum_crack)
										{
										case -1: Hx = -1; break;
										case 0: Hx = 1; break;
										case 1: Hx = 1; break;
										}

										//H = (1 / 2.)*(Hx - Hxi);
										H = Hx;

										return H;
									};

									std::function< Point<double>(Point<double> X)> bf = [current_element, id_base_DOF, HevisadeFunc](Point<double> X) -> Point<double>
									{
										Point<double> result;
										Point<double> base_bf = (*current_element->GetBasisFunctionInLocalID(id_base_DOF))(X);
										double heviside = HevisadeFunc(X);

										result = base_bf * heviside;
										return result;
									};
									std::function< Point<Point<double>>(Point<double> X)> derivative_bf = [current_element, id_base_DOF, HevisadeFunc](Point<double> X) -> Point<Point<double>>
									{
										Point<Point<double>> result;
										Point<Point<double>> base_derivative_of_bf = (*current_element->GetDerivativeOfBasisFunctionInLocalID(id_base_DOF))(X);
										double heviside = HevisadeFunc(X);

										result.x = base_derivative_of_bf.x * heviside;
										result.y = base_derivative_of_bf.y * heviside;
										result.z = base_derivative_of_bf.z * heviside;
										return result;
									};

									current_element->AppEndDOF(base_dof_count + i, bf, derivative_bf);
								}
							}*/

							for (int id_element = 0; id_element < this->vertexes[id_vertex].GetUpperElementCount(); id_element++)
							{
								auto current_element = this->GetElement(this->vertexes[id_vertex].GetUpperElement(id_element)->GetSelfGlobalId());
								int id_base_DOF = current_element->GetLocalID_forDOF(id_vertex);
								if (id_base_DOF == -1) {
									printf_s("We have problem with element[%d]. We can't find the vertex %d.\n", this->vertexes[id_vertex].GetUpperElement(id_element)->GetSelfGlobalId(), id_vertex);
								}

								std::function<double(Point<double>)> HevisadeFunc = [current_crack, current_element](Point<double> X)->double
								{
									double H;

									int Hx = 0, Hxi = 0;

									int signum_crack;
									int summ = 0;

									Point<double> P[3];
									if (current_element->SubGrid_for_integration.xyz.size() != 0)
									{
										P[0] = current_element->SubGrid_for_integration.xyz[4];
										P[1] = current_element->SubGrid_for_integration.xyz[5];
										P[2] = current_element->SubGrid_for_integration.xyz[6];
									}
									else {
										P[0] = current_crack->triangles[0].GetNode(0);
										P[1] = current_crack->triangles[0].GetNode(1);
										P[2] = current_crack->triangles[0].GetNode(2);
									}
									double A, B, C, D;
									math::GetPlaneEquation(P, A, B, C, D);
									signum_crack = math::GetSignum(A*X.x + B * X.y + C * X.z + D);

									switch (signum_crack)
									{
									case -1: Hx = -1; break;
									case 0: Hx = 1; break;
									case 1: Hx = 1; break;
									}

									//H = (1 / 2.)*(Hx - Hxi);
									H = Hx;

									return H;
								};

								std::function< Point<double>(Point<double> X)> bf = [current_element, id_base_DOF, HevisadeFunc](Point<double> X) -> Point<double>
								{
									Point<double> result;
									Point<double> base_bf = (*current_element->GetBasisFunctionInLocalID(id_base_DOF))(X);
									double heviside = HevisadeFunc(X);

									result = base_bf * heviside;
									return result;
								};
								std::function< Point<Point<double>>(Point<double> X)> derivative_bf = [current_element, id_base_DOF, HevisadeFunc](Point<double> X) -> Point<Point<double>>
								{
									Point<Point<double>> result;
									Point<Point<double>> base_derivative_of_bf = (*current_element->GetDerivativeOfBasisFunctionInLocalID(id_base_DOF))(X);
									double heviside = HevisadeFunc(X);

									result.x = base_derivative_of_bf.x * heviside;
									result.y = base_derivative_of_bf.y * heviside;
									result.z = base_derivative_of_bf.z * heviside;
									return result;
								};

								current_element->AppEndDOF(this->GetDOFsCount() - 1, bf, derivative_bf);

							}
						}
					}
				}
				
				//сингул€рные
				for (int id_front = 0; id_front < current_crack->fronts.size(); id_front++)
				{
					std::vector <int> id_extended_singular_vertexes;
					//via elements
					for (int i = 0; i < this->GetElementsCount(); i++)
					{
						auto current_element = this->GetElement(i);

						Point<double> node[6];
						int pattern_for_node[4] = { 0,1,2,3 };
						int cross_section_forNode[4], SUMM_forNode = 0, num_of_front_forNode[4];
						cross_section_forNode[0] = 0;
						cross_section_forNode[1] = 0;
						cross_section_forNode[2] = 0;
						cross_section_forNode[3] = 0;
						for (int id_segment = 0; id_segment < current_crack->fronts[id_front].size() && SUMM_forNode < 2; id_segment++)
						{
							for (int n = 0; n < 4; n++)
							{
								if (cross_section_forNode[n] == 0)
								{
									double len_left = math::SolveLengthVector(current_element->GetNode(pattern_for_node[n]),
										current_crack->xyz[current_crack->fronts[id_front][id_segment].id_left]);
									double len_right = math::SolveLengthVector(current_element->GetNode(pattern_for_node[n]),
										current_crack->xyz[current_crack->fronts[id_front][id_segment].id_right]);
									double len_full = math::SolveLengthVector(current_crack->xyz[current_crack->fronts[id_front][id_segment].id_left],
										current_crack->xyz[current_crack->fronts[id_front][id_segment].id_right]);
									if (math::IsEqual(len_left + len_right, len_full))
									{
										cross_section_forNode[n] = 1;
										num_of_front_forNode[n] = id_segment;
										SUMM_forNode++;
									}
								}
							}
						}

						Point<double> edge[6];
						int pattern_for_edge[6][2] = { { 0,1 },{ 0,2 },{ 0,3 },{ 1,2 },{ 1,3 },{ 2,3 } };
						int cross_section_forEdge[6], SUMM_forEdge = 0, num_of_front_forEdge[6];
						cross_section_forEdge[0] = 0;
						cross_section_forEdge[1] = 0;
						cross_section_forEdge[2] = 0;
						cross_section_forEdge[3] = 0;
						cross_section_forEdge[4] = 0;
						cross_section_forEdge[5] = 0;
						for (int id_segment = 0; id_segment < current_crack->fronts[id_front].size() && SUMM_forNode < 2 && SUMM_forEdge < 2; id_segment++)
						{
							for (int e = 0; e < 6; e++)
							{
								if (cross_section_forEdge[e] == 0)
								{
									Point<double> cross = math::GetIntersectionPointOfLineAndLine(
										current_element->GetNode(pattern_for_edge[e][0]),
										current_element->GetNode(pattern_for_edge[e][1]),
										current_crack->xyz[current_crack->fronts[id_front][id_segment].id_left],
										current_crack->xyz[current_crack->fronts[id_front][id_segment].id_right]);
									if (math::IsLineContainThePoint(current_element->GetNode(pattern_for_edge[e][0]), current_element->GetNode(pattern_for_edge[e][1]), cross) &&
										math::IsLineContainThePoint(current_crack->xyz[current_crack->fronts[id_front][id_segment].id_left], current_crack->xyz[current_crack->fronts[id_front][id_segment].id_right], cross) &&
										cross != current_element->GetNode(pattern_for_edge[e][0]) &&
										cross != current_element->GetNode(pattern_for_edge[e][1]))
									{
										cross_section_forEdge[e] = 1;
										num_of_front_forEdge[e] = id_segment;
										SUMM_forEdge++;
										edge[e] = cross;
									}
								}
							}
						}

						Point<double> face[6];
						int pattern_for_face[4][3] = { { 0,1,2 },{ 0,1,3 },{ 0,2,3 },{ 1,2,3 } };
						int pattern_for_face_via_edge[4][3] = { { 0,1,3 },{ 0,2,4 },{ 1,2,5 },{ 3,4,5 } };
						int cross_section_forFace[4], SUMM_forFace = 0, num_of_front_forFace[4];
						cross_section_forFace[0] = 0;
						cross_section_forFace[1] = 0;
						cross_section_forFace[2] = 0;
						cross_section_forFace[3] = 0;
						for (int id_segment = 0; id_segment < current_crack->fronts[id_front].size() && SUMM_forNode < 2 && SUMM_forEdge < 2 && SUMM_forFace < 2; id_segment++)
						{
							for (int f = 0; f < 4; f++)
							{
								if (cross_section_forFace[f] == 0)
								{
									Point<double> cross = math::GetIntersectionPointOfLineAndPlane(
										current_element->GetNode(pattern_for_face[f][0]),
										current_element->GetNode(pattern_for_face[f][1]),
										current_element->GetNode(pattern_for_face[f][2]),
										current_crack->xyz[current_crack->fronts[id_front][id_segment].id_left],
										current_crack->xyz[current_crack->fronts[id_front][id_segment].id_right]);

									double S = math::SolveSquareTriangle(current_element->GetNode(pattern_for_face[f][0]),
										current_element->GetNode(pattern_for_face[f][1]),
										current_element->GetNode(pattern_for_face[f][2]));
									double S0 = math::SolveSquareTriangle(cross,
										current_element->GetNode(pattern_for_face[f][1]),
										current_element->GetNode(pattern_for_face[f][2]));
									double S1 = math::SolveSquareTriangle(current_element->GetNode(pattern_for_face[f][0]),
										cross,
										current_element->GetNode(pattern_for_face[f][2]));
									double S2 = math::SolveSquareTriangle(current_element->GetNode(pattern_for_face[f][0]),
										current_element->GetNode(pattern_for_face[f][1]),
										cross);
									if (math::IsEqual(S0 + S1 + S2, S) &&
										!math::IsEqual(S0, 0.0) && !math::IsEqual(S1, 0.0) && !math::IsEqual(S2, 0.0))
									{
										cross_section_forFace[f] = 1;
										num_of_front_forFace[f] = id_segment;
										SUMM_forFace++;
									}
								}
							}
						}

						std::vector<int> _tmp_id_segment;
						if (SUMM_forFace != 0 || SUMM_forEdge != 0 || SUMM_forNode != 0)
						{
							for (int n = 0; n < 4; n++)
							{
								if (cross_section_forNode[n] != 0)
								{
									id_extended_singular_vertexes.push_back(current_element->GetIdNode(pattern_for_node[n]));
									_tmp_id_segment.push_back(id_front);
								}
							}
							for (int e = 0; e < 6; e++)
							{
								if (cross_section_forEdge[e] != 0)
								{
									id_extended_singular_vertexes.push_back(current_element->GetIdNode(pattern_for_edge[e][0]));
									id_extended_singular_vertexes.push_back(current_element->GetIdNode(pattern_for_edge[e][1]));
									_tmp_id_segment.push_back(id_front);
								}
							}
							for (int f = 0; f < 4; f++)
							{
								if (cross_section_forFace[f] != 0)
								{
									id_extended_singular_vertexes.push_back(current_element->GetIdNode(pattern_for_face[f][0]));
									id_extended_singular_vertexes.push_back(current_element->GetIdNode(pattern_for_face[f][1]));
									id_extended_singular_vertexes.push_back(current_element->GetIdNode(pattern_for_face[f][2]));
									_tmp_id_segment.push_back(id_front);
								}
							}

						}
						if (SUMM_forFace == 0 && SUMM_forEdge == 0 && SUMM_forNode == 0)
						{
							for (int id_segment = 0; id_segment < current_crack->fronts[id_front].size(); id_segment++)
							{
								if (current_element->IsContainThePoint(current_crack->xyz[current_crack->fronts[id_front][id_segment].id_left]) &&
									current_element->IsContainThePoint(current_crack->xyz[current_crack->fronts[id_front][id_segment].id_right]))
								{
									_tmp_id_segment.push_back(id_front);
								}
							}
						}

						if (_tmp_id_segment.size() != 0)
						{
							std::vector<int> aa;
							math::MakeQuickSort(_tmp_id_segment);
							math::MakeRemovalOfDuplication(_tmp_id_segment, aa);
						}
					}

					//via nodes
					for (int id_node = 0; id_node < this->GetVertexCount(); id_node++)
					{
						std::function< double(Point<double>, double&, Point<double>&, Point<double>&)> Qr = [current_crack, id_front](Point<double> X, double &r, Point<double> &X_on_line, Point<double> &normal_Line) -> double
						{
							double MinLen = 1E+15;
							Point<double>_X;
							Point<double>_X_dF_dt;
							int NumFr = -1, NumSegm = -1;
							double Q;

							for (int p = 0; p < current_crack->id_points_in_fronts[id_front].size() - 3; p++)
							{
								Point<double>P0 = current_crack->xyz[current_crack->id_points_in_fronts[id_front][p]];
								Point<double>P1 = current_crack->xyz[current_crack->id_points_in_fronts[id_front][p + 1]];
								Point<double>P2 = current_crack->xyz[current_crack->id_points_in_fronts[id_front][p + 2]];
								Point<double>P3 = current_crack->xyz[current_crack->id_points_in_fronts[id_front][p + 3]];
								if (p == 0 && current_crack->id_points_in_fronts[id_front][p] == current_crack->id_points_in_fronts[id_front][p + 1])
								{
									P0 = P1 + (P1 - P2);
								}
								if ((p + 3) == (current_crack->id_points_in_fronts[id_front].size() - 1) &&
									current_crack->id_points_in_fronts[id_front][p + 2] == current_crack->id_points_in_fronts[id_front][p + 3])
								{
									P3 = P2 + (P2 - P1);
								}

								Point<double>r3 = P1 * 3 - P2 * 3 + P3 + P0;
								Point<double>r2 = P1 * (-5) + P2 * 4 - P3;
								Point<double>r1 = P0 * (-1) + P2;
								Point<double>r0 = P1 * 2;
								r3 = P0 * (-1) + P1 * 3 - P2 * 3 + P3;
								r2 = P0 * 2 + P1 * (-5) + P2 * 4 - P3;
								r1 = P0 * (-1) + P2;
								r0 = P1 * 2;
								auto F = [&](double t) -> Point<double>
								{
									return (r3*pow(t, 3) + r2 * pow(t, 2) + r1 * t + r0) / 2.;
								};
								Point<double>f_0 = F(0);
								Point<double>f_0_5 = F(0.5);
								Point<double>f_1 = F(1);

								Point<double>_r2 = P1 * 9 - P2 * 9 + P3 * 3 + P0 * 3;
								Point<double>_r1 = P1 * (-10) + P2 * 8 - P3 * 2;
								Point<double>_r0 = P0 * (-1) + P2;
								_r2 = r3 * 3.;
								_r1 = r2 * 2.;
								_r0 = r1;
								auto dF_dt = [&](double t) -> Point<double>
								{
									return (_r2*pow(t, 2) + _r1 * t + _r0) / 2.;
								};
								Point<double>dF_dt_0 = dF_dt(0);
								Point<double>dF_dt_0_5 = dF_dt(0.5);
								Point<double>dF_dt_1 = dF_dt(1);

								double k5 = (r3 * _r2) * (-1 / 4.);
								double k4 = (r3 * _r1 + r2 * _r2) * (-1 / 4.);
								double k3 = (r3 * _r0 + r2 * _r1 + r1 * _r2) * (-1 / 4.);
								double k2 = (r2 * _r0 + r1 * _r1 + r0 * _r2) * (-1 / 4.) + (_r2 * X) * (1 / 2.);
								double k1 = (r1 * _r0 + r0 * _r1) * (-1 / 4.) + (_r1 * X) * (1 / 2.);
								double k0 = (r0 * _r0) * (-1 / 4.) + (_r0 * X) * (1 / 2.);
								auto Eq = [&](double t) -> double
								{
									double val = k5 * pow(t, 5) + k4 * pow(t, 4) + k3 * pow(t, 3) + k2 * pow(t, 2) + k1 * t + k0;
									return val;
								};
								std::vector<double> Res;
								math::SolveEquationViaBisectionMethod(0.0, 1.0, Eq, Res);

								for (int r = 0; r < Res.size(); r++)
								{
									Point<double>curr = F(Res[r]);
									double curr_len = math::SolveLengthVector(curr, X);
									if (curr_len < MinLen)
									{
										MinLen = curr_len;
										_X = curr;
										NumFr = id_front;
										NumSegm = p;
										_X_dF_dt = dF_dt(Res[r]);
									}
								}
							}


							if (NumFr != -1)
							{
								Point<double>Line[2];
								Line[0] = _X;
								Line[1] = _X + _X_dF_dt;
								Point<double>test;
								for (int tt = 0; tt < 3; tt++)
								{
									bool find = false;
									if (current_crack->fronts[NumFr][NumSegm].id_left == current_crack->triangles[current_crack->fronts[NumFr][NumSegm].id_base_triangles].GetIdNode(tt)
										|| current_crack->fronts[NumFr][NumSegm].id_right == current_crack->triangles[current_crack->fronts[NumFr][NumSegm].id_base_triangles].GetIdNode(tt))
									{
										find = true;
									}
									if (find == false)
									{
										int id_node = current_crack->triangles[current_crack->fronts[NumFr][NumSegm].id_base_triangles].GetIdNode(tt);
										test = current_crack->xyz[id_node];
										break;
									}
								}

								X_on_line = _X;
								DenseMatrix<double, double> test_matrix;
								normal_Line = math::MakeNormalForLine(Line, test, test_matrix);
								Q = math::SolveAngleBetweenV1AndV2(normal_Line, X - X_on_line);
								r = MinLen;
							}
							else {
								Q = 0;
								r = 0;

								Point<double>Line[2];
								//start&end segments
								int p[2] = { 0, (int)current_crack->id_points_in_fronts[id_front].size() - 4 };
								for (int pp = 0; pp < 2; pp++)
								{
									Point<double>P0 = current_crack->xyz[current_crack->id_points_in_fronts[id_front][p[pp]]];
									Point<double>P1 = current_crack->xyz[current_crack->id_points_in_fronts[id_front][p[pp] + 1]];
									Point<double>P2 = current_crack->xyz[current_crack->id_points_in_fronts[id_front][p[pp] + 2]];
									Point<double>P3 = current_crack->xyz[current_crack->id_points_in_fronts[id_front][p[pp] + 3]];
									if (p[pp] == 0 && current_crack->id_points_in_fronts[id_front][p[pp]] == current_crack->id_points_in_fronts[id_front][p[pp] + 1])
									{
										P0 = P1 + (P1 - P2);
									}
									if ((p[pp] + 3) == (current_crack->id_points_in_fronts[id_front].size() - 1) &&
										current_crack->id_points_in_fronts[id_front][p[pp] + 2] == current_crack->id_points_in_fronts[id_front][p[pp] + 3])
									{
										P3 = P2 + (P2 - P1);
									}

									Line[0] = P0;
									Line[1] = P1;
									Point<double>curr = math::MakeProgectionOfPointIntoLine(Line, X);
									double curr_len = math::SolveLengthVector(curr, X);
									if (curr_len < MinLen)
									{
										MinLen = curr_len;
										_X = curr;
										NumFr = id_front;
										NumSegm = p[pp];
										//_X_dF_dt = 0;
									}
								}

								Point<double>test;
								for (int tt = 0; tt < 3; tt++)
								{
									bool find = false;
									try {
										if (current_crack->fronts[NumFr][NumSegm].id_left == current_crack->triangles[current_crack->fronts[NumFr][NumSegm].id_base_triangles].GetIdNode(tt)
											|| current_crack->fronts[NumFr][NumSegm].id_right == current_crack->triangles[current_crack->fronts[NumFr][NumSegm].id_base_triangles].GetIdNode(tt))
										{
											find = true;
										}
									}
									catch (...)
									{
										printf_s("Number of front is %d\n", NumFr);
										printf_s("Number of segment is %d\n", NumSegm);

									};
									if (find == false)
									{
										test = current_crack->xyz[current_crack->triangles[current_crack->fronts[NumFr][NumSegm].id_base_triangles].GetIdNode(tt)];
										break;
									}
								}

								X_on_line = _X;
								DenseMatrix<double, double> test_matrix;
								normal_Line = math::MakeNormalForLine(Line, test, test_matrix);
								Q = math::SolveAngleBetweenV1AndV2(normal_Line, X - X_on_line);
								r = MinLen;

							};

							return Q;
						};
						Point<double> X_on_line, Normal;
						double r = 0;
						double Q = Qr(this->GetCoordinateViaID(id_node), r, X_on_line, Normal);

						if (r < current_crack->R) {
							id_extended_singular_vertexes.push_back(id_node);
						}
					}
				
					if (id_extended_singular_vertexes.size() != 0)
					{
						math::MakeQuickSort(id_extended_singular_vertexes);
						std::vector<int> _tmp;
						math::MakeRemovalOfDuplication(id_extended_singular_vertexes, _tmp);
						math::MakeCopyVector_A_into_B(_tmp, id_extended_singular_vertexes);

						for (int i = 0; i < id_extended_singular_vertexes.size(); i++)
						{
							int id_vertex = id_extended_singular_vertexes[i];
							this->AddDOF_ToEnd(id_vertex);
							this->AddDOF_ToEnd(id_vertex);
							this->AddDOF_ToEnd(id_vertex);
							this->AddDOF_ToEnd(id_vertex);

							//std::vector<int> update_tetrahedrons;
							//for (int id_edge = 0; id_edge < this->vertexes[id_vertex].GetUpperElementCount(); id_edge++)
							//{
							//	for (int id_face = 0; id_face < this->vertexes[id_vertex].GetUpperElement(id_edge)->GetUpperElementCount(); id_face++)
							//	{
							//		for (int id_element = 0; id_element < this->vertexes[id_vertex].GetUpperElement(id_edge)
							//			->GetUpperElement(id_face)
							//			->GetUpperElementCount(); id_element++)
							//		{
							//			update_tetrahedrons.push_back(this->vertexes[id_vertex].GetUpperElement(id_edge)
							//				->GetUpperElement(id_face)
							//				->GetUpperElement(id_element)->GetSelfGlobalId());
							//		}
							//	}
							//}
							//math::MakeQuickSort(update_tetrahedrons);
							//for (int id_element = 0; id_element < update_tetrahedrons.size(); id_element++)
							//{
							//	if (id_element == 0 || update_tetrahedrons[id_element] != update_tetrahedrons[id_element - 1])
							//	{
							//		auto current_element = this->GetElement(update_tetrahedrons[id_element]);
							for (int id_element = 0; id_element < this->vertexes[id_vertex].GetUpperElementCount(); id_element++)
							{
								auto current_element = this->GetElement(this->vertexes[id_vertex].GetUpperElement(id_element)->GetSelfGlobalId());
								current_element->SetIdDomain(current_element->GetIdDomain() + 1);

								if (current_element->GetIdDomain() > 4)
								{
									printf("Error in elem[%d]", current_element->GetSelfGlobalId());
								}

								int id_base_DOF = current_element->GetLocalID_forDOF(id_vertex);
								if (id_base_DOF == -1) {
									printf_s("We have problem with element[%d]. We can't find the vertex %d.\n", id_element, id_vertex);
								}

								std::function< double(Point<double>, double&, Point<double>&, Point<double>&)> Qr = [current_crack, id_front](Point<double> X, double &r, Point<double> &X_on_line, Point<double> &normal_Line) -> double
								{
									double MinLen = 1E+15;
									Point<double>_X;
									Point<double>_X_dF_dt;
									int NumFr = -1, NumSegm = -1;
									double Q;

									for (int p = 0; p < current_crack->id_points_in_fronts[id_front].size() - 3; p++)
									{
										Point<double>P0 = current_crack->xyz[current_crack->id_points_in_fronts[id_front][p]];
										Point<double>P1 = current_crack->xyz[current_crack->id_points_in_fronts[id_front][p + 1]];
										Point<double>P2 = current_crack->xyz[current_crack->id_points_in_fronts[id_front][p + 2]];
										Point<double>P3 = current_crack->xyz[current_crack->id_points_in_fronts[id_front][p + 3]];
										if (p == 0 && current_crack->id_points_in_fronts[id_front][p] == current_crack->id_points_in_fronts[id_front][p + 1])
										{
											P0 = P1 + (P1 - P2);
										}
										if ((p + 3) == (current_crack->id_points_in_fronts[id_front].size() - 1) &&
											current_crack->id_points_in_fronts[id_front][p + 2] == current_crack->id_points_in_fronts[id_front][p + 3])
										{
											P3 = P2 + (P2 - P1);
										}

										Point<double>r3 = P1 * 3 - P2 * 3 + P3 + P0;
										Point<double>r2 = P1 * (-5) + P2 * 4 - P3;
										Point<double>r1 = P0 * (-1) + P2;
										Point<double>r0 = P1 * 2;
										r3 = P0 * (-1) + P1 * 3 - P2 * 3 + P3;
										r2 = P0 * 2 + P1 * (-5) + P2 * 4 - P3;
										r1 = P0 * (-1) + P2;
										r0 = P1 * 2;
										auto F = [&](double t) -> Point<double>
										{
											return (r3*pow(t, 3) + r2 * pow(t, 2) + r1 * t + r0) / 2.;
										};
										Point<double>f_0 = F(0);
										Point<double>f_0_5 = F(0.5);
										Point<double>f_1 = F(1);

										Point<double>_r2 = P1 * 9 - P2 * 9 + P3 * 3 + P0 * 3;
										Point<double>_r1 = P1 * (-10) + P2 * 8 - P3 * 2;
										Point<double>_r0 = P0 * (-1) + P2;
										_r2 = r3 * 3.;
										_r1 = r2 * 2.;
										_r0 = r1;
										auto dF_dt = [&](double t) -> Point<double>
										{
											return (_r2*pow(t, 2) + _r1 * t + _r0) / 2.;
										};
										Point<double>dF_dt_0 = dF_dt(0);
										Point<double>dF_dt_0_5 = dF_dt(0.5);
										Point<double>dF_dt_1 = dF_dt(1);

										double k5 = (r3 * _r2) * (-1 / 4.);
										double k4 = (r3 * _r1 + r2 * _r2) * (-1 / 4.);
										double k3 = (r3 * _r0 + r2 * _r1 + r1 * _r2) * (-1 / 4.);
										double k2 = (r2 * _r0 + r1 * _r1 + r0 * _r2) * (-1 / 4.) + (_r2 * X) * (1 / 2.);
										double k1 = (r1 * _r0 + r0 * _r1) * (-1 / 4.) + (_r1 * X) * (1 / 2.);
										double k0 = (r0 * _r0) * (-1 / 4.) + (_r0 * X) * (1 / 2.);
										auto Eq = [&](double t) -> double
										{
											double val = k5 * pow(t, 5) + k4 * pow(t, 4) + k3 * pow(t, 3) + k2 * pow(t, 2) + k1 * t + k0;
											return val;
										};
										std::vector<double> Res;
										math::SolveEquationViaBisectionMethod(0.0, 1.0, Eq, Res);

										for (int r = 0; r < Res.size(); r++)
										{
											Point<double>curr = F(Res[r]);
											double curr_len = math::SolveLengthVector(curr, X);
											if (curr_len < MinLen)
											{
												MinLen = curr_len;
												_X = curr;
												NumFr = id_front;
												NumSegm = p;
												_X_dF_dt = dF_dt(Res[r]);
											}
										}
									}


									if (NumFr != -1)
									{
										Point<double>Line[2];
										Line[0] = _X;
										Line[1] = _X + _X_dF_dt;
										Point<double>test;
										for (int tt = 0; tt < 3; tt++)
										{
											bool find = false;
											if (current_crack->fronts[NumFr][NumSegm].id_left == current_crack->triangles[current_crack->fronts[NumFr][NumSegm].id_base_triangles].GetIdNode(tt)
												|| current_crack->fronts[NumFr][NumSegm].id_right == current_crack->triangles[current_crack->fronts[NumFr][NumSegm].id_base_triangles].GetIdNode(tt))
											{
												find = true;
											}
											if (find == false)
											{
												int id_node = current_crack->triangles[current_crack->fronts[NumFr][NumSegm].id_base_triangles].GetIdNode(tt);
												test = current_crack->xyz[id_node];
												break;
											}
										}

										X_on_line = _X;
										DenseMatrix<double, double> test_matrix;
										normal_Line = math::MakeNormalForLine(Line, test, test_matrix);
										Q = math::SolveAngleBetweenV1AndV2(normal_Line, X - X_on_line);
										r = MinLen;
									}
									else {
										Q = 0;
										r = 0;

										Point<double>Line[2];
										//start&end segments
										int p[2] = { 0, (int)current_crack->id_points_in_fronts[id_front].size() - 4 };
										for (int pp = 0; pp < 2; pp++)
										{
											Point<double>P0 = current_crack->xyz[current_crack->id_points_in_fronts[id_front][p[pp]]];
											Point<double>P1 = current_crack->xyz[current_crack->id_points_in_fronts[id_front][p[pp] + 1]];
											Point<double>P2 = current_crack->xyz[current_crack->id_points_in_fronts[id_front][p[pp] + 2]];
											Point<double>P3 = current_crack->xyz[current_crack->id_points_in_fronts[id_front][p[pp] + 3]];
											if (p[pp] == 0 && current_crack->id_points_in_fronts[id_front][p[pp]] == current_crack->id_points_in_fronts[id_front][p[pp] + 1])
											{
												P0 = P1 + (P1 - P2);
											}
											if ((p[pp] + 3) == (current_crack->id_points_in_fronts[id_front].size() - 1) &&
												current_crack->id_points_in_fronts[id_front][p[pp] + 2] == current_crack->id_points_in_fronts[id_front][p[pp] + 3])
											{
												P3 = P2 + (P2 - P1);
											}

											Line[0] = P0;
											Line[1] = P1;
											Point<double>curr = math::MakeProgectionOfPointIntoLine(Line, X);
											double curr_len = math::SolveLengthVector(curr, X);
											if (curr_len < MinLen)
											{
												MinLen = curr_len;
												_X = curr;
												NumFr = id_front;
												NumSegm = p[pp];
												//_X_dF_dt = 0;
											}
										}

										Point<double>test;
										for (int tt = 0; tt < 3; tt++)
										{
											bool find = false;
											try {
												if (current_crack->fronts[NumFr][NumSegm].id_left == current_crack->triangles[current_crack->fronts[NumFr][NumSegm].id_base_triangles].GetIdNode(tt)
													|| current_crack->fronts[NumFr][NumSegm].id_right == current_crack->triangles[current_crack->fronts[NumFr][NumSegm].id_base_triangles].GetIdNode(tt))
												{
													find = true;
												}
											}
											catch (...)
											{
												printf_s("Number of front is %d\n", NumFr);
												printf_s("Number of segment is %d\n", NumSegm);

											};
											if (find == false)
											{
												test = current_crack->xyz[current_crack->triangles[current_crack->fronts[NumFr][NumSegm].id_base_triangles].GetIdNode(tt)];
												break;
											}
										}

										X_on_line = _X;
										DenseMatrix<double, double> test_matrix;
										normal_Line = math::MakeNormalForLine(Line, test, test_matrix);
										Q = math::SolveAngleBetweenV1AndV2(normal_Line, X - X_on_line);
										r = MinLen;

									};

									return Q;
								}
								;
								for (int id_function = 0; id_function < 4; id_function++)
								{
									std::function< Point<double>(Point<double> X)> bf = [current_element, id_base_DOF, Qr, id_function](Point<double> X) -> Point<double>
									{
										Point<double> result;
										Point<double> base_bf = (*current_element->GetBasisFunctionInLocalID(id_base_DOF))(X);

										double val = 0;
										Point<double> X_on_line, normal_Line;
										double Q, r;
										Q = Qr(X, r, X_on_line, normal_Line);

										switch (id_function)
										{
										case 0:
											val = sqrt(r)*sin(Q / 2);
											break;
										case 1:
											val = sqrt(r)*cos(Q / 2);
											break;
										case 2:
											val = sqrt(r)*sin(Q / 2)*sin(Q);
											break;
										case 3:
											val = sqrt(r)*cos(Q / 2)*sin(Q);
											break;
										default:
											break;
										}

										result = base_bf * val;
										return result;
									};
									std::function< Point<Point<double>>(Point<double> X)> derivative_bf = [current_element, id_base_DOF, Qr, id_function](Point<double> X) -> Point<Point<double>>
									{
										Point<Point<double>> result;
										Point<double> base_bf = (*current_element->GetBasisFunctionInLocalID(id_base_DOF))(X);
										Point<Point<double>> base_derivative_of_bf = (*current_element->GetDerivativeOfBasisFunctionInLocalID(id_base_DOF))(X);

										double ksi = 0;
										Point<double> dksi;

										double val = 0;
										Point<double> X_on_line, normal_Line;
										double Q, r;
										Q = Qr(X, r, X_on_line, normal_Line);

										Point<double> _XX = X - X_on_line;
										Point<double> dr;
										dr.x = (_XX.x * 2) / (2. * r);
										dr.y = (_XX.y * 2) / (2. * r);
										dr.z = (_XX.z * 2) / (2. * r);
										double forQ = _XX * normal_Line / r;
										Point<double> dQ;
										dQ.x = (-1. / sqrt(1 - forQ * forQ)) * ((normal_Line.x)*r - _XX * normal_Line*dr.x) / (r*r);
										dQ.y = (-1. / sqrt(1 - forQ * forQ)) * ((normal_Line.y)*r - _XX * normal_Line*dr.y) / (r*r);
										dQ.z = (-1. / sqrt(1 - forQ * forQ)) * ((normal_Line.z)*r - _XX * normal_Line*dr.z) / (r*r);
										double du_dr, du_dQ;

										switch (id_function)
										{
										case 0:
											du_dr = sin(Q / 2.) / (2.*sqrt(r));
											du_dQ = sqrt(r)*cos(Q / 2.) / 2.;
											ksi = sqrt(r)*sin(Q / 2.);
											dksi.x = du_dr * dr.x + du_dQ * dQ.x;
											dksi.y = du_dr * dr.y + du_dQ * dQ.y;
											dksi.z = du_dr * dr.z + du_dQ * dQ.z;
											break;
										case 1:
											du_dr = cos(Q / 2.) / (2.*sqrt(r));
											du_dQ = -sqrt(r)*sin(Q / 2.) / 2.;
											ksi = sqrt(r)*cos(Q / 2.);
											dksi.x = du_dr * dr.x + du_dQ * dQ.x;
											dksi.y = du_dr * dr.y + du_dQ * dQ.y;
											dksi.z = du_dr * dr.z + du_dQ * dQ.z;
											break;
										case 2:
											du_dr = sin(Q / 2.)*sin(Q) / (2.*sqrt(r));
											du_dQ = sqrt(r)*(cos(Q / 2.)*sin(Q) / 2. + sin(Q / 2.)*cos(Q));
											ksi = sqrt(r)*sin(Q / 2.)*sin(Q);
											dksi.x = du_dr * dr.x + du_dQ * dQ.x;
											dksi.y = du_dr * dr.y + du_dQ * dQ.y;
											dksi.z = du_dr * dr.z + du_dQ * dQ.z;
											break;
										case 3:
											du_dr = cos(Q / 2.)*sin(Q) / (2.*sqrt(r));
											du_dQ = sqrt(r)*(-sin(Q / 2.)*sin(Q) / 2. + cos(Q / 2.)*cos(Q));
											ksi = sqrt(r)*cos(Q / 2.)*sin(Q);
											dksi.x = du_dr * dr.x + du_dQ * dQ.x;
											dksi.y = du_dr * dr.y + du_dQ * dQ.y;
											dksi.z = du_dr * dr.z + du_dQ * dQ.z;
											break;
										default:
											break;
										}

										result.x = base_derivative_of_bf.x * ksi + base_bf * dksi.x;
										result.y = base_derivative_of_bf.y * ksi + base_bf * dksi.y;
										result.z = base_derivative_of_bf.z * ksi + base_bf * dksi.z;

										return result;
									};
									std::function< Point<Point<double>>(Point<double> X)> derivative_bf_numerical = [current_element, id_base_DOF, Qr, id_function](Point<double> X) -> Point<Point<double>>
									{
										Point<Point<double>> result;
										Point<double> base_bf = (*current_element->GetBasisFunctionInLocalID(id_base_DOF))(X);
										Point<Point<double>> base_derivative_of_bf = (*current_element->GetDerivativeOfBasisFunctionInLocalID(id_base_DOF))(X);

										double ksi = 0;
										Point<double> dksi, dksi_new, _ksi;

										double val = 0;
										Point<double> X_on_line, normal_Line;
										double Q, r;
										Q = Qr(X, r, X_on_line, normal_Line);

										//for numerical derivative
										double delta = 0.01*abs(current_element->GetNode(0).z - current_element->GetNode(1).z);
										Point<double> _Q, _r;
										_Q.x = Qr(Point<double>(X.x + delta, X.y, X.z), _r.x, X_on_line, normal_Line);
										_Q.y = Qr(Point<double>(X.x, X.y + delta, X.z), _r.y, X_on_line, normal_Line);
										_Q.z = Qr(Point<double>(X.x, X.y, X.z + delta), _r.z, X_on_line, normal_Line);



										Point<double> _XX = X - X_on_line;
										Point<double> dr;
										dr.x = (_XX.x * 2) / (2. * r);
										dr.y = (_XX.y * 2) / (2. * r);
										dr.z = (_XX.z * 2) / (2. * r);
										double forQ = _XX * normal_Line / r;
										Point<double> dQ;
										dQ.x = (-1. / sqrt(1 - forQ * forQ)) * ((normal_Line.x)*r - _XX * normal_Line*dr.x) / (r*r);
										dQ.y = (-1. / sqrt(1 - forQ * forQ)) * ((normal_Line.y)*r - _XX * normal_Line*dr.y) / (r*r);
										dQ.z = (-1. / sqrt(1 - forQ * forQ)) * ((normal_Line.z)*r - _XX * normal_Line*dr.z) / (r*r);
										double du_dr, du_dQ;

										switch (id_function)
										{
										case 0:
											du_dr = sin(Q / 2.) / (2.*sqrt(r));
											du_dQ = sqrt(r)*cos(Q / 2.) / 2.;
											ksi = sqrt(r)*sin(Q / 2.);
											dksi.x = du_dr * dr.x + du_dQ * dQ.x;
											dksi.y = du_dr * dr.y + du_dQ * dQ.y;
											dksi.z = du_dr * dr.z + du_dQ * dQ.z;

											_ksi.x = sqrt(_r.x)*sin(_Q.x / 2.);
											_ksi.y = sqrt(_r.y)*sin(_Q.y / 2.);
											_ksi.z = sqrt(_r.z)*sin(_Q.z / 2.);
											dksi_new.x = (_ksi.x - ksi) / delta;
											dksi_new.y = (_ksi.y - ksi) / delta;
											dksi_new.z = (_ksi.z - ksi) / delta;
											break;
										case 1:
											du_dr = cos(Q / 2.) / (2.*sqrt(r));
											du_dQ = -sqrt(r)*sin(Q / 2.) / 2.;
											ksi = sqrt(r)*cos(Q / 2.);
											dksi.x = du_dr * dr.x + du_dQ * dQ.x;
											dksi.y = du_dr * dr.y + du_dQ * dQ.y;
											dksi.z = du_dr * dr.z + du_dQ * dQ.z;

											_ksi.x = sqrt(_r.x)*cos(_Q.x / 2.);
											_ksi.y = sqrt(_r.y)*cos(_Q.y / 2.);
											_ksi.z = sqrt(_r.z)*cos(_Q.z / 2.);
											dksi_new.x = (_ksi.x - ksi) / delta;
											dksi_new.y = (_ksi.y - ksi) / delta;
											dksi_new.z = (_ksi.z - ksi) / delta;
											break;
										case 2:
											du_dr = sin(Q / 2.)*sin(Q) / (2.*sqrt(r));
											du_dQ = sqrt(r)*(cos(Q / 2.)*sin(Q) / 2. + sin(Q / 2.)*cos(Q));
											ksi = sqrt(r)*sin(Q / 2.)*sin(Q);
											dksi.x = du_dr * dr.x + du_dQ * dQ.x;
											dksi.y = du_dr * dr.y + du_dQ * dQ.y;
											dksi.z = du_dr * dr.z + du_dQ * dQ.z;

											_ksi.x = sqrt(_r.x)*sin(_Q.x / 2.)*sin(_Q.x);
											_ksi.y = sqrt(_r.y)*sin(_Q.y / 2.)*sin(_Q.y);
											_ksi.z = sqrt(_r.z)*sin(_Q.z / 2.)*sin(_Q.z);
											dksi_new.x = (_ksi.x - ksi) / delta;
											dksi_new.y = (_ksi.y - ksi) / delta;
											dksi_new.z = (_ksi.z - ksi) / delta;
											break;
										case 3:
											du_dr = cos(Q / 2.)*sin(Q) / (2.*sqrt(r));
											du_dQ = sqrt(r)*(-sin(Q / 2.)*sin(Q) / 2. + cos(Q / 2.)*cos(Q));
											ksi = sqrt(r)*cos(Q / 2.)*sin(Q);
											dksi.x = du_dr * dr.x + du_dQ * dQ.x;
											dksi.y = du_dr * dr.y + du_dQ * dQ.y;
											dksi.z = du_dr * dr.z + du_dQ * dQ.z;

											_ksi.x = sqrt(_r.x)*cos(_Q.x / 2.)*sin(_Q.x);
											_ksi.y = sqrt(_r.y)*cos(_Q.y / 2.)*sin(_Q.y);
											_ksi.z = sqrt(_r.z)*cos(_Q.z / 2.)*sin(_Q.z);
											dksi_new.x = (_ksi.x - ksi) / delta;
											dksi_new.y = (_ksi.y - ksi) / delta;
											dksi_new.z = (_ksi.z - ksi) / delta;

											break;
										default:
											break;
										}

										result.x = base_derivative_of_bf.x * ksi + base_bf * dksi.x;
										result.y = base_derivative_of_bf.y * ksi + base_bf * dksi.y;
										result.z = base_derivative_of_bf.z * ksi + base_bf * dksi.z;

										result.x = base_derivative_of_bf.x * ksi + base_bf * dksi_new.x;
										result.y = base_derivative_of_bf.y * ksi + base_bf * dksi_new.y;
										result.z = base_derivative_of_bf.z * ksi + base_bf * dksi_new.z;

										return result;
									};


									current_element->AppEndDOF(this->GetDOFsCount() - 4 + id_function, bf, derivative_bf);
									
									
									//current_element->AppEndDOF(this->GetDOFsCount() - 4 + id_function, bf, derivative_bf_numerical);
									//current_element->SetIdDomain(2);
									//printf_s("\n%d\n", current_element->GetIdDomain());
									//current_element->AppEndDOF(base_dof_count + i*4 + id_function, bf, derivative_bf);
								}

							}
						}
					}
				}
			}

		}
		void CreateFunctionalSpaces()
		{
			this->SetDOFsCount((int)this->vertexes.size());
			for (int id_element = 0; id_element < this->GetElementsCount(); id_element++)
			{
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
			}

			//update boundary conditions
			for (int id_vertex = 0; id_vertex < this->boundary_vertexes.size(); id_vertex++)
			{
				auto boundary = &(this->boundary_vertexes[id_vertex]);
				std::function< Point<double>(Point<double> X)> bf = [boundary](Point<double> X) -> Point<double>
				{
					return Point<double> (0,0,0);
				};
				std::function< Point<Point<double>>(Point<double> X)> derivative_bf = [boundary](Point<double> X) -> Point<Point<double>>
				{
					Point<Point<double>> result;
					return result;
				};

				this->boundary_vertexes[id_vertex].ResizeDOF(1);
				this->boundary_vertexes[id_vertex].SetDOF(0, this->boundary_vertexes[id_vertex].GetIdNode(0), bf, derivative_bf);
			}
			for (int id_face = 0; id_face < this->boundary_faces.size(); id_face++)
			{
				auto boundary = &(this->boundary_faces[id_face]);
				
				int id_base_tetr = boundary->GetUpperElement(0)->GetSelfGlobalId();
				auto base_tetr = this->GetElement(id_base_tetr);
				boundary->ResizeDOF(boundary->GetNodesCount());
				for (int id_vertex = 0; id_vertex < boundary->GetNodesCount(); id_vertex++)
				{
					//int id_base_tetr = this->vertexes[boundary->GetIdNode(id_vertex)].GetUpperElement(0)->GetSelfGlobalId();
					//auto base_tetr = this->GetElement(id_base_tetr);
					int local_id_in_tetr = base_tetr->GetLocalID_forVertexID(boundary->GetIdNode(id_vertex));
					std::function< Point<double>(Point<double> X)> bf = *(base_tetr->GetBasisFunctionInLocalID(local_id_in_tetr));
					std::function< Point<Point<double>>(Point<double> X)> derivative_bf = *(base_tetr->GetDerivativeOfBasisFunctionInLocalID(local_id_in_tetr));
					boundary->SetDOF(id_vertex, boundary->GetIdNode(id_vertex), bf, derivative_bf);
				}
			}
		}
		void Initialization(math::SimpleGrid &base_grid, std::vector<std::vector<int>> first_boundary, std::vector<std::vector<int>> second_boundary)
		{
			try {
				math::MakeCopyVector_A_into_B(base_grid.xyz, *(this->GetCoordinates()));
				this->CreateXYZline();

				this->SetElementsCount((int)base_grid.nvtr.size());

				for (int id_element = 0; id_element < this->GetElementsCount(); id_element++)
				{
					auto new_element = this->GetElement(id_element);

					new_element->SetIdNodes(base_grid.nvtr[id_element]);
					for (int i = 0; i < new_element->GetNodesCount(); i++)
					{
						new_element->SetNode(i, this->GetPtrCoordinateViaID(new_element->GetIdNode(i)));
					}

					new_element->SetIdDomain(base_grid.nvkat[id_element]);
				}

				this->CreateQTree();


				//create topology
				CreateTopology();

				//1st boundary
				this->boundary_vertexes.resize(first_boundary.size());
				for (int id_vertex = 0; id_vertex < first_boundary.size(); id_vertex++)
				{
					this->boundary_vertexes[id_vertex].SetTypeBoundary(first_boundary[id_vertex][0]);
					this->boundary_vertexes[id_vertex].SetIdNode(0, first_boundary[id_vertex][1]);
					this->boundary_vertexes[id_vertex].SetNode(0, this->GetPtrCoordinateViaID(first_boundary[id_vertex][1]));
				}

				//2nd boundary

				CreateFunctionalSpaces();
				//CreateExtendedFunctionalSpace_ForCrack();
			}
			catch (const std::exception&)
			{
				printf_s("Error: MultiXFEM/MultiXFEM_Grid.h/void Initialization(geometry::Grid<geometry::Tetrahedron> &base_grid)\n");
			}
		}
		template <typename Dirichlet, typename Neumann>
		void Initialization(math::SimpleGrid &base_grid, std::vector<Dirichlet> &first_boundary, std::vector<Neumann> &second_boundary)
		{
			try {
				math::MakeCopyVector_A_into_B(base_grid.xyz, *(this->GetCoordinates()));
				this->CreateXYZline();

				this->SetElementsCount((int)base_grid.nvtr.size());

				for (int id_element = 0; id_element < this->GetElementsCount(); id_element++)
				{
					auto new_element = this->GetElement(id_element);

					new_element->SetIdNodes(base_grid.nvtr[id_element]);
					for (int i = 0; i < new_element->GetNodesCount(); i++)
					{
						new_element->SetNode(i, this->GetPtrCoordinateViaID(new_element->GetIdNode(i)));
					}

					new_element->SetIdDomain(base_grid.nvkat[id_element]);
				}

				this->CreateQTree();


				//create topology
				CreateTopology();

				//1st boundary
				{
					int N_1st = 0;
					for (int id_type = 0; id_type < first_boundary.size(); id_type++)
					{
						N_1st += first_boundary[id_type].id_vertexes.size();
					}
					this->boundary_vertexes.resize(N_1st);
					int id = 0;
					for (int id_type = 0; id_type < first_boundary.size(); id_type++)
					{
						for (int id_vertex = 0; id_vertex < first_boundary[id_type].id_vertexes.size(); id_vertex++)
						{
							this->boundary_vertexes[id].boundary_value = first_boundary[id_type].value;
							//Point<bool> _t;
							//auto val_old = first_boundary[id_type].value(_t);
							//auto val = this->boundary_vertexes[id].boundary_value(_t);
							//printf_s("b_v[%d]=(%.2lf, %.2lf, %.2lf) ->>> b_v_old[%d]=(%.2lf, %.2lf, %.2lf)\n", id, val.x, val.y, val.z, id, val_old.x, val_old.y, val_old.z);
							this->boundary_vertexes[id].SetIdNode(0, first_boundary[id_type].id_vertexes[id_vertex]);
							this->boundary_vertexes[id].SetNode(0, this->GetPtrCoordinateViaID(first_boundary[id_type].id_vertexes[id_vertex]));
							id++;
						}
					}
				}

				//2nd boundary
				{
					int N_2st = 0;
					for (int id_type = 0; id_type < second_boundary.size(); id_type++)
					{
						N_2st += second_boundary[id_type].id_vertexes_as_triangle.size();
					}
					this->boundary_faces.resize(N_2st);
					int id = 0;
					for (int id_type = 0; id_type < second_boundary.size(); id_type++)
					{
						for (int id_face = 0; id_face < second_boundary[id_type].id_vertexes_as_triangle.size(); id_face++)
						{
							int id_base_element = second_boundary[id_type].id_base_element[id_face];
							this->boundary_faces[id].SetUpperElementCount(1);
							this->boundary_faces[id].SetUpperElement(0, this->GetElement(id_base_element));

							if (second_boundary[id_type].values.size() == 0) //non individual values for friangles
							{
								this->boundary_faces[id].boundary_value = second_boundary[id_type].value;
							}
							else {
								this->boundary_faces[id].boundary_value = second_boundary[id_type].values[id_face];
							}
							
							for (int id_vertex = 0; id_vertex < this->boundary_faces[id].GetNodesCount(); id_vertex++)
							{
								this->boundary_faces[id].SetIdNode(id_vertex, second_boundary[id_type].id_vertexes_as_triangle[id_face][id_vertex]);
								this->boundary_faces[id].SetNode(id_vertex, this->GetPtrCoordinateViaID(second_boundary[id_type].id_vertexes_as_triangle[id_face][id_vertex]));
							}
							id++;
						}
					}
				}

				CreateFunctionalSpaces();
				CreateExtendedFunctionalSpace_ForCrack();

				/*for (int i = 0; i < this->GetElementsCount(); i++)
				{
					this->GetElement(i)->SetIdDomain(0);
				}*/
			}
			catch (const std::exception&)
			{
				printf_s("Error: MultiXFEM/MultiXFEM_Grid.h/void Initialization(geometry::Grid<geometry::Tetrahedron> &base_grid)\n");
			}
		}
		void Initialization(geometry::Grid<geometry::Tetrahedron> &base_grid, std::vector<std::vector<int>> &first_boundary, std::vector<std::vector<int>> &second_boundary)
		{
			try {
				math::MakeCopyVector_A_into_B(*(base_grid.GetCoordinates()), *(this->GetCoordinates()));
				this->CreateXYZline();

				this->SetElementsCount(base_grid.GetElementsCount());

				for (int id_element = 0; id_element < this->GetElementsCount(); id_element++)
				{
					auto new_element = this->GetElement(id_element);
					auto old_element = base_grid.GetElement(id_element);

					new_element->SetIdNodes(*old_element->GetIdNodes());
					for (int i = 0; i < new_element->GetNodesCount(); i++)
					{
						new_element->SetNode(i, this->GetPtrCoordinateViaID(new_element->GetIdNode(i)));
					}

				}

				this->CreateQTree();

				//create topology
				CreateTopology();

				//1st boundary
				this->boundary_vertexes.resize(first_boundary.size());
				for (int id_vertex = 0; id_vertex < first_boundary.size(); id_vertex++)
				{
					this->boundary_vertexes[id_vertex].SetTypeBoundary(first_boundary[id_vertex][0]);
					this->boundary_vertexes[id_vertex].SetIdNode(0, first_boundary[id_vertex][1]);
					this->boundary_vertexes[id_vertex].SetNode(0, this->GetPtrCoordinateViaID(first_boundary[id_vertex][1]));
				}

				//2nd boundary
			}
			catch (const std::exception&)
			{
				printf_s("Error: MultiXFEM/MultiXFEM_Grid.h/void Initialization(geometry::Grid<geometry::Tetrahedron> &base_grid)\n");
			}
		}
		template <typename Dirichlet, typename Neumann>
		void Initialization(geometry::Grid<geometry::Tetrahedron> &base_grid, std::vector<Dirichlet> &first_boundary, std::vector<Neumann> &second_boundary)
		{
			try {
				math::MakeCopyVector_A_into_B(*(base_grid.GetCoordinates()), *(this->GetCoordinates()));
				this->CreateXYZline();

				this->SetElementsCount(base_grid.GetElementsCount());

				for (int id_element = 0; id_element < this->GetElementsCount(); id_element++)
				{
					auto new_element = this->GetElement(id_element);
					auto old_element = base_grid.GetElement(id_element);

					new_element->SetIdNodes(*old_element->GetIdNodes());
					for (int i = 0; i < new_element->GetNodesCount(); i++)
					{
						new_element->SetNode(i, this->GetPtrCoordinateViaID(new_element->GetIdNode(i)));
					}

				}

				this->CreateQTree();

				//create topology
				CreateTopology();

				//1st boundary
				int N_1st = 0;
				for (int id_type = 0; id_type < first_boundary.size(); id_type++)
				{
					N_1st += first_boundary[id_type].id_vertexes.size();
				}
				this->boundary_vertexes.resize(N_1st);
				int id = 0;
				for (int id_type = 0; id_type < first_boundary.size(); id_type++)
				{
					for (int id_vertex = 0; id_vertex < first_boundary.size(); id_vertex++)
					{
						this->boundary_vertexes[id].boundary_value = first_boundary[id_type].value;
						this->boundary_vertexes[id].SetIdNode(0, first_boundary[id_type].id_vertexes[id_vertex]);
						this->boundary_vertexes[id].SetNode(0, this->GetPtrCoordinateViaID(first_boundary[id_type].id_vertexes[id_vertex]));
						id++;
					}
				}

				//2nd boundary
				int N_2st = 0;
				for (int id_type = 0; id_type < second_boundary.size(); id_type++)
				{
					N_2st += second_boundary[id_type].id_vertexes_as_triangle.size();
				}
				this->boundary_faces.resize(N_2st);
				int id = 0;
				for (int id_type = 0; id_type < second_boundary.size(); id_type++)
				{
					for (int id_face = 0; id_face < second_boundary.size(); id_face++)
					{
						this->boundary_faces[id].boundary_value = second_boundary[id_type].value;
						for (int id_vertex = 0; id_vertex < this->boundary_faces[id].GetNodesCount(); id_vertex++)
						{
							this->boundary_faces[id].SetIdNode(id_vertex, second_boundary[id_type].id_vertexes_as_triangle[id_face][id_vertex]);
							this->boundary_faces[id].SetNode(id_vertex, this->GetPtrCoordinateViaID(second_boundary[id_type].id_vertexes_as_triangle[id_face][id_vertex]));
						}
						id++;
					}
				}
			}
			catch (const std::exception&)
			{
				printf_s("Error: MultiXFEM/MultiXFEM_Grid.h/void Initialization(geometry::Grid<geometry::Tetrahedron> &base_grid)\n");
			}
		}

		template <typename Matrix>
		void CreationPortrait(Matrix &matrix)
		{
			try {
				printf_s("\n");
				std::vector<std::vector<int>> tmp_down_columns(this->GetDOFsCount()), tmp_up_columns(this->GetDOFsCount());
				std::vector<std::vector<int>> down_columns(this->GetDOFsCount()), up_columns(this->GetDOFsCount());
				for (int id_elem = 0; id_elem < this->GetElementsCount(); id_elem++)
				{
					if(id_elem % 1000 == 0)
						printf_s("Work with element[%d]\r", id_elem);
					auto element = this->GetElement(id_elem);
					for (int i = 0; i < element->GetDOFsCount(); i++)
					{
						int global_dof_i = element->GetDOFInLocalID(i);

						for (int j = i + 1; j < element->GetDOFsCount(); j++)
						{
							int global_dof_j = element->GetDOFInLocalID(j);
							if (global_dof_i > global_dof_j) //down elements of matrix
							{
								tmp_down_columns[global_dof_i].push_back(global_dof_j);
								tmp_up_columns[global_dof_j].push_back(global_dof_i);
							}
							else
							{
								tmp_down_columns[global_dof_j].push_back(global_dof_i);
								tmp_up_columns[global_dof_i].push_back(global_dof_j);
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

				matrix.Initialization(up_columns, down_columns);
			}
			catch (const std::exception&)
			{
				printf_s("Error: MultiXFEM/MultiXFEM_Grid.h/void CreationPortrait(Matrix matrix)\n");
			}
		}

		Point<double> GetSolutionInPoint(int id_element, Point<double> X, std::vector<Point<double>> &solution)
		{
			auto element = this->GetElement(id_element);
			Point<double> result;
			for (int j = 0; j < element->GetDOFsCount(); j++)
			{
				auto bf = (*element->GetBasisFunctionInLocalID(j))(X);
				auto value = solution[element->GetDOFInLocalID(j)];
				result.x += bf.x * value.x;
				result.y += bf.y * value.y;
				result.z += bf.z * value.z;
			}
			return result;
		}
		Point<Point<double>> GetDerevativeFromSolutionInPoint(int id_element, Point<double> X, std::vector<Point<double>> &solution)
		{
			auto element = this->GetElement(id_element);
			Point<Point<double>> result;
			for (int j = 0; j < element->GetDOFsCount(); j++)
			{
				auto bf = (*element->GetDerivativeOfBasisFunctionInLocalID(j))(X);
				auto value = solution[element->GetDOFInLocalID(j)];
				result.x += bf.x * value.x;
				result.y += bf.y * value.y;
				result.z += bf.z * value.z;
			}
			return result;
		}
		/// <summary>
			/// EPSILON
			/// </summary>
		Tensor2Rank3D GetStrainTensorFromSolutionInPoint(int id_element, Point<Point<double>> dU)
		{
			Tensor2Rank3D EPS;
			EPS.val[0][0] = dU.x.x;	EPS.val[0][1] = (dU.y.x + dU.x.y) / 2.;	EPS.val[0][2] = (dU.z.x + dU.x.z) / 2.;
			EPS.val[1][0] = EPS.val[0][1];	EPS.val[1][1] = dU.y.y;			EPS.val[1][2] = (dU.z.y + dU.y.z) / 2.;
			EPS.val[2][0] = EPS.val[0][2];	EPS.val[2][1] = EPS.val[1][2];	EPS.val[2][2] = dU.z.z;

			return EPS;
		}
		/// <summary>
		/// EPSILON
		/// </summary>
		Tensor2Rank3D GetStrainTensorFromSolutionInPoint(int id_element, Point<double> X, std::vector<Point<double>>& solution)
		{
			Point<Point<double>> dU = GetDerevativeFromSolutionInPoint(id_element, X, solution);

			return GetStrainTensorFromSolutionInPoint(id_element, dU);
		}
		/// <summary>
		/// SIGMA
		/// </summary>
		Tensor2Rank3D GetStressTensorFromSolutionInPoint(int id_element, Point<double> X, Tensor2Rank3D& StrainTensor)
		{
			Tensor2Rank3D SIG;
			auto _domain = this->GetDomain(this->GetElement(id_element)->GetIdDomain());
			SIG = _domain->forMech.solve_SIGMA(StrainTensor);

			return SIG;
		}
		/// <summary>
		/// SIGMA
		/// </summary>
		Tensor2Rank3D GetStressTensorFromSolutionInPoint(int id_element, Point<double> X, std::vector<Point<double>>& solution)
		{
			Tensor2Rank3D EPS = GetStrainTensorFromSolutionInPoint(id_element, X, solution);
			return GetStressTensorFromSolutionInPoint(id_element, X, EPS);
		}

		void printTecPlot3D(/*char *directory,*/ FILE *fdat, std::vector<std::vector<double>> &value, std::vector<std::vector<char>> name_value, char *name_zone)
		{
			int DEL_domain = 10;
			fprintf_s(fdat, "TITLE     = \"numerical\"\n");
			fprintf_s(fdat, "VARIABLES = \"x\"\n \"y\"\n \"z\"\n");
			for (int i = 0; i < value.size(); i++)
			{
				fprintf_s(fdat, " \"");
				for (int j = 0; j < value[i].size(); j++)
				{
					if (name_value[i][j] != '\0')
					{
						fprintf_s(fdat, "%c", name_value[i][j]);
					}
					else
					{
						break;
					}
				}
				fprintf_s(fdat, "\"\n");
			}

			/*fprintf_s(fdat, " \"sigma_xx\"\n");
			fprintf_s(fdat, " \"sigma_yy\"\n");
			fprintf_s(fdat, " \"sigma_zz\"\n");
			fprintf_s(fdat, " \"eps_xx\"\n");
			fprintf_s(fdat, " \"eps_yy\"\n");
			fprintf_s(fdat, " \"eps_zz\"\n");*/

			int num_elem = this->GetElementsCount();
			for (int i = 0; i < this->GetElementsCount(); i++)
			{
				if (this->GetElement(i)->GetIdDomain() == DEL_domain) num_elem--;
			}
			fprintf_s(fdat, "ZONE T=\"%s\"\n", name_zone);
			fprintf_s(fdat, " N=%d,  E=%d, F=FEBLOCK ET=Tetrahedron \n", this->GetVertexCount(), num_elem);
			fprintf_s(fdat, " VARLOCATION=(NODAL NODAL NODAL");
			for (int i = 0; i < value.size(); i++)
			{
				fprintf_s(fdat, " CELLCENTERED");
			}
			fprintf_s(fdat, ")\n");

			for (int i = 0; i < this->GetVertexCount(); i++)
				fprintf_s(fdat, "%.10e\n", this->GetCoordinateViaID(i).x);
			fprintf_s(fdat, "\n");
			for (int i = 0; i < this->GetVertexCount(); i++)
				fprintf_s(fdat, "%.10e\n", this->GetCoordinateViaID(i).y);
			fprintf_s(fdat, "\n");
			for (int i = 0; i < this->GetVertexCount(); i++)
				fprintf_s(fdat, "%.10e\n", this->GetCoordinateViaID(i).z);
			fprintf_s(fdat, "\n");

			for (int i = 0; i < value.size(); i++)
			{
				for (int j = 0; j < this->GetElementsCount(); j++)
				{
					if (this->GetElement(i)->GetIdDomain() != DEL_domain)
					{
						fprintf_s(fdat, "%.10lf\n", value[i][j]);
					}
				}
				fprintf_s(fdat, "\n");
			}

			for (int i = 0; i < this->GetElementsCount(); i++)
			{
				if (this->GetElement(i)->GetIdDomain() != DEL_domain)
				{
					for (int j = 0; j < this->GetElement(i)->GetNodesCount(); j++)
						fprintf_s(fdat, "%d ", this->GetElement(i)->GetIdNode(j) + 1);
				}
				fprintf_s(fdat, "\n");
			}

			fclose(fdat);
		}
	
		~Grid()
		{
			DOFs_count = 0;
			std::vector<geometry::Crack> v_cr;
			std::vector<geometry::Crack>(v_cr).swap(this->cracks);
			std::vector<BoundaryVertex_forMech> v_b1;
			std::vector<BoundaryVertex_forMech>(v_b1).swap(this->boundary_vertexes);
			std::vector<BoundaryFace_forMech> v_b2;
			std::vector<BoundaryFace_forMech>(v_b2).swap(this->boundary_faces);
			std::vector<int> v_p;
			std::vector<int>(v_p).swap(this->accordance_DOF_and_vertex);
			std::vector<topology::Vertex<topology::lower::EmptyElement, topology::upper::Tetrahedron>> v_v;
			std::vector<topology::Vertex<topology::lower::EmptyElement, topology::upper::Tetrahedron>>(v_v).swap(this->vertexes);

			this->DeleteGeometriGrid();
		}
private:

	};
}