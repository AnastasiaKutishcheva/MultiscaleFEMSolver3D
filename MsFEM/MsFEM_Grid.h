#pragma once
#include "../Library/GeometryGrid.h"
#include "../Library/GeometryShape.h"
#include "../Library/TopologyShape.h"
#include "../Library/FunctionalShape.h"
#include "../Library/DenseMatrix.h"
#include "../Library/Tensor.h"
#include "../Library/Math.h"
#include "../FEM/FEM.h"
#include <stdio.h>
#include <vector>
#include <iostream>

namespace MsFEM {
	class FiniteElement_forMech_Poly :
		public geometry::Polyhedron,
		//public topology::Tetrahedron<topology::lower::Triangle, topology::upper::EmptyElement>,
		public topology::Polyhedron<topology::lower::Polygon, topology::upper::EmptyElement>,
		public functional::Shape<int, Point<double>>
	{
	public:
		math::SimpleGrid SubGrid_for_integration;
		math::SimpleGrid SubGrid_for_integr_byTriangles;

		FEM::Grid_forMech self_grid;
		std::vector< std::vector<Point<double>>> self_basis_functions;

		char self_direction[1000];

		FiniteElement_forMech_Poly() { return; };
		~FiniteElement_forMech_Poly() { return; };

		void SolveLocalMatrix(DenseMatrix<Tensor2Rank3D, Point<double>>& local_matix, std::function<std::vector<std::vector<double>>(Point<double>)>& koefD)
		{
			try {
				local_matix.SetSize(this->GetDOFsCount());

				this->SetIntegrationLaw(4);

				/*if (this->GetIdDomain() == 2)
					printf_s("\nThis elem[%d] is singular\n", this->GetSelfGlobalId());*/

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
							double cohesive_koef = -7.4e+6;
							math::ResizeVector(mult_res, 3, 3);
							mult_res[0][0] = BF_I_inX.x * BF_J_inX.x * cohesive_koef;
							mult_res[0][1] = BF_I_inX.x * BF_J_inX.y * cohesive_koef;
							mult_res[0][2] = BF_I_inX.x * BF_J_inX.z * cohesive_koef;

							mult_res[1][0] = BF_I_inX.y * BF_J_inX.x * cohesive_koef;
							mult_res[1][1] = BF_I_inX.y * BF_J_inX.y * cohesive_koef;
							mult_res[1][2] = BF_I_inX.y * BF_J_inX.z * cohesive_koef;

							mult_res[2][0] = BF_I_inX.z * BF_J_inX.x * cohesive_koef;
							mult_res[2][1] = BF_I_inX.z * BF_J_inX.y * cohesive_koef;
							mult_res[2][2] = BF_I_inX.z * BF_J_inX.z * cohesive_koef;

							result = mult_res;
							//result = 1.0;
							return result;
						};

						double V = this->GetVolume();
						local_matix.A[_bf_local_id_I][_bf_local_id_J] = this->SolveIntegral(StiffnessMatrix);
						local_matix.A[_bf_local_id_J][_bf_local_id_I] = this->SolveIntegral(StiffnessMatrix).T();
					}
				}
				/*for (int _bf_local_id_I = 0; _bf_local_id_I < this->GetDOFsCount(); _bf_local_id_I++)
				{
					for (int _bf_local_id_J = _bf_local_id_I + 1; _bf_local_id_J < this->GetDOFsCount(); _bf_local_id_J++)
					{
						local_matix.A[_bf_local_id_I][_bf_local_id_J] = local_matix.A[_bf_local_id_J][_bf_local_id_I].T();
					}
					for (int i = 0; i < 3; i++)
					{
						for (int j = i; j < 3; j++)
						{
							local_matix.A[_bf_local_id_I][_bf_local_id_I].val[i][j] = local_matix.A[_bf_local_id_I][_bf_local_id_I].val[j][i];
						}
					}
				}*/
			}
			catch (const std::exception&)
			{
				printf_s("Error: XFEM/MultiXFEM_Grid.h/void SolveLocalMatrix(DenseMatrix<Tensor2Rank3D> &local_matix, std::function<std::vector<std::vector<double>>(Point)> &koefD)\n");
			}
		}
		void SolveLocalCohesiveMatrix(DenseMatrix<Tensor2Rank3D, Point<double>>& local_matix, std::function<double(Point<double>)>& koefCohesive)
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
							mult_res[0][0] = BF_I_inX.x * BF_J_inX.x * koefCohesive(X);
							mult_res[0][1] = BF_I_inX.x * BF_J_inX.y * koefCohesive(X);
							mult_res[0][2] = BF_I_inX.x * BF_J_inX.z * koefCohesive(X);

							mult_res[1][0] = BF_I_inX.y * BF_J_inX.x * koefCohesive(X);
							mult_res[1][1] = BF_I_inX.y * BF_J_inX.y * koefCohesive(X);
							mult_res[1][2] = BF_I_inX.y * BF_J_inX.z * koefCohesive(X);

							mult_res[2][0] = BF_I_inX.z * BF_J_inX.x * koefCohesive(X);
							mult_res[2][1] = BF_I_inX.z * BF_J_inX.y * koefCohesive(X);
							mult_res[2][2] = BF_I_inX.z * BF_J_inX.z * koefCohesive(X);

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
	class BoundaryVertex_forMech_Poly :
		public geometry::Vertex,
		//public topology::Vertex<topology::lower::EmptyElement, topology::upper::Segment>,
		public topology::Vertex<topology::lower::EmptyElement, topology::upper::Tetrahedron>,
		public functional::Shape<int, Point<double>>
	{
	public:
		std::function<Point<double>(Point<bool>&)> boundary_value;

		BoundaryVertex_forMech_Poly() { return; };
		~BoundaryVertex_forMech_Poly() { return; };

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
	class BoundaryFace_forMech_Poly :
		public geometry::Polygon,
		public topology::Polygon<topology::lower::Polyline, topology::upper::Polyhedron>,
		public functional::Shape<int, Point<double>>
	{
	public:
		std::function<Point<double>(Point<double>)> boundary_value;

		BoundaryFace_forMech_Poly() { return; };
		~BoundaryFace_forMech_Poly() { return; };

		void SetTypeBoundary(int type)
		{
			this->type_boundary = type;
		}
		int GetTypeBoundary()
		{
			return this->type_boundary;
		}

		void SolveLocalBoundaryVector(std::vector<Point<double>>& local_vector, std::function<Point<double>(Point<double>)>& boundary_value)
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
	class Grid_forMech : public geometry::Grid<FiniteElement_forMech_Poly>
	{
		int DOFs_count;
		struct Face: //2D elements
			public topology::Polygon<topology::lower::Polyline, topology::upper::Polyhedron>,
			public geometry::Polygon,
			public functional::Shape<int, Point<double>>
		{
		public:	
			Face() {};
			~Face() {};
		};
		std::vector<Face> faces;
		struct Edge: //1D elements
			public topology::Polyline<topology::lower::Vertex, topology::upper::Polygon>,
			public geometry::Segment, 
			public functional::Shape<int, Point<double>>
		{
		public:
			Edge() {};
			~Edge() {};
		};
		std::vector<Edge> edges;
		struct Vertex: //0D elements
			public topology::Vertex<topology::lower::EmptyElement, topology::upper::Polyline>,
			public geometry::Vertex,			
			public functional::Shape<int, Point<double>>
		{
		public:
			Vertex() {};
			~Vertex() {};
		};
		std::vector<Vertex> vertexes;
		 
	public:
		std::vector<BoundaryVertex_forMech_Poly> boundary_vertexes;
		std::vector<BoundaryFace_forMech_Poly> boundary_faces;

		bool is_print_logFile;

		std::vector<int> accordance_DOF_and_vertex;

		Grid_forMech()
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

		void CreateTopology()
		{
			try {
				//create elements
				for (int id_element = 0; id_element < this->GetElementsCount(); id_element++)
				{
					this->GetElement(id_element)->SetSelfGlobalId(id_element);
					this->GetElement(id_element)->SetLowerElementCount(this->GetElement(id_element)->GetCountFaces());
				}
				//create faces topology
				{
					struct face_temp
					{
						std::vector<int> num;
						int element, local_in_element;
					public:
						face_temp()
						{
							element = 0; local_in_element = 0;
						}
						face_temp(int element, int local_in_element,
							std::vector<int> &n)
						{
							this->element = element; this->local_in_element = local_in_element;
							math::MakeCopyVector_A_into_B(n, this->num);
						}

						void operator= (face_temp A)
						{
							element = A.element;
							local_in_element = A.local_in_element;
							this->num.resize(A.num.size());
							for (int i = 0; i < A.num.size(); i++)
							{
								this->num[i] = A.num[i];
							};
						}
						bool operator< (face_temp A)
						{
							for (int i = 0; i < num.size(); i++)
							{
								if (num[i] < A.num[i]) return true;
								if (num[i] > A.num[i]) return false;
							}
							return false;
						}

						bool operator>(face_temp A)
						{
							for (int i = 0; i < num.size(); i++)
							{
								if (num[i] > A.num[i]) return true;
								if (num[i] < A.num[i]) return false;
							}
							return false;
						}
						bool operator!= (face_temp A)
						{
							for (int i = 0; i < num.size(); i++)
							{
								if (num[i] != A.num[i]) return true;
							}
							return false;
						}
					};
					struct face
					{
						std::vector<int> num;
						int self_global_id;

						std::vector<int> element, local_in_element;
					public:
						face(face_temp A)
						{
							element.push_back(A.element);
							local_in_element.push_back(A.local_in_element);

							num.resize(A.num.size());
							for (int i = 0; i < num.size(); i++)
								num[i] = A.num[i];
						}
						face()
						{
						}
						void set_elem(face_temp A)
						{
							element.push_back(A.element);
							local_in_element.push_back(A.local_in_element);

						}
					};
					auto MakeRemovalOfDuplication_f = [](std::vector<face_temp>& A, std::vector<face>& B) -> void
					{
						if (A.size() != 0)
						{
							int k = 0, j = 0;

							B.push_back(face());
							B[k] = A[0];
							k++;

							for (int i = 1; i < A.size(); i++)
							{
								if (math::GetConfluence(A[j].num, A[i].num).size() != A[j].num.size()) //не повтор
								{
									j = i;
									B.push_back(face());
									B[k] = A[i];
									k++;
								}
								else {
									B[k-1].local_in_element.push_back(A[i].local_in_element);
									B[k-1].element.push_back(A[i].element);
								}
							}
						}
					};

					std::vector <face_temp> f_temp;
					std::vector <face> faces;
					topology::lower::Polyhedron tmp_3D;
					topology::lower::Polygon tmp_2D;
					
					for (int id_element = 0; id_element < this->GetElementsCount(); id_element++)
					{
						auto element3D = this->GetElement(id_element);

						for (int id_face = 0; id_face < element3D->GetCountFaces(); id_face++)
						{
							int _id_elem = id_element;
							int _id_face_in_elem = id_face;
							std::vector<int> id_in_face;
							element3D->GetGlobalIdFace(id_face, id_in_face);
							//math::MakeQuickSort(id_in_face);
							f_temp.push_back(face_temp(_id_elem, _id_face_in_elem, id_in_face));
						}
					}
					math::MakeQuickSort(f_temp);
					MakeRemovalOfDuplication_f(f_temp, faces);

					this->faces.resize(faces.size());
					for (int id_face = 0; id_face < faces.size(); id_face++)
					{
						this->faces[id_face].SetUpperElementCount((int)faces[id_face].element.size());
						this->faces[id_face].SetSelfGlobalId(id_face);
									
						//add Upper and Lower Elements
						for (int id_volume = 0; id_volume < faces[id_face].element.size(); id_volume++)
						{
							this->faces[id_face].SetUpperElement(id_volume, this->GetElement(faces[id_face].element[id_volume]));
							this->GetElement(faces[id_face].element[id_volume])->SetLowerElement(faces[id_face].local_in_element[id_volume], &this->faces[id_face]);
						}

						//save geometry
						std::vector<Point<double>*> P;
						this->GetPtrCoordinateViaID(faces[id_face].num, P);
						this->faces[id_face].SetGeometry(faces[id_face].num, P);
					}
				}

				//create edge topology
				if(true)
				{
					struct edge_temp
					{
						int num1, num2;
						int element, local_in_element;
						int face, local_in_face;
					public:
						edge_temp()
						{
							element = 0; local_in_element = 0;
							face = 0; local_in_face = 0;
							num1 = 0;
							num2 = 0;
						}
						edge_temp(int element, int local_in_element,
							int face, int local_in_face,
							std::vector<int>& n)
						{
							this->element = element; this->local_in_element = local_in_element;
							this->face = face; this->local_in_face = local_in_face;

							num1 = n[0];
							num2 = n[1];
							if (n[0] > n[1])
							{
								num1 = n[1];
								num2 = n[0];
							}
						}

						void operator= (edge_temp A)
						{
							this->element = A.element; this->local_in_element = A.local_in_element;
							this->face = A.face; this->local_in_face = A.local_in_face;

							num1 = A.num1;
							num2 = A.num2;
						}
						bool operator< (edge_temp A)
						{
							if (num1 < A.num1) return true;
							if (num1 > A.num1) return false;
							if (num1 == A.num1 && num2 < A.num2) return true;
							return false;
						}

						bool operator>(edge_temp A)
						{
							if (num1 > A.num1) return true;
							if (num1 < A.num1) return false;
							if (num1 == A.num1 && num2 > A.num2) return true;
							return false;
						}
						bool operator!= (edge_temp A)
						{
							if (num1 != A.num1 || num2 != A.num2) return true;
							return false;
						}

						bool operator == (edge_temp A)
						{
							if (num1 == A.num1 && num2 == A.num2) return true;

							return false;
						}
					};
					struct edge
					{
						int num1, num2;
						int self_global_id;

						std::vector<int> element, local_in_element;
						std::vector<int> face, local_in_face;
					public:
						edge(edge_temp A)
						{
							//element.push_back(A.element);
							//local_in_element.push_back(A.local_in_element);
							//face.push_back(A.face);
							//local_in_face.push_back(A.local_in_face);
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
							//face.push_back(A.face);
							//local_in_face.push_back(A.local_in_face);
							face.push_back(A.face);
							local_in_face.push_back(A.local_in_face);

						}
					};
					auto MakeRemovalOfDuplication_e = [](std::vector<edge_temp>& A, std::vector<edge>& B) -> void
					{
						if (A.size() != 0)
						{
							int k = 0, j = 0;

							B.push_back(edge());
							B[k] = A[0];
							k++;

							for (int i = 1; i < A.size(); i++)
							{
								if (A[j].num1 != A[i].num1 || A[j].num2 != A[i].num2) //не повтор
								{
									j = i;
									B.push_back(edge());
									B[k] = A[i];
									k++;
								}
								else {
									B[k-1].local_in_face.push_back(A[i].local_in_face);
									B[k-1].face.push_back(A[i].face);
								}
							}
						}
					};

					std::vector <edge_temp> e_temp;
					std::vector <edge> edges_vector;

					for (int id_face = 0; id_face < this->faces.size(); id_face++)
					{
						auto element2D = &this->faces[id_face];
						element2D->SetLowerElementCount(element2D->GetNodesCount());

						int local_id_faces_in_elem_0 = -1;
						auto element3D = this->GetElement(element2D->GetUpperElement(0)->GetSelfGlobalId());
						for(int id = 0; id < element3D->GetLowerElementCount(); id++)
						{
							if (element2D->GetSelfGlobalId() == element3D->GetLowerElement(id)->GetSelfGlobalId())
							{
								local_id_faces_in_elem_0 = id;
								break;
							}
						}

						std::vector<int> id_nodes;
						for (int id_edge = 0; id_edge < element2D->GetLowerElementCount(); id_edge++)
						{
							int _id_elem = -1;
							int _id_face_in_elem = -1;
							int _id_face = id_face;
							int _id_edge_in_face = id_edge;
							element3D->GetGlobalIdEdge(local_id_faces_in_elem_0, id_edge, id_nodes);
							e_temp.push_back(edge_temp(
								_id_elem, _id_face_in_elem,
								_id_face, _id_edge_in_face,
								id_nodes));
						}
					}
					math::MakeQuickSort(e_temp);
					MakeRemovalOfDuplication_e(e_temp, edges_vector);

					//add Upper and Lower Elements
					this->edges.resize(edges_vector.size());
					for (int id_edge = 0; id_edge < edges_vector.size(); id_edge++)
					{
						this->edges[id_edge].SetUpperElementCount((int)edges_vector[id_edge].face.size());
						this->edges[id_edge].SetSelfGlobalId(id_edge);
						
						for (int id_face = 0; id_face < edges_vector[id_edge].face.size(); id_face++)
						{
							this->edges[id_edge].SetUpperElement(id_face, &this->faces[edges_vector[id_edge].face[id_face]]);
							this->faces[edges_vector[id_edge].face[id_face]].SetLowerElement(edges_vector[id_edge].local_in_face[id_face], &this->edges[id_edge]);
						}

						//save geometry
						std::vector<int> id_nodes(2);
						id_nodes[0] = edges_vector[id_edge].num1;
						id_nodes[1] = edges_vector[id_edge].num2;
						std::vector<Point<double>*> P;
						this->GetPtrCoordinateViaID(id_nodes, P);
						this->edges[id_edge].SetGeometry(id_nodes, P);
					}
				}

				//create vertexes topology
				if(true)
				{
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
					auto MakeRemovalOfDuplication_n = [](std::vector<node_temp>& A, std::vector<node>& B) -> void
					{
						if (A.size() != 0)
						{
							int k = 0, j = 0;

							B.push_back(node());
							B[k] = A[0];
							k++;

							for (int i = 1; i < A.size(); i++)
							{
								if (A[j].num != A[i].num) //не повтор
								{
									j = i;
									B.push_back(node());
									B[k] = A[i];
									k++;
								}
								else {
									B[k - 1].local_in_edge.push_back(A[i].local_in_edge);
									B[k - 1].edge.push_back(A[i].edge);
								}
							}
						}
					};

					std::vector <node_temp> n_temp;
					std::vector <node> nodes_vector;

					std::vector <int> id_nodes_in_edge(2);

					for (int id_edge = 0; id_edge < this->edges.size(); id_edge++)
					{
						auto element1D = &this->edges[id_edge];
						auto element2D = element1D->GetUpperElement(0);
						auto element3D = element2D->GetUpperElement(0);

						for(int id_vertex = 0; id_vertex < element1D->GetLowerElementCount(); id_vertex++)
						{
							int _id_elem = -1;
							int _id_face_in_elem = -1;
							int _id_face = -1;
							int _id_edge_in_face = -1;
							int _id_edge = id_edge;
							int _id_vertex_in_edge = id_vertex;
							n_temp.push_back(node_temp(_id_elem, _id_face_in_elem,
								_id_face, _id_edge_in_face,
								_id_edge, _id_vertex_in_edge,
								this->edges[id_edge].GetIdNode(id_vertex)));
						}
					}
					math::MakeQuickSort(n_temp);
					MakeRemovalOfDuplication_n(n_temp, nodes_vector);

					this->vertexes.resize(nodes_vector.size());
					//add Upper and Lower Elements
					for (int id_node = 0; id_node < nodes_vector.size(); id_node++)
					{
						this->vertexes[id_node].SetUpperElementCount((int)nodes_vector[id_node].edge.size());
						this->vertexes[id_node].SetSelfGlobalId(id_node);
						
						for (int id_edge = 0; id_edge < nodes_vector[id_node].edge.size(); id_edge++)
						{
							this->vertexes[id_node].SetUpperElement(id_edge, &this->edges[nodes_vector[id_node].edge[id_edge]]);
							this->edges[nodes_vector[id_node].edge[id_edge]].SetLowerElement(nodes_vector[id_node].local_in_edge[id_edge], &this->vertexes[id_node]);
						}

						//save geometry
						std::vector<int> id_nodes(1);
						id_nodes[0] = nodes_vector[id_node].num;
						std::vector<Point<double>*> P;
						this->GetPtrCoordinateViaID(id_nodes, P);
						this->vertexes[id_node].SetGeometry(id_nodes, P);
					}
				}

				//sort local_edges for faces
				for (int id_face = 0; id_face < faces.size(); id_face++)
				{
					std::vector<std::vector<int>> old_positions(this->faces[id_face].GetLowerElementCount());
					std::vector<int> new_positions(this->faces[id_face].GetLowerElementCount());
					for (int i = 0; i < old_positions.size(); i++)
					{
						old_positions[i].resize(2);
						old_positions[i][0] = this->faces[id_face].GetLowerElement(i)->GetSelfGlobalId();
						old_positions[i][1] = 0;
					}
					new_positions[0] = old_positions[0][0];
					old_positions[0][1] = 1;
					int id_target_vertex = this->edges[new_positions[0]].GetLowerElement(1)->GetSelfGlobalId();
					for (int i = 1; i < new_positions.size(); i++)
					{
						for (int j = 1; j < old_positions.size(); j++)
						{
							if (old_positions[j][1] == 0)
							{
								if (this->edges[old_positions[j][0]].GetLowerElement(0)->GetSelfGlobalId() == id_target_vertex)
								{
									id_target_vertex = this->edges[old_positions[j][0]].GetLowerElement(1)->GetSelfGlobalId();
									new_positions[i] = old_positions[j][0];
									old_positions[j][1] = 1;
									break;
								}
								else if (this->edges[old_positions[j][0]].GetLowerElement(1)->GetSelfGlobalId() == id_target_vertex) {
									id_target_vertex = this->edges[old_positions[j][0]].GetLowerElement(0)->GetSelfGlobalId();
									new_positions[i] = old_positions[j][0];
									old_positions[j][1] = 1;
									break;
								}
							}
						}
					}

					for (int i = 0; i < new_positions.size(); i++)
					{
						this->faces[id_face].SetLowerElement(i, &this->edges[new_positions[i]]);
					}
				}

				
			}
			catch (const std::exception&)
			{
				printf_s("Error: MultiXFEM/MultiXFEM_Grid.h/void CreateTopology()\n");
			}
		}
		void CreateFunctionalSpaces(char* fine_mesh_dir)
		{
			//create basis functions
			this->SetDOFsCount((int)this->vertexes.size() + this->edges.size() + this->faces.size() + this->GetElementsCount());
			for (int id_element = 0; id_element < this->GetElementsCount(); id_element++)
			{
				auto element = this->GetElement(id_element);
				sprintf_s(element->self_direction, sizeof(element->self_direction), "%s/Elem_%d/SubElem_0", fine_mesh_dir, id_element);
								
				std::vector<int> id_self_edges;
				{
					std::vector<int> id_self_edges_tmp;
					for (int id_face_local = 0; id_face_local < element->GetLowerElementCount(); id_face_local++)
					{
						int id_face = element->GetLowerElement(id_face_local)->GetSelfGlobalId();
						for (int id_edge_local = 0; id_edge_local < this->faces[id_face].GetLowerElementCount(); id_edge_local++)
						{
							int id_edge = this->faces[id_face].GetLowerElement(id_edge_local)->GetSelfGlobalId();
							id_self_edges_tmp.push_back(id_edge);
						}
					}
					math::MakeQuickSort(id_self_edges_tmp);
					math::MakeRemovalOfDuplication(id_self_edges_tmp, id_self_edges);
				}
				std::vector<std::vector<int>> boundary_vertexes_in_macro(id_self_edges.size());

				int DOF_count = element->GetNodesCount(); //vertexes
				DOF_count += 1; //self volume
				DOF_count += element->GetLowerElementCount(); //faces
				DOF_count += id_self_edges.size(); //edges
				element->ResizeDOF(DOF_count);
				element->self_basis_functions.resize(DOF_count);
					
				//read sub meshes
				math::SimpleGrid base_grid_simple;
				base_grid_simple.ReadFromNVTR(element->self_direction, 4);
				std::vector<math::SimpleGrid> faces_grid_simple(element->GetLowerElementCount());
				std::vector<math::SimpleGrid> edges_grid_simple(id_self_edges.size());
				{
					for (int id_face_local = 0; id_face_local < element->GetLowerElementCount(); id_face_local++)
					{
						char name_in_file[1000];
						sprintf_s(name_in_file, "%s/Face_id%d.dat", element->self_direction, id_face_local);
						faces_grid_simple[id_face_local].ReadFromSalomeDat(name_in_file, 2, this->GetVertexCount());
					}

					for (int id_edge_local = 0; id_edge_local < id_self_edges.size(); id_edge_local++)
					{
						int id_edge = id_self_edges[id_edge_local];
						std::vector<int> neigbor_faces_local_in_elem;
						for (int id_face_local = 0; id_face_local < this->edges[id_edge].GetUpperElementCount() && neigbor_faces_local_in_elem.size() < 2; id_face_local++)
						{
							for (int id_element_local = 0; id_element_local < this->edges[id_edge].GetUpperElement(id_face_local)->GetUpperElementCount(); id_element_local++)
							{
								int id_element = this->edges[id_edge].GetUpperElement(id_face_local)->GetUpperElement(id_element_local)->GetSelfGlobalId();
								if (id_element == element->GetSelfGlobalId())
								{
									int global_id = this->edges[id_edge].GetUpperElement(id_face_local)->GetSelfGlobalId();
									int local_id = -1;
									for (int i = 0; i < element->GetLowerElementCount(); i++)
									{
										if (element->GetLowerElement(i)->GetSelfGlobalId() == global_id)
										{
											local_id = i;
											break;
										}
									}
									neigbor_faces_local_in_elem.push_back(local_id);
									break;
								}
							}
						}

						for (int id_elem_f1 = 0; id_elem_f1 < faces_grid_simple[neigbor_faces_local_in_elem[0]].nvtr.size(); id_elem_f1++)
						{
							for (int id_elem_f2 = 0; id_elem_f2 < faces_grid_simple[neigbor_faces_local_in_elem[1]].nvtr.size(); id_elem_f2++)
							{
								auto segment = math::GetConfluence(faces_grid_simple[neigbor_faces_local_in_elem[0]].nvtr[id_elem_f1], faces_grid_simple[neigbor_faces_local_in_elem[1]].nvtr[id_elem_f2]);

								if (segment.size() == 2)
								{
									std::vector<Point<double>> X(2);
									X[0] = base_grid_simple.xyz[segment[0]];
									X[1] = base_grid_simple.xyz[segment[1]];

									if (X[1] < X[0])
									{
										Point<double> tmp = X[0];
										X[0] = X[1];
										X[1] = tmp;

										int tmp2 = segment[0];
										segment[0] = segment[1];
										segment[1] = tmp2;
									}

									edges_grid_simple[id_edge_local].nvtr.push_back(segment);
									edges_grid_simple[id_edge_local].nvkat.push_back(faces_grid_simple[neigbor_faces_local_in_elem[0]].nvkat[id_elem_f1]);
									break;
								}
								if (segment.size() > 2)
								{
									printf_s("\nError in create segment!!!\n");
									//sleep(100000);
								}
							}
						}

						//create map for edge
						std::vector<int> vertex_in_edge;
						vertex_in_edge.reserve(edges_grid_simple[id_edge_local].nvtr.size() * 2);
						for (int i = 0; i < edges_grid_simple[id_edge_local].nvtr.size(); i++)
						{
							vertex_in_edge.push_back(edges_grid_simple[id_edge_local].nvtr[i][0]);
							vertex_in_edge.push_back(edges_grid_simple[id_edge_local].nvtr[i][1]);
						}
						math::MakeQuickSort(vertex_in_edge);
						math::MakeRemovalOfDuplication(vertex_in_edge, edges_grid_simple[id_edge_local].vertex_map);
						edges_grid_simple[id_edge_local].xyz.resize(edges_grid_simple[id_edge_local].vertex_map.size());
						for (int i = 0; i < edges_grid_simple[id_edge_local].xyz.size(); i++)
						{
							edges_grid_simple[id_edge_local].xyz[i] = base_grid_simple.xyz[edges_grid_simple[id_edge_local].vertex_map[i]];
						}
					}

					//renumeration of the NVTR into self
					for (int id_face_local = 0; id_face_local < faces_grid_simple.size(); id_face_local++)
					{
						for (int i = 0; i < faces_grid_simple[id_face_local].nvtr.size(); i++)
						{
							for (int j = 0; j < faces_grid_simple[id_face_local].nvtr[i].size(); j++)
							{
								faces_grid_simple[id_face_local].nvtr[i][j] = faces_grid_simple[id_face_local].TransferIdGlobalIntoSelf(faces_grid_simple[id_face_local].nvtr[i][j]);
							}
						}
					}
					for (int id_edge_local = 0; id_edge_local < edges_grid_simple.size(); id_edge_local++)
					{
						for (int i = 0; i < edges_grid_simple[id_edge_local].nvtr.size(); i++)
						{
							for (int j = 0; j < edges_grid_simple[id_edge_local].nvtr[i].size(); j++)
							{
								edges_grid_simple[id_edge_local].nvtr[i][j] = edges_grid_simple[id_edge_local].TransferIdGlobalIntoSelf(edges_grid_simple[id_edge_local].nvtr[i][j]);
							}
						}
					}
				}

				//create basis funtions
				bool is_solve = false;
				
				for (int id_bf = 0; id_bf < element->GetDOFsCount(); id_bf++)
				{
					FILE* fin;
					char name_in[1000];
					sprintf_s(name_in, "%s/BasisFunction_%d.txt", element->self_direction, id_bf);
					fopen_s(&fin, name_in, "r");
					if (fin == NULL)
					{
						is_solve = true;
						break;
					}
					else {
						fclose(fin);
					}
				}
				if (is_solve == false)
				{
					for (int i = 0; i < this->GetDomainsCount(); i++)
						element->self_grid.AddDomain(*this->GetDomain(i));
					struct _Dirichlet {
						int global_id, in_elem_id;
						std::vector<int> id_vertexes;
						std::vector<Point<double>> xyz;
						std::function<Point<double>(Point<bool>&, int id_vertex)> value;
					};
					std::vector<_Dirichlet> first_boundaries_empty;
					struct _Neumann {
						/*double value;
						Point<double> vector;*/
						std::function<Point<double>(Point<double>)> value;
						std::vector<std::function<Point<double>(Point<double>)>> values;
						std::vector<std::vector<int>> id_vertexes_as_triangle;
						std::vector<int> id_base_element;
					};
					std::vector<_Neumann> second_boundary_empty;

					element->self_grid.Initialization(base_grid_simple, first_boundaries_empty, second_boundary_empty);

					for (int id_bf = 0; id_bf < element->GetDOFsCount(); id_bf++)
					{
						FILE* fin;
						char name_in[1000];
						sprintf_s(name_in, "%s/BasisFunction_%d.txt", element->self_direction, id_bf);
						fopen_s(&fin, name_in, "r");

						element->self_basis_functions[id_bf].resize(element->self_grid.GetDOFsCount());

						for (int i = 0; i < element->self_grid.GetDOFsCount(); i++)
						{
							fscanf_s(fin, "%lf %lf %lf", &element->self_basis_functions[id_bf][i].x, &element->self_basis_functions[id_bf][i].y, &element->self_basis_functions[id_bf][i].z);
						}

						fclose(fin);
					}
				}
				if(is_solve)
				{
					//solve edges problems
					std::vector<FEM::Grid_1D_forMech> solver_grid_edges(id_self_edges.size()); //output [num_edge]
					std::vector<std::vector<std::vector<Point<double>>>> Solution_edges(id_self_edges.size()); //output [num_edge][num_BF][num_DOF]
					std::vector<std::vector<int>> boundary_vertexes_local(id_self_edges.size());
					for (int id_edge_local = 0; id_edge_local < id_self_edges.size(); id_edge_local++)
					{
						Solution_edges[id_edge_local].resize(3);

						for (int i = 0; i < this->GetDomainsCount(); i++)
							solver_grid_edges[id_edge_local].AddDomain(*this->GetDomain(i));


						CSSD_Matrix<Tensor2Rank3D, Point<double>> Base_stiffness_matrix;
						Base_stiffness_matrix.print_logs = this->is_print_logFile;

						//create Boundary values and vertexes
						struct _Dirichlet {
							std::vector<int> id_vertexes;
							std::function<Point<double>(Point<bool>&)> value;
						};
						std::vector<int> elements_for_vertex(edges_grid_simple[id_edge_local].xyz.size());
						for (int i = 0; i < edges_grid_simple[id_edge_local].nvtr.size(); i++)
						{
							elements_for_vertex[edges_grid_simple[id_edge_local].nvtr[i][0]]++;
							elements_for_vertex[edges_grid_simple[id_edge_local].nvtr[i][1]]++;
						}
						for (int i = 0; i < elements_for_vertex.size(); i++)
						{
							if (elements_for_vertex[i] == 1)
							{
								boundary_vertexes_local[id_edge_local].push_back(edges_grid_simple[id_edge_local].TransferIdSelfIntoGlobal(i));
							}
						}
						std::vector<std::vector<int>> vertex_in_boundaries(2);
						std::vector<_Dirichlet> first_boundaries(2);
						first_boundaries[0].id_vertexes.push_back(edges_grid_simple[id_edge_local].TransferIdGlobalIntoSelf(boundary_vertexes_local[id_edge_local][0]));
						first_boundaries[1].id_vertexes.push_back(edges_grid_simple[id_edge_local].TransferIdGlobalIntoSelf(boundary_vertexes_local[id_edge_local][1]));
						boundary_vertexes_in_macro[id_edge_local].resize(2);
						if (math::IsEqual(edges_grid_simple[id_edge_local].xyz[first_boundaries[0].id_vertexes[0]], this->edges[id_self_edges[id_edge_local]].GetNode(0)))
						{
							boundary_vertexes_in_macro[id_edge_local][0] = this->edges[id_self_edges[id_edge_local]].GetIdNode(0);
							boundary_vertexes_in_macro[id_edge_local][1] = this->edges[id_self_edges[id_edge_local]].GetIdNode(1);
						}
						else {
							boundary_vertexes_in_macro[id_edge_local][1] = this->edges[id_self_edges[id_edge_local]].GetIdNode(0);
							boundary_vertexes_in_macro[id_edge_local][0] = this->edges[id_self_edges[id_edge_local]].GetIdNode(1);
						}

						//First function (linear)
						printf_s("create linear function for edge %d(local: %d) - bf linear 0\n", id_self_edges[id_edge_local], id_edge_local);
						{
							first_boundaries[0].value = [](Point<bool>& arg)->Point<double>
							{
								arg.x = true;
								arg.y = true;
								arg.z = true;
								Point<double> result;
								result.x = 1;
								result.y = 1;
								result.z = 1;
								return result;
							};
							first_boundaries[1].value = [](Point<bool>& arg)->Point<double>
							{
								arg.x = true;
								arg.y = true;
								arg.z = true;
								Point<double> result;
								result.x = 0;
								result.y = 0;
								result.z = 0;
								return result;
							};

							struct _Neumann {
								/*double value;
								Point<double> vector;*/
								std::function<Point<double>(Point<double>)> value;
								std::vector<std::function<Point<double>(Point<double>)>> values;
								std::vector<std::vector<int>> id_vertexes_as_triangle;
								std::vector<int> id_base_element;
							};
							std::vector<_Neumann> second_boundary_empty;

							//right_side
							std::function<Point<double>(bool& is_solve, int id_element, Point<double> X)> sourse = [](bool& is_solve, int id_element, Point<double> X) -> Point<double>
							{
								is_solve = false;
								return Point<double>(0, 0, 0);
							};

							FEM::FEM_1D_forElasticDeformation(
								false,
								1e-15,
								edges_grid_simple[id_edge_local], //input
								first_boundaries, //input
								second_boundary_empty, //input
								sourse,
								element->self_direction, //output
								solver_grid_edges[id_edge_local], //output
								Solution_edges[id_edge_local][0], //output
								Base_stiffness_matrix //output
							);
						}

						//Second function (linear)
						printf_s("create linear function for edge %d(local: %d) - bf linear 1\n", id_self_edges[id_edge_local], id_edge_local);
						{
							first_boundaries[0].value = [](Point<bool>& arg) -> Point<double>
							{
								arg.x = true;
								arg.y = true;
								arg.z = true;
								Point<double> result;
								result.x = 0;
								result.y = 0;
								result.z = 0;
								return result;
							};
							first_boundaries[1].value = [](Point<bool>& arg) -> Point<double>
							{
								arg.x = true;
								arg.y = true;
								arg.z = true;
								Point<double> result;
								result.x = 1;
								result.y = 1;
								result.z = 1;
								return result;
							};

							struct _Neumann {
								/*double value;
								Point<double> vector;*/
								std::function<Point<double>(Point<double>)> value;
								std::vector<std::function<Point<double>(Point<double>)>> values;
								std::vector<std::vector<int>> id_vertexes_as_triangle;
								std::vector<int> id_base_element;
							};
							std::vector<_Neumann> second_boundary_empty;

							//right_side
							std::function<Point<double>(bool& is_solve, int id_element, Point<double> X)> sourse = [](bool& is_solve, int id_element, Point<double> X) -> Point<double>
							{
								is_solve = false;
								return Point<double>(0, 0, 0);
							};

							FEM::FEM_1D_forElasticDeformation(
								false,
								1e-15,
								edges_grid_simple[id_edge_local], //input
								first_boundaries, //input
								second_boundary_empty, //input
								sourse,
								element->self_direction, //output
								solver_grid_edges[id_edge_local], //output
								Solution_edges[id_edge_local][1], //output
								Base_stiffness_matrix //output
							);
						}

						//Third function (bubble)
						printf_s("create bubble function for edge %d(local: %d) - bf bubble\n", id_self_edges[id_edge_local], id_edge_local);
						{
							first_boundaries[0].value = [](Point<bool>& arg) -> Point<double>
							{
								arg.x = true;
								arg.y = true;
								arg.z = true;
								return Point<double>(0, 0, 0);
							};
							first_boundaries[1].value = [](Point<bool>& arg) -> Point<double>
							{
								arg.x = true;
								arg.y = true;
								arg.z = true;
								return Point<double>(0, 0, 0);
							};

							struct _Neumann {
								/*double value;
								Point<double> vector;*/
								std::function<Point<double>(Point<double>)> value;
								std::vector<std::function<Point<double>(Point<double>)>> values;
								std::vector<std::vector<int>> id_vertexes_as_triangle;
								std::vector<int> id_base_element;
							};
							std::vector<_Neumann> second_boundary_empty;

							auto BaseEdge = &this->edges[id_self_edges[id_edge_local]];

							//right_side
							std::function<Point<double>(bool& is_solve, int id_element, Point<double> X)> sourse = [&solver_grid_edges, BaseEdge, id_edge_local, this](bool& is_solve, int id_element, Point<double> X) -> Point<double>
							{
								is_solve = true;
								Point<double> result;
								Point<double> self_result;

								Point<double> X0_self = BaseEdge->GetSelfNode(0);
								Point<double> X1_self = BaseEdge->GetSelfNode(1);
								double size_elem = math::SolveLengthVector(X0_self, X1_self);

								double lambda, mu;
								this->GetDomain(solver_grid_edges[id_edge_local].GetElement(id_element)->GetIdDomain())->forMech.GetLameCoefficients(lambda, mu);

								auto basis = *BaseEdge->GetSelfReverseBasis();
								double laplasian_phi = -8 / (size_elem * size_elem) * (basis[0][0] * basis[0][0] + basis[0][1] * basis[0][1] + basis[0][2] * basis[0][2]);
								Point<double> DD_phi;
								DD_phi.x = -8 / (size_elem * size_elem) * (basis[0][0] * basis[0][0] + basis[0][0] * basis[0][1] + basis[0][0] * basis[0][2]);
								DD_phi.y = -8 / (size_elem * size_elem) * (basis[0][0] * basis[0][1] + basis[0][1] * basis[0][1] + basis[0][1] * basis[0][2]);
								DD_phi.z = -8 / (size_elem * size_elem) * (basis[0][0] * basis[0][2] + basis[0][2] * basis[0][1] + basis[0][2] * basis[0][2]);

								result = DD_phi * (lambda + mu) + Point<double>(laplasian_phi, laplasian_phi, laplasian_phi) * mu;

								return  result * (-1.0);
							};

							FEM::FEM_1D_forElasticDeformation(
								false,
								1e-15,
								edges_grid_simple[id_edge_local], //input
								first_boundaries, //input
								second_boundary_empty, //input
								sourse,
								element->self_direction, //output
								solver_grid_edges[id_edge_local], //output
								Solution_edges[id_edge_local][2], //output
								Base_stiffness_matrix //output
							);
						}
					}
				
					//solve faces problems
					std::vector<FEM::Grid_2D_forMech> solver_grid_faces(faces_grid_simple.size()); //output [num_face]
					std::vector<std::vector<std::vector<Point<double>>>> Solution_faces(faces_grid_simple.size()); //output [num_face][num_BF][num_DOF]
					std::vector<std::vector<int>> target_id_for_BF(solver_grid_faces.size());
					for (int id_face_local = 0; id_face_local < solver_grid_faces.size(); id_face_local++)
					{
						int id_face = element->GetLowerElement(id_face_local)->GetSelfGlobalId();
						Solution_faces[id_face_local].resize(faces[id_face].GetLowerElementCount() * 2 + 1);

						for (int i = 0; i < this->GetDomainsCount(); i++)
							solver_grid_faces[id_face_local].AddDomain(*this->GetDomain(i));

						CSSD_Matrix<Tensor2Rank3D, Point<double>> Base_stiffness_matrix;
						Base_stiffness_matrix.print_logs = this->is_print_logFile;

						//create Boundary values and vertexes
						printf_s("create Boundary values and vertexes for face %d(local: %d)\n", id_face, id_face_local);
						struct _Dirichlet {
							int global_id, in_elem_id;
							std::vector<int> id_vertexes;
							std::vector<Point<double>> xyz;
							std::function<Point<double>(Point<bool>&, int id_vertex)> value;
						};
						std::vector<_Dirichlet> first_boundaries(faces[id_face].GetLowerElementCount());
						for (int id_edge_local = 0; id_edge_local < first_boundaries.size(); id_edge_local++)
						{
							first_boundaries[id_edge_local].global_id = faces[id_face].GetLowerElement(id_edge_local)->GetSelfGlobalId();
							first_boundaries[id_edge_local].in_elem_id = math::GetPositionInSortVector(id_self_edges, first_boundaries[id_edge_local].global_id);

							for (int i = 0; i < edges_grid_simple[first_boundaries[id_edge_local].in_elem_id].xyz.size(); i++)
							{
								int global_id_vertex = edges_grid_simple[first_boundaries[id_edge_local].in_elem_id].TransferIdSelfIntoGlobal(i);
								first_boundaries[id_edge_local].id_vertexes.push_back(faces_grid_simple[id_face_local].TransferIdGlobalIntoSelf(global_id_vertex));
								//first_boundaries[id_edge_local].xyz.push_back(faces_grid_simple[id_face_local].xyz[first_boundaries[id_edge_local].id_vertexes[first_boundaries[id_edge_local].id_vertexes.size() - 1]]);
							}
						}

						//linear function
						for (int id_bf = 0; id_bf < faces[id_face].GetLowerElementCount(); id_bf++)
						{
							printf_s("create linear functions for face %d(local: %d) - bf linear %d\n", id_face, id_face_local, id_bf);
							for (int i = 0; i < first_boundaries.size(); i++)
							{
								first_boundaries[i].value = [](Point<bool>& arg, int id) -> Point<double>
								{
									arg.x = true;
									arg.y = true;
									arg.z = true;
									return Point<double>(0, 0, 0);
								};
							}

							int left_edge_id_global = this->faces[id_face].GetLowerElement((id_bf - 1) < 0 ? faces[id_face].GetLowerElementCount() - 1 : id_bf - 1)->GetSelfGlobalId();
							int left_edge_id = math::GetPositionInSortVector(id_self_edges, left_edge_id_global);

							int right_edge_id_global = this->faces[id_face].GetLowerElement(id_bf)->GetSelfGlobalId();
							int right_edge_id = math::GetPositionInSortVector(id_self_edges, right_edge_id_global);

							int target_vertex_global = math::GetConfluence(boundary_vertexes_in_macro[left_edge_id], boundary_vertexes_in_macro[right_edge_id])[0];
							target_id_for_BF[id_face_local].push_back(target_vertex_global);

							first_boundaries[id_bf].value = [&Solution_edges, &boundary_vertexes_in_macro, &edges_grid_simple, &faces_grid_simple, id_face_local, right_edge_id, target_vertex_global](Point<bool>& arg, int id_vertex) -> Point<double>
							{
								Point<double> result;
								arg.x = true;
								arg.y = true;
								arg.z = true;

								int id_bf_in_edge;
								if (target_vertex_global == boundary_vertexes_in_macro[right_edge_id][0]) id_bf_in_edge = 0;
								else id_bf_in_edge = 1;

								int id_vertex_global = faces_grid_simple[id_face_local].TransferIdSelfIntoGlobal(id_vertex);

								result = Solution_edges[right_edge_id][id_bf_in_edge][edges_grid_simple[right_edge_id].TransferIdGlobalIntoSelf(id_vertex_global)];

								return result;
							};
							first_boundaries[(id_bf - 1) < 0 ? faces[id_face].GetLowerElementCount() - 1 : id_bf - 1].value = [&Solution_edges, &boundary_vertexes_in_macro, &edges_grid_simple, &faces_grid_simple, id_face_local, left_edge_id, target_vertex_global](Point<bool>& arg, int id_vertex) -> Point<double>
							{
								Point<double> result;
								arg.x = true;
								arg.y = true;
								arg.z = true;

								int id_bf_in_edge;
								if (target_vertex_global == boundary_vertexes_in_macro[left_edge_id][0]) id_bf_in_edge = 0;
								else id_bf_in_edge = 1;

								int id_vertex_global = faces_grid_simple[id_face_local].TransferIdSelfIntoGlobal(id_vertex);

								result = Solution_edges[left_edge_id][id_bf_in_edge][edges_grid_simple[left_edge_id].TransferIdGlobalIntoSelf(id_vertex_global)];

								return result;
							};

							struct _Neumann {
								/*double value;
								Point<double> vector;*/
								std::function<Point<double>(Point<double>)> value;
								std::vector<std::function<Point<double>(Point<double>)>> values;
								std::vector<std::vector<int>> id_vertexes_as_triangle;
								std::vector<int> id_base_element;
							};
							std::vector<_Neumann> second_boundary_empty;

							//right_side
							std::function<Point<double>(bool& is_solve, int id_element, Point<double> X)> sourse = [](bool& is_solve, int id_element, Point<double> X) -> Point<double>
							{
								is_solve = false;
								return Point<double>(0, 0, 0);
							};

							FEM::FEM_2D_forElasticDeformation(
								false,
								1e-15,
								faces_grid_simple[id_face_local], //input
								first_boundaries, //input
								second_boundary_empty, //input
								sourse,
								element->self_direction, //output
								solver_grid_faces[id_face_local], //output
								Solution_faces[id_face_local][id_bf], //output
								Base_stiffness_matrix //output
							);

							if(this->is_print_logFile)
							{
								printf_s("Print the mech result into .dat file... \n");

								FILE* fout_tech;
								char name_u_tech[5000];
								sprintf_s(name_u_tech, "%s/Face_%d(global%d)_linear_bf%d.dat", element->self_direction, id_face_local, id_face, id_bf);
								fopen_s(&fout_tech, name_u_tech, "w");
								char name_in_file[1000];
								sprintf_s(name_in_file, "Face_%d(global%d)_linear_bf%d", id_face_local, id_face, id_bf);
								std::vector<std::vector<char>> name_value(4);
								char name_v_tmp[6][100];
								sprintf_s(name_v_tmp[0], "Ux");
								sprintf_s(name_v_tmp[1], "Uy");
								sprintf_s(name_v_tmp[2], "Uz");
								sprintf_s(name_v_tmp[3], "Material");
								for (int i = 0; i < name_value.size(); i++)
								{
									name_value[i].resize(100);
									for (int j = 0; j < name_value[i].size(); j++)
									{
										name_value[i][j] = name_v_tmp[i][j];
									}
								}
								std::vector<std::vector<double>> value(4);
								value[0].resize(solver_grid_faces[id_face_local].GetVertexCount());
								value[1].resize(solver_grid_faces[id_face_local].GetVertexCount());
								value[2].resize(solver_grid_faces[id_face_local].GetVertexCount());
								value[3].resize(solver_grid_faces[id_face_local].GetElementsCount());

								for (int i = 0; i < solver_grid_faces[id_face_local].GetVertexCount(); i++)
								{
									value[0][i] = Solution_faces[id_face_local][id_bf][i].x;
									value[1][i] = Solution_faces[id_face_local][id_bf][i].y;
									value[2][i] = Solution_faces[id_face_local][id_bf][i].z;
								}
								for (int i = 0; i < solver_grid_faces[id_face_local].GetElementsCount(); i++)
								{
									auto element = solver_grid_faces[id_face_local].GetElement(i);
									value[3][i] = element->GetIdDomain();

								}
								solver_grid_faces[id_face_local].printTecPlot3D_DiffDomains(fout_tech, value, name_value, name_in_file);
								fclose(fout_tech);
							}
						}

						//bubble functions by edges
						Point<double> centr_of_polygon = this->faces[id_face].GetWeightCentr();
						for (int id_bf = 0; id_bf < faces[id_face].GetLowerElementCount(); id_bf++)
						{
							printf_s("create bubble functions by edges for face %d(local: %d) - bf edge_bubble %d\n", id_face, id_face_local, id_bf);

							for (int i = 0; i < first_boundaries.size(); i++)
							{
								first_boundaries[i].value = [](Point<bool>& arg, int id) -> Point<double>
								{
									arg.x = true;
									arg.y = true;
									arg.z = true;
									return Point<double>(0, 0, 0);
								};
							}

							int target_edge_id_global = this->faces[id_face].GetLowerElement(id_bf)->GetSelfGlobalId();
							int target_edge_id = math::GetPositionInSortVector(id_self_edges, target_edge_id_global);

							target_id_for_BF[id_face_local].push_back(target_edge_id_global);

							first_boundaries[id_bf].value = [&Solution_edges, &edges_grid_simple, &faces_grid_simple, id_face_local, target_edge_id](Point<bool>& arg, int id_vertex) -> Point<double>
							{
								Point<double> result;
								arg.x = true;
								arg.y = true;
								arg.z = true;

								int id_vertex_global = faces_grid_simple[id_face_local].TransferIdSelfIntoGlobal(id_vertex);

								result = Solution_edges[target_edge_id][2][edges_grid_simple[target_edge_id].TransferIdGlobalIntoSelf(id_vertex_global)];

								return result;
							};
							struct _Neumann {
								/*double value;
								Point<double> vector;*/
								std::function<Point<double>(Point<double>)> value;
								std::vector<std::function<Point<double>(Point<double>)>> values;
								std::vector<std::vector<int>> id_vertexes_as_triangle;
								std::vector<int> id_base_element;
							};
							std::vector<_Neumann> second_boundary_empty;


							int id_edge_local = id_bf;
							Point<double> Edge[2];
							/*Edge[0] = base_grid_simple.xyz[boundary_vertexes[target_edge_id][0]];
							Edge[1] = base_grid_simple.xyz[boundary_vertexes[target_edge_id][1]];*/
							Edge[0] = this->edges[target_edge_id_global].GetNode(0);
							Edge[1] = this->edges[target_edge_id_global].GetNode(1);

							//делаем проекцию всех узлов на плоскость (Edge[0], Edge[1], centr_of_polygon)
							Point<double> plane[3];
							plane[0] = Edge[0];
							plane[1] = Edge[1];
							plane[2] = centr_of_polygon;
							//разворачиваем базис
							geometry::Triangle test_triag;
							test_triag.SetNode(0, &plane[0]);
							test_triag.SetNode(1, &plane[1]);
							test_triag.SetNode(2, &plane[2]);
							test_triag.CreateBasis();
							Point<double> Edge_self[2];
							Edge_self[0] = test_triag.MakeTransferXyzIntoSelf(Edge[0]);
							Edge_self[1] = test_triag.MakeTransferXyzIntoSelf(Edge[1]);
							//делаем поправочный вектор
							//Point<double> correction_vector = test_triag.MakeTransferXyzIntoSelf(plane[0]) * -1.;
							Point<double> correction_vector;

							std::vector<Point<double>> vertex_in_plane(this->faces[id_face].GetNodesCount());
							for (int i = 0; i < vertex_in_plane.size(); i++)
							{
								//vertex_in_plane[i] = math::MakeProgectionOfPointIntoLine(plane, this->faces[id_face].GetNode(i));
								vertex_in_plane[i] = test_triag.MakeTransferXyzIntoSelf(this->faces[id_face].GetNode(i)) + correction_vector;
							}

							//строим прямоугольник
							Point<double> p_min = vertex_in_plane[0];
							Point<double> p_max = vertex_in_plane[0];
							for (int i = 0; i < vertex_in_plane.size(); i++)
							{
								//Point<double> p_in_self = test_triag.MakeTransferXyzIntoSelf(vertex_in_plane[i])+ correction_vector;
								if (vertex_in_plane[i].x < p_min.x) p_min.x = vertex_in_plane[i].x;
								if (vertex_in_plane[i].y < p_min.y) p_min.y = vertex_in_plane[i].y;
								if (vertex_in_plane[i].z < p_min.z) p_min.z = vertex_in_plane[i].z;
								if (vertex_in_plane[i].x > p_max.x) p_max.x = vertex_in_plane[i].x;
								if (vertex_in_plane[i].y > p_max.y) p_max.y = vertex_in_plane[i].y;
								if (vertex_in_plane[i].z > p_max.z) p_max.z = vertex_in_plane[i].z;
							}

							//right_side
							std::function<Point<double>(bool& is_solve, int id_element, Point<double> X)> sourse = [&Edge_self, &test_triag, &p_min, &p_max, &correction_vector, &solver_grid_faces, &id_face_local](bool& is_solve, int id_element, Point<double> X) -> Point<double>
							{
								is_solve = true;
								Point<double> result;

								Point<double> X_self = test_triag.MakeTransferXyzIntoSelf(X) + correction_vector;
								double hx = abs(p_max.x - p_min.x);
								double hy = abs(p_max.y - p_min.y);
								double lambda, mu;
								solver_grid_faces[id_face_local].GetDomain(solver_grid_faces[id_face_local].GetElement(id_element)->GetIdDomain())->forMech.GetLameCoefficients(lambda, mu);
								auto A = *test_triag.GetSelfReverseBasis();

								double ksi_x = (X_self.x - p_min.x) / hx;
								double ksi_y = (X_self.y - p_min.y) / hy;

								double dphi_dx;
								double dphi_dy;
								double dphi_dxdx;
								double dphi_dydy;
								double dphi_dxdy;
								if (math::IsEqual(Edge_self[0].y, p_min.y))
								{
									dphi_dx = (-4 * ksi_x * ksi_x + 4 * ksi_x) * -1.0 / hy;
									dphi_dy = (1 - ksi_y) * (-8 / hx * ksi_x + 4 / hx);
									dphi_dxdx = (1 - ksi_y) * (-8 / (hx * hx));
									dphi_dydy = 0;
									dphi_dxdy = (-1.0 / hy) * (-8 / hx * ksi_x + 4 / hx);
								}
								else {
									dphi_dx = (-4 * ksi_x * ksi_x + 4 * ksi_x) * 1.0 / hy;
									dphi_dy = (ksi_y) * (-8 / hx * ksi_x + 4 / hx);
									dphi_dxdx = (ksi_y) * (-8 / (hx * hx));
									dphi_dydy = 0;
									dphi_dxdy = (1.0 / hy) * (-8 / hx * ksi_x + 4 / hx);
								}

								double dphi_dXdX = A[0][0] * A[0][0] * dphi_dxdx + 2 * A[0][0] * A[1][0] * dphi_dxdy + A[1][0] * A[1][0] * dphi_dydy;
								double dphi_dYdY = A[0][1] * A[0][1] * dphi_dxdx + 2 * A[0][1] * A[1][1] * dphi_dxdy + A[1][1] * A[1][1] * dphi_dydy;
								double dphi_dZdZ = A[0][2] * A[0][2] * dphi_dxdx + 2 * A[0][2] * A[1][2] * dphi_dxdy + A[1][2] * A[1][2] * dphi_dydy;
								double dphi_dXdY = A[0][0] * A[0][1] * dphi_dxdx + (A[0][0] * A[1][1] + A[0][1] * A[1][0]) * dphi_dxdy + A[1][0] * A[1][1] * dphi_dydy;
								double dphi_dXdZ = A[0][0] * A[0][2] * dphi_dxdx + (A[0][0] * A[1][2] + A[0][2] * A[1][0]) * dphi_dxdy + A[1][0] * A[1][2] * dphi_dydy;
								double dphi_dYdZ = A[0][1] * A[0][2] * dphi_dxdx + (A[0][1] * A[1][2] + A[0][2] * A[1][1]) * dphi_dxdy + A[1][1] * A[1][2] * dphi_dydy;

								double laplasian_phi = dphi_dXdX + dphi_dYdY + dphi_dZdZ;
								Point<double> DD_phi;
								DD_phi.x = dphi_dXdX + dphi_dXdY + dphi_dXdZ;
								DD_phi.y = dphi_dXdY + dphi_dYdY + dphi_dYdZ;
								DD_phi.z = dphi_dXdZ + dphi_dYdZ + dphi_dZdZ;

								result = DD_phi * (lambda + mu) + Point<double>(laplasian_phi, laplasian_phi, laplasian_phi) * mu;

								return result * -1.;
							};
							std::function<Point<double>(bool& is_solve, int id_element, Point<double> X)> sourse_radial = [&Edge_self, &test_triag, p_min, p_max, correction_vector, &solver_grid_faces, id_face_local](bool& is_solve, int id_element, Point<double> X) -> Point<double>
							{
								is_solve = true;
								Point<double> result;

								Point<double> X_self = test_triag.MakeTransferXyzIntoSelf(X);
								double hx = abs(p_max.x - p_min.x);
								double hy = abs(p_max.y - p_min.y);
								double lambda, mu;
								solver_grid_faces[id_face_local].GetDomain(solver_grid_faces[id_face_local].GetElement(id_element)->GetIdDomain())->forMech.GetLameCoefficients(lambda, mu);
								auto A = *test_triag.GetSelfReverseBasis();

								Point<double> Edge_centr = (Edge_self[0] + Edge_self[1]) / 2.;
								double r = math::SolveLengthVector(X_self, Edge_centr);
								double hr = math::SolveLengthVector(Edge_self[0], Edge_centr);
								if (r > hr) return result;

								double ksi_r = r / hr;
								double phi = -1 * ksi_r * ksi_r + 1;
								double dr_dx = (X_self.x - Edge_centr.x) / hr;
								double dr_dy = (X_self.y - Edge_centr.y) / hr;
								double dr_dz = (X_self.z - Edge_centr.z) / hr;
								double dr_dxdx = ((X_self.y - Edge_centr.y) * (X_self.y - Edge_centr.y) + (X_self.z - Edge_centr.z) * (X_self.z - Edge_centr.z)) / (hr * hr * hr);
								double dr_dydy = ((X_self.x - Edge_centr.x) * (X_self.x - Edge_centr.x) + (X_self.z - Edge_centr.z) * (X_self.z - Edge_centr.z)) / (hr * hr * hr);
								double dr_dzdz = ((X_self.x - Edge_centr.x) * (X_self.x - Edge_centr.x) + (X_self.y - Edge_centr.y) * (X_self.y - Edge_centr.y)) / (hr * hr * hr);
								double dr_dxdy = -(X_self.x - Edge_centr.x) * (X_self.y - Edge_centr.y) / (hr * hr * hr);
								double dr_dxdz = -(X_self.x - Edge_centr.x) * (X_self.z - Edge_centr.z) / (hr * hr * hr);
								double dr_dydz = -(X_self.y - Edge_centr.y) * (X_self.z - Edge_centr.z) / (hr * hr * hr);
								double dphi_dr = -2 * r / (hr * hr);
								double dphi_drr = -2 / (hr * hr);
								

								double dphi_dx = dphi_dr * dr_dx;
								double dphi_dy = dphi_dr * dr_dy;
								double dphi_dz = dphi_dr * dr_dz;
								double dphi_dxdx = dphi_drr * dr_dx * dr_dx + dphi_dr * dr_dxdx;
								double dphi_dydy = dphi_drr * dr_dy * dr_dy + dphi_dr * dr_dydy;
								double dphi_dzdz = dphi_drr * dr_dz * dr_dz + dphi_dr * dr_dzdz;
								double dphi_dxdy = dphi_drr * dr_dx * dr_dy + dphi_dr * dr_dxdy;
								double dphi_dxdz = dphi_drr * dr_dx * dr_dz + dphi_dr * dr_dxdz;
								double dphi_dydz = dphi_drr * dr_dy * dr_dz + dphi_dr * dr_dydz;
								

								double dphi_dXdX = A[0][0] * A[0][0] * dphi_dxdx
									+ A[1][0] * A[1][0] * dphi_dydy
									+ A[2][0] * A[2][0] * dphi_dzdz
									+ 2 * A[0][0] * A[1][0] * dphi_dxdy
									+ 2 * A[0][0] * A[2][0] * dphi_dxdz
									+ 2 * A[1][0] * A[2][0] * dphi_dydz;
								double dphi_dYdY = A[0][1] * A[0][1] * dphi_dxdx
									+ A[1][1] * A[1][1] * dphi_dydy
									+ A[2][1] * A[2][1] * dphi_dzdz
									+ 2 * A[0][1] * A[1][1] * dphi_dxdy
									+ 2 * A[0][1] * A[2][1] * dphi_dxdz
									+ 2 * A[1][1] * A[2][1] * dphi_dydz;
								double dphi_dZdZ = A[0][2] * A[0][2] * dphi_dxdx
									+ A[1][2] * A[1][2] * dphi_dydy
									+ A[2][2] * A[2][2] * dphi_dzdz
									+ 2 * A[0][2] * A[1][2] * dphi_dxdy
									+ 2 * A[0][2] * A[2][2] * dphi_dxdz
									+ 2 * A[1][2] * A[2][2] * dphi_dydz;
								double dphi_dXdY = A[0][0] * A[0][1] * dphi_dxdx
									+ A[1][0] * A[1][1] * dphi_dydy
									+ A[2][0] * A[2][1] * dphi_dzdz
									+ (A[0][0] * A[1][1] + A[0][1] * A[1][0]) * dphi_dxdy
									+ (A[0][0] * A[2][1] + A[0][1] * A[2][0]) * dphi_dxdz
									+ (A[1][0] * A[2][1] + A[1][1] * A[2][0]) * dphi_dydz;
								double dphi_dXdZ = A[0][0] * A[0][2] * dphi_dxdx
									+ A[1][0] * A[1][2] * dphi_dydy
									+ A[2][0] * A[2][2] * dphi_dzdz
									+ (A[0][0] * A[1][2] + A[0][2] * A[1][0]) * dphi_dxdy
									+ (A[0][0] * A[2][2] + A[0][2] * A[2][0]) * dphi_dxdz
									+ (A[1][0] * A[2][2] + A[1][2] * A[2][0]) * dphi_dydz;
								double dphi_dYdZ = A[0][1] * A[0][2] * dphi_dxdx
									+ A[1][1] * A[1][2] * dphi_dydy
									+ A[2][1] * A[2][2] * dphi_dzdz
									+ (A[0][1] * A[1][2] + A[0][2] * A[1][1]) * dphi_dxdy
									+ (A[0][1] * A[2][2] + A[0][2] * A[2][1]) * dphi_dxdz
									+ (A[1][1] * A[2][2] + A[1][2] * A[2][1]) * dphi_dydz;

								double laplasian_phi = dphi_dXdX + dphi_dYdY + dphi_dZdZ;
								Point<double> DD_phi;
								DD_phi.x = dphi_dXdX + dphi_dXdY + dphi_dXdZ;
								DD_phi.y = dphi_dXdY + dphi_dYdY + dphi_dYdZ;
								DD_phi.z = dphi_dXdZ + dphi_dYdZ + dphi_dZdZ;

								result = DD_phi * (lambda + mu) + Point<double>(laplasian_phi, laplasian_phi, laplasian_phi) * mu;

								return result * -1.;
							};

							FEM::FEM_2D_forElasticDeformation(
								false,
								1e-15,
								faces_grid_simple[id_face_local], //input
								first_boundaries, //input
								second_boundary_empty, //input
								sourse,
								element->self_direction, //output
								solver_grid_faces[id_face_local], //output
								Solution_faces[id_face_local][id_bf + faces[id_face].GetLowerElementCount()], //output
								Base_stiffness_matrix //output
							);

							/*for (int i = 0; i < solver_grid_faces[id_face_local].GetElementsCount(); i++)
							{
								for (int j = 0; j < solver_grid_faces[id_face_local].GetElement(i)->GetDOFsCount(); j++)
								{
									Point<double> X = solver_grid_faces[id_face_local].GetElement(i)->GetNode(j);
									Point<double> X_self = test_triag.MakeTransferXyzIntoSelf(X) + correction_vector;

									double hx = abs(p_max.x - p_min.x);
									double hy = abs(p_max.y - p_min.y);
									double ksi_x = (X_self.x - p_min.x) / hx;
									double ksi_y = (X_self.y - p_min.y) / hy;

									double phi;
									if (math::IsEqual(Edge_self[0].y, p_min.y))
									{
										phi = (-4 * ksi_x * ksi_x + 4 * ksi_x) * (1 - ksi_y);
									}
									else {
										phi = (-4 * ksi_x * ksi_x + 4 * ksi_x) * (ksi_y);
									}


									Solution_faces[id_face_local][id_bf + faces[id_face].GetLowerElementCount()][solver_grid_faces[id_face_local].GetElement(i)->GetDOFInLocalID(j)].x = phi;
								}
							}*/

							if(this->is_print_logFile)
							{
								printf_s("Print the mech result into .dat file... ");

								FILE* fout_tech;
								char name_u_tech[5000];
								sprintf_s(name_u_tech, "%s/Face_%d(global%d)_bubble_by_edge_bf%d.dat", element->self_direction, id_face_local, id_face, id_bf + faces[id_face].GetLowerElementCount());
								fopen_s(&fout_tech, name_u_tech, "w");
								char name_in_file[1000];
								sprintf_s(name_in_file, "Face_%d(global%d)_bubble_by_edge_bf%d", id_face_local, id_face, id_bf + faces[id_face].GetLowerElementCount());
								std::vector<std::vector<char>> name_value(4);
								char name_v_tmp[6][100];
								sprintf_s(name_v_tmp[0], "Ux");
								sprintf_s(name_v_tmp[1], "Uy");
								sprintf_s(name_v_tmp[2], "Uz");
								sprintf_s(name_v_tmp[3], "Material");
								for (int i = 0; i < name_value.size(); i++)
								{
									name_value[i].resize(100);
									for (int j = 0; j < name_value[i].size(); j++)
									{
										name_value[i][j] = name_v_tmp[i][j];
									}
								}
								std::vector<std::vector<double>> value(4);
								value[0].resize(solver_grid_faces[id_face_local].GetVertexCount());
								value[1].resize(solver_grid_faces[id_face_local].GetVertexCount());
								value[2].resize(solver_grid_faces[id_face_local].GetVertexCount());
								value[3].resize(solver_grid_faces[id_face_local].GetElementsCount());

								for (int i = 0; i < solver_grid_faces[id_face_local].GetVertexCount(); i++)
								{
									value[0][i] = Solution_faces[id_face_local][id_bf + faces[id_face].GetLowerElementCount()][i].x;
									value[1][i] = Solution_faces[id_face_local][id_bf + faces[id_face].GetLowerElementCount()][i].y;
									value[2][i] = Solution_faces[id_face_local][id_bf + faces[id_face].GetLowerElementCount()][i].z;
								}
								for (int i = 0; i < solver_grid_faces[id_face_local].GetElementsCount(); i++)
								{
									auto element = solver_grid_faces[id_face_local].GetElement(i);
									value[3][i] = element->GetIdDomain();

								}
								solver_grid_faces[id_face_local].printTecPlot3D_DiffDomains(fout_tech, value, name_value, name_in_file);
								fclose(fout_tech);
							}
						}

						//bubble function
						printf_s("create bubble function for face %d(local: %d) - bf bubble\n", id_face, id_face_local);
						if (true)
						{
							for (int i = 0; i < first_boundaries.size(); i++)
							{
								first_boundaries[i].value = [](Point<bool>& arg, int id) -> Point<double>
								{
									arg.x = true;
									arg.y = true;
									arg.z = true;
									return Point<double>(0, 0, 0);
								};
							}
							struct _Neumann {
								/*double value;
								Point<double> vector;*/
								std::function<Point<double>(Point<double>)> value;
								std::vector<std::function<Point<double>(Point<double>)>> values;
								std::vector<std::vector<int>> id_vertexes_as_triangle;
								std::vector<int> id_base_element;
							};
							std::vector<_Neumann> second_boundary_empty;

							int target_edge_id_global = this->faces[id_face].GetLowerElement(0)->GetSelfGlobalId();
							int target_edge_id = math::GetPositionInSortVector(id_self_edges, target_edge_id_global);

							Point<double> BaseEdge[2];
							/*BaseEdge[0] = this->GetCoordinateViaID(boundary_vertexes[target_edge_id][0]);
							BaseEdge[1] = this->GetCoordinateViaID(boundary_vertexes[target_edge_id][1]);*/
							BaseEdge[0] = this->edges[target_edge_id_global].GetNode(0);
							BaseEdge[1] = this->edges[target_edge_id_global].GetNode(1);

							//делаем проекцию всех узлов на плоскость (BaseEdge[0], BaseEdge[1], centr_of_polygon)
							Point<double> plane[3];
							plane[0] = BaseEdge[0];
							plane[1] = BaseEdge[1];
							plane[2] = centr_of_polygon;
							std::vector<Point<double>> vertex_in_plane(this->faces[id_face].GetNodesCount());
							DenseMatrix<double, double> _matrix;
							for (int i = 0; i < vertex_in_plane.size(); i++)
								//vertex_in_plane[i] = math::MakeProgectionOfPointIntoPlane(plane, this->faces[id_face].GetNode(i), _matrix);
								vertex_in_plane[i] = this->faces[id_face].GetNode(i);
							//разворачиваем базис
							geometry::Triangle test_triag;
							test_triag.SetNode(0, &plane[0]);
							test_triag.SetNode(1, &plane[1]);
							test_triag.SetNode(2, &plane[2]);
							test_triag.CreateBasis();
							//строим прямоугольник
							Point<double> p_min = test_triag.MakeTransferXyzIntoSelf(vertex_in_plane[0]);
							Point<double> p_max = test_triag.MakeTransferXyzIntoSelf(vertex_in_plane[0]);
							for (int i = 0; i < vertex_in_plane.size(); i++)
							{
								Point<double> p_in_self = test_triag.MakeTransferXyzIntoSelf(vertex_in_plane[i]);
								if (p_in_self.x < p_min.x) p_min.x = p_in_self.x;
								if (p_in_self.y < p_min.y) p_min.y = p_in_self.y;
								if (p_in_self.z < p_min.z) p_min.z = p_in_self.z;
								if (p_in_self.x > p_max.x) p_max.x = p_in_self.x;
								if (p_in_self.y > p_max.y) p_max.y = p_in_self.y;
								if (p_in_self.z > p_max.z) p_max.z = p_in_self.z;
							}
							//нельзя переводить прямоугольник в базовые координаты. придется всегда работать через проекцию?
							//попробуем таки по новым
							/*Point<double> p_min_xyz = test_triag.MakeTransferSelfIntoXyz(p_min);
							Point<double> p_max_xyz = test_triag.MakeTransferSelfIntoXyz(p_max);*/

							//right_side
							std::function<Point<double>(bool& is_solve, int id_element, Point<double> X)> sourse = [&test_triag, &plane, &p_min, &p_max, &solver_grid_faces, &id_face_local](bool& is_solve, int id_element, Point<double> X) -> Point<double>
							{
								is_solve = true;
								Point<double> result;
								double f;

								Point<double> X_in_plane = test_triag.MakeTransferXyzIntoSelf(X);
								auto A = *test_triag.GetSelfReverseBasis();

								double hx = abs(p_max.x - p_min.x);
								double hy = abs(p_max.y - p_min.y);
								double ksi_x = (X_in_plane.x - p_min.x) / hx;
								double ksi_y = (X_in_plane.y - p_min.y) / hy;
								double dphi_dx = (-8 * ksi_x / hx + 4 / hx) * (-4 * ksi_y * ksi_y + 4 * ksi_y);
								double dphi_dy = (-4 * ksi_x * ksi_x + 4 * ksi_x) * (-8 / hy * ksi_y + 4 / hy);
								double dphi_dxdx = (-8 * ksi_x / (hx * hx)) * (-4 * ksi_y * ksi_y + 4 * ksi_y);
								double dphi_dydy = (-4 * ksi_x * ksi_x + 4 * ksi_x) * (-8 / (hy * hy));
								double dphi_dxdy = (-8 * ksi_x / hx + 4 / hx) * (-8 / hy * ksi_y + 4 / hy);

								double dphi_dXdX = A[0][0] * A[0][0] * dphi_dxdx + 2 * A[0][0] * A[1][0] * dphi_dxdy + A[1][0] * A[1][0] * dphi_dydy;
								double dphi_dYdY = A[0][1] * A[0][1] * dphi_dxdx + 2 * A[0][1] * A[1][1] * dphi_dxdy + A[1][1] * A[1][1] * dphi_dydy;
								double dphi_dZdZ = A[0][2] * A[0][2] * dphi_dxdx + 2 * A[0][2] * A[1][2] * dphi_dxdy + A[1][2] * A[1][2] * dphi_dydy;
								double dphi_dXdY = A[0][0] * A[0][1] * dphi_dxdx + (A[0][0] * A[1][1] + A[0][1] * A[1][0]) * dphi_dxdy + A[1][0] * A[1][1] * dphi_dydy;
								double dphi_dXdZ = A[0][0] * A[0][2] * dphi_dxdx + (A[0][0] * A[1][2] + A[0][2] * A[1][0]) * dphi_dxdy + A[1][0] * A[1][2] * dphi_dydy;
								double dphi_dYdZ = A[0][1] * A[0][2] * dphi_dxdx + (A[0][1] * A[1][2] + A[0][2] * A[1][1]) * dphi_dxdy + A[1][1] * A[1][2] * dphi_dydy;

								double laplasian_phi = dphi_dXdX + dphi_dYdY + dphi_dZdZ;
								Point<double> DD_phi;
								DD_phi.x = dphi_dXdX + dphi_dXdY + dphi_dXdZ;
								DD_phi.y = dphi_dXdY + dphi_dYdY + dphi_dYdZ;
								DD_phi.z = dphi_dXdZ + dphi_dYdZ + dphi_dZdZ;

								double lambda, mu;
								solver_grid_faces[id_face_local].GetDomain(solver_grid_faces[id_face_local].GetElement(id_element)->GetIdDomain())->forMech.GetLameCoefficients(lambda, mu);
								result = DD_phi * (lambda + mu) + Point<double>(laplasian_phi, laplasian_phi, laplasian_phi) * mu;

								return result * -1.;
							};
							std::function<Point<double>(bool& is_solve, int id_element, Point<double> X)> sourse_radial = [centr_of_polygon, &vertex_in_plane, &solver_grid_faces, id_face_local](bool& is_solve, int id_element, Point<double> X) -> Point<double>
							{
								is_solve = true;
								Point<double> result;


								double lambda, mu;
								solver_grid_faces[id_face_local].GetDomain(solver_grid_faces[id_face_local].GetElement(id_element)->GetIdDomain())->forMech.GetLameCoefficients(lambda, mu);
								
								double min_len = math::SolveLengthVector(vertex_in_plane[0], centr_of_polygon);
								for (int i = 1; i < vertex_in_plane.size(); i++)
								{
									double len  = math::SolveLengthVector(vertex_in_plane[i], centr_of_polygon);
									if (len < min_len)
									{
										min_len = len;
									}
								}
								for (int i = 0; i < vertex_in_plane.size(); i++)
								{
									int next = (i + 1) >= vertex_in_plane.size() ? 0 : (i + 1);
									double len = math::SolveLengthPointAndLine(vertex_in_plane[i], vertex_in_plane[next], centr_of_polygon);
									if (len < min_len)
									{
										min_len = len;
									}
								}

								double r = math::SolveLengthVector(X, centr_of_polygon);
								double hr = min_len;
								if (r > hr) return result;

								double ksi_r = r / hr;
								double phi = -1 * ksi_r * ksi_r + 1;
								double dr_dx = (X.x - centr_of_polygon.x) / hr;
								double dr_dy = (X.y - centr_of_polygon.y) / hr;
								double dr_dz = (X.z - centr_of_polygon.z) / hr;
								double dr_dxdx = ((X.y - centr_of_polygon.y) * (X.y - centr_of_polygon.y) + (X.z - centr_of_polygon.z) * (X.z - centr_of_polygon.z)) / (hr * hr * hr);
								double dr_dydy = ((X.x - centr_of_polygon.x) * (X.x - centr_of_polygon.x) + (X.z - centr_of_polygon.z) * (X.z - centr_of_polygon.z)) / (hr * hr * hr);
								double dr_dzdz = ((X.x - centr_of_polygon.x) * (X.x - centr_of_polygon.x) + (X.y - centr_of_polygon.y) * (X.y - centr_of_polygon.y)) / (hr * hr * hr);
								double dr_dxdy = -(X.x - centr_of_polygon.x) * (X.y - centr_of_polygon.y) / (hr * hr * hr);
								double dr_dxdz = -(X.x - centr_of_polygon.x) * (X.z - centr_of_polygon.z) / (hr * hr * hr);
								double dr_dydz = -(X.y - centr_of_polygon.y) * (X.z - centr_of_polygon.z) / (hr * hr * hr);
								double dphi_dr = -2 * r / (hr * hr);
								double dphi_drr = -2 / (hr * hr);


								double dphi_dx = dphi_dr * dr_dx;
								double dphi_dy = dphi_dr * dr_dy;
								double dphi_dz = dphi_dr * dr_dz;
								double dphi_dxdx = dphi_drr * dr_dx * dr_dx + dphi_dr * dr_dxdx;
								double dphi_dydy = dphi_drr * dr_dy * dr_dy + dphi_dr * dr_dydy;
								double dphi_dzdz = dphi_drr * dr_dz * dr_dz + dphi_dr * dr_dzdz;
								double dphi_dxdy = dphi_drr * dr_dx * dr_dy + dphi_dr * dr_dxdy;
								double dphi_dxdz = dphi_drr * dr_dx * dr_dz + dphi_dr * dr_dxdz;
								double dphi_dydz = dphi_drr * dr_dy * dr_dz + dphi_dr * dr_dydz;


								double laplasian_phi = dphi_dxdx + dphi_dydy + dphi_dzdz;
								Point<double> DD_phi;
								DD_phi.x = dphi_dxdx + dphi_dxdy + dphi_dxdz;
								DD_phi.y = dphi_dxdy + dphi_dydy + dphi_dydz;
								DD_phi.z = dphi_dxdz + dphi_dydz + dphi_dzdz;

								result = DD_phi * (lambda + mu) + Point<double>(laplasian_phi, laplasian_phi, laplasian_phi) * mu;

								return result * -1.;
							};



							FEM::FEM_2D_forElasticDeformation(
								false,
								1e-15,
								faces_grid_simple[id_face_local], //input
								first_boundaries, //input
								second_boundary_empty, //input
								sourse,
								element->self_direction, //output
								solver_grid_faces[id_face_local], //output
								Solution_faces[id_face_local][0 + 2 * faces[id_face].GetLowerElementCount()], //output
								Base_stiffness_matrix //output
							);

							if(this->is_print_logFile){
								printf_s("Print the mech result into .dat file... ");

								FILE* fout_tech;
								char name_u_tech[5000];
								sprintf_s(name_u_tech, "%s/Face_%d(global%d)_bubble_bf%d.dat", element->self_direction, id_face_local, id_face, 0 + 2 * faces[id_face].GetLowerElementCount());
								fopen_s(&fout_tech, name_u_tech, "w");
								char name_in_file[1000];
								sprintf_s(name_in_file, "Face_%d(global%d)_bubble_bf%d", id_face_local, id_face, 0 + 2 * faces[id_face].GetLowerElementCount());
								std::vector<std::vector<char>> name_value(4);
								char name_v_tmp[6][100];
								sprintf_s(name_v_tmp[0], "Ux");
								sprintf_s(name_v_tmp[1], "Uy");
								sprintf_s(name_v_tmp[2], "Uz");
								sprintf_s(name_v_tmp[3], "Material");
								for (int i = 0; i < name_value.size(); i++)
								{
									name_value[i].resize(100);
									for (int j = 0; j < name_value[i].size(); j++)
									{
										name_value[i][j] = name_v_tmp[i][j];
									}
								}
								std::vector<std::vector<double>> value(4);
								value[0].resize(solver_grid_faces[id_face_local].GetVertexCount());
								value[1].resize(solver_grid_faces[id_face_local].GetVertexCount());
								value[2].resize(solver_grid_faces[id_face_local].GetVertexCount());
								value[3].resize(solver_grid_faces[id_face_local].GetElementsCount());

								for (int i = 0; i < solver_grid_faces[id_face_local].GetVertexCount(); i++)
								{
									value[0][i] = Solution_faces[id_face_local][0 + 2 * faces[id_face].GetLowerElementCount()][i].x;
									value[1][i] = Solution_faces[id_face_local][0 + 2 * faces[id_face].GetLowerElementCount()][i].y;
									value[2][i] = Solution_faces[id_face_local][0 + 2 * faces[id_face].GetLowerElementCount()][i].z;
								}
								for (int i = 0; i < solver_grid_faces[id_face_local].GetElementsCount(); i++)
								{
									auto element = solver_grid_faces[id_face_local].GetElement(i);
									value[3][i] = element->GetIdDomain();

								}
								solver_grid_faces[id_face_local].printTecPlot3D_DiffDomains(fout_tech, value, name_value, name_in_file);
								fclose(fout_tech);
							}
						}
					}
				
					//solve volume problems
					if (true)
					{
						CSSD_Matrix<Tensor2Rank3D, Point<double>> Base_stiffness_matrix;
						Base_stiffness_matrix.print_logs = this->is_print_logFile;

						for (int i = 0; i < this->GetDomainsCount(); i++)
							element->self_grid.AddDomain(*this->GetDomain(i));

						//create Boundary values and vertexes
						printf_s("create Boundary values and vertexes for volume\n");
						struct _Dirichlet {
							int global_id, in_elem_id;
							std::vector<int> id_vertexes;
							std::vector<Point<double>> xyz;
							std::function<Point<double>(Point<bool>&, int id_vertex)> value;
						};
						std::vector<_Dirichlet> first_boundaries(element->GetLowerElementCount());
						for (int id_face_local = 0; id_face_local < first_boundaries.size(); id_face_local++)
						{
							first_boundaries[id_face_local].global_id = element->GetLowerElement(id_face_local)->GetSelfGlobalId();
							first_boundaries[id_face_local].in_elem_id = id_face_local;

							for (int i = 0; i < faces_grid_simple[id_face_local].xyz.size(); i++)
							{
								int global_id_vertex = faces_grid_simple[id_face_local].TransferIdSelfIntoGlobal(i);
								first_boundaries[id_face_local].id_vertexes.push_back(global_id_vertex);
								first_boundaries[id_face_local].xyz.push_back(base_grid_simple.xyz[global_id_vertex]);
							}
						}
						struct _Neumann {
							/*double value;
							Point<double> vector;*/
							std::function<Point<double>(Point<double>)> value;
							std::vector<std::function<Point<double>(Point<double>)>> values;
							std::vector<std::vector<int>> id_vertexes_as_triangle;
							std::vector<int> id_base_element;
						};
						std::vector<_Neumann> second_boundary_empty;

						//linear problems
						for (int id_bf = 0; id_bf < element->GetNodesCount(); id_bf++)
						{
							printf_s("create linear function for volume %d - bf linear %d\n", id_element, id_bf);

							int target_vertex_id_global = element->GetIdNode(id_bf);

							for (int i = 0; i < first_boundaries.size(); i++)
							{
								first_boundaries[i].value = [](Point<bool>& arg, int id) -> Point<double>
								{
									arg.x = true;
									arg.y = true;
									arg.z = true;
									return Point<double>(0, 0, 0);
								};
							}

							//find the target faces and create boundary conditions
							for (int id_face_local = 0; id_face_local < element->GetLowerElementCount(); id_face_local++)
							{
								for (int i = 0; i < element->GetLowerElement(id_face_local)->GetLowerElementCount(); i++)
								{
									if (target_id_for_BF[id_face_local][i] == target_vertex_id_global)
									{
										target_vertex_id_global *= 1;
										first_boundaries[id_face_local].value = [&Solution_faces, id_face_local, i, &faces_grid_simple](Point<bool>& arg, int id_vertex)->Point <double>
										{
											arg.x = true; arg.y = true; arg.z = true;
											Point<double> result = Solution_faces[id_face_local][i][faces_grid_simple[id_face_local].TransferIdGlobalIntoSelf(id_vertex)];

											return result;
										};
										break;
									}
								}
							}

							//right_side
							std::function<Point<double>(bool& is_solve, int id_element, Point<double> X)> sourse = [](bool& is_solve, int id_element, Point<double> X) -> Point<double>
							{
								is_solve = false;
								return Point<double>(0, 0, 0);
							};

							FEM::FEM_3D_forElasticDeformation(
								false,
								1.e-15,
								base_grid_simple, //input
								first_boundaries, //input
								second_boundary_empty, //input
								sourse,
								element->self_direction, //output
								element->self_grid, //output
								element->self_basis_functions[id_bf], //output
								Base_stiffness_matrix //output
							);

							if(this->is_print_logFile){
								printf_s("Print the mech result into .dat file... ");

								FILE* fout_tech;
								char name_u_tech[1000];
								sprintf_s(name_u_tech, "%s/Volume_linear_bf%d.dat", element->self_direction, id_bf);
								fopen_s(&fout_tech, name_u_tech, "w");
								char name_in_file[1000];
								sprintf_s(name_in_file, "Volume_linear_bf%d", id_bf);
								std::vector<std::vector<char>> name_value(4);
								char name_v_tmp[6][100];
								sprintf_s(name_v_tmp[0], "Ux");
								sprintf_s(name_v_tmp[1], "Uy");
								sprintf_s(name_v_tmp[2], "Uz");
								sprintf_s(name_v_tmp[3], "Material");
								for (int i = 0; i < name_value.size(); i++)
								{
									name_value[i].resize(100);
									for (int j = 0; j < name_value[i].size(); j++)
									{
										name_value[i][j] = name_v_tmp[i][j];
									}
								}
								std::vector<std::vector<double>> value(4);
								value[0].resize(element->self_grid.GetVertexCount());
								value[1].resize(element->self_grid.GetVertexCount());
								value[2].resize(element->self_grid.GetVertexCount());
								value[3].resize(element->self_grid.GetElementsCount());

								for (int i = 0; i < element->self_grid.GetVertexCount(); i++)
								{
									value[0][i] = element->self_basis_functions[id_bf][i].x;
									value[1][i] = element->self_basis_functions[id_bf][i].y;
									value[2][i] = element->self_basis_functions[id_bf][i].z;
								}
								for (int i = 0; i < element->self_grid.GetElementsCount(); i++)
								{
									auto el = element->self_grid.GetElement(i);
									value[3][i] = el->GetIdDomain();

								}
								element->self_grid.printTecPlot3D_DiffDomains(fout_tech, value, name_value, name_in_file);
								fclose(fout_tech);
							}
						}

						//bubble problems by edges
						for (int id_bf = 0; id_bf < id_self_edges.size(); id_bf++)
						{
							printf_s("create edges bubble function for volume %d - bf edge_bubble %d\n", id_element, id_bf);

							int target_edges_id_global = id_self_edges[id_bf];

							for (int i = 0; i < first_boundaries.size(); i++)
							{
								first_boundaries[i].value = [](Point<bool>& arg, int id) -> Point<double>
								{
									arg.x = true;
									arg.y = true;
									arg.z = true;
									return Point<double>(0, 0, 0);
								};
							}

							//find the target faces and create boundary conditions
							int id_base_face_global;
							for (int id_face_local = 0; id_face_local < element->GetLowerElementCount(); id_face_local++)
							{
								for (int i = 0; i < element->GetLowerElement(id_face_local)->GetLowerElementCount(); i++)
								{
									if (target_id_for_BF[id_face_local][i + element->GetLowerElement(id_face_local)->GetLowerElementCount()] == target_edges_id_global)
									{
										id_base_face_global = element->GetLowerElement(id_face_local)->GetSelfGlobalId();
										int id_solution_in_face = i + faces[id_base_face_global].GetLowerElementCount();
										first_boundaries[id_face_local].value = [&Solution_faces, id_face_local, id_solution_in_face, &faces_grid_simple](Point<bool>& arg, int id_vertex)->Point <double>
										{
											arg.x = true; arg.y = true; arg.z = true;

											Point<double> result = Solution_faces[id_face_local][id_solution_in_face][faces_grid_simple[id_face_local].TransferIdGlobalIntoSelf(id_vertex)];

											return result;
										};
									}
								}
							}

							Point<double> base_face_centr = this->faces[id_base_face_global].GetWeightCentr();
							Point<double> target_edge[2];
							target_edge[0] = this->edges[target_edges_id_global].GetNode(0);
							target_edge[1] = this->edges[target_edges_id_global].GetNode(1);
							geometry::Triangle test_triang;
							test_triang.SetNode(0, &target_edge[0]);
							test_triang.SetNode(1, &target_edge[1]);
							test_triang.SetNode(2, &base_face_centr);
							test_triang.CreateBasis();
							std::vector<Point<double>> new_element_nodes(element->GetNodesCount());
							for (int i = 0; i < new_element_nodes.size(); i++)
							{
								new_element_nodes[i] = test_triang.MakeTransferXyzIntoSelf(element->GetNode(i));
							}
							Point<double> new_down_parallelepiped = new_element_nodes[0];
							Point<double> new_up_parallelepiped = new_element_nodes[0];
							for (int i = 1; i < new_element_nodes.size(); i++)
							{
								if (new_down_parallelepiped.x > new_element_nodes[i].x) new_down_parallelepiped.x = new_element_nodes[i].x;
								if (new_down_parallelepiped.y > new_element_nodes[i].y) new_down_parallelepiped.y = new_element_nodes[i].y;
								if (new_down_parallelepiped.z > new_element_nodes[i].z) new_down_parallelepiped.z = new_element_nodes[i].z;
								if (new_up_parallelepiped.x < new_element_nodes[i].x) new_up_parallelepiped.x = new_element_nodes[i].x;
								if (new_up_parallelepiped.y < new_element_nodes[i].y) new_up_parallelepiped.y = new_element_nodes[i].y;
								if (new_up_parallelepiped.z < new_element_nodes[i].z) new_up_parallelepiped.z = new_element_nodes[i].z;
							}
							Point<double> size_parallelepiped = new_up_parallelepiped - new_down_parallelepiped;


							//right_side
							std::function<Point<double>(bool& is_solve, int id_element, Point<double> X)> sourse = [&target_edge, &size_parallelepiped, &test_triang, &new_down_parallelepiped, &new_up_parallelepiped, &element](bool& is_solve, int id_element, Point<double> X) -> Point<double>
							{
								is_solve = true;
								Point<double> result;
								/*Point<double> new_X = test_triang.MakeTransferXyzIntoSelf(X);
								result.x = -8.0 / pow(size_parallelepiped.x, 2) * ((new_up_parallelepiped.y - new_X.y) / size_parallelepiped.y) * ((new_up_parallelepiped.z - new_X.z) / size_parallelepiped.z);
								result.y = result.x;
								result.z = result.x;*/


								double hx = size_parallelepiped.x;
								double hy = size_parallelepiped.y;
								double hz = size_parallelepiped.z;
								double lambda, mu;
								element->self_grid.GetDomain(element->self_grid.GetElement(id_element)->GetIdDomain())->forMech.GetLameCoefficients(lambda, mu);
								auto A = *test_triang.GetSelfReverseBasis();

								Point<double> X_self = test_triang.MakeTransferXyzIntoSelf(X);
								Point<double> Edge_self[2];
								Edge_self[0] = test_triang.MakeTransferXyzIntoSelf(target_edge[0]);
								Edge_self[1] = test_triang.MakeTransferXyzIntoSelf(target_edge[1]);

								double ksi_x = (X_self.x - new_down_parallelepiped.x) / hx;
								double ksi_y = (X_self.y - new_down_parallelepiped.y) / hy;
								double ksi_z = (X_self.z - new_down_parallelepiped.z) / hz;

								double phi_x, phi_y, phi_z;
								double dphi_x, dphi_y, dphi_z;
								double ddphi_x, ddphi_y, ddphi_z;
								double dphi_dxdx;
								double dphi_dydy;
								double dphi_dzdz;
								double dphi_dxdy;
								double dphi_dxdz;
								double dphi_dydz;
								phi_x = -4 * ksi_x * ksi_x + 4 * ksi_x;
								dphi_x = (-8 / hx * ksi_x + 4 / hx);
								ddphi_x = -8 / (hx * hx);
								if (math::IsEqual(Edge_self[0].y, new_down_parallelepiped.y))
								{
									phi_y = 1 - ksi_y;
									dphi_y = -1 / hy;
									ddphi_y = 0;
								}
								else {
									phi_y = ksi_y;
									dphi_y = 1 / hy;
									ddphi_y = 0;
								}
								if (math::IsEqual(Edge_self[0].z, new_down_parallelepiped.z))
								{
									phi_z = 1 - ksi_z;
									dphi_z = -1 / hz;
									ddphi_z = 0;
								}
								else {
									phi_z = ksi_z;
									dphi_z = 1 / hz;
									ddphi_z = 0;
								}

								dphi_dxdx = ddphi_x * phi_y * phi_z;
								dphi_dxdy = dphi_x * dphi_y * phi_z;
								dphi_dxdz = dphi_x * phi_y * dphi_z;
								dphi_dydy = phi_x * ddphi_y * phi_z;
								dphi_dydz = phi_x * dphi_y * dphi_z;
								dphi_dzdz = phi_x * phi_y * ddphi_z;

								double dphi_dXdX = A[0][0] * A[0][0] * dphi_dxdx
									+ A[1][0] * A[1][0] * dphi_dydy
									+ A[2][0] * A[2][0] * dphi_dzdz
									+ 2 * A[0][0] * A[1][0] * dphi_dxdy
									+ 2 * A[0][0] * A[2][0] * dphi_dxdz
									+ 2 * A[1][0] * A[2][0] * dphi_dydz;
								double dphi_dYdY = A[0][1] * A[0][1] * dphi_dxdx
									+ A[1][1] * A[1][1] * dphi_dydy
									+ A[2][1] * A[2][1] * dphi_dzdz
									+ 2 * A[0][1] * A[1][1] * dphi_dxdy
									+ 2 * A[0][1] * A[2][1] * dphi_dxdz
									+ 2 * A[1][1] * A[2][1] * dphi_dydz;
								double dphi_dZdZ = A[0][2] * A[0][2] * dphi_dxdx
									+ A[1][2] * A[1][2] * dphi_dydy
									+ A[2][2] * A[2][2] * dphi_dzdz
									+ 2 * A[0][2] * A[1][2] * dphi_dxdy
									+ 2 * A[0][2] * A[2][2] * dphi_dxdz
									+ 2 * A[1][2] * A[2][2] * dphi_dydz;
								double dphi_dXdY = A[0][0] * A[0][1] * dphi_dxdx
									+ A[1][0] * A[1][1] * dphi_dydy
									+ A[2][0] * A[2][1] * dphi_dzdz
									+ (A[0][0] * A[1][1] + A[0][1] * A[1][0]) * dphi_dxdy
									+ (A[0][0] * A[2][1] + A[0][1] * A[2][0]) * dphi_dxdz
									+ (A[1][0] * A[2][1] + A[1][1] * A[2][0]) * dphi_dydz;
								double dphi_dXdZ = A[0][0] * A[0][2] * dphi_dxdx
									+ A[1][0] * A[1][2] * dphi_dydy
									+ A[2][0] * A[2][2] * dphi_dzdz
									+ (A[0][0] * A[1][2] + A[0][2] * A[1][0]) * dphi_dxdy
									+ (A[0][0] * A[2][2] + A[0][2] * A[2][0]) * dphi_dxdz
									+ (A[1][0] * A[2][2] + A[1][2] * A[2][0]) * dphi_dydz;
								double dphi_dYdZ = A[0][1] * A[0][2] * dphi_dxdx
									+ A[1][1] * A[1][2] * dphi_dydy
									+ A[2][1] * A[2][2] * dphi_dzdz
									+ (A[0][1] * A[1][2] + A[0][2] * A[1][1]) * dphi_dxdy
									+ (A[0][1] * A[2][2] + A[0][2] * A[2][1]) * dphi_dxdz
									+ (A[1][1] * A[2][2] + A[1][2] * A[2][1]) * dphi_dydz;

								double laplasian_phi = dphi_dXdX + dphi_dYdY + dphi_dZdZ;
								Point<double> DD_phi;
								DD_phi.x = dphi_dXdX + dphi_dXdY + dphi_dXdZ;
								DD_phi.y = dphi_dXdY + dphi_dYdY + dphi_dYdZ;
								DD_phi.z = dphi_dXdZ + dphi_dYdZ + dphi_dZdZ;

								result = DD_phi * (lambda + mu) + Point<double>(laplasian_phi, laplasian_phi, laplasian_phi) * mu;

								return result * (-1.0);
							};

							FEM::FEM_3D_forElasticDeformation(
								false,
								1.e-15,
								base_grid_simple, //input
								first_boundaries, //input
								second_boundary_empty, //input
								sourse,
								element->self_direction, //output
								element->self_grid, //output
								element->self_basis_functions[id_bf + element->GetNodesCount()], //output
								Base_stiffness_matrix //output
							);

							if (this->is_print_logFile) {
								printf_s("Print the mech result into .dat file... ");

								FILE* fout_tech;
								char name_u_tech[1000];
								sprintf_s(name_u_tech, "%s/Volume_bubble_by_edge_bf%d.dat", element->self_direction, id_bf + element->GetNodesCount());
								fopen_s(&fout_tech, name_u_tech, "w");
								char name_in_file[1000];
								sprintf_s(name_in_file, "Volume_bubble_by_edge_bf%d", id_bf + element->GetNodesCount());
								std::vector<std::vector<char>> name_value(4);
								char name_v_tmp[6][100];
								sprintf_s(name_v_tmp[0], "Ux");
								sprintf_s(name_v_tmp[1], "Uy");
								sprintf_s(name_v_tmp[2], "Uz");
								sprintf_s(name_v_tmp[3], "Material");
								for (int i = 0; i < name_value.size(); i++)
								{
									name_value[i].resize(100);
									for (int j = 0; j < name_value[i].size(); j++)
									{
										name_value[i][j] = name_v_tmp[i][j];
									}
								}
								std::vector<std::vector<double>> value(4);
								value[0].resize(element->self_grid.GetVertexCount());
								value[1].resize(element->self_grid.GetVertexCount());
								value[2].resize(element->self_grid.GetVertexCount());
								value[3].resize(element->self_grid.GetElementsCount());

								for (int i = 0; i < element->self_grid.GetVertexCount(); i++)
								{
									value[0][i] = element->self_basis_functions[id_bf + element->GetNodesCount()][i].x;
									value[1][i] = element->self_basis_functions[id_bf + element->GetNodesCount()][i].y;
									value[2][i] = element->self_basis_functions[id_bf + element->GetNodesCount()][i].z;
								}
								for (int i = 0; i < element->self_grid.GetElementsCount(); i++)
								{
									auto el = element->self_grid.GetElement(i);
									value[3][i] = el->GetIdDomain();

								}
								element->self_grid.printTecPlot3D_DiffDomains(fout_tech, value, name_value, name_in_file);
								fclose(fout_tech);
							}
						}

						//bubble problems by faces
						for (int id_bf = 0; id_bf < element->GetLowerElementCount(); id_bf++)
						{
							printf_s("create faces bubble function for volume %d - bf face_bubble %d\n", id_element, id_bf);

							int target_faces_id_global = element->GetLowerElement(id_bf)->GetSelfGlobalId();

							for (int i = 0; i < first_boundaries.size(); i++)
							{
								first_boundaries[i].value = [](Point<bool>& arg, int id) -> Point<double>
								{
									arg.x = true;
									arg.y = true;
									arg.z = true;
									return Point<double>(0, 0, 0);
								};
							}

							//find create boundary conditions for target face
							int id_solution_in_face = Solution_faces[id_bf].size() - 1;
							first_boundaries[id_bf].value = [&Solution_faces, id_bf, id_solution_in_face, &faces_grid_simple](Point<bool>& arg, int id_vertex)->Point <double>
							{
								arg.x = true; arg.y = true; arg.z = true;
								Point<double> result = Solution_faces[id_bf][id_solution_in_face][faces_grid_simple[id_bf].TransferIdGlobalIntoSelf(id_vertex)];

								return result;
							};

							Point<double> base_face_centr = this->faces[target_faces_id_global].GetWeightCentr();
							int id_target_edge = this->faces[target_faces_id_global].GetLowerElement(0)->GetSelfGlobalId();
							Point<double> target_edge[2];
							target_edge[0] = this->edges[id_target_edge].GetNode(0);
							target_edge[1] = this->edges[id_target_edge].GetNode(1);
							geometry::Triangle test_triang;
							test_triang.SetNode(0, &target_edge[0]);
							test_triang.SetNode(1, &target_edge[1]);
							test_triang.SetNode(2, &base_face_centr);
							test_triang.CreateBasis();
							std::vector<Point<double>> new_element_nodes(element->GetNodesCount());
							for (int i = 0; i < new_element_nodes.size(); i++)
							{
								new_element_nodes[i] = test_triang.MakeTransferXyzIntoSelf(element->GetNode(i));
							}
							Point<double> new_down_parallelepiped = new_element_nodes[0];
							Point<double> new_up_parallelepiped = new_element_nodes[0];
							for (int i = 1; i < new_element_nodes.size(); i++)
							{
								if (new_down_parallelepiped.x > new_element_nodes[i].x) new_down_parallelepiped.x = new_element_nodes[i].x;
								if (new_down_parallelepiped.y > new_element_nodes[i].y) new_down_parallelepiped.y = new_element_nodes[i].y;
								if (new_down_parallelepiped.z > new_element_nodes[i].z) new_down_parallelepiped.z = new_element_nodes[i].z;
								if (new_up_parallelepiped.x < new_element_nodes[i].x) new_up_parallelepiped.x = new_element_nodes[i].x;
								if (new_up_parallelepiped.y < new_element_nodes[i].y) new_up_parallelepiped.y = new_element_nodes[i].y;
								if (new_up_parallelepiped.z < new_element_nodes[i].z) new_up_parallelepiped.z = new_element_nodes[i].z;
							}
							Point<double> size_parallelepiped = new_up_parallelepiped - new_down_parallelepiped;


							//right_side
							std::function<Point<double>(bool& is_solve, int id_element, Point<double> X)> sourse = [&size_parallelepiped, &new_down_parallelepiped, &element, &target_edge, &test_triang, &new_up_parallelepiped](bool& is_solve, int id_element, Point<double> X) -> Point<double>
							{
								is_solve = true;
								Point<double> result;

								double hx = size_parallelepiped.x;
								double hy = size_parallelepiped.y;
								double hz = size_parallelepiped.z;
								double lambda, mu;
								element->self_grid.GetDomain(element->self_grid.GetElement(id_element)->GetIdDomain())->forMech.GetLameCoefficients(lambda, mu);
								auto A = *test_triang.GetSelfReverseBasis();

								Point<double> X_self = test_triang.MakeTransferXyzIntoSelf(X);
								Point<double> Edge_self[2];
								Edge_self[0] = test_triang.MakeTransferXyzIntoSelf(target_edge[0]);
								Edge_self[1] = test_triang.MakeTransferXyzIntoSelf(target_edge[1]);

								double ksi_x = (X_self.x - new_down_parallelepiped.x) / hx;
								double ksi_y = (X_self.y - new_down_parallelepiped.y) / hy;
								double ksi_z = (X_self.z - new_down_parallelepiped.z) / hz;

								double phi_x, phi_y, phi_z;
								double dphi_x, dphi_y, dphi_z;
								double ddphi_x, ddphi_y, ddphi_z;
								double dphi_dxdx;
								double dphi_dydy;
								double dphi_dzdz;
								double dphi_dxdy;
								double dphi_dxdz;
								double dphi_dydz;
								phi_x = -4 * ksi_x * ksi_x + 4 * ksi_x;
								dphi_x = (-8 / hx * ksi_x + 4 / hx);
								ddphi_x = -8 / (hx * hx);
								phi_y = -4 * ksi_y * ksi_y + 4 * ksi_y;
								dphi_y = (-8 / hy * ksi_y + 4 / hy);
								ddphi_y = -8 / (hy * hy);
								if (math::IsEqual(Edge_self[0].z, new_down_parallelepiped.z))
								{
									phi_z = 1 - ksi_z;
									dphi_z = -1 / hz;
									ddphi_z = 0;
								}
								else {
									phi_z = ksi_z;
									dphi_z = 1 / hz;
									ddphi_z = 0;
								}

								dphi_dxdx = ddphi_x * phi_y * phi_z;
								dphi_dxdy = dphi_x * dphi_y * phi_z;
								dphi_dxdz = dphi_x * phi_y * dphi_z;
								dphi_dydy = phi_x * ddphi_y * phi_z;
								dphi_dydz = phi_x * dphi_y * dphi_z;
								dphi_dzdz = phi_x * phi_y * ddphi_z;

								double dphi_dXdX = A[0][0] * A[0][0] * dphi_dxdx
									+ A[1][0] * A[1][0] * dphi_dydy
									+ A[2][0] * A[2][0] * dphi_dzdz
									+ 2 * A[0][0] * A[1][0] * dphi_dxdy
									+ 2 * A[0][0] * A[2][0] * dphi_dxdz
									+ 2 * A[1][0] * A[2][0] * dphi_dydz;
								double dphi_dYdY = A[0][1] * A[0][1] * dphi_dxdx
									+ A[1][1] * A[1][1] * dphi_dydy
									+ A[2][1] * A[2][1] * dphi_dzdz
									+ 2 * A[0][1] * A[1][1] * dphi_dxdy
									+ 2 * A[0][1] * A[2][1] * dphi_dxdz
									+ 2 * A[1][1] * A[2][1] * dphi_dydz;
								double dphi_dZdZ = A[0][2] * A[0][2] * dphi_dxdx
									+ A[1][2] * A[1][2] * dphi_dydy
									+ A[2][2] * A[2][2] * dphi_dzdz
									+ 2 * A[0][2] * A[1][2] * dphi_dxdy
									+ 2 * A[0][2] * A[2][2] * dphi_dxdz
									+ 2 * A[1][2] * A[2][2] * dphi_dydz;
								double dphi_dXdY = A[0][0] * A[0][1] * dphi_dxdx
									+ A[1][0] * A[1][1] * dphi_dydy
									+ A[2][0] * A[2][1] * dphi_dzdz
									+ (A[0][0] * A[1][1] + A[0][1] * A[1][0]) * dphi_dxdy
									+ (A[0][0] * A[2][1] + A[0][1] * A[2][0]) * dphi_dxdz
									+ (A[1][0] * A[2][1] + A[1][1] * A[2][0]) * dphi_dydz;
								double dphi_dXdZ = A[0][0] * A[0][2] * dphi_dxdx
									+ A[1][0] * A[1][2] * dphi_dydy
									+ A[2][0] * A[2][2] * dphi_dzdz
									+ (A[0][0] * A[1][2] + A[0][2] * A[1][0]) * dphi_dxdy
									+ (A[0][0] * A[2][2] + A[0][2] * A[2][0]) * dphi_dxdz
									+ (A[1][0] * A[2][2] + A[1][2] * A[2][0]) * dphi_dydz;
								double dphi_dYdZ = A[0][1] * A[0][2] * dphi_dxdx
									+ A[1][1] * A[1][2] * dphi_dydy
									+ A[2][1] * A[2][2] * dphi_dzdz
									+ (A[0][1] * A[1][2] + A[0][2] * A[1][1]) * dphi_dxdy
									+ (A[0][1] * A[2][2] + A[0][2] * A[2][1]) * dphi_dxdz
									+ (A[1][1] * A[2][2] + A[1][2] * A[2][1]) * dphi_dydz;

								double laplasian_phi = dphi_dXdX + dphi_dYdY + dphi_dZdZ;
								Point<double> DD_phi;
								DD_phi.x = dphi_dXdX + dphi_dXdY + dphi_dXdZ;
								DD_phi.y = dphi_dXdY + dphi_dYdY + dphi_dYdZ;
								DD_phi.z = dphi_dXdZ + dphi_dYdZ + dphi_dZdZ;

								result = DD_phi * (lambda + mu) + Point<double>(laplasian_phi, laplasian_phi, laplasian_phi) * mu;

								return result * (-1.0);
							};

							FEM::FEM_3D_forElasticDeformation(
								false,
								1.e-15,
								base_grid_simple, //input
								first_boundaries, //input
								second_boundary_empty, //input
								sourse,
								element->self_direction, //output
								element->self_grid, //output
								element->self_basis_functions[id_bf + element->GetNodesCount() + id_self_edges.size()], //output
								Base_stiffness_matrix //output
							);

							if (this->is_print_logFile) {
								printf_s("Print the mech result into .dat file... ");

								FILE* fout_tech;
								char name_u_tech[1000];
								sprintf_s(name_u_tech, "%s/Volume_bubble_by_face_bf%d.dat", element->self_direction, id_bf + element->GetNodesCount() + (int)id_self_edges.size());
								fopen_s(&fout_tech, name_u_tech, "w");
								char name_in_file[1000];
								sprintf_s(name_in_file, "Volume_bubble_by_face_bf%d", id_bf + element->GetNodesCount() + (int)id_self_edges.size());
								std::vector<std::vector<char>> name_value(4);
								char name_v_tmp[6][100];
								sprintf_s(name_v_tmp[0], "Ux");
								sprintf_s(name_v_tmp[1], "Uy");
								sprintf_s(name_v_tmp[2], "Uz");
								sprintf_s(name_v_tmp[3], "Material");
								for (int i = 0; i < name_value.size(); i++)
								{
									name_value[i].resize(100);
									for (int j = 0; j < name_value[i].size(); j++)
									{
										name_value[i][j] = name_v_tmp[i][j];
									}
								}
								std::vector<std::vector<double>> value(4);
								value[0].resize(element->self_grid.GetVertexCount());
								value[1].resize(element->self_grid.GetVertexCount());
								value[2].resize(element->self_grid.GetVertexCount());
								value[3].resize(element->self_grid.GetElementsCount());

								for (int i = 0; i < element->self_grid.GetVertexCount(); i++)
								{
									value[0][i] = element->self_basis_functions[id_bf + element->GetNodesCount() + id_self_edges.size()][i].x;
									value[1][i] = element->self_basis_functions[id_bf + element->GetNodesCount() + id_self_edges.size()][i].y;
									value[2][i] = element->self_basis_functions[id_bf + element->GetNodesCount() + id_self_edges.size()][i].z;
								}
								for (int i = 0; i < element->self_grid.GetElementsCount(); i++)
								{
									auto el = element->self_grid.GetElement(i);
									value[3][i] = el->GetIdDomain();

								}
								element->self_grid.printTecPlot3D_DiffDomains(fout_tech, value, name_value, name_in_file);
								fclose(fout_tech);
							}
						}

						//bubble problems by volume
						if (true) {
							printf_s("create bubble function for volume %d - bf bubble\n", id_element);

							for (int i = 0; i < first_boundaries.size(); i++)
							{
								first_boundaries[i].value = [](Point<bool>& arg, int id) -> Point<double>
								{
									arg.x = true;
									arg.y = true;
									arg.z = true;
									return Point<double>(0, 0, 0);
								};
							}

							Point<double> down_parallelepiped = element->GetNode(0);
							Point<double> up_parallelepiped = element->GetNode(0);
							for (int i = 1; i < element->GetNodesCount(); i++)
							{
								if (down_parallelepiped.x > element->GetNode(i).x) down_parallelepiped.x = element->GetNode(i).x;
								if (down_parallelepiped.y > element->GetNode(i).y) down_parallelepiped.y = element->GetNode(i).y;
								if (down_parallelepiped.z > element->GetNode(i).z) down_parallelepiped.z = element->GetNode(i).z;
								if (up_parallelepiped.x < element->GetNode(i).x) up_parallelepiped.x = element->GetNode(i).x;
								if (up_parallelepiped.y < element->GetNode(i).y) up_parallelepiped.y = element->GetNode(i).y;
								if (up_parallelepiped.z < element->GetNode(i).z) up_parallelepiped.z = element->GetNode(i).z;
							}
							Point<double> size_parallelepiped = up_parallelepiped - down_parallelepiped;


							//right_side
							std::function<Point<double>(bool& is_solve, int id_element, Point<double> X)> sourse = [&size_parallelepiped, &up_parallelepiped, &down_parallelepiped, &element](bool& is_solve, int id_element, Point<double> X) -> Point<double>
							{
								is_solve = true;
								Point<double> result;

								double hx = size_parallelepiped.x;
								double hy = size_parallelepiped.y;
								double hz = size_parallelepiped.z;
								double lambda, mu;
								element->self_grid.GetDomain(element->self_grid.GetElement(id_element)->GetIdDomain())->forMech.GetLameCoefficients(lambda, mu);

								double ksi_x = (X.x - down_parallelepiped.x) / hx;
								double ksi_y = (X.y - down_parallelepiped.y) / hy;
								double ksi_z = (X.z - down_parallelepiped.z) / hz;

								double phi_x, phi_y, phi_z;
								double dphi_x, dphi_y, dphi_z;
								double ddphi_x, ddphi_y, ddphi_z;
								double dphi_dxdx;
								double dphi_dydy;
								double dphi_dzdz;
								double dphi_dxdy;
								double dphi_dxdz;
								double dphi_dydz;
								phi_x = -4 * ksi_x * ksi_x + 4 * ksi_x;
								dphi_x = (-8 / hx * ksi_x + 4 / hx);
								ddphi_x = -8 / (hx * hx);
								phi_y = -4 * ksi_y * ksi_y + 4 * ksi_y;
								dphi_y = (-8 / hy * ksi_y + 4 / hy);
								ddphi_y = -8 / (hy * hy);
								phi_z = -4 * ksi_z * ksi_z + 4 * ksi_z;
								dphi_z = (-8 / hz * ksi_z + 4 / hz);
								ddphi_z = -8 / (hz * hz);

								dphi_dxdx = ddphi_x * phi_y * phi_z;
								dphi_dxdy = dphi_x * dphi_y * phi_z;
								dphi_dxdz = dphi_x * phi_y * dphi_z;
								dphi_dydy = phi_x * ddphi_y * phi_z;
								dphi_dydz = phi_x * dphi_y * dphi_z;
								dphi_dzdz = phi_x * phi_y * ddphi_z;


								double laplasian_phi = dphi_dxdx + dphi_dydy + dphi_dzdz;
								Point<double> DD_phi;
								DD_phi.x = dphi_dxdx + dphi_dxdy + dphi_dxdz;
								DD_phi.y = dphi_dxdy + dphi_dydy + dphi_dydz;
								DD_phi.z = dphi_dxdz + dphi_dydz + dphi_dzdz;

								result = DD_phi * (lambda + mu) + Point<double>(laplasian_phi, laplasian_phi, laplasian_phi) * mu;

								return result * (-1.0);
							};

							FEM::FEM_3D_forElasticDeformation(
								false,
								1.e-15,
								base_grid_simple, //input
								first_boundaries, //input
								second_boundary_empty, //input
								sourse,
								element->self_direction, //output
								element->self_grid, //output
								element->self_basis_functions[element->GetNodesCount() + id_self_edges.size() + element->GetLowerElementCount()], //output
								Base_stiffness_matrix //output
							);

							if (this->is_print_logFile) {
								printf_s("Print the mech result into .dat file... ");

								FILE* fout_tech;
								char name_u_tech[5000];
								sprintf_s(name_u_tech, "%s/Volume_bubble_bf%d.dat", element->self_direction, element->GetNodesCount() + id_self_edges.size() + element->GetLowerElementCount());
								fopen_s(&fout_tech, name_u_tech, "w");
								char name_in_file[1000];
								sprintf_s(name_in_file, "Volume_bubble_bf%d", element->GetNodesCount() + (int)id_self_edges.size() + element->GetLowerElementCount());
								std::vector<std::vector<char>> name_value(4);
								char name_v_tmp[6][100];
								sprintf_s(name_v_tmp[0], "Ux");
								sprintf_s(name_v_tmp[1], "Uy");
								sprintf_s(name_v_tmp[2], "Uz");
								sprintf_s(name_v_tmp[3], "Material");
								for (int i = 0; i < name_value.size(); i++)
								{
									name_value[i].resize(100);
									for (int j = 0; j < name_value[i].size(); j++)
									{
										name_value[i][j] = name_v_tmp[i][j];
									}
								}
								std::vector<std::vector<double>> value(4);
								value[0].resize(element->self_grid.GetVertexCount());
								value[1].resize(element->self_grid.GetVertexCount());
								value[2].resize(element->self_grid.GetVertexCount());
								value[3].resize(element->self_grid.GetElementsCount());

								for (int i = 0; i < element->self_grid.GetVertexCount(); i++)
								{
									value[0][i] = element->self_basis_functions[element->GetNodesCount() + id_self_edges.size() + element->GetLowerElementCount()][i].x;
									value[1][i] = element->self_basis_functions[element->GetNodesCount() + id_self_edges.size() + element->GetLowerElementCount()][i].y;
									value[2][i] = element->self_basis_functions[element->GetNodesCount() + id_self_edges.size() + element->GetLowerElementCount()][i].z;
								}
								for (int i = 0; i < element->self_grid.GetElementsCount(); i++)
								{
									auto el = element->self_grid.GetElement(i);
									value[3][i] = el->GetIdDomain();

								}
								element->self_grid.printTecPlot3D_DiffDomains(fout_tech, value, name_value, name_in_file);
								fclose(fout_tech);
							}
						}
					}

					//save solutions
					for (int id_bf = 0; id_bf < element->self_basis_functions.size(); id_bf++)
					{
						FILE* fout;
						char name_out[1000];
						sprintf_s(name_out, "%s/BasisFunction_%d.txt", element->self_direction, id_bf);
						fopen_s(&fout, name_out, "w");
						for (int i = 0; i < element->self_basis_functions[id_bf].size(); i++)
						{
							fprintf_s(fout, "%.15e %.15e %.15e\n", element->self_basis_functions[id_bf][i].x, element->self_basis_functions[id_bf][i].y, element->self_basis_functions[id_bf][i].z);
						}
						fclose(fout);
					}
				}

				//create MACRO basis functions
				for (int _id_bf = 0; _id_bf < element->self_basis_functions.size(); _id_bf++)
				{
					std::function< Point<double>(Point<double> X)> bf = [_id_bf, element](Point<double> X) -> Point<double>
					{
						Point<double> result;

						double len;
						int id_micro_elem = element->self_grid.GetNearestElementID(X, len);
						result = element->self_grid.GetSolutionInPoint(id_micro_elem, X, element->self_basis_functions[_id_bf]);

						/*if (_id_bf <= 8)
						{
							result = element->self_grid.GetSolutionInPoint(id_micro_elem, X, element->self_basis_functions[_id_bf]);
						}
						else {
							auto f = [&](double x, double xmax, double xmin) -> double*
							{
								double phi[3];
								double x_loc = (x - xmin) / (xmax - xmin);
								phi[0] = (1 - x_loc);
								phi[1] = (x_loc);
								phi[2] = -4 * x_loc * x_loc + 4 * x_loc;

								return phi;
							};

							int numeration[3][27] = {
							{ 0,1,0,1,0,1,0,1, 2,0,0,1,1,2,0,1,2,0,1,2, 2,2,0,2,1,2, 2 },
							{ 0,0,1,1,0,0,1,1, 0,2,0,2,0,1,1,1,0,2,2,1, 2,2,2,1,2,0, 2 },
							{ 0,0,0,0,1,1,1,1, 0,0,2,0,2,0,2,2,1,1,1,1, 0,1,2,2,2,2, 2 } };

							result.x = f(X.x, element->GetNode(7).x, element->GetNode(0).x)[numeration[0][_id_bf]] *
								f(X.y, element->GetNode(7).y, element->GetNode(0).y)[numeration[1][_id_bf]] *
								f(X.z, element->GetNode(7).z, element->GetNode(0).z)[numeration[2][_id_bf]];
						}*/

						result.y = result.z;
						result.x = result.z;

						return result;
					};
					std::function< Point<Point<double>>(Point<double> X)> derivative_bf = [element, _id_bf](Point<double> X) -> Point<Point<double>>
					{
						Point<Point<double>> result;

						double len;
						int id_micro_elem = element->self_grid.GetNearestElementID(X, len);
						result = element->self_grid.GetDerevativeFromSolutionInPoint(id_micro_elem, X, element->self_basis_functions[_id_bf]);

						/*if (_id_bf <= 8)
						{
							result = element->self_grid.GetDerevativeFromSolutionInPoint(id_micro_elem, X, element->self_basis_functions[_id_bf]);
						}
						else {
							auto f = [&](double x, double xmax, double xmin) -> double*
							{
								double phi[3];
								double x_loc = (x - xmin) / (xmax - xmin);
								phi[0] = (1 - x_loc);
								phi[1] = (x_loc);
								phi[2] = -4 * x_loc * x_loc + 4 * x_loc;

								return phi;
							};
							auto df = [&](double x, double xmax, double xmin) -> double*
							{
								double phi[3];
								double dx_loc = 1 / (xmax - xmin);
								double x_loc = (x - xmin) / (xmax - xmin);
								phi[0] = -dx_loc;
								phi[1] = +dx_loc;
								phi[2] = -8 * x_loc * dx_loc + 4;
								return phi;
							};

							int numeration[3][27] = {
								{ 0,1,0,1,0,1,0,1, 2,0,0,1,1,2,0,1,2,0,1,2, 2,2,0,2,1,2, 2 },
								{ 0,0,1,1,0,0,1,1, 0,2,0,2,0,1,1,1,0,2,2,1, 2,2,2,1,2,0, 2 },
								{ 0,0,0,0,1,1,1,1, 0,0,2,0,2,0,2,2,1,1,1,1, 0,1,2,2,2,2, 2 } };

							double resx, resy, resz;
							resx = df(X.x, element->GetNode(7).x, element->GetNode(0).x)[numeration[0][_id_bf]];
							resy = f(X.y, element->GetNode(7).y, element->GetNode(0).y)[numeration[1][_id_bf]];
							resz = f(X.z, element->GetNode(7).z, element->GetNode(0).z)[numeration[2][_id_bf]];
							result.x.x = resx * resy * resz;

							resx = f(X.x, element->GetNode(7).x, element->GetNode(0).x)[numeration[0][_id_bf]];
							resy = df(X.y, element->GetNode(7).y, element->GetNode(0).y)[numeration[1][_id_bf]];
							resz = f(X.z, element->GetNode(7).z, element->GetNode(0).z)[numeration[2][_id_bf]];
							result.x.y = resx * resy * resz;

							resx = f(X.x, element->GetNode(7).x, element->GetNode(0).x)[numeration[0][_id_bf]];
							resy = df(X.y, element->GetNode(7).y, element->GetNode(0).y)[numeration[1][_id_bf]];
							resz = f(X.z, element->GetNode(7).z, element->GetNode(0).z)[numeration[2][_id_bf]];
							result.x.z = resx * resy * resz;
						}*/

						result.y = result.z;
						result.x = result.z;

						return result;
					};

					int target_global_id;
					if (_id_bf < element->GetNodesCount())
						target_global_id = element->GetIdNode(_id_bf);
					else if (_id_bf < id_self_edges.size() + element->GetNodesCount())
						target_global_id = this->edges[id_self_edges[_id_bf - element->GetNodesCount()]].GetSelfGlobalId() + this->GetVertexCount();
					else if (_id_bf < id_self_edges.size() + element->GetNodesCount() + element->GetLowerElementCount())
						target_global_id = element->GetLowerElement(_id_bf - element->GetNodesCount() - id_self_edges.size())->GetSelfGlobalId() + this->GetVertexCount() + this->edges.size();
					else
						target_global_id = element->GetSelfGlobalId() + this->GetVertexCount() + this->edges.size() + this->faces.size();
					element->SetDOF(_id_bf, target_global_id, bf, derivative_bf);
				}
			}

			for (int id_face = 0; id_face < this->faces.size(); id_face++)
			{
				this->faces[id_face].ResizeDOF(this->faces[id_face].GetLowerElementCount() + this->faces[id_face].GetNodesCount() + 1);
				for (int id_vertex_local = 0; id_vertex_local < this->faces[id_face].GetNodesCount(); id_vertex_local++)
				{
					int id_vertex_global = this->faces[id_face].GetIdNode(id_vertex_local);
					auto local_bf = this->GetElement(this->faces[id_face].GetUpperElement(0)->GetSelfGlobalId())->GetBasisFunctionInGlobalID(id_vertex_global);
					auto local_Dbf = this->GetElement(this->faces[id_face].GetUpperElement(0)->GetSelfGlobalId())->GetDerivativeOfBasisFunctionInGlobalID(id_vertex_global);
					this->faces[id_face].SetDOF(id_vertex_local, id_vertex_global, *local_bf, *local_Dbf);
				}
				for (int id_edge_local = 0; id_edge_local < this->faces[id_face].GetLowerElementCount(); id_edge_local++)
				{
					int id_edge_global = this->faces[id_face].GetLowerElement(id_edge_local)->GetSelfGlobalId();
					auto local_bf = this->GetElement(this->faces[id_face].GetUpperElement(0)->GetSelfGlobalId())->GetBasisFunctionInGlobalID(this->vertexes.size() + id_edge_global);
					auto local_Dbf = this->GetElement(this->faces[id_face].GetUpperElement(0)->GetSelfGlobalId())->GetDerivativeOfBasisFunctionInGlobalID(this->vertexes.size() + id_edge_global);
					this->faces[id_face].SetDOF(id_edge_local + this->faces[id_face].GetNodesCount(), this->vertexes.size() + id_edge_global, *local_bf, *local_Dbf);
				}

				auto local_bf = this->GetElement(this->faces[id_face].GetUpperElement(0)->GetSelfGlobalId())->GetBasisFunctionInGlobalID(this->vertexes.size() + this->edges.size() + id_face);
				auto local_Dbf = this->GetElement(this->faces[id_face].GetUpperElement(0)->GetSelfGlobalId())->GetDerivativeOfBasisFunctionInGlobalID(this->vertexes.size() + +this->edges.size() + id_face);
				this->faces[id_face].SetDOF(this->faces[id_face].GetLowerElementCount() + this->faces[id_face].GetNodesCount(), this->vertexes.size() + this->edges.size() + id_face, *local_bf, *local_Dbf);
			}
		}
		template <typename Dirichlet, typename Neumann>
		void CreateBoundaryConditions(std::vector<Dirichlet>& first_boundary, std::vector<Neumann>& second_boundary)
		{
			//1st boundary
			//по смыслу правильно, но не работает
			if (false) 
			{
				for (int id_type = 1; id_type < first_boundary.size(); id_type++)
				{
					std::vector<int> global_dof_id;
					std::vector<int> target_faces;
					for (int i = 0; i < first_boundary[id_type].id_vertexes.size(); i++)
					{
						global_dof_id.push_back(first_boundary[id_type].id_vertexes[i]);
					}
					for (int id_edge = 0; id_edge < this->edges.size(); id_edge++)
					{
						int id_edge_0 = this->edges[id_edge].GetLowerElement(0)->GetSelfGlobalId();
						int id_edge_1 = this->edges[id_edge].GetLowerElement(1)->GetSelfGlobalId();

						bool is_find_0 = false;
						bool is_find_1 = false;
						for (int i = 0; i < first_boundary[id_type].id_vertexes.size() && (is_find_0 == false || is_find_1 == false); i++)
						{
							int boundary_vert = first_boundary[id_type].id_vertexes[i];
							if (boundary_vert == id_edge_0)
							{
								is_find_0 = true;
							}
							if (boundary_vert == id_edge_1)
							{
								is_find_1 = true;
							}
						}
						if (is_find_0 && is_find_1)
						{
							global_dof_id.push_back(id_edge + this->vertexes.size());
						}
					}
					for (int id_face = 0; id_face < this->faces.size(); id_face++)
					{
						std::vector<int> face(this->faces[id_face].GetNodesCount());
						for (int i = 0; i < face.size(); i++)
						{
							face[i] = this->faces[id_face].GetIdNode(i);
						}
						auto find = math::GetConfluence(face, first_boundary[id_type].id_vertexes);
						if (find.size() == face.size())
						{
							target_faces.push_back(id_face);
							global_dof_id.push_back(id_face + this->edges.size() + this->vertexes.size());
						}
					}
					math::MakeQuickSort(global_dof_id);

					DenseMatrix_Tensor Matrix(global_dof_id.size());
					
					for (int id_face_local = 0; id_face_local < target_faces.size(); id_face_local++)
					{
						auto face = &(this->faces[target_faces[id_face_local]]);

						face->SetIntegrationLaw(3);
						
						DenseMatrix<Tensor2Rank3D, Point<double>> local_matix(face->GetDOFsCount());
						for (int I = 0; I < face->GetDOFsCount(); I++)
						{
							//auto basis_function_I = face->GetBasisFunctionInLocalID(I);
							int id_elem = face->GetUpperElement(0)->GetSelfGlobalId();
							auto basis_function_I = this->GetElement(id_elem)->GetBasisFunctionInGlobalID(face->GetDOFInLocalID(I));
							for (int J = 0; J < face->GetDOFsCount(); J++)
							{
								auto basis_function_J = this->GetElement(face->GetUpperElement(0)->GetSelfGlobalId())->GetBasisFunctionInGlobalID(face->GetDOFInLocalID(J));

								std::function<Tensor2Rank3D(Point<double>)> MassMatrix = [&basis_function_I, &basis_function_J]
								(Point<double> X) -> Tensor2Rank3D
								{
									Tensor2Rank3D result;

									Point<double> BF_I_inX = (*basis_function_I)(X);
									Point<double> BF_J_inX = (*basis_function_J)(X);

									math::GetTensorMult(BF_J_inX, BF_I_inX, result);
									result.SetValue(BF_I_inX.x * BF_J_inX.x);

									return result;
								};
								std::function<double(Point<double>)> Volume = []
								(Point<double> X) -> double
								{
									return 1.0;
								};

								double V = face->GetVolume();
								double V_solve = face->SolveIntegral(Volume);
								local_matix.A[I][J] = face->SolveIntegral(MassMatrix);
							}
							std::function<Point<double>(Point<double>)> RightSide = [&basis_function_I, &first_boundary, &id_type]
							(Point<double> X) -> Point<double>
							{
								Point<double> result;

								Point<bool> is_in;
								Point<double> BF_I_inX = (*basis_function_I)(X);
								Point<double> F = Point<double>(1, 1, 1);//first_boundary[id_type].value(is_in);

								result.x = BF_I_inX.x * F.x;
								result.y = BF_I_inX.x * F.y;
								result.z = BF_I_inX.x * F.z;

								return result;
							};
							local_matix.F[I] = face->SolveIntegral(RightSide);
						}

						for (int I = 0; I < local_matix.GetSize(); I++)
						{
							int id = face->GetDOFInLocalID(I);
							int i_global = math::GetPositionInSortVector(global_dof_id, id);
							if (i_global != -1)
							{
								for (int J = 0; J < face->GetDOFsCount(); J++)
								{
									id = face->GetDOFInLocalID(J);
									int j_global = math::GetPositionInSortVector(global_dof_id, id);
									if (j_global != -1)
									{
										Matrix.A[i_global][j_global] += local_matix.A[I][J];
									}
								}
								Matrix.F[i_global] += local_matix.F[I];
							}
						}
					}

					//Matrix.Gauss_LU();
					Matrix.Gauss();

					for (int i = 0; i < global_dof_id.size(); i++)
					{
						MsFEM::BoundaryVertex_forMech_Poly temp;
						temp.boundary_value = [Matrix, i, &first_boundary, id_type](Point<bool>& is_enter) -> Point<double>
						{ 
							first_boundary[id_type].value(is_enter);
							return Point<double>(Matrix.X[i].x, Matrix.X[i].x, Matrix.X[i].x);
							//return Matrix.X[i];
						};
						
						temp.ResizeDOF(1);
						if (global_dof_id[i] < this->vertexes.size())
						{
							temp.SetIdNode(0, global_dof_id[i]);
							temp.SetNode(0, this->GetPtrCoordinateViaID(global_dof_id[i]));
							int target_edge_id = this->vertexes[global_dof_id[i]].GetUpperElement(0)->GetSelfGlobalId();
							int target_face_id = this->edges[target_edge_id].GetUpperElement(0)->GetSelfGlobalId();
							int target_elem_id = this->faces[target_face_id].GetUpperElement(0)->GetSelfGlobalId();
							temp.SetDOF(0, global_dof_id[i], *this->GetElement(target_elem_id)->GetBasisFunctionInGlobalID(global_dof_id[i]), *this->GetElement(target_elem_id)->GetDerivativeOfBasisFunctionInGlobalID(global_dof_id[i]));
						}
						else if (global_dof_id[i] < (this->vertexes.size() + this->edges.size()))
						{
							int target_face_id = this->edges[global_dof_id[i] - this->vertexes.size()].GetUpperElement(0)->GetSelfGlobalId();
							int target_elem_id = this->faces[target_face_id].GetUpperElement(0)->GetSelfGlobalId();
							temp.SetDOF(0, global_dof_id[i], *this->GetElement(target_elem_id)->GetBasisFunctionInGlobalID(global_dof_id[i]), *this->GetElement(target_elem_id)->GetDerivativeOfBasisFunctionInGlobalID(global_dof_id[i]));
						}
						else
						{
							int id_face = global_dof_id[i] - this->vertexes.size() - this->edges.size();
							int target_elem_id = this->faces[id_face].GetUpperElement(0)->GetSelfGlobalId();
							temp.SetDOF(0, global_dof_id[i], *this->GetElement(target_elem_id)->GetBasisFunctionInGlobalID(global_dof_id[i]), *this->GetElement(target_elem_id)->GetDerivativeOfBasisFunctionInGlobalID(global_dof_id[i]));
						}
						this->boundary_vertexes.push_back(temp);
					}

					std::vector<Point<double>> test_points(10);
					test_points[0] = Point<double>(0, 0, 2);
					test_points[1] = Point<double>(2, 0, 2);
					test_points[2] = Point<double>(0, 2, 2);
					test_points[3] = Point<double>(2, 2, 2);
					test_points[4] = Point<double>(1, 0, 2);
					test_points[5] = Point<double>(1, 2, 2);
					test_points[6] = Point<double>(0, 1, 2);
					test_points[7] = Point<double>(2, 1, 2);
					test_points[8] = Point<double>(1, 1, 2);
					test_points[9] = Point<double>(0.5, 0.5, 2);
					std::vector<Point<double>> U(test_points.size());
					Point<bool> test;
					for (int i = 0; i < test_points.size(); i++)
					{
						for (int j = 0; j < global_dof_id.size(); j++)
						{
							Point<double> bf_j = (*this->boundary_vertexes[j].GetBasisFunctionInLocalID(0))(test_points[i]);
							Point<double> val = this->boundary_vertexes[j].boundary_value(test);
							U[i].x += bf_j.x * val.x;
							U[i].y += bf_j.x * val.y;
							U[i].z += bf_j.x * val.z;
						}
					}
				}
			}
			//все компоненты считаем по х-базису
			if (true)
			{
				for (int id_type = 0; id_type < first_boundary.size(); id_type++)
				{
					std::vector<int> global_dof_id;
					std::vector<Point<double>> values;
					std::vector<int> target_faces;
					for (int i = 0; i < first_boundary[id_type].id_vertexes.size(); i++)
					{
						global_dof_id.push_back(first_boundary[id_type].id_vertexes[i]);
					}
					for (int id_edge = 0; id_edge < this->edges.size(); id_edge++)
					{
						int id_edge_0 = this->edges[id_edge].GetLowerElement(0)->GetSelfGlobalId();
						int id_edge_1 = this->edges[id_edge].GetLowerElement(1)->GetSelfGlobalId();

						bool is_find_0 = false;
						bool is_find_1 = false;
						for (int i = 0; i < first_boundary[id_type].id_vertexes.size() && (is_find_0 == false || is_find_1 == false); i++)
						{
							int boundary_vert = first_boundary[id_type].id_vertexes[i];
							if (boundary_vert == id_edge_0)
							{
								is_find_0 = true;
							}
							if (boundary_vert == id_edge_1)
							{
								is_find_1 = true;
							}
						}
						if (is_find_0 && is_find_1)
						{
							global_dof_id.push_back(id_edge + this->vertexes.size());
						}
					}
					for (int id_face = 0; id_face < this->faces.size(); id_face++)
					{
						std::vector<int> face(this->faces[id_face].GetNodesCount());
						for (int i = 0; i < face.size(); i++)
						{
							face[i] = this->faces[id_face].GetIdNode(i);
						}
						auto find = math::GetConfluence(face, first_boundary[id_type].id_vertexes);
						if (find.size() == face.size())
						{
							target_faces.push_back(id_face);
							global_dof_id.push_back(id_face + this->edges.size() + this->vertexes.size());
						}
					}
					math::MakeQuickSort(global_dof_id);

					DenseMatrix<double, double> Matrix(global_dof_id.size());
					values.resize(global_dof_id.size());

					//собрали матрицу массы
					for (int id_face_local = 0; id_face_local < target_faces.size(); id_face_local++)
					{
						auto face = &(this->faces[target_faces[id_face_local]]);

						face->SetIntegrationLaw(3);

						DenseMatrix<double, double> local_matix(face->GetDOFsCount());
						for (int I = 0; I < face->GetDOFsCount(); I++)
						{
							//auto basis_function_I = face->GetBasisFunctionInLocalID(I);
							int id_elem = face->GetUpperElement(0)->GetSelfGlobalId();
							auto basis_function_I = this->GetElement(id_elem)->GetBasisFunctionInGlobalID(face->GetDOFInLocalID(I));
							for (int J = 0; J < face->GetDOFsCount(); J++)
							{
								auto basis_function_J = this->GetElement(face->GetUpperElement(0)->GetSelfGlobalId())->GetBasisFunctionInGlobalID(face->GetDOFInLocalID(J));

								std::function<double(Point<double>)> MassMatrix = [&basis_function_I, &basis_function_J]
								(Point<double> X) -> double
								{
									Point<double> BF_I_inX = (*basis_function_I)(X);
									Point<double> BF_J_inX = (*basis_function_J)(X);

									return BF_I_inX.x * BF_J_inX.x;
								};
								std::function<double(Point<double>)> Volume = []
								(Point<double> X) -> double
								{
									return 1.0;
								};

								double V = face->GetVolume();
								double V_solve = face->SolveIntegral(Volume);
								local_matix.A[I][J] = face->SolveIntegral(MassMatrix);
							}
						}

						for (int I = 0; I < local_matix.GetSize(); I++)
						{
							int id = face->GetDOFInLocalID(I);
							int i_global = math::GetPositionInSortVector(global_dof_id, id);
							if (i_global != -1)
							{
								for (int J = 0; J < face->GetDOFsCount(); J++)
								{
									id = face->GetDOFInLocalID(J);
									int j_global = math::GetPositionInSortVector(global_dof_id, id);
									if (j_global != -1)
									{
										Matrix.A[i_global][j_global] += local_matix.A[I][J];
									}
								}
							}
						}
					}

					//ищем покомпонентные разложения по базису
					for (int xx = 0; xx < 3; xx++)
					{
						//правая часть
						for (int i = 0; i < Matrix.F.size(); i++)
						{
							Matrix.F[i] = 0.0;
						}
						for (int id_face_local = 0; id_face_local < target_faces.size(); id_face_local++)
						{
							auto face = &(this->faces[target_faces[id_face_local]]);

							face->SetIntegrationLaw(3);

							DenseMatrix<double, double> local_matix(face->GetDOFsCount());
							for (int I = 0; I < face->GetDOFsCount(); I++)
							{
								//auto basis_function_I = face->GetBasisFunctionInLocalID(I);
								int id_elem = face->GetUpperElement(0)->GetSelfGlobalId();
								auto basis_function_I = this->GetElement(id_elem)->GetBasisFunctionInGlobalID(face->GetDOFInLocalID(I));

								std::function<double(Point<double>)> RightSide = [&basis_function_I, &first_boundary, &id_type, xx]
								(Point<double> X) -> double
								{
									double result;

									Point<bool> is_in;
									Point<double> BF_I_inX = (*basis_function_I)(X);
									Point<double> F = first_boundary[id_type].value(is_in);

									switch (xx)
									{
									case 0: result = F.x * BF_I_inX.x; break;
									case 1: result = F.y * BF_I_inX.x; break;
									case 2: result = F.z * BF_I_inX.x; break;
									default:
										break;
									}

									return result;
								};
								local_matix.F[I] = face->SolveIntegral(RightSide);
							}

							for (int I = 0; I < local_matix.GetSize(); I++)
							{
								int id = face->GetDOFInLocalID(I);
								int i_global = math::GetPositionInSortVector(global_dof_id, id);
								if (i_global != -1)
								{
									Matrix.F[i_global] += local_matix.F[I];
								}
							}
						}

						//Matrix.Gauss_LU();
						Matrix.Gauss();

						for (int i = 0; i < values.size(); i++)
						{
							switch (xx)
							{
							case 0: values[i].x = Matrix.X[i]; break;
							case 1: values[i].y = Matrix.X[i]; break;
							case 2: values[i].z = Matrix.X[i]; break;
							default:
								break;
							}

						}
					}

					for (int i = 0; i < global_dof_id.size(); i++)
					{
						MsFEM::BoundaryVertex_forMech_Poly temp;
						temp.boundary_value = [values, i, &first_boundary, id_type](Point<bool>& is_enter) -> Point<double>
						{
							first_boundary[id_type].value(is_enter);
							return values[i];
							//return Matrix.X[i];
						};

						temp.ResizeDOF(1);
						if (global_dof_id[i] < this->vertexes.size())
						{
							temp.SetIdNode(0, global_dof_id[i]);
							temp.SetNode(0, this->GetPtrCoordinateViaID(global_dof_id[i]));
							int target_edge_id = this->vertexes[global_dof_id[i]].GetUpperElement(0)->GetSelfGlobalId();
							int target_face_id = this->edges[target_edge_id].GetUpperElement(0)->GetSelfGlobalId();
							int target_elem_id = this->faces[target_face_id].GetUpperElement(0)->GetSelfGlobalId();
							temp.SetDOF(0, global_dof_id[i], *this->GetElement(target_elem_id)->GetBasisFunctionInGlobalID(global_dof_id[i]), *this->GetElement(target_elem_id)->GetDerivativeOfBasisFunctionInGlobalID(global_dof_id[i]));
						}
						else if (global_dof_id[i] < (this->vertexes.size() + this->edges.size()))
						{
							int target_face_id = this->edges[global_dof_id[i] - this->vertexes.size()].GetUpperElement(0)->GetSelfGlobalId();
							int target_elem_id = this->faces[target_face_id].GetUpperElement(0)->GetSelfGlobalId();
							temp.SetDOF(0, global_dof_id[i], *this->GetElement(target_elem_id)->GetBasisFunctionInGlobalID(global_dof_id[i]), *this->GetElement(target_elem_id)->GetDerivativeOfBasisFunctionInGlobalID(global_dof_id[i]));
						}
						else
						{
							int id_face = global_dof_id[i] - this->vertexes.size() - this->edges.size();
							int target_elem_id = this->faces[id_face].GetUpperElement(0)->GetSelfGlobalId();
							temp.SetDOF(0, global_dof_id[i], *this->GetElement(target_elem_id)->GetBasisFunctionInGlobalID(global_dof_id[i]), *this->GetElement(target_elem_id)->GetDerivativeOfBasisFunctionInGlobalID(global_dof_id[i]));
						}
						this->boundary_vertexes.push_back(temp);
					}

					/*std::vector<Point<double>> test_points(10);
					test_points[0] = Point<double>(0, 0, 2);
					test_points[1] = Point<double>(2, 0, 2);
					test_points[2] = Point<double>(0, 2, 2);
					test_points[3] = Point<double>(2, 2, 2);
					test_points[4] = Point<double>(1, 0, 2);
					test_points[5] = Point<double>(1, 2, 2);
					test_points[6] = Point<double>(0, 1, 2);
					test_points[7] = Point<double>(2, 1, 2);
					test_points[8] = Point<double>(1, 1, 2);
					test_points[9] = Point<double>(0.5, 0.5, 2);
					std::vector<Point<double>> U(test_points.size());
					Point<bool> test;
					for (int i = 0; i < test_points.size(); i++)
					{
						for (int j = 0; j < global_dof_id.size(); j++)
						{
							Point<double> bf_j = (*this->boundary_vertexes[j].GetBasisFunctionInLocalID(0))(test_points[i]);
							Point<double> val = this->boundary_vertexes[j].boundary_value(test);
							U[i].x += bf_j.x * val.x;
							U[i].y += bf_j.y * val.y;
							U[i].z += bf_j.z * val.z;
						}
					}*/
				}
			}
		}
		template <typename Dirichlet, typename Neumann>
		void Initialization(math::SimpleGrid& base_grid, std::vector<Dirichlet>& first_boundary, std::vector<Neumann>& second_boundary, char *fine_mesh_dir)
		{
			try {

				//create geometry space
				math::MakeCopyVector_A_into_B(base_grid.xyz, *(this->GetCoordinates()));
				this->CreateXYZline();
				this->SetElementsCount((int)base_grid.nvtr.size());
				for (int id_element = 0; id_element < this->GetElementsCount(); id_element++)
				{
					auto new_element = this->GetElement(id_element);

					std::vector<int> nodes_global(base_grid.nvtr[id_element][0]);
					std::vector<Point<double>*> P(base_grid.nvtr[id_element][0]);
					int ii = 1;
					for (int n = 0; n < nodes_global.size(); n++)
					{
						nodes_global[n] = base_grid.nvtr[id_element][ii];
						ii++;
						P[n] = this->GetPtrCoordinateViaID(nodes_global[n]);
					}
					std::vector<std::vector<int>> nodes_faces(base_grid.nvtr[id_element][ii]);
					ii++;
					for (int f = 0; f < nodes_faces.size(); f++)
					{
						nodes_faces[f].resize(base_grid.nvtr[id_element][ii]);
						ii++;
						for (int nf = 0; nf < nodes_faces[f].size(); nf++)
						{
							nodes_faces[f][nf] = base_grid.nvtr[id_element][ii];
							ii++;
						}
					}

					new_element->SetGeometry(nodes_global, nodes_faces, P);
					new_element->SetIdDomain(base_grid.nvkat[id_element]);
				}
				this->CreateQTree();

				//create topology
				CreateTopology();
					
				//create functional spaces
				CreateFunctionalSpaces(fine_mesh_dir);

				CreateBoundaryConditions(first_boundary, second_boundary);
			}
			catch (const std::exception&)
			{
				printf_s("Error: MultiXFEM/MultiXFEM_Grid.h/void Initialization(geometry::Grid<geometry::Tetrahedron> &base_grid)\n");
			}
		}
		
		template <typename Matrix>
		void CreationPortrait(Matrix& matrix)
		{
			try {
				printf_s("\n");
				std::vector<std::vector<int>> tmp_down_columns(this->GetDOFsCount()), tmp_up_columns(this->GetDOFsCount());
				std::vector<std::vector<int>> down_columns(this->GetDOFsCount()), up_columns(this->GetDOFsCount());
				for (int id_elem = 0; id_elem < this->GetElementsCount(); id_elem++)
				{
					if (id_elem % 100 == 0)
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
					if (id_string % 1 == 0)
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

		Point<double> GetSolutionInPoint(int id_element, Point<double> X, std::vector<Point<double>>& solution)
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
		Point<Point<double>> GetDerevativeFromSolutionInPoint(int id_element, Point<double> X, std::vector<Point<double>>& solution)
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

		void printTecPlot3D(/*char *directory,*/ FILE* fdat, std::vector<std::vector<double>>& value, std::vector<std::vector<char>> name_value, char* name_zone)
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

		~Grid_forMech()
		{
			DOFs_count = 0;
			std::vector<BoundaryVertex_forMech_Poly> v_b1;
			std::vector<BoundaryVertex_forMech_Poly>(v_b1).swap(this->boundary_vertexes);
			std::vector<BoundaryFace_forMech_Poly> v_b2;
			std::vector<BoundaryFace_forMech_Poly>(v_b2).swap(this->boundary_faces);
			std::vector<int> v_p;
			std::vector<int>(v_p).swap(this->accordance_DOF_and_vertex);
			std::vector<Vertex> v_v;
			std::vector<Vertex>(v_v).swap(this->vertexes);
			std::vector<Edge> v_e;
			std::vector<Edge>(v_e).swap(this->edges);
			std::vector<Face> v_f;
			std::vector<Face>(v_f).swap(this->faces);

			this->DeleteGeometriGrid();
		}
	private:

	};

	class FiniteElement_forMech_Poly_Order1 :
		public geometry::Polyhedron,
		public topology::Polyhedron<topology::lower::Polygon, topology::upper::EmptyElement>,
		public functional::Shape<int, Point<double>>
	{
	public:
		math::SimpleGrid SubGrid_for_integration;
		math::SimpleGrid SubGrid_for_integr_byTriangles;

		FEM::Grid_forMech self_grid;
		std::vector< std::vector<Point<double>>> self_basis_functions;

		char self_direction[1000];

		FiniteElement_forMech_Poly_Order1() { return; };
		~FiniteElement_forMech_Poly_Order1() { return; };

		void SolveLocalMatrix(DenseMatrix<Tensor2Rank3D, Point<double>>& local_matix, std::function<std::vector<std::vector<double>>(Point<double>)>& koefD)
		{
			try {
				local_matix.SetSize(this->GetDOFsCount());

				this->SetIntegrationLaw(4);

				/*if (this->GetIdDomain() == 2)
					printf_s("\nThis elem[%d] is singular\n", this->GetSelfGlobalId());*/

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

							std::vector<std::vector<double>> mult_res;
							double cohesive_koef = -7.4e+6;
							math::ResizeVector(mult_res, 3, 3);
							mult_res[0][0] = BF_I_inX.x * BF_J_inX.x * cohesive_koef;
							mult_res[1][1] = BF_I_inX.y * BF_J_inX.y * cohesive_koef;
							mult_res[2][2] = BF_I_inX.z * BF_J_inX.z * cohesive_koef;

							result = mult_res;
							//result = 1.0;
							return result;
						};

						double V = this->GetVolume();
						local_matix.A[_bf_local_id_I][_bf_local_id_J] = this->SolveIntegral(StiffnessMatrix);
						local_matix.A[_bf_local_id_J][_bf_local_id_I] = local_matix.A[_bf_local_id_I][_bf_local_id_J].T();
					}
				}
				/*for (int _bf_local_id_I = 0; _bf_local_id_I < this->GetDOFsCount(); _bf_local_id_I++)
				{
					for (int _bf_local_id_J = _bf_local_id_I + 1; _bf_local_id_J < this->GetDOFsCount(); _bf_local_id_J++)
					{
						local_matix.A[_bf_local_id_I][_bf_local_id_J] = local_matix.A[_bf_local_id_J][_bf_local_id_I].T();
					}
					for (int i = 0; i < 3; i++)
					{
						for (int j = i; j < 3; j++)
						{
							local_matix.A[_bf_local_id_I][_bf_local_id_I].val[i][j] = local_matix.A[_bf_local_id_I][_bf_local_id_I].val[j][i];
						}
					}
				}*/

				for (int i = 0; i < local_matix.A.size(); i++)
				{
					Tensor2Rank3D tmp;
					for (int j = 0; j < local_matix.A[i].size(); j++)
					{
						tmp += local_matix.A[i][j];
					}
				}
			}
			catch (const char* err)
			{
				std::cerr << "Error (SolveLocalMatrix): " << err << std::endl;
				int a;
				scanf_s("%d", &a);
			}
		}
		void SolveLocalCohesiveMatrix(DenseMatrix<Tensor2Rank3D, Point<double>>& local_matix, std::function<double(Point<double>)>& koefCohesive)
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
							mult_res[0][0] = BF_I_inX.x * BF_J_inX.x * koefCohesive(X);
							mult_res[1][1] = BF_I_inX.y * BF_J_inX.y * koefCohesive(X);
							mult_res[2][2] = BF_I_inX.z * BF_J_inX.z * koefCohesive(X);

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
	class BoundaryVertex_forMech_Poly_Order1 :
		public geometry::Vertex,
		//public topology::Vertex<topology::lower::EmptyElement, topology::upper::Segment>,
		public topology::Vertex<topology::lower::EmptyElement, topology::upper::Polyhedron>,
		public functional::Shape<int, Point<double>>
	{
	public:
		std::function<Point<double>(Point<bool>&)> boundary_value;

		BoundaryVertex_forMech_Poly_Order1() { return; };
		~BoundaryVertex_forMech_Poly_Order1() { return; };

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
	class BoundaryFace_forMech_Poly_Order1:
		public geometry::Polygon,
		public topology::Polygon<topology::lower::Vertex, topology::upper::Polyhedron>,
		public functional::Shape<int, Point<double>>
	{
	public:
		std::function<Point<double>(Point<double>)> boundary_value;

		BoundaryFace_forMech_Poly_Order1() { return; };
		~BoundaryFace_forMech_Poly_Order1() { return; };

		void SetTypeBoundary(int type)
		{
			this->type_boundary = type;
		}
		int GetTypeBoundary()
		{
			return this->type_boundary;
		}

		void SolveLocalBoundaryVector(std::vector<Point<double>>& local_vector, std::function<Point<double>(Point<double>)>& boundary_value)
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
	class Grid_forMech_Poly_Order1 : public geometry::Grid<FiniteElement_forMech_Poly_Order1>
	{
		int DOFs_count;
		struct Face : //2D elements
			public topology::Polygon<topology::lower::Polyline, topology::upper::Polyhedron>,
			public geometry::Polygon,
			public functional::Shape<int, Point<double>>
		{
		public:
			Face() {};
			~Face() {};
		};
		std::vector<Face> faces;
		struct Edge : //1D elements
			public topology::Polyline<topology::lower::Vertex, topology::upper::Polygon>,
			public geometry::Segment,
			public functional::Shape<int, Point<double>>
		{
		public:
			Edge() {};
			~Edge() {};
		};
		std::vector<Edge> edges;
		struct Vertex : //0D elements
			public topology::Vertex<topology::lower::EmptyElement, topology::upper::Polyline>,
			public geometry::Vertex,
			public functional::Shape<int, Point<double>>
		{
		public:
			Vertex() {};
			~Vertex() {};
		};
		std::vector<Vertex> vertexes;

	public:
		std::vector<BoundaryVertex_forMech_Poly_Order1> boundary_vertexes;
		std::vector<BoundaryFace_forMech_Poly_Order1> boundary_faces;

		bool is_print_logFile;

		std::vector<int> accordance_DOF_and_vertex;

		Grid_forMech_Poly_Order1()
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

		void CreateTopology()
		{
			try {
				//create elements
				for (int id_element = 0; id_element < this->GetElementsCount(); id_element++)
				{
					this->GetElement(id_element)->SetSelfGlobalId(id_element);
					this->GetElement(id_element)->SetLowerElementCount(this->GetElement(id_element)->GetCountFaces());
				}
				//create faces topology
				{
					struct face_temp
					{
						std::vector<int> num;
						int element, local_in_element;
					public:
						face_temp()
						{
							element = 0; local_in_element = 0;
						}
						face_temp(int element, int local_in_element,
							std::vector<int>& n)
						{
							this->element = element; this->local_in_element = local_in_element;
							math::MakeCopyVector_A_into_B(n, this->num);
						}

						void operator= (face_temp A)
						{
							element = A.element;
							local_in_element = A.local_in_element;
							this->num.resize(A.num.size());
							for (int i = 0; i < A.num.size(); i++)
							{
								this->num[i] = A.num[i];
							};
						}
						bool operator< (face_temp A)
						{
							for (int i = 0; i < num.size(); i++)
							{
								if (num[i] < A.num[i]) return true;
								if (num[i] > A.num[i]) return false;
							}
							return false;
						}

						bool operator>(face_temp A)
						{
							for (int i = 0; i < num.size(); i++)
							{
								if (num[i] > A.num[i]) return true;
								if (num[i] < A.num[i]) return false;
							}
							return false;
						}
						bool operator!= (face_temp A)
						{
							for (int i = 0; i < num.size(); i++)
							{
								if (num[i] != A.num[i]) return true;
							}
							return false;
						}
					};
					struct face
					{
						std::vector<int> num;
						int self_global_id;

						std::vector<int> element, local_in_element;
					public:
						face(face_temp A)
						{
							element.push_back(A.element);
							local_in_element.push_back(A.local_in_element);

							num.resize(A.num.size());
							for (int i = 0; i < num.size(); i++)
								num[i] = A.num[i];
						}
						face()
						{
						}
						void set_elem(face_temp A)
						{
							element.push_back(A.element);
							local_in_element.push_back(A.local_in_element);

						}
					};
					auto MakeRemovalOfDuplication_f = [](std::vector<face_temp>& A, std::vector<face>& B) -> void
					{
						if (A.size() != 0)
						{
							int k = 0, j = 0;

							B.push_back(face());
							B[k] = A[0];
							k++;

							for (int i = 1; i < A.size(); i++)
							{
								if (math::GetConfluence(A[j].num, A[i].num).size() != A[j].num.size()) //не повтор
								{
									j = i;
									B.push_back(face());
									B[k] = A[i];
									k++;
								}
								else {
									B[k - 1].local_in_element.push_back(A[i].local_in_element);
									B[k - 1].element.push_back(A[i].element);
								}
							}
						}
					};

					std::vector <face_temp> f_temp;
					std::vector <face> faces;
					topology::lower::Polyhedron tmp_3D;
					topology::lower::Polygon tmp_2D;

					for (int id_element = 0; id_element < this->GetElementsCount(); id_element++)
					{
						auto element3D = this->GetElement(id_element);

						for (int id_face = 0; id_face < element3D->GetCountFaces(); id_face++)
						{
							int _id_elem = id_element;
							int _id_face_in_elem = id_face;
							std::vector<int> id_in_face;
							element3D->GetGlobalIdFace(id_face, id_in_face);
							//math::MakeQuickSort(id_in_face);
							f_temp.push_back(face_temp(_id_elem, _id_face_in_elem, id_in_face));
						}
					}
					math::MakeQuickSort(f_temp);
					MakeRemovalOfDuplication_f(f_temp, faces);

					this->faces.resize(faces.size());
					for (int id_face = 0; id_face < faces.size(); id_face++)
					{
						this->faces[id_face].SetUpperElementCount((int)faces[id_face].element.size());
						this->faces[id_face].SetSelfGlobalId(id_face);

						//add Upper and Lower Elements
						for (int id_volume = 0; id_volume < faces[id_face].element.size(); id_volume++)
						{
							this->faces[id_face].SetUpperElement(id_volume, this->GetElement(faces[id_face].element[id_volume]));
							this->GetElement(faces[id_face].element[id_volume])->SetLowerElement(faces[id_face].local_in_element[id_volume], &this->faces[id_face]);
						}

						//save geometry
						std::vector<Point<double>*> P;
						this->GetPtrCoordinateViaID(faces[id_face].num, P);
						this->faces[id_face].SetGeometry(faces[id_face].num, P);
					}
				}

				//create edge topology
				if (true)
				{
					struct edge_temp
					{
						int num1, num2;
						int element, local_in_element;
						int face, local_in_face;
					public:
						edge_temp()
						{
							element = 0; local_in_element = 0;
							face = 0; local_in_face = 0;
							num1 = 0;
							num2 = 0;
						}
						edge_temp(int element, int local_in_element,
							int face, int local_in_face,
							std::vector<int>& n)
						{
							this->element = element; this->local_in_element = local_in_element;
							this->face = face; this->local_in_face = local_in_face;

							num1 = n[0];
							num2 = n[1];
							if (n[0] > n[1])
							{
								num1 = n[1];
								num2 = n[0];
							}
						}

						void operator= (edge_temp A)
						{
							this->element = A.element; this->local_in_element = A.local_in_element;
							this->face = A.face; this->local_in_face = A.local_in_face;

							num1 = A.num1;
							num2 = A.num2;
						}
						bool operator< (edge_temp A)
						{
							if (num1 < A.num1) return true;
							if (num1 > A.num1) return false;
							if (num1 == A.num1 && num2 < A.num2) return true;
							return false;
						}

						bool operator>(edge_temp A)
						{
							if (num1 > A.num1) return true;
							if (num1 < A.num1) return false;
							if (num1 == A.num1 && num2 > A.num2) return true;
							return false;
						}
						bool operator!= (edge_temp A)
						{
							if (num1 != A.num1 || num2 != A.num2) return true;
							return false;
						}

						bool operator == (edge_temp A)
						{
							if (num1 == A.num1 && num2 == A.num2) return true;

							return false;
						}
					};
					struct edge
					{
						int num1, num2;
						int self_global_id;

						std::vector<int> element, local_in_element;
						std::vector<int> face, local_in_face;
					public:
						edge(edge_temp A)
						{
							//element.push_back(A.element);
							//local_in_element.push_back(A.local_in_element);
							//face.push_back(A.face);
							//local_in_face.push_back(A.local_in_face);
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
							//face.push_back(A.face);
							//local_in_face.push_back(A.local_in_face);
							face.push_back(A.face);
							local_in_face.push_back(A.local_in_face);

						}
					};
					auto MakeRemovalOfDuplication_e = [](std::vector<edge_temp>& A, std::vector<edge>& B) -> void
					{
						if (A.size() != 0)
						{
							int k = 0, j = 0;

							B.push_back(edge());
							B[k] = A[0];
							k++;

							for (int i = 1; i < A.size(); i++)
							{
								if (A[j].num1 != A[i].num1 || A[j].num2 != A[i].num2) //не повтор
								{
									j = i;
									B.push_back(edge());
									B[k] = A[i];
									k++;
								}
								else {
									B[k - 1].local_in_face.push_back(A[i].local_in_face);
									B[k - 1].face.push_back(A[i].face);
								}
							}
						}
					};

					std::vector <edge_temp> e_temp;
					std::vector <edge> edges_vector;

					for (int id_face = 0; id_face < this->faces.size(); id_face++)
					{
						auto element2D = &this->faces[id_face];
						element2D->SetLowerElementCount(element2D->GetNodesCount());

						int local_id_faces_in_elem_0 = -1;
						auto element3D = this->GetElement(element2D->GetUpperElement(0)->GetSelfGlobalId());
						for (int id = 0; id < element3D->GetLowerElementCount(); id++)
						{
							if (element2D->GetSelfGlobalId() == element3D->GetLowerElement(id)->GetSelfGlobalId())
							{
								local_id_faces_in_elem_0 = id;
								break;
							}
						}

						std::vector<int> id_nodes;
						for (int id_edge = 0; id_edge < element2D->GetLowerElementCount(); id_edge++)
						{
							int _id_elem = -1;
							int _id_face_in_elem = -1;
							int _id_face = id_face;
							int _id_edge_in_face = id_edge;
							element3D->GetGlobalIdEdge(local_id_faces_in_elem_0, id_edge, id_nodes);
							e_temp.push_back(edge_temp(
								_id_elem, _id_face_in_elem,
								_id_face, _id_edge_in_face,
								id_nodes));
						}
					}
					math::MakeQuickSort(e_temp);
					MakeRemovalOfDuplication_e(e_temp, edges_vector);

					//add Upper and Lower Elements
					this->edges.resize(edges_vector.size());
					for (int id_edge = 0; id_edge < edges_vector.size(); id_edge++)
					{
						this->edges[id_edge].SetUpperElementCount((int)edges_vector[id_edge].face.size());
						this->edges[id_edge].SetSelfGlobalId(id_edge);

						for (int id_face = 0; id_face < edges_vector[id_edge].face.size(); id_face++)
						{
							this->edges[id_edge].SetUpperElement(id_face, &this->faces[edges_vector[id_edge].face[id_face]]);
							this->faces[edges_vector[id_edge].face[id_face]].SetLowerElement(edges_vector[id_edge].local_in_face[id_face], &this->edges[id_edge]);
						}

						//save geometry
						std::vector<int> id_nodes(2);
						id_nodes[0] = edges_vector[id_edge].num1;
						id_nodes[1] = edges_vector[id_edge].num2;
						std::vector<Point<double>*> P;
						this->GetPtrCoordinateViaID(id_nodes, P);
						this->edges[id_edge].SetGeometry(id_nodes, P);
					}
				}

				//create vertexes topology
				if (true)
				{
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
					auto MakeRemovalOfDuplication_n = [](std::vector<node_temp>& A, std::vector<node>& B) -> void
					{
						if (A.size() != 0)
						{
							int k = 0, j = 0;

							B.push_back(node());
							B[k] = A[0];
							k++;

							for (int i = 1; i < A.size(); i++)
							{
								if (A[j].num != A[i].num) //не повтор
								{
									j = i;
									B.push_back(node());
									B[k] = A[i];
									k++;
								}
								else {
									B[k - 1].local_in_edge.push_back(A[i].local_in_edge);
									B[k - 1].edge.push_back(A[i].edge);
								}
							}
						}
					};

					std::vector <node_temp> n_temp;
					std::vector <node> nodes_vector;

					std::vector <int> id_nodes_in_edge(2);

					for (int id_edge = 0; id_edge < this->edges.size(); id_edge++)
					{
						auto element1D = &this->edges[id_edge];
						auto element2D = element1D->GetUpperElement(0);
						auto element3D = element2D->GetUpperElement(0);

						for (int id_vertex = 0; id_vertex < element1D->GetLowerElementCount(); id_vertex++)
						{
							int _id_elem = -1;
							int _id_face_in_elem = -1;
							int _id_face = -1;
							int _id_edge_in_face = -1;
							int _id_edge = id_edge;
							int _id_vertex_in_edge = id_vertex;
							n_temp.push_back(node_temp(_id_elem, _id_face_in_elem,
								_id_face, _id_edge_in_face,
								_id_edge, _id_vertex_in_edge,
								this->edges[id_edge].GetIdNode(id_vertex)));
						}
					}
					math::MakeQuickSort(n_temp);
					MakeRemovalOfDuplication_n(n_temp, nodes_vector);

					this->vertexes.resize(nodes_vector.size());
					//add Upper and Lower Elements
					for (int id_node = 0; id_node < nodes_vector.size(); id_node++)
					{
						this->vertexes[id_node].SetUpperElementCount((int)nodes_vector[id_node].edge.size());
						this->vertexes[id_node].SetSelfGlobalId(id_node);

						for (int id_edge = 0; id_edge < nodes_vector[id_node].edge.size(); id_edge++)
						{
							this->vertexes[id_node].SetUpperElement(id_edge, &this->edges[nodes_vector[id_node].edge[id_edge]]);
							this->edges[nodes_vector[id_node].edge[id_edge]].SetLowerElement(nodes_vector[id_node].local_in_edge[id_edge], &this->vertexes[id_node]);
						}

						//save geometry
						std::vector<int> id_nodes(1);
						id_nodes[0] = nodes_vector[id_node].num;
						std::vector<Point<double>*> P;
						this->GetPtrCoordinateViaID(id_nodes, P);
						this->vertexes[id_node].SetGeometry(id_nodes, P);
					}
				}

				//sort local_edges for faces
				for (int id_face = 0; id_face < faces.size(); id_face++)
				{
					std::vector<std::vector<int>> old_positions(this->faces[id_face].GetLowerElementCount());
					std::vector<int> new_positions(this->faces[id_face].GetLowerElementCount());
					for (int i = 0; i < old_positions.size(); i++)
					{
						old_positions[i].resize(2);
						old_positions[i][0] = this->faces[id_face].GetLowerElement(i)->GetSelfGlobalId();
						old_positions[i][1] = 0;
					}
					new_positions[0] = old_positions[0][0];
					old_positions[0][1] = 1;
					int id_target_vertex = this->edges[new_positions[0]].GetLowerElement(1)->GetSelfGlobalId();
					for (int i = 1; i < new_positions.size(); i++)
					{
						for (int j = 1; j < old_positions.size(); j++)
						{
							if (old_positions[j][1] == 0)
							{
								if (this->edges[old_positions[j][0]].GetLowerElement(0)->GetSelfGlobalId() == id_target_vertex)
								{
									id_target_vertex = this->edges[old_positions[j][0]].GetLowerElement(1)->GetSelfGlobalId();
									new_positions[i] = old_positions[j][0];
									old_positions[j][1] = 1;
									break;
								}
								else if (this->edges[old_positions[j][0]].GetLowerElement(1)->GetSelfGlobalId() == id_target_vertex) {
									id_target_vertex = this->edges[old_positions[j][0]].GetLowerElement(0)->GetSelfGlobalId();
									new_positions[i] = old_positions[j][0];
									old_positions[j][1] = 1;
									break;
								}
							}
						}
					}

					for (int i = 0; i < new_positions.size(); i++)
					{
						this->faces[id_face].SetLowerElement(i, &this->edges[new_positions[i]]);
					}
				}


			}
			catch (const std::exception&)
			{
				printf_s("Error: MultiXFEM/MultiXFEM_Grid.h/void CreateTopology()\n");
			}
		}
		void CreateFunctionalSpaces(char* fine_mesh_dir)
		{
			try {
				//create basis functions
				this->SetDOFsCount((int)this->vertexes.size());
				omp_set_num_threads(math::NUM_THREADS);
#pragma omp parallel for
				for (int id_element = 0; id_element < this->GetElementsCount(); id_element++)
				{
					printf_s("Create basis functions for element[%d]                  \r", id_element);

					auto element = this->GetElement(id_element);
					sprintf_s(element->self_direction, sizeof(element->self_direction), "%s/Elem_%d/SubElem_0", fine_mesh_dir, id_element);

					std::vector<int> id_self_edges;
					{
						std::vector<int> id_self_edges_tmp;
						for (int id_face_local = 0; id_face_local < element->GetLowerElementCount(); id_face_local++)
						{
							int id_face = (element->GetLowerElement(id_face_local))->GetSelfGlobalId();
							//printf_s("id_face(global) = %d\r", id_face);
							for (int id_edge_local = 0; id_edge_local < this->faces[id_face].GetLowerElementCount(); id_edge_local++)
							{
								int id_edge = this->faces[id_face].GetLowerElement(id_edge_local)->GetSelfGlobalId();
								id_self_edges_tmp.push_back(id_edge);
							}
						}
						math::MakeQuickSort(id_self_edges_tmp);
						math::MakeRemovalOfDuplication(id_self_edges_tmp, id_self_edges);
					}
					std::vector<std::vector<int>> boundary_vertexes_in_macro(id_self_edges.size());

					int DOF_count = element->GetNodesCount(); //vertexes
					element->ResizeDOF(DOF_count);
					element->self_basis_functions.resize(DOF_count);

					//read sub meshes
					math::SimpleGrid base_grid_simple;
					if (!base_grid_simple.ReadFromNVTR(element->self_direction, 4))
					{
						printf_s("\nThere is a problem with self-NVTR for element[%d]\n", id_element);
						int a;
						scanf_s("%d", &a);
					}
					std::vector<math::SimpleGrid> faces_grid_simple(element->GetLowerElementCount());
					std::vector<math::SimpleGrid> edges_grid_simple(id_self_edges.size());
					{
						for (int id_face_local = 0; id_face_local < element->GetLowerElementCount(); id_face_local++)
						{
							char name_in_file[1000];
							sprintf_s(name_in_file, "%s/Face_id%d.dat", element->self_direction, id_face_local);
							faces_grid_simple[id_face_local].ReadFromSalomeDat(name_in_file, 2, this->GetVertexCount());
						}

						for (int id_edge_local = 0; id_edge_local < id_self_edges.size(); id_edge_local++)
						{
							int id_edge = id_self_edges[id_edge_local];
							std::vector<int> neigbor_faces_local_in_elem;
							for (int id_face_local = 0; id_face_local < this->edges[id_edge].GetUpperElementCount() && neigbor_faces_local_in_elem.size() < 2; id_face_local++)
							{
								for (int id_element_local = 0; id_element_local < this->edges[id_edge].GetUpperElement(id_face_local)->GetUpperElementCount(); id_element_local++)
								{
									int id_element = this->edges[id_edge].GetUpperElement(id_face_local)->GetUpperElement(id_element_local)->GetSelfGlobalId();
									if (id_element == element->GetSelfGlobalId())
									{
										int global_id = this->edges[id_edge].GetUpperElement(id_face_local)->GetSelfGlobalId();
										int local_id = -1;
										for (int i = 0; i < element->GetLowerElementCount(); i++)
										{
											if (element->GetLowerElement(i)->GetSelfGlobalId() == global_id)
											{
												local_id = i;
												break;
											}
										}
										neigbor_faces_local_in_elem.push_back(local_id);
										break;
									}
								}
							}

							for (int id_elem_f1 = 0; id_elem_f1 < faces_grid_simple[neigbor_faces_local_in_elem[0]].nvtr.size(); id_elem_f1++)
							{
								for (int id_elem_f2 = 0; id_elem_f2 < faces_grid_simple[neigbor_faces_local_in_elem[1]].nvtr.size(); id_elem_f2++)
								{
									auto segment = math::GetConfluence(faces_grid_simple[neigbor_faces_local_in_elem[0]].nvtr[id_elem_f1], faces_grid_simple[neigbor_faces_local_in_elem[1]].nvtr[id_elem_f2]);

									if (segment.size() == 2)
									{
										std::vector<Point<double>> X(2);
										X[0] = base_grid_simple.xyz[segment[0]];
										X[1] = base_grid_simple.xyz[segment[1]];

										if (X[1] < X[0])
										{
											Point<double> tmp = X[0];
											X[0] = X[1];
											X[1] = tmp;

											int tmp2 = segment[0];
											segment[0] = segment[1];
											segment[1] = tmp2;
										}

										edges_grid_simple[id_edge_local].nvtr.push_back(segment);
										edges_grid_simple[id_edge_local].nvkat.push_back(faces_grid_simple[neigbor_faces_local_in_elem[0]].nvkat[id_elem_f1]);
										break;
									}
									if (segment.size() > 2)
									{
										printf_s("\nError in create segment!!!\n");
										//sleep(100000);
									}
								}
							}

							//create map for edge
							std::vector<int> vertex_in_edge;
							vertex_in_edge.reserve(edges_grid_simple[id_edge_local].nvtr.size() * 2);
							for (int i = 0; i < edges_grid_simple[id_edge_local].nvtr.size(); i++)
							{
								vertex_in_edge.push_back(edges_grid_simple[id_edge_local].nvtr[i][0]);
								vertex_in_edge.push_back(edges_grid_simple[id_edge_local].nvtr[i][1]);
							}
							math::MakeQuickSort(vertex_in_edge);
							math::MakeRemovalOfDuplication(vertex_in_edge, edges_grid_simple[id_edge_local].vertex_map);
							edges_grid_simple[id_edge_local].xyz.resize(edges_grid_simple[id_edge_local].vertex_map.size());
							for (int i = 0; i < edges_grid_simple[id_edge_local].xyz.size(); i++)
							{
								edges_grid_simple[id_edge_local].xyz[i] = base_grid_simple.xyz[edges_grid_simple[id_edge_local].vertex_map[i]];
							}
						}

						//renumeration of the NVTR into self
						for (int id_face_local = 0; id_face_local < faces_grid_simple.size(); id_face_local++)
						{
							for (int i = 0; i < faces_grid_simple[id_face_local].nvtr.size(); i++)
							{
								for (int j = 0; j < faces_grid_simple[id_face_local].nvtr[i].size(); j++)
								{
									faces_grid_simple[id_face_local].nvtr[i][j] = faces_grid_simple[id_face_local].TransferIdGlobalIntoSelf(faces_grid_simple[id_face_local].nvtr[i][j]);
								}
							}
						}
						for (int id_edge_local = 0; id_edge_local < edges_grid_simple.size(); id_edge_local++)
						{
							for (int i = 0; i < edges_grid_simple[id_edge_local].nvtr.size(); i++)
							{
								for (int j = 0; j < edges_grid_simple[id_edge_local].nvtr[i].size(); j++)
								{
									edges_grid_simple[id_edge_local].nvtr[i][j] = edges_grid_simple[id_edge_local].TransferIdGlobalIntoSelf(edges_grid_simple[id_edge_local].nvtr[i][j]);
								}
							}
						}
					}

					//create basis funtions
					bool is_solve = false;

					for (int id_bf = 0; id_bf < element->GetDOFsCount(); id_bf++)
					{
						FILE* fin;
						char name_in[2000];
						sprintf_s(name_in, "%s/BasisFunction_%d.txt", element->self_direction, id_bf);
						fopen_s(&fin, name_in, "r");
						if (fin == NULL)
						{
							is_solve = true;
							break;
						}
						else {
							fclose(fin);
						}
					}
					if (is_solve == false)
					{
						for (int i = 0; i < this->GetDomainsCount(); i++)
							element->self_grid.AddDomain(*this->GetDomain(i));
						struct _Dirichlet {
							int global_id, in_elem_id;
							std::vector<int> id_vertexes;
							std::vector<Point<double>> xyz;
							std::function<Point<double>(Point<bool>&, int id_vertex)> value;
						};
						std::vector<_Dirichlet> first_boundaries_empty;
						struct _Neumann {
							/*double value;
							Point<double> vector;*/
							std::function<Point<double>(Point<double>)> value;
							std::vector<std::function<Point<double>(Point<double>)>> values;
							std::vector<std::vector<int>> id_vertexes_as_triangle;
							std::vector<int> id_base_element;
						};
						std::vector<_Neumann> second_boundary_empty;

						element->self_grid.Initialization(base_grid_simple, first_boundaries_empty, second_boundary_empty);

						for (int id_bf = 0; id_bf < element->GetDOFsCount(); id_bf++)
						{
							FILE* fin;
							char name_in[1000];
							sprintf_s(name_in, "%s/BasisFunction_%d.txt", element->self_direction, id_bf);
							fopen_s(&fin, name_in, "r");

							element->self_basis_functions[id_bf].resize(element->self_grid.GetDOFsCount());

							for (int i = 0; i < element->self_grid.GetDOFsCount(); i++)
							{
								fscanf_s(fin, "%lf %lf %lf", &element->self_basis_functions[id_bf][i].x, &element->self_basis_functions[id_bf][i].y, &element->self_basis_functions[id_bf][i].z);
							}

							fclose(fin);
						}
					}
					if (is_solve == true)
					{
						//solve edges problems
						std::vector<FEM::Grid_1D_forMech> solver_grid_edges(id_self_edges.size()); //output [num_edge]
						std::vector<std::vector<std::vector<Point<double>>>> Solution_edges(id_self_edges.size()); //output [num_edge][num_BF][num_DOF]
						std::vector<std::vector<int>> boundary_vertexes_local(id_self_edges.size());
						for (int id_edge_local = 0; id_edge_local < id_self_edges.size(); id_edge_local++)
						{
							Solution_edges[id_edge_local].resize(2);

							for (int i = 0; i < this->GetDomainsCount(); i++)
								solver_grid_edges[id_edge_local].AddDomain(*this->GetDomain(i));


							CSSD_Matrix<Tensor2Rank3D, Point<double>> Base_stiffness_matrix;
							Base_stiffness_matrix.print_logs = this->is_print_logFile;

							//create Boundary values and vertexes
							struct _Dirichlet {
								std::vector<int> id_vertexes;
								std::function<Point<double>(Point<bool>&)> value;
							};
							std::vector<int> elements_for_vertex(edges_grid_simple[id_edge_local].xyz.size());
							for (int i = 0; i < edges_grid_simple[id_edge_local].nvtr.size(); i++)
							{
								elements_for_vertex[edges_grid_simple[id_edge_local].nvtr[i][0]]++;
								elements_for_vertex[edges_grid_simple[id_edge_local].nvtr[i][1]]++;
							}
							for (int i = 0; i < elements_for_vertex.size(); i++)
							{
								if (elements_for_vertex[i] == 1)
								{
									boundary_vertexes_local[id_edge_local].push_back(edges_grid_simple[id_edge_local].TransferIdSelfIntoGlobal(i));
								}
							}
							std::vector<std::vector<int>> vertex_in_boundaries(2);
							std::vector<_Dirichlet> first_boundaries(2);
							first_boundaries[0].id_vertexes.push_back(edges_grid_simple[id_edge_local].TransferIdGlobalIntoSelf(boundary_vertexes_local[id_edge_local][0]));
							first_boundaries[1].id_vertexes.push_back(edges_grid_simple[id_edge_local].TransferIdGlobalIntoSelf(boundary_vertexes_local[id_edge_local][1]));
							boundary_vertexes_in_macro[id_edge_local].resize(2);
							if (math::IsEqual(edges_grid_simple[id_edge_local].xyz[first_boundaries[0].id_vertexes[0]], this->edges[id_self_edges[id_edge_local]].GetNode(0)))
							{
								boundary_vertexes_in_macro[id_edge_local][0] = this->edges[id_self_edges[id_edge_local]].GetIdNode(0);
								boundary_vertexes_in_macro[id_edge_local][1] = this->edges[id_self_edges[id_edge_local]].GetIdNode(1);
							}
							else {
								boundary_vertexes_in_macro[id_edge_local][1] = this->edges[id_self_edges[id_edge_local]].GetIdNode(0);
								boundary_vertexes_in_macro[id_edge_local][0] = this->edges[id_self_edges[id_edge_local]].GetIdNode(1);
							}

							//First function (linear)
							printf_s("create linear function for edge %d(local: %d) - bf linear 0\n", id_self_edges[id_edge_local], id_edge_local);
							{
								first_boundaries[0].value = [](Point<bool>& arg)->Point<double>
								{
									arg.x = true;
									arg.y = true;
									arg.z = true;
									Point<double> result;
									result.x = 1;
									result.y = 1;
									result.z = 1;
									return result;
								};
								first_boundaries[1].value = [](Point<bool>& arg)->Point<double>
								{
									arg.x = true;
									arg.y = true;
									arg.z = true;
									Point<double> result;
									result.x = 0;
									result.y = 0;
									result.z = 0;
									return result;
								};

								struct _Neumann {
									/*double value;
									Point<double> vector;*/
									std::function<Point<double>(Point<double>)> value;
									std::vector<std::function<Point<double>(Point<double>)>> values;
									std::vector<std::vector<int>> id_vertexes_as_triangle;
									std::vector<int> id_base_element;
								};
								std::vector<_Neumann> second_boundary_empty;

								//right_side
								std::function<Point<double>(bool& is_solve, int id_element, Point<double> X)> sourse = [](bool& is_solve, int id_element, Point<double> X) -> Point<double>
								{
									is_solve = false;
									return Point<double>(0, 0, 0);
								};

								FEM::FEM_1D_forElasticDeformation(
									false,
									1e-15,
									edges_grid_simple[id_edge_local], //input
									first_boundaries, //input
									second_boundary_empty, //input
									sourse,
									element->self_direction, //output
									solver_grid_edges[id_edge_local], //output
									Solution_edges[id_edge_local][0], //output
									Base_stiffness_matrix //output
								);
							}

							//Second function (linear)
							printf_s("create linear function for edge %d(local: %d) - bf linear 1\n", id_self_edges[id_edge_local], id_edge_local);
							{
								first_boundaries[0].value = [](Point<bool>& arg) -> Point<double>
								{
									arg.x = true;
									arg.y = true;
									arg.z = true;
									Point<double> result;
									result.x = 0;
									result.y = 0;
									result.z = 0;
									return result;
								};
								first_boundaries[1].value = [](Point<bool>& arg) -> Point<double>
								{
									arg.x = true;
									arg.y = true;
									arg.z = true;
									Point<double> result;
									result.x = 1;
									result.y = 1;
									result.z = 1;
									return result;
								};

								struct _Neumann {
									/*double value;
									Point<double> vector;*/
									std::function<Point<double>(Point<double>)> value;
									std::vector<std::function<Point<double>(Point<double>)>> values;
									std::vector<std::vector<int>> id_vertexes_as_triangle;
									std::vector<int> id_base_element;
								};
								std::vector<_Neumann> second_boundary_empty;

								//right_side
								std::function<Point<double>(bool& is_solve, int id_element, Point<double> X)> sourse = [](bool& is_solve, int id_element, Point<double> X) -> Point<double>
								{
									is_solve = false;
									return Point<double>(0, 0, 0);
								};

								FEM::FEM_1D_forElasticDeformation(
									false,
									1e-15,
									edges_grid_simple[id_edge_local], //input
									first_boundaries, //input
									second_boundary_empty, //input
									sourse,
									element->self_direction, //output
									solver_grid_edges[id_edge_local], //output
									Solution_edges[id_edge_local][1], //output
									Base_stiffness_matrix //output
								);
							}

							////Выравниевание на х (костыль!)
							/*for (int i = 0; i < Solution_edges[id_edge_local][0].size(); i++)
							{
								Solution_edges[id_edge_local][0][i].y = Solution_edges[id_edge_local][0][i].x;
								Solution_edges[id_edge_local][0][i].z = Solution_edges[id_edge_local][0][i].x;

								Solution_edges[id_edge_local][1][i].y = Solution_edges[id_edge_local][1][i].x;
								Solution_edges[id_edge_local][1][i].z = Solution_edges[id_edge_local][1][i].x;
							}*/
						}

						//solve faces problems
						std::vector<FEM::Grid_2D_forMech> solver_grid_faces(faces_grid_simple.size()); //output [num_face]
						std::vector<std::vector<std::vector<Point<double>>>> Solution_faces(faces_grid_simple.size()); //output [num_face][num_BF][num_DOF]
						std::vector<std::vector<int>> target_id_for_BF(solver_grid_faces.size());
						for (int id_face_local = 0; id_face_local < solver_grid_faces.size(); id_face_local++)
						{
							int id_face = element->GetLowerElement(id_face_local)->GetSelfGlobalId();
							Solution_faces[id_face_local].resize(faces[id_face].GetLowerElementCount());

							for (int i = 0; i < this->GetDomainsCount(); i++)
								solver_grid_faces[id_face_local].AddDomain(*this->GetDomain(i));

							CSSD_Matrix<Tensor2Rank3D, Point<double>> Base_stiffness_matrix;
							Base_stiffness_matrix.print_logs = this->is_print_logFile;

							//create Boundary values and vertexes
							printf_s("create Boundary values and vertexes for face %d(local: %d)\n", id_face, id_face_local);
							struct _Dirichlet {
								int global_id, in_elem_id;
								std::vector<int> id_vertexes;
								std::vector<Point<double>> xyz;
								std::function<Point<double>(Point<bool>&, int id_vertex)> value;
							};
							std::vector<_Dirichlet> first_boundaries(faces[id_face].GetLowerElementCount());
							for (int id_edge_local = 0; id_edge_local < first_boundaries.size(); id_edge_local++)
							{
								first_boundaries[id_edge_local].global_id = faces[id_face].GetLowerElement(id_edge_local)->GetSelfGlobalId();
								first_boundaries[id_edge_local].in_elem_id = math::GetPositionInSortVector(id_self_edges, first_boundaries[id_edge_local].global_id);

								for (int i = 0; i < edges_grid_simple[first_boundaries[id_edge_local].in_elem_id].xyz.size(); i++)
								{
									int global_id_vertex = edges_grid_simple[first_boundaries[id_edge_local].in_elem_id].TransferIdSelfIntoGlobal(i);
									first_boundaries[id_edge_local].id_vertexes.push_back(faces_grid_simple[id_face_local].TransferIdGlobalIntoSelf(global_id_vertex));
									//first_boundaries[id_edge_local].xyz.push_back(faces_grid_simple[id_face_local].xyz[first_boundaries[id_edge_local].id_vertexes[first_boundaries[id_edge_local].id_vertexes.size() - 1]]);
								}
							}

							//linear function
							for (int id_bf = 0; id_bf < faces[id_face].GetLowerElementCount(); id_bf++)
							{
								printf_s("create linear functions for face %d(local: %d) - bf linear %d\n", id_face, id_face_local, id_bf);
								for (int i = 0; i < first_boundaries.size(); i++)
								{
									first_boundaries[i].value = [](Point<bool>& arg, int id) -> Point<double>
									{
										arg.x = true;
										arg.y = true;
										arg.z = true;
										return Point<double>(0, 0, 0);
									};
								}

								int left_edge_id_global = this->faces[id_face].GetLowerElement((id_bf - 1) < 0 ? faces[id_face].GetLowerElementCount() - 1 : id_bf - 1)->GetSelfGlobalId();
								int left_edge_id = math::GetPositionInSortVector(id_self_edges, left_edge_id_global);

								int right_edge_id_global = this->faces[id_face].GetLowerElement(id_bf)->GetSelfGlobalId();
								int right_edge_id = math::GetPositionInSortVector(id_self_edges, right_edge_id_global);

								int target_vertex_global = math::GetConfluence(boundary_vertexes_in_macro[left_edge_id], boundary_vertexes_in_macro[right_edge_id])[0];
								target_id_for_BF[id_face_local].push_back(target_vertex_global);

								first_boundaries[id_bf].value = [&Solution_edges, &boundary_vertexes_in_macro, &edges_grid_simple, &faces_grid_simple, id_face_local, right_edge_id, target_vertex_global](Point<bool>& arg, int id_vertex) -> Point<double>
								{
									Point<double> result;
									arg.x = true;
									arg.y = true;
									arg.z = true;

									int id_bf_in_edge;
									if (target_vertex_global == boundary_vertexes_in_macro[right_edge_id][0]) id_bf_in_edge = 0;
									else id_bf_in_edge = 1;

									int id_vertex_global = faces_grid_simple[id_face_local].TransferIdSelfIntoGlobal(id_vertex);

									result = Solution_edges[right_edge_id][id_bf_in_edge][edges_grid_simple[right_edge_id].TransferIdGlobalIntoSelf(id_vertex_global)];
									return result;
								};
								first_boundaries[(id_bf - 1) < 0 ? faces[id_face].GetLowerElementCount() - 1 : id_bf - 1].value = [&Solution_edges, &boundary_vertexes_in_macro, &edges_grid_simple, &faces_grid_simple, id_face_local, left_edge_id, target_vertex_global](Point<bool>& arg, int id_vertex) -> Point<double>
								{
									Point<double> result;
									arg.x = true;
									arg.y = true;
									arg.z = true;

									int id_bf_in_edge;
									if (target_vertex_global == boundary_vertexes_in_macro[left_edge_id][0]) id_bf_in_edge = 0;
									else id_bf_in_edge = 1;

									int id_vertex_global = faces_grid_simple[id_face_local].TransferIdSelfIntoGlobal(id_vertex);

									result = Solution_edges[left_edge_id][id_bf_in_edge][edges_grid_simple[left_edge_id].TransferIdGlobalIntoSelf(id_vertex_global)];
									return result;
								};

								struct _Neumann {
									/*double value;
									Point<double> vector;*/
									std::function<Point<double>(Point<double>)> value;
									std::vector<std::function<Point<double>(Point<double>)>> values;
									std::vector<std::vector<int>> id_vertexes_as_triangle;
									std::vector<int> id_base_element;
								};
								std::vector<_Neumann> second_boundary_empty;

								//right_side
								std::function<Point<double>(bool& is_solve, int id_element, Point<double> X)> sourse = [](bool& is_solve, int id_element, Point<double> X) -> Point<double>
								{
									is_solve = false;
									return Point<double>(0, 0, 0);
								};

								FEM::FEM_2D_forElasticDeformation(
									false,
									1e-15,
									faces_grid_simple[id_face_local], //input
									first_boundaries, //input
									second_boundary_empty, //input
									sourse,
									element->self_direction, //output
									solver_grid_faces[id_face_local], //output
									Solution_faces[id_face_local][id_bf], //output
									Base_stiffness_matrix //output
								);

								if (element->GetSelfGlobalId() == 0)
								{
									printf_s("Print the mech result into .dat file... \n");

									FILE* fout_tech;
									char name_u_tech[5000];
									sprintf_s(name_u_tech, "%s/Face_%d(global%d)_linear_bf%d(%d).dat", element->self_direction, id_face_local, id_face, id_bf, target_vertex_global);
									fopen_s(&fout_tech, name_u_tech, "w");
									char name_in_file[1000];
									sprintf_s(name_in_file, "Face_%d(global%d)_linear_bf%d(%d)", id_face_local, id_face, id_bf, target_vertex_global);
									std::vector<std::vector<char>> name_value(4);
									char name_v_tmp[6][100];
									sprintf_s(name_v_tmp[0], "Ux");
									sprintf_s(name_v_tmp[1], "Uy");
									sprintf_s(name_v_tmp[2], "Uz");
									sprintf_s(name_v_tmp[3], "Material");
									for (int i = 0; i < name_value.size(); i++)
									{
										name_value[i].resize(100);
										for (int j = 0; j < name_value[i].size(); j++)
										{
											name_value[i][j] = name_v_tmp[i][j];
										}
									}
									std::vector<std::vector<double>> value(4);
									value[0].resize(solver_grid_faces[id_face_local].GetVertexCount());
									value[1].resize(solver_grid_faces[id_face_local].GetVertexCount());
									value[2].resize(solver_grid_faces[id_face_local].GetVertexCount());
									value[3].resize(solver_grid_faces[id_face_local].GetElementsCount());

									for (int i = 0; i < solver_grid_faces[id_face_local].GetVertexCount(); i++)
									{
										value[0][i] = Solution_faces[id_face_local][id_bf][i].x;
										value[1][i] = Solution_faces[id_face_local][id_bf][i].y;
										value[2][i] = Solution_faces[id_face_local][id_bf][i].z;
									}
									for (int i = 0; i < solver_grid_faces[id_face_local].GetElementsCount(); i++)
									{
										auto element = solver_grid_faces[id_face_local].GetElement(i);
										value[3][i] = element->GetIdDomain();

									}
									solver_grid_faces[id_face_local].printTecPlot3D_DiffDomains(fout_tech, value, name_value, name_in_file);
									fclose(fout_tech);
								}

								////Выравниевание на х (костыль!)
								//for (int i = 0; i < Solution_faces[id_face_local][id_bf].size(); i++)
								//{
								//	Solution_faces[id_face_local][id_bf][i].y = Solution_faces[id_face_local][id_bf][i].x;
								//	Solution_faces[id_face_local][id_bf][i].z = Solution_faces[id_face_local][id_bf][i].x;
								//}
							}
						}

						//solve volume problems
						if (true)
						{
							CSSD_Matrix<Tensor2Rank3D, Point<double>> Base_stiffness_matrix;
							Base_stiffness_matrix.print_logs = this->is_print_logFile;

							for (int i = 0; i < this->GetDomainsCount(); i++)
								element->self_grid.AddDomain(*this->GetDomain(i));

							//create Boundary values and vertexes
							printf_s("create Boundary values and vertexes for volume\n");
							struct _Dirichlet {
								int global_id, in_elem_id;
								std::vector<int> id_vertexes;
								std::vector<Point<double>> xyz;
								std::function<Point<double>(Point<bool>&, int id_vertex)> value;
							};
							std::vector<_Dirichlet> first_boundaries(element->GetLowerElementCount());
							for (int id_face_local = 0; id_face_local < first_boundaries.size(); id_face_local++)
							{
								first_boundaries[id_face_local].global_id = element->GetLowerElement(id_face_local)->GetSelfGlobalId();
								first_boundaries[id_face_local].in_elem_id = id_face_local;

								for (int i = 0; i < faces_grid_simple[id_face_local].xyz.size(); i++)
								{
									int global_id_vertex = faces_grid_simple[id_face_local].TransferIdSelfIntoGlobal(i);
									first_boundaries[id_face_local].id_vertexes.push_back(global_id_vertex);
									first_boundaries[id_face_local].xyz.push_back(base_grid_simple.xyz[global_id_vertex]);
								}
							}
							struct _Neumann {
								/*double value;
								Point<double> vector;*/
								std::function<Point<double>(Point<double>)> value;
								std::vector<std::function<Point<double>(Point<double>)>> values;
								std::vector<std::vector<int>> id_vertexes_as_triangle;
								std::vector<int> id_base_element;
							};
							std::vector<_Neumann> second_boundary_empty;

							//linear problems
							for (int id_bf = 0; id_bf < element->GetNodesCount(); id_bf++)
							{
								printf_s("create linear function for volume %d - bf linear %d\n", id_element, id_bf);

								int target_vertex_id_global = element->GetIdNode(id_bf);

								for (int i = 0; i < first_boundaries.size(); i++)
								{
									first_boundaries[i].value = [](Point<bool>& arg, int id) -> Point<double>
									{
										arg.x = true;
										arg.y = true;
										arg.z = true;
										return Point<double>(0, 0, 0);
									};
								}

								//find the target faces and create boundary conditions
								for (int id_face_local = 0; id_face_local < element->GetLowerElementCount(); id_face_local++)
								{
									for (int i = 0; i < element->GetLowerElement(id_face_local)->GetLowerElementCount(); i++)
									{
										if (target_id_for_BF[id_face_local][i] == target_vertex_id_global)
										{
											target_vertex_id_global *= 1;
											first_boundaries[id_face_local].value = [&Solution_faces, id_face_local, i, &faces_grid_simple, &element, id_bf](Point<bool>& arg, int id_vertex)->Point <double>
											{
												arg.x = true; arg.y = true; arg.z = true;
												Point<double> result = Solution_faces[id_face_local][i][faces_grid_simple[id_face_local].TransferIdGlobalIntoSelf(id_vertex)];

												//{
												//	auto f = [&](double x, double xmax, double xmin) -> double*
												//	{
												//		double phi[3];
												//		double x_loc = (x - xmin) / (xmax - xmin);
												//		phi[0] = (1 - x_loc);
												//		phi[1] = (x_loc);
												//		phi[2] = -4 * x_loc * x_loc + 4 * x_loc;

												//		return phi;
												//	};

												//	//int numeration[3][27] = {
												//	//{ 0,1,0,1,0,1,0,1, 2,0,0,1,1,2,0,1,2,0,1,2, 2,2,0,2,1,2, 2 },
												//	//{ 0,0,1,1,0,0,1,1, 0,2,0,2,0,1,1,1,0,2,2,1, 2,2,2,1,2,0, 2 },
												//	//{ 0,0,0,0,1,1,1,1, 0,0,2,0,2,0,2,2,1,1,1,1, 0,1,2,2,2,2, 2 } };
												//	int numeration[3][27] = {
												//		{ 0,0,0,0,1,1,1,1, 2,0,0,1,1,2,0,1,2,0,1,2, 2,2,0,2,1,2, 2 },
												//		{ 0,0,1,1,0,0,1,1, 0,2,0,2,0,1,1,1,0,2,2,1, 2,2,2,1,2,0, 2 },
												//		{ 0,1,0,1,0,1,0,1, 0,0,2,0,2,0,2,2,1,1,1,1, 0,1,2,2,2,2, 2 } };

												//	Point<double> X = faces_grid_simple[id_face_local].xyz[faces_grid_simple[id_face_local].TransferIdGlobalIntoSelf(id_vertex)];

												//	result.z = f(X.x, element->GetNode(7).x, element->GetNode(0).x)[numeration[0][id_bf]] *
												//		f(X.y, element->GetNode(7).y, element->GetNode(0).y)[numeration[1][id_bf]] *
												//		f(X.z, element->GetNode(7).z, element->GetNode(0).z)[numeration[2][id_bf]];
												//	result.x = result.z;
												//	result.y = result.z;
												//}

												return result;
											};
											break;
										}
									}
								}

								//right_side
								std::function<Point<double>(bool& is_solve, int id_element, Point<double> X)> sourse = [](bool& is_solve, int id_element, Point<double> X) -> Point<double>
								{
									is_solve = false;
									return Point<double>(0, 0, 0);
								};

								FEM::FEM_3D_forElasticDeformation(
									false,
									5.e-15,
									base_grid_simple, //input
									first_boundaries, //input
									second_boundary_empty, //input
									sourse,
									element->self_direction, //output
									element->self_grid, //output
									element->self_basis_functions[id_bf], //output
									Base_stiffness_matrix //output
								);

								if (this->is_print_logFile) {
									printf_s("Print the mech result into .dat file... ");

									FILE* fout_tech;
									char name_u_tech[1000];
									sprintf_s(name_u_tech, "%s/Volume_linear_bf%d.dat", element->self_direction, id_bf);
									fopen_s(&fout_tech, name_u_tech, "w");
									char name_in_file[1000];
									sprintf_s(name_in_file, "Volume_linear_bf%d", id_bf);
									std::vector<std::vector<char>> name_value(4);
									char name_v_tmp[6][100];
									sprintf_s(name_v_tmp[0], "Ux");
									sprintf_s(name_v_tmp[1], "Uy");
									sprintf_s(name_v_tmp[2], "Uz");
									sprintf_s(name_v_tmp[3], "Material");
									for (int i = 0; i < name_value.size(); i++)
									{
										name_value[i].resize(100);
										for (int j = 0; j < name_value[i].size(); j++)
										{
											name_value[i][j] = name_v_tmp[i][j];
										}
									}
									std::vector<std::vector<double>> value(4);
									value[0].resize(element->self_grid.GetVertexCount());
									value[1].resize(element->self_grid.GetVertexCount());
									value[2].resize(element->self_grid.GetVertexCount());
									value[3].resize(element->self_grid.GetElementsCount());

									for (int i = 0; i < element->self_grid.GetVertexCount(); i++)
									{
										value[0][i] = element->self_basis_functions[id_bf][i].x;
										value[1][i] = element->self_basis_functions[id_bf][i].y;
										value[2][i] = element->self_basis_functions[id_bf][i].z;
									}
									for (int i = 0; i < element->self_grid.GetElementsCount(); i++)
									{
										auto el = element->self_grid.GetElement(i);
										value[3][i] = el->GetIdDomain();

									}
									element->self_grid.printTecPlot3D_DiffDomains(fout_tech, value, name_value, name_in_file);
									fclose(fout_tech);
								}
							}
						}

						//save solutions
						for (int id_bf = 0; id_bf < element->self_basis_functions.size(); id_bf++)
						{
							FILE* fout;
							char name_out[1000];
							sprintf_s(name_out, "%s/BasisFunction_%d.txt", element->self_direction, id_bf);
							fopen_s(&fout, name_out, "w");
							for (int i = 0; i < element->self_basis_functions[id_bf].size(); i++)
							{
								fprintf_s(fout, "%.15e %.15e %.15e\n", element->self_basis_functions[id_bf][i].x, element->self_basis_functions[id_bf][i].y, element->self_basis_functions[id_bf][i].z);
							}
							fclose(fout);

							//output into techplot
							if (/*element->GetSelfGlobalId() == 0 || element->GetSelfGlobalId() == 1*/id_bf == 0) {
								FILE* fout_tech;
								char name_u_tech[5000];
								sprintf_s(name_u_tech, "%s/Techplot_BasisFunction_%d(%d).dat", element->self_direction, id_bf, element->GetIdNode(id_bf));
								fopen_s(&fout_tech, name_u_tech, "w");
								char name_in_file[1000];
								sprintf_s(name_in_file, "BasisFunction_%d(%d)", id_bf, element->GetIdNode(id_bf));
								std::vector<std::vector<char>> name_value(4);
								char name_v_tmp[4][100];
								sprintf_s(name_v_tmp[0], "Ux");
								sprintf_s(name_v_tmp[1], "Uy");
								sprintf_s(name_v_tmp[2], "Uz");
								sprintf_s(name_v_tmp[3], "Material");
								for (int i = 0; i < name_value.size(); i++)
								{
									name_value[i].resize(100);
									for (int j = 0; j < name_value[i].size(); j++)
									{
										name_value[i][j] = name_v_tmp[i][j];
									}
								}
								std::vector<std::vector<double>> value(name_value.size());
								value[0].resize(element->self_grid.GetVertexCount());
								value[1].resize(element->self_grid.GetVertexCount());
								value[2].resize(element->self_grid.GetVertexCount());
								value[3].resize(element->self_grid.GetElementsCount());
								for (int i = 0; i < element->self_grid.GetVertexCount(); i++)
								{
									value[0][i] = element->self_basis_functions[id_bf][i].x;
									value[1][i] = element->self_basis_functions[id_bf][i].y;
									value[2][i] = element->self_basis_functions[id_bf][i].z;
								}
								for (int i = 0; i < element->self_grid.GetElementsCount(); i++)
								{
									value[3][i] = element->self_grid.GetElement(i)->GetIdDomain();
								}

								element->self_grid.printTecPlot3D_DiffDomains(fout_tech, value, name_value, name_in_file);
								fclose(fout_tech);
							}
						}
					}

					//create MACRO basis functions
					for (int _id_bf = 0; _id_bf < element->self_basis_functions.size(); _id_bf++)
					{
						std::function< Point<double>(Point<double> X)> bf = [_id_bf, element](Point<double> X) -> Point<double>
						{
							Point<double> result;

							double len;
							int id_micro_elem = element->self_grid.GetNearestElementID(X, len);
							result = element->self_grid.GetSolutionInPoint(id_micro_elem, X, element->self_basis_functions[_id_bf]);

							/*if (_id_bf <= 8)
							{
								result = element->self_grid.GetSolutionInPoint(id_micro_elem, X, element->self_basis_functions[_id_bf]);
							}
							else {
								auto f = [&](double x, double xmax, double xmin) -> double*
								{
									double phi[3];
									double x_loc = (x - xmin) / (xmax - xmin);
									phi[0] = (1 - x_loc);
									phi[1] = (x_loc);
									phi[2] = -4 * x_loc * x_loc + 4 * x_loc;

									return phi;
								};

								//int numeration[3][27] = {
								//{ 0,1,0,1,0,1,0,1, 2,0,0,1,1,2,0,1,2,0,1,2, 2,2,0,2,1,2, 2 },
								//{ 0,0,1,1,0,0,1,1, 0,2,0,2,0,1,1,1,0,2,2,1, 2,2,2,1,2,0, 2 },
								//{ 0,0,0,0,1,1,1,1, 0,0,2,0,2,0,2,2,1,1,1,1, 0,1,2,2,2,2, 2 } };
								int numeration[3][27] = {
									{ 0,0,0,0,1,1,1,1, 2,0,0,1,1,2,0,1,2,0,1,2, 2,2,0,2,1,2, 2 },
									{ 0,0,1,1,0,0,1,1, 0,2,0,2,0,1,1,1,0,2,2,1, 2,2,2,1,2,0, 2 },
									{ 0,1,0,1,0,1,0,1, 0,0,2,0,2,0,2,2,1,1,1,1, 0,1,2,2,2,2, 2 } };

								result.z = f(X.x, element->GetNode(7).x, element->GetNode(0).x)[numeration[0][_id_bf]] *
									f(X.y, element->GetNode(7).y, element->GetNode(0).y)[numeration[1][_id_bf]] *
									f(X.z, element->GetNode(7).z, element->GetNode(0).z)[numeration[2][_id_bf]];
							}

							result.y = result.z;
							result.x = result.z;*/

							return result;
						};
						std::function< Point<Point<double>>(Point<double> X)> derivative_bf = [element, _id_bf](Point<double> X) -> Point<Point<double>>
						{
							Point<Point<double>> result;

							double len;
							int id_micro_elem = element->self_grid.GetNearestElementID(X, len);
							result = element->self_grid.GetDerevativeFromSolutionInPoint(id_micro_elem, X, element->self_basis_functions[_id_bf]);

							/*if (_id_bf <= 8)
							{
								result = element->self_grid.GetDerevativeFromSolutionInPoint(id_micro_elem, X, element->self_basis_functions[_id_bf]);
							}
							else {
								auto f = [&](double x, double xmax, double xmin) -> double*
								{
									double phi[3];
									double x_loc = (x - xmin) / (xmax - xmin);
									phi[0] = (1 - x_loc);
									phi[1] = (x_loc);
									phi[2] = -4 * x_loc * x_loc + 4 * x_loc;

									return phi;
								};
								auto df = [&](double x, double xmax, double xmin) -> double*
								{
									double phi[3];
									double dx_loc = 1 / (xmax - xmin);
									double x_loc = (x - xmin) / (xmax - xmin);
									phi[0] = -dx_loc;
									phi[1] = +dx_loc;
									phi[2] = -8 * x_loc * dx_loc + 4;
									return phi;
								};

								//int numeration[3][27] = {
								//	{ 0,1,0,1,0,1,0,1, 2,0,0,1,1,2,0,1,2,0,1,2, 2,2,0,2,1,2, 2 },
								//	{ 0,0,1,1,0,0,1,1, 0,2,0,2,0,1,1,1,0,2,2,1, 2,2,2,1,2,0, 2 },
								//	{ 0,0,0,0,1,1,1,1, 0,0,2,0,2,0,2,2,1,1,1,1, 0,1,2,2,2,2, 2 } };
								int numeration[3][27] = {
									{ 0,0,0,0,1,1,1,1, 2,0,0,1,1,2,0,1,2,0,1,2, 2,2,0,2,1,2, 2 },
									{ 0,0,1,1,0,0,1,1, 0,2,0,2,0,1,1,1,0,2,2,1, 2,2,2,1,2,0, 2 },
									{ 0,1,0,1,0,1,0,1, 0,0,2,0,2,0,2,2,1,1,1,1, 0,1,2,2,2,2, 2 } };

								double resx, resy, resz;
								resx = df(X.x, element->GetNode(7).x, element->GetNode(0).x)[numeration[0][_id_bf]];
								resy = f(X.y, element->GetNode(7).y, element->GetNode(0).y)[numeration[1][_id_bf]];
								resz = f(X.z, element->GetNode(7).z, element->GetNode(0).z)[numeration[2][_id_bf]];
								result.z.x = resx * resy * resz;

								resx = f(X.x, element->GetNode(7).x, element->GetNode(0).x)[numeration[0][_id_bf]];
								resy = df(X.y, element->GetNode(7).y, element->GetNode(0).y)[numeration[1][_id_bf]];
								resz = f(X.z, element->GetNode(7).z, element->GetNode(0).z)[numeration[2][_id_bf]];
								result.z.y = resx * resy * resz;

								resx = f(X.x, element->GetNode(7).x, element->GetNode(0).x)[numeration[0][_id_bf]];
								resy = f(X.y, element->GetNode(7).y, element->GetNode(0).y)[numeration[1][_id_bf]];
								resz = df(X.z, element->GetNode(7).z, element->GetNode(0).z)[numeration[2][_id_bf]];
								result.z.z = resx * resy * resz;
							}

							result.y = result.z;
							result.x = result.z;*/

							return result;
						};

						int target_global_id;
						target_global_id = element->GetIdNode(_id_bf);
						element->SetDOF(_id_bf, target_global_id, bf, derivative_bf);
					}
				}

				for (int id_face = 0; id_face < this->faces.size(); id_face++)
				{
					this->faces[id_face].ResizeDOF(this->faces[id_face].GetLowerElementCount() + this->faces[id_face].GetNodesCount() + 1);
					for (int id_vertex_local = 0; id_vertex_local < this->faces[id_face].GetNodesCount(); id_vertex_local++)
					{
						int id_vertex_global = this->faces[id_face].GetIdNode(id_vertex_local);
						auto local_bf = this->GetElement(this->faces[id_face].GetUpperElement(0)->GetSelfGlobalId())->GetBasisFunctionInGlobalID(id_vertex_global);
						auto local_Dbf = this->GetElement(this->faces[id_face].GetUpperElement(0)->GetSelfGlobalId())->GetDerivativeOfBasisFunctionInGlobalID(id_vertex_global);
						this->faces[id_face].SetDOF(id_vertex_local, id_vertex_global, *local_bf, *local_Dbf);
					}
				}
			}
			catch (const char* err)
			{
				std::cerr << "Error (Grid_forMech_Poly_Order1.CreateFunctionalSpaces): " << err << std::endl;
			}
		}
		template <typename Dirichlet, typename Neumann>
		void CreateBoundaryConditions(std::vector<Dirichlet>& first_boundary, std::vector<Neumann>& second_boundary)
		{
			//1st boundary
			for (int id_type = 0; id_type < first_boundary.size(); id_type++)
			{
				for (int i = 0; i < first_boundary[id_type].id_vertexes.size(); i++)
				{
					int global_dof_id = first_boundary[id_type].id_vertexes[i];

					MsFEM::BoundaryVertex_forMech_Poly_Order1 temp;
					temp.boundary_value = first_boundary[id_type].value;

					temp.ResizeDOF(1);

					temp.SetIdNode(0, global_dof_id);
					temp.SetNode(0, this->GetPtrCoordinateViaID(global_dof_id));
					int target_edge_id = this->vertexes[global_dof_id].GetUpperElement(0)->GetSelfGlobalId();
					int target_face_id = this->edges[target_edge_id].GetUpperElement(0)->GetSelfGlobalId();
					int target_elem_id = this->faces[target_face_id].GetUpperElement(0)->GetSelfGlobalId();
					temp.SetDOF(0, global_dof_id, *this->GetElement(target_elem_id)->GetBasisFunctionInGlobalID(global_dof_id), *this->GetElement(target_elem_id)->GetDerivativeOfBasisFunctionInGlobalID(global_dof_id));

					this->boundary_vertexes.push_back(temp);
				}
			}
		}
		template <typename Dirichlet, typename Neumann>
		void Initialization(math::SimpleGrid& base_grid, std::vector<Dirichlet>& first_boundary, std::vector<Neumann>& second_boundary, char* fine_mesh_dir)
		{
			try {

				//create geometry space
				math::MakeCopyVector_A_into_B(base_grid.xyz, *(this->GetCoordinates()));
				this->CreateXYZline();
				this->SetElementsCount((int)base_grid.nvtr.size());
				for (int id_element = 0; id_element < this->GetElementsCount(); id_element++)
				{
					auto new_element = this->GetElement(id_element);

					std::vector<int> nodes_global(base_grid.nvtr[id_element][0]);
					std::vector<Point<double>*> P(base_grid.nvtr[id_element][0]);
					int ii = 1;
					for (int n = 0; n < nodes_global.size(); n++)
					{
						nodes_global[n] = base_grid.nvtr[id_element][ii];
						ii++;
						P[n] = this->GetPtrCoordinateViaID(nodes_global[n]);
					}
					std::vector<std::vector<int>> nodes_faces(base_grid.nvtr[id_element][ii]);
					ii++;
					for (int f = 0; f < nodes_faces.size(); f++)
					{
						nodes_faces[f].resize(base_grid.nvtr[id_element][ii]);
						ii++;
						for (int nf = 0; nf < nodes_faces[f].size(); nf++)
						{
							nodes_faces[f][nf] = base_grid.nvtr[id_element][ii];
							ii++;
						}
					}

					new_element->SetGeometry(nodes_global, nodes_faces, P);
					new_element->SetIdDomain(base_grid.nvkat[id_element]);
				}
				this->CreateQTree();

				//create topology
				CreateTopology();

				//create functional spaces
				CreateFunctionalSpaces(fine_mesh_dir);

				CreateBoundaryConditions(first_boundary, second_boundary);
			}
			catch (const std::exception&)
			{
				printf_s("Error: MultiXFEM/MultiXFEM_Grid.h/void Initialization(geometry::Grid<geometry::Tetrahedron> &base_grid)\n");
			}
		}

		template <typename Matrix>
		void CreationPortrait(Matrix& matrix)
		{
			try {
				printf_s("\n");
				std::vector<std::vector<int>> tmp_down_columns(this->GetDOFsCount()), tmp_up_columns(this->GetDOFsCount());
				std::vector<std::vector<int>> down_columns(this->GetDOFsCount()), up_columns(this->GetDOFsCount());
				for (int id_elem = 0; id_elem < this->GetElementsCount(); id_elem++)
				{
					if (id_elem % 100 == 0)
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
					if (id_string % 1 == 0)
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

		Point<double> GetSolutionInPoint(int id_element, Point<double> X, std::vector<Point<double>>& solution)
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
		Point<Point<double>> GetDerevativeFromSolutionInPoint(int id_element, Point<double> X, std::vector<Point<double>>& solution)
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

		void printTecPlot3D(/*char *directory,*/ FILE* fdat, std::vector<std::vector<double>>& value, std::vector<std::vector<char>> name_value, char* name_zone)
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

		~Grid_forMech_Poly_Order1()
		{
			DOFs_count = 0;
			std::vector<BoundaryVertex_forMech_Poly_Order1> v_b1;
			std::vector<BoundaryVertex_forMech_Poly_Order1>(v_b1).swap(this->boundary_vertexes);
			std::vector<BoundaryFace_forMech_Poly_Order1> v_b2;
			std::vector<BoundaryFace_forMech_Poly_Order1>(v_b2).swap(this->boundary_faces);
			std::vector<int> v_p;
			std::vector<int>(v_p).swap(this->accordance_DOF_and_vertex);
			std::vector<Vertex> v_v;
			std::vector<Vertex>(v_v).swap(this->vertexes);
			std::vector<Edge> v_e;
			std::vector<Edge>(v_e).swap(this->edges);
			std::vector<Face> v_f;
			std::vector<Face>(v_f).swap(this->faces);

			this->DeleteGeometriGrid();
		}
	private:

	};
}