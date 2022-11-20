#pragma once
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

namespace FEM{
	class FiniteElement_forScal :
		public geometry::Tetrahedron,
		//public topology::Tetrahedron<topology::lower::Triangle, topology::upper::EmptyElement>,
		public topology::Tetrahedron<topology::lower::Vertex, topology::upper::EmptyElement>,
		public functional::Shape<int, double>
	{
	public:

		FiniteElement_forScal() { return; };
		~FiniteElement_forScal() { return; };

		void SolveLocalMatrix(DenseMatrix<double, double> &local_matix, 
			std::function<double(int,Point<double>)> &koefStiffness, 
			std::function<double(int,Point<double>)> &koefMass, 
			std::function<double(int,Point<double>)> &F)
		{
			try {
				local_matix.SetSize(this->GetDOFsCount());
				int global_id = this->GetSelfGlobalId();

				this->SetIntegrationLaw(1);

				for (int _bf_local_id_I = 0; _bf_local_id_I < this->GetDOFsCount(); _bf_local_id_I++)
				{
					auto D_basis_function_I = this->GetDerivativeOfBasisFunctionInLocalID(_bf_local_id_I);
					auto basis_function_I = this->GetBasisFunctionInLocalID(_bf_local_id_I);
					
					std::function<double(Point<double>)> RightSide = [&global_id, &_bf_local_id_I, &basis_function_I, &F]
					(Point<double> X) -> double
					{
						double result;
						auto right_side = F(global_id, X);
						double BF_I_inX = (*basis_function_I)(X);

						result = right_side * BF_I_inX;
						//result = 1.0;
						return result;
					};
					local_matix.F[_bf_local_id_I] = this->SolveIntegral(RightSide);

					for (int _bf_local_id_J = 0; _bf_local_id_J < this->GetDOFsCount(); _bf_local_id_J++)
					{
						auto D_basis_function_J = this->GetDerivativeOfBasisFunctionInLocalID(_bf_local_id_J);
						auto basis_function_J = this->GetBasisFunctionInLocalID(_bf_local_id_J);

						Point<double> o = this->GetWeightCentr();
						auto val = (*this->GetDerivativeOfBasisFunctionInLocalID(0))(o);

						std::function<double(Point<double>)> StiffnessMatrix = [&global_id, &_bf_local_id_I, &_bf_local_id_J, &koefStiffness, &D_basis_function_I, &D_basis_function_J]
						(Point<double> X) -> double
						{
							double result;
							auto S = koefStiffness(global_id, X);
							Point<double> derevativeBF_I_inX = (*D_basis_function_I)(X);
							Point<double> derevativeBF_J_inX = (*D_basis_function_J)(X);

							result = S * (
								derevativeBF_I_inX.x*derevativeBF_J_inX.x +
								derevativeBF_I_inX.y*derevativeBF_J_inX.y +
								derevativeBF_I_inX.z*derevativeBF_J_inX.z);
							//result = 1.0;
							return result;
						};
						std::function<double(Point<double>)> MassMatrix = [&global_id, &_bf_local_id_I, &_bf_local_id_J, &basis_function_I, &basis_function_J, &koefMass]
						(Point<double> X) -> double
						{
							double result;
							auto M = koefMass(global_id, X);
							double BF_I_inX = (*basis_function_I)(X);
							double BF_J_inX = (*basis_function_J)(X);

							result = M * BF_I_inX * BF_J_inX;
							//result = 1.0;
							return result;
						};
						

						double V = this->GetVolume();
						double S = this->SolveIntegral(StiffnessMatrix);
						double M = this->SolveIntegral(MassMatrix);
						local_matix.A[_bf_local_id_I][_bf_local_id_J] = S + M;
					}
				}
			}
			catch (const std::exception&)
			{
				printf_s("Error: XFEM/MultiXFEM_Grid.h/void SolveLocalMatrix(DenseMatrix<Tensor2Rank3D> &local_matix, std::function<std::vector<std::vector<double>>(Point)> &koefD)\n");
			}
		}
		void SolveLocalMatrix(DenseMatrix<double, double>& local_matix,
			std::function<double(int, Point<double>)>& koefStiffness)
		{
			try {
				local_matix.SetSize(this->GetDOFsCount());
				int global_id = this->GetSelfGlobalId();

				this->SetIntegrationLaw(1);

				for (int _bf_local_id_I = 0; _bf_local_id_I < this->GetDOFsCount(); _bf_local_id_I++)
				{
					auto D_basis_function_I = this->GetDerivativeOfBasisFunctionInLocalID(_bf_local_id_I);
					auto basis_function_I = this->GetBasisFunctionInLocalID(_bf_local_id_I);
										
					for (int _bf_local_id_J = 0; _bf_local_id_J < this->GetDOFsCount(); _bf_local_id_J++)
					{
						auto D_basis_function_J = this->GetDerivativeOfBasisFunctionInLocalID(_bf_local_id_J);
						auto basis_function_J = this->GetBasisFunctionInLocalID(_bf_local_id_J);

						Point<double> o = this->GetWeightCentr();
						auto val = (*this->GetDerivativeOfBasisFunctionInLocalID(0))(o);

						std::function<double(Point<double>)> StiffnessMatrix = [&global_id, &_bf_local_id_I, &_bf_local_id_J, &koefStiffness, &D_basis_function_I, &D_basis_function_J]
						(Point<double> X) -> double
						{
							double result;
							auto S = koefStiffness(global_id, X);
							Point<double> derevativeBF_I_inX = (*D_basis_function_I)(X);
							Point<double> derevativeBF_J_inX = (*D_basis_function_J)(X);

							result = S * (
								derevativeBF_I_inX.x * derevativeBF_J_inX.x +
								derevativeBF_I_inX.y * derevativeBF_J_inX.y +
								derevativeBF_I_inX.z * derevativeBF_J_inX.z);
							//result = 1.0;
							return result;
						};


						double V = this->GetVolume();
						double S = this->SolveIntegral(StiffnessMatrix);
						local_matix.A[_bf_local_id_I][_bf_local_id_J] = S;
					}
				}
			}
			catch (const std::exception&)
			{
				printf_s("Error: XFEM/MultiXFEM_Grid.h/void SolveLocalMatrix(DenseMatrix<Tensor2Rank3D> &local_matix, std::function<std::vector<std::vector<double>>(Point)> &koefD)\n");
			}
		}
		void SolveLocalMatrix(DenseMatrix<double, double>& local_matix,
			std::function<Tensor2Rank3D(int, Point<double>)>& koefStiffness)
		{
			try {
				local_matix.SetSize(this->GetDOFsCount());
				int global_id = this->GetSelfGlobalId();

				this->SetIntegrationLaw(1);

				for (int _bf_local_id_I = 0; _bf_local_id_I < this->GetDOFsCount(); _bf_local_id_I++)
				{
					auto D_basis_function_I = this->GetDerivativeOfBasisFunctionInLocalID(_bf_local_id_I);
					auto basis_function_I = this->GetBasisFunctionInLocalID(_bf_local_id_I);

					for (int _bf_local_id_J = 0; _bf_local_id_J < this->GetDOFsCount(); _bf_local_id_J++)
					{
						auto D_basis_function_J = this->GetDerivativeOfBasisFunctionInLocalID(_bf_local_id_J);
						auto basis_function_J = this->GetBasisFunctionInLocalID(_bf_local_id_J);

						Point<double> o = this->GetWeightCentr();
						auto val = (*this->GetDerivativeOfBasisFunctionInLocalID(0))(o);

						std::function<double(Point<double>)> StiffnessMatrix = [&global_id, &_bf_local_id_I, &_bf_local_id_J, &koefStiffness, &D_basis_function_I, &D_basis_function_J]
						(Point<double> X) -> double
						{
							double result;
							auto S = koefStiffness(global_id, X);
							Point<double> derevativeBF_I_inX = (*D_basis_function_I)(X);
							Point<double> derevativeBF_J_inX = (*D_basis_function_J)(X);

							//auto S_gradU = S * derevativeBF_J_inX;

							result = derevativeBF_I_inX * (S * derevativeBF_J_inX);
							//result = 1.0;
							return result;
						};


						double V = this->GetVolume();
						double S = this->SolveIntegral(StiffnessMatrix);
						local_matix.A[_bf_local_id_I][_bf_local_id_J] = S;
					}
				}
			}
			catch (const std::exception&)
			{
				printf_s("Error: XFEM/MultiXFEM_Grid.h/void SolveLocalMatrix(DenseMatrix<Tensor2Rank3D> &local_matix, std::function<std::vector<std::vector<double>>(Point)> &koefD)\n");
			}
		}
	};
	class BoundaryVertex_forScal :
		public geometry::Vertex,
		//public topology::Vertex<topology::lower::EmptyElement, topology::upper::Segment>,
		public topology::Vertex<topology::lower::EmptyElement, topology::upper::Tetrahedron>,
		public functional::Shape<int, double>
	{
	public:
		std::function<double(Point<double> &X)> boundary_value;

		BoundaryVertex_forScal() { return; };
		~BoundaryVertex_forScal() { return; };

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
	class BoundaryFace_forScal :
		public geometry::Triangle,
		public topology::Triangle<topology::lower::Vertex, topology::upper::Tetrahedron>,
		public functional::Shape<int, double>
	{
	public:
		std::function<double(Point<double> &X)> boundary_value;

		BoundaryFace_forScal() { return; };
		~BoundaryFace_forScal() { return; };

		void SetTypeBoundary(int type)
		{
			this->type_boundary = type;
		}
		int GetTypeBoundary()
		{
			return this->type_boundary;
		}

		void SolveLocalBoundaryVector(std::vector<double> &local_vector, std::function<double(Point<double> &X)> &boundary_value)
		{
			try {
				/*local_vector.resize(this->GetDOFsCount());
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

						return result;
					};

					double V = this->GetVolume();
					local_vector[_bf_local_id_I] = this->SolveIntegral(ValueInPoint);

				}*/
			}
			catch (const std::exception&)
			{
				printf_s("Error: XFEM/MultiXFEM_Grid.h/void SolveLocalMatrix(DenseMatrix<Tensor2Rank3D> &local_matix, std::function<std::vector<std::vector<double>>(Point)> &koefD)\n");
			}
		}

	private:
		int type_boundary;
	};
	class Grid_forScal : public geometry::Grid<FiniteElement_forScal>
	{
		int DOFs_count;
		std::vector<topology::Vertex<topology::lower::EmptyElement, topology::upper::Tetrahedron>> vertexes; //Vertexes

	public:
		std::vector<geometry::Crack> cracks;
		std::vector<BoundaryVertex_forScal> boundary_vertexes;
		std::vector<BoundaryFace_forScal> boundary_faces;

		std::vector<int> accordance_DOF_and_vertex;

		Grid_forScal()
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
						if (A[j] != A[i]) //םו ןמגעמנ
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
					std::function< double(Point<double> X)> bf = [element, _id_vertex](Point<double> X) -> double
					{
						double result;
						result = element->alpha[_id_vertex][0] + element->alpha[_id_vertex][1] * X.x + element->alpha[_id_vertex][2] * X.y + element->alpha[_id_vertex][3] * X.z;
						return result;
					};
					std::function< Point<double>(Point<double> X)> derivative_bf = [element, _id_vertex](Point<double> X) -> Point<double>
					{
						Point<double> result;

						result.x = element->alpha[_id_vertex][1];
						result.y = element->alpha[_id_vertex][2];
						result.z = element->alpha[_id_vertex][3];

						return result;
					};
					element->SetDOF(_id_vertex, element->GetIdNode(_id_vertex), bf, derivative_bf);
				}
			}

			//update boundary conditions
			for (int id_vertex = 0; id_vertex < this->boundary_vertexes.size(); id_vertex++)
			{
				auto boundary = &(this->boundary_vertexes[id_vertex]);
				std::function< double(Point<double> X)> bf = [boundary](Point<double> X) -> double
				{
					return 0;
				};
				std::function< Point<double>(Point<double> X)> derivative_bf = [boundary](Point<double> X) -> Point<double>
				{
					Point<double> result;
					return result;
				};

				this->boundary_vertexes[id_vertex].ResizeDOF(1);
				this->boundary_vertexes[id_vertex].SetDOF(0, this->boundary_vertexes[id_vertex].GetIdNode(0), bf, derivative_bf);
			}
			//for (int id_face = 0; id_face < this->boundary_faces.size(); id_face++)
			//{
			//	auto boundary = &(this->boundary_faces[id_face]);

			//	int id_base_tetr = boundary->GetUpperElement(0)->GetSelfGlobalId();
			//	auto base_tetr = this->GetElement(id_base_tetr);
			//	boundary->ResizeDOF(boundary->GetNodesCount());
			//	for (int id_vertex = 0; id_vertex < boundary->GetNodesCount(); id_vertex++)
			//	{
			//		//int id_base_tetr = this->vertexes[boundary->GetIdNode(id_vertex)].GetUpperElement(0)->GetSelfGlobalId();
			//		//auto base_tetr = this->GetElement(id_base_tetr);
			//		int local_id_in_tetr = base_tetr->GetLocalID_forVertexID(boundary->GetIdNode(id_vertex));
			//		std::function< Point<double>(Point<double> X)> bf = *(base_tetr->GetBasisFunctionInLocalID(local_id_in_tetr));
			//		std::function< Point<Point<double>>(Point<double> X)> derivative_bf = *(base_tetr->GetDerivativeOfBasisFunctionInLocalID(local_id_in_tetr));
			//		boundary->SetDOF(id_vertex, boundary->GetIdNode(id_vertex), bf, derivative_bf);
			//	}
			//}
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
				/*{
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

							this->boundary_faces[id].boundary_value = second_boundary[id_type].value;
							for (int id_vertex = 0; id_vertex < this->boundary_faces[id].GetNodesCount(); id_vertex++)
							{
								this->boundary_faces[id].SetIdNode(id_vertex, second_boundary[id_type].id_vertexes_as_triangle[id_face][id_vertex]);
								this->boundary_faces[id].SetNode(id_vertex, this->GetPtrCoordinateViaID(second_boundary[id_type].id_vertexes_as_triangle[id_face][id_vertex]));
							}
							id++;
						}
					}
				}*/

				CreateFunctionalSpaces();

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
					if (id_elem % 1000 == 0)
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

		double GetSolutionInPoint(int id_element, Point<double> X, std::vector<double> &solution)
		{
			auto element = this->GetElement(id_element);
			double result=0;
			for (int j = 0; j < element->GetDOFsCount(); j++)
			{
				auto bf = (*element->GetBasisFunctionInLocalID(j))(X);
				auto value = solution[element->GetDOFInLocalID(j)];
				result += bf * value;
			}
			return result;
		}
		Point<double> GetDerevativeFromSolutionInPoint(int id_element, Point<double> X, std::vector<double> &solution)
		{
			auto element = this->GetElement(id_element);
			Point<double> result;
			for (int j = 0; j < element->GetDOFsCount(); j++)
			{
				auto bf = (*element->GetDerivativeOfBasisFunctionInLocalID(j))(X);
				auto value = solution[element->GetDOFInLocalID(j)];
				result += bf * value;
			}
			return result;
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
		void printTecPlot3D_DiffDomains(/*char *directory,*/ FILE* fdat, std::vector<std::vector<double>>& value, std::vector<std::vector<char>> name_value, char* name_zone)
		{
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

			std::vector<int> in_domain(this->GetDomainsCount());
			for (int i = 0; i < this->GetElementsCount(); i++)
			{
				in_domain[this->GetElement(i)->GetIdDomain()]++;
			}

			for (int id_domain = 0; id_domain < this->GetDomainsCount(); id_domain++)
			{
				if (in_domain[id_domain] != 0)
				{
					fprintf_s(fdat, "ZONE T=\"%s_%d\"\n", name_zone, id_domain);
					fprintf_s(fdat, " N=%d,  E=%d, F=FEBLOCK ET=Tetrahedron \n", this->GetVertexCount(), in_domain[id_domain]);
					fprintf_s(fdat, " VARLOCATION=(NODAL NODAL NODAL");
					for (int i = 0; i < value.size(); i++)
					{
						if (value[i].size() == this->GetElementsCount()) fprintf_s(fdat, " CELLCENTERED");
						else fprintf_s(fdat, " NODAL");
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
						if (value[i].size() == this->GetElementsCount())
						{
							for (int j = 0; j < this->GetElementsCount(); j++)
								if (this->GetElement(j)->GetIdDomain() == id_domain)
									fprintf_s(fdat, "%.10e\n", value[i][j]);
						}
						else {
							for (int j = 0; j < value[i].size(); j++)
								fprintf_s(fdat, "%.10e\n", value[i][j]);
						}
						fprintf_s(fdat, "\n");
					}

					for (int i = 0; i < this->GetElementsCount(); i++)
					{
						if (this->GetElement(i)->GetIdDomain() == id_domain)
						{
							for (int j = 0; j < this->GetElement(i)->GetNodesCount(); j++)
								fprintf_s(fdat, "%d ", this->GetElement(i)->GetIdNode(j) + 1);
						}
						fprintf_s(fdat, "\n");
					}
				}
			}
			fclose(fdat);
		}


		~Grid_forScal()
		{
			DOFs_count = 0;
			std::vector<geometry::Crack> v_cr;
			std::vector<geometry::Crack>(v_cr).swap(this->cracks);
			std::vector<BoundaryVertex_forScal> v_b1;
			std::vector<BoundaryVertex_forScal>(v_b1).swap(this->boundary_vertexes);
			std::vector<BoundaryFace_forScal> v_b2;
			std::vector<BoundaryFace_forScal>(v_b2).swap(this->boundary_faces);
			std::vector<int> v_p;
			std::vector<int>(v_p).swap(this->accordance_DOF_and_vertex);
			std::vector<topology::Vertex<topology::lower::EmptyElement, topology::upper::Tetrahedron>> v_v;
			std::vector<topology::Vertex<topology::lower::EmptyElement, topology::upper::Tetrahedron>>(v_v).swap(this->vertexes);

			this->DeleteGeometriGrid();
		}
	private:

	};

	class FiniteElement_forScal_OrderBF2 :
		public geometry::Tetrahedron,
		//public topology::Tetrahedron<topology::lower::Triangle, topology::upper::EmptyElement>,
		public topology::Tetrahedron<topology::lower::Segment, topology::upper::EmptyElement>,
		public functional::Shape<int, double>
	{
	public:

		FiniteElement_forScal_OrderBF2() { return; };
		~FiniteElement_forScal_OrderBF2() { return; };

		void SolveLocalMatrix(DenseMatrix<double, double>& local_matix,
			std::function<double(int, Point<double>)>& koefStiffness,
			std::function<double(int, Point<double>)>& koefMass,
			std::function<double(int, Point<double>)>& F)
		{
			try {
				local_matix.SetSize(this->GetDOFsCount());
				int global_id = this->GetSelfGlobalId();

				this->SetIntegrationLaw(4);

				for (int _bf_local_id_I = 0; _bf_local_id_I < this->GetDOFsCount(); _bf_local_id_I++)
				{
					auto D_basis_function_I = this->GetDerivativeOfBasisFunctionInLocalID(_bf_local_id_I);
					auto basis_function_I = this->GetBasisFunctionInLocalID(_bf_local_id_I);

					std::function<double(Point<double>)> RightSide = [&global_id, &_bf_local_id_I, &basis_function_I, &F]
					(Point<double> X) -> double
					{
						double result;
						auto right_side = F(global_id, X);
						double BF_I_inX = (*basis_function_I)(X);

						result = right_side * BF_I_inX;
						//result = 1.0;
						return result;
					};
					local_matix.F[_bf_local_id_I] = this->SolveIntegral(RightSide);

					for (int _bf_local_id_J = _bf_local_id_I; _bf_local_id_J < this->GetDOFsCount(); _bf_local_id_J++)
					{
						auto D_basis_function_J = this->GetDerivativeOfBasisFunctionInLocalID(_bf_local_id_J);
						auto basis_function_J = this->GetBasisFunctionInLocalID(_bf_local_id_J);

						Point<double> o = this->GetWeightCentr();
						auto val = (*this->GetDerivativeOfBasisFunctionInLocalID(0))(o);

						std::function<double(Point<double>)> StiffnessMatrix = [&global_id, &_bf_local_id_I, &_bf_local_id_J, &koefStiffness, &D_basis_function_I, &D_basis_function_J]
						(Point<double> X) -> double
						{
							double result;
							auto S = koefStiffness(global_id, X);
							Point<double> derevativeBF_I_inX = (*D_basis_function_I)(X);
							Point<double> derevativeBF_J_inX = (*D_basis_function_J)(X);

							result = S * (
								derevativeBF_I_inX.x * derevativeBF_J_inX.x +
								derevativeBF_I_inX.y * derevativeBF_J_inX.y +
								derevativeBF_I_inX.z * derevativeBF_J_inX.z);
							//result = 1.0;
							return result;
						};
						std::function<double(Point<double>)> MassMatrix = [&global_id, &_bf_local_id_I, &_bf_local_id_J, &basis_function_I, &basis_function_J, &koefMass]
						(Point<double> X) -> double
						{
							double result;
							auto M = koefMass(global_id, X);
							double BF_I_inX = (*basis_function_I)(X);
							double BF_J_inX = (*basis_function_J)(X);

							result = M * BF_I_inX * BF_J_inX;
							//result = 1.0;
							return result;
						};


						double V = this->GetVolume();
						double S = this->SolveIntegral(StiffnessMatrix);
						double M = this->SolveIntegral(MassMatrix);
						local_matix.A[_bf_local_id_I][_bf_local_id_J] = S + M;
						local_matix.A[_bf_local_id_J][_bf_local_id_I] = S + M;
					}
				}
			}
			catch (const std::exception&)
			{
				printf_s("Error: FEM/FEM_Grid.h/void SolveLocalMatrix(..)\n");
			}
		}
		void SolveLocalMatrix(DenseMatrix<double, double>& local_matix,
			std::function<double(int, Point<double>)>& koefStiffness,
			std::function<double(int, Point<double>)>& F)
		{
			try {
				local_matix.SetSize(this->GetDOFsCount());
				int global_id = this->GetSelfGlobalId();

				this->SetIntegrationLaw(4);

				for (int _bf_local_id_I = 0; _bf_local_id_I < this->GetDOFsCount(); _bf_local_id_I++)
				{
					auto D_basis_function_I = this->GetDerivativeOfBasisFunctionInLocalID(_bf_local_id_I);
					auto basis_function_I = this->GetBasisFunctionInLocalID(_bf_local_id_I);

					std::function<double(Point<double>)> RightSide = [&global_id, &_bf_local_id_I, &basis_function_I, &F]
					(Point<double> X) -> double
					{
						double result;
						auto right_side = F(global_id, X);
						double BF_I_inX = (*basis_function_I)(X);

						result = right_side * BF_I_inX;
						//result = 1.0;
						return result;
					};
					local_matix.F[_bf_local_id_I] = this->SolveIntegral(RightSide);

					for (int _bf_local_id_J = _bf_local_id_I; _bf_local_id_J < this->GetDOFsCount(); _bf_local_id_J++)
					{
						auto D_basis_function_J = this->GetDerivativeOfBasisFunctionInLocalID(_bf_local_id_J);
						auto basis_function_J = this->GetBasisFunctionInLocalID(_bf_local_id_J);

						Point<double> o = this->GetWeightCentr();
						auto val = (*this->GetDerivativeOfBasisFunctionInLocalID(0))(o);

						std::function<double(Point<double>)> StiffnessMatrix = [&global_id, &_bf_local_id_I, &_bf_local_id_J, &koefStiffness, &D_basis_function_I, &D_basis_function_J]
						(Point<double> X) -> double
						{
							double result;
							auto S = koefStiffness(global_id, X);
							Point<double> derevativeBF_I_inX = (*D_basis_function_I)(X);
							Point<double> derevativeBF_J_inX = (*D_basis_function_J)(X);

							result = S * (
								derevativeBF_I_inX.x * derevativeBF_J_inX.x +
								derevativeBF_I_inX.y * derevativeBF_J_inX.y +
								derevativeBF_I_inX.z * derevativeBF_J_inX.z);
							//result = 1.0;
							return result;
						};
						
						double V = this->GetVolume();
						double S = this->SolveIntegral(StiffnessMatrix);
						local_matix.A[_bf_local_id_I][_bf_local_id_J] = S;
						local_matix.A[_bf_local_id_J][_bf_local_id_I] = S;
					}
				}
			}
			catch (const std::exception&)
			{
				printf_s("Error: FEM/FEM_Grid.h/void SolveLocalMatrix(..)\n");
			}
		}
		void SolveLocalMatrix(DenseMatrix<double, double>& local_matix,
			std::function<double(int, Point<double>)>& koefStiffness)
		{
			try {
				local_matix.SetSize(this->GetDOFsCount());
				int global_id = this->GetSelfGlobalId();

				this->SetIntegrationLaw(4);

				for (int _bf_local_id_I = 0; _bf_local_id_I < this->GetDOFsCount(); _bf_local_id_I++)
				{
					auto D_basis_function_I = this->GetDerivativeOfBasisFunctionInLocalID(_bf_local_id_I);
					auto basis_function_I = this->GetBasisFunctionInLocalID(_bf_local_id_I);

					for (int _bf_local_id_J = _bf_local_id_I; _bf_local_id_J < this->GetDOFsCount(); _bf_local_id_J++)
					{
						auto D_basis_function_J = this->GetDerivativeOfBasisFunctionInLocalID(_bf_local_id_J);
						auto basis_function_J = this->GetBasisFunctionInLocalID(_bf_local_id_J);

						Point<double> o = this->GetWeightCentr();
						auto val = (*this->GetDerivativeOfBasisFunctionInLocalID(0))(o);

						std::function<double(Point<double>)> StiffnessMatrix = [&global_id, &_bf_local_id_I, &_bf_local_id_J, &koefStiffness, &D_basis_function_I, &D_basis_function_J]
						(Point<double> X) -> double
						{
							double result;
							auto S = koefStiffness(global_id, X);
							Point<double> derevativeBF_I_inX = (*D_basis_function_I)(X);
							Point<double> derevativeBF_J_inX = (*D_basis_function_J)(X);

							result = S * (
								derevativeBF_I_inX.x * derevativeBF_J_inX.x +
								derevativeBF_I_inX.y * derevativeBF_J_inX.y +
								derevativeBF_I_inX.z * derevativeBF_J_inX.z);
							//result = 1.0;
							return result;
						};

						double V = this->GetVolume();
						double S = this->SolveIntegral(StiffnessMatrix);
						local_matix.A[_bf_local_id_I][_bf_local_id_J] = S;
						local_matix.A[_bf_local_id_J][_bf_local_id_I] = S;
					}
				}
			}
			catch (const std::exception&)
			{
				printf_s("Error: FEM/FEM_Grid.h/void SolveLocalMatrix(..)\n");
			}
		}
	};
	class BoundaryVertex_forScal_OrderBF2 :
		public geometry::Vertex,
		//public topology::Vertex<topology::lower::EmptyElement, topology::upper::Segment>,
		public topology::Vertex<topology::lower::EmptyElement, topology::upper::Tetrahedron>,
		public functional::Shape<int, double>
	{
	public:
		std::function<double(Point<double>& X)> boundary_value;

		BoundaryVertex_forScal_OrderBF2() { return; };
		~BoundaryVertex_forScal_OrderBF2() { return; };

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
	class BoundaryEdges_forScal_OrderBF2 :
		public geometry::Vertex,
		//public topology::Vertex<topology::lower::EmptyElement, topology::upper::Segment>,
		public topology::Segment<topology::lower::EmptyElement, topology::upper::Tetrahedron>,
		public functional::Shape<int, double>
	{
	public:
		std::function<double(Point<double>& X)> boundary_value;

		BoundaryEdges_forScal_OrderBF2() { return; };
		~BoundaryEdges_forScal_OrderBF2() { return; };

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
	class BoundaryFace_forScal_OrderBF2 :
		public geometry::Triangle,
		public topology::Triangle<topology::lower::Vertex, topology::upper::Tetrahedron>,
		public functional::Shape<int, double>
	{
	public:
		std::function<double(Point<double>& X)> boundary_value;

		BoundaryFace_forScal_OrderBF2() { return; };
		~BoundaryFace_forScal_OrderBF2() { return; };

		void SetTypeBoundary(int type)
		{
			this->type_boundary = type;
		}
		int GetTypeBoundary()
		{
			return this->type_boundary;
		}

		void SolveLocalBoundaryVector(std::vector<double>& local_vector, std::function<double(Point<double>& X)>& boundary_value)
		{
			try {
				/*local_vector.resize(this->GetDOFsCount());
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

						return result;
					};

					double V = this->GetVolume();
					local_vector[_bf_local_id_I] = this->SolveIntegral(ValueInPoint);

				}*/
			}
			catch (const std::exception&)
			{
				printf_s("Error: XFEM/MultiXFEM_Grid.h/void SolveLocalMatrix(DenseMatrix<Tensor2Rank3D> &local_matix, std::function<std::vector<std::vector<double>>(Point)> &koefD)\n");
			}
		}

	private:
		int type_boundary;
	};
	class Grid_forScal_OrderBF2 : public geometry::Grid<FiniteElement_forScal_OrderBF2>
	{
		int DOFs_count;
		struct Edge : //1D elements
			public topology::Segment<topology::lower::Vertex, topology::upper::Tetrahedron>,
			public geometry::Segment,
			public functional::Shape<int, Point<double>>
		{
		public:
			Point<double> geo_centr;

			Edge() {};
			~Edge() {};
		};
		std::vector<Edge> edges;
		struct Vertex : //0D elements
			public topology::Vertex<topology::lower::EmptyElement, topology::upper::SegmentToTetrahedron>,
			public geometry::Vertex,
			public functional::Shape<int, Point<double>>
		{
		public:
			Vertex() {};
			~Vertex() {};
		};
		std::vector<Vertex> vertexes;

	public:
		std::vector<BoundaryVertex_forScal_OrderBF2> boundary_vertexes;
		std::vector<BoundaryEdges_forScal_OrderBF2> boundary_edges;
		std::vector<BoundaryFace_forScal_OrderBF2> boundary_faces;

		bool is_print_logFile;

		std::vector<int> accordance_DOF_and_vertex;

		Grid_forScal_OrderBF2()
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
				topology::lower::TetrahedronToEdges test_tetrahedron;

				//create elements
				for (int id_element = 0; id_element < this->GetElementsCount(); id_element++)
				{
					this->GetElement(id_element)->SetSelfGlobalId(id_element);
					this->GetElement(id_element)->SetLowerElementCount(test_tetrahedron.GetLowerElementCount()); //in tetrahedron 6 edges
				}
				//create edge topology
				if (true)
				{
					struct edge_temp
					{
						int num1, num2;
						int element, local_in_element;
					public:
						edge_temp()
						{
							element = 0; local_in_element = 0;
							num1 = 0;
							num2 = 0;
						}
						edge_temp(int element, int local_in_element,
							std::vector<int>& n)
						{
							this->element = element; this->local_in_element = local_in_element;

							num1 = n[0];
							num2 = n[1];
							if (n[0] > n[1])
							{
								num1 = n[1];
								num2 = n[0];
							}
						}
						edge_temp(int element, int local_in_element, int n0, int n1)
						{
							this->element = element; this->local_in_element = local_in_element;

							num1 = n0;
							num2 = n1;
							if (n0 > n1)
							{
								num1 = n1;
								num2 = n0;
							}
						}

						void operator= (edge_temp A)
						{
							this->element = A.element; this->local_in_element = A.local_in_element;

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
					public:
						edge(edge_temp A)
						{
							element.push_back(A.element);
							local_in_element.push_back(A.local_in_element);

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
							element.push_back(A.element);
							local_in_element.push_back(A.local_in_element);

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
								if (A[j].num1 != A[i].num1 || A[j].num2 != A[i].num2) //םו ןמגעמנ
								{
									j = i;
									B.push_back(edge());
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

					std::vector <edge_temp> e_temp;
					std::vector <edge> edges_vector;
					std::vector<std::vector<int>> pattern;
					test_tetrahedron.GetLowerElementPatternInLocal(pattern);

					for (int id_tetr = 0; id_tetr < this->GetElementsCount(); id_tetr++)
					{
						for(int id_edge = 0; id_edge < pattern.size(); id_edge++)
						{
							int _id_elem = id_tetr;
							int _id_edge_in_elem = id_edge;
							e_temp.push_back(edge_temp(
								_id_elem, _id_edge_in_elem,
								this->GetElement(id_tetr)->GetIdNode(pattern[id_edge][0]),
								this->GetElement(id_tetr)->GetIdNode(pattern[id_edge][1])));
						}
					}
					math::MakeQuickSort(e_temp);
					MakeRemovalOfDuplication_e(e_temp, edges_vector);

					//add Upper and Lower Elements
					this->edges.resize(edges_vector.size());
					for (int id_edge = 0; id_edge < edges_vector.size(); id_edge++)
					{
						this->edges[id_edge].SetUpperElementCount((int)edges_vector[id_edge].element.size());
						this->edges[id_edge].SetSelfGlobalId(id_edge);

						for (int id_tetr = 0; id_tetr < edges_vector[id_edge].element.size(); id_tetr++)
						{
							this->edges[id_edge].SetUpperElement(id_tetr, this->GetElement(edges_vector[id_edge].element[id_tetr]));
							this->GetElement(edges_vector[id_edge].element[id_tetr])->SetLowerElement(edges_vector[id_edge].local_in_element[id_tetr], &this->edges[id_edge]);
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
						int edge, local_in_edge;
					public:
						node_temp()
						{
							edge = 0; local_in_edge = 0;
							num = 0;
						}
						node_temp(int edge, int local_in_edge, int n1)
						{
							this->edge = edge; this->local_in_edge = local_in_edge;

							num = n1;
						}

						void operator= (node_temp A)
						{
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

						std::vector<int> edge, local_in_edge;
					public:
						node(node_temp A)
						{
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
								if (A[j].num != A[i].num) //םו ןמגעמנ
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

						for (int id_vertex = 0; id_vertex < element1D->GetLowerElementCount(); id_vertex++)
						{
							int _id_edge = id_edge;
							int _id_vertex_in_edge = id_vertex;
							n_temp.push_back(node_temp(
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

			}
			catch (const std::exception&)
			{
				printf_s("Error: FEM/FEM_Grid.h/Grid_forScal_OrderBF2::CreateTopology()\n");
			}
		}

		void CreateFunctionalSpaces()
		{
			this->SetDOFsCount((int)this->vertexes.size() + (int)this->edges.size());
			for (int id_element = 0; id_element < this->GetElementsCount(); id_element++)
			{
				auto element = this->GetElement(id_element);
				element->ResizeDOF(element->GetNodesCount() + element->GetLowerElementCount()); //nodes+edges = 2 Order
				element->SolveAlphaMatrix();

				//nodes functions
				for (int _id_vertex = 0; _id_vertex < element->GetNodesCount(); _id_vertex++)
				{
					std::function< double(Point<double> X)> bf = [element, _id_vertex](Point<double> X) -> double
					{
						double result;
						double Li = element->alpha[_id_vertex][0] + element->alpha[_id_vertex][1] * X.x + element->alpha[_id_vertex][2] * X.y + element->alpha[_id_vertex][3] * X.z;
						result = Li * (2 * Li - 1);
						return result;
					};
					std::function< Point<double>(Point<double> X)> derivative_bf = [element, _id_vertex](Point<double> X) -> Point<double>
					{
						Point<double> result;
						double Li = element->alpha[_id_vertex][0] + element->alpha[_id_vertex][1] * X.x + element->alpha[_id_vertex][2] * X.y + element->alpha[_id_vertex][3] * X.z;
						double Li_dX = element->alpha[_id_vertex][1];
						double Li_dY = element->alpha[_id_vertex][2];
						double Li_dZ = element->alpha[_id_vertex][3];

						result.x = Li_dX * (2 * Li - 1) + Li * 2 * Li_dX;
						result.y = Li_dY * (2 * Li - 1) + Li * 2 * Li_dY;
						result.z = Li_dZ * (2 * Li - 1) + Li * 2 * Li_dZ;

						return result;
					};
					element->SetDOF(_id_vertex, element->GetIdNode(_id_vertex), bf, derivative_bf);
				}

				//edges functions
				
				for (int _local_id_edge = 0; _local_id_edge < element->GetLowerElementCount(); _local_id_edge++)
				{
					auto target_edge = &this->edges[element->GetLowerElement(_local_id_edge)->GetSelfGlobalId()];
					int _local_id_node_0 = element->GetLocalID_forVertexID(target_edge->GetIdNode(0));
					int _local_id_node_1 = element->GetLocalID_forVertexID(target_edge->GetIdNode(1));

					std::function< double(Point<double> X)> bf = [element, _local_id_node_0, _local_id_node_1](Point<double> X) -> double
					{
						double result;
						double L0 = element->alpha[_local_id_node_0][0] 
							+ element->alpha[_local_id_node_0][1] * X.x 
							+ element->alpha[_local_id_node_0][2] * X.y 
							+ element->alpha[_local_id_node_0][3] * X.z;
						double L1 = element->alpha[_local_id_node_1][0]
							+ element->alpha[_local_id_node_1][1] * X.x
							+ element->alpha[_local_id_node_1][2] * X.y
							+ element->alpha[_local_id_node_1][3] * X.z;

						result = 4 * L0 * L1;

						return result;
					};
					std::function< Point<double>(Point<double> X)> derivative_bf = [element, _local_id_node_0, _local_id_node_1](Point<double> X) -> Point<double>
					{
						Point<double> result;
						double L0 = element->alpha[_local_id_node_0][0]
							+ element->alpha[_local_id_node_0][1] * X.x
							+ element->alpha[_local_id_node_0][2] * X.y
							+ element->alpha[_local_id_node_0][3] * X.z;
						double L1 = element->alpha[_local_id_node_1][0]
							+ element->alpha[_local_id_node_1][1] * X.x
							+ element->alpha[_local_id_node_1][2] * X.y
							+ element->alpha[_local_id_node_1][3] * X.z;
						double L0_dX = element->alpha[_local_id_node_0][1];
						double L0_dY = element->alpha[_local_id_node_0][2];
						double L0_dZ = element->alpha[_local_id_node_0][3];
						double L1_dX = element->alpha[_local_id_node_1][1];
						double L1_dY = element->alpha[_local_id_node_1][2];
						double L1_dZ = element->alpha[_local_id_node_1][3];

						result.x = 4 * L0_dX * L1 + 4 * L0 * L1_dX;
						result.y = 4 * L0_dY * L1 + 4 * L0 * L1_dY;
						result.z = 4 * L0_dZ * L1 + 4 * L0 * L1_dZ;

						return result;
					};
					
					element->SetDOF(element->GetNodesCount() + _local_id_edge,
						vertexes.size() + element->GetLowerElement(_local_id_edge)->GetSelfGlobalId(), bf, derivative_bf);
				}
			}
		}
		template <typename Dirichlet, typename Neumann>
		void CreateBoundaryConditions(std::vector<Dirichlet>& first_boundary, std::vector<Neumann>& second_boundary)
		{
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
					std::vector<int> edges_on_boundaries(this->edges.size());
					//base vertexes (1 Order)
					for (int id_vertex = 0; id_vertex < first_boundary[id_type].id_vertexes.size(); id_vertex++)
					{
						int global_id_vertex = first_boundary[id_type].id_vertexes[id_vertex];
						this->boundary_vertexes[id].boundary_value = first_boundary[id_type].value;
						this->boundary_vertexes[id].SetIdNode(0, global_id_vertex);
						this->boundary_vertexes[id].SetNode(0, this->GetPtrCoordinateViaID(global_id_vertex));
						int id_tetr = this->vertexes[global_id_vertex].GetUpperElement(0)->GetUpperElement(0)->GetSelfGlobalId();
						this->boundary_vertexes[id].SetDOF(global_id_vertex, 
							*this->GetElement(id_tetr)->GetBasisFunctionInGlobalID(global_id_vertex),
							*this->GetElement(id_tetr)->GetDerivativeOfBasisFunctionInGlobalID(global_id_vertex));
						id++;

						//check additional edges-DOF for current boundary condition
						auto target_vertex = &(this->vertexes[global_id_vertex]);
						for (int i = 0; i < target_vertex->GetUpperElementCount(); i++)
						{
							auto target_edge = target_vertex->GetUpperElement(i);
							edges_on_boundaries[target_edge->GetSelfGlobalId()]++;
						}
					}

					//check additional edges-DOF for current boundary condition
					for (int id_edge = 0; id_edge < edges_on_boundaries.size(); id_edge++)
					{
						//all vertexes of edge are in boundary condition
						if (edges_on_boundaries[id_edge] >= 2)
						{
							int id_DOF = this->vertexes.size() + id_edge;
							int id_left_vertex = this->edges[id_edge].GetLowerElement(0)->GetSelfGlobalId();
							int id_right_vertex = this->edges[id_edge].GetLowerElement(0)->GetSelfGlobalId();
							this->edges[id_edge].geo_centr = (vertexes[id_left_vertex].GetNode(0) + vertexes[id_right_vertex].GetNode(0)) / 2.0;
							
							BoundaryVertex_forScal_OrderBF2 new_BV;// = new BoundaryVertex_forScal_OrderBF2();
							new_BV.boundary_value = first_boundary[id_type].value;
							new_BV.SetIdNode(0, id_DOF);
							new_BV.SetNode(0, &this->edges[id_edge].geo_centr); //ךמסעûכü!!!
							int id_tetr = this->edges[id_edge].GetUpperElement(0)->GetSelfGlobalId();
							new_BV.SetDOF(id_DOF,
								*this->GetElement(id_tetr)->GetBasisFunctionInGlobalID(id_DOF),
								*this->GetElement(id_tetr)->GetDerivativeOfBasisFunctionInGlobalID(id_DOF));

							this->boundary_vertexes.push_back(new_BV);
						}
					}
				}
			}

			//2nd boundary
			/*{
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

						this->boundary_faces[id].boundary_value = second_boundary[id_type].value;
						for (int id_vertex = 0; id_vertex < this->boundary_faces[id].GetNodesCount(); id_vertex++)
						{
							this->boundary_faces[id].SetIdNode(id_vertex, second_boundary[id_type].id_vertexes_as_triangle[id_face][id_vertex]);
							this->boundary_faces[id].SetNode(id_vertex, this->GetPtrCoordinateViaID(second_boundary[id_type].id_vertexes_as_triangle[id_face][id_vertex]));
						}
						id++;
					}
				}
			}*/

		}
		template <typename Dirichlet, typename Neumann>
		void Initialization(math::SimpleGrid& base_grid, std::vector<Dirichlet>& first_boundary, std::vector<Neumann>& second_boundary)
		{
			try {
				//create geometry space
				math::MakeCopyVector_A_into_B(base_grid.xyz, *(this->GetCoordinates()));
				this->CreateXYZline();
				this->SetElementsCount((int)base_grid.nvtr.size());
				for (int id_element = 0; id_element < this->GetElementsCount(); id_element++)
				{
					auto new_element = this->GetElement(id_element);

					std::vector<Point<double>*> P(base_grid.nvtr[id_element].size());
					for (int n = 0; n < P.size(); n++)
					{
						P[n] = this->GetPtrCoordinateViaID(base_grid.nvtr[id_element][n]);
					}
					new_element->SetGeometry(base_grid.nvtr[id_element], P);
					new_element->SetIdDomain(base_grid.nvkat[id_element]);
				}
				this->CreateQTree();

				//create topology
				CreateTopology();

				//create functional spaces
				CreateFunctionalSpaces();

				//create boundary conditions (geometry&topology&function)
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

		double GetSolutionInPoint(int id_element, Point<double> X, std::vector<double>& solution)
		{
			auto element = this->GetElement(id_element);
			double result = 0;
			for (int j = 0; j < element->GetDOFsCount(); j++)
			{
				auto bf = (*element->GetBasisFunctionInLocalID(j))(X);
				auto value = solution[element->GetDOFInLocalID(j)];
				result += bf * value;
			}
			return result;
		}
		Point<double> GetDerevativeFromSolutionInPoint(int id_element, Point<double> X, std::vector<double>& solution)
		{
			auto element = this->GetElement(id_element);
			Point<double> result;
			for (int j = 0; j < element->GetDOFsCount(); j++)
			{
				auto bf = (*element->GetDerivativeOfBasisFunctionInLocalID(j))(X);
				auto value = solution[element->GetDOFInLocalID(j)];
				result += bf * value;
			}
			return result;
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
		void printTecPlot3D_DiffDomains(/*char *directory,*/ FILE* fdat, std::vector<std::vector<double>>& value, std::vector<std::vector<char>> name_value, char* name_zone)
		{
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

			std::vector<int> in_domain(this->GetDomainsCount());
			for (int i = 0; i < this->GetElementsCount(); i++)
			{
				in_domain[this->GetElement(i)->GetIdDomain()]++;
			}

			for (int id_domain = 0; id_domain < this->GetDomainsCount(); id_domain++)
			{
				if (in_domain[id_domain] != 0)
				{
					fprintf_s(fdat, "ZONE T=\"%s_%d\"\n", name_zone, id_domain);
					fprintf_s(fdat, " N=%d,  E=%d, F=FEBLOCK ET=Tetrahedron \n", this->GetVertexCount(), in_domain[id_domain]);
					fprintf_s(fdat, " VARLOCATION=(NODAL NODAL NODAL");
					for (int i = 0; i < value.size(); i++)
					{
						if (value[i].size() == this->GetElementsCount()) fprintf_s(fdat, " CELLCENTERED");
						else fprintf_s(fdat, " NODAL");
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
						if (value[i].size() == this->GetElementsCount())
						{
							for (int j = 0; j < this->GetElementsCount(); j++)
								if (this->GetElement(j)->GetIdDomain() == id_domain)
									fprintf_s(fdat, "%.10e\n", value[i][j]);
						}
						else {
							for (int j = 0; j < value[i].size(); j++)
								fprintf_s(fdat, "%.10e\n", value[i][j]);
						}
						fprintf_s(fdat, "\n");
					}

					for (int i = 0; i < this->GetElementsCount(); i++)
					{
						if (this->GetElement(i)->GetIdDomain() == id_domain)
						{
							for (int j = 0; j < this->GetElement(i)->GetNodesCount(); j++)
								fprintf_s(fdat, "%d ", this->GetElement(i)->GetIdNode(j) + 1);
						}
						fprintf_s(fdat, "\n");
					}
				}
			}
			fclose(fdat);
		}


		~Grid_forScal_OrderBF2()
		{
			DOFs_count = 0;
			std::vector<BoundaryVertex_forScal_OrderBF2> v_b1;
			std::vector<BoundaryVertex_forScal_OrderBF2>(v_b1).swap(this->boundary_vertexes);
			std::vector<BoundaryFace_forScal_OrderBF2> v_b2;
			std::vector<BoundaryFace_forScal_OrderBF2>(v_b2).swap(this->boundary_faces);
			std::vector<int> v_p;
			std::vector<int>(v_p).swap(this->accordance_DOF_and_vertex);
			std::vector<Vertex> v_v;
			std::vector<Vertex>(v_v).swap(this->vertexes);

			this->DeleteGeometriGrid();
		}
	private:

	};
	   	
	class FiniteElement_forMech :
		public geometry::Tetrahedron,
		public topology::Tetrahedron<topology::lower::Vertex, topology::upper::EmptyElement>,
		public functional::Shape<int, Point<double>>
	{
	public:
		FiniteElement_forMech() { return; };
		~FiniteElement_forMech() { return; };

		//
		double v_mech;
		double eps_mech;

		void SolveLocalMatrix(DenseMatrix<Tensor2Rank3D, Point<double>> &local_matix, std::function<std::vector<std::vector<double>>(Point<double>)> &koefD)
		{
			try {
				local_matix.SetSize(this->GetDOFsCount());

				this->SetIntegrationLaw(4);

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

						if (true)
						{
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
								grad_right[4][2] = derevativeBF_I_inX.z.y;
								grad_right[4][1] = derevativeBF_I_inX.y.z;
								grad_right[5][2] = derevativeBF_I_inX.z.x;
								grad_right[5][0] = derevativeBF_I_inX.x.z;
								grad_right[3][1] = derevativeBF_I_inX.y.x;
								grad_right[3][0] = derevativeBF_I_inX.x.y;

								gradT_left[0][0] = derevativeBF_J_inX.x.x;
								gradT_left[1][1] = derevativeBF_J_inX.y.y;
								gradT_left[2][2] = derevativeBF_J_inX.z.z;
								gradT_left[2][4] = derevativeBF_J_inX.z.y;
								gradT_left[1][4] = derevativeBF_J_inX.y.z;
								gradT_left[2][5] = derevativeBF_J_inX.z.x;
								gradT_left[0][5] = derevativeBF_J_inX.x.z;
								gradT_left[1][3] = derevativeBF_J_inX.y.x;
								gradT_left[0][3] = derevativeBF_J_inX.x.y;

								std::vector<std::vector<double>> mult_gradT_D;
								math::MultMatrixMatrix(gradT_left, D, mult_gradT_D);
								std::vector<std::vector<double>> mult_gradT_D_grad;
								math::MultMatrixMatrix(mult_gradT_D, grad_right, mult_gradT_D_grad);

								result = mult_gradT_D_grad;
								//result = 1.0;
								return result;
							};

							double V = this->GetVolume();
							local_matix.A[_bf_local_id_I][_bf_local_id_J] = this->SolveIntegral(StiffnessMatrix);
							local_matix.A[_bf_local_id_J][_bf_local_id_I] = local_matix.A[_bf_local_id_I][_bf_local_id_J].T();
						}
						if (false)
						{
							eps_mech = 3e+7;
							v_mech = 0.3;
							double a = v_mech / (1 - v_mech);
							double b = (1 - 2 * v_mech) / (2 * (1 - v_mech));
							double k = eps_mech * (1 - v_mech) / ((1 + v_mech) * (1 - 2 * v_mech));

							std::function<double(Point<double>)> G_00 = [&D_basis_function_I, &D_basis_function_J, a, b, k](Point<double> X)->double
							{
								double res;
								Point<Point<double>> derevativeBF_I_inX = (*D_basis_function_I)(X);
								Point<Point<double>> derevativeBF_J_inX = (*D_basis_function_J)(X);
								res = derevativeBF_J_inX.x.x * k * derevativeBF_I_inX.x.x + derevativeBF_J_inX.x.y * k * b * derevativeBF_I_inX.x.y + derevativeBF_J_inX.x.z * k * b * derevativeBF_I_inX.x.z;
								return  res; 
							};
							local_matix.A[_bf_local_id_I][_bf_local_id_J].val[0][0] = this->SolveIntegral(G_00);
							std::function<double(Point<double>)> G_01 = [&D_basis_function_I, &D_basis_function_J, a, b, k](Point<double> X)->double
							{
								double res;
								Point<Point<double>> derevativeBF_I_inX = (*D_basis_function_I)(X);
								Point<Point<double>> derevativeBF_J_inX = (*D_basis_function_J)(X);
								res = derevativeBF_J_inX.x.x * k * a * derevativeBF_I_inX.x.y + derevativeBF_J_inX.x.y * k * b * derevativeBF_I_inX.x.x;;
								return  res;
							};
							local_matix.A[_bf_local_id_I][_bf_local_id_J].val[0][1] = this->SolveIntegral(G_01);
							std::function<double(Point<double>)> G_02 = [&D_basis_function_I, &D_basis_function_J, a, b, k](Point<double> X)->double
							{
								double res;
								Point<Point<double>> derevativeBF_I_inX = (*D_basis_function_I)(X);
								Point<Point<double>> derevativeBF_J_inX = (*D_basis_function_J)(X);
								res = derevativeBF_J_inX.x.x * k * a * derevativeBF_I_inX.x.z + derevativeBF_J_inX.x.z * k * b * derevativeBF_I_inX.x.x;
								return  res;
							};
							local_matix.A[_bf_local_id_I][_bf_local_id_J].val[0][2] = this->SolveIntegral(G_02);
							std::function<double(Point<double>)> G_10 = [&D_basis_function_I, &D_basis_function_J, a, b, k](Point<double> X)->double
							{
								double res;
								Point<Point<double>> derevativeBF_I_inX = (*D_basis_function_I)(X);
								Point<Point<double>> derevativeBF_J_inX = (*D_basis_function_J)(X);
								res = derevativeBF_J_inX.x.y * k * a * derevativeBF_I_inX.x.x + derevativeBF_J_inX.x.x * k * b * derevativeBF_I_inX.x.y;
								return  res;
							};
							local_matix.A[_bf_local_id_I][_bf_local_id_J].val[1][0] = this->SolveIntegral(G_10);
							std::function<double(Point<double>)> G_11 = [&D_basis_function_I, &D_basis_function_J, a, b, k](Point<double> X)->double
							{
								double res;
								Point<Point<double>> derevativeBF_I_inX = (*D_basis_function_I)(X);
								Point<Point<double>> derevativeBF_J_inX = (*D_basis_function_J)(X);
								res = derevativeBF_J_inX.x.y * k * derevativeBF_I_inX.x.y + derevativeBF_J_inX.x.x * k * b * derevativeBF_I_inX.x.x + derevativeBF_J_inX.x.z * k * b * derevativeBF_I_inX.x.z;
								return  res;
							};
							local_matix.A[_bf_local_id_I][_bf_local_id_J].val[1][1] = this->SolveIntegral(G_11);
							std::function<double(Point<double>)> G_12 = [&D_basis_function_I, &D_basis_function_J, a, b, k](Point<double> X)->double
							{
								double res;
								Point<Point<double>> derevativeBF_I_inX = (*D_basis_function_I)(X);
								Point<Point<double>> derevativeBF_J_inX = (*D_basis_function_J)(X);
								return derevativeBF_J_inX.x.y * k * a * derevativeBF_I_inX.x.z + derevativeBF_J_inX.x.z * k * b * derevativeBF_I_inX.x.y;
							};
							local_matix.A[_bf_local_id_I][_bf_local_id_J].val[1][2] = this->SolveIntegral(G_12);
							std::function<double(Point<double>)> G_20 = [&D_basis_function_I, &D_basis_function_J, a, b, k](Point<double> X)->double
							{
								double res;
								Point<Point<double>> derevativeBF_I_inX = (*D_basis_function_I)(X);
								Point<Point<double>> derevativeBF_J_inX = (*D_basis_function_J)(X);
								return derevativeBF_J_inX.x.z * k * a * derevativeBF_I_inX.x.x + derevativeBF_J_inX.x.x * k * b * derevativeBF_I_inX.x.z;
							};
							local_matix.A[_bf_local_id_I][_bf_local_id_J].val[2][0] = this->SolveIntegral(G_20);
							std::function<double(Point<double>)> G_21 = [&D_basis_function_I, &D_basis_function_J, a, b, k](Point<double> X)->double
							{
								double res;
								Point<Point<double>> derevativeBF_I_inX = (*D_basis_function_I)(X);
								Point<Point<double>> derevativeBF_J_inX = (*D_basis_function_J)(X);
								return derevativeBF_J_inX.x.z * k * a * derevativeBF_I_inX.x.y + derevativeBF_J_inX.x.y * k * b * derevativeBF_I_inX.x.z;
							};
							local_matix.A[_bf_local_id_I][_bf_local_id_J].val[2][1] = this->SolveIntegral(G_21);
							std::function<double(Point<double>)> G_22 = [&D_basis_function_I, &D_basis_function_J, a, b, k](Point<double> X)->double
							{
								double res;
								Point<Point<double>> derevativeBF_I_inX = (*D_basis_function_I)(X);
								Point<Point<double>> derevativeBF_J_inX = (*D_basis_function_J)(X);
								return derevativeBF_J_inX.x.z * k * derevativeBF_I_inX.x.z + derevativeBF_J_inX.x.y * k * b * derevativeBF_I_inX.x.y + derevativeBF_J_inX.x.x * k * b * derevativeBF_I_inX.x.x;
							};
							local_matix.A[_bf_local_id_I][_bf_local_id_J].val[2][2] = this->SolveIntegral(G_22);

						}
					}
				}

				/*for (int _bf_local_id_I = 0; _bf_local_id_I < this->GetDOFsCount(); _bf_local_id_I++)
				{
					for (int _bf_local_id_J = _bf_local_id_I+1; _bf_local_id_J < this->GetDOFsCount(); _bf_local_id_J++)
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
		void SolveLocalMatrix(DenseMatrix<Tensor2Rank3D, Point<double>>& local_matix, 
			std::function<std::vector<std::vector<double>>(Point<double>)>& koef_forStiffnessMatrix,
			std::function<double(Point<double>)>& koef_forMassMatrix,
			std::function<Point<double>(Point<double>)>& koef_forRightSide)
		{
			try {
				local_matix.SetSize(this->GetDOFsCount());
				std::vector<std::vector<Tensor2Rank3D>> MassMatrix;
				math::ResizeVector(MassMatrix, this->GetDOFsCount(), this->GetDOFsCount());

				this->SetIntegrationLaw(4);

				for (int _bf_local_id_I = 0; _bf_local_id_I < this->GetDOFsCount(); _bf_local_id_I++)
				{
					Tensor2Rank3D summ;
					Point<double> summ_vect;
					auto D_basis_function_I = this->GetDerivativeOfBasisFunctionInLocalID(_bf_local_id_I);
					auto basis_function_I = this->GetBasisFunctionInLocalID(_bf_local_id_I);
					for (int _bf_local_id_J = _bf_local_id_I; _bf_local_id_J < this->GetDOFsCount(); _bf_local_id_J++)
					{
						auto D_basis_function_J = this->GetDerivativeOfBasisFunctionInLocalID(_bf_local_id_J);
						auto basis_function_J = this->GetBasisFunctionInLocalID(_bf_local_id_J);

						Point<double> o = this->GetWeightCentr();
						auto val = (*this->GetDerivativeOfBasisFunctionInLocalID(0))(o);

						std::function<Tensor2Rank3D(Point<double>)> StiffnessMatrix = [&_bf_local_id_I, &_bf_local_id_J, &koef_forStiffnessMatrix, &D_basis_function_I, &D_basis_function_J]
						(Point<double> X) -> Tensor2Rank3D
						{
							Tensor2Rank3D result;

							std::vector<std::vector<double>> grad_right, gradT_left;
							math::ResizeVector(grad_right, 6, 3);
							math::ResizeVector(gradT_left, 3, 6);
							auto D = koef_forStiffnessMatrix(X);
							Point<Point<double>> derevativeBF_I_inX = (*D_basis_function_I)(X);
							Point<Point<double>> derevativeBF_J_inX = (*D_basis_function_J)(X);

							grad_right[0][0] = derevativeBF_I_inX.x.x;
							grad_right[1][1] = derevativeBF_I_inX.y.y;
							grad_right[2][2] = derevativeBF_I_inX.z.z;
							grad_right[4][2] = derevativeBF_I_inX.z.y;
							grad_right[4][1] = derevativeBF_I_inX.y.z;
							grad_right[5][2] = derevativeBF_I_inX.z.x;
							grad_right[5][0] = derevativeBF_I_inX.x.z;
							grad_right[3][1] = derevativeBF_I_inX.y.x;
							grad_right[3][0] = derevativeBF_I_inX.x.y;

							gradT_left[0][0] = derevativeBF_J_inX.x.x;
							gradT_left[1][1] = derevativeBF_J_inX.y.y;
							gradT_left[2][2] = derevativeBF_J_inX.z.z;
							gradT_left[2][4] = derevativeBF_J_inX.z.y;
							gradT_left[1][4] = derevativeBF_J_inX.y.z;
							gradT_left[2][5] = derevativeBF_J_inX.z.x;
							gradT_left[0][5] = derevativeBF_J_inX.x.z;
							gradT_left[1][3] = derevativeBF_J_inX.y.x;
							gradT_left[0][3] = derevativeBF_J_inX.x.y;

							std::vector<std::vector<double>> mult_gradT_D;
							math::MultMatrixMatrix(gradT_left, D, mult_gradT_D);
							std::vector<std::vector<double>> mult_gradT_D_grad;
							math::MultMatrixMatrix(mult_gradT_D, grad_right, mult_gradT_D_grad);

							result = mult_gradT_D_grad;
							//result = 1.0;
							return result;
						};
						std::function<Tensor2Rank3D(Point<double>)> MassMatrix_simple = [&_bf_local_id_I, &_bf_local_id_J, &koef_forMassMatrix, &basis_function_I, &basis_function_J]
						(Point<double> X) -> Tensor2Rank3D
						{
							Tensor2Rank3D result;

							Point<double> BF_I_inX = (*basis_function_I)(X);
							Point<double> BF_J_inX = (*basis_function_J)(X);

							//result = BF_I_inX * BF_J_inX;
							/*result.val[0][0] = BF_I_inX.x * BF_J_inX.x;
							result.val[0][1] = BF_I_inX.x * BF_J_inX.y;
							result.val[0][2] = BF_I_inX.x * BF_J_inX.z;

							result.val[1][0] = BF_I_inX.y * BF_J_inX.x;
							result.val[1][1] = BF_I_inX.y * BF_J_inX.y;
							result.val[1][2] = BF_I_inX.y * BF_J_inX.z;

							result.val[2][0] = BF_I_inX.z * BF_J_inX.x;
							result.val[2][1] = BF_I_inX.z * BF_J_inX.y;
							result.val[2][2] = BF_I_inX.z * BF_J_inX.z;*/

							result.val[0][0] = BF_I_inX.x * BF_J_inX.x;
							result.val[0][1] = 0;
							result.val[0][2] = 0;

							result.val[1][0] = 0;
							result.val[1][1] = BF_I_inX.y * BF_J_inX.y;
							result.val[1][2] = 0;

							result.val[2][0] = 0;
							result.val[2][1] = 0;
							result.val[2][2] = BF_I_inX.z * BF_J_inX.z;

							result = result * koef_forMassMatrix(X);
							return result;
						};

						//double V = this->GetVolume();
						local_matix.A[_bf_local_id_I][_bf_local_id_J] = this->SolveIntegral(StiffnessMatrix);
						//local_matix.A[_bf_local_id_I][_bf_local_id_J] += this->SolveIntegral(MassMatrix_simple);
						MassMatrix[_bf_local_id_I][_bf_local_id_J] += this->SolveIntegral(MassMatrix_simple);
						local_matix.A[_bf_local_id_I][_bf_local_id_J] += MassMatrix[_bf_local_id_I][_bf_local_id_J];
						local_matix.A[_bf_local_id_J][_bf_local_id_I] = local_matix.A[_bf_local_id_I][_bf_local_id_J].T();

						summ += local_matix.A[_bf_local_id_I][_bf_local_id_J];
					}

					std::function<Point<double>(Point<double>)> RightSide = [&basis_function_I, &koef_forRightSide]
					(Point<double> X) -> Point<double>
					{
						//Tensor2Rank3D result;
						Point<double> result;

						std::vector<std::vector<double>> f, phi;
						math::ResizeVector(f, 3, 1);
						math::ResizeVector(phi, 1, 3);
						Point<double> BF_I_inX = (*basis_function_I)(X);
						Point<double> F_inX = koef_forRightSide(X);

						f[0][0] = F_inX.x;
						f[1][0] = F_inX.y;
						f[2][0] = F_inX.z;

						phi[0][0] = BF_I_inX.x;
						phi[0][1] = BF_I_inX.y;
						phi[0][2] = BF_I_inX.z;

						result = Point<double>(F_inX.x * BF_I_inX.x, F_inX.y * BF_I_inX.y, F_inX.z * BF_I_inX.z);
						return result;
					};

					double V = this->GetVolume();
					local_matix.F[_bf_local_id_I] = this->SolveIntegral(RightSide);
					/*local_matix.F[_bf_local_id_I] = Point<double>(0, 0, 0);
					for (int j = 0; j < this->GetDOFsCount(); j++)
					{
						Point<double> X = this->GetNode(j);
						Point<double> F_inX = koef_forRightSide(X);
						local_matix.F[_bf_local_id_I] += MassMatrix[_bf_local_id_I][j] * F_inX;
					}*/

					summ_vect.x = summ.val[0][0] + summ.val[0][1] + summ.val[0][2];
					summ_vect.y = summ.val[1][0] + summ.val[1][1] + summ.val[1][2];
					summ_vect.z = summ.val[2][0] + summ.val[2][1] + summ.val[2][2];
				}
			}
			catch (const std::exception&)
			{
				printf_s("Error: XFEM/MultiXFEM_Grid.h/void SolveLocalMatrix(DenseMatrix<Tensor2Rank3D> &local_matix, std::function<std::vector<std::vector<double>>(Point)> &koefD)\n");
			}
		}
		void SolveRightSide(std::vector<Point<double>> &local_right_side, std::function<Point<double>(Point<double> X)>& sourse)
		{
			try {
				local_right_side.resize(this->GetDOFsCount());

				this->SetIntegrationLaw(4);

				for (int _bf_local_id_I = 0; _bf_local_id_I < this->GetDOFsCount(); _bf_local_id_I++)
				{
					auto basis_function_I = this->GetBasisFunctionInLocalID(_bf_local_id_I);

					Point<double> o = this->GetWeightCentr();
					auto val = (*this->GetDerivativeOfBasisFunctionInLocalID(0))(o);


					std::function<Point<double>(Point<double>)> MassMatrix = [&basis_function_I, &sourse]
					(Point<double> X) -> Point<double>
					{
						//Tensor2Rank3D result;
						Point<double> result;

						std::vector<std::vector<double>> f, phi;
						math::ResizeVector(f, 3, 1);
						math::ResizeVector(phi, 1, 3);
						Point<double> BF_I_inX = (*basis_function_I)(X);
						Point<double> F_inX = sourse(X);

						f[0][0] = F_inX.x;
						f[1][0] = F_inX.y;
						f[2][0] = F_inX.z;

						phi[0][0] = BF_I_inX.x;
						phi[0][1] = BF_I_inX.y;
						phi[0][2] = BF_I_inX.z;

						/*std::vector<std::vector<double>> mult_f_phi;
						math::MultMatrixMatrix(f, phi, mult_f_phi);*/

						//result = mult_f_phi;
						result = Point<double>(F_inX.x * BF_I_inX.x, F_inX.y * BF_I_inX.y, F_inX.z * BF_I_inX.z);
						return result;
					};

					double V = this->GetVolume();
					local_right_side[_bf_local_id_I] = this->SolveIntegral(MassMatrix);

				}
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
		std::function<Point<double>(Point<bool>&, int id_vertex)> boundary_value;
		int id_type;

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
		int id_type;

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
						/*func.x = 0.3333333333333333333;
						func.y = 0.3333333333333333333;
						func.z = 0.3333333333333333333;*/

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
	class Grid_forMech : public geometry::Grid<FiniteElement_forMech>
	{
		int DOFs_count;
		std::vector<topology::Vertex<topology::lower::EmptyElement, topology::upper::Tetrahedron>> vertexes; //Vertexes

	public:
		std::vector<geometry::Crack> cracks;
		std::vector<BoundaryVertex_forMech> boundary_vertexes;
		std::vector<BoundaryFace_forMech> boundary_faces;

		std::vector<int> accordance_DOF_and_vertex;

		Grid_forMech()
		{
			DOFs_count = 0;
		};

		topology::Vertex<topology::lower::EmptyElement, topology::upper::Tetrahedron>* GetTopologyVertex(int id)
		{
			return &(vertexes[id]);
		}

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
						if (A[j] != A[i]) //םו ןמגעמנ
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
					return Point<double>(0, 0, 0);
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
							this->boundary_vertexes[id].id_type = id_type;
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
							this->boundary_faces[id].id_type = id_type;

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

				for (int i = 0; i < this->GetElementsCount(); i++)
				{
					this->GetElement(i)->SetIdDomain(base_grid.nvkat[i]);
				}
			}
			catch (const std::exception&)
			{
				printf_s("Error: MultiXFEM/MultiXFEM_Grid.h/void Initialization(geometry::Grid<geometry::Tetrahedron> &base_grid)\n");
			}
		}
		template <typename Dirichlet, typename Neumann>
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
					if (id_elem % 1000 == 0)
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
		Tensor2Rank3D GetStrainTensorFromSolutionInPoint(int id_element, Point<double> X, Point<Point<double>> dU)
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
		Tensor2Rank3D GetStrainTensorFromSolutionInPoint(int id_element, Point<double> X, std::vector<Point<double>> &solution)
		{
			Point<Point<double>> dU = GetDerevativeFromSolutionInPoint(id_element, X, solution);

			return GetStrainTensorFromSolutionInPoint(id_element, X, dU);
		}
		/// <summary>
		/// SIGMA
		/// </summary>
		Tensor2Rank3D GetStressTensorFromSolutionInPoint(int id_element, Point<double> X, Tensor2Rank3D &StrainTensor)
		{
			Tensor2Rank3D SIG;
			auto _domain = this->GetDomain(this->GetElement(id_element)->GetIdDomain());
			SIG = _domain->forMech.solve_SIGMA(StrainTensor);

			return SIG;
		}
		/// <summary>
		/// SIGMA
		/// </summary>
		Tensor2Rank3D GetStressTensorFromSolutionInPoint(int id_element, Point<double> X, std::vector<Point<double>> &solution)
		{
			Tensor2Rank3D EPS = GetStrainTensorFromSolutionInPoint(id_element, X, solution);
			return GetStressTensorFromSolutionInPoint(id_element, X, EPS);
		}
		/// <summary>
		/// SIGMA
		/// </summary>
		Tensor2Rank3D GetStressTensorFromSolutionInPoint(int id_element, Point<double> X, Point<Point<double>>& dU)
		{
			Tensor2Rank3D EPS = GetStrainTensorFromSolutionInPoint(id_element, X, dU);
			return GetStressTensorFromSolutionInPoint(id_element, X, EPS);
		}
		double GetVonMisesStress(Tensor2Rank3D &Strees)
		{
			double mises_stress = sqrt((Strees.val[0][0] - Strees.val[1][1]) * (Strees.val[0][0] - Strees.val[1][1])
				+ (Strees.val[1][1] - Strees.val[2][2]) * (Strees.val[1][1] - Strees.val[2][2])
				+ (Strees.val[0][0] - Strees.val[2][2]) * (Strees.val[0][0] - Strees.val[2][2])
				+ 6 * (Strees.val[0][1] * Strees.val[1][0] + Strees.val[0][2] * Strees.val[2][0] + Strees.val[1][2] * Strees.val[2][1]) / 2.0);
			
			return mises_stress;
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

		void printTecPlot3D_DiffDomains(/*char *directory,*/ FILE* fdat, std::vector<std::vector<double>>& value, std::vector<std::vector<char>> name_value, char* name_zone)
		{
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

			std::vector<int> in_domain(this->GetDomainsCount());
			for (int i = 0; i < this->GetElementsCount(); i++)
			{
				in_domain[this->GetElement(i)->GetIdDomain()]++;
			}

			for (int id_domain = 0; id_domain < this->GetDomainsCount(); id_domain++)
			{
				if (in_domain[id_domain] != 0)
				{
					fprintf_s(fdat, "ZONE T=\"%s_%d\"\n", name_zone, id_domain);
					fprintf_s(fdat, " N=%d,  E=%d, F=FEBLOCK ET=Tetrahedron \n", this->GetVertexCount(), in_domain[id_domain]);
					fprintf_s(fdat, " VARLOCATION=(NODAL NODAL NODAL");
					for (int i = 0; i < value.size(); i++)
					{
						if (value[i].size() == this->GetElementsCount()) fprintf_s(fdat, " CELLCENTERED");
						else fprintf_s(fdat, " NODAL");
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
						if (value[i].size() == this->GetElementsCount())
						{
							for (int j = 0; j < this->GetElementsCount(); j++)
								if (this->GetElement(j)->GetIdDomain() == id_domain)
									fprintf_s(fdat, "%.10e\n", value[i][j]);
						}
						else {
							for (int j = 0; j < value[i].size(); j++)
								fprintf_s(fdat, "%.10e\n", value[i][j]);
						}
						fprintf_s(fdat, "\n");
					}

					for (int i = 0; i < this->GetElementsCount(); i++)
					{
						if (this->GetElement(i)->GetIdDomain() == id_domain)
						{
							for (int j = 0; j < this->GetElement(i)->GetNodesCount(); j++)
								fprintf_s(fdat, "%d ", this->GetElement(i)->GetIdNode(j) + 1);
						}
						fprintf_s(fdat, "\n");
					}
				}
			}
			fclose(fdat);
		}

		~Grid_forMech()
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

	class FiniteElement_1D_forMech :
		public geometry::Segment,
		public topology::Segment<topology::lower::Vertex, topology::upper::EmptyElement>,
		public functional::Shape<int, Point<double>>
	{
	public:
		FiniteElement_1D_forMech() { return; };
		~FiniteElement_1D_forMech() { return; };

		void SolveLocalMatrix(DenseMatrix<Tensor2Rank3D, Point<double>>& local_matix, std::function<std::vector<std::vector<double>>(Point<double>)>& koefD)
		{
			try {
				local_matix.SetSize(this->GetDOFsCount());

				this->SetIntegrationLaw(4);

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
							grad_right[4][2] = derevativeBF_I_inX.z.y;
							grad_right[4][1] = derevativeBF_I_inX.y.z;
							grad_right[5][2] = derevativeBF_I_inX.z.x;
							grad_right[5][0] = derevativeBF_I_inX.x.z;
							grad_right[3][1] = derevativeBF_I_inX.y.x;
							grad_right[3][0] = derevativeBF_I_inX.x.y;

							gradT_left[0][0] = derevativeBF_J_inX.x.x;
							gradT_left[1][1] = derevativeBF_J_inX.y.y;
							gradT_left[2][2] = derevativeBF_J_inX.z.z;
							gradT_left[2][4] = derevativeBF_J_inX.z.y;
							gradT_left[1][4] = derevativeBF_J_inX.y.z;
							gradT_left[2][5] = derevativeBF_J_inX.z.x;
							gradT_left[0][5] = derevativeBF_J_inX.x.z;
							gradT_left[1][3] = derevativeBF_J_inX.y.x;
							gradT_left[0][3] = derevativeBF_J_inX.x.y;

							std::vector<std::vector<double>> mult_gradT_D;
							math::MultMatrixMatrix(gradT_left, D, mult_gradT_D);
							std::vector<std::vector<double>> mult_gradT_D_grad;
							math::MultMatrixMatrix(mult_gradT_D, grad_right, mult_gradT_D_grad);

							result = mult_gradT_D_grad;
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
			}
			catch (const std::exception&)
			{
				printf_s("Error: XFEM/MultiXFEM_Grid.h/void SolveLocalMatrix(DenseMatrix<Tensor2Rank3D> &local_matix, std::function<std::vector<std::vector<double>>(Point)> &koefD)\n");
			}
		}

		void SolveRightSide(std::vector<Point<double>>& local_vector, std::function<Point<double>(Point<double>)>& sourse)
		{
			try {
				local_vector.resize(this->GetDOFsCount());

				this->SetIntegrationLaw(4);

				for (int _bf_local_id_I = 0; _bf_local_id_I < this->GetDOFsCount(); _bf_local_id_I++)
				{
					auto basis_function_I = this->GetBasisFunctionInLocalID(_bf_local_id_I);

					Point<double> o = this->GetWeightCentr();
					auto val = (*this->GetDerivativeOfBasisFunctionInLocalID(0))(o);

					std::function<Point<double>(Point<double>)> MassMatrix = [&basis_function_I, &sourse]
					(Point<double> X) -> Point<double>
					{
						//Tensor2Rank3D result;
						Point<double> result;

						std::vector<std::vector<double>> f, phi;
						math::ResizeVector(f, 3, 1);
						math::ResizeVector(phi, 1, 3);
						Point<double> BF_I_inX = (*basis_function_I)(X);
						Point<double> F_inX = sourse(X);

						f[0][0] = F_inX.x;
						f[1][0] = F_inX.y;
						f[2][0] = F_inX.z;

						phi[0][0] = BF_I_inX.x;
						phi[0][1] = BF_I_inX.y;
						phi[0][2] = BF_I_inX.z;

						/*std::vector<std::vector<double>> mult_f_phi;
						math::MultMatrixMatrix(f, phi, mult_f_phi);*/

						//result = mult_f_phi;
						result = Point<double>(F_inX.x * BF_I_inX.x, F_inX.y * BF_I_inX.y, F_inX.z * BF_I_inX.z);
						//result = Point<double>(1,1,1);
						return result;
					};

					double V = this->GetVolume();
					local_vector[_bf_local_id_I] = this->SolveIntegral(MassMatrix);

				}
			}
			catch (const std::exception&)
			{
				printf_s("Error: XFEM/MultiXFEM_Grid.h/void SolveLocalMatrix(DenseMatrix<Tensor2Rank3D> &local_matix, std::function<std::vector<std::vector<double>>(Point)> &koefD)\n");
			}
		}
	};
	class BoundaryVertex_1D_forMech :
		public geometry::Vertex,
		public topology::Vertex<topology::lower::EmptyElement, topology::upper::SegmentUpper>,
		public functional::Shape<int, Point<double>>
	{
	public:
		std::function<Point<double>(Point<bool>&)> boundary_value;

		BoundaryVertex_1D_forMech() { return; };
		~BoundaryVertex_1D_forMech() { return; };

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
	class Grid_1D_forMech : public geometry::Grid<FiniteElement_1D_forMech>
	{
		int DOFs_count;
		std::vector<topology::Vertex<topology::lower::EmptyElement, topology::upper::SegmentUpper>> vertexes; //Vertexes

	public:
		std::vector<BoundaryVertex_1D_forMech> boundary_vertexes;

		std::vector<int> accordance_DOF_and_vertex;

		Grid_1D_forMech()
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
				topology::lower::Segment tmp_1D;
				topology::lower::Vertex tmp_0D;
				std::vector<std::vector<int>> node_pattern;
				tmp_1D.GetLowerElementPatternInLocal(node_pattern);
				n_temp.reserve(node_pattern.size() * this->GetElementsCount());

				for (int id_element = 0; id_element < this->GetElementsCount(); id_element++)
				{
					auto element3D = this->GetElement(id_element);

					for (int id_vertex = 0; id_vertex < tmp_1D.GetLowerElementCount(); id_vertex++)
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
				math::MakeRemovalOfDuplication(n_temp, nodes_vector);

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

		void CreateFunctionalSpaces()
		{
			this->SetDOFsCount((int)this->vertexes.size());
			for (int id_element = 0; id_element < this->GetElementsCount(); id_element++)
			{
				auto element = this->GetElement(id_element);
				element->ResizeDOF(element->GetNodesCount());
				element->CreateBasis();

				for (int _id_vertex = 0; _id_vertex < element->GetNodesCount(); _id_vertex++)
				{
					std::function< Point<double>(Point<double> X)> bf = [element, _id_vertex](Point<double> X) -> Point<double>
					{
						Point<double> result;
						double self_result;
						Point<double> X_in_self = element->MakeTransferXyzIntoSelf(X);
						Point<double> X0_self = element->GetSelfNode(0);
						Point<double> X1_self = element->GetSelfNode(1);
						if (X0_self < X1_self)
						{
							switch (_id_vertex)
							{
							case 0:
								self_result = math::SolveLengthVector(X1_self, X_in_self) / math::SolveLengthVector(X1_self, X0_self);
								break;
							case 1:
								self_result = math::SolveLengthVector(X0_self, X_in_self) / math::SolveLengthVector(X1_self, X0_self);
								break;
							}
						}
						else {
							switch (_id_vertex)
							{
							case 1:
								self_result = math::SolveLengthVector(X1_self, X_in_self) / math::SolveLengthVector(X1_self, X0_self);
								break;
							case 0:
								self_result = math::SolveLengthVector(X0_self, X_in_self) / math::SolveLengthVector(X1_self, X0_self);
								break;
							}
						}
						result.x = self_result;
						result.y = self_result;
						result.z = self_result;

						return result;
					};
					std::function< Point<Point<double>>(Point<double> X)> derivative_bf = [element, _id_vertex](Point<double> X) -> Point<Point<double>>
					{
						Point<Point<double>> result;
						double self_result;
						Point<double> X_in_self = element->MakeTransferXyzIntoSelf(X);
						Point<double> X0_self = element->GetSelfNode(0);
						Point<double> X1_self = element->GetSelfNode(1);
						double size_elem = math::SolveLengthVector(X0_self, X1_self);

						if (X0_self < X1_self)
						{
							switch (_id_vertex)
							{
							case 0: self_result = -1. / size_elem;
								break;
							case 1: self_result = 1. / size_elem;
								break;
							}
						}
						else
						{
							switch (_id_vertex)
							{
							case 1: self_result = -1. / size_elem;
								break;
							case 0: self_result = 1. / size_elem;
								break;
							}
						}

						result.x.x = (*element->GetSelfReverseBasis())[0][0] * self_result + (*element->GetSelfReverseBasis())[1][0] * 0.0 + (*element->GetSelfReverseBasis())[2][0] * 0.0;
						result.x.y = (*element->GetSelfReverseBasis())[0][1] * self_result + (*element->GetSelfReverseBasis())[1][1] * 0.0 + (*element->GetSelfReverseBasis())[2][1] * 0.0;
						result.x.z = (*element->GetSelfReverseBasis())[0][2] * self_result + (*element->GetSelfReverseBasis())[1][2] * 0.0 + (*element->GetSelfReverseBasis())[2][2] * 0.0;
						/*result.x.x = (*element->GetSelfReverseBasis())[0][0] * self_result + (*element->GetSelfReverseBasis())[0][1] * 0.0 + (*element->GetSelfReverseBasis())[0][2] * 0.0;
						result.x.y = (*element->GetSelfReverseBasis())[1][0] * self_result + (*element->GetSelfReverseBasis())[1][1] * 0.0 + (*element->GetSelfReverseBasis())[1][2] * 0.0;
						result.x.z = (*element->GetSelfReverseBasis())[2][0] * self_result + (*element->GetSelfReverseBasis())[2][1] * 0.0 + (*element->GetSelfReverseBasis())[2][2] * 0.0;*/


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
					return Point<double>(0, 0, 0);
				};
				std::function< Point<Point<double>>(Point<double> X)> derivative_bf = [boundary](Point<double> X) -> Point<Point<double>>
				{
					Point<Point<double>> result;
					return result;
				};

				this->boundary_vertexes[id_vertex].ResizeDOF(1);
				this->boundary_vertexes[id_vertex].SetDOF(0, this->boundary_vertexes[id_vertex].GetIdNode(0), bf, derivative_bf);
			}
		}


		template <typename Dirichlet, typename Neumann>
		void Initialization(math::SimpleGrid& base_grid, std::vector<Dirichlet>& first_boundary, std::vector<Neumann>& second_boundary)
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


				CreateFunctionalSpaces();

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

		template <typename Matrix>
		void CreationPortrait(Matrix& matrix)
		{
			try {
				printf_s("\n");
				std::vector<std::vector<int>> tmp_down_columns(this->GetDOFsCount()), tmp_up_columns(this->GetDOFsCount());
				std::vector<std::vector<int>> down_columns(this->GetDOFsCount()), up_columns(this->GetDOFsCount());
				for (int id_elem = 0; id_elem < this->GetElementsCount(); id_elem++)
				{
					if (id_elem % 1000 == 0)
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

		~Grid_1D_forMech()
		{
			DOFs_count = 0;
			std::vector<BoundaryVertex_1D_forMech> v_b1;
			std::vector<BoundaryVertex_1D_forMech>(v_b1).swap(this->boundary_vertexes);
			std::vector<int> v_p;
			std::vector<int>(v_p).swap(this->accordance_DOF_and_vertex);
			std::vector<topology::Vertex<topology::lower::EmptyElement, topology::upper::SegmentUpper>> v_v;
			std::vector<topology::Vertex<topology::lower::EmptyElement, topology::upper::SegmentUpper>>(v_v).swap(this->vertexes);

			this->DeleteGeometriGrid();
		}
	private:

	};

	class FiniteElement_2D_forMech :
		public geometry::Triangle,
		public topology::Triangle<topology::lower::Vertex, topology::upper::EmptyElement>,
		public functional::Shape<int, Point<double>>
	{
	public:
		FiniteElement_2D_forMech() { return; };
		~FiniteElement_2D_forMech() { return; };

		void SolveLocalMatrix(DenseMatrix<Tensor2Rank3D, Point<double>>& local_matix, std::function<std::vector<std::vector<double>>(Point<double>)>& koefD)
		{
			try {
				local_matix.SetSize(this->GetDOFsCount());

				this->SetIntegrationLaw(3);

				for (int _bf_local_id_I = 0; _bf_local_id_I < this->GetDOFsCount(); _bf_local_id_I++)
				{
					for (int _bf_local_id_J = 0; _bf_local_id_J < this->GetDOFsCount(); _bf_local_id_J++)
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

							grad_right[3][0] = derevativeBF_I_inX.x.y;
							grad_right[3][1] = derevativeBF_I_inX.y.x;
							grad_right[4][2] = derevativeBF_I_inX.z.y;
							grad_right[4][1] = derevativeBF_I_inX.y.z;
							grad_right[5][2] = derevativeBF_I_inX.z.x;
							grad_right[5][0] = derevativeBF_I_inX.x.z;
							

							gradT_left[0][0] = derevativeBF_J_inX.x.x;
							gradT_left[1][1] = derevativeBF_J_inX.y.y;
							gradT_left[2][2] = derevativeBF_J_inX.z.z;

							gradT_left[2][4] = derevativeBF_J_inX.z.y;
							gradT_left[1][4] = derevativeBF_J_inX.y.z;
							gradT_left[2][5] = derevativeBF_J_inX.z.x;
							gradT_left[0][5] = derevativeBF_J_inX.x.z;
							gradT_left[1][3] = derevativeBF_J_inX.y.x;
							gradT_left[0][3] = derevativeBF_J_inX.x.y;

							std::vector<std::vector<double>> mult_gradT_D;
							math::MultMatrixMatrix(gradT_left, D, mult_gradT_D);
							std::vector<std::vector<double>> mult_gradT_D_grad;
							math::MultMatrixMatrix(mult_gradT_D, grad_right, mult_gradT_D_grad);

							result = mult_gradT_D_grad;
							//result = 1.0;
							return result;
						};

						double V = this->GetVolume();
						local_matix.A[_bf_local_id_I][_bf_local_id_J] = this->SolveIntegral(StiffnessMatrix);
						//local_matix.A[_bf_local_id_J][_bf_local_id_I] = local_matix.A[_bf_local_id_I][_bf_local_id_J].T();
					
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
			catch (const std::exception&)
			{
				printf_s("Error: XFEM/MultiXFEM_Grid.h/void SolveLocalMatrix(DenseMatrix<Tensor2Rank3D> &local_matix, std::function<std::vector<std::vector<double>>(Point)> &koefD)\n");
			}
		}
		void SolveRightSide(DenseMatrix<Tensor2Rank3D, Point<double>>& local_matix, std::function<Point<double>(Point<double>)>& sourse)
		{
			try {
				local_matix.SetSize(this->GetDOFsCount());

				this->SetIntegrationLaw(4);

				for (int _bf_local_id_I = 0; _bf_local_id_I < this->GetDOFsCount(); _bf_local_id_I++)
				{
					auto basis_function_I = this->GetBasisFunctionInLocalID(_bf_local_id_I);

					Point<double> o = this->GetWeightCentr();
					auto val = (*this->GetDerivativeOfBasisFunctionInLocalID(0))(o);

					std::function<Point<double>(Point<double>)> MassMatrix = [&basis_function_I, &sourse]
					(Point<double> X) -> Point<double>
					{
						//Tensor2Rank3D result;
						Point<double> result;

						std::vector<std::vector<double>> f, phi;
						math::ResizeVector(f, 3, 1);
						math::ResizeVector(phi, 1, 3);
						Point<double> BF_I_inX = (*basis_function_I)(X);
						Point<double> F_inX = sourse(X);

						f[0][0] = F_inX.x;
						f[1][0] = F_inX.y;
						f[2][0] = F_inX.z;

						phi[0][0] = BF_I_inX.x;
						phi[0][1] = BF_I_inX.y;
						phi[0][2] = BF_I_inX.z;

						/*std::vector<std::vector<double>> mult_f_phi;
						math::MultMatrixMatrix(f, phi, mult_f_phi);*/

						//result = mult_f_phi;
						result = Point<double>(F_inX.x * BF_I_inX.x, F_inX.y * BF_I_inX.y, F_inX.z * BF_I_inX.z);
						return result;
					};

					double V = this->GetVolume();
					local_matix.F[_bf_local_id_I] = this->SolveIntegral(MassMatrix);

				}
			}
			catch (const std::exception&)
			{
				printf_s("Error: XFEM/MultiXFEM_Grid.h/void SolveLocalMatrix(DenseMatrix<Tensor2Rank3D> &local_matix, std::function<std::vector<std::vector<double>>(Point)> &koefD)\n");
			}
		}
		void SolveRightSide(std::vector<Point<double>>& local_vector, std::function<Point<double>(Point<double>)>& sourse)
		{
			try {
				local_vector.resize(this->GetDOFsCount());

				this->SetIntegrationLaw(4);

				for (int _bf_local_id_I = 0; _bf_local_id_I < this->GetDOFsCount(); _bf_local_id_I++)
				{
					auto basis_function_I = this->GetBasisFunctionInLocalID(_bf_local_id_I);

					Point<double> o = this->GetWeightCentr();
					auto val = (*this->GetDerivativeOfBasisFunctionInLocalID(0))(o);

					std::function<Point<double>(Point<double>)> MassMatrix = [&basis_function_I, &sourse]
					(Point<double> X) -> Point<double>
					{
						//Tensor2Rank3D result;
						Point<double> result;

						std::vector<std::vector<double>> f, phi;
						math::ResizeVector(f, 3, 1);
						math::ResizeVector(phi, 1, 3);
						Point<double> BF_I_inX = (*basis_function_I)(X);
						Point<double> F_inX = sourse(X);

						f[0][0] = F_inX.x;
						f[1][0] = F_inX.y;
						f[2][0] = F_inX.z;

						phi[0][0] = BF_I_inX.x;
						phi[0][1] = BF_I_inX.y;
						phi[0][2] = BF_I_inX.z;

						/*std::vector<std::vector<double>> mult_f_phi;
						math::MultMatrixMatrix(f, phi, mult_f_phi);*/

						//result = mult_f_phi;
						result = Point<double>(F_inX.x * BF_I_inX.x, F_inX.y * BF_I_inX.y, F_inX.z * BF_I_inX.z);
						return result;
					};

					double V = this->GetVolume();
					local_vector[_bf_local_id_I] = this->SolveIntegral(MassMatrix);

				}
			}
			catch (const std::exception&)
			{
				printf_s("Error: XFEM/MultiXFEM_Grid.h/void SolveLocalMatrix(DenseMatrix<Tensor2Rank3D> &local_matix, std::function<std::vector<std::vector<double>>(Point)> &koefD)\n");
			}
		}
	};
	class BoundaryVertex_2D_forMech :
		public geometry::Vertex,
		public topology::Vertex<topology::lower::EmptyElement, topology::upper::TriangleUpper>,
		public functional::Shape<int, Point<double>>
	{
	public:
		std::function<Point<double>(Point<bool>&, int id_vertex)> boundary_value;

		BoundaryVertex_2D_forMech() { return; };
		~BoundaryVertex_2D_forMech() { return; };

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
	class Grid_2D_forMech : public geometry::Grid<FiniteElement_2D_forMech>
	{
		int DOFs_count;
		std::vector<topology::Vertex<topology::lower::EmptyElement, topology::upper::TriangleUpper>> vertexes; //Vertexes

	public:
		std::vector<BoundaryVertex_2D_forMech> boundary_vertexes;

		std::vector<int> accordance_DOF_and_vertex;

		Grid_2D_forMech()
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
							if (A[j].num != A[i].num) //םו ןמגעמנ
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
				topology::lower::TriangleToVertex tmp_2D;
				topology::lower::Vertex tmp_0D;
				std::vector<std::vector<int>> node_pattern;
				tmp_2D.GetLowerElementPatternInLocal(node_pattern);
				n_temp.reserve(node_pattern.size() * this->GetElementsCount());

				for (int id_element = 0; id_element < this->GetElementsCount(); id_element++)
				{
					auto element2D = this->GetElement(id_element);

					for (int id_vertex = 0; id_vertex < tmp_2D.GetLowerElementCount(); id_vertex++)
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
							element2D->GetIdNode(node_pattern[id_vertex][0])));
					}
				}

				math::MakeQuickSort(n_temp);
				MakeRemovalOfDuplication_n(n_temp, nodes_vector);

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

		void CreateFunctionalSpaces()
		{
			this->SetDOFsCount((int)this->vertexes.size());
			for (int id_element = 0; id_element < this->GetElementsCount(); id_element++)
			{
				auto element = this->GetElement(id_element);
				element->ResizeDOF(element->GetNodesCount());
				element->CreateBasis();

				for (int _id_vertex = 0; _id_vertex < element->GetNodesCount(); _id_vertex++)
				{
					std::function< Point<double>(Point<double> X)> bf = [element, _id_vertex](Point<double> X) -> Point<double>
					{
						Point<double> result;
						double self_result;
						Point<double> X_in_self = element->MakeTransferXyzIntoSelf(X);

						element->SolveAlphaMatrix();

						self_result = element->alpha[_id_vertex][0] + element->alpha[_id_vertex][1] * X_in_self.x + element->alpha[_id_vertex][2] * X_in_self.y;

						result.x = self_result;
						result.y = self_result;
						result.z = self_result;

						return result;
					};
					std::function< Point<Point<double>>(Point<double> X)> derivative_bf = [element, _id_vertex](Point<double> X) -> Point<Point<double>>
					{
						Point<Point<double>> result;
						element->SolveAlphaMatrix();

						result.x.x = (*element->GetSelfReverseBasis())[0][0] * element->alpha[_id_vertex][1] + (*element->GetSelfReverseBasis())[1][0] * element->alpha[_id_vertex][2];
						result.x.y = (*element->GetSelfReverseBasis())[0][1] * element->alpha[_id_vertex][1] + (*element->GetSelfReverseBasis())[1][1] * element->alpha[_id_vertex][2];
						result.x.z = (*element->GetSelfReverseBasis())[0][2] * element->alpha[_id_vertex][1] + (*element->GetSelfReverseBasis())[1][2] * element->alpha[_id_vertex][2];

						/*result.x.x = (*element->GetSelfBasis())[0][0] * element->alpha[_id_vertex][1] + (*element->GetSelfBasis())[0][1] * element->alpha[_id_vertex][2];
						result.x.y = (*element->GetSelfBasis())[1][0] * element->alpha[_id_vertex][1] + (*element->GetSelfBasis())[1][1] * element->alpha[_id_vertex][2];
						result.x.z = (*element->GetSelfBasis())[2][0] * element->alpha[_id_vertex][1] + (*element->GetSelfBasis())[2][1] * element->alpha[_id_vertex][2];*/

						/*result.x.x = (*element->GetSelfReverseBasis())[0][0] * element->alpha[_id_vertex][1] + (*element->GetSelfReverseBasis())[0][1] * element->alpha[_id_vertex][2] + (*element->GetSelfReverseBasis())[0][2] * 0.0;
						result.x.y = (*element->GetSelfReverseBasis())[1][0] * element->alpha[_id_vertex][1] + (*element->GetSelfReverseBasis())[1][1] * element->alpha[_id_vertex][2] + (*element->GetSelfReverseBasis())[1][2] * 0.0;
						result.x.z = (*element->GetSelfReverseBasis())[2][0] * element->alpha[_id_vertex][1] + (*element->GetSelfReverseBasis())[2][1] * element->alpha[_id_vertex][2] + (*element->GetSelfReverseBasis())[2][2] * 0.0;*/


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
					return Point<double>(0, 0, 0);
				};
				std::function< Point<Point<double>>(Point<double> X)> derivative_bf = [boundary](Point<double> X) -> Point<Point<double>>
				{
					Point<Point<double>> result;
					return result;
				};

				this->boundary_vertexes[id_vertex].ResizeDOF(1);
				this->boundary_vertexes[id_vertex].SetDOF(0, this->boundary_vertexes[id_vertex].GetIdNode(0), bf, derivative_bf);
			}
		}


		template <typename Dirichlet, typename Neumann>
		void Initialization(math::SimpleGrid& base_grid, std::vector<Dirichlet>& first_boundary, std::vector<Neumann>& second_boundary)
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


				CreateFunctionalSpaces();

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

		template <typename Matrix>
		void CreationPortrait(Matrix& matrix)
		{
			try {
				printf_s("\n");
				std::vector<std::vector<int>> tmp_down_columns(this->GetDOFsCount()), tmp_up_columns(this->GetDOFsCount());
				std::vector<std::vector<int>> down_columns(this->GetDOFsCount()), up_columns(this->GetDOFsCount());
				for (int id_elem = 0; id_elem < this->GetElementsCount(); id_elem++)
				{
					if (id_elem % 1000 == 0)
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
		void printTecPlot3D_DiffDomains(/*char *directory,*/ FILE* fdat, std::vector<std::vector<double>>& value, std::vector<std::vector<char>> name_value, char* name_zone)
		{
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

			std::vector<int> in_domain(this->GetDomainsCount());
			for (int i = 0; i < this->GetElementsCount(); i++)
			{
				in_domain[this->GetElement(i)->GetIdDomain()]++;
			}

			for (int id_domain = 0; id_domain < this->GetDomainsCount(); id_domain++)
			{
				if (in_domain[id_domain] != 0)
				{
					fprintf_s(fdat, "ZONE T=\"%s_%d\"\n", name_zone, id_domain);
					fprintf_s(fdat, " N=%d,  E=%d, F=FEBLOCK ET=Triangle \n", this->GetVertexCount(), in_domain[id_domain]);
					fprintf_s(fdat, " VARLOCATION=(NODAL NODAL NODAL");
					for (int i = 0; i < value.size(); i++)
					{
						if (value[i].size() == this->GetElementsCount()) fprintf_s(fdat, " CELLCENTERED");
						else fprintf_s(fdat, " NODAL");
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
						if (value[i].size() == this->GetElementsCount())
						{
							for (int j = 0; j < this->GetElementsCount(); j++)
								if (this->GetElement(j)->GetIdDomain() == id_domain)
									fprintf_s(fdat, "%.10lf\n", value[i][j]);
						}
						else {
							for (int j = 0; j < value[i].size(); j++)
								fprintf_s(fdat, "%.10lf\n", value[i][j]);
						}
						fprintf_s(fdat, "\n");
					}

					for (int i = 0; i < this->GetElementsCount(); i++)
					{
						if (this->GetElement(i)->GetIdDomain() == id_domain)
						{
							for (int j = 0; j < this->GetElement(i)->GetNodesCount(); j++)
								fprintf_s(fdat, "%d ", this->GetElement(i)->GetIdNode(j) + 1);
						}
						fprintf_s(fdat, "\n");
					}
				}
			}
			fclose(fdat);
		}

		~Grid_2D_forMech()
		{
			DOFs_count = 0;
			std::vector<BoundaryVertex_2D_forMech> v_b1;
			std::vector<BoundaryVertex_2D_forMech>(v_b1).swap(this->boundary_vertexes);
			std::vector<int> v_p;
			std::vector<int>(v_p).swap(this->accordance_DOF_and_vertex);
			std::vector<topology::Vertex<topology::lower::EmptyElement, topology::upper::TriangleUpper>> v_v;
			std::vector<topology::Vertex<topology::lower::EmptyElement, topology::upper::TriangleUpper>>(v_v).swap(this->vertexes);

			this->DeleteGeometriGrid();
		}
	private:

	};

	class FiniteElement_forMech_Order2 :
		public geometry::Tetrahedron,
		public topology::Tetrahedron<topology::lower::Segment, topology::upper::EmptyElement>,
		public functional::Shape<int, Point<double>>
	{
	public:
		FiniteElement_forMech_Order2() { return; };
		~FiniteElement_forMech_Order2() { return; };

		void SolveLocalMatrix(DenseMatrix<Tensor2Rank3D, Point<double>> &local_matix, std::function<std::vector<std::vector<double>>(Point<double>)> &koefD)
		{
			try {
				local_matix.SetSize(this->GetDOFsCount());

				this->SetIntegrationLaw(4);

				for (int _bf_local_id_I = 0; _bf_local_id_I < this->GetDOFsCount(); _bf_local_id_I++)
				{
					for (int _bf_local_id_J = 0; _bf_local_id_J < this->GetDOFsCount(); _bf_local_id_J++)
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
							grad_right[4][2] = derevativeBF_I_inX.z.y;
							grad_right[4][1] = derevativeBF_I_inX.y.z;
							grad_right[5][2] = derevativeBF_I_inX.z.x;
							grad_right[5][0] = derevativeBF_I_inX.x.z;
							grad_right[3][1] = derevativeBF_I_inX.y.x;
							grad_right[3][0] = derevativeBF_I_inX.x.y;

							gradT_left[0][0] = derevativeBF_J_inX.x.x;
							gradT_left[1][1] = derevativeBF_J_inX.y.y;
							gradT_left[2][2] = derevativeBF_J_inX.z.z;
							gradT_left[2][4] = derevativeBF_J_inX.z.y;
							gradT_left[1][4] = derevativeBF_J_inX.y.z;
							gradT_left[2][5] = derevativeBF_J_inX.z.x;
							gradT_left[0][5] = derevativeBF_J_inX.x.z;
							gradT_left[1][3] = derevativeBF_J_inX.y.x;
							gradT_left[0][3] = derevativeBF_J_inX.x.y;

							std::vector<std::vector<double>> mult_gradT_D;
							math::MultMatrixMatrix(gradT_left, D, mult_gradT_D);
							std::vector<std::vector<double>> mult_gradT_D_grad;
							math::MultMatrixMatrix(mult_gradT_D, grad_right, mult_gradT_D_grad);

							result = mult_gradT_D_grad;
							//result = 1.0;
							return result;
						};

						double V = this->GetVolume();
						local_matix.A[_bf_local_id_I][_bf_local_id_J] = this->SolveIntegral(StiffnessMatrix);
					}
				}
			}
			catch (const std::exception&)
			{
				printf_s("Error: XFEM/MultiXFEM_Grid.h/void SolveLocalMatrix(DenseMatrix<Tensor2Rank3D> &local_matix, std::function<std::vector<std::vector<double>>(Point)> &koefD)\n");
			}
		}
		void SolveLocalMatrix(DenseMatrix<Tensor2Rank3D, Point<double>>& local_matix,
			std::function<std::vector<std::vector<double>>(Point<double>)>& koef_forStiffnessMatrix,
			std::function<double(Point<double>)>& koef_forMassMatrix,
			std::function<Point<double>(Point<double>)>& koef_forRightSide)
		{
			try {
				local_matix.SetSize(this->GetDOFsCount());
				std::vector<std::vector<Tensor2Rank3D>> MassMatrix;
				math::ResizeVector(MassMatrix, this->GetDOFsCount(), this->GetDOFsCount());
				

				this->SetIntegrationLaw(4);

				for (int _bf_local_id_I = 0; _bf_local_id_I < this->GetDOFsCount(); _bf_local_id_I++)
				{
					//Tensor2Rank3D summ;
					//Point<double> summ_vect;
					auto D_basis_function_I = this->GetDerivativeOfBasisFunctionInLocalID(_bf_local_id_I);
					auto basis_function_I = this->GetBasisFunctionInLocalID(_bf_local_id_I);
					for (int _bf_local_id_J = 0; _bf_local_id_J < this->GetDOFsCount(); _bf_local_id_J++)
					{
						auto D_basis_function_J = this->GetDerivativeOfBasisFunctionInLocalID(_bf_local_id_J);
						auto basis_function_J = this->GetBasisFunctionInLocalID(_bf_local_id_J);

						Point<double> o = this->GetWeightCentr();
						auto val = (*this->GetDerivativeOfBasisFunctionInLocalID(0))(o);

						std::function<Tensor2Rank3D(Point<double>)> StiffnessMatrix = [&_bf_local_id_I, &_bf_local_id_J, &koef_forStiffnessMatrix, &D_basis_function_I, &D_basis_function_J]
						(Point<double> X) -> Tensor2Rank3D
						{
							Tensor2Rank3D result;

							std::vector<std::vector<double>> grad_right, gradT_left;
							math::ResizeVector(grad_right, 6, 3);
							math::ResizeVector(gradT_left, 3, 6);
							auto D = koef_forStiffnessMatrix(X);
							Point<Point<double>> derevativeBF_I_inX = (*D_basis_function_I)(X);
							Point<Point<double>> derevativeBF_J_inX = (*D_basis_function_J)(X);

							grad_right[0][0] = derevativeBF_I_inX.x.x;
							grad_right[1][1] = derevativeBF_I_inX.y.y;
							grad_right[2][2] = derevativeBF_I_inX.z.z;
							grad_right[4][2] = derevativeBF_I_inX.z.y;
							grad_right[4][1] = derevativeBF_I_inX.y.z;
							grad_right[5][2] = derevativeBF_I_inX.z.x;
							grad_right[5][0] = derevativeBF_I_inX.x.z;
							grad_right[3][1] = derevativeBF_I_inX.y.x;
							grad_right[3][0] = derevativeBF_I_inX.x.y;

							gradT_left[0][0] = derevativeBF_J_inX.x.x;
							gradT_left[1][1] = derevativeBF_J_inX.y.y;
							gradT_left[2][2] = derevativeBF_J_inX.z.z;
							gradT_left[2][4] = derevativeBF_J_inX.z.y;
							gradT_left[1][4] = derevativeBF_J_inX.y.z;
							gradT_left[2][5] = derevativeBF_J_inX.z.x;
							gradT_left[0][5] = derevativeBF_J_inX.x.z;
							gradT_left[1][3] = derevativeBF_J_inX.y.x;
							gradT_left[0][3] = derevativeBF_J_inX.x.y;

							std::vector<std::vector<double>> mult_gradT_D;
							math::MultMatrixMatrix(gradT_left, D, mult_gradT_D);
							std::vector<std::vector<double>> mult_gradT_D_grad;
							math::MultMatrixMatrix(mult_gradT_D, grad_right, mult_gradT_D_grad);

							result = mult_gradT_D_grad;
							//result = 1.0;
							return result;
						};
						std::function<Tensor2Rank3D(Point<double>)> MassMatrix_simple = [&_bf_local_id_I, &_bf_local_id_J, &koef_forMassMatrix, &basis_function_I, &basis_function_J]
						(Point<double> X) -> Tensor2Rank3D
						{
							Tensor2Rank3D result;

							Point<double> BF_I_inX = (*basis_function_I)(X);
							Point<double> BF_J_inX = (*basis_function_J)(X);

							//result = BF_I_inX * BF_J_inX;
							/*result.val[0][0] = BF_I_inX.x * BF_J_inX.x;
							result.val[0][1] = BF_I_inX.x * BF_J_inX.y;
							result.val[0][2] = BF_I_inX.x * BF_J_inX.z;

							result.val[1][0] = BF_I_inX.y * BF_J_inX.x;
							result.val[1][1] = BF_I_inX.y * BF_J_inX.y;
							result.val[1][2] = BF_I_inX.y * BF_J_inX.z;

							result.val[2][0] = BF_I_inX.z * BF_J_inX.x;
							result.val[2][1] = BF_I_inX.z * BF_J_inX.y;
							result.val[2][2] = BF_I_inX.z * BF_J_inX.z;*/

							/*result.val[0][0] = BF_I_inX.x * BF_J_inX.x;
							result.val[0][1] = result.val[0][0];
							result.val[0][2] = result.val[0][0];

							result.val[1][0] = result.val[0][0];
							result.val[1][1] = result.val[0][0];
							result.val[1][2] = result.val[0][0];

							result.val[2][0] = result.val[0][0];
							result.val[2][1] = result.val[0][0];
							result.val[2][2] = result.val[0][0];*/

							result.val[0][0] = BF_I_inX.x * BF_J_inX.x;
							result.val[0][1] = 0;
							result.val[0][2] = 0;

							result.val[1][0] = 0;
							result.val[1][1] = result.val[0][0];
							result.val[1][2] = 0;

							result.val[2][0] = 0;
							result.val[2][1] = 0;
							result.val[2][2] = result.val[0][0];

							result = result * koef_forMassMatrix(X);
							return result;
						};

						//double V = this->GetVolume();
						local_matix.A[_bf_local_id_I][_bf_local_id_J] = this->SolveIntegral(StiffnessMatrix);
						//local_matix.A[_bf_local_id_I][_bf_local_id_J] += this->SolveIntegral(MassMatrix_simple);
						MassMatrix[_bf_local_id_I][_bf_local_id_J] += this->SolveIntegral(MassMatrix_simple);
						local_matix.A[_bf_local_id_I][_bf_local_id_J] += MassMatrix[_bf_local_id_I][_bf_local_id_J];

						//summ += local_matix.A[_bf_local_id_I][_bf_local_id_J];
					}

					std::function<Point<double>(Point<double>)> RightSide = [&basis_function_I, &koef_forRightSide]
					(Point<double> X) -> Point<double>
					{
						//Tensor2Rank3D result;
						Point<double> result;

						std::vector<std::vector<double>> f, phi;
						math::ResizeVector(f, 3, 1);
						math::ResizeVector(phi, 1, 3);
						Point<double> BF_I_inX = (*basis_function_I)(X);
						Point<double> F_inX = koef_forRightSide(X);

						f[0][0] = F_inX.x;
						f[1][0] = F_inX.y;
						f[2][0] = F_inX.z;

						phi[0][0] = BF_I_inX.x;
						phi[0][1] = BF_I_inX.y;
						phi[0][2] = BF_I_inX.z;

						result = Point<double>(F_inX.x * BF_I_inX.x, F_inX.y * BF_I_inX.y, F_inX.z * BF_I_inX.z);
						return result;
					};

					double V = this->GetVolume();
					local_matix.F[_bf_local_id_I] = this->SolveIntegral(RightSide);
					//local_matix.F[_bf_local_id_I] = Point<double>(0, 0, 0);
					//for (int j = 0; j < this->GetDOFsCount(); j++)
					//{
					//	Point<double> X = Point<double>(0, 0, 0); //this->GetNode(j);
					//	Point<double> F_inX = koef_forRightSide(X);
					//	local_matix.F[_bf_local_id_I] += MassMatrix[_bf_local_id_I][j] * F_inX;
					//}

					
					/*summ_vect.x = summ.val[0][0] + summ.val[0][1] + summ.val[0][2];
					summ_vect.y = summ.val[1][0] + summ.val[1][1] + summ.val[1][2];
					summ_vect.z = summ.val[2][0] + summ.val[2][1] + summ.val[2][2];*/
				}
			}
			catch (const std::exception&)
			{
				printf_s("Error: XFEM/MultiXFEM_Grid.h/void SolveLocalMatrix(DenseMatrix<Tensor2Rank3D> &local_matix, std::function<std::vector<std::vector<double>>(Point)> &koefD)\n");
			}
		}
		void SolveRightSide(std::vector<Point<double>>& local_right_side, std::function<Point<double>(Point<double> X)>& sourse)
		{
			try {
				local_right_side.resize(this->GetDOFsCount());

				this->SetIntegrationLaw(4);

				for (int _bf_local_id_I = 0; _bf_local_id_I < this->GetDOFsCount(); _bf_local_id_I++)
				{
					auto basis_function_I = this->GetBasisFunctionInLocalID(_bf_local_id_I);

					Point<double> o = this->GetWeightCentr();
					auto val = (*this->GetDerivativeOfBasisFunctionInLocalID(0))(o);


					std::function<Point<double>(Point<double>)> MassMatrix = [&basis_function_I, &sourse]
					(Point<double> X) -> Point<double>
					{
						//Tensor2Rank3D result;
						Point<double> result;

						std::vector<std::vector<double>> f, phi;
						math::ResizeVector(f, 3, 1);
						math::ResizeVector(phi, 1, 3);
						Point<double> BF_I_inX = (*basis_function_I)(X);
						Point<double> F_inX = sourse(X);

						f[0][0] = F_inX.x;
						f[1][0] = F_inX.y;
						f[2][0] = F_inX.z;

						phi[0][0] = BF_I_inX.x;
						phi[0][1] = BF_I_inX.y;
						phi[0][2] = BF_I_inX.z;

						/*std::vector<std::vector<double>> mult_f_phi;
						math::MultMatrixMatrix(f, phi, mult_f_phi);*/

						//result = mult_f_phi;
						result = Point<double>(F_inX.x * BF_I_inX.x, F_inX.y * BF_I_inX.y, F_inX.z * BF_I_inX.z);
						return result;
					};

					double V = this->GetVolume();
					local_right_side[_bf_local_id_I] = this->SolveIntegral(MassMatrix);

				}
			}
			catch (const std::exception&)
			{
				printf_s("Error: XFEM/MultiXFEM_Grid.h/void SolveLocalMatrix(DenseMatrix<Tensor2Rank3D> &local_matix, std::function<std::vector<std::vector<double>>(Point)> &koefD)\n");
			}
		}
	};
	class BoundaryVertex_forMech_Order2 :
		public geometry::Vertex,
		//public topology::Vertex<topology::lower::EmptyElement, topology::upper::Segment>,
		public topology::Vertex<topology::lower::EmptyElement, topology::upper::Tetrahedron>,
		public functional::Shape<int, Point<double>>
	{
	public:
		int id_type;
		std::function<Point<double>(Point<bool>&, int)> boundary_value;

		BoundaryVertex_forMech_Order2() { return; };
		~BoundaryVertex_forMech_Order2() { return; };

		void SetTypeBoundary(int type)
		{
			this->type_boundary = type;
		}
		int GetTypeBoundary()
		{
			return this->type_boundary;
		}

		BoundaryVertex_forMech_Order2 operator = (const BoundaryVertex_forMech_Order2& A)
		{
			this->boundary_value = A.boundary_value;
			this->id_type = A.id_type;

			this->ResizeDOF(1);
			this->SetIdNode(0, A.id_nodes[0]);
			this->SetNode(0, A.nodes[0]);
			this->SetDOF(0, A.dofs[0], A.basis_functions[0], A.derivative_of_basis_functions[0]);

			return *this;
		}
		bool operator > (const BoundaryVertex_forMech_Order2& A)
		{
			if (this->GetDOFInLocalID(0) > A.dofs[0])
			{
				return true;
			}
			return false;
		}
		bool operator < (const BoundaryVertex_forMech_Order2& A)
		{
			if (this->GetDOFInLocalID(0) < A.dofs[0])
			{
				return true;
			}
			return false;
		}

	private:
		int type_boundary;

	};
	class BoundaryEdge_forMech_Order2 :
		public geometry::Vertex,
		//public topology::Vertex<topology::lower::EmptyElement, topology::upper::Segment>,
		public topology::Vertex<topology::lower::EmptyElement, topology::upper::Tetrahedron>,
		public functional::Shape<int, Point<double>>
	{
	public:
		int id_type;
		std::function<Point<double>(Point<bool>&)> boundary_value;

		BoundaryEdge_forMech_Order2() { return; };
		~BoundaryEdge_forMech_Order2() { return; };

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
	class BoundaryFace_forMech_Order2 :
		public geometry::Triangle,
		public topology::Triangle<topology::lower::Segment, topology::upper::Tetrahedron>,
		public functional::Shape<int, Point<double>>
	{
	public:
		int id_type;
		std::function<Point<double>(Point<double>)> boundary_value;

		BoundaryFace_forMech_Order2() { return; };
		~BoundaryFace_forMech_Order2() { return; };

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
	class Grid_forMech_Order2 : public geometry::Grid<FiniteElement_forMech_Order2>
	{
		int DOFs_count;
		struct Edge : //1D elements
			public topology::Segment<topology::lower::Vertex, topology::upper::Tetrahedron>,
			public geometry::Segment,
			public functional::Shape<int, Point<double>>
		{
		public:
			Edge() {};
			~Edge() {};
		};
		std::vector<Edge> edges;
		struct Vertex : //0D elements
			public topology::Vertex<topology::lower::EmptyElement, topology::upper::SegmentToTetrahedron>,
			public geometry::Vertex,
			public functional::Shape<int, Point<double>>
		{
		public:
			Vertex() {};
			~Vertex() {};
		};
		std::vector<Vertex> vertexes;

	public:
		std::vector<geometry::Crack> cracks;
		std::vector<BoundaryVertex_forMech_Order2> boundary_vertexes;
		std::vector<BoundaryFace_forMech_Order2> boundary_faces;

		std::vector<int> accordance_DOF_and_vertex;
		std::vector<int> accordance_DOF_and_edge;

		Grid_forMech_Order2()
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
		void SetDOFsCount(int count_vertexes, int count_edges)///!!!
		{
			this->DOFs_count = count_vertexes + count_edges;
			this->accordance_DOF_and_vertex.resize(count_vertexes);
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
					this->GetElement(id_element)->SetLowerElementCount(6); //Lower elements is edges
				}
				//create edge topology
				if (true)
				{
					struct edge_temp
					{
						int num1, num2;
						int element, local_in_element;
					public:
						edge_temp()
						{
							element = 0; local_in_element = 0;
							num1 = 0;
							num2 = 0;
						}
						edge_temp(int element, int local_in_element,
							std::vector<int>& n)
						{
							this->element = element; this->local_in_element = local_in_element;
							
							num1 = n[0];
							num2 = n[1];
							if (n[0] > n[1])
							{
								num1 = n[1];
								num2 = n[0];
							}
						}
						edge_temp(int element, int local_in_element,
							int n0, int n1)
						{
							this->element = element; this->local_in_element = local_in_element;

							num1 = n0;
							num2 = n1;
							if (n0 > n1)
							{
								num1 = n1;
								num2 = n0;
							}
						}

						void operator= (edge_temp A)
						{
							this->element = A.element; this->local_in_element = A.local_in_element;
							
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
					public:
						edge(edge_temp A)
						{
							element.push_back(A.element);
							local_in_element.push_back(A.local_in_element);

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
							element.push_back(A.element);
							local_in_element.push_back(A.local_in_element);
							//face.push_back(A.face);
							//local_in_face.push_back(A.local_in_face);

						}
					};
					

					std::vector <edge_temp> e_temp;
					std::vector <edge> edges_vector;
					topology::lower::TetrahedronToEdges tmp_3D;
					topology::lower::Segment tmp_1D;
					std::vector<std::vector<int>> edge_pattern;
					tmp_3D.GetLowerElementPatternInLocal(edge_pattern);
					e_temp.reserve(edge_pattern.size()* this->GetElementsCount());

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
								if (A[j].num1 != A[i].num1 || A[j].num2 != A[i].num2) //םו ןמגעמנ
								{
									j = i;
									B.push_back(edge());
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
					for (int id_element = 0; id_element < this->GetElementsCount(); id_element++)
					{
						auto element3D = this->GetElement(id_element);

						for (int id_edge = 0; id_edge < tmp_3D.GetLowerElementCount(); id_edge++)
						{
							int n1, n2;
							n1 = element3D->GetIdNode(edge_pattern[id_edge][0]);
							n2 = element3D->GetIdNode(edge_pattern[id_edge][1]);
							if (n2 < n1)
							{
								n2 = element3D->GetIdNode(edge_pattern[id_edge][0]);
								n1 = element3D->GetIdNode(edge_pattern[id_edge][1]);
							}
							e_temp.push_back(edge_temp(id_element, id_edge, n1, n2));
						}
					}
					math::MakeQuickSort(e_temp);
					MakeRemovalOfDuplication_e(e_temp, edges_vector);

					//add Upper and Lower Elements
					this->edges.resize(edges_vector.size());
					for (int id_edge = 0; id_edge < edges_vector.size(); id_edge++)
					{
						this->edges[id_edge].SetUpperElementCount((int)edges_vector[id_edge].element.size());
						this->edges[id_edge].SetSelfGlobalId(id_edge);

						for (int id_volume = 0; id_volume < edges_vector[id_edge].element.size(); id_volume++)
						{
							this->edges[id_edge].SetUpperElement(id_volume, this->GetElement(edges_vector[id_edge].element[id_volume]));
							this->GetElement(edges_vector[id_edge].element[id_volume])->SetLowerElement(edges_vector[id_edge].local_in_element[id_volume], &this->edges[id_edge]);
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
							if (A[j] != A[i]) //םו ןמגעמנ
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

					for (int id_edge = 0; id_edge < this->edges.size(); id_edge++)
					{
						auto element1D = &this->edges[id_edge];

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
					localRemovalOfDuplication_n(n_temp, nodes_vector);

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

				
			}
			catch (const std::exception&)
			{
				printf_s("Error: MultiXFEM/MultiXFEM_Grid.h/void CreateTopology()\n");
			}
		}

		void CreateFunctionalSpaces()
		{
			this->SetDOFsCount((int)this->vertexes.size(), (int)this->edges.size());
			for (int id_element = 0; id_element < this->GetElementsCount(); id_element++)
			{
				auto element = this->GetElement(id_element);
				element->ResizeDOF(element->GetNodesCount() + element->GetLowerElementCount());
				element->SolveAlphaMatrix();

				for (int _id_vertex = 0; _id_vertex < element->GetNodesCount(); _id_vertex++)
				{
					std::function< Point<double>(Point<double> X)> bf = [element, _id_vertex](Point<double> X) -> Point<double>
					{
						Point<double> result;

						double L1 = element->alpha[0][0] + element->alpha[0][1] * X.x + element->alpha[0][2] * X.y + element->alpha[0][3] * X.z;
						double L2 = element->alpha[1][0] + element->alpha[1][1] * X.x + element->alpha[1][2] * X.y + element->alpha[1][3] * X.z;
						double L3 = element->alpha[2][0] + element->alpha[2][1] * X.x + element->alpha[2][2] * X.y + element->alpha[2][3] * X.z;
						double L4 = element->alpha[3][0] + element->alpha[3][1] * X.x + element->alpha[3][2] * X.y + element->alpha[3][3] * X.z;

						switch (_id_vertex)
						{
						case 0: result.x = L1 * (2 * L1 - 1); break;
						case 1: result.x = L2 * (2 * L2 - 1); break;
						case 2: result.x = L3 * (2 * L3 - 1); break;
						case 3: result.x = L4 * (2 * L4 - 1); break;
						default: result.x = 0;
							break;
						}

						result.y = result.x;
						result.z = result.x;

						return result;
					};
					std::function< Point<Point<double>>(Point<double> X)> derivative_bf = [element, _id_vertex](Point<double> X) -> Point<Point<double>>
					{
						Point<Point<double>> result;

						double L1 = element->alpha[0][0] + element->alpha[0][1] * X.x + element->alpha[0][2] * X.y + element->alpha[0][3] * X.z;
						double L2 = element->alpha[1][0] + element->alpha[1][1] * X.x + element->alpha[1][2] * X.y + element->alpha[1][3] * X.z;
						double L3 = element->alpha[2][0] + element->alpha[2][1] * X.x + element->alpha[2][2] * X.y + element->alpha[2][3] * X.z;
						double L4 = element->alpha[3][0] + element->alpha[3][1] * X.x + element->alpha[3][2] * X.y + element->alpha[3][3] * X.z;
						
						double L1_dx = element->alpha[0][1];
						double L2_dx = element->alpha[1][1];
						double L3_dx = element->alpha[2][1];
						double L4_dx = element->alpha[3][1];

						double L1_dy = element->alpha[0][2];
						double L2_dy = element->alpha[1][2];
						double L3_dy = element->alpha[2][2];
						double L4_dy = element->alpha[3][2];

						double L1_dz = element->alpha[0][3];
						double L2_dz = element->alpha[1][3];
						double L3_dz = element->alpha[2][3];
						double L4_dz = element->alpha[3][3];

						switch (_id_vertex)
						{
						case 0: 
							result.x.x = L1_dx * (2 * L1 - 1) + L1 * 2 * L1_dx;
							result.x.y = L1_dy * (2 * L1 - 1) + L1 * 2 * L1_dy;
							result.x.z = L1_dz * (2 * L1 - 1) + L1 * 2 * L1_dz;
							break;
						case 1: 
							result.x.x = L2_dx * (2 * L2 - 1) + L2 * 2 * L2_dx;
							result.x.y = L2_dy * (2 * L2 - 1) + L2 * 2 * L2_dy;
							result.x.z = L2_dz * (2 * L2 - 1) + L2 * 2 * L2_dz;
							break;
						case 2:
							result.x.x = L3_dx * (2 * L3 - 1) + L3 * 2 * L3_dx;
							result.x.y = L3_dy * (2 * L3 - 1) + L3 * 2 * L3_dy;
							result.x.z = L3_dz * (2 * L3 - 1) + L3 * 2 * L3_dz;
							break;
						case 3: 
							result.x.x = L4_dx * (2 * L4 - 1) + L4 * 2 * L4_dx;
							result.x.y = L4_dy * (2 * L4 - 1) + L4 * 2 * L4_dy;
							result.x.z = L4_dz * (2 * L4 - 1) + L4 * 2 * L4_dz;
							break;
						default:
							result.x.x = 0;
							result.x.y = 0;
							result.x.z = 0;
							break;
						}

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
				for (int _id_edge = 0; _id_edge < element->GetLowerElementCount(); _id_edge++)
				{
					auto target_edge = &this->edges[element->GetLowerElement(_id_edge)->GetSelfGlobalId()];
					std::function< Point<double>(Point<double> X)> bf = [element, _id_edge](Point<double> X) -> Point<double>
					{
						Point<double> result;

						double L1 = element->alpha[0][0] + element->alpha[0][1] * X.x + element->alpha[0][2] * X.y + element->alpha[0][3] * X.z;
						double L2 = element->alpha[1][0] + element->alpha[1][1] * X.x + element->alpha[1][2] * X.y + element->alpha[1][3] * X.z;
						double L3 = element->alpha[2][0] + element->alpha[2][1] * X.x + element->alpha[2][2] * X.y + element->alpha[2][3] * X.z;
						double L4 = element->alpha[3][0] + element->alpha[3][1] * X.x + element->alpha[3][2] * X.y + element->alpha[3][3] * X.z;

						switch (_id_edge)
						{
						case 0: result.x = 4 * L1 * L2;	    break;
						case 1: result.x = 4 * L1 * L3;	    break;
						case 2: result.x = 4 * L1 * L4;	    break;
						case 3: result.x = 4 * L2 * L3;	    break;
						case 4: result.x = 4 * L2 * L4;	    break;
						case 5: result.x = 4 * L3 * L4;	    break;
						default: result.x = 0;
							break;
						}

						result.y = result.x;
						result.z = result.x;

						return result;
					};
					std::function< Point<Point<double>>(Point<double> X)> derivative_bf = [element, _id_edge](Point<double> X) -> Point<Point<double>>
					{
						Point<Point<double>> result;

						double L1 = element->alpha[0][0] + element->alpha[0][1] * X.x + element->alpha[0][2] * X.y + element->alpha[0][3] * X.z;
						double L2 = element->alpha[1][0] + element->alpha[1][1] * X.x + element->alpha[1][2] * X.y + element->alpha[1][3] * X.z;
						double L3 = element->alpha[2][0] + element->alpha[2][1] * X.x + element->alpha[2][2] * X.y + element->alpha[2][3] * X.z;
						double L4 = element->alpha[3][0] + element->alpha[3][1] * X.x + element->alpha[3][2] * X.y + element->alpha[3][3] * X.z;
						
						double L1_dx = element->alpha[0][1];
						double L2_dx = element->alpha[1][1];
						double L3_dx = element->alpha[2][1];
						double L4_dx = element->alpha[3][1];

						double L1_dy = element->alpha[0][2];
						double L2_dy = element->alpha[1][2];
						double L3_dy = element->alpha[2][2];
						double L4_dy = element->alpha[3][2];

						double L1_dz = element->alpha[0][3];
						double L2_dz = element->alpha[1][3];
						double L3_dz = element->alpha[2][3];
						double L4_dz = element->alpha[3][3];

						switch (_id_edge)
						{
						case 0:
							result.x.x = 4 * L1_dx * L2 + 4 * L1 * L2_dx;
							result.x.y = 4 * L1_dy * L2 + 4 * L1 * L2_dy;
							result.x.z = 4 * L1_dz * L2 + 4 * L1 * L2_dz;
							break;
						case 1:
							result.x.x = 4 * L1_dx * L3 + 4 * L1 * L3_dx;
							result.x.y = 4 * L1_dy * L3 + 4 * L1 * L3_dy;
							result.x.z = 4 * L1_dz * L3 + 4 * L1 * L3_dz;
							break;
						case 2:
							result.x.x = 4 * L1_dx * L4 + 4 * L1 * L4_dx;
							result.x.y = 4 * L1_dy * L4 + 4 * L1 * L4_dy;
							result.x.z = 4 * L1_dz * L4 + 4 * L1 * L4_dz;
							break;
						case 3:
							result.x.x = 4 * L2_dx * L3 + 4 * L2 * L3_dx;
							result.x.y = 4 * L2_dy * L3 + 4 * L2 * L3_dy;
							result.x.z = 4 * L2_dz * L3 + 4 * L2 * L3_dz;
							break;
						case 4:
							result.x.x = 4 * L2_dx * L4 + 4 * L2 * L4_dx;
							result.x.y = 4 * L2_dy * L4 + 4 * L2 * L4_dy;
							result.x.z = 4 * L2_dz * L4 + 4 * L2 * L4_dz;
							break;
						case 5:
							result.x.x = 4 * L3_dx * L4 + 4 * L3 * L4_dx;
							result.x.y = 4 * L3_dy * L4 + 4 * L3 * L4_dy;
							result.x.z = 4 * L3_dz * L4 + 4 * L3 * L4_dz;
							break;
						default:
							result.x.x = 0;
							result.x.y = 0;
							result.x.z = 0;
							break;
						}

						result.y.x = result.x.x;
						result.y.y = result.x.y;
						result.y.z = result.x.z;
						result.z.x = result.x.x;
						result.z.y = result.x.y;
						result.z.z = result.x.z;

						return result;
					};
					element->SetDOF(element->GetNodesCount() + _id_edge, this->GetVertexCount() + target_edge->GetSelfGlobalId(), bf, derivative_bf);
				}
			}

		}

		template <typename Dirichlet, typename Neumann>
		void CreateBoundaryConditions(std::vector<Dirichlet>& first_boundary, std::vector<Neumann>& second_boundary)
		{
			//1st boundary
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
					math::MakeQuickSort(global_dof_id);

					for (int i = 0; i < global_dof_id.size(); i++)
					{
						FEM::BoundaryVertex_forMech_Order2 temp;

						temp.boundary_value = first_boundary[id_type].value;
						temp.id_type = id_type;

						temp.ResizeDOF(1);
						if (global_dof_id[i] < this->vertexes.size())
						{
							temp.SetIdNode(0, global_dof_id[i]);
							temp.SetNode(0, this->GetPtrCoordinateViaID(global_dof_id[i]));
							int target_edge_id = this->vertexes[global_dof_id[i]].GetUpperElement(0)->GetSelfGlobalId();
							int target_elem_id = this->edges[target_edge_id].GetUpperElement(0)->GetSelfGlobalId();
							temp.SetDOF(0, global_dof_id[i], *this->GetElement(target_elem_id)->GetBasisFunctionInGlobalID(global_dof_id[i]), *this->GetElement(target_elem_id)->GetDerivativeOfBasisFunctionInGlobalID(global_dof_id[i]));
						}
						else
						{
							int target_elem_id = this->edges[global_dof_id[i] - this->vertexes.size()].GetUpperElement(0)->GetSelfGlobalId();
							temp.SetDOF(0, global_dof_id[i], *this->GetElement(target_elem_id)->GetBasisFunctionInGlobalID(global_dof_id[i]), *this->GetElement(target_elem_id)->GetDerivativeOfBasisFunctionInGlobalID(global_dof_id[i]));
						}
						this->boundary_vertexes.push_back(temp);
					}
				}
				math::MakeQuickSort(this->boundary_vertexes);
			}

			//2nd boundary
			if(true) {
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
						this->boundary_faces[id].SetLowerElementCount(3);
						this->boundary_faces[id].ResizeDOF(this->boundary_faces[id].GetLowerElementCount() + this->boundary_faces[id].GetNodesCount());
						this->boundary_faces[id].id_type = id_type;

						if (second_boundary[id_type].values.size() == 0) //non individual values for friangles
						{
							this->boundary_faces[id].boundary_value = second_boundary[id_type].value;
						}
						else {
							this->boundary_faces[id].boundary_value = second_boundary[id_type].values[id_face];
						}

						//DOF from vertexes
						for (int id_vertex = 0; id_vertex < this->boundary_faces[id].GetNodesCount(); id_vertex++)
						{
							int global_id_DOF = second_boundary[id_type].id_vertexes_as_triangle[id_face][id_vertex];
							this->boundary_faces[id].SetIdNode(id_vertex, second_boundary[id_type].id_vertexes_as_triangle[id_face][id_vertex]);
							this->boundary_faces[id].SetNode(id_vertex, this->GetPtrCoordinateViaID(second_boundary[id_type].id_vertexes_as_triangle[id_face][id_vertex]));
							this->boundary_faces[id].SetDOF(id_vertex, global_id_DOF,
								*this->GetElement(id_base_element)->GetBasisFunctionInGlobalID(global_id_DOF),
								*this->GetElement(id_base_element)->GetDerivativeOfBasisFunctionInGlobalID(global_id_DOF));
						}

						//DOF from edges
						topology::lower::Triangle test_triang;
						std::vector<std::vector<int>> pattern;
						test_triang.GetLowerElementPatternInLocal(pattern);
						for (int id_self_edge = 0; id_self_edge < pattern.size(); id_self_edge++)
						{
							for (int local_id_in_elem = 0; local_id_in_elem < this->GetElement(id_base_element)->GetLowerElementCount(); local_id_in_elem++)
							{
								auto target_edge = &this->edges[this->GetElement(id_base_element)->GetLowerElement(local_id_in_elem)->GetSelfGlobalId()];
								if ((this->boundary_faces[id].GetIdNode(pattern[id_self_edge][0]) == target_edge->GetIdNode(0) &&
									this->boundary_faces[id].GetIdNode(pattern[id_self_edge][1]) == target_edge->GetIdNode(1)) ||
									(this->boundary_faces[id].GetIdNode(pattern[id_self_edge][0]) == target_edge->GetIdNode(1) &&
										this->boundary_faces[id].GetIdNode(pattern[id_self_edge][1]) == target_edge->GetIdNode(0)))
								{
									this->boundary_faces[id].SetLowerElement(id_self_edge, target_edge);
									break;
								}
							}
						}
						for (int i = 0; i < this->boundary_faces[id].GetLowerElementCount(); i++)
						{
							auto lower_edge = &this->edges[this->boundary_faces[id].GetLowerElement(i)->GetSelfGlobalId()];
							int global_id_DOF = lower_edge->GetSelfGlobalId() + this->vertexes.size(); ///?????
							this->boundary_faces[id].SetDOF(this->boundary_faces[id].GetNodesCount() + i, global_id_DOF,
								*this->GetElement(id_base_element)->GetBasisFunctionInGlobalID(global_id_DOF),
								*this->GetElement(id_base_element)->GetDerivativeOfBasisFunctionInGlobalID(global_id_DOF));
						}
						id++;
					}
				}
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

				CreateFunctionalSpaces();

				CreateBoundaryConditions(first_boundary, second_boundary);

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

		void printTecPlot3D_DiffDomains(/*char *directory,*/ FILE* fdat, std::vector<std::vector<double>>& value, std::vector<std::vector<char>> name_value, char* name_zone)
		{
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

			std::vector<int> in_domain(this->GetDomainsCount());
			for (int i = 0; i < this->GetElementsCount(); i++)
			{
				in_domain[this->GetElement(i)->GetIdDomain()]++;
			}

			for (int id_domain = 0; id_domain < this->GetDomainsCount(); id_domain++)
			{
				if (in_domain[id_domain] != 0)
				{
					fprintf_s(fdat, "ZONE T=\"%s_%d\"\n", name_zone, id_domain);
					fprintf_s(fdat, " N=%d,  E=%d, F=FEBLOCK ET=Tetrahedron \n", this->GetVertexCount(), in_domain[id_domain]);
					fprintf_s(fdat, " VARLOCATION=(NODAL NODAL NODAL");
					for (int i = 0; i < value.size(); i++)
					{
						if (value[i].size() == this->GetElementsCount()) fprintf_s(fdat, " CELLCENTERED");
						else fprintf_s(fdat, " NODAL");
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
						if (value[i].size() == this->GetElementsCount())
						{
							for (int j = 0; j < this->GetElementsCount(); j++)
								if (this->GetElement(j)->GetIdDomain() == id_domain)
									fprintf_s(fdat, "%.10e\n", value[i][j]);
						}
						else {
							for (int j = 0; j < value[i].size(); j++)
								fprintf_s(fdat, "%.10e\n", value[i][j]);
						}
						fprintf_s(fdat, "\n");
					}

					for (int i = 0; i < this->GetElementsCount(); i++)
					{
						if (this->GetElement(i)->GetIdDomain() == id_domain)
						{
							for (int j = 0; j < this->GetElement(i)->GetNodesCount(); j++)
								fprintf_s(fdat, "%d ", this->GetElement(i)->GetIdNode(j) + 1);
						}
						fprintf_s(fdat, "\n");
					}
				}
			}
			fclose(fdat);
		}


		~Grid_forMech_Order2()
		{
			DOFs_count = 0;
			std::vector<geometry::Crack> v_cr;
			std::vector<geometry::Crack>(v_cr).swap(this->cracks);
			std::vector<BoundaryVertex_forMech_Order2> v_b1;
			std::vector<BoundaryVertex_forMech_Order2>(v_b1).swap(this->boundary_vertexes);
			std::vector<BoundaryFace_forMech_Order2> v_b2;
			std::vector<BoundaryFace_forMech_Order2>(v_b2).swap(this->boundary_faces);
			std::vector<int> v_p;
			std::vector<int>(v_p).swap(this->accordance_DOF_and_vertex);
			std::vector<Vertex> v_v;
			std::vector<Vertex>(v_v).swap(this->vertexes);
			std::vector<Edge> v_e;
			std::vector<Edge>(v_e).swap(this->edges);

			this->DeleteGeometriGrid();
		}
	private:

	};
}