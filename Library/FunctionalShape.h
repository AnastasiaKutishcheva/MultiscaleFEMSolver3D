#pragma once
#include <vector>
#include <functional>
#include "Point.h"

namespace functional
{
	template <typename TypeDOF, typename TypeBF>
	class Shape
	{
	public:
		int dofs_count;
		Shape()
		{
			dofs_count = -1;
		};
		~Shape()
		{
			std::vector<TypeDOF> v_dofs;
			std::vector<TypeDOF>(v_dofs).swap(this->dofs);
			std::vector<std::function<TypeBF(Point<double>)>> v_basis_functions;
			std::vector<std::function<TypeBF(Point<double>)>>(v_basis_functions).swap(this->basis_functions);
			std::vector<std::function<Point<TypeBF>(Point<double>)>> v_d_basis_functions;
			std::vector<std::function<Point<TypeBF>(Point<double>)>>(v_d_basis_functions).swap(this->derivative_of_basis_functions);
		};

		int GetDOFsCount()
		{
			if (dofs_count > 0)
				return dofs_count;
			return (int)dofs.size();
		}
		std::vector<TypeDOF>* GetElementDOFs()
		{
			return &(this->dofs);
		}
		TypeDOF GetDOFInLocalID(int local_DOF_id)
		{
			try {
				return this->dofs[local_DOF_id];
			}
			catch (const std::exception&)
			{
				printf_s("Error: Library/FunctionalShape.h/TypeDOF GetDOFInLocalID(int local_DOF_id)\n");
			}
		}
		int GetLocalID_forDOF(TypeDOF global_DOF_id)
		{
			int result = -1;
			for (int i = 0; i < this->dofs.size(); i++)
			{
				if (this->dofs[i] == global_DOF_id)
				{
					result = i;
					break;
				}
			}
			return result;
		}

		void SetDOF(TypeDOF global_DOF_id, std::function<TypeBF(Point<double>)> &basis_function, std::function<Point<TypeBF>(Point<double>)> &derivative_of_basis_function)
		{
			try {
				this->dofs_count++;
				this->dofs.push_back(global_DOF_id);
				this->basis_functions.push_back(basis_function);
				this->derivative_of_basis_functions.push_back(derivative_of_basis_function);
			}
			catch (const std::exception&)
			{
				printf_s("Error: Library/FunctionalShape.h/void SetDOF(TypeDOF global_DOF_id, std::function<TypeBF(Point<double>)>* basis_function, std::function<Point<TypeBF>(Point<double>)>* derivative_of_basis_function)\n");
			}
		}
		void SetDOF(int local_DOF_id, TypeDOF global_DOF_id, std::function<TypeBF(Point<double>)> basis_function, std::function<Point<TypeBF>(Point<double>)> derivative_of_basis_function)
		{
			try {
				/*for (int i = (int)this->dofs.size()-1; i < local_DOF_id; i++)
				{
					TypeDOF _tmp_dof;
					this->dofs.push_back(_tmp_dof);
					this->basis_functions.push_back(NULL);
					this->derivative_of_basis_functions.push_back(NULL);
				}*/
				this->dofs[local_DOF_id] = global_DOF_id;
				this->basis_functions[local_DOF_id] = basis_function;
				this->derivative_of_basis_functions[local_DOF_id] = derivative_of_basis_function;
			}
			catch (const std::exception&)
			{
				printf_s("Error: Library/FunctionalShape.h/void SetDOF(int local_DOF_id, TypeDOF global_DOF_id, std::function<TypeBF(Point<double>)>* basis_function, std::function<Point<TypeBF>(Point<double>)>* derivative_of_basis_function)\n");
			}
		}
		void AppEndDOF(TypeDOF global_DOF_id, std::function<TypeBF(Point<double>)> &basis_function, std::function<Point<TypeBF>(Point<double>)> &derivative_of_basis_function)
		{
			try {
				this->dofs_count++;
				this->dofs.push_back(global_DOF_id);
				this->basis_functions.push_back(basis_function);
				this->derivative_of_basis_functions.push_back(derivative_of_basis_function);
			}
			catch (const std::exception&)
			{
				printf_s("Error: Library/FunctionalShape.h/void SetDOF(int local_DOF_id, TypeDOF global_DOF_id, std::function<TypeBF(Point<double>)>* basis_function, std::function<Point<TypeBF>(Point<double>)>* derivative_of_basis_function)\n");
			}
		}
		

		void ResizeDOF(int size)
		{
			try {
				this->dofs_count = size;
				this->dofs.resize(size);
				this->basis_functions.resize(size);
				this->derivative_of_basis_functions.resize(size);
			}
			catch (const std::exception&)
			{
				printf_s("Error: Library/FunctionalShape.h/void ResizeDOF(int size)\n");
			}
		}

		std::function<TypeBF(Point<double>)>* GetBasisFunctionInLocalID(int local_id)
		{
			try {
				return &basis_functions[local_id];
			}
			catch (const std::exception&)
			{
				printf_s("Error: Library/FunctionalShape.h/std::function<TypeBF(Point)>* GetBasisFunctionInLocalID(int local_id)\n");
			}
		}
		std::function<TypeBF(Point<double>)>* GetBasisFunctionInGlobalID(int global_id)
		{
			try {
				int local_id = -1;
				for(local_id = 0; local_id < this->dofs.size(); local_id++)
				{
					if (global_id == this->dofs[local_id]) break;
				}
				return &basis_functions[local_id];
			}
			catch (const std::exception&)
			{
				printf_s("Error: Library/FunctionalShape.h/std::function<TypeBF(Point)>* GetBasisFunctionInLocalID(int local_id)\n");
			}
		}
		//( (../dx), (.../dy), (.../dz) )
		std::function<Point<TypeBF>(Point<double>)>* GetDerivativeOfBasisFunctionInLocalID(int local_id)
		{
			try {
				return &derivative_of_basis_functions[local_id];
			}
			catch (const std::exception&)
			{
				printf_s("Error: Library/FunctionalShape.h/std::function<Point<TypeBF>(Point<double>)>* GetDerivativeOfBasisFunctionInLocalID(int local_id)\n");
			}
		}
		std::function<Point<TypeBF>(Point<double>)>* GetDerivativeOfBasisFunctionInGlobalID(int global_id)
		{
			try {
				int local_id = -1;
				for (local_id = 0; local_id < this->dofs.size(); local_id++)
				{
					if (global_id == this->dofs[local_id]) break;
				}
				return &derivative_of_basis_functions[local_id];
			}
			catch (const std::exception&)
			{
				printf_s("Error: Library/FunctionalShape.h/std::function<TypeBF(Point)>* GetBasisFunctionInLocalID(int local_id)\n");
			}
		}

		
	//private:
		std::vector<TypeDOF> dofs; //global degrees of freedom in local numeration
		std::vector<std::function<TypeBF(Point<double>)>> basis_functions;
		std::vector<std::function<Point<TypeBF>(Point<double>)>> derivative_of_basis_functions;
	};
}