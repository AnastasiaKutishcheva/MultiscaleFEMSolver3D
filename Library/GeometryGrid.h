#include "Point.h"
#include "QTree.h"
#include "Domain.h"

#pragma once
namespace geometry
{
template <class element> class Grid{
private:
	QTree<Grid> *q_tree;
	bool is_null_volume_in_grid;
	double volume;

	std::vector<element> elements;
	std::vector<Domain> domains;
	std::vector <double> X, Y, Z;
	std::vector<Point<double>> xyz;

	double GetPartVolume(Point<double> min, Point<double> max)
	{
		double part_volume = 0;
		int find_inside_points = 0;
		Point<double> middle = (max + min) / 2.0;

		if (this->GetElementID(min) >= 0) find_inside_points++;
		if (this->GetElementID(max) >= 0) find_inside_points++;
		if (this->GetElementID(Point<double>(max.x, min.y, min.z)) >= 0) find_inside_points++;
		if (this->GetElementID(Point<double>(min.x, max.y, min.z)) >= 0) find_inside_points++;
		if (this->GetElementID(Point<double>(max.x, max.y, min.z)) >= 0) find_inside_points++;
		if (this->GetElementID(Point<double>(min.x, min.y, max.z)) >= 0) find_inside_points++;
		if (this->GetElementID(Point<double>(max.x, min.y, max.z)) >= 0) find_inside_points++;
		if (this->GetElementID(Point<double>(min.x, max.y, max.z)) >= 0) find_inside_points++;
		if (find_inside_points == 0 && this->GetElementID(middle) >= 0) find_inside_points = -1;

		switch (find_inside_points)
		{
		case 0: 
			part_volume = 0; break;
		case 8: 
			part_volume = (max.x - min.x) * (max.y - min.y) * (max.z - min.z); break;

		default:
			if ((max.z - min.z) > (this->Z[this->Z.size() - 1] - this->Z[0]) * 1e-2)
			{
				part_volume += GetPartVolume(Point<double>(min.x, min.y, min.z), Point<double>(middle.x, middle.y, middle.z));
				part_volume += GetPartVolume(Point<double>(middle.x, min.y, min.z), Point<double>(max.x, middle.y, middle.z));
				part_volume += GetPartVolume(Point<double>(min.x, middle.y, min.z), Point<double>(middle.x, max.y, middle.z));
				part_volume += GetPartVolume(Point<double>(middle.x, middle.y, min.z), Point<double>(max.x, max.y, middle.z));

				part_volume += GetPartVolume(Point<double>(min.x, min.y, middle.z), Point<double>(middle.x, middle.y, max.z));
				part_volume += GetPartVolume(Point<double>(middle.x, min.y, middle.z), Point<double>(max.x, middle.y, max.z));
				part_volume += GetPartVolume(Point<double>(min.x, middle.y, middle.z), Point<double>(middle.x, max.y, max.z));
				part_volume += GetPartVolume(Point<double>(middle.x, middle.y, middle.z), Point<double>(max.x, max.y, max.z));
			}
			else {
				part_volume = (max.x - min.x) * (max.y - min.y) * (max.z - min.z) * find_inside_points / 8.0;
				if (find_inside_points == -1) part_volume = 0;
			}
			break;
		}

		return part_volume;
	}

public:
	
	Grid()
	{
	};
	void CreateQTree()
	{
		this->q_tree = new QTree<Grid>(*this);
	}
	void CreateXYZline()
	{
		std::vector <double> Xtemp, Ytemp, Ztemp;
		Xtemp.resize(xyz.size());
		Ytemp.resize(xyz.size());
		Ztemp.resize(xyz.size());
		for (int i = 0; i < xyz.size(); i++)
		{
			Xtemp[i] = xyz[i].x;
			Ytemp[i] = xyz[i].y;
			Ztemp[i] = xyz[i].z;
		}
		math::MakeQuickSort(Xtemp, 0, (int)Xtemp.size() - 1);
		math::MakeQuickSort(Ytemp, 0, (int)Ytemp.size() - 1);
		math::MakeQuickSort(Ztemp, 0, (int)Ztemp.size() - 1);
		math::MakeRemovalOfDuplication(Xtemp, X);
		math::MakeRemovalOfDuplication(Ytemp, Y);
		math::MakeRemovalOfDuplication(Ztemp, Z);
	}
	void UpdateCoordinates(std::vector<Point<double>> &new_xyz) {
		math::MakeCopyVector_A_into_B(new_xyz, this->xyz);
		this->CreateXYZline();
		this->CreateQTree();

		for (int i = 0; i < this->GetElementsCount(); i++)
		{
			for (int j = 0; j < this->elements[i].GetNodesCount(); j++)
			{
				this->elements[i].SetNode(j, &(this->xyz[this->elements[i].GetIdNodes(j)]));
			}
			this->elements[i].SolveAlphaMatrix();
		}
	}
	void MoveCoordinates(std::vector<Point<double>>& move_vectors) 
	{
		for (int i = 0; i < this->xyz.size(); i++)
		{
			this->xyz[i] += move_vectors[i];
		}
		this->CreateXYZline();
		this->CreateQTree();

		for (int i = 0; i < this->GetElementsCount(); i++)
		{
			for (int j = 0; j < this->elements[i].GetNodesCount(); j++)
			{
				this->elements[i].SetNode(j, &(this->xyz[this->elements[i].GetIdNode(j)]));
			}
			this->elements[i].SolveAlphaMatrix();
		}
	}
	void ReMoveCoordinates(std::vector<Point<double>>& move_vectors)
	{
		for (int i = 0; i < this->xyz.size(); i++)
		{
			this->xyz[i] -= move_vectors[i];
		}
		this->CreateXYZline();
		this->CreateQTree();

		for (int i = 0; i < this->GetElementsCount(); i++)
		{
			for (int j = 0; j < this->elements[i].GetNodesCount(); j++)
			{
				this->elements[i].SetNode(j, &(this->xyz[this->elements[i].GetIdNode(j)]));
			}
			this->elements[i].SolveAlphaMatrix();
		}
	}

	
	double GetFullVolumeByPart()
	{
		if (this->X.size() == 0) CreateXYZline();

		return this->GetPartVolume(this->GetMinCoordinate(), this->GetMaxCoordinate());
	}


	Point<double> GetCoordinateViaID(int id)
	{
		try
		{
			return this->xyz[id];
		}
		catch (const std::exception&)
		{
			printf_s("Error: Grid.h/GetCoordinateViaID(int id)\n");
		}
	}
	Point<double>* GetPtrCoordinateViaID(int id)
	{
		try
		{
			return &(this->xyz[id]);
		}
		catch (const std::exception&)
		{
			printf_s("Error: Grid.h/Point<double>* GetPtrCoordinateViaID(int id)\n");
		}
	}
	void GetPtrCoordinateViaID(std::vector<int> &id, std::vector<Point<double>*> &P)
	{
		try
		{
			P.resize(id.size());
			for (int i = 0; i < id.size(); i++)
			{
				P[i] = &(this->xyz[id[i]]);
			}
		}
		catch (const std::exception&)
		{
			printf_s("Error: Grid.h/void GetPtrCoordinateViaID(std::vector<int> &id, std::vector<Point<double>*> &P)\n");
		}
	}
	int GetVertexCount()
	{
		return (int)this->xyz.size();
	};
	std::vector<Point<double>>* GetCoordinates()
	{
		try
		{
			return &(this->xyz);
		}
		catch (const std::exception&)
		{
			printf_s("Error: Library/GeometryGrid.h/std::vector<Point<double>>* GetCoordinates()\n");
		}
	}
	Point <double> GetMinCoordinate()
	{
		if (this->X.size() == 0 || this->Y.size() == 0 || this->Z.size() == 0)
		{
			this->CreateXYZline();
		}
		return Point<double>(this->X[0], this->Y[0], this->Z[0]);
	}
	Point <double> GetMaxCoordinate()
	{
		if (this->X.size() == 0 || this->Y.size() == 0 || this->Z.size() == 0)
		{
			this->CreateXYZline();
		}
		return Point<double>(this->X[this->X.size() - 1], this->Y[this->Y.size() - 1], this->Z[this->Z.size() - 1]);
	}
	void MoveTheVertex(int id, Point<double> vector)
	{
		this->xyz[id] += vector;
	}

	void AddDomain()
	{
		this->domains.push_back(Domain());
	}
	void AddDomain(Domain &new_domain)
	{
		this->domains.push_back(new_domain);
	}
	Domain* GetDomain(int id)
	{
		try
		{
			return &(this->domains[id]);
		}
		catch (const std::exception&)
		{
			printf_s("Error: Library/GeometryGrid.h/Domain* GetDomain(int id)\n");
		}
	}
	int GetDomainsCount()
	{
		return this->domains.size();
	}

	int GetElementsCount()
	{
		return (int)this->elements.size();
	}
	element* GetElement(int id)
	{
		try
		{
			return &(this->elements[id]);
		}
		catch (const std::exception&)
		{
			printf_s("Error: Library/GeometryGrid.h/element* GetElement(int id)\n");
		}
	}
	void ReservationElementsCount(int count)
	{
		this->elements.Reserv(count);
	}
	void SetElementsCount(int count)
	{
		this->elements.resize(count);
	}
	double GetVolume()
	{

		try
		{
			if (volume < 0 || math::IsEqual(volume, 0.0))
			{
				bool is_null_volume_in_grid = false;
				if (this->elements.size() != 0)
				{
					volume = 0;
					for(int el = 0; el < this->elements.size(); el++)
					{
						volume += this->elements[el].GetVolume();
						if (this->elements[el].GetVolume() <= 0)
						{
							is_null_volume_in_grid = true;
							volume = -1;
							break;
						}
					}
				}
			}
			return volume;
		}
		catch (const std::exception&)
		{
			printf_s("Error: Grid.h/GetCoordinateViaID(int id)\n");
		}
	};
	
	int GetElementID( Point<double> A )
	{
		try
		{
			int num = -1;
			if (this->q_tree == NULL)
			{
				this->CreateQTree();
			}
			if (this->q_tree->IsElemInQTree(A, num, *this))
				return num;
			else return -1;
		}
		catch (const std::exception&)
		{
			printf_s("Error: Grid.h/ GetElementID( Point A )\n");
		}
		
	}
	int GetNearestElementID(Point<double> A, double &lenght)
	{
		try
		{
			int num = -1;
			if (this->q_tree == NULL)
			{
				this->CreateQTree();
			}
			if (this->q_tree->IsNearestElemInQTree(A, num, lenght, *this))
				return num;
			else return -1 * num;
		}
		catch (const std::exception&)
		{
			
		}
	}

	void Zoom(double coef)
	{
		try
		{
			for (int i = 0; i < this->xyz.size(); i++)
			{
				this->xyz[i] *= coef;
			}
			for (int i = 0; i < this->X.size(); i++)
				X[i] *= coef;
			for (int i = 0; i < this->Y.size(); i++)
				Y[i] *= coef;
			for (int i = 0; i < this->Z.size(); i++)
				Z[i] *= coef;

			for (int i = 0; i < this->elements.size(); i++)
			{
				this->elements[i].UpdateOfGeometricalProperties();
			}
		}
		catch (const std::exception&)
		{
			printf_s("Error: Grid.h/ Zoom(double coef)\n");
		}
		
	}

	~Grid()
	{
		std::vector<element> v_el;
		std::vector<element>(v_el).swap(this->elements);
		std::vector<Domain> v_dom;
		std::vector<Domain>(v_dom).swap(this->domains);
		std::vector<Point<double>> v_p;
		std::vector<Point<double>>(v_p).swap(this->xyz);
	}

	protected:
	void DeleteGeometriGrid()
	{
		this->~Grid();
	}
};
}