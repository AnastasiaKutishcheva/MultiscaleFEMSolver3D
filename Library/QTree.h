#pragma once
#include "Point.h"
#include <vector>

template <class Grid> 
class QTree {
private:
	Point<double> bottom, top;
	std::vector<int> elems;
	QTree *next[8];
	int level;

	void CreateBranch(QTree *box, Grid &grid)
	{
		int maxk = (int)log(grid.GetElementsCount())*10;
		if (maxk == 0) maxk = grid.GetElementsCount();
		if (this->level < 20 && box->elems.size() > maxk)
		{
			Point<double> Center = (box->top + box->bottom) / 2;
			int next_level = this->level + 1;

			box->next[0] = new QTree(next_level, box->bottom, Center, box, grid);
			box->next[1] = new QTree(next_level, Point<double>(Center.x, box->bottom.y, box->bottom.z),
				Point<double>(box->top.x, Center.y, Center.z),
				box, grid);
			box->next[2] = new QTree(next_level, Point<double>(box->bottom.x, Center.y, box->bottom.z),
				Point<double>(Center.x, box->top.y, Center.z),
				box, grid);
			box->next[3] = new QTree(next_level, Point<double>(Center.x, Center.y, box->bottom.z),
				Point<double>(box->top.x, box->top.y, Center.z),
				box, grid);

			box->next[4] = new QTree(next_level, Point<double>(box->bottom.x, box->bottom.y, Center.z),
				Point<double>(Center.x, Center.y, box->top.z),
				box, grid);
			box->next[5] = new QTree(next_level, Point<double>(Center.x, box->bottom.y, Center.z),
				Point<double>(box->top.x, Center.y, box->top.z),
				box, grid);
			box->next[6] = new QTree(next_level, Point<double>(box->bottom.x, Center.y, Center.z),
				Point<double>(Center.x, box->top.y, box->top.z),
				box, grid);
			box->next[7] = new QTree(next_level, Point<double>(Center.x, Center.y, Center.z),
				Point<double>(box->top.x, box->top.y, box->top.z),
				box, grid);

			std::vector<int> v;
			std::vector<int>(v).swap(box->elems);
		}
		else {
			return;
		}
	}
	bool IsElemInQTree(Point<double> A, QTree *box, int &num, Grid &grid)
	{
		auto next = box->next[0];
		int size = box->elems.size();
		//if (box->next[0] != NULL || box->elems.size() == 0)//значит еще не лист
		if (box->next[0] != NULL && box->elems.size() == 0)//значит еще не лист
		{
			//определ€ем в каком сегменте лежит
			for (int i = 0; i < 8; i++)
			{
				if (box->next[i]->IsInclude(A))
				{
					//printf( "ENTER segm %d\t", i );
					IsElemInQTree(A, box->next[i], num, grid);
					break;
				}
			}
		}
		else { //уже лист
			for (int i = 0; i < box->elems.size(); i++)
				if (grid.GetElement(box->elems[i])->IsContainThePoint(A))
				{
					num = box->elems[i];
					return true;
				}
			return false; // так и не нашли элемента
		}
		return true;
	}
	bool IsNearestElemInQTree(Point<double> A, QTree *box, int &num, double &len, Grid &grid)
	{
		auto next = box->next[0];
		int size = box->elems.size();
		//if (box->next[0] != NULL || box->elems.size() == 0)//значит еще не лист
		if (box->next[0] != NULL && box->elems.size() == 0)//значит еще не лист
		{
			//определ€ем в каком сегменте лежит
			for (int i = 0; i < 8; i++)
			{
				if (box->next[i]->IsInclude(A))
				{
					//printf( "ENTER segm %d\t", i );
					//SearchInTree(A, box->next[i], num);
					IsNearestElemInQTree(A, box->next[i], num, len, grid);
					break;
				}
			}
		}
		else { //уже лист
			for (int i = 0; i < box->elems.size(); i++)
			{
				if (grid.GetElement(box->elems[i])->IsContainThePoint(A))
				{
					num = box->elems[i];
					len = 0.0;
					return true;
				}
			}
			double discrepancy = 0;
			if (box->elems.size() != 0)
			{
				grid.GetElement(box->elems[0])->IsContainThePoint(A, len);
				num = box->elems[0];
				for (int i = 1; i < box->elems.size(); i++)
				{
					grid.GetElement(box->elems[i])->IsContainThePoint(A, discrepancy);
					if (discrepancy < len)
					{
						len = discrepancy;
						num = box->elems[i];
					}
				}
			}
			return false; // так и не нашли элемента
		}
		return true;
	}
public:
	QTree()
	{
		try
		{
			level = 0;
			for (int i = 0; i < 8; i++) this->next[i] = NULL;
		}
		catch (const std::exception&)
		{
			printf_s("Error: QTree.h/QTree()\n");
		}
	}
	QTree(Grid &grid)
	{
		try {
			this->elems.resize(grid.GetElementsCount());
			for (int i = 0; i < grid.GetElementsCount(); i++)
				this->elems[i] = i;
						
			this->bottom = grid.GetMinCoordinate();
			this->top = grid.GetMaxCoordinate();

			//this->prev = NULL;
			for (int i = 0; i < 8; i++) this->next[i] = NULL;
			this->level = 0;

			CreateBranch(this, grid);
		}
		catch (const std::exception&)
		{
			printf_s("Error: QTree.h/QTree(Grid &grid)\n");
		}
	}
	QTree(int level, Point<double> X1, Point<double> X2, QTree *prev, Grid &grid)
	{
		try {
			//this->prev = prev;
			for (int i = 0; i < 8; i++) this->next[i] = NULL;
			this->level = level;

			if (X1 < X2)
			{
				this->bottom = X1;
				this->top = X2;
			}
			else {
				this->bottom = X2;
				this->top = X1;
			}

			for (int i = 0; i < prev->elems.size(); i++)
			{
				if (grid.GetElement(prev->elems[i])->IsPartOfElemInBox(this->bottom, this->top))
					this->elems.push_back(prev->elems[i]);
			}
			CreateBranch(this, grid);
		}
		catch (const std::exception&)
		{
			printf_s("Error: QTree.h/QTree(Point<double> X1, Point<double> X2, Tree *prev, Grid &grid)\n");
		}
	}
	
	bool IsElemInQTree(Point<double> A, int &num, Grid &grid)
	{
		try
		{
			return this->IsElemInQTree(A, this, num, grid);
		}
		catch (const std::exception&)
		{
			printf_s("Error: QTree.h/IsElemInQTree\n");
		}
		
	}
	bool IsNearestElemInQTree(Point<double> A, int &num, double &lenght, Grid &grid)
	{
		try
		{
			return this->IsNearestElemInQTree(A, this, num, lenght, grid);
		}
		catch (const std::exception&)
		{
			printf_s("Error: QTree.h/IsNearestElemInQTree\n");
		}
		
	}
	

	bool IsInclude(Point<double> X)
	{
		if (this->bottom.x <= X.x && X.x <= this->top.x &&
			this->bottom.y <= X.y && X.y <= this->top.y &&
			this->bottom.z <= X.z && X.z <= this->top.z)
			return true;
		return false;
	}
};