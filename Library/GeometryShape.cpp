#pragma once
#include "GeometryShape.h"
#include "Math.h"
#include <vector>



//==================Vertex============================================
geometry::Vertex::Vertex()
{
	try
	{
	}
	catch (const std::exception&)
	{
		printf_s("Error: GeometryShape.cpp/geometry::Vertex::Vertex()");
	}
}
geometry::Vertex::Vertex(std::vector<Point<double> *> &P)
{
	try
	{
		this->UpdateNodes(std::vector<Point<double> *> P);
	}
	catch (const std::exception&)
	{
		printf_s("Error: GeometryShape.cpp/geometry::Triangle::Triangle(std::vector<Point<double> *> &P)");
	}
}
double geometry::Vertex::GetVolume()
{
	try
	{
		return 0;
	}
	catch (const std::exception&)
	{
		printf_s("Error: GeometryShape.cpp/geometry::Triangle::GetVolume()");
	}
}

bool geometry::Vertex::IsContainThePoint(Point<double> A)
{
	try
	{
		return false;
	}
	catch (const std::exception&)
	{
		printf_s("Error: GeometryShape.cpp/geometry::Tetrahedron::IsContainThePoint(Point<double> A)");
	}
}

//==================Triangle============================================
geometry::Tetrahedron::Tetrahedron()
{
	try
	{
	}
	catch (const std::exception&)
	{
		printf_s("Error: GeometryShape.cpp/geometry::Tetrahedron::Tetrahedron()");
	}
}
geometry::Triangle::Triangle(std::vector<Point<double> *> &P)
{
	try
	{
		this->UpdateNodes(std::vector<Point<double> *> P);
	}
	catch (const std::exception&)
	{
		printf_s("Error: GeometryShape.cpp/geometry::Triangle::Triangle(std::vector<Point<double> *> &P)");
	}
}
double geometry::Triangle::GetVolume()
{
	try
	{
		Point<double> _nodes[4];
		for (int i = 0; i < 4; i++)
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
		printf_s("Error: GeometryShape.cpp/geometry::Triangle::GetVolume()");
	}
}

bool geometry::Triangle::IsContainThePoint(Point<double> A)
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
		P[2] = X;
		V1 = geometry::Triangle(P[0], P[1], P[2]).Square();
		P[0] = this->GetNode(0);
		P[1] = X;
		P[2] = this->GetNode(2);
		V2 = geometry::Triangle(P[0], P[1], P[2]).Square();
		P[0] = X;
		P[1] = this->GetNode(1);
		P[2] = this->GetNode(2);
		V3 = geometry::Triangle(P[0], P[1], P[2]).Square();

		double summ = V1 + V2 + V3;
		double res = round(abs(summ - V) / V, 5);
		if (res <= 1E-4)
			return true;
		return false;
	}
	catch (const std::exception&)
	{
		printf_s("Error: GeometryShape.cpp/geometry::Tetrahedron::IsContainThePoint(Point<double> A)");
	}
}

//==================Tetrahedron============================================
geometry::Tetrahedron::Tetrahedron()
{
	try
	{
	}
	catch (const std::exception&)
	{
		printf_s("Error: GeometryShape.cpp/geometry::Tetrahedron::Tetrahedron(std::vector<Point<double> *> &P)");
	}
}
geometry::Tetrahedron::Tetrahedron(std::vector<Point<double> *> &P)
{
	try
	{
		this->UpdateNodes(std::vector<Point<double> *> P);
	}
	catch (const std::exception&)
	{
		printf_s("Error: GeometryShape.cpp/geometry::Tetrahedron::Tetrahedron(std::vector<Point<double> *> &P)");
	}
}
double geometry::Tetrahedron::GetVolume()
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
		printf_s("Error: GeometryShape.cpp/geometry::Tetrahedron::GetVolume()");
	}
}
double geometry::Tetrahedron::GetVolume(Point<double>[4] _P)
{
	try
	{
		detD_1 = _P[1].z*_P[2].y*_P[0].x + _P[3].z*_P[0].x*_P[1].y + _P[0].x*_P[2].z*_P[3].y - _P[2].z*_P[0].x*_P[1].y - _P[0].x*_P[3].z*_P[2].y - _P[1].z*_P[3].y*_P[0].x
			- _P[3].z*_P[1].x*_P[0].y + _P[0].z*_P[3].y*_P[1].x - _P[0].z*_P[2].y*_P[1].x + _P[2].z*_P[1].x*_P[0].y + _P[3].z*_P[1].x*_P[2].y - _P[3].z*_P[1].y*_P[2].x + _P[3].z*_P[0].y*_P[2].x
			- _P[1].z*_P[3].x*_P[2].y - _P[1].z*_P[0].y*_P[2].x + _P[0].z*_P[3].x*_P[2].y - _P[0].z*_P[3].x*_P[1].y - _P[2].z*_P[1].x*_P[3].y + _P[2].z*_P[1].y*_P[3].x - _P[2].z*_P[0].y*_P[3].x
			+ _P[1].z*_P[2].x*_P[3].y + _P[1].z*_P[0].y*_P[3].x - _P[0].z*_P[2].x*_P[3].y + _P[0].z*_P[2].x*_P[1].y;
		detD_1 = abs(detD_1) / 6.;

		return detD_1;
	}
	catch (const std::exception&)
	{
		printf_s("Error: GeometryShape.cpp/geometry::Tetrahedron::GetVolume(Point[4] _P)");
	}
}

bool geometry::Tetrahedron::IsContainThePoint(Point<double> A)
{
	try
	{
		Point<double> P[4];
		double V = GetVolume(), V1, V2, V3, V4;

		P[0] = this->GetNode(0);
		P[1] = this->GetNode(1);
		P[2] = this->GetNode(2);
		P[3] = A;
		V1 = Volume(P);
		P[0] = this->GetNode(0);
		P[1] = this->GetNode(1);
		P[2] = A;
		P[3] = this->GetNode(3);
		V2 = Volume(P);
		P[0] = this->GetNode(0);
		P[1] = A;
		P[2] = this->GetNode(2);
		P[3] = this->GetNode(3);
		V3 = Volume(P);
		P[0] = A;
		P[1] = this->GetNode(1);
		P[2] = this->GetNode(2);
		P[3] = this->GetNode(3);
		V4 = Volume(P);

		double summ = V1 + V2 + V3 + V4;
		double res = round(abs(summ - V) / V, 5);
		if (res <= 1E-15)
			return true;
		return false;
	}
	catch (const std::exception&)
	{
		printf_s("Error: GeometryShape.cpp/geometry::Tetrahedron::IsContainThePoint(Point<double> A)");
	}
}

//==================Crack============================================
//geometry::Crack::Crack(char in_file_name[1000])
//{
//	this->Input(in_file_name);
}
//bool geometry::Crack::Input(char in_file_name[1000])
//{
//	try {
//		FILE *f;
//		fopen_s(&f, name, "r");
//
//		if (f == NULL) return false;
//				
//		int num_nodes;
//		fscanf_s(f, "%d", &num_nodes);
//		this->xyz.resize(num_nodes);
//		for (int i = 0; i < num_nodes; i++)
//		{
//			fscanf_s(f, "%lf %lf %lf", &this->xyz[i].x, &this->xyz[i].y, &this->xyz[i].z);
//		}
//
//		int num_planes;
//		fscanf_s(f, "%d", &num_planes);
//		this->triangles.resize(num_planes);
//		this->triangles_in_front.resize(num_planes);
//		for (int i = 0; i < num_planes; i++)
//		{
//			this->triangles_in_front[i] = false;
//			fscanf_s(f, "%d %d %d", 
//				&this->triangles[i].id_node[0],
//				&this->triangles[i].id_node[1], 
//				&this->triangles[i].id_node[2]);
//			this->triangles[i].node[0] = &(this->triangles.xyz[grid.geom_defects.cracks[num_cr].tr[i].id_node[0]]);
//			this->triangles[i].node[1] = &(this->triangles.xyz[grid.geom_defects.cracks[num_cr].tr[i].id_node[1]]);
//			this->triangles[i].node[2] = &(this->triangles.xyz[grid.geom_defects.cracks[num_cr].tr[i].id_node[2]]);
//		}
//
//		int num_fronts;
//		fscanf_s(f, "%d", &num_fronts);
//		this->fronts.resize(num_fronts);
//		this->points_in_fronts.resize(num_fronts);
//		for (int i = 0; i < num_fronts; i++)
//		{
//			int num_elem_in_front;
//			fscanf_s(f, "%d", &num_elem_in_front);
//			this->fronts[i].resize(num_elem_in_front);
//			for (int j = 0; j < num_elem_in_front; j++)
//			{
//				fscanf_s(f, "%d %d %d", &this->fronts[i][j].id_base_triangles,
//					&this->fronts[i][j].id_left, &this->fronts[i][j].id_right);
//
//				this->triangles_in_front[this->fronts[i][j].id_base_triangles] = true;
//			}
//
//			//closed front
//			if (this->fronts[i][0].id_left == this->fronts[i][num_elem_in_front - 1].id_right)
//			{
//				this->points_in_fronts[i].push_back(this->fronts[i][num_elem_in_front - 1].id_left);
//				for (int j = 0; j < num_elem_in_front; j++)
//				{
//					this->points_in_fronts[i].push_back(this->fronts[i][j].id_left);
//				}
//				this->points_in_fronts[i].push_back(this->fronts[i][num_elem_in_front - 1].id_right);
//				this->points_in_fronts[i].push_back(this->fronts[i][0].id_right);
//			}
//			else {
//				this->points_in_fronts[i].push_back(this->fronts[i][0].id_left);
//				for (int j = 0; j < num_elem_in_front; j++)
//				{
//					this->points_in_fronts[i].push_back(this->fronts[i][j].id_left);
//				}
//				this->points_in_fronts[i].push_back(this->fronts[i][num_elem_in_front - 1].id_right);
//				this->points_in_fronts[i].push_back(this->fronts[i][num_elem_in_front - 1].id_right);
//			}
//		}
//
//		fclose(f);
//		return true;
//	}
//	catch (const std::exception&)
//	{
//		printf_s("Error: GeometryShape.cpp/bool geometry::Crack::Input(char in_file_name[1000])");
//	}
//}
