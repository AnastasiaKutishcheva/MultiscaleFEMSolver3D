#pragma once
#include <vector>
#include "Point.h"

namespace topology
{
	namespace lower
	{
		template <class LowerElement>
		class Shape
		{
		public:
			Shape()
			{
				this->id_self = -1;
			};
			~Shape()
			{
				std::vector<LowerElement*> v_lower_elements;
				std::vector<LowerElement*>(v_lower_elements).swap(this->lower_elements);
			};

			void SetLowerElement(LowerElement *element)
			{
				try
				{
					this->lower_elements.push_back(element);
				}
				catch (const std::exception&)
				{
					printf_s("Error: TopologyShape.h/topology::Shape::SetLowerElement(LowerElement *element)\n");
				}
			}
			void SetLowerElement(int local_id, LowerElement *element)
			{
				try
				{
					this->lower_elements[local_id] = element;
				}
				catch (const std::exception&)
				{
					printf_s("Error: TopologyShape.h/topology::Shape::SetLowerElement(LowerElement *element)\n");
				}
			}
			void SetLowerElements(std::vector<LowerElement*> &element)
			{
				try
				{
					this->lower_elements.resize(element.size());
					for (int i = 0; i < element.size(); i++)
					{
						this->lower_elements[i] = element[i];
					}
				}
				catch (const std::exception&)
				{
					printf_s("Error: TopologyShape.h/topology::Shape::SetLowerElements(std::vector<LowerElement*> &element)\n");
				}
			}
			LowerElement* GetLowerElement(int number)
			{
				try
				{
					return this->lower_elements[number];
				}
				catch (const std::exception&)
				{
					printf_s("Error: TopologyShape.h/topology::Shape::GetLowerElement(int number)\n");
				}
			}
			int GetLowerElementCount()
			{
				try
				{
					return (int)this->lower_elements.size();
				}
				catch (const std::exception&)
				{
					printf_s("Error: TopologyShape.h/topology::Shape::GetLowerElementCount()\n");
				}
			}
			void SetLowerElementCount(int count)
			{
				try
				{
					this->lower_elements.resize(count);
				}
				catch (const std::exception&)
				{
					printf_s("Error: Library/TopologyShape.h/void SetLowerElementCount(int count)\n");
				}
			}

			//pattern[id lower element in local][id nodes in lower element in local]
			virtual void GetLowerElementPatternInLocal(std::vector<std::vector<int>> &local_pattern) = 0;

			void SetSelfGlobalId(int id)
			{
				this->id_self = id;
			}
			int GetSelfGlobalId()
			{
				return this->id_self;
			}

		private:
			std::vector<LowerElement*> lower_elements;
			int id_self; //in global grid numeration
		};

		class EmptyElement : public Shape <EmptyElement>
		{
		public:
			EmptyElement() {};
			~EmptyElement() {};

			void GetLowerElementPatternInLocal(std::vector<std::vector<int>> &local_pattern)
			{
				return;
			}
		};

		class Vertex : public Shape <EmptyElement>
			//class Vertex : public Shape <EmptyElement, Segment>
		{
		public:
			Vertex() {};
			~Vertex() {};

			void GetLowerElementPatternInLocal(std::vector<std::vector<int>> &local_pattern)
			{
				try
				{
					local_pattern.resize(1);
					local_pattern[0].resize(1);
					local_pattern[0][0] = 0;
					return;
				}
				catch (const std::exception&)
				{
					printf_s("Error: TopologyShape.h/topology::Node::GetLowerElementPatternInLocal(std::vector<std::vector<int>> &local_pattern)\n");
				}
			}
		};

		class Segment : public Shape <Vertex>
		{
		public:
			Segment()
			{
				this->SetLowerElementCount(2);
			};
			~Segment() {};

			void GetLowerElementPatternInLocal(std::vector<std::vector<int>> &local_pattern)
			{
				try
				{
					local_pattern.resize(2);
					local_pattern[0].resize(1);
					local_pattern[1].resize(1);
					local_pattern[0][0] = 0;
					local_pattern[1][0] = 1;
					return;
				}
				catch (const std::exception&)
				{
					printf_s("Error: TopologyShape.h/topology::Segment::GetLowerElementPatternInLocal(std::vector<std::vector<int>> &local_pattern)\n");
				}
			}
		};

		class Polyline : public Shape <Vertex>
		{
		public:
			Polyline() {
				this->SetLowerElementCount(2);
			};
			~Polyline() {};

			void GetLowerElementPatternInLocal(std::vector<std::vector<int>> &local_pattern)
			{
				try
				{
					local_pattern.resize(2);
					local_pattern[0].resize(1);
					local_pattern[1].resize(1);
					local_pattern[0][0] = 0;
					local_pattern[1][0] = 1;
					return;
				}
				catch (const std::exception&)
				{
					printf_s("Error: TopologyShape.h/topology::Polyline::GetLowerElementPatternInLocal(std::vector<std::vector<int>> &local_pattern)\n");
				}
			}
		};

		class Triangle : public Shape <Segment>
		{
		public:
			Triangle() {
				this->SetLowerElementCount(3);
			};
			~Triangle() {};

			void GetLowerElementPatternInLocal(std::vector<std::vector<int>> &local_pattern)
			{
				try
				{
					local_pattern.resize(3);
					local_pattern[0].resize(2);
					local_pattern[1].resize(2);
					local_pattern[2].resize(2);
					local_pattern[0][0] = 0; local_pattern[0][1] = 1;
					local_pattern[1][0] = 0; local_pattern[1][1] = 2;
					local_pattern[2][0] = 1; local_pattern[2][1] = 2;
					return;
				}
				catch (const std::exception&)
				{
					printf_s("Error: TopologyShape.h/topology::Triangle::GetLowerElementPatternInLocal(std::vector<std::vector<int>> &local_pattern)\n");
				}
			}
		};

		class TriangleToVertex : public Shape <Vertex>
		{
		public:
			TriangleToVertex() {
				this->SetLowerElementCount(3);
			};
			~TriangleToVertex() {};

			void GetLowerElementPatternInLocal(std::vector<std::vector<int>>& local_pattern)
			{
				try
				{
					local_pattern.resize(3);
					local_pattern[0].resize(1);
					local_pattern[1].resize(1);
					local_pattern[2].resize(1);
					local_pattern[0][0] = 0;
					local_pattern[1][0] = 1;
					local_pattern[2][0] = 2;
					return;
				}
				catch (const std::exception&)
				{
					printf_s("Error: TopologyShape.h/topology::Tetrahedron::GetLowerElementPatternInLocal(std::vector<std::vector<int>> &local_pattern)\n");
				}
			}
		};

		class Polygon : public Shape <Polyline>
		{
		public:
			Polygon() {};
			~Polygon() {};

			void GetLowerElementPatternInLocal(std::vector<std::vector<int>> &local_pattern)
			{
				try
				{
					local_pattern.resize(this->GetLowerElementCount());
					for (int i = 0; i < local_pattern.size(); i++)
					{
						local_pattern[i].resize(2);
						local_pattern[i][0] = i; local_pattern[i][1] = i + 1;
					}
					local_pattern[local_pattern.size() - 1][1] = 0;
					return;
				}
				catch (const std::exception&)
				{
					printf_s("Error: TopologyShape.h/topology::Polygon::GetLowerElementPatternInLocal(std::vector<std::vector<int>> &local_pattern)\n");
				}
			}
		};

		class Tetrahedron : public Shape <Triangle>
		{
		public:
			Tetrahedron() {
				this->SetLowerElementCount(4);
			};
			~Tetrahedron() {};

			void GetLowerElementPatternInLocal(std::vector<std::vector<int>> &local_pattern)
			{
				try
				{
					local_pattern.resize(4);
					local_pattern[0].resize(3);
					local_pattern[1].resize(3);
					local_pattern[2].resize(3);
					local_pattern[3].resize(3);
					local_pattern[0][0] = 0; local_pattern[0][1] = 1; local_pattern[0][2] = 2;
					local_pattern[1][0] = 0; local_pattern[1][1] = 1; local_pattern[1][2] = 3;
					local_pattern[2][0] = 0; local_pattern[2][1] = 2; local_pattern[2][2] = 3;
					local_pattern[3][0] = 1; local_pattern[3][1] = 2; local_pattern[3][2] = 3;
					return;
				}
				catch (const std::exception&)
				{
					printf_s("Error: TopologyShape.h/topology::Tetrahedron::GetLowerElementPatternInLocal(std::vector<std::vector<int>> &local_pattern)\n");
				}
			}
		};

		class TetrahedronToVertex : public Shape <Vertex>
		{
		public:
			TetrahedronToVertex() {
				this->SetLowerElementCount(4);
			};
			~TetrahedronToVertex() {};

			void GetLowerElementPatternInLocal(std::vector<std::vector<int>> &local_pattern)
			{
				try
				{
					local_pattern.resize(4);
					local_pattern[0].resize(1);
					local_pattern[1].resize(1);
					local_pattern[2].resize(1);
					local_pattern[3].resize(1);
					local_pattern[0][0] = 0;
					local_pattern[1][0] = 1;
					local_pattern[2][0] = 2;
					local_pattern[3][0] = 3;
					return;
				}
				catch (const std::exception&)
				{
					printf_s("Error: TopologyShape.h/topology::Tetrahedron::GetLowerElementPatternInLocal(std::vector<std::vector<int>> &local_pattern)\n");
				}
			}
		};

		class TetrahedronToEdges : public Shape <Segment>
		{
		public:
			TetrahedronToEdges() {
				this->SetLowerElementCount(6);
			};
			~TetrahedronToEdges() {};

			void GetLowerElementPatternInLocal(std::vector<std::vector<int>> &local_pattern)
			{
				try
				{
					local_pattern.resize(6);
					local_pattern[0].resize(2);
					local_pattern[1].resize(2);
					local_pattern[2].resize(2);
					local_pattern[3].resize(2);
					local_pattern[4].resize(2);
					local_pattern[5].resize(2);
					local_pattern[0][0] = 0; local_pattern[0][1] = 1;
					local_pattern[1][0] = 0; local_pattern[1][1] = 2;
					local_pattern[2][0] = 0; local_pattern[2][1] = 3;
					local_pattern[3][0] = 1; local_pattern[3][1] = 2;
					local_pattern[4][0] = 1; local_pattern[4][1] = 3;
					local_pattern[5][0] = 2; local_pattern[5][1] = 3;

					return;
				}
				catch (const std::exception&)
				{
					printf_s("Error: TopologyShape.h/topology::Tetrahedron::GetLowerElementPatternInLocal(std::vector<std::vector<int>> &local_pattern)\n");
				}
			}
		};

		class Polyhedron : public Shape <Polygon>
		{
		public:
			Polyhedron() {};
			~Polyhedron() {};

			//This realiation isn't right!
			void GetLowerElementPatternInLocal(std::vector<std::vector<int>> &local_pattern)
			{
				try
				{
					local_pattern.resize(this->GetLowerElementCount());
					for (int i = 0; i < local_pattern.size(); i++)
					{
						auto lower_elem = this->GetLowerElement(i);
						local_pattern[i].resize(lower_elem->GetLowerElementCount());
						for (int j = 0; j < local_pattern[i].size(); j++)
						{
							local_pattern[i][j] = i; local_pattern[i][1] = i + 1;
						}
					}
					local_pattern[local_pattern.size() - 1][1] = 0;
					return;
				}
				catch (const std::exception&)
				{
					printf_s("Error: TopologyShape.h/topology::Polyhedron::GetLowerElementPatternInLocal(std::vector<std::vector<int>> &local_pattern)\n");
				}
			}
		};
	}
	namespace upper {
		template <class UpperElement>
		class Shape
		{
		public:
			Shape()
			{
				this->id_self = -1;
			};
			~Shape()
			{
				std::vector<UpperElement*> v_upper_elements;
				std::vector<UpperElement*>(v_upper_elements).swap(this->upper_elements);
			};

			
			void SetUpperElement(int local_id, UpperElement *element)
			{
				try
				{
					this->upper_elements[local_id] = element;
				}
				catch (const std::exception&)
				{
					printf_s("Error: TopologyShape.h/void SetUpperElement(int local_id, UpperElement *element)\n");
				}
			}
			void SetUpperElement(UpperElement *element)
			{
				try
				{
					this->upper_elements.push_back(element);
				}
				catch (const std::exception&)
				{
					printf_s("Error: TopologyShape.h/topology::Shape::SetUpperElement(UpperElement *element)\n");
				}
			}
			void SetUpperElements(std::vector<UpperElement*> &element)
			{
				try
				{
					this->upper_elements.resize(element.size());
					for (int i = 0; i < element.size(); i++)
					{
						this->upper_elements[i] = element[i];
					}
				}
				catch (const std::exception&)
				{
					printf_s("Error: TopologyShape.h/topology::Shape::SetUpperElements(std::vector<UpperElement*> &element)\n");
				}
			}
			UpperElement* GetUpperElement(int number)
			{
				try
				{
					return this->upper_elements[number];
				}
				catch (const std::exception&)
				{
					printf_s("Error: TopologyShape.h/topology::Shape::GetUpperElement(int number)\n");
				}
			}
			int GetUpperElementCount()
			{
				try
				{
					return (int)this->upper_elements.size();
				}
				catch (const std::exception&)
				{
					printf_s("Error: TopologyShape.h/topology::Shape::GetUpperElementCount()\n");
				}
			}
			void SetUpperElementCount(int count)
			{
				try
				{
					this->upper_elements.resize(count);
				}
				catch (const std::exception&)
				{
					printf_s("Error: Library/TopologyShape.h/void SetUpperElementCount(int count)\n");
				}
			}
			
			void SetSelfGlobalId(int id)
			{
				this->id_self = id;
			}
			int GetSelfGlobalId()
			{
				return this->id_self;
			}

		private:
			std::vector<UpperElement*> upper_elements;
			int id_self; //in global grid numeration
		};

		class EmptyElement : public Shape <EmptyElement>
		{
		public:
			EmptyElement() {};
			~EmptyElement() {};

		};
		
		class Polyhedron : public Shape <EmptyElement>
		{
		public:
			Polyhedron() {};
			~Polyhedron() {};
		};

		class Tetrahedron : public Shape <EmptyElement>
		{
		public:
			Tetrahedron() {
			};
			~Tetrahedron() {};

		};

		class Polygon : public Shape <Polyhedron>
		{
		public:
			Polygon() {};
			~Polygon() {};
		};

		class Triangle : public Shape <Tetrahedron>
		{
		public:
			Triangle() {
			};
			~Triangle() {};

		};

		class TriangleUpper : public Shape <EmptyElement>
		{
		public:
			TriangleUpper() {
			};
			~TriangleUpper() {};

		};

		class Polyline : public Shape <Polygon>
		{
		public:
			Polyline() {
			};
			~Polyline() {};
		};

		class Segment : public Shape <Triangle>
		{
		public:
			Segment()
			{
			};
			~Segment() {};

		};
		class SegmentToTetrahedron : public Shape <Tetrahedron>
		{
		public:
			SegmentToTetrahedron()
			{
			};
			~SegmentToTetrahedron() {};

		};
		class SegmentUpper : public Shape <EmptyElement>
		{
		public:
			SegmentUpper()
			{
			};
			~SegmentUpper() {};

		};

		class Vertex : public Shape <Segment>
		{
		public:
			Vertex() {};
			~Vertex() {};
		};

		class VertexToTetrahedron : public Shape <Tetrahedron>
		{
		public:
			VertexToTetrahedron() {};
			~VertexToTetrahedron() {};
		};
	}
	
	template <class LowerElement, class UpperElement>
	class Shape : public LowerElement, public UpperElement
	{
	public:
		void SetSelfGlobalId(int id)
		{
			LowerElement::SetSelfGlobalId(id);
			UpperElement::SetSelfGlobalId(id);
			return;
		}
		int GetSelfGlobalId()
		{
			return LowerElement::GetSelfGlobalId();
		}
		
	};

	template <class LowerElement, class UpperElement>
	class Polyhedron : public Shape <LowerElement, UpperElement>
	{
	public:
		Polyhedron()
		{
		}; 
		Polyhedron(int num_of_low_elements)
		{
			this->SetLowerElementCount(num_of_low_elements);
		};
	};


	template <class LowerElement, class UpperElement>
	class Tetrahedron : public Shape <LowerElement, UpperElement>
	{
	public:
		Tetrahedron()
		{
			this->SetLowerElementCount(4);
		};
	};

	template <class LowerElement, class UpperElement>
	class Polygon : public Shape <LowerElement, UpperElement>
	{
	public:
		Polygon()
		{
		};
		Polygon(int num_of_low_elements)
		{
			this->SetLowerElementCount(num_of_low_elements);
		};
	};

	template <class LowerElement, class UpperElement>
	class Triangle : public Shape <LowerElement, UpperElement>
	{
	public:
		Triangle()
		{
			this->SetLowerElementCount(3);
		};
	};

	template <class LowerElement, class UpperElement>
	class Polyline : public Shape <LowerElement, UpperElement>
	{
	public:
		Polyline()
		{
			this->SetLowerElementCount(2);
		};
	};

	template <class LowerElement, class UpperElement>
	class Segment : public Shape <LowerElement, UpperElement>
	{
	public:
		Segment()
		{
			this->SetLowerElementCount(2);
		};
	};

	template <class LowerElement, class UpperElement>
	class Vertex : public Shape <LowerElement, UpperElement>
	{
	public:
		Vertex()
		{
			this->SetLowerElementCount(0);
		};
	};

	template <class LowerElement, class UpperElement>
	class EmptyElement : public Shape <LowerElement, UpperElement>
	{
	};
}