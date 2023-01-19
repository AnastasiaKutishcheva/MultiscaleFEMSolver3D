#include "../FEM/FEM.h"
#include <Windows.h>
#include <iostream>

struct integrate_for_rpho_eff
{
	Point<double> integr_EE_full, integr_JJ_full;
	std::vector<Point<double>> integr_in_domain_EE, integr_in_domain_JJ;
	std::vector<double> V_in_domain;
};

template <typename Grid>
double EffectiveResistance_new_LIGHT(Grid& grid,
	std::vector<double>& U, 
	double& norm_Jsqrt, 
	double& rpho_eff_summ, 
	double& GLOBAL_SIGMA_hmm, 
	integrate_for_rpho_eff& result_of_integration)
{
	//вычисляем сразу ||J|| и ||gardU||
	double norm_J = 0;
	norm_Jsqrt = 0;
	double norm_gradU = 0;

	GLOBAL_SIGMA_hmm = 0;
	double norm_gradUsqrt_on_points = 0;
	double norm_Jsqrt_on_points = 0;


	double V_global = 0;
	for (int elem_coarse = 0; elem_coarse < grid.GetElementsCount(); elem_coarse++)
	{
		auto element = grid.GetElement(elem_coarse);
		V_global += element->GetVolume();
		result_of_integration.V_in_domain[element->GetIdDomain()] += element->GetVolume();
	}

	std::vector<double> macro_norm_J_master;
	std::vector<double> macro_norm_Jsqrt_master;
	std::vector<double> macro_norm_gradU_master;
	std::vector<double> macro_norm_Jsqrt_on_points_master;
	std::vector<double> macro_norm_gradUsqrt_on_points_master;
	for (int elem_coarse = 0; elem_coarse < grid.GetElementsCount(); elem_coarse++)
	{
		auto element = grid.GetElement(elem_coarse);
		
		std::vector<double> value_in_local_nodes(element->GetDOFsCount());
		for (int i = 0; i < element->GetDOFsCount(); i++)
			value_in_local_nodes[i] = U[element->GetDOFInLocalID(i)];

		if (true)
		{
			std::vector<Point<double>> x;
			std::vector<double> w;
			element->SetIntegrationLaw(4);
			element->GetIntegrationProperties(w, x);
			for (int ii = 0; ii < w.size(); ii++)
			{
				std::function<double(int, Point<double>)> StiffnessCoef = [&](int elem, Point<double> x) ->double
				{
					return grid.GetDomain(element->GetIdDomain())->forElectrical.sigma;
				};
				Point<double> E, J;
				for (int k = 0; k < element->GetDOFsCount(); k++)
				{
					auto BF_d = (*element->GetDerivativeOfBasisFunctionInLocalID(k))(x[ii]);
					E.x += value_in_local_nodes[k] * BF_d.x;
					E.y += value_in_local_nodes[k] * BF_d.y;
					E.z += value_in_local_nodes[k] * BF_d.z;
				}
				double sigma;
				sigma = StiffnessCoef(elem_coarse, x[ii]);
				J = Point<double>(sigma * E.x, sigma * E.y, sigma * E.z);

				result_of_integration.integr_EE_full.x += w[ii] * E.x * E.x;
				result_of_integration.integr_EE_full.y += w[ii] * E.y * E.y;
				result_of_integration.integr_EE_full.z += w[ii] * E.z * E.z;
				result_of_integration.integr_JJ_full.x += w[ii] * J.x * J.x;
				result_of_integration.integr_JJ_full.y += w[ii] * J.y * J.y;
				result_of_integration.integr_JJ_full.z += w[ii] * J.z * J.z;

				result_of_integration.integr_in_domain_EE[element->GetIdDomain()].x += w[ii] * E.x * E.x;
				result_of_integration.integr_in_domain_EE[element->GetIdDomain()].y += w[ii] * E.y * E.y;
				result_of_integration.integr_in_domain_EE[element->GetIdDomain()].z += w[ii] * E.z * E.z;

				result_of_integration.integr_in_domain_JJ[element->GetIdDomain()].x += w[ii] * J.x * J.x;
				result_of_integration.integr_in_domain_JJ[element->GetIdDomain()].y += w[ii] * J.y * J.y;
				result_of_integration.integr_in_domain_JJ[element->GetIdDomain()].z += w[ii] * J.z * J.z;

				norm_Jsqrt += w[ii] * sqrt(J.x * J.x + J.y * J.y + J.z * J.z);
				norm_J += w[ii] * (J.x * J.x + J.y * J.y + J.z * J.z);
				norm_gradU += w[ii] * (E.x * E.x + E.y * E.y + E.z * E.z);
				//norm_Jsqrt_on_points += tmpJsqrt_on_points;
				//norm_gradUsqrt_on_points += tmpUsqrt_on_points;
			};
		}
	}

	norm_J = sqrt(abs(norm_J));
	printf("\n----------->\nnorm_J = %.15lf\n<-----------\n", norm_J);
	norm_gradU = sqrt(abs(norm_gradU));
	printf("\n----------->\nnorm_gradU = %.15lf\n<-----------\n", norm_gradU);

	printf("\n----------->\nI = %.15lf\n<-----------\n", norm_Jsqrt);

	rpho_eff_summ = norm_gradUsqrt_on_points / norm_Jsqrt_on_points;

	return (norm_gradU / norm_J);

}

/// <summary>
/// 
/// </summary>
/// <typeparam name="Grid">type of solver grid (usualy FEM::Grid_forScal)</typeparam>
/// <param name="typeXYZ"> 0 - X, 1 - Y, 2 - Z</param>
/// <param name="typeGEO"> 0 - cylinder, 1 - box</param>
/// <param name="SectionPositions">(0;1) - position of solver sections</param>
/// <param name="solver_grid">full solver grid</param>
/// <param name="U">solution of the direct problem</param>
/// <param name="I_in_sections">total current in cross section</param>
/// <returns>effective resistivity (Om*m)</returns>
template <typename Grid>
double EffectiveResistance_ByСrossSection(int typeXYZ, int typeGEO, double deltaU, std::vector<double> &SectionPositions, Grid &solver_grid, std::vector<double> &U, std::vector<double> &I_in_sections)
{
	//create cross section
	math::SimpleGrid section_grid;
	Point<double> section_centr;
	double sizeX, sizeY;
	auto Min = solver_grid.GetMinCoordinate();
	auto Max = solver_grid.GetMaxCoordinate();

	//получение параметров
	switch (typeXYZ)
	{
	case 0:
		section_centr.x = (Min.y + Max.y) / 2.0;
		section_centr.y = (Min.z + Max.z) / 2.0;
		section_centr.z = (Min.x);
		sizeX = Max.y - Min.y;
		sizeY = Max.z - Min.z;
		break;

	case 1:
		section_centr.x = (Min.x + Max.x) / 2.0;
		section_centr.y = (Min.z + Max.z) / 2.0;
		section_centr.z = Min.y;
		sizeX = Max.x - Min.x;
		sizeY = Max.z - Min.z;
		break;

	case 2:
		section_centr.x = (Min.x + Max.x) / 2.0;
		section_centr.y = (Min.y + Max.y) / 2.0;
		section_centr.z = Min.z;
		sizeX = Max.x - Min.x;
		sizeY = Max.y - Min.y;
		break;
	default:
		break;
	}
	//построение сетки
	switch (typeGEO)
	{
	case 0:
		section_grid.CreateGrid2D_circle(section_centr, sizeX/2.0, sizeY/2.0, 40);
		break;
	case 1:
		section_grid.CreateGrid2D_rectangle(
			Point<double>(section_centr.x - sizeX / 2.0, section_centr.y - sizeY / 2.0, section_centr.z),
			Point<double>(section_centr.x + sizeX / 2.0, section_centr.y + sizeY / 2.0, section_centr.z),
			20, 20);
		break;
	}
	//разворот сетки
	switch (typeXYZ)
	{
	case 0:
		for (int i = 0; i < section_grid.xyz.size(); i++)
		{
			Point<double> rotation_point;
			rotation_point.x = section_grid.xyz[i].z;
			rotation_point.y = section_grid.xyz[i].x;
			rotation_point.z = section_grid.xyz[i].y;

			section_grid.xyz[i] = rotation_point;
		}
		break;

	case 1:
		section_centr.x = (Min.x + Max.x) / 2.0;
		section_centr.y = (Min.z + Max.z) / 2.0;
		section_centr.z = Min.y;
		sizeX = Max.x - Min.x;
		sizeY = Max.z - Min.z;
		for (int i = 0; i < section_grid.xyz.size(); i++)
		{
			Point<double> rotation_point;
			rotation_point.x = section_grid.xyz[i].x;
			rotation_point.y = section_grid.xyz[i].z;
			rotation_point.z = section_grid.xyz[i].y;

			section_grid.xyz[i] = rotation_point;
		}
		break;
	}


	//считаем интеграллы в сечениях
	double I_middle = 0;
	double S = 0; //площадь сечения
	I_in_sections.resize(SectionPositions.size());
	omp_set_num_threads(math::NUM_THREADS);
#pragma omp parallel for schedule(dynamic)
	for (int id_section = 0; id_section < SectionPositions.size(); id_section++)
	{
		//подвижка сетки на нужную позицию
		double new_x;
		switch (typeXYZ)
		{
		case 0:
			new_x = Min.x + (Max.x - Min.x) * SectionPositions[id_section];
			for (int i = 0; i < section_grid.xyz.size(); i++)
			{
				section_grid.xyz[i].x = new_x;
			}

			break;

		case 1:
			new_x = Min.y + (Max.y - Min.y) * SectionPositions[id_section];
			for (int i = 0; i < section_grid.xyz.size(); i++)
			{
				section_grid.xyz[i].y = new_x;
			}

			break;

		case 2:
			new_x = Min.z + (Max.z - Min.z) * SectionPositions[id_section];
			for (int i = 0; i < section_grid.xyz.size(); i++)
			{
				section_grid.xyz[i].z = new_x;
			}
			break;
		}

		for (int id_tr = 0; id_tr < section_grid.nvtr.size(); id_tr++)
		{
			geometry::Triangle triangle(
				section_grid.xyz[section_grid.nvtr[id_tr][0]],
				section_grid.xyz[section_grid.nvtr[id_tr][1]],
				section_grid.xyz[section_grid.nvtr[id_tr][2]]);

			std::function<double(Point<double>)> I = [&solver_grid, &U, typeXYZ](Point<double> x) ->double
			{
				double res = 0;

				double len;
				int id_3Delem = solver_grid.GetNearestElementID(x, len);

				if (id_3Delem >= 0)
				{
					Point<double> E = solver_grid.GetDerevativeFromSolutionInPoint(id_3Delem, x, U);
					double sigma = solver_grid.GetDomain(solver_grid.GetElement(id_3Delem)->GetIdDomain())->forElectrical.sigma;

					switch (typeXYZ)
					{
					case 0:
						res = sigma * E.x;
						break;
					case 1:
						res = sigma * E.y;
						break;
					case 2:
						res = sigma * E.z;
						break;
					}
				}
				return res;
			};

			triangle.SetIntegrationLaw(3);
			I_in_sections[id_section] += triangle.SolveIntegral(I);

			if (id_section == 0)
				S += triangle.GetVolume();
		}

		//проверка сеток для сечений
		FILE* fout_tech;
		char name_u_tech[5000];
		sprintf_s(name_u_tech, "./section_%d.dat", id_section);
		fopen_s(&fout_tech, name_u_tech, "w");
		char name_in_file[1000];
		sprintf_s(name_in_file, "section_%d", id_section);
		section_grid.printTecPlot3D(fout_tech, name_in_file);
	}

	double L; //расстояние между электродами
	switch (typeXYZ)
	{
	case 0:
		L = Max.x - Min.x;
		break;
	case 1:
		L = Max.y - Min.y;
		break;
	case 2:
		L = Max.z - Min.z;
		break;
	}
	for (int id_section = 0; id_section < SectionPositions.size(); id_section++)
		I_middle += I_in_sections[id_section] / I_in_sections.size();

	double rpho_eff = S * deltaU / (I_middle * L);

	return rpho_eff;
}
template <typename Grid>
double EffectiveResistance_ByСrossSection_forTensorDomains(int typeXYZ, int typeGEO, double deltaU, std::vector<double>& SectionPositions, Grid& solver_grid, std::vector<double>& U, std::vector<double>& I_in_sections)
{
	//create cross section
	math::SimpleGrid section_grid;
	Point<double> section_centr;
	double sizeX, sizeY;
	auto Min = solver_grid.GetMinCoordinate();
	auto Max = solver_grid.GetMaxCoordinate();

	//получение параметров
	switch (typeXYZ)
	{
	case 0:
		section_centr.x = (Min.y + Max.y) / 2.0;
		section_centr.y = (Min.z + Max.z) / 2.0;
		section_centr.z = (Min.x);
		sizeX = Max.y - Min.y;
		sizeY = Max.z - Min.z;
		break;

	case 1:
		section_centr.x = (Min.x + Max.x) / 2.0;
		section_centr.y = (Min.z + Max.z) / 2.0;
		section_centr.z = Min.y;
		sizeX = Max.x - Min.x;
		sizeY = Max.z - Min.z;
		break;

	case 2:
		section_centr.x = (Min.x + Max.x) / 2.0;
		section_centr.y = (Min.y + Max.y) / 2.0;
		section_centr.z = Min.z;
		sizeX = Max.x - Min.x;
		sizeY = Max.y - Min.y;
		break;
	default:
		break;
	}
	//построение сетки
	switch (typeGEO)
	{
	case 0:
		section_grid.CreateGrid2D_circle(section_centr, sizeX / 2.0, sizeY / 2.0, 40);
		break;
	case 1:
		section_grid.CreateGrid2D_rectangle(
			Point<double>(section_centr.x - sizeX / 2.0, section_centr.y - sizeY / 2.0, section_centr.z),
			Point<double>(section_centr.x + sizeX / 2.0, section_centr.y + sizeY / 2.0, section_centr.z),
			20, 20);
		break;
	}
	//разворот сетки
	switch (typeXYZ)
	{
	case 0:
		for (int i = 0; i < section_grid.xyz.size(); i++)
		{
			Point<double> rotation_point;
			rotation_point.x = section_grid.xyz[i].z;
			rotation_point.y = section_grid.xyz[i].x;
			rotation_point.z = section_grid.xyz[i].y;

			section_grid.xyz[i] = rotation_point;
		}
		break;

	case 1:
		section_centr.x = (Min.x + Max.x) / 2.0;
		section_centr.y = (Min.z + Max.z) / 2.0;
		section_centr.z = Min.y;
		sizeX = Max.x - Min.x;
		sizeY = Max.z - Min.z;
		for (int i = 0; i < section_grid.xyz.size(); i++)
		{
			Point<double> rotation_point;
			rotation_point.x = section_grid.xyz[i].x;
			rotation_point.y = section_grid.xyz[i].z;
			rotation_point.z = section_grid.xyz[i].y;

			section_grid.xyz[i] = rotation_point;
		}
		break;
	}


	//считаем интеграллы в сечениях
	double I_middle = 0;
	double S = 0; //площадь сечения
	I_in_sections.resize(SectionPositions.size());
	omp_set_num_threads(math::NUM_THREADS);
#pragma omp parallel for schedule(dynamic)
	for (int id_section = 0; id_section < SectionPositions.size(); id_section++)
	{
		//подвижка сетки на нужную позицию
		double new_x;
		switch (typeXYZ)
		{
		case 0:
			new_x = Min.x + (Max.x - Min.x) * SectionPositions[id_section];
			for (int i = 0; i < section_grid.xyz.size(); i++)
			{
				section_grid.xyz[i].x = new_x;
			}

			break;

		case 1:
			new_x = Min.y + (Max.y - Min.y) * SectionPositions[id_section];
			for (int i = 0; i < section_grid.xyz.size(); i++)
			{
				section_grid.xyz[i].y = new_x;
			}

			break;

		case 2:
			new_x = Min.z + (Max.z - Min.z) * SectionPositions[id_section];
			for (int i = 0; i < section_grid.xyz.size(); i++)
			{
				section_grid.xyz[i].z = new_x;
			}
			break;
		}

		for (int id_tr = 0; id_tr < section_grid.nvtr.size(); id_tr++)
		{
			geometry::Triangle triangle(
				section_grid.xyz[section_grid.nvtr[id_tr][0]],
				section_grid.xyz[section_grid.nvtr[id_tr][1]],
				section_grid.xyz[section_grid.nvtr[id_tr][2]]);

			std::function<double(Point<double>)> I = [&solver_grid, &U, typeXYZ](Point<double> x) ->double
			{
				double res;

				double len;
				int id_3Delem = solver_grid.GetNearestElementID(x, len);

				Point<double> E = solver_grid.GetDerevativeFromSolutionInPoint(id_3Delem, x, U);
				Tensor2Rank3D sigma = solver_grid.GetDomain(solver_grid.GetElement(id_3Delem)->GetIdDomain())->forElectrical.sigma_tensor;

				Point<double> J = sigma * E;
				switch (typeXYZ)
				{
				case 0:
					res = J.x;
					break;
				case 1:
					res = J.y;
					break;
				case 2:
					res = J.z;
					break;
				}

				return res;
			};

			triangle.SetIntegrationLaw(3);
			I_in_sections[id_section] += triangle.SolveIntegral(I);

			if (id_section == 0)
				S += triangle.GetVolume();
		}

		//проверка сеток для сечений
		FILE* fout_tech;
		char name_u_tech[5000];
		sprintf_s(name_u_tech, "./section_%d.dat", id_section);
		fopen_s(&fout_tech, name_u_tech, "w");
		char name_in_file[1000];
		sprintf_s(name_in_file, "section_%d", id_section);
		section_grid.printTecPlot3D(fout_tech, name_in_file);
	}

	double L; //расстояние между электродами
	switch (typeXYZ)
	{
	case 0:
		L = Max.x - Min.x;
		break;
	case 1:
		L = Max.y - Min.y;
		break;
	case 2:
		L = Max.z - Min.z;
		break;
	}
	for (int id_section = 0; id_section < SectionPositions.size(); id_section++)
		I_middle += I_in_sections[id_section] / I_in_sections.size();

	double rpho_eff = S * deltaU / (I_middle * L);

	return rpho_eff;
}

template <typename Grid>
double EffectiveResistance_ByBoundary(int typeXYZ,
	std::vector<std::vector<std::vector<int>>>& boundaries,
	std::vector<int>& target_types,
	double deltaU,
	Grid& solver_grid,
	std::vector<double>& U,
	std::vector<double>& I_in_sections,
	std::vector<double>& S_in_sections)
{
	//считаем интеграллы в сечениях
	double I_middle = 0;
	S_in_sections.resize(2); //площадь сечения
	I_in_sections.resize(2);
	auto Min = solver_grid.GetMinCoordinate();
	auto Max = solver_grid.GetMaxCoordinate();

	auto SolveI = [&solver_grid, &U, typeXYZ](std::vector<std::vector<int>>& boundary) -> double
	{
		double result_I = 0;
		for (int id_tr = 0; id_tr < boundary.size(); id_tr++)
		{
			if (boundary[id_tr][0] >= 0)
			{
				geometry::Triangle triangle(
					solver_grid.GetCoordinateViaID(boundary[id_tr][1]),
					solver_grid.GetCoordinateViaID(boundary[id_tr][2]),
					solver_grid.GetCoordinateViaID(boundary[id_tr][3]));

				int id_base_tetr = boundary[id_tr][0];

				std::function<double(Point<double>)> I = [&solver_grid, &U, typeXYZ, id_base_tetr](Point<double> x) ->double
				{
					double res;

					Point<double> E = solver_grid.GetDerevativeFromSolutionInPoint(id_base_tetr, x, U);
					double sigma = solver_grid.GetDomain(solver_grid.GetElement(id_base_tetr)->GetIdDomain())->forElectrical.sigma;

					switch (typeXYZ)
					{
					case 0:
						res = sigma * E.x;
						break;
					case 1:
						res = sigma * E.y;
						break;
					case 2:
						res = sigma * E.z;
						break;
					}

					return res;
				};

				triangle.SetIntegrationLaw(1);
				result_I += triangle.SolveIntegral(I);
			}
		}
		return result_I;
	};
	auto SolveS = [&solver_grid](std::vector<std::vector<int>>& boundary) -> double
	{
		double result_S = 0;
		for (int id_tr = 0; id_tr < boundary.size(); id_tr++)
		{
			if (boundary[id_tr][0] >= 0)
			{
				geometry::Triangle triangle(
					solver_grid.GetCoordinateViaID(boundary[id_tr][1]),
					solver_grid.GetCoordinateViaID(boundary[id_tr][2]),
					solver_grid.GetCoordinateViaID(boundary[id_tr][3]));
				result_S += triangle.GetVolume();
			}
		}
		return result_S;
	};

	S_in_sections[0] = SolveS(boundaries[target_types[0]]);
	I_in_sections[0] = SolveI(boundaries[target_types[0]]);
	S_in_sections[1] = SolveS(boundaries[target_types[1]]);
	I_in_sections[1] = SolveI(boundaries[target_types[1]]);


	double L; //расстояние между электродами
	switch (typeXYZ)
	{
	case 0:
		L = Max.x - Min.x;
		break;
	case 1:
		L = Max.y - Min.y;
		break;
	case 2:
		L = Max.z - Min.z;
		break;
	}
	
	double rpho_eff_0 = S_in_sections[0] * deltaU / (I_in_sections[0] * L);
	double rpho_eff_1 = S_in_sections[1] * deltaU / (I_in_sections[1] * L);

	return (rpho_eff_0 + rpho_eff_1)/2.0;
}
template <typename Grid>
double EffectiveResistance_ByBoundary_forTensorDomains(int typeXYZ,
	std::vector<std::vector<std::vector<int>>>& boundaries,
	std::vector<int>& target_types,
	double deltaU,
	Grid& solver_grid,
	std::vector<double>& U,
	std::vector<double>& I_in_sections,
	std::vector<double>& S_in_sections)
{
	//считаем интеграллы в сечениях
	double I_middle = 0;
	S_in_sections.resize(2); //площадь сечения
	I_in_sections.resize(2);
	auto Min = solver_grid.GetMinCoordinate();
	auto Max = solver_grid.GetMaxCoordinate();

	auto SolveI = [&solver_grid, &U, typeXYZ](std::vector<std::vector<int>>& boundary) -> double
	{
		double result_I = 0;
		for (int id_tr = 0; id_tr < boundary.size(); id_tr++)
		{
			geometry::Triangle triangle(
				solver_grid.GetCoordinateViaID(boundary[id_tr][1]),
				solver_grid.GetCoordinateViaID(boundary[id_tr][2]),
				solver_grid.GetCoordinateViaID(boundary[id_tr][3]));

			int id_base_tetr = boundary[id_tr][0];

			std::function<double(Point<double>)> I = [&solver_grid, &U, typeXYZ, id_base_tetr](Point<double> x) ->double
			{
				double res;

				Point<double> E = solver_grid.GetDerevativeFromSolutionInPoint(id_base_tetr, x, U);
				Tensor2Rank3D sigma = solver_grid.GetDomain(solver_grid.GetElement(id_base_tetr)->GetIdDomain())->forElectrical.sigma_tensor;
				Point<double> J = sigma * E;

				switch (typeXYZ)
				{
				case 0:
					res = J.x;
					break;
				case 1:
					res = J.y;
					break;
				case 2:
					res = J.z;
					break;
				}

				return res;
			};

			triangle.SetIntegrationLaw(3);
			result_I += triangle.SolveIntegral(I);
		}
		return result_I;
	};
	auto SolveS = [&solver_grid](std::vector<std::vector<int>>& boundary) -> double
	{
		double result_S = 0;
		for (int id_tr = 0; id_tr < boundary.size(); id_tr++)
		{
			geometry::Triangle triangle(
				solver_grid.GetCoordinateViaID(boundary[id_tr][1]),
				solver_grid.GetCoordinateViaID(boundary[id_tr][2]),
				solver_grid.GetCoordinateViaID(boundary[id_tr][3]));
			result_S += triangle.GetVolume();
		}
		return result_S;
	};

	S_in_sections[0] = SolveS(boundaries[target_types[0]]);
	I_in_sections[0] = SolveI(boundaries[target_types[0]]);
	S_in_sections[1] = SolveS(boundaries[target_types[1]]);
	I_in_sections[1] = SolveI(boundaries[target_types[1]]);


	double L; //расстояние между электродами
	switch (typeXYZ)
	{
	case 0:
		L = Max.x - Min.x;
		break;
	case 1:
		L = Max.y - Min.y;
		break;
	case 2:
		L = Max.z - Min.z;
		break;
	}

	double rpho_eff_0 = S_in_sections[0] * deltaU / (I_in_sections[0] * L);
	double rpho_eff_1 = S_in_sections[1] * deltaU / (I_in_sections[1] * L);

	return (rpho_eff_0 + rpho_eff_1) / 2.0;
}



void EffectiveSigmaZ_3Dmsh_ManyObj_WithBoundaryTag_1DOF(char* properties_file)
{
	math::SimpleGrid geo_grid; //input
	FEM::Grid_forScal solver_grid; //output
	std::vector<double> Solution; //output

	clock_t t_after = clock();
	double start = omp_get_wtime();

	printf_s("\n====================Scalar Electrical Problem======================\n");
	//char properties_file[1000] = { "D:/testing23062022/testing1408/new_FEM/m0/properties.txt" };
	//char properties_file[1000] = { "D:/testing23062022/testing1408/new_FEM/simple/properties.txt" };
	//char properties_file[1000] = { "properties.txt" };
	printf_s("The properties file contains the information:\n==========================================\n");
	printf_s("\tDo you need to print to a log file? 0-false/1-true\n");
	printf_s("\tinput: Mesh directory: box_200x200x200/90el\n");
	printf_s("\t<Critical residual value for SLAE solution>\n");
	printf_s("\t<Type effective properties (X/Y/Z)>\n");
	printf_s("\t<Type cross section for solve effective properties (circle/rectangle)>\n");
	printf_s("\t<properties: Boundary values>\n");
	printf_s("\t<properties: Materials>\n");
	printf_s("\t<output: Result directory: box_200x200x200/90el/result>\n==========================================\n");


	//printf_s("Enter the name of the properties file: ");
	//scanf_s("%s", &properties_file);

	FILE* f_properties;
	fopen_s(&f_properties, properties_file, "r");
	if (f_properties == NULL)
	{
		printf_s("\nError in properties file\n");
	}
	bool is_print_logFile = false;
	char mesh_directory[1000];
	char cracks_directory[1000];
	char base_result_directory[1000];


	{
		char _line[1000];
		math::ReadNonEmptyLine_forNumbers(f_properties, _line);
		std::vector<int> val_i;
		math::ParserStringToVectorInt(_line, val_i, " ");

		if (val_i[0] == 1) is_print_logFile = true;

		math::ReadNonEmptyLine_forNumbers(f_properties, _line);
		math::ParserStringToVectorInt(_line, val_i, " ");
		math::NUM_THREADS = val_i[0];
	}

	int num_in_FE;
	{
		char _line[1000];
		math::ReadNonEmptyLine_forNumbers(f_properties, _line);
		std::vector<int> val_i;
		math::ParserStringToVectorInt(_line, val_i, " ");

		num_in_FE = val_i[0];
	}

	math::ReadNonEmptyLine(f_properties, mesh_directory);
	std::vector<std::vector<std::vector<int>>> boundary_faces;
	switch (num_in_FE)
	{
	case 9:
		geo_grid.ReadFromMSH_v2208(mesh_directory, 1, 4, 4, boundary_faces);
		break;
	case 10:
		geo_grid.ReadFromMSH(mesh_directory, 1, 4, 4, boundary_faces);
		break;
	default:
		break;
	}

	double critical_residual;
	{
		char _line[1000];
		math::ReadNonEmptyLine_forNumbers(f_properties, _line);
		std::vector<float> val;
		math::ParserStringToVectorFloat(_line, val, " ");
		critical_residual = val[0];
	}

	int typeXYZ_effective = 0, typeGEO_section = 0;
	std::vector<double> SectionPosition;
	{
		char _line[1000];
		math::ReadNonEmptyLine(f_properties, _line);
		for (int i = 0; i < sizeof(_line); i++)
			if (_line[i] == ' ') { _line[i] = '\0'; break; }

		if (!strcmp(_line, "X"))
			typeXYZ_effective = 0;
		if (!strcmp(_line, "Y"))
			typeXYZ_effective = 1;
		if (!strcmp(_line, "Z"))
			typeXYZ_effective = 2;

		math::ReadNonEmptyLine(f_properties, _line);
		for (int i = 0; i < sizeof(_line); i++)
			if (_line[i] == ' ') { _line[i] = '\0'; break; }

		if (!strcmp(_line, "circle"))
			typeGEO_section = 0;
		if (!strcmp(_line, "rectangle"))
			typeGEO_section = 1;

		math::ReadNonEmptyLine_forNumbers(f_properties, _line);
		std::vector<float> val_f;
		math::ParserStringToVectorFloat(_line, val_f, " ");
		SectionPosition.resize(val_f.size());
		for (int i = 0; i < val_f.size(); i++)
			SectionPosition[i] = (double)val_f[i];
	}

	//read Boundary values and vertexes
	struct _Dirichlet {
		/*Point<double> value_const;
		Point<bool> is_condition;*/
		std::vector<int> id_vertexes;
		std::function<double(Point<double>)> value;
	};
	std::vector<int> id_boundaries;
	std::vector<_Dirichlet> first_boundaries;
	{
		int Nb;

		{
			char _line[1000];
			math::ReadNonEmptyLine_forNumbers(f_properties, _line);
			std::vector<int> val;
			math::ParserStringToVectorInt(_line, val, " ");
			Nb = val[0];
		}

		first_boundaries.resize(Nb);
		id_boundaries.resize(Nb);
		char line[1000];
		for (int i = 0; i < Nb; i++)
		{
			double value;
			bool is_condition;

			math::ReadNonEmptyLine_forNumbers(f_properties, line);

			if (line[0] != '#')
			{
				id_boundaries[i] = math::ParserCharToInt(line[0]);
				int jj = 0;
				int curr_i = 0;
				char _tmp_line[1000];
				for (int j = 1; j < 1000; j++)
				{
					switch (line[j])
					{
					case '\0': j = 1000; break;
					case '\n': j = 1000; break;
					case '#': j = 1000; break;
					case '*':
					{
						switch (jj)
						{
						case 0: is_condition = false; break;
						default:
							break;
						}
						jj++;
						break;
					}
					case ' ':
					{
						if (curr_i > 0)
						{
							std::vector<float> _val;
							_tmp_line[curr_i] = '\0';
							math::ParserStringToVectorFloat(_tmp_line, _val, " *");
							curr_i = 0;

							switch (jj)
							{
							case 0: is_condition = true; value = _val[0]; break;
							default:
								break;
							}
							jj++;
						}
						break;
					}
					case '\t': break;
					default:
					{
						_tmp_line[curr_i] = line[j];
						curr_i++;
					}
					break;
					}
				}
			}

			first_boundaries[i].value = [value](Point<double> x)->double
			{
				return value;
			};
		}

		if (Nb != 0)
		{
			for (int id_type = 0; id_type < Nb; id_type++)
			{
				std::vector<int> tmp_vert;
				for (int i = 0; i < boundary_faces[id_boundaries[id_type]].size(); i++)
				{
					for (int j = 1; j < boundary_faces[id_boundaries[id_type]][i].size(); j++)
					{
						tmp_vert.push_back(boundary_faces[id_boundaries[id_type]][i][j]);
					}
				}
				math::MakeQuickSort(tmp_vert);
				math::MakeRemovalOfDuplication(tmp_vert, first_boundaries[id_type].id_vertexes);
			}
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
	std::vector<_Neumann> second_boundary;
	{
		int Nb;
		{
			char _line[1000];
			math::ReadNonEmptyLine_forNumbers(f_properties, _line);
			std::vector<int> val;
			math::ParserStringToVectorInt(_line, val, " ");
			Nb = val[0];
		}
		second_boundary.resize(Nb);
		std::vector<int> id_boundaries(Nb);
		char line[1000];
		std::vector<bool> is_individual_values(Nb);
		std::vector<double> individual_values(Nb);
		for (int i = 0; i < Nb; i++)
		{
			is_individual_values[i] = false;
			double value;
			Point<double> vector;

			math::ReadNonEmptyLine_forNumbers(f_properties, line);
			id_boundaries[i] = math::ParserCharToInt(line[0]);
			int jj = 0;
			int curr_i = 0;
			char _tmp_line[1000];
			for (int j = 1; j < 1000; j++)
			{
				if ((curr_i > 0) && (line[j] == '\n' || line[j] == '\0' || line[j] == '#'))
				{
					std::vector<float> val;
					_tmp_line[curr_i] = '\0';
					math::ParserStringToVectorFloat(_tmp_line, val, " ");
					curr_i = 0;

					value = val[0];

					if (val.size() == 4)
					{
						vector.x = val[1];
						vector.y = val[2];
						vector.z = val[3];

						individual_values[i] = value;
						if (math::IsEqual(math::SolveLengthVector(vector), 0.0))
						{
							is_individual_values[i] = true;
						}
					}
					break;
				}
				else
				{
					if (line[j] != '\t')
					{
						_tmp_line[curr_i] = line[j];
						curr_i++;
					}
				}
			}

			if (is_individual_values[i] == false)
			{
				second_boundary[i].value = [value, vector](Point<double> X) -> Point<double>
				{
					Point <double> res;
					res.x = vector.x * value;
					res.y = vector.y * value;
					res.z = vector.z * value;
					return res;
				};
			}
		}

		if (Nb != 0)
		{
			for (int id_type = 0; id_type < Nb; id_type++)
			{
				second_boundary[id_type].id_vertexes_as_triangle.resize(boundary_faces[id_boundaries[id_type]].size());
				for (int id_triang = 0; id_triang < boundary_faces[id_boundaries[id_type]].size(); id_triang++)
				{
					printf_s("Create boundary condition %d (triangle %d/%d)\r", id_boundaries[id_type], id_triang, boundary_faces[id_boundaries[id_type]].size());
					int test_vertex = -1;
					int base_elem = boundary_faces[id_boundaries[id_type]][id_triang][0];

					for (int i = 1; i < boundary_faces[id_boundaries[id_type]][id_triang].size(); i++)
						second_boundary[id_type].id_vertexes_as_triangle[id_triang].push_back(boundary_faces[id_boundaries[id_type]][id_triang][i]);
					second_boundary[id_type].id_base_element.push_back(base_elem);

					//по собственным нормалям
					if (is_individual_values[id_type] == true)
					{
						double A, B, C, D;
						Point<double> test_vector = geo_grid.xyz[test_vertex];
						std::vector<Point<double>> vertexes(3);
						vertexes[0] = geo_grid.xyz[second_boundary[id_type].id_vertexes_as_triangle[id_triang][0]];
						vertexes[1] = geo_grid.xyz[second_boundary[id_type].id_vertexes_as_triangle[id_triang][1]];
						vertexes[2] = geo_grid.xyz[second_boundary[id_type].id_vertexes_as_triangle[id_triang][2]];
						bool reverse = false;
						math::GetPlaneEquation(vertexes, A, B, C, D);

						if (Point<double>(A, B, C) * test_vector > 0)
						{
							int t = second_boundary[id_type].id_vertexes_as_triangle[id_triang][1];
							second_boundary[id_type].id_vertexes_as_triangle[id_triang][1] = second_boundary[id_type].id_vertexes_as_triangle[id_triang][2];
							second_boundary[id_type].id_vertexes_as_triangle[id_triang][2] = t;

							vertexes[0] = geo_grid.xyz[second_boundary[id_type].id_vertexes_as_triangle[id_triang][0]];
							vertexes[1] = geo_grid.xyz[second_boundary[id_type].id_vertexes_as_triangle[id_triang][1]];
							vertexes[2] = geo_grid.xyz[second_boundary[id_type].id_vertexes_as_triangle[id_triang][2]];
						}
						Point<double> normal;
						double d;
						math::GetPlaneEquation(vertexes, normal.x, normal.y, normal.z, d);
						normal /= math::SolveLengthVector(normal);

						double _val = individual_values[id_type];
						std::function<Point<double>(Point<double>)> curr_val = [_val, normal](Point<double> X) -> Point<double>
						{
							Point <double> res;
							res.x = normal.x * _val;
							res.y = normal.y * _val;
							res.z = normal.z * _val;
							return res;
						};
						second_boundary[id_type].values.push_back(curr_val);
					}
				}
			}
		}
	}

	//read Materials
	struct _material
	{
		double _sigma;

		_material(double _sigma)
		{
			this->_sigma = _sigma;
		}
	};
	std::vector<_material> _materials;
	int N_domain;
	{
		char _line[1000];
		math::ReadNonEmptyLine_forNumbers(f_properties, _line);
		std::vector<int> val;
		math::ParserStringToVectorInt(_line, val, " ");
		N_domain = val[0];
	}
	for (int id_domain = 0; id_domain < N_domain; id_domain++)
	{
		double _sigma;
		{
			char _line[1000];
			math::ReadNonEmptyLine_forNumbers(f_properties, _line);
			std::vector<float> val;
			math::ParserStringToVectorFloat(_line, val, " ");
			_sigma = val[0];
		}
		_materials.push_back(_material(_sigma));
	}

	math::ReadNonEmptyLine(f_properties, base_result_directory);
	fclose(f_properties);

	wchar_t _tmp_wc[1000];
	math::Char_To_Wchar_t(base_result_directory, _tmp_wc, 500);
	CreateDirectory((LPCTSTR)_tmp_wc, NULL);

	{
		printf_s("\n=============================*===============================\n");

		for (int i = 0; i < _materials.size(); i++)
		{
			solver_grid.AddDomain();
			auto domain = solver_grid.GetDomain(i);
			domain->forElectrical.sigma = _materials[i]._sigma;
		}

		std::function<double(int, Point<double>)> StiffnessCoef = [&](int elem, Point<double> X)->double
		{
			//return (20E-7)*(20E-7);

			return solver_grid.GetDomain(solver_grid.GetElement(elem)->GetIdDomain())->forElectrical.sigma;
		};

		FEM::FEM_forBaseEliptic(
			is_print_logFile,
			critical_residual,
			geo_grid, //input
			first_boundaries, //input
			second_boundary, //input
			StiffnessCoef, //input
			base_result_directory, //output
			solver_grid, //output
			Solution //output
		);

		printf("Solve effective resistans... \n");
		double R = 0, R_new = 0, norm_Jsqrt = 0, rpho_eff_summ, GLOBAL_SIGMA_hmm;
		integrate_for_rpho_eff result_of_integration;
		result_of_integration.integr_in_domain_EE.resize(solver_grid.GetDomainsCount());
		result_of_integration.integr_in_domain_JJ.resize(solver_grid.GetDomainsCount());
		result_of_integration.V_in_domain.resize(solver_grid.GetDomainsCount());
		R_new = EffectiveResistance_new_LIGHT(solver_grid, Solution, norm_Jsqrt, rpho_eff_summ, GLOBAL_SIGMA_hmm, result_of_integration);

		double deltaU = abs(first_boundaries[0].value(Point<double>(0, 0, 0)) - first_boundaries[1].value(Point<double>(0, 0, 0)));
		std::vector<double> I_in_sections;
		double R_sections = EffectiveResistance_ByСrossSection(typeXYZ_effective, typeGEO_section, deltaU, SectionPosition, solver_grid, Solution, I_in_sections);

		std::vector<double> I_in_bases, S_in_bases;
		double R_bases = EffectiveResistance_ByBoundary(
			typeXYZ_effective,
			boundary_faces,
			id_boundaries,
			deltaU,
			solver_grid,
			Solution,
			I_in_bases, S_in_bases);
		printf("\tcomplit\n");

		clock_t t_before = clock();

		printf("TIME = %lf sec\n", (t_before - t_after) * 1.0 / CLK_TCK);
		double end = omp_get_wtime();
		printf("start = %.16g\nend = %.16g\ndiff = %.16g\n", start, end, end - start);

		FILE* fout;
		FILE* fout_new;
		char name_out[1000];
		sprintf_s(name_out, "%s/effective_sigma.txt", base_result_directory);
		fopen_s(&fout, name_out, "a");

		for (int i = 0; i < solver_grid.GetDomainsCount(); i++)
			fprintf_s(fout, "sigma%d = %.5e\trpho%d = %.5e\n",
				i, solver_grid.GetDomain(i)->forElectrical.sigma,
				i, 1 / solver_grid.GetDomain(i)->forElectrical.sigma);

		char type;
		switch (typeXYZ_effective)
		{
		case 0: type = 'X'; break;
		case 1: type = 'Y'; break;
		case 2: type = 'Z'; break;
		}
		fprintf_s(fout, "By volume: sigma(%c)_eff = %.15lf\trpho(%c)_eff = %.15lf\n", type, 1. / R_new, type, R_new);
		fprintf_s(fout, "By sections: sigma(%c)_eff = %.15lf\trpho(%c)_eff = %.15lf\n", type, 1. / R_sections, type, R_sections);
		fprintf_s(fout, "By bases: sigma(%c)_eff = %.15lf\trpho(%c)_eff = %.15lf\n", type, 1. / R_bases, type, R_bases);

		fprintf_s(fout, "\nIntegration results:\n");
		fprintf_s(fout, "\tCross sections:\n");
		fprintf_s(fout, "\t\tPositions:");
		for (int i = 0; i < SectionPosition.size(); i++)
		{
			fprintf_s(fout, "\t%.2lf", SectionPosition[i]);
		}
		fprintf_s(fout, "\n\t\tIntegral:");
		for (int i = 0; i < SectionPosition.size(); i++)
		{
			fprintf_s(fout, "\t%.4e", I_in_sections[i]);
		}
		fprintf_s(fout, "\n\n");

		fprintf_s(fout, "\tBases:\n");
		fprintf_s(fout, "\n\t\tArea:");
		for (int i = 0; i < I_in_bases.size(); i++)
		{
			fprintf_s(fout, "\t%.4e", S_in_bases[i]);
		}
		fprintf_s(fout, "\n\t\tIntegral:");
		for (int i = 0; i < I_in_bases.size(); i++)
		{
			fprintf_s(fout, "\t%.4e", I_in_bases[i]);
		}
		fprintf_s(fout, "\n\n");

		fprintf_s(fout, "\t(E.x*E.x, E.y*E.y, E.z*E.z) in all domain: %.2e %.2e %.2e\n", result_of_integration.integr_EE_full.x, result_of_integration.integr_EE_full.y, result_of_integration.integr_EE_full.z);
		fprintf_s(fout, "\t(J.x*J.x, J.y*J.y, J.z*J.z) in all domain: %.2e %.2e %.2e\n", result_of_integration.integr_JJ_full.x, result_of_integration.integr_JJ_full.y, result_of_integration.integr_JJ_full.z);
		fprintf_s(fout, "\tsqrt(E.x*E.x + E.y*E.y + E.z*E.z) in all domain: %.2e\n", sqrt(result_of_integration.integr_EE_full.x + result_of_integration.integr_EE_full.y + result_of_integration.integr_EE_full.z));
		fprintf_s(fout, "\tsqrt(J.x*J.x + J.y*J.y + J.z*J.z) in all domain: %.2e\n", sqrt(result_of_integration.integr_JJ_full.x + result_of_integration.integr_JJ_full.y + result_of_integration.integr_JJ_full.z));

		for (int i = 0; i < solver_grid.GetDomainsCount(); i++)
		{
			fprintf_s(fout, "\tIn %d domain:\n", i);
			fprintf_s(fout, "\tV = %.2e\n", result_of_integration.V_in_domain[i]);
			fprintf_s(fout, "\t\t(E.x*E.x, E.y*E.y, E.z*E.z): %.2e %.2e %.2e\n", result_of_integration.integr_in_domain_EE[i].x, result_of_integration.integr_in_domain_EE[i].y, result_of_integration.integr_in_domain_EE[i].z);
			fprintf_s(fout, "\t\t(J.x*J.x, J.y*J.y, J.z*J.z): %.2e %.2e %.2e\n", result_of_integration.integr_in_domain_JJ[i].x, result_of_integration.integr_in_domain_JJ[i].y, result_of_integration.integr_in_domain_JJ[i].z);
		}

		double V = 0;
		std::vector<double> Vin(solver_grid.GetDomainsCount());
		int qwe = 0;
		double conc = 0.0;

		for (int i = 0; i < solver_grid.GetElementsCount(); i++)
		{
			auto element = solver_grid.GetElement(i);
			Vin[element->GetIdDomain()] += element->GetVolume();
			V += element->GetVolume();
		}
		for (int i = 0; i < Vin.size(); i++)
		{
			printf("concentration(id%d) = %.5lf %%\n", i, (Vin[i] / V) * 100.);
			fprintf_s(fout, "concentration(id%d) = %.5lf %%\n", i, (Vin[i] / V) * 100.);
		}

		fprintf_s(fout, "TIME = %lf sec\n\n", (t_before - t_after) * 1.0 / CLK_TCK);
		fclose(fout);
		printf_s("\t complite\n");

		//output solution
		printf_s("Print the mech result into .dat file... ");
		{
			FILE* fout_tech;
			char name_u_tech[5000];
			sprintf_s(name_u_tech, "%s/U.dat", base_result_directory);
			fopen_s(&fout_tech, name_u_tech, "w");
			char name_in_file[1000];
			sprintf_s(name_in_file, "ScalarPotential");
			std::vector<std::vector<char>> name_value(5);
			char name_v_tmp[5][100];
			sprintf_s(name_v_tmp[0], "U");
			sprintf_s(name_v_tmp[1], "Jx");
			sprintf_s(name_v_tmp[2], "Jy");
			sprintf_s(name_v_tmp[3], "Jz");
			sprintf_s(name_v_tmp[4], "material");
			for (int i = 0; i < name_value.size(); i++)
			{
				name_value[i].resize(100);
				for (int j = 0; j < name_value[i].size(); j++)
				{
					name_value[i][j] = name_v_tmp[i][j];
				}
			}
			std::vector<std::vector<double>> value(5);
			value[0].resize(solver_grid.GetElementsCount());
			value[1].resize(solver_grid.GetElementsCount());
			value[2].resize(solver_grid.GetElementsCount());
			value[3].resize(solver_grid.GetElementsCount());
			value[4].resize(solver_grid.GetElementsCount());
			double sigma_inv_max = 0;
			int elem_sigma_max = 0;
			for (int i = 0; i < solver_grid.GetElementsCount(); i++)
			{
				auto element = solver_grid.GetElement(i);

				Point<double> Centr = element->GetWeightCentr();

				double U = solver_grid.GetSolutionInPoint(i, Centr, Solution);
				Point<double> dU = solver_grid.GetDerevativeFromSolutionInPoint(i, Centr, Solution);
				double sigma = solver_grid.GetDomain(element->GetIdDomain())->forElectrical.sigma;

				value[0][i] = U;
				value[1][i] = -sigma * dU.x;
				value[2][i] = -sigma * dU.y;
				value[3][i] = -sigma * dU.z;
				value[4][i] = element->GetIdDomain();
			}
			//solver_grid.printTecPlot3D(fout_tech, value, name_value, name_in_file);
			solver_grid.printTecPlot3D_DiffDomains(fout_tech, value, name_value, name_in_file);

			fclose(fout_tech);
		}
		printf_s("\t complite\n");



		solver_grid.~Grid_forScal();
	}
}
void EffectiveSigmaZ_3Dmsh_ManyObj_WithBoundaryTag_2DOF()
{
	math::SimpleGrid geo_grid; //input
	FEM::Grid_forScal_OrderBF2 solver_grid; //output
	std::vector<double> Solution; //output

	clock_t t_after = clock();
	double start = omp_get_wtime();

	printf_s("\n====================Scalar Electrical Problem======================\n");
	//char properties_file[1000] = { "D:/testing23062022/testing1408/new_FEM/m0/properties.txt" };
	//char properties_file[1000] = { "D:/testing23062022/testing1408/new_FEM/simple/properties.txt" };
	char properties_file[1000] = { "properties.txt" };
	printf_s("The properties file contains the information:\n==========================================\n");
	printf_s("\tDo you need to print to a log file? 0-false/1-true\n");
	printf_s("\tinput: Mesh directory: box_200x200x200/90el\n");
	printf_s("\t<Critical residual value for SLAE solution>\n");
	printf_s("\t<Type effective properties (X/Y/Z)>\n");
	printf_s("\t<Type cross section for solve effective properties (circle/rectangle)>\n");
	printf_s("\t<properties: Boundary values>\n");
	printf_s("\t<properties: Materials>\n");
	printf_s("\t<output: Result directory: box_200x200x200/90el/result>\n==========================================\n");


	printf_s("Enter the name of the properties file: ");
	//scanf_s("%s", &properties_file);

	FILE* f_properties;
	fopen_s(&f_properties, properties_file, "r");
	if (f_properties == NULL)
	{
		printf_s("\nError in properties file\n");
	}
	bool is_print_logFile = false;
	char mesh_directory[1000];
	char cracks_directory[1000];
	char base_result_directory[1000];


	int _flag;
	fscanf_s(f_properties, "%d", &_flag);
	if (_flag == 1) is_print_logFile = true;

	math::ReadNonEmptyLine(f_properties, mesh_directory);
	std::vector<std::vector<std::vector<int>>> boundary_faces;
	geo_grid.ReadFromMSH_v2208(mesh_directory, 1, 4, 4, boundary_faces);

	double critical_residual;
	{
		char _line[1000];
		math::ReadNonEmptyLine_forNumbers(f_properties, _line);
		std::vector<float> val;
		math::ParserStringToVectorFloat(_line, val, " ");
		critical_residual = val[0];
	}

	int typeXYZ_effective = 0, typeGEO_section = 0;
	std::vector<double> SectionPosition;
	{
		char _line[1000];
		math::ReadNonEmptyLine(f_properties, _line);
		for (int i = 0; i < sizeof(_line); i++)
			if (_line[i] == ' ') { _line[i] = '\0'; break; }

		if (!strcmp(_line, "X"))
			typeXYZ_effective = 0;
		if (!strcmp(_line, "Y"))
			typeXYZ_effective = 1;
		if (!strcmp(_line, "Z"))
			typeXYZ_effective = 2;

		math::ReadNonEmptyLine(f_properties, _line);
		for (int i = 0; i < sizeof(_line); i++)
			if (_line[i] == ' ') { _line[i] = '\0'; break; }

		if (!strcmp(_line, "circle"))
			typeGEO_section = 0;
		if (!strcmp(_line, "rectangle"))
			typeGEO_section = 1;

		math::ReadNonEmptyLine_forNumbers(f_properties, _line);
		std::vector<float> val_f;
		math::ParserStringToVectorFloat(_line, val_f, " ");
		SectionPosition.resize(val_f.size());
		for (int i = 0; i < val_f.size(); i++)
			SectionPosition[i] = (double)val_f[i];
	}

	//read Boundary values and vertexes
	struct _Dirichlet {
		/*Point<double> value_const;
		Point<bool> is_condition;*/
		std::vector<int> id_vertexes;
		std::function<double(Point<double>)> value;
	};
	std::vector<_Dirichlet> first_boundaries;
	std::vector<int> id_boundaries;
	{
		int Nb;
		{
			char _line[1000];
			math::ReadNonEmptyLine_forNumbers(f_properties, _line);
			std::vector<int> val;
			math::ParserStringToVectorInt(_line, val, " ");
			Nb = val[0];
		}

		first_boundaries.resize(Nb);
		id_boundaries.resize(Nb);
		char line[1000];
		for (int i = 0; i < Nb; i++)
		{
			double value;
			bool is_condition;

			math::ReadNonEmptyLine_forNumbers(f_properties, line);

			if (line[0] != '#')
			{
				id_boundaries[i] = math::ParserCharToInt(line[0]);
				int jj = 0;
				int curr_i = 0;
				char _tmp_line[1000];
				for (int j = 1; j < 1000; j++)
				{
					switch (line[j])
					{
					case '\0': j = 1000; break;
					case '\n': j = 1000; break;
					case '#': j = 1000; break;
					case '*':
					{
						switch (jj)
						{
						case 0: is_condition = false; break;
						default:
							break;
						}
						jj++;
						break;
					}
					case ' ':
					{
						if (curr_i > 0)
						{
							std::vector<float> _val;
							_tmp_line[curr_i] = '\0';
							math::ParserStringToVectorFloat(_tmp_line, _val, " *");
							curr_i = 0;

							switch (jj)
							{
							case 0: is_condition = true; value = _val[0]; break;
							default:
								break;
							}
							jj++;
						}
						break;
					}
					case '\t': break;
					default:
					{
						_tmp_line[curr_i] = line[j];
						curr_i++;
					}
					break;
					}
				}
			}

			first_boundaries[i].value = [value](Point<double> x)->double
			{
				return value;
			};
		}

		if (Nb != 0)
		{
			for (int id_type = 0; id_type < Nb; id_type++)
			{
				std::vector<int> tmp_vert;
				for (int i = 0; i < boundary_faces[id_boundaries[id_type]].size(); i++)
				{
					for (int j = 1; j < boundary_faces[id_boundaries[id_type]][i].size(); j++)
					{
						tmp_vert.push_back(boundary_faces[id_boundaries[id_type]][i][j]);
					}
				}
				math::MakeQuickSort(tmp_vert);
				math::MakeRemovalOfDuplication(tmp_vert, first_boundaries[id_type].id_vertexes);
			}
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
	std::vector<_Neumann> second_boundary;
	{
		int Nb;
		{
			char _line[1000];
			math::ReadNonEmptyLine_forNumbers(f_properties, _line);
			std::vector<int> val;
			math::ParserStringToVectorInt(_line, val, " ");
			Nb = val[0];
		}
		second_boundary.resize(Nb);
		std::vector<int> id_boundaries(Nb);
		char line[1000];
		std::vector<bool> is_individual_values(Nb);
		std::vector<double> individual_values(Nb);
		for (int i = 0; i < Nb; i++)
		{
			is_individual_values[i] = false;
			double value;
			Point<double> vector;

			math::ReadNonEmptyLine_forNumbers(f_properties, line);
			id_boundaries[i] = math::ParserCharToInt(line[0]);
			int jj = 0;
			int curr_i = 0;
			char _tmp_line[1000];
			for (int j = 1; j < 1000; j++)
			{
				if ((curr_i > 0) && (line[j] == '\n' || line[j] == '\0' || line[j] == '#'))
				{
					std::vector<float> val;
					_tmp_line[curr_i] = '\0';
					math::ParserStringToVectorFloat(_tmp_line, val, " ");
					curr_i = 0;

					value = val[0];
					vector.x = val[1];
					vector.y = val[2];
					vector.z = val[3];

					individual_values[i] = value;
					if (math::IsEqual(math::SolveLengthVector(vector), 0.0))
					{
						is_individual_values[i] = true;
					}

					break;
				}
				else
				{
					if (line[j] != '\t')
					{
						_tmp_line[curr_i] = line[j];
						curr_i++;
					}
				}
			}

			if (is_individual_values[i] == false)
			{
				second_boundary[i].value = [value, vector](Point<double> X) -> Point<double>
				{
					Point <double> res;
					res.x = vector.x * value;
					res.y = vector.y * value;
					res.z = vector.z * value;
					return res;
				};
			}
		}

		if (Nb != 0)
		{
			for (int id_type = 0; id_type < Nb; id_type++)
			{
				second_boundary[id_type].id_vertexes_as_triangle.resize(boundary_faces[id_boundaries[id_type]].size());
				for (int id_triang = 0; id_triang < boundary_faces[id_boundaries[id_type]].size(); id_triang++)
				{
					printf_s("Create boundary condition %d (triangle %d/%d)\r", id_boundaries[id_type], id_triang, boundary_faces[id_boundaries[id_type]].size());
					int test_vertex = -1;
					int base_elem = boundary_faces[id_boundaries[id_type]][id_triang][0];

					for (int i = 1; i < boundary_faces[id_boundaries[id_type]][id_triang].size(); i++)
						second_boundary[id_type].id_vertexes_as_triangle[id_triang].push_back(boundary_faces[id_boundaries[id_type]][id_triang][i]);
					second_boundary[id_type].id_base_element.push_back(base_elem);

					//по собственным нормалям
					if (is_individual_values[id_type] == true)
					{
						double A, B, C, D;
						Point<double> test_vector = geo_grid.xyz[test_vertex];
						std::vector<Point<double>> vertexes(3);
						vertexes[0] = geo_grid.xyz[second_boundary[id_type].id_vertexes_as_triangle[id_triang][0]];
						vertexes[1] = geo_grid.xyz[second_boundary[id_type].id_vertexes_as_triangle[id_triang][1]];
						vertexes[2] = geo_grid.xyz[second_boundary[id_type].id_vertexes_as_triangle[id_triang][2]];
						bool reverse = false;
						math::GetPlaneEquation(vertexes, A, B, C, D);

						if (Point<double>(A, B, C) * test_vector > 0)
						{
							int t = second_boundary[id_type].id_vertexes_as_triangle[id_triang][1];
							second_boundary[id_type].id_vertexes_as_triangle[id_triang][1] = second_boundary[id_type].id_vertexes_as_triangle[id_triang][2];
							second_boundary[id_type].id_vertexes_as_triangle[id_triang][2] = t;

							vertexes[0] = geo_grid.xyz[second_boundary[id_type].id_vertexes_as_triangle[id_triang][0]];
							vertexes[1] = geo_grid.xyz[second_boundary[id_type].id_vertexes_as_triangle[id_triang][1]];
							vertexes[2] = geo_grid.xyz[second_boundary[id_type].id_vertexes_as_triangle[id_triang][2]];
						}
						Point<double> normal;
						double d;
						math::GetPlaneEquation(vertexes, normal.x, normal.y, normal.z, d);
						normal /= math::SolveLengthVector(normal);

						double _val = individual_values[id_type];
						std::function<Point<double>(Point<double>)> curr_val = [_val, normal](Point<double> X) -> Point<double>
						{
							Point <double> res;
							res.x = normal.x * _val;
							res.y = normal.y * _val;
							res.z = normal.z * _val;
							return res;
						};
						second_boundary[id_type].values.push_back(curr_val);
					}
				}
			}
		}
	}

	//read Materials
	struct _material
	{
		double _sigma;

		_material(double _sigma)
		{
			this->_sigma = _sigma;
		}
	};
	std::vector<_material> _materials;
	int N_domain;
	{
		char _line[1000];
		math::ReadNonEmptyLine_forNumbers(f_properties, _line);
		std::vector<int> val;
		math::ParserStringToVectorInt(_line, val, " ");
		N_domain = val[0];
	}
	for (int id_domain = 0; id_domain < N_domain; id_domain++)
	{
		double _sigma;
		{
			char _line[1000];
			math::ReadNonEmptyLine_forNumbers(f_properties, _line);
			std::vector<float> val;
			math::ParserStringToVectorFloat(_line, val, " ");
			_sigma = val[0];
		}
		_materials.push_back(_material(_sigma));
	}

	math::ReadNonEmptyLine(f_properties, base_result_directory);
	fclose(f_properties);

	wchar_t _tmp_wc[1000];
	math::Char_To_Wchar_t(base_result_directory, _tmp_wc, 500);
	CreateDirectory((LPCTSTR)_tmp_wc, NULL);

	{
		printf_s("\n=============================*===============================\n");

		for (int i = 0; i < _materials.size(); i++)
		{
			solver_grid.AddDomain();
			auto domain = solver_grid.GetDomain(i);
			domain->forElectrical.sigma = _materials[i]._sigma;
		}

		std::function<double(int, Point<double>)> StiffnessCoef = [&](int elem, Point<double> X)->double
		{
			//return (20E-7)*(20E-7);

			return solver_grid.GetDomain(solver_grid.GetElement(elem)->GetIdDomain())->forElectrical.sigma;
		};

		FEM::FEM_forBaseEliptic_OrderBF2(
			is_print_logFile,
			critical_residual,
			geo_grid, //input
			first_boundaries, //input
			second_boundary, //input
			StiffnessCoef, //input
			base_result_directory, //output
			solver_grid, //output
			Solution //output
		);

		printf("Solve effective resistans... \n");
		double R = 0, R_new = 0, norm_Jsqrt = 0, rpho_eff_summ, GLOBAL_SIGMA_hmm;
		integrate_for_rpho_eff result_of_integration;
		result_of_integration.integr_in_domain_EE.resize(solver_grid.GetDomainsCount());
		result_of_integration.integr_in_domain_JJ.resize(solver_grid.GetDomainsCount());
		result_of_integration.V_in_domain.resize(solver_grid.GetDomainsCount());
		R_new = EffectiveResistance_new_LIGHT(solver_grid, Solution, norm_Jsqrt, rpho_eff_summ, GLOBAL_SIGMA_hmm, result_of_integration);

		double deltaU = abs(first_boundaries[0].value(Point<double>(0, 0, 0)) - first_boundaries[1].value(Point<double>(0, 0, 0)));
		std::vector<double> I_in_sections;
		double R_sections = EffectiveResistance_ByСrossSection(typeXYZ_effective, typeGEO_section, deltaU, SectionPosition, solver_grid, Solution, I_in_sections);
		
		std::vector<double> I_in_bases, S_in_bases;
		double R_bases = EffectiveResistance_ByBoundary(
			typeXYZ_effective,
			boundary_faces,
			id_boundaries,
			deltaU,
			solver_grid,
			Solution,
			I_in_bases, S_in_bases);
		printf("\tcomplit\n");

		clock_t t_before = clock();

		printf("TIME = %lf sec\n", (t_before - t_after) * 1.0 / CLK_TCK);
		double end = omp_get_wtime();
		printf("start = %.16g\nend = %.16g\ndiff = %.16g\n", start, end, end - start);

		FILE* fout;
		FILE* fout_new;
		char name_out[1000];
		sprintf_s(name_out, "%s/effective_sigma.txt", base_result_directory);
		fopen_s(&fout, name_out, "a");

		for (int i = 0; i < solver_grid.GetDomainsCount(); i++)
			fprintf_s(fout, "sigma%d = %.5e\trpho%d = %.5e\n",
				i, solver_grid.GetDomain(i)->forElectrical.sigma,
				i, 1 / solver_grid.GetDomain(i)->forElectrical.sigma);

		char type;
		switch (typeXYZ_effective)
		{
		case 0: type = 'X'; break;
		case 1: type = 'Y'; break;
		case 2: type = 'Z'; break;
		}
		fprintf_s(fout, "By volume: sigma(%c)_eff = %.15lf\trpho(%c)_eff = %.15lf\n", type, 1. / R_new, type, R_new);
		fprintf_s(fout, "By sections: sigma(%c)_eff = %.15lf\trpho(%c)_eff = %.15lf\n", type, 1. / R_sections, type, R_sections);
		fprintf_s(fout, "By bases: sigma(%c)_eff = %.15lf\trpho(%c)_eff = %.15lf\n", type, 1. / R_bases, type, R_bases);

		fprintf_s(fout, "\nIntegration results:\n");
		fprintf_s(fout, "\tCross sections:\n");
		fprintf_s(fout, "\t\tPositions:");
		for (int i = 0; i < SectionPosition.size(); i++)
		{
			fprintf_s(fout, "\t%.2lf", SectionPosition[i]);
		}
		fprintf_s(fout, "\n\t\tIntegral:");
		for (int i = 0; i < SectionPosition.size(); i++)
		{
			fprintf_s(fout, "\t%.2e", I_in_sections[i]);
		}
		fprintf_s(fout, "\n\n");

		fprintf_s(fout, "\tBases:\n");
		fprintf_s(fout, "\n\t\tArea:");
		for (int i = 0; i < I_in_bases.size(); i++)
		{
			fprintf_s(fout, "\t%.2e", S_in_bases[i]);
		}
		fprintf_s(fout, "\n\t\tIntegral:");
		for (int i = 0; i < I_in_bases.size(); i++)
		{
			fprintf_s(fout, "\t%.2e", I_in_bases[i]);
		}
		fprintf_s(fout, "\n\n");
		fprintf_s(fout, "\t(E.x*E.x, E.y*E.y, E.z*E.z) in all domain: %.2e %.2e %.2e\n", result_of_integration.integr_EE_full.x, result_of_integration.integr_EE_full.y, result_of_integration.integr_EE_full.z);
		fprintf_s(fout, "\t(J.x*J.x, J.y*J.y, J.z*J.z) in all domain: %.2e %.2e %.2e\n", result_of_integration.integr_JJ_full.x, result_of_integration.integr_JJ_full.y, result_of_integration.integr_JJ_full.z);
		fprintf_s(fout, "\tsqrt(E.x*E.x + E.y*E.y + E.z*E.z) in all domain: %.2e\n", sqrt(result_of_integration.integr_EE_full.x + result_of_integration.integr_EE_full.y + result_of_integration.integr_EE_full.z));
		fprintf_s(fout, "\tsqrt(J.x*J.x + J.y*J.y + J.z*J.z) in all domain: %.2e\n", sqrt(result_of_integration.integr_JJ_full.x + result_of_integration.integr_JJ_full.y + result_of_integration.integr_JJ_full.z));

		for (int i = 0; i < solver_grid.GetDomainsCount(); i++)
		{
			fprintf_s(fout, "\tIn %d domain:\n", i);
			fprintf_s(fout, "\tV = %.2e\n", result_of_integration.V_in_domain[i]);
			fprintf_s(fout, "\t\t(E.x*E.x, E.y*E.y, E.z*E.z): %.2e %.2e %.2e\n", result_of_integration.integr_in_domain_EE[i].x, result_of_integration.integr_in_domain_EE[i].y, result_of_integration.integr_in_domain_EE[i].z);
			fprintf_s(fout, "\t\t(J.x*J.x, J.y*J.y, J.z*J.z): %.2e %.2e %.2e\n", result_of_integration.integr_in_domain_JJ[i].x, result_of_integration.integr_in_domain_JJ[i].y, result_of_integration.integr_in_domain_JJ[i].z);
		}

		double V = 0;
		std::vector<double> Vin(solver_grid.GetDomainsCount());
		int qwe = 0;
		double conc = 0.0;

		for (int i = 0; i < solver_grid.GetElementsCount(); i++)
		{
			auto element = solver_grid.GetElement(i);
			Vin[element->GetIdDomain()] += element->GetVolume();
			V += element->GetVolume();
		}
		for (int i = 0; i < Vin.size(); i++)
		{
			printf("concentration(id%d) = %.5lf %%\n", i, (Vin[i] / V) * 100.);
			fprintf_s(fout, "concentration(id%d) = %.5lf %%\n", i, (Vin[i] / V) * 100.);
		}

		fprintf_s(fout, "TIME = %lf sec\n\n", (t_before - t_after) * 1.0 / CLK_TCK);
		fclose(fout);
		printf_s("\t complite\n");

		//output solution
		printf_s("Print the mech result into .dat file... ");
		{
			FILE* fout_tech;
			char name_u_tech[5000];
			sprintf_s(name_u_tech, "%s/U.dat", base_result_directory);
			fopen_s(&fout_tech, name_u_tech, "w");
			char name_in_file[1000];
			sprintf_s(name_in_file, "ScalarPotential");
			std::vector<std::vector<char>> name_value(5);
			char name_v_tmp[5][100];
			sprintf_s(name_v_tmp[0], "U");
			sprintf_s(name_v_tmp[1], "Jx");
			sprintf_s(name_v_tmp[2], "Jy");
			sprintf_s(name_v_tmp[3], "Jz");
			sprintf_s(name_v_tmp[4], "material");
			for (int i = 0; i < name_value.size(); i++)
			{
				name_value[i].resize(100);
				for (int j = 0; j < name_value[i].size(); j++)
				{
					name_value[i][j] = name_v_tmp[i][j];
				}
			}
			std::vector<std::vector<double>> value(5);
			value[0].resize(solver_grid.GetElementsCount());
			value[1].resize(solver_grid.GetElementsCount());
			value[2].resize(solver_grid.GetElementsCount());
			value[3].resize(solver_grid.GetElementsCount());
			value[4].resize(solver_grid.GetElementsCount());
			double sigma_inv_max = 0;
			int elem_sigma_max = 0;
			for (int i = 0; i < solver_grid.GetElementsCount(); i++)
			{
				auto element = solver_grid.GetElement(i);

				Point<double> Centr = element->GetWeightCentr();

				double U = solver_grid.GetSolutionInPoint(i, Centr, Solution);
				Point<double> dU = solver_grid.GetDerevativeFromSolutionInPoint(i, Centr, Solution);
				double sigma = solver_grid.GetDomain(element->GetIdDomain())->forElectrical.sigma;

				value[0][i] = U;
				value[1][i] = -sigma * dU.x;
				value[2][i] = -sigma * dU.y;
				value[3][i] = -sigma * dU.z;
				value[4][i] = element->GetIdDomain();
			}
			//solver_grid.printTecPlot3D(fout_tech, value, name_value, name_in_file);
			solver_grid.printTecPlot3D_DiffDomains(fout_tech, value, name_value, name_in_file);

			fclose(fout_tech);
		}
		printf_s("\t complite\n");



		solver_grid.~Grid_forScal_OrderBF2();
	}
}
void EffectiveSigmaZ_3Dmsh_ManyObj_WithBoundaryTag_1DOF()
{
	char properties_file[1000] = { "properties.txt" };
	EffectiveSigmaZ_3Dmsh_ManyObj_WithBoundaryTag_1DOF(properties_file);
}

void ElectroStaticTensor_ManyObj_1DOF(char* properties_file)
{
	math::SimpleGrid geo_grid; //input
	FEM::Grid_forScal solver_grid; //output
	std::vector<double> Solution; //output

	clock_t t_after = clock();
	double start = omp_get_wtime();

	printf_s("\n====================Scalar Electrical Problem======================\n");
	//char properties_file[1000] = { "D:/testing23062022/testing1408/new_FEM/m0/properties.txt" };
	//char properties_file[1000] = { "D:/testing23062022/testing1408/new_FEM/simple/properties.txt" };
	//char properties_file[1000] = { "properties.txt" };
	printf_s("The properties file contains the information:\n==========================================\n");
	printf_s("\tDo you need to print to a log file? 0-false/1-true\n");
	printf_s("\tinput: Mesh directory: box_200x200x200/90el\n");
	printf_s("\t<Critical residual value for SLAE solution>\n");
	printf_s("\t<Type effective properties (X/Y/Z)>\n");
	printf_s("\t<Type cross section for solve effective properties (circle/rectangle)>\n");
	printf_s("\t<properties: Boundary values>\n");
	printf_s("\t<properties: Materials>\n");
	printf_s("\t<output: Result directory: box_200x200x200/90el/result>\n==========================================\n");


	printf_s("Enter the name of the properties file: ");
	//scanf_s("%s", &properties_file);

	FILE* f_properties;
	fopen_s(&f_properties, properties_file, "r");
	if (f_properties == NULL)
	{
		printf_s("\nError in properties file\n");
	}
	bool is_print_logFile = false;
	char mesh_directory[1000];
	char cracks_directory[1000];
	char base_result_directory[1000];


	int _flag;
	fscanf_s(f_properties, "%d", &_flag);
	if (_flag == 1) is_print_logFile = true;

	math::ReadNonEmptyLine(f_properties, mesh_directory);
	std::vector<std::vector<std::vector<int>>> boundary_faces;
	geo_grid.ReadFromMSH_v2208(mesh_directory, 1, 4, 4, boundary_faces);

	double critical_residual;
	{
		char _line[1000];
		math::ReadNonEmptyLine_forNumbers(f_properties, _line);
		std::vector<float> val;
		math::ParserStringToVectorFloat(_line, val, " ");
		critical_residual = val[0];
	}

	int typeXYZ_effective = 0, typeGEO_section = 0;
	std::vector<double> SectionPosition;
	{
		char _line[1000];
		math::ReadNonEmptyLine(f_properties, _line);
		for (int i = 0; i < sizeof(_line); i++)
			if (_line[i] == ' ') { _line[i] = '\0'; break; }

		if (!strcmp(_line, "X"))
			typeXYZ_effective = 0;
		if (!strcmp(_line, "Y"))
			typeXYZ_effective = 1;
		if (!strcmp(_line, "Z"))
			typeXYZ_effective = 2;

		math::ReadNonEmptyLine(f_properties, _line);
		for (int i = 0; i < sizeof(_line); i++)
			if (_line[i] == ' ') { _line[i] = '\0'; break; }

		if (!strcmp(_line, "circle"))
			typeGEO_section = 0;
		if (!strcmp(_line, "rectangle"))
			typeGEO_section = 1;

		math::ReadNonEmptyLine_forNumbers(f_properties, _line);
		std::vector<float> val_f;
		math::ParserStringToVectorFloat(_line, val_f, " ");
		SectionPosition.resize(val_f.size());
		for (int i = 0; i < val_f.size(); i++)
			SectionPosition[i] = (double)val_f[i];
	}

	//read Boundary values and vertexes
	struct _Dirichlet {
		/*Point<double> value_const;
		Point<bool> is_condition;*/
		std::vector<int> id_vertexes;
		std::function<double(Point<double>)> value;
	};
	std::vector<int> id_boundaries;
	std::vector<_Dirichlet> first_boundaries;
	{
		int Nb;

		{
			char _line[1000];
			math::ReadNonEmptyLine_forNumbers(f_properties, _line);
			std::vector<int> val;
			math::ParserStringToVectorInt(_line, val, " ");
			Nb = val[0];
		}

		first_boundaries.resize(Nb);
		id_boundaries.resize(Nb);
		char line[1000];
		for (int i = 0; i < Nb; i++)
		{
			double value;
			bool is_condition;

			math::ReadNonEmptyLine_forNumbers(f_properties, line);

			if (line[0] != '#')
			{
				id_boundaries[i] = math::ParserCharToInt(line[0]);
				int jj = 0;
				int curr_i = 0;
				char _tmp_line[1000];
				for (int j = 1; j < 1000; j++)
				{
					switch (line[j])
					{
					case '\0': j = 1000; break;
					case '\n': j = 1000; break;
					case '#': j = 1000; break;
					case '*':
					{
						switch (jj)
						{
						case 0: is_condition = false; break;
						default:
							break;
						}
						jj++;
						break;
					}
					case ' ':
					{
						if (curr_i > 0)
						{
							std::vector<float> _val;
							_tmp_line[curr_i] = '\0';
							math::ParserStringToVectorFloat(_tmp_line, _val, " *");
							curr_i = 0;

							switch (jj)
							{
							case 0: is_condition = true; value = _val[0]; break;
							default:
								break;
							}
							jj++;
						}
						break;
					}
					case '\t': break;
					default:
					{
						_tmp_line[curr_i] = line[j];
						curr_i++;
					}
					break;
					}
				}
			}

			first_boundaries[i].value = [value](Point<double> x)->double
			{
				return value;
			};
		}

		if (Nb != 0)
		{
			for (int id_type = 0; id_type < Nb; id_type++)
			{
				std::vector<int> tmp_vert;
				for (int i = 0; i < boundary_faces[id_boundaries[id_type]].size(); i++)
				{
					for (int j = 1; j < boundary_faces[id_boundaries[id_type]][i].size(); j++)
					{
						tmp_vert.push_back(boundary_faces[id_boundaries[id_type]][i][j]);
					}
				}
				math::MakeQuickSort(tmp_vert);
				math::MakeRemovalOfDuplication(tmp_vert, first_boundaries[id_type].id_vertexes);
			}
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
	std::vector<_Neumann> second_boundary;
	{
		int Nb;
		{
			char _line[1000];
			math::ReadNonEmptyLine_forNumbers(f_properties, _line);
			std::vector<int> val;
			math::ParserStringToVectorInt(_line, val, " ");
			Nb = val[0];
		}
		second_boundary.resize(Nb);
		std::vector<int> id_boundaries(Nb);
		char line[1000];
		std::vector<bool> is_individual_values(Nb);
		std::vector<double> individual_values(Nb);
		for (int i = 0; i < Nb; i++)
		{
			is_individual_values[i] = false;
			double value;
			Point<double> vector;

			math::ReadNonEmptyLine_forNumbers(f_properties, line);
			id_boundaries[i] = math::ParserCharToInt(line[0]);
			int jj = 0;
			int curr_i = 0;
			char _tmp_line[1000];
			for (int j = 1; j < 1000; j++)
			{
				if ((curr_i > 0) && (line[j] == '\n' || line[j] == '\0' || line[j] == '#'))
				{
					std::vector<float> val;
					_tmp_line[curr_i] = '\0';
					math::ParserStringToVectorFloat(_tmp_line, val, " ");
					curr_i = 0;

					value = val[0];
					vector.x = val[1];
					vector.y = val[2];
					vector.z = val[3];

					individual_values[i] = value;
					if (math::IsEqual(math::SolveLengthVector(vector), 0.0))
					{
						is_individual_values[i] = true;
					}

					break;
				}
				else
				{
					if (line[j] != '\t')
					{
						_tmp_line[curr_i] = line[j];
						curr_i++;
					}
				}
			}

			if (is_individual_values[i] == false)
			{
				second_boundary[i].value = [value, vector](Point<double> X) -> Point<double>
				{
					Point <double> res;
					res.x = vector.x * value;
					res.y = vector.y * value;
					res.z = vector.z * value;
					return res;
				};
			}
		}

		if (Nb != 0)
		{
			for (int id_type = 0; id_type < Nb; id_type++)
			{
				second_boundary[id_type].id_vertexes_as_triangle.resize(boundary_faces[id_boundaries[id_type]].size());
				for (int id_triang = 0; id_triang < boundary_faces[id_boundaries[id_type]].size(); id_triang++)
				{
					printf_s("Create boundary condition %d (triangle %d/%d)\r", id_boundaries[id_type], id_triang, boundary_faces[id_boundaries[id_type]].size());
					int test_vertex = -1;
					int base_elem = boundary_faces[id_boundaries[id_type]][id_triang][0];

					for (int i = 1; i < boundary_faces[id_boundaries[id_type]][id_triang].size(); i++)
						second_boundary[id_type].id_vertexes_as_triangle[id_triang].push_back(boundary_faces[id_boundaries[id_type]][id_triang][i]);
					second_boundary[id_type].id_base_element.push_back(base_elem);

					//по собственным нормалям
					if (is_individual_values[id_type] == true)
					{
						double A, B, C, D;
						Point<double> test_vector = geo_grid.xyz[test_vertex];
						std::vector<Point<double>> vertexes(3);
						vertexes[0] = geo_grid.xyz[second_boundary[id_type].id_vertexes_as_triangle[id_triang][0]];
						vertexes[1] = geo_grid.xyz[second_boundary[id_type].id_vertexes_as_triangle[id_triang][1]];
						vertexes[2] = geo_grid.xyz[second_boundary[id_type].id_vertexes_as_triangle[id_triang][2]];
						bool reverse = false;
						math::GetPlaneEquation(vertexes, A, B, C, D);

						if (Point<double>(A, B, C) * test_vector > 0)
						{
							int t = second_boundary[id_type].id_vertexes_as_triangle[id_triang][1];
							second_boundary[id_type].id_vertexes_as_triangle[id_triang][1] = second_boundary[id_type].id_vertexes_as_triangle[id_triang][2];
							second_boundary[id_type].id_vertexes_as_triangle[id_triang][2] = t;

							vertexes[0] = geo_grid.xyz[second_boundary[id_type].id_vertexes_as_triangle[id_triang][0]];
							vertexes[1] = geo_grid.xyz[second_boundary[id_type].id_vertexes_as_triangle[id_triang][1]];
							vertexes[2] = geo_grid.xyz[second_boundary[id_type].id_vertexes_as_triangle[id_triang][2]];
						}
						Point<double> normal;
						double d;
						math::GetPlaneEquation(vertexes, normal.x, normal.y, normal.z, d);
						normal /= math::SolveLengthVector(normal);

						double _val = individual_values[id_type];
						std::function<Point<double>(Point<double>)> curr_val = [_val, normal](Point<double> X) -> Point<double>
						{
							Point <double> res;
							res.x = normal.x * _val;
							res.y = normal.y * _val;
							res.z = normal.z * _val;
							return res;
						};
						second_boundary[id_type].values.push_back(curr_val);
					}
				}
			}
		}
	}

	//read Materials
	struct _material
	{
		Tensor2Rank3D _sigma_tensor;
		//double _sigma;

		_material(std::vector<double> _sigma)
		{
			switch (_sigma.size())
			{
			case 1: this->_sigma_tensor.SetValueDiag(_sigma[0]); break;
			case 6: this->_sigma_tensor.SetValue_from_xx_yy_zz_xy_xz_yz(_sigma); break;
			default:
				break;
			}
		}
	};
	std::vector<_material> _materials;
	int N_domain;
	{
		char _line[1000];
		math::ReadNonEmptyLine_forNumbers(f_properties, _line);
		std::vector<int> val;
		math::ParserStringToVectorInt(_line, val, " ");
		N_domain = val[0];
	}
	for (int id_domain = 0; id_domain < N_domain; id_domain++)
	{
		std::vector<double> _sigma;
		{
			char _line[1000];
			math::ReadNonEmptyLine_forNumbers(f_properties, _line);
			std::vector<double> val;
			math::ParserStringToVectorDouble(_line, val, " ");
			_sigma = val;
		}
		_materials.push_back(_material(_sigma));
	}

	math::ReadNonEmptyLine(f_properties, base_result_directory);
	fclose(f_properties);

	wchar_t _tmp_wc[1000];
	math::Char_To_Wchar_t(base_result_directory, _tmp_wc, 500);
	CreateDirectory((LPCTSTR)_tmp_wc, NULL);

	{
		printf_s("\n=============================*===============================\n");

		for (int i = 0; i < _materials.size(); i++)
		{
			solver_grid.AddDomain();
			auto domain = solver_grid.GetDomain(i);
			domain->forElectrical.sigma_tensor = _materials[i]._sigma_tensor;
		}

		std::function<Tensor2Rank3D(int, Point<double>)> StiffnessCoef = [&](int elem, Point<double> X)->Tensor2Rank3D
		{
			return solver_grid.GetDomain(solver_grid.GetElement(elem)->GetIdDomain())->forElectrical.sigma_tensor;
		};

		FEM::FEM_forBaseEliptic_tensor(
			is_print_logFile,
			critical_residual,
			geo_grid, //input
			first_boundaries, //input
			second_boundary, //input
			StiffnessCoef, //input
			base_result_directory, //output
			solver_grid, //output
			Solution //output
		);

		printf("Solve effective resistans... \n");
		double R = 0, R_new = 0, norm_Jsqrt = 0, rpho_eff_summ, GLOBAL_SIGMA_hmm;
		double deltaU = abs(first_boundaries[0].value(Point<double>(0, 0, 0)) - first_boundaries[1].value(Point<double>(0, 0, 0)));
		std::vector<double> I_in_sections;
		double R_sections = EffectiveResistance_ByСrossSection_forTensorDomains(typeXYZ_effective, typeGEO_section, deltaU, SectionPosition, solver_grid, Solution, I_in_sections);

		std::vector<double> I_in_bases, S_in_bases;
		double R_bases = EffectiveResistance_ByBoundary_forTensorDomains(
			typeXYZ_effective,
			boundary_faces,
			id_boundaries,
			deltaU,
			solver_grid,
			Solution,
			I_in_bases, S_in_bases);
		printf("\tcomplit\n");

		clock_t t_before = clock();

		printf("TIME = %lf sec\n", (t_before - t_after) * 1.0 / CLK_TCK);
		double end = omp_get_wtime();
		printf("start = %.16g\nend = %.16g\ndiff = %.16g\n", start, end, end - start);

		FILE* fout;
		FILE* fout_new;
		char name_out[1000];
		sprintf_s(name_out, "%s/effective_sigma.txt", base_result_directory);
		fopen_s(&fout, name_out, "a");

		for (int i = 0; i < solver_grid.GetDomainsCount(); i++)
			fprintf_s(fout, "sigma%d = %s\n",
				i, solver_grid.GetDomain(i)->forElectrical.sigma_tensor.print_xx_yy_zz_xy_xz_yz());

		char type;
		switch (typeXYZ_effective)
		{
		case 0: type = 'X'; break;
		case 1: type = 'Y'; break;
		case 2: type = 'Z'; break;
		}
		fprintf_s(fout, "By sections: sigma(%c)_eff = %.15lf\trpho(%c)_eff = %.15lf\n", type, 1. / R_sections, type, R_sections);
		fprintf_s(fout, "By bases: sigma(%c)_eff = %.15lf\trpho(%c)_eff = %.15lf\n", type, 1. / R_bases, type, R_bases);

		fprintf_s(fout, "\nIntegration results:\n");
		fprintf_s(fout, "\tCross sections:\n");
		fprintf_s(fout, "\t\tPositions:");
		for (int i = 0; i < SectionPosition.size(); i++)
		{
			fprintf_s(fout, "\t%.2lf", SectionPosition[i]);
		}
		fprintf_s(fout, "\n\t\tIntegral:");
		for (int i = 0; i < SectionPosition.size(); i++)
		{
			fprintf_s(fout, "\t%.4e", I_in_sections[i]);
		}
		fprintf_s(fout, "\n\n");

		fprintf_s(fout, "\tBases:\n");
		fprintf_s(fout, "\n\t\tArea:");
		for (int i = 0; i < I_in_bases.size(); i++)
		{
			fprintf_s(fout, "\t%.4e", S_in_bases[i]);
		}
		fprintf_s(fout, "\n\t\tIntegral:");
		for (int i = 0; i < I_in_bases.size(); i++)
		{
			fprintf_s(fout, "\t%.4e", I_in_bases[i]);
		}
		fprintf_s(fout, "\n\n");

		double V = 0;
		std::vector<double> Vin(solver_grid.GetDomainsCount());
		int qwe = 0;
		double conc = 0.0;

		for (int i = 0; i < solver_grid.GetElementsCount(); i++)
		{
			auto element = solver_grid.GetElement(i);
			Vin[element->GetIdDomain()] += element->GetVolume();
			V += element->GetVolume();
		}
		for (int i = 0; i < Vin.size(); i++)
		{
			printf("concentration(id%d) = %.5lf %%\n", i, (Vin[i] / V) * 100.);
			fprintf_s(fout, "concentration(id%d) = %.5lf %%\n", i, (Vin[i] / V) * 100.);
		}

		fprintf_s(fout, "TIME = %lf sec\n\n", (t_before - t_after) * 1.0 / CLK_TCK);
		fclose(fout);
		printf_s("\t complite\n");

		//output solution
		printf_s("Print the mech result into .dat file... ");
		{
			FILE* fout_tech;
			char name_u_tech[5000];
			sprintf_s(name_u_tech, "%s/U.dat", base_result_directory);
			fopen_s(&fout_tech, name_u_tech, "w");
			char name_in_file[1000];
			sprintf_s(name_in_file, "ScalarPotential");
			std::vector<std::vector<char>> name_value(5);
			char name_v_tmp[5][100];
			sprintf_s(name_v_tmp[0], "U");
			sprintf_s(name_v_tmp[1], "Jx");
			sprintf_s(name_v_tmp[2], "Jy");
			sprintf_s(name_v_tmp[3], "Jz");
			sprintf_s(name_v_tmp[4], "material");
			for (int i = 0; i < name_value.size(); i++)
			{
				name_value[i].resize(100);
				for (int j = 0; j < name_value[i].size(); j++)
				{
					name_value[i][j] = name_v_tmp[i][j];
				}
			}
			std::vector<std::vector<double>> value(5);
			value[0].resize(solver_grid.GetElementsCount());
			value[1].resize(solver_grid.GetElementsCount());
			value[2].resize(solver_grid.GetElementsCount());
			value[3].resize(solver_grid.GetElementsCount());
			value[4].resize(solver_grid.GetElementsCount());
			double sigma_inv_max = 0;
			int elem_sigma_max = 0;
			for (int i = 0; i < solver_grid.GetElementsCount(); i++)
			{
				auto element = solver_grid.GetElement(i);

				Point<double> Centr = element->GetWeightCentr();

				double U = solver_grid.GetSolutionInPoint(i, Centr, Solution);
				Point<double> dU = solver_grid.GetDerevativeFromSolutionInPoint(i, Centr, Solution);
				Tensor2Rank3D sigma = solver_grid.GetDomain(element->GetIdDomain())->forElectrical.sigma_tensor;
				Point<double> J = sigma * dU;

				value[0][i] = U;
				value[1][i] = -J.x;
				value[2][i] = -J.y;
				value[3][i] = -J.z;
				value[4][i] = element->GetIdDomain();
			}
			//solver_grid.printTecPlot3D(fout_tech, value, name_value, name_in_file);
			solver_grid.printTecPlot3D_DiffDomains(fout_tech, value, name_value, name_in_file);

			fclose(fout_tech);
		}
		printf_s("\t complite\n");



		solver_grid.~Grid_forScal();
	}
}

void main()
{
	//EffectiveSigmaZ_3Dmsh_ManyObj_WithBoundaryTag_1DOF();
	//EffectiveSigmaZ_3Dmsh_ManyObj_WithBoundaryTag_2DOF();

	/*math::SimpleGrid circle;
	CreateGrid2D_circle(Point<double>(100, 100, 100), 50, 50, 10, circle);

	FILE* fout_tech;
	char name_u_tech[5000];
	sprintf_s(name_u_tech, "./circle.dat");
	fopen_s(&fout_tech, name_u_tech, "w");
	char name_in_file[1000];
	sprintf_s(name_in_file, "circle");

	circle.printTecPlot3D(fout_tech, name_in_file);*/

	/*math::SimpleGrid rectangle;
	rectangle.CreateGrid2D_rectangle(Point<double>(50, 50, 50), Point<double>(100, 100, 100), 10, 50);
	max(0, 1);

	FILE* fout_tech;
	char name_u_tech[5000];
	sprintf_s(name_u_tech, "./rectangle.dat");
	fopen_s(&fout_tech, name_u_tech, "w");
	char name_in_file[1000];
	sprintf_s(name_in_file, "rectangle");

	rectangle.printTecPlot3D(fout_tech, name_in_file);*/

	//char base_name[1000] = { "D:/testing23062022/testing19102022/param_for_solver" };
	char base_name[1000] = { "./param_for_solver" };
	//char base_name[1000] = {"D:/testing23062022/sand_init/param_for_solver"};
	char properties_file[1000];
	int I = 0;
	sprintf_s(properties_file, sizeof(properties_file), "%s_%d.txt", base_name, I);
	FILE* fparam;
	fopen_s(&fparam, properties_file, "r");
	while (fparam != NULL)
	{
		fclose(fparam);
		EffectiveSigmaZ_3Dmsh_ManyObj_WithBoundaryTag_1DOF(properties_file);
		//ElectroStaticTensor_ManyObj_1DOF(properties_file);
		I++;
		sprintf_s(properties_file, sizeof(properties_file), "%s_%d.txt", base_name, I);
		fopen_s(&fparam, properties_file, "r");
	}

	int a;
	scanf_s("%d", &a);
}