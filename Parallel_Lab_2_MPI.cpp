#include <mpi.h>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <omp.h>
#include <sstream>
#include <time.h>
#include <vector>
#include <string>
#include <algorithm>
#include <numeric>
#include <math.h>

#define MaxOmpThreads 8

bool isThroughLinkDownOfNode(int nodeId, int countOfRows, int countOfColumns, int nodeIndexOfRow, int nodeIndexOfColumn, int processId,
	int loop, int k1) {

	if ((nodeIndexOfColumn == (countOfColumns - 1)) || (nodeIndexOfRow == (countOfRows - 1))) {
		//std::cout << "Process " << processId << ", NodeId " << nodeId << " row: " << nodeIndexOfRow <<
		//	", column: " << nodeIndexOfColumn << ", marked: " << false << std::endl;
		return false;
	}

	int currentFixedId = nodeId - nodeIndexOfRow;
	bool result = ((currentFixedId % loop) != 0) && ((currentFixedId % loop) > k1 - 1) && ((currentFixedId % loop) < loop);

	//std::cout << "Process " << processId << ", NodeId " << nodeId << " row: " << nodeIndexOfRow <<
	//	", column: " << nodeIndexOfColumn << ", marked: " << result << std::endl;

	return result;
}

void createNodesOfGraph(int processId, int Px, int Py, int ib, int ie, int jb, int je, int countOfColumns, int countOfRows, int k1,
	int k2, std::vector<int>& IA, std::vector<int>& JA, std::map<int, std::vector<int>>& resultGraph, int& NLocal,
	int& NOwn, std::vector<int>& Part, std::vector<int>& G2L) {

	int loop = k1 + k2;
	std::map<int, int> haloGraph;

#pragma omp parallel 
	{
#pragma omp for
		for (int j = jb - 1; j <= je + 1; j++)
		{
			bool isHaloJ = false;

			if (j < 0 || j == countOfRows) {
				continue;
			}
			else {
				isHaloJ = j == (jb - 1) || j == (je + 1);
			}


			for (int i = ib - 1; i <= ie + 1; i++)
			{
				bool isHaloI = false;
				if (i < 0 || i == countOfColumns)
				{
					continue;
				}
				else
				{
					isHaloI = i == (ib - 1) || i == (ie + 1);
				}


				std::vector<int> nodesNeighbors;

				int nodeId = i + j * countOfColumns;
				int nodeIndexOfColumn = nodeId % countOfColumns;
				int nodeIndexOfRow = nodeId / countOfColumns;

				bool upExist = (nodeId - countOfColumns - 1) >= 0;
				bool leftNodeExist = (nodeIndexOfColumn > 0);
				bool rightNodeExist = (nodeIndexOfColumn < countOfColumns - 1);
				bool downNodeExist = (nodeIndexOfRow < countOfRows - 1);

				bool throughLinkDown = isThroughLinkDownOfNode(nodeId, countOfRows, countOfColumns, nodeIndexOfRow, nodeIndexOfColumn, processId, loop, k1);
				bool throughLinkUp = upExist && leftNodeExist && isThroughLinkDownOfNode(nodeId - countOfColumns - 1, countOfRows, countOfColumns, nodeIndexOfRow - 1, nodeIndexOfColumn - 1, processId, loop, k1);

				bool upNodeExist = (nodeIndexOfRow > 0);

				// Добавление боковых сверху
				if (throughLinkUp)
				{
					int throwingNodeId = nodeId - countOfColumns - 1;
					nodesNeighbors.push_back(throwingNodeId);
				}

				// Добавляем верхнюю вершину в соседи
				if (upNodeExist)
				{
					nodesNeighbors.push_back(nodeId - countOfColumns);
				}

				// Добавляем левую вершину в соседи
				if (leftNodeExist)
				{
					nodesNeighbors.push_back(nodeId - 1);
				}

				// Добавляем текущую вершину в соседи
				nodesNeighbors.push_back(nodeId);

				// Добавляем правую вершину в соседи
				if (rightNodeExist)
				{
					nodesNeighbors.push_back(nodeId + 1);
				}

				// Добавляем нижнюю вершину в соседи
				if (downNodeExist)
				{
					nodesNeighbors.push_back(nodeId + countOfColumns);
				}

				// Добавление боковых снизу
				if (throughLinkDown)
				{
					int throwingNodeId = nodeId + countOfColumns + 1;
					nodesNeighbors.push_back(throwingNodeId);
				}
				//

				// Угловые точки, не входящие в гало
				// XOR
				if (isHaloI || isHaloJ)
				{
					if (isHaloI ^ isHaloJ)
					{
						int partProcessId = processId;
						if (j == jb - 1) {
							partProcessId -= Px;
						}
						if (j == je + 1)
						{
							partProcessId += Px;
						}
						if (i == ib - 1)
						{
							partProcessId -= 1;
						}
						if (i == ie + 1)
						{
							partProcessId += 1;
						}

						haloGraph.insert(std::pair<int, int>(nodeId, partProcessId));
					}
				}
				else
				{
					resultGraph.insert(std::pair<int, std::vector<int>>(nodeId, nodesNeighbors));
				}
			}
		}
	}

	NOwn = resultGraph.size();
	NLocal = NOwn + haloGraph.size();
	IA.resize(NOwn);
	IA.push_back(0);

	int i = 0;
	// Заполнение портретов
	for (auto tempPair : resultGraph)
	{
		std::vector<int> nodesNeighbors = tempPair.second;
		IA.at(i + 1) = IA.at(i) + nodesNeighbors.size();

		for (int element : nodesNeighbors)
		{
			JA.push_back(element);
		}

		G2L.push_back(tempPair.first);
		Part.push_back(processId);

		i++;
	}

	for (auto tempPair : haloGraph)
	{
		G2L.push_back(tempPair.first);
		Part.push_back(tempPair.second);
	}


	//std::cout << "Process " << processId << " NO: " << NO << ", NH: " << NH << std::endl;

	/*if (processId == 1)
	{
		std::cout << "Own" << std::endl;

		for (auto temp : resultGraph) {
			std::cout << temp.first << std::endl;
		}

		std::cout << "Halo" << std::endl;

		for (auto temp : haloGraph) {
			std::cout << temp.first << std::endl;
		}
	}*/

	int s = 1;

}


void calculateSubAreaIndexes(int processId, int Px, int Py, int& ib, int& ie, int& jb, int& je, int Nx, int Ny) {

	std::vector<int> offsetX, offsetY;
	offsetX.resize(Py);
	offsetY.resize(Px);

	int processXId = processId % Py;
	int processYId = processId / Py;

	bool includeHalo = false;

	// По текущей строке
	if (offsetX.at(processYId) == 0)
	{
		offsetX.at(processYId) = Nx % Px;
		MPI_Bcast(&offsetX.at(processYId), 1, MPI_INT, 0, MPI_COMM_WORLD);
	}

	int countToOneProcessX = Nx / Px;
	bool remainderX = (processXId < offsetX.at(processYId));

	ib = processXId * countToOneProcessX + ((remainderX) ? processXId : 0);

	if (!remainderX) {
		ib += offsetX.at(processYId);
	}

	ie = ib + countToOneProcessX - ((remainderX) ? 0 : 1);

	if (includeHalo) {
		if (ib > 0) {
			ib--;
		}
		if (ie < Nx - 1)
		{
			ie++;
		}
	}

	// По текущему столбцу
	if (offsetY.at(processXId) == 0)
	{
		offsetY.at(processXId) = Ny % Py;
		MPI_Bcast(&offsetY.at(processXId), 1, MPI_INT, 0, MPI_COMM_WORLD);
	}
	int countToOneProcessY = Ny / Py;
	bool remainderY = (processYId < offsetY.at(processXId));

	jb = processYId * countToOneProcessY + ((remainderY) ? processYId : 0);

	if (!remainderY) {
		jb += offsetY.at(processXId);
	}

	je = jb + countToOneProcessY - ((remainderY) ? 0 : 1);

	if (includeHalo) {
		if (jb > 0) {
			jb--;
		}
		if (je < Ny - 1) {
			je++;
		}
	}
}

std::string vectorToString(std::vector<int> intVector)
{
	std::stringstream result;
	result << "[";

	for (int i = 0; i < intVector.size(); i++)
	{

		if (i == intVector.size() - 1)
		{
			result << std::to_string(intVector.at(i));
		}
		else
		{
			result << std::to_string(intVector.at(i)) << ", ";
		}
	}

	result << "]";

	return result.str();
}

std::string createAdjacencyList(std::map<int, std::vector<int>> graph)
{
	std::string result;

	for (std::pair<int, std::vector<int>> nodePair : graph)
	{
		std::stringstream resultString;
		resultString << "[" << nodePair.first << "] -> {";

		std::vector<int> neighbors = nodePair.second;
		for (int i = 0; i < neighbors.size(); i++)
		{
			int neighborID = neighbors[i];

			if (i == neighbors.size() - 1)
			{
				resultString << std::to_string(neighborID);
			}
			else
			{
				resultString << std::to_string(neighborID) << ", ";
			}
		}

		resultString << "}" << std::endl;
		result += resultString.str();
	}

	return result;
}

void printResult(std::map<int, std::vector<int>> graph, std::vector<int> IA, std::vector<int> JA, double seconds, int countOfNodes,
	int processId, std::vector<int> G2L, std::vector<int> Part)
{
	std::cout << "Process " << processId << ": N (size of matrix) - " << countOfNodes << " nodes" << std::endl;

	std::string IAstring = vectorToString(IA);
	std::cout << "Process " << processId << ": IA: " << IAstring << std::endl;

	std::string JAstring = vectorToString(JA);
	std::cout << "Process " << processId << ": JA: " << JAstring << std::endl;

	std::string adjancecyString = createAdjacencyList(graph);
	std::cout << "Process " << processId << ": Adjacency list: " << std::endl;
	std::cout << adjancecyString << std::endl;

	std::string G2LString = vectorToString(G2L);
	std::cout << "Process " << processId << ": G2L: " << G2LString << std::endl;

	std::string PartString = vectorToString(Part);
	std::cout << "Process " << processId << ": Part: " << PartString << std::endl;

	//std::cout << "Execution time: " << seconds << " s." << std::endl;
	//std::cout << "Time per 1 node: " << seconds / countOfNodes << " s." << std::endl;
}

// Проверка, что Px * Py = P, указанному при запуске программы
bool isCountOfProcessesValid(int Px, int Py, int countOfProcesses, int processId) {
	if (countOfProcesses != Px * Py) {
		if (processId == 0) {
			std::cout << "The number of MPI processes does not match the mesh decomposition parameters. (Количество MPI процессов не совпадает с параметрами декомпозции сетки.)";
			//MPI_Abort(MPI_COMM_WORLD, -1);
		}
		return false;
	}

	return true;
}

// Инициализация стартовых переменных
bool initVariables(int argc, char* argv[], int& Nx, int& Ny, int& k1, int& k2, int& Px, int& Py, double& tol,
	std::vector<int>& arguments, bool& isPrint) {
	// Первый параметр - ссылка на сборку
	for (int i = 1; i < argc; i++)
	{
		try
		{
			if (i == 7)
			{
				tol = std::stod(argv[i]);
			}
			else
			{
				arguments.push_back(atoi(argv[i]));
			}
		}
		catch (const std::exception)
		{
			std::cout << "Incorrect entry of launch parameters (Некорректный ввод параметров запуска)" << std::endl;
		}
	}

	if (arguments.empty() || arguments.size() < 7)
	{
		std::cout << "Enter 8 elements of the graph portrait (Введите 8 элементов портрета графа) (Nx, Ny, k1, k2, Px, Py, tol, isPrint)!" << std::endl;
		return false;
	}

	Nx = arguments[0] + 1;
	Ny = arguments[1] + 1;
	k1 = arguments[2];
	k2 = arguments[3];
	Py = arguments[4];
	Px = arguments[5];

	if (arguments.size() > 5)
	{
		isPrint = (arguments[6] == 1) ? true : false;
	}

	return true;
}

int main(int argc, char* argv[])
{
	setlocale(LC_ALL, "Russian");

	omp_set_num_threads(MaxOmpThreads);
	int Nx, Ny, Py, Px, k1, k2, NLocal, NOwn;
	int countOfProcesses, processId;
	// Точность решения
	double tol;
	bool isPrint;
	std::vector<int> arguments, IA, JA, offsetX, offsetY, Part, G2L;
	std::map<int, std::vector<int>> resultGraph;
	// Диапазоны индексов подобласти каждого из процесса
	int ib, ie, jb, je;

	if (!initVariables(argc, argv, Nx, Ny, k1, k2, Px, Py, tol, arguments, isPrint)) {
		return 0;
	}

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &processId);
	MPI_Comm_size(MPI_COMM_WORLD, &countOfProcesses);

	if (!isCountOfProcessesValid(Px, Py, countOfProcesses, processId)) {
		return 0;
	}

	calculateSubAreaIndexes(processId, Px, Py, ib, ie, jb, je, Nx, Ny);
	
	double start = omp_get_wtime();
	createNodesOfGraph(processId, Px, Py, ib, ie, jb, je, Nx, Ny, k1, k2, IA, JA, resultGraph, NLocal, NOwn, Part, G2L);
	double end = omp_get_wtime();
	double seconds = end - start;
	std::cout << "Process " << processId << ": Total nodes: " << NOwn << " elem." << std::endl;
	std::cout << "Process " << processId << ": First stage time: " << seconds << " s." << std::endl;
	if (isPrint)
	{
		printResult(resultGraph, IA, JA, seconds, NOwn, processId, G2L, Part);
	}

	printf("Process = %d, ib: %d, ie: %d\, jb: %d, je: %d\n", processId, ib, ie, jb, je);

	MPI_Finalize();

	return 0;
}
