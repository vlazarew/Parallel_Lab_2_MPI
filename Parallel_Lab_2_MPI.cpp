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


#pragma region Этап 1. Генерация портрета

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

void createNodesOfGraph(int processId, int processesColumns, int processesRows, int columnIndexBegin, int columnIndexEnd, int rowIndexBegin, int rowIndexEnd, int countOfColumns, int countOfRows, int k1,
	int k2, std::vector<int>& IA, std::vector<int>& JA, std::map<int, std::vector<int>>& resultGraph, int& NLocal,
	int& NOwn, std::vector<int>& Part, std::vector<int>& G2L, std::vector<int>& L2G) {

	int loop = k1 + k2;
	std::map<int, int> haloGraph;
	int localIndex = 0;

#pragma omp parallel 
	{
#pragma omp for
		for (int j = rowIndexBegin - 1; j <= rowIndexEnd + 1; j++)
		{
			bool isHaloJ = false;

			if (j < 0 || j == countOfRows) {
				continue;
			}
			else {
				isHaloJ = j == (rowIndexBegin - 1) || j == (rowIndexEnd + 1);
			}


			for (int i = columnIndexBegin - 1; i <= columnIndexEnd + 1; i++)
			{
				bool isHaloI = false;
				if (i < 0 || i == countOfColumns)
				{
					continue;
				}
				else
				{
					isHaloI = i == (columnIndexBegin - 1) || i == (columnIndexEnd + 1);
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
						if (j == rowIndexBegin - 1) {
							partProcessId -= processesColumns;
						}
						if (j == rowIndexEnd + 1)
						{
							partProcessId += processesColumns;
						}
						if (i == columnIndexBegin - 1)
						{
							partProcessId -= 1;
						}
						if (i == columnIndexEnd + 1)
						{
							partProcessId += 1;
						}

#pragma omp critical
						haloGraph.insert(std::pair<int, int>(nodeId, partProcessId));
					}
				}
				else
				{
#pragma omp critical
					resultGraph.insert(std::pair<int, std::vector<int>>(nodeId, nodesNeighbors));
				}
			}
		}
	}

	G2L.resize(countOfColumns * countOfRows);
#pragma omp parallel
	{

#pragma omp for
		for (int i = 0; i < G2L.size(); i++)
		{
			G2L.at(i) = -1;
		}
	}

	NOwn = resultGraph.size();
	NLocal = NOwn + haloGraph.size();
	IA.resize(NOwn);
	IA.push_back(0);

	int i = 0;
	for (auto tempPair : resultGraph) {
		L2G.push_back(tempPair.first);
		Part.push_back(processId);
		G2L.at(L2G.at(i)) = i;
		i++;
	}

	for (auto tempPair : haloGraph)
	{
		L2G.push_back(tempPair.first);
		Part.push_back(tempPair.second);
		G2L.at(L2G.at(i)) = i;
		i++;
	}

	i = 0;
	// Заполнение портретов
	for (auto tempPair : resultGraph)
	{
		int counterIA = 0;
		for (int element : tempPair.second)
		{
			int localIndexOfNode = G2L.at(element);
			if (localIndexOfNode <= NLocal && localIndexOfNode >= 0)
			{
				JA.push_back(localIndexOfNode);
				counterIA++;
			}
		}

		IA.at(i + 1) = IA.at(i) + counterIA;
		i++;
	}
}

void calculateSubAreaIndexes(int processId, int processesColumns, int processesRows, int& columnIndexBegin, int& columnIndexEnd,
	int& rowIndexBegin, int& rowIndexEnd, int countOfColumns, int countOfRows) {

	std::vector<int> offsetColumns, offsetRows;
	offsetColumns.resize(processesColumns);
	offsetRows.resize(processesRows);

	int processColumnId = processId % processesColumns;
	int processRowId = processId / processesColumns;

	bool includeHalo = false;

	// По текущей строке
	if (offsetColumns.at(processColumnId) == 0)
	{
		offsetColumns.at(processColumnId) = countOfColumns % processesColumns;
		MPI_Bcast(&offsetColumns.at(processColumnId), 1, MPI_INT, 0, MPI_COMM_WORLD);
	}

	int countToOneProcessX = countOfColumns / processesColumns;
	bool remainderX = (processColumnId < offsetColumns.at(processColumnId));

	columnIndexBegin = processColumnId * countToOneProcessX + ((remainderX) ? processColumnId : 0);

	if (!remainderX) {
		columnIndexBegin += offsetColumns.at(processColumnId);
	}

	columnIndexEnd = columnIndexBegin + countToOneProcessX - ((remainderX) ? 0 : 1);

	/*if (includeHalo) {
		if (ib > 0) {
			ib--;
		}
		if (ie < Nx - 1)
		{
			ie++;
		}
	}*/

	// По текущему столбцу
	if (offsetRows.at(processRowId) == 0)
	{
		offsetRows.at(processRowId) = countOfRows % processesRows;
		MPI_Bcast(&offsetRows.at(processRowId), 1, MPI_INT, 0, MPI_COMM_WORLD);
	}
	int countToOneProcessY = countOfRows / processesRows;
	bool remainderY = (processRowId < offsetRows.at(processRowId));

	rowIndexBegin = processRowId * countToOneProcessY + ((remainderY) ? processRowId : 0);

	if (!remainderY) {
		rowIndexBegin += offsetRows.at(processRowId);
	}

	rowIndexEnd = rowIndexBegin + countToOneProcessY - ((remainderY) ? 0 : 1);

	/*if (includeHalo) {
		if (jb > 0) {
			jb--;
		}
		if (je < Ny - 1) {
			je++;
		}
	}*/
}

std::string vectorToString(std::vector<int> intVector, std::vector<int> L2G, bool isGlobal = false)
{
	std::stringstream result;
	result << "[";

	for (int i = 0; i < intVector.size(); i++)
	{
		int currenntElem;
		if (isGlobal)
		{
			currenntElem = L2G.at(intVector.at(i));
		}
		else
		{
			currenntElem = intVector.at(i);
		}

		if (i == intVector.size() - 1)
		{
			result << std::to_string(currenntElem);
		}
		else
		{
			result << std::to_string(currenntElem) << ", ";
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

std::string printResult(std::map<int, std::vector<int>> graph, std::vector<int> IA, std::vector<int> JA, double seconds,
	int countOfNodes, int processId, std::vector<int> G2L, std::vector<int> L2G, std::vector<int> Part)
{
	std::stringstream result;

	result << "Process " << processId << ": N (size of matrix) - " << countOfNodes << " nodes" << std::endl;

	std::string IAstring = vectorToString(IA, L2G);
	result << "Process " << processId << ": IA: " << IAstring << std::endl;

	std::string JAstring = vectorToString(JA, L2G);
	result << "Process " << processId << ": JA: " << JAstring << std::endl;

	std::string JAGlobalString = vectorToString(JA, L2G, true);
	result << "Process " << processId << ": JA (Global): " << JAGlobalString << std::endl;

	std::string adjancecyString = createAdjacencyList(graph);
	result << "Process " << processId << ": Adjacency list: " << std::endl << adjancecyString << std::endl;

	std::string G2LString = vectorToString(G2L, L2G);
	result << "Process " << processId << ": G2L: " << G2LString << std::endl;

	std::string L2GString = vectorToString(L2G, L2G);
	result << "Process " << processId << ": L2G: " << L2GString << std::endl;

	std::string PartString = vectorToString(Part, L2G);
	result << "Process " << processId << ": Part: " << PartString << std::endl;

	//std::cout << "Execution time: " << seconds << " s." << std::endl;
	//std::cout << "Time per 1 node: " << seconds / countOfNodes << " s." << std::endl;

	return result.str();
}

// Проверка, что Px * Py = P, указанному при запуске программы
bool isCountOfProcessesValid(int processesColumns, int processesRows, int countOfProcesses, int processId) {
	if (countOfProcesses != processesColumns * processesRows) {
		if (processId == 0) {
			std::cout << "The number of MPI processes does not match the mesh decomposition parameters. (Количество MPI процессов не совпадает с параметрами декомпозции сетки.)";
			//MPI_Abort(MPI_COMM_WORLD, -1);
		}
		return false;
	}

	return true;
}

// Инициализация стартовых переменных
bool initVariables(int argc, char* argv[], int& countOfColumns, int& countOfRows, int& k1, int& k2, int& processesColumns,
	int& processesRows, double& tol, std::vector<int>& arguments, bool& isPrint) {
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

	countOfColumns = arguments[0] + 1;
	countOfRows = arguments[1] + 1;
	k1 = arguments[2];
	k2 = arguments[3];
	processesColumns = arguments[4];
	processesRows = arguments[5];

	if (arguments.size() > 5)
	{
		isPrint = (arguments[6] == 1) ? true : false;
	}

	return true;
}

#pragma endregion

#pragma region Этап 2. Генерация СЛАУ

void makeSLAE(std::vector<int> IA, std::vector<int> JA, std::vector<double>& A, std::vector<double>& b, int NOwn, int NLocal,
	std::vector<int> L2G)
{
	// i - номер строки, j - номер столбца
#pragma omp parallel
	{
#pragma omp for
		for (int i = 0; i < NOwn; i++)
		{
			double rowSum = 0;
			int diagonalIndex = 0;
			int startIndex = IA.at(i);
			int endIndex = IA.at(i + 1);

			for (int j = startIndex; j < endIndex; ++j)
			{
				int indexOfNode = JA.at(j);
				if (i == indexOfNode)
				{
					diagonalIndex = j;
					continue;
				}
				else
				{
					double value = cos(i * L2G.at(indexOfNode) + i + L2G.at(indexOfNode));
					A.at(j) = value;
					rowSum = rowSum + std::fabs(value);
				}
			}
			A.at(diagonalIndex) = 1.234 * rowSum;
			b.at(i) = sin(L2G.at(i));
		}
	}
}

std::string printSLAE(std::vector<double> A, std::vector<double> b, std::vector<int> IA, int NOwn, int NLocal, int processId)
{
	std::stringstream resultString;

	resultString << "Process " << processId << ": Matrix and right-hand side coefficients: " << std::endl;

	int countOfLinks = 0;
	int currentIndex = 0;
	for (int i = 1; i <= NOwn; i++)
	{
		countOfLinks = IA.at(i) - IA.at(i - 1);

		resultString << std::to_string(i - 1) << ". [";
		for (int j = 0; j < countOfLinks; j++)
		{
			resultString << std::fixed << std::showpoint << std::setprecision(5) << A.at(currentIndex);
			currentIndex++;
			if (j != countOfLinks - 1)
			{
				resultString << ", ";
			}
		}

		resultString << "]";

		resultString << " = " << std::fixed << std::showpoint << std::setprecision(5) << b.at(i - 1) << std::endl;
	}

	return resultString.str();
}


#pragma endregion

#pragma region Этап 3. Построение схемы обменов

void createCom(int NLocal, int NOwn, std::vector<int> IA, std::vector<int> JA, std::vector<int> Part, std::vector<int> L2G,
	std::vector<int> G2L, std::vector<int>& Neighbours, std::vector<int>& SendOffet, std::vector<int>& RecvOffset,
	std::vector<int>& Send, std::vector<int>& Recv, int countOfProcesses) {

	std::vector<std::vector<int>> SendToProcess;
	std::vector<std::vector<int>> RecvFromProcess;
	SendToProcess.resize(countOfProcesses);
	RecvFromProcess.resize(countOfProcesses);

	for (int i = 0; i < NOwn; i++)
	{
		int startIndex = IA.at(i);
		int endIndex = IA.at(i + 1);

		for (int j = startIndex; j < endIndex; ++j)
		{
			int indexOfNode = JA.at(j);

			// Значит попали в вершину, которая от другого процесса
			if (indexOfNode >= NOwn && indexOfNode < NLocal)
			{
				int haloProcess = Part.at(indexOfNode);
				SendToProcess.at(haloProcess).push_back(i);
				RecvFromProcess.at(haloProcess).push_back(indexOfNode);
			}

		}
	}

	for (int haloProcessId = 0; haloProcessId < countOfProcesses; haloProcessId++)
	{
		std::vector<int> SendToProcessHalo = SendToProcess.at(haloProcessId);
		std::vector<int> RecvFromProcessHalo = RecvFromProcess.at(haloProcessId);

		if (SendToProcessHalo.empty() && RecvFromProcessHalo.empty())
		{
			continue;
		}

		Neighbours.push_back(haloProcessId);

		// SendToProcessHalo.size() и RecvFromProcessHalo.size() должны быть равны
		for (int i = 0; i < SendToProcessHalo.size(); i++) {
			SendToProcessHalo.at(i) = L2G.at(SendToProcessHalo.at(i));
			RecvFromProcessHalo.at(i) = L2G.at(RecvFromProcessHalo.at(i));
		}

		// Сортируем по возрастанию глобальных номеров
		std::sort(SendToProcessHalo.begin(), SendToProcessHalo.end());
		std::sort(RecvFromProcessHalo.begin(), RecvFromProcessHalo.end());

		// Удаляем возможные дубликаты
		SendToProcessHalo.resize(std::unique(SendToProcessHalo.begin(), SendToProcessHalo.end()) - SendToProcessHalo.begin());
		RecvFromProcessHalo.resize(std::unique(RecvFromProcessHalo.begin(), RecvFromProcessHalo.end()) - RecvFromProcessHalo.begin());

		// Конвертируем обратно в локальные
		for (int i = 0; i < SendToProcessHalo.size(); i++) {
			SendToProcessHalo.at(i) = G2L.at(SendToProcessHalo.at(i));
			RecvFromProcessHalo.at(i) = G2L.at(RecvFromProcessHalo.at(i));
		}

		// Добавляем все в вектора Send и Recv
		for (int i = 0; i < SendToProcessHalo.size(); i++) {
			Send.push_back(SendToProcessHalo.at(i));
			Recv.push_back(RecvFromProcessHalo.at(i));
		}

		SendOffet.push_back(Send.size());
		RecvOffset.push_back(Recv.size());
	}

}

std::string printCom(std::vector<int> Neighbours, std::vector<int> SendOffet, std::vector<int> RecvOffset,
	std::vector<int> Send, std::vector<int> Recv, int processId, std::vector<int> L2G) {

	std::stringstream result;

	std::string NeighboursString = vectorToString(Neighbours, L2G);
	result << "Process " << processId << ": Neighbours: " << NeighboursString << std::endl;

	std::string SendOffsetString = vectorToString(SendOffet, L2G);
	result << "Process " << processId << ": SendOffset: " << SendOffsetString << std::endl;

	std::string RecvOffsetString = vectorToString(RecvOffset, L2G);
	result << "Process " << processId << ": RecvOffset: " << RecvOffsetString << std::endl;

	std::string SendString = vectorToString(Send, L2G);
	result << "Process " << processId << ": Send: " << SendString << std::endl;

	std::string SendStringGlobal = vectorToString(Send, L2G, true);
	result << "Process " << processId << ": Send (Global): " << SendStringGlobal << std::endl;

	std::string RecvString = vectorToString(Recv, L2G);
	result << "Process " << processId << ": Recv: " << RecvString << std::endl;

	std::string RecvStringGlobal = vectorToString(Recv, L2G, true);
	result << "Process " << processId << ": Recv (Global): " << RecvStringGlobal << std::endl;

	return result.str();
}

#pragma endregion

#pragma region Этап 4. Решение СЛАУ

double scalar(std::vector<double> x1, std::vector<double> x2, double& allTime, int& countOfCalls)
{
	double start = omp_get_wtime();

	double result = 0;
	std::vector<double> resultVector;
	std::vector<double> mpiResult;
	int vectorSize = x1.size();
	resultVector.resize(vectorSize);

#pragma omp parallel
	{
#pragma omp for
		for (int i = 0; i < vectorSize; i++)
		{
			resultVector.at(i) = x1.at(i) * x2.at(i);
		}
	}

	double processResult = 0;
	for (int i = 0; i < vectorSize; i++) {
		processResult += resultVector.at(i);
	}

	MPI_Allreduce(&processResult, &result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	double end = omp_get_wtime();
	allTime += (end - start);
	countOfCalls++;

	return result;
}

// Расчет L2 нормы вектора
double normalizeVector(std::vector<double> x)
{
	double result = 0;

#pragma omp parallel
	{
#pragma omp for
		for (int i = 0; i < x.size(); i++)
		{
			result += x.at(i) * x.at(i);
		}
	}

	return std::sqrt(result);
}

std::vector<double> spMV(int countOfNodes, std::vector<int> IA, std::vector<int> JA, std::vector<double> A, std::vector<double> x,
	double& allTime, int& countOfCalls, std::vector<int> Send, std::vector<int> Recv, std::vector<int> Neighbours, std::vector<int> SendOffset, std::vector<int> RecvOffset, int processId)
{
	double start = omp_get_wtime();

	std::vector<double> xSend(countOfNodes);
	std::map<int, double> xRecv;
	int total = Send.size() + Recv.size();
	std::vector<MPI_Request> request(total);

	//std::cout << "req" << request.size() << std::endl;
	//std::cout << "send" << Send.size() << std::endl;

	//if (processId == 0) {
	//std::cout << "process = " << processId << "neighbours size = " << Neighbours.size() << std::endl;
	//std::cout << "process = " << processId << "sendoffset size = " << SendOffset.size() << std::endl;
	for (int i = 0; i < Neighbours.size(); i++) {
		//std::cout << "process = " << processId << "sendoffset(i) = " << SendOffset.at(i) << "sendoffset(i+1) = " << SendOffset.at(i + 1) << std::endl;
		for (int j = SendOffset.at(i); j < SendOffset.at(i + 1); j++)
		{
			//std::cout << "process = " << processId << "Send.at(j) = " << Send.at(j) << std::endl;
			xSend.at(Send.at(j)) = x.at(Send.at(j));
			//for (int i = 0; i < xSend.size(); i++) {
			//	std::cout << xSend.at(Send.at(j)) << std::endl;
			//}
			//for (int i = 0; i < SendOffset.size(); i++) {
			//	std::cout << "offset"<< SendOffset.at(i) << std::endl;
			//std::cout << "govno" << xSend.at(Send.at(j)) << std::endl;
			//std::cout << "process = " << processId << "i = " << i << ", j = " << j << std::endl;
			//}
			MPI_Isend(&xSend.at(Send.at(j)), 1, MPI_DOUBLE, Neighbours.at(i), 0, MPI_COMM_WORLD, &request.at(j));
		}
		//std::cout << "process = " << processId << "recvoffset  i= " << SendOffset.at(i) << "recvoffset i+1 = " << SendOffset.at(i + 1) << std::endl;
		for (int j = RecvOffset.at(i); j < RecvOffset.at(i + 1); j++)
		{
			//std::cout << "process = " << processId << "Send.at(j) = " << Send.at(j) << std::endl;
			//xRecv.at(Recv.at(j)) = 0;
			xRecv[Recv.at(j)] = 0;
			MPI_Irecv(&xRecv.at(Recv.at(j)), 1, MPI_DOUBLE, Neighbours.at(i), 0, MPI_COMM_WORLD, &request.at(Send.size() + j));
		}
	}
	//}

	MPI_Waitall(total, request.data(), MPI_STATUSES_IGNORE);
	std::vector<double> result(countOfNodes);

	for (int i = 0; i < countOfNodes; i++)
	{
		result.at(i) = 0;
		int end = IA.at(std::min(i + 1, countOfNodes));
		for (int j = IA.at(i); j < end; j++)
		{
			if (JA.at(j) < countOfNodes)
			{
				result.at(i) += A.at(j) * x.at(JA.at(j));
			}
			else {
				result.at(i) += A.at(j) * xRecv.at(JA.at(j));
			}
		}
	}

	/*if (x.size() == 0)
	{
		x.resize(countOfNodes);
	}

#pragma omp parallel
	{
#pragma omp for
		for (int i = 0; i < countOfNodes; i++)
		{
			for (int k = IA.at(i); k < IA.at(i + 1); k++)
			{
				if (JA.at(k) < countOfNodes)
				{
					result.at(i) += A.at(k) * x.at(JA.at(k));
				}
			}
		}
	}*/

	double end = omp_get_wtime();
	allTime += (end - start);
	countOfCalls++;

	return result;
}

std::vector<double> linearCombination(std::vector<double> x1, std::vector<double> x2, double a1, double a2,
	double& allTime, int& countOfCalls)
{
	double start = omp_get_wtime();

	std::vector<double> result;
	int vectorSize = x1.size();
	result.resize(vectorSize);

#pragma omp parallel
	{
#pragma omp for
		for (int i = 0; i < vectorSize; i++)
		{
			result.at(i) = x1.at(i) * a1 + x2.at(i) * a2;
		}
	}

	double end = omp_get_wtime();
	allTime += (end - start);
	countOfCalls++;

	return result;
}

void createMMatrixFromA(int countOfNodes, std::vector<int> IA, std::vector<int> JA, std::vector<double> A,
	std::vector<int>& IAM, std::vector<int>& JAM, std::vector<double>& AM)
{
	IAM.resize(countOfNodes + 1);
	AM.resize(countOfNodes);
	JAM.resize(countOfNodes);
	IAM.at(0) = 0;

#pragma omp parallel
	{
#pragma omp for
		for (int i = 0; i < countOfNodes; i++)
		{
			int startIndex = IA.at(i);
			int endIndex = IA.at(i + 1);

			for (int j = startIndex; j < endIndex; ++j)
			{
				int indexOfNode = JA.at(j);
				if (i == indexOfNode)
				{
					AM.at(i) = A.at(j);
					IAM.at(i + 1) = i + 1;
					JAM.at(i) = i;
				}
			}
		}
	}
}

void reverseMMatrix(std::vector<double>& AM)
{
#pragma omp parallel
	{
#pragma omp for
		for (int i = 0; i < AM.size(); i++)
		{
			AM.at(i) = 1 / AM.at(i);
		}
	}
}

std::string solveSLAE(std::vector<int> IA, std::vector<int> JA, std::vector<double> A, std::vector<double> b, int NOwn, int NLocal,
	double tol, std::vector<double>& xRes, int& n, double& res, std::vector<int> Neighbours, std::vector<int> SendOffet,
	std::vector<int> RecvOffset, std::vector<int> Send, std::vector<int> Recv, int processId)
{
	std::stringstream result;

	//if (processId == 0) {
	double allTimeSpMV = 0;
	double allTimeLinear = 0;
	double allTimeScalar = 0;

	int countOfCallsSpMV = 0;
	int countOfCallsLinear = 0;
	int countOfCallsScalar = 0;

	double normB = normalizeVector(b);

	std::vector<double> zCurr(NOwn);
	std::vector<double> pPrev;
	std::vector<double> pCurr;
	std::vector<double> qCurr;
	double poPrev;
	double poCurr;
	double bettaCurr;
	double alphaCurr;

	std::vector<double> xPrev;
	// Вектор с невязкой
	std::vector<double> rPrev;
	std::vector<double> rCurr;
	xPrev.resize(NOwn);


	std::vector<double> AX0 = spMV(NOwn, IA, JA, A, xPrev, allTimeSpMV, countOfCallsSpMV, Send, Recv, Neighbours, SendOffet, RecvOffset, processId);
	rPrev = linearCombination(b, AX0, 1, -1, allTimeLinear, countOfCallsLinear);
	//rPrev = b;

	bool convergence = false;
	int k = 1;
	std::vector<int> IAM;
	std::vector<int> JAM;
	std::vector<double> AM(NOwn);
	createMMatrixFromA(NOwn, IA, JA, A, IAM, JAM, AM);
	reverseMMatrix(AM);

	/*for (int i = 0; i < NOwn; i += 1) {
		//solution.x[i] = 0;
		//r[i] = area.b[i];
		int d = IA.at(i);
		while (JA.at(d) != i) {
			d += 1;
		}
		AM.at(i) = 1.0 / A.at(d);
	}*/


	do
	{
		if (k == NOwn) {
			convergence = true;
			continue;
		}

		/*for (int i = 0; i < NOwn; i += 1) {
			zCurr.at(i) = AM.at(i) * rPrev.at(i);
		}*/

		zCurr = spMV(NOwn, IAM, JAM, AM, rPrev, allTimeSpMV, countOfCallsSpMV, Send, Recv, Neighbours, SendOffet, RecvOffset, processId);
		poCurr = scalar(rPrev, zCurr, allTimeScalar, countOfCallsScalar);

		if (k == 1)
		{
			pCurr = zCurr;
		}
		else
		{
			bettaCurr = poCurr / poPrev;
			pCurr = linearCombination(zCurr, pPrev, 1, bettaCurr, allTimeLinear, countOfCallsLinear);
		}

		qCurr = spMV(NOwn, IA, JA, A, pCurr, allTimeSpMV, countOfCallsSpMV, Send, Recv, Neighbours, SendOffet, RecvOffset, processId);
		alphaCurr = poCurr / scalar(pCurr, qCurr, allTimeScalar, countOfCallsScalar);
		xRes = linearCombination(xPrev, pCurr, 1, alphaCurr, allTimeLinear, countOfCallsLinear);

		rCurr = linearCombination(rPrev, qCurr, 1, -alphaCurr, allTimeLinear, countOfCallsLinear);

		res = normalizeVector(rCurr);
		result << "Step " << k << " ||b - Ax|| = " << res << " po = " << poCurr << std::endl;
		if (poCurr < tol || k >= 15000)
		{
			convergence = true;
		}
		else
		{
			k++;
			poPrev = poCurr;
			pPrev = pCurr;
			xPrev = xRes;
			rPrev = rCurr;
		}
	} while (!convergence);

	result << "Average time scalar: " << (allTimeScalar / countOfCallsScalar) << " s." << std::endl;
	result << "Average time linear: " << (allTimeLinear / countOfCallsLinear) << " s." << std::endl;
	result << "Average time SpMV: " << (allTimeSpMV / countOfCallsSpMV) << " s." << std::endl;
	result << "Total time scalar: " << (allTimeScalar) << " s." << std::endl;
	result << "Total time linear: " << (allTimeLinear) << " s." << std::endl;
	result << "Total time SpMV: " << (allTimeSpMV) << " s." << std::endl;

	n = k;
	//}
	return result.str();
}

std::string printSolveVector(double res, int n, std::vector<double>& x)
{
	std::stringstream result;

	result << "x = [";
	for (int i = 0; i < x.size(); i++)
	{
		if (i == x.size() - 1)
		{
			result << std::fixed << std::showpoint << std::setprecision(5) << x.at(i) << " ]" << std::endl;
		}
		else
		{
			result << std::fixed << std::showpoint << std::setprecision(5) << x.at(i) << ", ";
		}
	}

	result << "Misclosure rate = " << res << std::endl;
	result << "Count of steps = " << n << std::endl;

	return result.str();
}

#pragma endregion

int main(int argc, char* argv[])
{
	setlocale(LC_ALL, "Russian");

	omp_set_num_threads(MaxOmpThreads);
	// Nx = countOfColumns
	// Ny = countOfRows
	// Px = processesColumns
	// Py = processesRows
	int countOfColumns, countOfRows, processesColumns, processesRows, k1, k2, NLocal, NOwn;
	int countOfProcesses, processId;
	// Точность решения
	double tol;
	bool isPrint;
	std::vector<int> arguments, IA, JA, Part, G2L, L2G;
	std::map<int, std::vector<int>> resultGraph;
	// Диапазоны индексов подобласти каждого из процесса
	// ib = columnIndexBegin
	// ie = columnIndexEnd
	// jb = rowIndexBegin
	// je = rowIndexEnd
	int columnIndexBegin, columnIndexEnd, rowIndexBegin, rowIndexEnd;

	std::stringstream debuginfo;

#pragma region Первый этап

#pragma endregion

	if (!initVariables(argc, argv, countOfColumns, countOfRows, k1, k2, processesColumns, processesRows, tol, arguments, isPrint)) {
		return 0;
	}

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &processId);
	MPI_Comm_size(MPI_COMM_WORLD, &countOfProcesses);

	if (!isCountOfProcessesValid(processesColumns, processesRows, countOfProcesses, processId)) {
		return 0;
	}

	calculateSubAreaIndexes(processId, processesColumns, processesRows, columnIndexBegin, columnIndexEnd, rowIndexBegin, rowIndexEnd,
		countOfColumns, countOfRows);

	double start = omp_get_wtime();
	createNodesOfGraph(processId, processesColumns, processesRows, columnIndexBegin, columnIndexEnd, rowIndexBegin, rowIndexEnd,
		countOfColumns, countOfRows, k1, k2, IA, JA, resultGraph, NLocal, NOwn, Part, G2L, L2G);
	double end = omp_get_wtime();
	double seconds = end - start;
	debuginfo << "Process " << processId << ": Total nodes: " << NOwn << " elem." << std::endl;
	debuginfo << "Process " << processId << ": First stage time: " << seconds << " s." << std::endl;
	//printf("Process = %d, columnIndexBegin: %d, columnIndexEnd: %d\, rowIndexBegin: %d, rowIndexEnd: %d\n", processId, columnIndexBegin,
	//	columnIndexEnd, rowIndexBegin, rowIndexEnd);
#pragma endregion

#pragma region Второй этап
	// Вектор ненулевых коэффициентов матрицы
	std::vector<double> A;
	A.resize(IA.at(NOwn));
	// Вектор правой части
	std::vector<double> b;
	b.resize(NOwn);

	double startSecondStage = omp_get_wtime();
	makeSLAE(IA, JA, A, b, NOwn, NLocal, L2G);
	double endSecondStage = omp_get_wtime();
	double endSecondStagePrint = endSecondStage - startSecondStage;
	debuginfo << "Process " << processId << ": Second stage time: " << endSecondStagePrint << " s." << std::endl;
	//std::cout << "Second stage time per 1 elem: " << endSecondStagePrint / NOwn << " s." << std::endl << std::endl;
#pragma endregion

#pragma region Третий этап
	std::vector<int> Neighbours;
	std::vector<int> SendOffet;
	std::vector<int> RecvOffset;
	std::vector<int> Send;
	std::vector<int> Recv;

	SendOffet.push_back(0);
	RecvOffset.push_back(0);

	double startThirdStage = omp_get_wtime();
	createCom(NLocal, NOwn, IA, JA, Part, L2G, G2L, Neighbours, SendOffet, RecvOffset, Send, Recv,
		countOfProcesses);
	double endThirdStage = omp_get_wtime();
	double endThirdStagePrint = endThirdStage - startThirdStage;
	debuginfo << "Process " << processId << ": Third stage time: " << endThirdStagePrint << " s." << std::endl;
#pragma endregion

#pragma region Четвертый этап
	// Вектор решения
	std::vector<double> xRes;
	// Количество итераций
	int n = 0;
	// L2 норма невязки
	double res = 0;

	double startFourthStage = omp_get_wtime();
	debuginfo << solveSLAE(IA, JA, A, b, NOwn, NLocal, tol, xRes, n, res, Neighbours, SendOffet, RecvOffset, Send, Recv, processId) <<
		std::endl;
	double endFourthStage = omp_get_wtime();
	double endFourthStagePrint = endFourthStage - startFourthStage;
	debuginfo << "Process " << processId << ": Fourth stage time: " << endFourthStagePrint << " s." << std::endl;
	//std::cout << "Third stage time per 1 elem (Время третьего этапа на 1 элемент): " << endThirdStagePrint / countOfNodes << " s." <<
	//	std::endl << std::endl;

#pragma endregion

	if (isPrint)
	{
		debuginfo << printResult(resultGraph, IA, JA, seconds, NOwn, processId, G2L, L2G, Part) << std::endl;
		debuginfo << printSLAE(A, b, IA, NOwn, NLocal, processId) << std::endl;
		debuginfo << printCom(Neighbours, SendOffet, RecvOffset, Send, Recv, processId, L2G) << std::endl;
		debuginfo << printSolveVector(res, n, xRes) << std::endl;
	}

	std::cout << debuginfo.str() << std::endl << "     ////////////////////////////////////     " << std::endl << std::endl;

	MPI_Finalize();

	return 0;
}
