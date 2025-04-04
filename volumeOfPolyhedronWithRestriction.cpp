/*
 * MaoHuPi (c) 2025
 * 邊長約束下的凸多面體之體積計算
 *
 * input:
0
10e-5
10
1
4
6
0 1 5
1 2 6
2 0 7
0 3 6
1 3 7
2 3 5

 * output:
19.478
 */

#include <iostream>
#include <vector>
#include <array>

#include <cmath>
#include <cstdlib>
#include <ctime>

using namespace std;

double distance(array<double, 3> p1, array<double, 3> p2)
{
	return sqrt(
		pow(p1.data()[0] - p2.data()[0], 2) +
		pow(p1.data()[1] - p2.data()[1], 2) +
		pow(p1.data()[2] - p2.data()[2], 2));
}
array<double, 3> add(array<double, 3> p1, array<double, 3> p2)
{
	return {
		p1.data()[0] + p2.data()[0],
		p1.data()[1] + p2.data()[1],
		p1.data()[2] + p2.data()[2]};
}
array<double, 3> sub(array<double, 3> p1, array<double, 3> p2)
{
	return {
		p1.data()[0] - p2.data()[0],
		p1.data()[1] - p2.data()[1],
		p1.data()[2] - p2.data()[2]};
}
array<double, 3> mul(double r, array<double, 3> p)
{
	return {
		r * p.data()[0],
		r * p.data()[1],
		r * p.data()[2]};
}
array<double, 3> cross(array<double, 3> p1, array<double, 3> p2)
{
	return {
		p1.data()[1] * p2.data()[2] - p1.data()[2] * p2.data()[1],
		p1.data()[2] * p2.data()[0] - p1.data()[0] * p2.data()[2],
		p1.data()[0] * p2.data()[1] - p1.data()[1] * p2.data()[0]};
}
double dot(array<double, 3> p1, array<double, 3> p2)
{
	return p1.data()[0] * p2.data()[0] + p1.data()[1] * p2.data()[1] + p1.data()[2] * p2.data()[2];
}
int sign(double n)
{
	// return (n > 0) - (n < 0);
	if (n > 0)
		return 1;
	else if (n < 0)
		return -1;
	else
		return 0;
}

void printPoint(array<double, 3> point)
{
	cout << '(' << point[0] << ',' << point[1] << ',' << point[2] << ')' << endl;
}
void printPoints(vector<array<double, 3>> points)
{
	cout << '[';
	for (int i = 0; i < points.size(); i++)
	{
		cout << '(' << points[i][0] << ',' << points[i][1] << ',' << points[i][2] << ')';
		if (i < points.size() - 1)
			cout << ',';
	}
	cout << ']' << endl;
}

vector<array<double, 3>> applyRestriction(vector<array<double, 3>> pointData, vector<array<double, 3>> restrictionData, double errorThreshold)
{
	srand(time(NULL));
	for (int t = 0; t < 500; t++)
	{
		double error = 0;
		for (array<double, 3> restriction : restrictionData)
		{
			array<double, 3> p1 = pointData[(int)restriction[0]],
							 p2 = pointData[(int)restriction[1]],
							 o = {0, 0, 0};
			double currentDistance = distance(p1, p2),
				   targetDistance = restriction[2];

			array<double, 3> vecP1P2 = sub(p2, p1);
			double lenP1P2 = distance(o, vecP1P2),
				   moveLength = (currentDistance - targetDistance) / 2;

			if (lenP1P2 == 0)
			{
				vecP1P2[0] = rand();
				vecP1P2[1] = rand();
				vecP1P2[2] = rand();
				lenP1P2 = distance(o, vecP1P2);
			}
			vecP1P2 = mul(moveLength / lenP1P2, vecP1P2);
			pointData[(int)restriction[0]] = add(p1, vecP1P2);
			pointData[(int)restriction[1]] = sub(p2, vecP1P2);
			error += pow((currentDistance - targetDistance), 2);
		}
		if (error < errorThreshold)
			break;
	}
	return pointData;
}

double getVolume(vector<array<double, 3>> pointData, bool logSamplePoints)
{
	array<double, 3> centerPoint;
	for (array<double, 3> point : pointData)
	{
		centerPoint = add(centerPoint, point);
	}
	centerPoint = mul(1.0 / 3.0, centerPoint);

	for (int i = 0; i < pointData.size(); i++)
	{
		pointData[i] = sub(pointData[i], centerPoint);
	}

	int xMin = 0;
	int xMax = 0;
	int yMin = 0;
	int yMax = 0;
	int zMin = 0;
	int zMax = 0;
	for (array<double, 3> point : pointData)
	{
		if (point[0] < xMin)
			xMin = point[0];
		if (point[0] > xMax)
			xMax = point[0];
		if (point[1] < yMin)
			yMin = point[1];
		if (point[1] > yMax)
			yMax = point[1];
		if (point[2] < zMin)
			zMin = point[2];
		if (point[2] > zMax)
			zMax = point[2];
	}
	xMin -= 1;
	xMax += 1;
	yMin -= 1;
	yMax += 1;
	zMin -= 1;
	zMax += 1;
	double volumeSum = 0;
	vector<array<double, 3>> sample;
	for (int x = xMin; x <= xMax; x++)
		for (int y = yMin; y <= yMax; y++)
			for (int z = zMin; z <= zMax; z++)
			{
				// x = 100; y = 100; z = 100; //
				array<double, 3> samplePoint = {(double)x, (double)y, (double)z};
				vector<array<double, 3>> vecToPoints(pointData.size());
				for (int i = 0; i < pointData.size(); i++)
				{
					vecToPoints[i] = sub(pointData[i], samplePoint);
				}
				int sameDirectionForVecDataCount = 0;
				int compositionCount = 0;
				for (int p1 = 0; p1 < vecToPoints.size(); p1++)
					for (int p2 = p1; p2 < vecToPoints.size(); p2++)
						for (int p3 = p2; p3 < vecToPoints.size(); p3++)
						{
							if (p1 == p2 || p2 == p3 || p3 == p1)
								continue;
							compositionCount++;
							array<array<double, 3>, 3> vecData = {vecToPoints[p1], vecToPoints[p2], vecToPoints[p3]};
							array<double, 3> vecSum = cross(sub(vecData[1], vecData[0]), sub(vecData[2], vecData[0]));
							int sameDirectionCount = 0;
							for (int i = 0; i < vecToPoints.size(); i++)
							{
								sameDirectionCount += sign(dot(vecToPoints[i], vecSum));
							}
							sameDirectionCount = abs(sameDirectionCount);
							if (sameDirectionCount < pointData.size())
							{
								sameDirectionForVecDataCount++;
							}
						}
				if (sameDirectionForVecDataCount == compositionCount)
				{
					volumeSum++;
					if (logSamplePoints)
						sample.push_back(samplePoint);
				}
			}
	if (logSamplePoints)
	{
		printPoints(pointData);
		printPoints(sample);
	}
	return volumeSum;
}

int main()
{
	const bool customizeMode = true;

	// 是否印出擬合點與取樣點資料
	bool logSamplePoints = false;
	// 擬合約束條件時，結束擬合所採用的最大誤差值
	double errorThreshold = 10e-5;
	// 採樣所用之網格細分數（以二次伸縮進行實作）
	double resolution = 15;
	// 取樣幾次後取平均
	int sampleTimes = 1;
	if (customizeMode)
	{
		cout << "logSamplePoints(0 or 1, 0): ";
		cin >> logSamplePoints;
		cout << "errorThreshold(double, 10e-5): ";
		cin >> errorThreshold;
		cout << "resolution(double, 15): ";
		cin >> resolution;
		cout << "sampleTimes(int, 1): ";
		cin >> sampleTimes;
	}

	// 節點數目
	int pointNum = 4;
	if (customizeMode)
	{
		cout << "pointNum(int, 4): ";
		cin >> pointNum;
	}
	// 約束條件
	vector<array<double, 3>> restrictionData;
	if (customizeMode)
	{
		int restrictionNum;
		cout << "restrictionData(lineCount and data with [point 1 index(int)] [point 2 index(int)] [target distance(double)]): " << endl;
		cin >> restrictionNum;
		int p1Index, p2Index;
		double targetDistance;
		array<double, 3> restriction;
		for (int i = 0; i < restrictionNum; i++)
		{
			cin >> p1Index >> p2Index >> targetDistance;
			restriction = {(double)p1Index, (double)p2Index, targetDistance * resolution};
			restrictionData.push_back(restriction);
		}
	} else {
		double a = 5, b = 6, c = 7;
		a *= resolution;
		b *= resolution;
		c *= resolution;
		restrictionData = {
			{0, 1, a}, {1, 2, b}, {2, 0, c}, {0, 3, b}, {1, 3, c}, {2, 3, a}};
	}

	// 取樣迴圈
	double volumeAv = 0;
	for (int i = 0; i < sampleTimes; i++)
	{
		vector<std::array<double, 3>> pointData(pointNum);
		for (int j = 0; j < pointNum; j++)
		{
			pointData[j] = {0, 0, 0};
		}
		pointData = applyRestriction(pointData, restrictionData, errorThreshold);
		double a = getVolume(pointData, logSamplePoints) / pow(resolution, 3);
		volumeAv += a;
	}
	volumeAv /= sampleTimes;

	// 印出體積平均值
	cout << volumeAv << endl;
	return 0;
}