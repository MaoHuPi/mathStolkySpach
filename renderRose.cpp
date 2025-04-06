/*
 * 2025 (c) MaoHuPi
 * renderRose.cpp
 * global x: right, y: front, z: up
 * the camera is rotated 180 degree by default to match the right side of rose surface
 * 
 * 輸出紀錄：1920 * 1080 * 400張 (.ppm) => 4.08GB
 */

/* include */
#include <iostream>
#include <string>
#include <cmath>
#include <functional>
#include <array>
#include <map>
#include <tuple>
#include <sstream>

/* namespace */
using namespace std;

/* define and constant */
#define PI 3.14159265358979323846
// const double OUTPUT_SCALE = 0.1;
const double OUTPUT_SCALE = 1;
const int SCREEN_WIDTH = (int)(1920 * OUTPUT_SCALE);
const int SCREEN_HEIGHT = (int)(1080 * OUTPUT_SCALE);
const double SCREEN_DEPTH = (int)(960 * OUTPUT_SCALE);
const double DELTA_UNIT = 0.00001;
const bool USE_LIGHT = true;
const bool USE_FOG = true;
const double FOG_RADIUS = 100;
const int OUTPUT_FRAMES = 400;
const bool BOXES_ACCELERATION = false;

double BALL_Y = -12;

/* operation */
// scalar
double sign(double n)
{
	return (n > 0) - (n < 0);
}
// vector
void printVector(array<double, 3> vec)
{
	double *v = vec.data();
	cout << '(' << v[0] << ',' << v[1] << ',' << v[2] << ')' << endl;
}
double vectorDistance(array<double, 3> A, array<double, 3> B)
{
	double *a = A.data();
	double *b = B.data();
	return sqrt(pow(a[0] - b[0], 2) + pow(a[1] - b[1], 2) + pow(a[2] - b[2], 2));
}
double vectorLength(array<double, 3> vec)
{
	double *v = vec.data();
	return sqrt(pow(v[0], 2) + pow(v[1], 2) + pow(v[2], 2));
}
array<double, 3> vectorCopy(array<double, 3> vec)
{
	double *v = vec.data();
	return {v[0], v[1], v[2]};
}
array<double, 3> vectorAddition(array<double, 3> A, array<double, 3> B)
{
	double *a = A.data();
	double *b = B.data();
	return {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
}
array<double, 3> vectorScale(double r, array<double, 3> vec)
{
	double *v = vec.data();
	return {r * v[0], r * v[1], r * v[2]};
}
double vectorDot(array<double, 3> A, array<double, 3> B)
{
	double *a = A.data();
	double *b = B.data();
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
array<double, 3> vectorCross(array<double, 3> A, array<double, 3> B)
{
	double *a = A.data();
	double *b = B.data();
	return {
		a[1] * b[2] - a[2] * b[1],
		a[2] * b[0] - a[0] * b[2],
		a[0] * b[1] - a[1] * b[0]};
}
array<double, 4> vector2quaternion(array<double, 3> vec)
{
	double *v = vec.data();
	return {0, v[0], v[1], v[2]};
}
array<int, 3> vectorD2I(array<double, 3> vec)
{
	double *v = vec.data();
	return {(int)v[0], (int)v[1], (int)v[2]};
}
array<double, 3> vectorI2D(array<int, 3> vec)
{
	int *v = vec.data();
	return {(double)v[0], (double)v[1], (double)v[2]};
}
// quaternion
array<double, 3> quaternion2vector(array<double, 4> qua)
{
	double *q = qua.data();
	return {q[1], q[2], q[3]};
}
array<double, 4> quaternionConjugate(array<double, 4> qua)
{
	double *q = qua.data();
	return {q[0], -q[1], -q[2], -q[3]};
}
array<double, 4> quaternionMultiply(array<double, 4> A, array<double, 4> B)
{
	double *a = A.data();
	double *b = B.data();
	return {
		a[0] * b[0] - a[1] * b[1] - a[2] * b[2] - a[3] * b[3],
		a[0] * b[1] + a[1] * b[0] + a[2] * b[3] - a[3] * b[2],
		a[0] * b[2] + a[2] * b[0] - a[1] * b[3] + a[3] * b[1],
		a[0] * b[3] + a[3] * b[0] + a[1] * b[2] - a[2] * b[1]};
}
array<double, 4> rotationQuaternion(array<double, 3> axis, double angle)
{
	double *a = axis.data();
	double sineHalfAngle = sin(angle / 2.0);
	return {
		cos(angle / 2.0),
		a[0] * sineHalfAngle,
		a[1] * sineHalfAngle,
		a[2] * sineHalfAngle};
}
array<double, 4> mat3_normalized_to_quat_fast(array<array<double, 3>, 3> param_mat)
{
	// https://github.com/blender/blender/blob/756538b4a117cb51a15e848fa6170143b6aafcd8/source/blender/blenlib/intern/math_rotation.c#L272
	/* Method outlined by Mike Day, ref: https://math.stackexchange.com/a/3183435/220949
	 * with an additional `sqrtf(..)` for higher precision result.
	 * Removing the `sqrt` causes tests to fail unless the precision is set to 1e-6 or larger. */
	double q[4] = {0, 0, 0, 0};
	double mat[3][3];
	mat[0][0] = param_mat[0].data()[0];
	mat[0][1] = param_mat[0].data()[1];
	mat[0][2] = param_mat[0].data()[2];
	mat[1][0] = param_mat[1].data()[0];
	mat[1][1] = param_mat[1].data()[1];
	mat[1][2] = param_mat[1].data()[2];
	mat[2][0] = param_mat[2].data()[0];
	mat[2][1] = param_mat[2].data()[1];
	mat[2][2] = param_mat[2].data()[2];
	if (mat[2][2] < 0.0)
	{
		if (mat[0][0] > mat[1][1])
		{
			const double trace = 1.0 + mat[0][0] - mat[1][1] - mat[2][2];
			double s = 2.0 * sqrt(trace);
			if (mat[1][2] < mat[2][1])
			{
				s = -s;
			}
			q[1] = 0.25 * s;
			s = 1.0 / s;
			q[0] = (mat[1][2] - mat[2][1]) * s;
			q[2] = (mat[0][1] + mat[1][0]) * s;
			q[3] = (mat[2][0] + mat[0][2]) * s;
			if ((trace == 1.0) && (q[0] == 0.0 && q[2] == 0.0 && q[3] == 0.0))
			{
				q[1] = 1.0;
			}
		}
		else
		{
			const double trace = 1.0 - mat[0][0] + mat[1][1] - mat[2][2];
			double s = 2.0 * sqrt(trace);
			if (mat[2][0] < mat[0][2])
			{
				s = -s;
			}
			q[2] = 0.25 * s;
			s = 1.0 / s;
			q[0] = (mat[2][0] - mat[0][2]) * s;
			q[1] = (mat[0][1] + mat[1][0]) * s;
			q[3] = (mat[1][2] + mat[2][1]) * s;
			if ((trace == 1.0) && (q[0] == 0.0 && q[1] == 0.0 && q[3] == 0.0))
			{
				q[2] = 1.0;
			}
		}
	}
	else
	{
		if (mat[0][0] < -mat[1][1])
		{
			const double trace = 1.0 - mat[0][0] - mat[1][1] + mat[2][2];
			double s = 2.0 * sqrt(trace);
			if (mat[0][1] < mat[1][0])
			{
				s = -s;
			}
			q[3] = 0.25 * s;
			s = 1.0 / s;
			q[0] = (mat[0][1] - mat[1][0]) * s;
			q[1] = (mat[2][0] + mat[0][2]) * s;
			q[2] = (mat[1][2] + mat[2][1]) * s;
			if ((trace == 1.0) && (q[0] == 0.0 && q[1] == 0.0 && q[2] == 0.0))
			{
				q[3] = 1.0;
			}
		}
		else
		{
			/* NOTE(@campbellbarton): A zero matrix will fall through to this block,
			 * needed so a zero scaled matrices to return a quaternion without rotation, see: T101848. */
			const double trace = 1.0 + mat[0][0] + mat[1][1] + mat[2][2];
			double s = 2.0 * sqrt(trace);
			q[0] = 0.25 * s;
			s = 1.0 / s;
			q[1] = (mat[1][2] - mat[2][1]) * s;
			q[2] = (mat[2][0] - mat[0][2]) * s;
			q[3] = (mat[0][1] - mat[1][0]) * s;
			if ((trace == 1.0) && (q[1] == 0.0 && q[2] == 0.0 && q[3] == 0.0))
			{
				q[0] = 1.0;
			}
		}
	}

	return {q[0], q[1], q[2], q[3]};
}

/* object */
// implicit object
class implicitObject
{
public:
	string id;
	function<double(array<double, 3>)> implicit;
	function<array<int, 3>(array<double, 3>)> render;
	vector<array<array<double, 3>, 8>> updateBoxes;
	/*
	 * 7|\--------/|1
	 *  | 5\----/3 |
	 *  |  |    |  |
	 *  | 6/----\4 |
	 * 8|/--------\|2
	 */
	vector<array<array<double, 2>, 8>> updateBoxesPoly;
	vector<array<double, 4>> updateBoxesRect;
	bool updateFlag;
	bool updateBoxesCheck;
	implicitObject()
		: id(""), implicit(nullptr), render(nullptr), updateBoxes({}), updateBoxesPoly({}), updateBoxesRect({}), updateFlag(true), updateBoxesCheck(true) {}
	implicitObject(string param_id, function<double(array<double, 3>)> param_implicit, function<array<int, 3>(array<double, 3>)> param_render, vector<array<array<double, 3>, 8>> param_updateBoxes)
	{
		id = param_id;
		implicit = param_implicit;
		render = param_render;
		updateBoxes = param_updateBoxes;
		updateBoxesPoly = {};
		updateBoxesRect = {};
		updateFlag = true;
		updateBoxesCheck = BOXES_ACCELERATION ? param_updateBoxes.size() > 0 : false;
	}
};

/* specific function */
// array<double, 3> cameraVector2globalVector(array<double, 3> pos, array<double, 4> rot, array<double, 3> vec)
// {
// 	return vectorAddition(
// 		pos,
// 		quaternion2vector(
// 			quaternionMultiply(
// 				quaternionMultiply(rot, vector2quaternion(vec)),
// 				quaternionConjugate(rot))));
// }
array<double, 2> rose_global2plain(array<double, 3> pos)
{
	double *p = pos.data();
	double l = vectorLength(pos);
	double tanL = tan(l);
	double cotL = 1 / tanL;
	double tanLSquared = pow(tanL, 2);
	double cotLSquared = pow(cotL, 2);
	double diractionLength = l * sqrt((tanLSquared + cotLSquared) * (tanLSquared + cotLSquared + 1));
	return {
		l,
		atan2(sqrt(tanLSquared + cotLSquared + 1) * (p[1] * cotL - p[2] * tanL), p[0] * (tanLSquared + cotLSquared) - p[1] * tanL - p[2] * cotL)};
}
array<double, 3> rose_plain2global(array<double, 2> point)
{
	double *p = point.data();
	double l = p[0], theta = p[1];
	double tanLSquared = pow(tan(l), 2);
	double cotLSquared = 1 / tanLSquared;
	double c1 = l * sqrt((tanLSquared + cotLSquared) * (tanLSquared + cotLSquared + 1));
	double c2 = sqrt(tanLSquared + cotLSquared + 1);
	double c3 = cos(theta) / (tanLSquared + cotLSquared + 1);
	return {c1 * c3,
			c1 / (tan(l) + pow(1 / tan(l), 3)) * (cotLSquared * sin(theta) / c2 - c3),
			c1 / (pow(tan(l), 3) + 1 / tan(l)) * (-tanLSquared * sin(theta) / c2 - c3)};
}
const array<int, 3> empty_color()
{
	// return {135, 206, 235};
	return {0, 0, 0};
}
array<double, 3> light;
array<int, 3> renderByNormal(array<double, 3> pos, array<double, 3> normal, array<double, 3> color)
{
	if (USE_LIGHT)
	{
		normal = vectorScale(1 / vectorLength(normal), normal);
		if (vectorDot(pos, normal) > 0)
			normal = vectorScale(-1, normal);
		double lightness = max(min((-vectorDot(light, normal) + 1) / 2.0, 1.0), 0.0);
		// double lightness = max(min(-vectorDot(light, normal), 1.0), 0.0);
		color = vectorScale(lightness, color);
	}
	if (USE_FOG)
	{
		double fogginess = max(min(vectorLength(pos) / FOG_RADIUS, 1.0), 0.0);
		color = vectorAddition(
			vectorScale(1.0 - fogginess, color),
			vectorScale(fogginess, vectorI2D(empty_color())));
	}
	return vectorD2I(color);
}

/* program entrance */
int main()
{
	/* define implicit objects */
	// floor
	implicitObject obj_floor("floor", [](array<double, 3> pos) -> double
							 { return pos.data()[2] - (-10); }, [](array<double, 3> pos) -> array<int, 3>
							 { return renderByNormal(pos, {0, 0, 1}, {255, 255, 255}); }, {{{{100, -100, 0.001}, {100, -100, -0.001}, {100, -0.001, 0.001}, {100, -0.001, -0.001}, {-100, -0.001, 0.001}, {-100, -0.001, -0.001}, {-100, -100, 0.001}, {-100, -100, 0.001}}}});
	// ball
	function<double(array<double, 3>)> ball_implicit = [](array<double, 3> pos) -> double
	{
		double *p = pos.data();
		double oriRes = pow(p[0], 2) + pow(p[1] - BALL_Y, 2) + pow(p[2], 2) - pow(5, 2);
		// return oriRes;
		return sign(oriRes) * sqrt(abs(oriRes));
	};
	function<array<int, 3>(array<double, 3>)> ball_render = [](array<double, 3> pos) -> array<int, 3>
	{
		array<double, 3> nVec = vectorAddition(pos, {0, -BALL_Y, 0});
		nVec = vectorScale(1 / vectorLength(nVec), nVec);
		return renderByNormal(pos, nVec, {255, 255, 255});
	};
	implicitObject obj_ball("ball", ball_implicit, ball_render, {{{{5, BALL_Y - 5, 5}, {5, BALL_Y - 5, -5}, {5, BALL_Y + 5, 5}, {5, BALL_Y + 5, -5}, {-5, BALL_Y + 5, 5}, {-5, BALL_Y + 5, -5}, {-5, BALL_Y - 5, 5}, {-5, BALL_Y - 5, -5}}}});
	// rose
	function<double(array<double, 3>)> rose_implicit = [](array<double, 3> pos) -> double
	{
		array<double, 3> roseCenter = {0, 0, 0};
		pos = vectorAddition(pos, vectorScale(-1, roseCenter));
		double tanL = tan(vectorLength(pos) * 1.0);
		double *p = pos.data();
		double oriRes = p[0] + tanL * p[1] + (1 / tanL) * p[2];
		// return tanh(0.5 * oriRes);
		return 0.05 * oriRes / (1 + abs(oriRes));
	};
	function<array<int, 3>(array<double, 3>)> rose_render = [](array<double, 3> pos) -> array<int, 3>
	{
		double *p = pos.data();
		array<double, 2> pos2 = rose_global2plain(pos);
		double *p2 = pos2.data();
		p2[0] = abs(p2[0]);
		while (p2[1] < 0)
			p2[1] += 2 * PI;
		array<double, 3> rePos = rose_plain2global({p2[0], p2[1]});
		array<double, 3> nVec = vectorCross(
			vectorAddition(rose_plain2global({p2[0] + DELTA_UNIT, p2[1]}), vectorScale(-1, rePos)),
			vectorAddition(vectorScale(((int)floor(p2[0] / (PI / 2)) % 2 == 0 ? 1 : -1), rose_plain2global({p2[0], p2[1] + DELTA_UNIT})), vectorScale(-1, rePos)));
		//
		array<double, 3> color = (vector<array<double, 3>>{
			{255, 0, 0},
			{255, 0, 0},
			{255, 0, 0},
			{255, 255, 255}})[(((int)(sqrt(p2[0]) / PI * 4000) % 8 <= 6) * 2 + ((int)(p2[1] / PI * 720) % 8 <= 6))];
		return renderByNormal(pos, nVec, color);
	};
	implicitObject obj_rose("rose", rose_implicit, rose_render, {});

	/* scene setup */
	light = {1, -1, -1};
	light = vectorScale(1 / vectorLength(light), light);
	vector<implicitObject> implicitObjectList = {obj_rose};
	// vector<implicitObject> implicitObjectList = {obj_floor, obj_ball};
	map<string, implicitObject> implicitObjectMap;
	for (implicitObject obj : implicitObjectList)
	{
		implicitObjectMap[obj.id] = obj;
	}
	array<double, 3> cameraPos = {0, 0, 0.01};
	array<double, 4> cameraRot = rotationQuaternion({0, 0, 1}, PI);

	/* image render loop */
	for (int nth = 0; nth < OUTPUT_FRAMES; nth++)
	{
		double r = 0.001 + nth * (2 * PI) / OUTPUT_FRAMES;
		double theta = PI;
		cameraPos = rose_plain2global({r, theta});
		array<double, 3> tanR = vectorAddition(rose_plain2global({r + DELTA_UNIT, theta}), vectorScale(-1, cameraPos));
		tanR = vectorScale(1 / vectorLength(tanR), tanR);
		array<double, 3> tanTheta = vectorAddition(rose_plain2global({r, theta + DELTA_UNIT}), vectorScale(-1, cameraPos));
		tanTheta = vectorScale(1 / vectorLength(tanTheta) * ((int)floor(r / (PI / 2)) % 2 == 0 ? 1 : -1), tanTheta);
		array<double, 3> normal = vectorCross(tanR, tanTheta);
		normal = vectorScale(1 / vectorLength(tanTheta), normal); // 因為徑向、角度方向的兩切向量不一定垂直，故該兩單位長度向量之外積結果不一定為單位向量
		tanTheta = vectorCross(normal, tanR);					  // 因為不是正交坐標系，故角度方向的切向量需做修正，方能保證正交
		cameraRot = mat3_normalized_to_quat_fast({{{-tanTheta.data()[0], tanR[0], normal[0]}, {-tanTheta.data()[1], tanR[1], normal[1]}, {-tanTheta.data()[2], tanR[2], normal[2]}}});
		cameraRot = quaternionConjugate(cameraRot); // 旋轉四元數似乎是反向的？待理解
		light = quaternion2vector(quaternionMultiply(quaternionMultiply(cameraRot, vector2quaternion({0, 1, 0})), quaternionConjugate(cameraRot)));
		cameraPos = vectorAddition(cameraPos, vectorScale(0.3, normal));
		// cameraRot = quaternionMultiply(cameraRot, rotationQuaternion({1, 0, 0}, -PI / 2)); // 測試是否能正確地「看地板」，且地面向畫面下方移動
		// printVector(cameraPos);

		array<double, 4> cameraRotInverse = quaternionConjugate(cameraRot);
		if (BOXES_ACCELERATION)
		{
			for (implicitObject &obj : implicitObjectList)
			{
				vector<array<array<double, 2>, 8>> updateBoxesPoly = {};
				vector<array<double, 4>> updateBoxesRect = {};
				if (obj.updateBoxes.size() > 0)
				{
					int rectOutOfScreenCount = 0;
					for (array<array<double, 3>, 8> box : obj.updateBoxes)
					{
						array<array<double, 2>, 8> boxPoly = {};
						for (int i = 0; i < 8; i++)
						{
							array<double, 3> vertex = box[i];
							double *relativeVertexPosition = quaternion2vector(quaternionMultiply(quaternionMultiply(cameraRotInverse, vector2quaternion(vectorAddition(vertex, vectorScale(-1, cameraPos)))), cameraRot)).data();
							// 使用方錐投影，將 relativeVertexPosition 投影至 screen 上，並加入 boxPoly。
							// t * relativeVertexPosition = {x, cameraDepth, z}
							// => t * relativeVertexPosition[1] = cameraDepth
							// => t = cameraDepth / relativeVertexPosition[1]
							// => x = relativeVertexPosition[0] * t
							// => z = relativeVertexPosition[2] * t
							// => i = cameraHeight / 2 - relativeVertexPosition[2] * t
							// => j = cameraWidth / 2 + relativeVertexPosition[0] * t
							double t = SCREEN_DEPTH / relativeVertexPosition[1];
							boxPoly[i] = array<double, 2>{SCREEN_HEIGHT / 2 - relativeVertexPosition[2] * t,
														  SCREEN_WIDTH / 2 + relativeVertexPosition[0] * t};
						}
						updateBoxesPoly.push_back(boxPoly);
						double iMax = 0, iMin = SCREEN_HEIGHT, jMax = 0, jMin = SCREEN_WIDTH;
						for (array<double, 2> screenVertex : boxPoly)
						{
							double *v = screenVertex.data();
							if (v[0] > iMax)
								iMax = v[0];
							if (v[0] < iMin)
								iMin = v[0];
							if (v[1] > jMax)
								jMax = v[1];
							if (v[1] < jMin)
								jMin = v[1];
						}
						if (iMax == 0 || iMin == SCREEN_HEIGHT || jMax == 0 || jMin == SCREEN_WIDTH)
							rectOutOfScreenCount++;
						// updateBoxesRect.push_back({iMin, jMin, iMax, jMax});
						else
							updateBoxesRect.push_back({iMin, jMin, iMax, jMax});
					}
					obj.updateFlag = rectOutOfScreenCount < obj.updateBoxes.size();
					obj.updateBoxesCheck = obj.updateFlag;
				}
				else
				{
					obj.updateFlag = true;
					obj.updateBoxesCheck = false;
					updateBoxesPoly = {};
					updateBoxesRect = {};
				}
				// cout << '['
				// 	 << '(' << updateBoxesPoly[0][0].data()[0] << ',' << updateBoxesPoly[0][0].data()[1] << ')' << ','
				// 	 << '(' << updateBoxesPoly[0][1].data()[0] << ',' << updateBoxesPoly[0][1].data()[1] << ')' << ','
				// 	 << '(' << updateBoxesPoly[0][2].data()[0] << ',' << updateBoxesPoly[0][2].data()[1] << ')' << ','
				// 	 << '(' << updateBoxesPoly[0][3].data()[0] << ',' << updateBoxesPoly[0][3].data()[1] << ')' << ','
				// 	 << '(' << updateBoxesPoly[0][4].data()[0] << ',' << updateBoxesPoly[0][4].data()[1] << ')' << ','
				// 	 << '(' << updateBoxesPoly[0][5].data()[0] << ',' << updateBoxesPoly[0][5].data()[1] << ')' << ','
				// 	 << '(' << updateBoxesPoly[0][6].data()[0] << ',' << updateBoxesPoly[0][6].data()[1] << ')' << ','
				// 	 << '(' << updateBoxesPoly[0][7].data()[0] << ',' << updateBoxesPoly[0][7].data()[1] << ')' << ']' << endl;
				// cout << "x=" << updateBoxesRect[0].data()[0] << '\n'
				// 	 << "y=" << updateBoxesRect[0].data()[1] << '\n'
				// 	 << "x=" << updateBoxesRect[0].data()[2] << '\n'
				// 	 << "y=" << updateBoxesRect[0].data()[3] << endl;
				obj.updateBoxesPoly = updateBoxesPoly;
				obj.updateBoxesRect = updateBoxesRect;
			}
		}

		vector<vector<array<int, 3>>> renderFrame(SCREEN_HEIGHT, vector<array<int, 3>>(SCREEN_WIDTH, {0, 0, 0}));
		for (int i = 0; i < SCREEN_HEIGHT; i++)
			for (int j = 0; j < SCREEN_WIDTH; j++)
			{
				array<double, 3> cameraVec = {
					(double)j - (double)SCREEN_WIDTH / 2.0,
					SCREEN_DEPTH,
					(double)SCREEN_HEIGHT / 2.0 - (double)i};
				array<double, 3> diractionVec = quaternion2vector(quaternionMultiply(
					quaternionMultiply(cameraRot, vector2quaternion(cameraVec)),
					quaternionConjugate(cameraRot)));
				diractionVec = vectorScale(1 / vectorLength(diractionVec), diractionVec);
				array<double, 3> globalVec = cameraPos;
				double minDistance = INFINITY;
				string targetId = "";
				for (implicitObject obj : implicitObjectList)
				{
					if (obj.updateFlag == false)
						continue;
					if (BOXES_ACCELERATION && obj.updateBoxesCheck)
					{
						bool skipFlag = true;
						for (array<double, 4> rect : obj.updateBoxesRect)
						{
							double *r = rect.data();
							if (r[0] < i && r[1] < j && r[2] > i && r[3] > j)
							{
								skipFlag = false;
								break;
							}
						}
						if (skipFlag)
							continue;

						skipFlag = false;
						for (array<array<double, 2>, 8> poly : obj.updateBoxesPoly)
						{
							vector<array<array<double, 2>, 2>> edges = {
								{poly[0], poly[1]},
								{poly[2], poly[3]},
								{poly[4], poly[5]},
								{poly[6], poly[7]},
								{poly[0], poly[6]},
								{poly[1], poly[7]},
								{poly[2], poly[4]},
								{poly[3], poly[5]},
								{poly[0], poly[2]},
								{poly[1], poly[3]},
								{poly[6], poly[4]},
								{poly[7], poly[5]}};
							for (array<array<double, 2>, 2> edge : edges)
							{
								double *p1 = edge[0].data();
								double *p2 = edge[1].data();
								double edgeNormal[2] = {p2[1] - p1[1], -(p2[0] - p1[0])};
								int sameSide;
								for (array<double, 2> vertex : poly)
								{
									double *v = vertex.data();
									sameSide += (int)sign(edgeNormal[0] * (v[0] - i) + edgeNormal[1] * (v[1] - j));
								}
								if (abs(sameSide) == 8)
								{
									skipFlag = true;
									break;
								}
							}
							if (skipFlag)
								break;
						}
						if (skipFlag)
							continue;
					}
					// TODO: 利用 updateBoxesPoly 判斷此像素 (i, j) 是否需要測距，若否則跳過
					array<double, 3> tracingVec = vectorCopy(globalVec);
					const double resSign = sign(obj.implicit(tracingVec));
					// for (int k = 0; k < 100; k++)
					for (int k = 0; k < 400; k++)
					{
						double res = obj.implicit(tracingVec) / resSign;
						tracingVec = vectorAddition(tracingVec, vectorScale(res, diractionVec));
						if (isnan(res) || isinf(res))
							break;
						// else if (res < 0.01)
						else if (res < 0.001)
						{
							double objDistance = vectorLength(vectorAddition(tracingVec, vectorScale(-1, globalVec)));
							if (objDistance < minDistance)
							{
								minDistance = objDistance;
								targetId = obj.id;
							}
							break;
						}
					}
				}
				array<int, 3> color;
				if (isnan(minDistance) || isinf(minDistance))
				{
					color = empty_color();
				}
				else
				{
					/* distance */
					// int gray = min(max((int)(minDistance), 0), 255);
					// color = {gray, gray, gray};

					/* object color, lighting, and fog effect */
					color = implicitObjectMap[targetId].render(vectorAddition(vectorScale(minDistance, diractionVec), globalVec));
				}
				int *c = color.data();
				renderFrame[i][j] = color;
			}

		/* save image as ppm format */
		ostringstream colorData;
		for (size_t i = 0; i < SCREEN_HEIGHT; ++i)
			for (size_t j = 0; j < SCREEN_WIDTH; ++j)
				for (size_t k = 0; k < 3; ++k)
				{
					colorData << renderFrame[i][j][k] << ' ';
				}
		string imageData = "P3\n" + to_string(SCREEN_WIDTH) + " " + to_string(SCREEN_HEIGHT) + "\n255\n" + colorData.str();
		FILE *fp = fopen(("out/renderRose_" + to_string(nth) + ".ppm").c_str(), "w");
		fwrite(imageData.c_str(), 1, imageData.size(), fp);
		fclose(fp);
		cout << nth * 100 / OUTPUT_FRAMES << '%' << endl;
		// cameraPos = vectorAddition(cameraPos, {0, 0, 0.05});
		// BALL_Y += 0.05;
	}

	return 0;
}