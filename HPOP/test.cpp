//------------------------------------------------------------------------------
//
// High Precision Orbit Propagator
//
//
// Last modified:
//
//   2000/03/04  OMO  Final version (1st edition)
//   2005/04/14  OMO  Final version (2nd reprint)
//
// (c) 1999-2008  O. Montenbruck, E. Gill and Meysam Mahooti
//
//------------------------------------------------------------------------------
//#include <stdc++.h>
//#include <conio.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <time.h>
#include <vector>
#include "GNU_iomanip.h"

#include "Polygon_Relation.h"
#include "SAT_Const.h"
#include "SAT_DE.h"
#include "SAT_Force.h"
#include "SAT_RefSys.h"
#include "SAT_Time.h"
#include "SAT_VecMat.h"
#include "APC_Planets.h"
#include "eopspw.h"

using namespace std;

//------------------------------------------------------------------------------
//
// Global types and data
//
//------------------------------------------------------------------------------
Matrix cnm(361, 361), snm(361, 361);
const double R_ref = 6378.1363e3;   // Earth's radius [m]; GGM03S
const double GM_ref = 398600.4415e9; // [m^3/s^2]; GGM03S
spwdata spwarr[spwsize];
eopdata eoparr[eopsize];
int dat;
char interp = 'l';
double jd, jdeopstart, mfme, dut1, lod, xp, yp, ddpsi, ddeps, dx, dy, x, y, s,
deltapsi, deltaeps;
double jdspwstart, f107a, f107, f107bar, ap, avgap, kp, sumkp, aparr[8], kparr[8];

// Record for passing global data between Deriv and the calling program
struct AuxParam {
	double  Mjd0_TT;
	double  Area_drag, Area_solar, mass, CR, CD;
	int     n, m;
	bool    Sun, Moon, SRad, Drag;
};

//------------------------------------------------------------------------------
//
// Accel
//
// Purpose:
//
//   Computes the acceleration of an Earth orbiting satellite due to
//    - the Earth's harmonic gravity field,
//    - the gravitational perturbations of the Sun and Moon
//    - the solar radiation pressure and
//    - the atmospheric drag
//
// Input/Output:
//
//   Mjd_TT      Terrestrial Time (Modified Julian Date)
//   r           Satellite position vector in the ICRF/EME2000 system
//   v           Satellite velocity vector in the ICRF/EME2000 system
//   Area_drag   Cross-section
//   Area_solar  Cross-section
//   mass        Spacecraft mass
//   CR          Radiation pressure coefficient
//   CD          Drag coefficient
//   <return>    Acceleration (a=d^2r/dt^2) in the ICRF/EME2000 system
//
//------------------------------------------------------------------------------
Vector Accel(double Mjd_TT, const Vector& r, const Vector& v, double Area_drag,
	double Area_solar, double mass, double CR, double CD, int n, int m,
	bool FlagSun, bool FlagMoon, bool FlagSRad, bool FlagDrag)
{
	double Mjd_UTC = 0.0;
	double Mjd_UT1;
	double T1;     // Julian cent. since J2000
	Vector a(3), r_Sun(3), r_Moon(3);
	Matrix P(3, 3), N(3, 3), T(3, 3), E(3, 3);
	char fluxtype, f81type, inputtype;

	// Acceleration due to harmonic gravity field
	Mjd_UTC = Mjd_TT - IERS::TT_UTC(Mjd_UTC) / 86400.0;
	Mjd_UT1 = Mjd_UTC + IERS::UT1_UTC(Mjd_UTC) / 86400.0;

	jd = Mjd_UTC + 2400000.5;
	mfme = 1440.0*(Mjd_UTC - floor(Mjd_UTC));

	findeopparam(jd, mfme, interp, eoparr, jdeopstart, dut1, dat, lod, xp, yp,
		ddpsi, ddeps, dx, dy, x, y, s, deltapsi, deltaeps);

	fluxtype = 'o';
	f81type = 'c';
	inputtype = 'a';
	findatmosparam(jd, mfme, interp, fluxtype, f81type, inputtype, spwarr,
		jdspwstart, f107, f107bar, ap, avgap, aparr, kp, sumkp, kparr);

	IERS::Set(dut1, -dat, xp, yp);

	P = PrecMatrix(MJD_J2000, Mjd_TT);
	N = NutMatrix(Mjd_TT);
	T = N*P;
	E = GHAMatrix(Mjd_UT1) * T;
	// E = PoleMatrix(Mjd_UTC) * GHAMatrix(Mjd_UT1) * T;

	a = AccelHarmonic(r, E, GM_ref, R_ref, cnm, snm, n, m);

	// Luni-solar perturbations
	//  r_Sun  = Sun(Mjd_TT);
	//  r_Moon = Moon(Mjd_TT);
	T1 = (Mjd_TT - MJD_J2000) / 36525.0;
	r_Sun = AU*Transp(EclMatrix(Mjd_TT)*P)*SunPos(T1);
	r_Moon = Transp(EclMatrix(Mjd_TT)*P)*MoonPos(T1);

	if (FlagSun)  a += AccelPointMass(r, r_Sun, GM_Sun);
	if (FlagMoon) a += AccelPointMass(r, r_Moon, GM_Moon);

	// Solar radiation pressure
	if (FlagSRad) a += AccelSolrad(r, r_Sun, Area_solar, mass, CR, P_Sol, AU);

	// Atmospheric drag
	if (FlagDrag) a += AccelDrag(Mjd_TT, r, v, T, E, Area_drag, mass, CD);

	// Acceleration
	return a;
}

//------------------------------------------------------------------------------
//
// Deriv
//
// Purpose:
//
//   Computes the derivative of the state vector
//
// Note:
//
//   pAux is expected to point to a variable of type AuxDataRecord, which is
//   used to communicate with the other program sections and to hold data
//   between subsequent calls of this function
//
//------------------------------------------------------------------------------
void Deriv(double t, const Vector& y, Vector& yp, void* pAux)
{
	// Pointer to auxiliary data record
	AuxParam* p = static_cast<AuxParam*>(pAux);

	// Time
	double  Mjd_TT = (*p).Mjd0_TT + t / 86400.0;

	// State vector components
	Vector r = y.slice(0, 2);
	Vector v = y.slice(3, 5);

	// Acceleration
	Vector a(3);

	a = Accel(Mjd_TT, r, v, (*p).Area_drag, (*p).Area_solar, (*p).mass, (*p).CR, (*p).CD,
		(*p).n, (*p).m, (*p).Sun, (*p).Moon, (*p).SRad, (*p).Drag);

	// State vector derivative
	yp = Stack(v, a);
};

//------------------------------------------------------------------------------
//
// Ephemeris computation
//
//------------------------------------------------------------------------------
void Ephemeris(const Vector& Y0, int N_Step, double Step, AuxParam p, Vector Eph[])
{
	int       i;
	double    t, t_end;
	double    relerr, abserr;        // Accuracy requirements
	DE        Orb(Deriv, 6, &p);      // Object for integrating the eq. of motion
	Vector    Y(6);

	relerr = 1.0e-13;
	abserr = 1.0e-6;
	t = 0.0;
	Y = Y0;
	Orb.Init(t, relerr, abserr);
	for (i = 0; i <= N_Step; i++) {
		t_end = Step*i;
		Orb.Integ(t_end, Y);
		Eph[i] = Y;
	};
}

//------------------------------------------------------------------------------
//
// AccessWindows computation
//
//------------------------------------------------------------------------------
int AccessWindowsPoint(Vector SatEci[], Vector SatEcef[], Geodetic Target, Vector Coe[], double Fov[], const int N_Step, const int MaxNumOrbit, double Sway, double Mjd_UTC, double ** AccessTime)
{
	double CentralTime(Vector SatEcef[], int WindowTime[], Geodetic Target, int* CenTime);
	int n = 0; //记录轨道圈数,+1
	int num = 0; //记录时间窗口个数
	int(*OrbitCircle)[2] = new int[MaxNumOrbit][2]; //动态分配轨道圈数个数，最大20圈,每圈存储长度为2的开始和结束时刻序号大小
	double thev, theh;
	double mjdtt, mjdut1, mjdutc;
	int N = N_Step + 1; ////总的轨道个数
	//矩形传感器视场角
	thev = Fov[0] * Rad; //对应STK的Vertical，垂直飞行方向，
	theh = Fov[1] * Rad; //对应STK的Horizontal，平行飞行方向
	// 利用读入的coe矩阵生成轨道根数u，i，om，用于Eci2Ev函数
	double *Arglat = new double[N];
	double *Inclination = new double[N];
	double *Ascend = new double[N];
	int Syear, Smonth, Sday, Shour, Smin, Fyear, Fmonth, Fday, Fhour, Fmin;
	double Ssec, Fsec;
	Matrix pp(3, 3), nn(3, 3), theta(3, 3), ee(3, 3), po(3, 3); //eci2ecef的矩阵
	for (int i = 0; i < N; i++)
	{
		Arglat[i] = fmod(Coe[i](4) + Coe[i](5), pi2*Deg); //纬度幅角 = 真近点角 + 近地点幅角，这里单位用了deg
		Inclination[i] = Coe[i](2);
		Ascend[i] = Coe[i](3);
	}
	// 根据纬度幅角划分一天内轨道运行时段,将一天时间T个数据按纬度幅角为90°分为n个时间段
	OrbitCircle[n][0] = 0; //第一圈开始时刻m的值
	for (int i = 1; i < N; i++)
	{
		if (90 - Arglat[i] > 0.0 && 90 - Arglat[i + 1] <= 0.0)
		{
			OrbitCircle[n][1] = i;
			OrbitCircle[n + 1][0] = i + 1;
			n = n + 1;
		}
	}
	OrbitCircle[n][1] = N;
	//搜索满足可见性的轨道时刻
	double P = 0.0 *Rad; //设置传感器姿态角
	double Y = 0.0 *Rad;
	double R = Sway *Rad;
	Vector TargetEci(3), DeltaEci(3), DeltaOrbit(3), DeltaBody(3), DeltaBodyXOZ(3), DeltaBodyYOZ(3);
	Vector PlaneXOZ(0, 0, 1), PlaneYOZ(0, 0, 1);
	double H, MaxAngle, IncludeAngle, AngleXOZ, AngleYOZ;
	int StartPoint, FinishPoint;
	int step = 15; //粗搜索步长为step * 1秒
	int(*TimeWindow)[4] = new int[n + 1][4]; //动态分配时间窗口数量
	int WindowTime[2];
	int CenTime;
	int Cyear, Cmonth, Cday, Chour, Cmin;
	double Angle, Csec;
	Vector TargetEcef(3);
	TargetEcef = Target.Position(R_Earth, f_Earth); //将目标点经纬高转换为地固坐标XYZ, 单位km
	TargetEcef = operator / (TargetEcef, 1000);
	//根据矩形视场的半角判断点目标可见性,分成n个时间段进行二分法查找
	for (int j = 0; j < n + 1; j++) //轨道圈数=n+1
	{
		int FLAG = 0;
		TimeWindow[j][0] = OrbitCircle[j][0]; //开始时刻的范围
		TimeWindow[j][1] = OrbitCircle[j][0] + step;
		TimeWindow[j][2] = OrbitCircle[j][0]; //结束时刻的范围
		TimeWindow[j][3] = OrbitCircle[j][0] + step;
		int flag = 0; //判断是否找到开始时刻
		for (int i = OrbitCircle[j][0] + step; i <= OrbitCircle[j][1]; i = i + step) //间隔1min，粗搜索, 默认i = 1时不可见
		{
			//判断卫星和目标是否同侧
			H = Norm(SatEcef[i].slice(0, 2)) * 1000; //此时卫星与地心的距离,m
			MaxAngle = acos(R_Earth / H); //卫星与地球的最大夹角
			IncludeAngle = Ang(SatEcef[i].slice(0, 2), TargetEcef);
			//cout << MaxAngle << "," << IncludeAngle << endl;
			//如果卫星和目标是否同侧，则根据卫星与目标夹角范围进行时间窗口计算
			if (IncludeAngle <= MaxAngle) //卫星与目标在同侧
			{
				//计算时间和转换矩阵
				mjdtt = Mjd_UTC + i / 86400.0 + IERS::TT_UTC(Mjd_UTC + i / 86400.0) / 86400.0;
				mjdutc = Mjd_UTC + i / 86400.0;
				mjdut1 = Mjd_UTC + i / 86400.0 + IERS::UT1_UTC(Mjd_UTC + i / 86400.0) / 86400.0;
				po = PoleMatrix(mjdutc);
				pp = PrecMatrix(MJD_J2000, mjdtt);
				nn = NutMatrix(mjdtt);
				theta = GHAMatrix(mjdut1);
				ee = po*theta*nn*pp;
				//计算卫星和地面目标的地固坐标到惯性坐标（R）
				TargetEci = Inv(ee)*TargetEcef;
				DeltaEci = TargetEci - SatEci[i].slice(0, 2); //计算J2000到卫星轨道坐标系的矢量差,km
				//cout << DeltaEci(0) << "," << DeltaEci(1) << "," << DeltaEci(2) << endl;
				//利用轨道根数（Arglat, Inclination, Ascend）, 调用ECIOrbitMatrix函数进行转换得到轨道坐标系系下的下矢量差
				DeltaOrbit = ECIOrbitMatrix(Arglat[i] * Rad, Inclination[i] * Rad, Ascend[i] * Rad)*DeltaEci;
				//利用卫星姿态角（Roll, Pitch, Yaw）, 调用BodyOrbitMatrix函数进行转换得到传感器（卫星本体）姿态下矢量差
				DeltaBody = BodyOrbitMatrix(R, P, Y)*DeltaOrbit;
				DeltaBodyXOZ(0) = DeltaBody(0);
				DeltaBodyXOZ(1) = 0;
				DeltaBodyXOZ(2) = DeltaBody(2);
				DeltaBodyYOZ(0) = 0;
				DeltaBodyYOZ(1) = DeltaBody(1);
				DeltaBodyYOZ(2) = DeltaBody(2);
				//根据矢量差分别求xOz和yOz夹角
				AngleXOZ = Ang(DeltaBodyXOZ, PlaneXOZ);
				AngleYOZ = Ang(DeltaBodyYOZ, PlaneYOZ);
				//将卫星与目标的夹角与视场半角进行比较，记录时间窗口
				if (AngleXOZ <= theh && AngleYOZ <= thev)
				{
					TimeWindow[j][0] = TimeWindow[j][0];
					TimeWindow[j][1] = TimeWindow[j][0] + step;
					flag = 1;
					TimeWindow[j][2] = i;
					TimeWindow[j][3] = i + step;
				}
				else
				{
					if (flag != 1) //如果flag!= 1；说明先不可见
					{
						TimeWindow[j][0] = i;
						TimeWindow[j][1] = i + step;
						TimeWindow[j][2] = i;
						TimeWindow[j][3] = i + step;
					}
					else //如果flag == 1；说明先可见后不可见，找到了开始时刻和结束时刻的范围
					{
						while (TimeWindow[j][1] - TimeWindow[j][0] > 1) //二分法查找具体开始时刻，精度小于1秒停止查找
						{
							int temps = round(TimeWindow[j][0] + (TimeWindow[j][1] - TimeWindow[j][0]) / 2); //中间时刻
							//计算时间和转换矩阵
							mjdtt = Mjd_UTC + temps / 86400.0 + IERS::TT_UTC(Mjd_UTC + temps / 86400.0) / 86400.0;
							mjdutc = Mjd_UTC + temps / 86400.0;
							mjdut1 = Mjd_UTC + temps / 86400.0 + IERS::UT1_UTC(Mjd_UTC + temps / 86400.0) / 86400.0;
							po = PoleMatrix(mjdutc);
							pp = PrecMatrix(MJD_J2000, mjdtt);
							nn = NutMatrix(mjdtt);
							theta = GHAMatrix(mjdut1);
							ee = po*theta*nn*pp;
							//计算卫星和地面目标的地固坐标到惯性坐标（R）
							TargetEci = Inv(ee) * TargetEcef;
							DeltaEci = TargetEci - SatEci[temps].slice(0, 2);
							DeltaOrbit = ECIOrbitMatrix(Arglat[temps] * Rad, Inclination[temps] * Rad, Ascend[temps] * Rad)*DeltaEci;
							//利用卫星姿态角（Roll, Pitch, Yaw）, 调用BodyOrbitMatrix函数进行转换得到传感器（卫星本体）姿态下矢量差
							DeltaBody = BodyOrbitMatrix(R, P, Y)*DeltaOrbit;
							DeltaBodyXOZ(0) = DeltaBody(0);
							DeltaBodyXOZ(1) = 0;
							DeltaBodyXOZ(2) = DeltaBody(2);
							DeltaBodyYOZ(0) = 0;
							DeltaBodyYOZ(1) = DeltaBody(1);
							DeltaBodyYOZ(2) = DeltaBody(2);
							//根据矢量差分别求xOz和yOz夹角
							AngleXOZ = Ang(DeltaBodyXOZ, PlaneXOZ);
							AngleYOZ = Ang(DeltaBodyYOZ, PlaneYOZ);
							if (AngleXOZ <= theh && AngleYOZ <= thev)
								TimeWindow[j][1] = temps;
							else
								TimeWindow[j][0] = temps;
						}
						while (TimeWindow[j][3] - TimeWindow[j][2] > 1) //二分法查找具体结束时刻，精度小于1秒停止查找
						{
							int tempf = round(TimeWindow[j][2] + (TimeWindow[j][3] - TimeWindow[j][2]) / 2); //中间时刻
							//计算时间和转换矩阵
							mjdtt = Mjd_UTC + tempf / 86400.0 + IERS::TT_UTC(Mjd_UTC + tempf / 86400.0) / 86400.0;
							mjdutc = Mjd_UTC + tempf / 86400.0;
							mjdut1 = Mjd_UTC + tempf / 86400.0 + IERS::UT1_UTC(Mjd_UTC + tempf / 86400.0) / 86400.0;
							po = PoleMatrix(mjdutc);
							pp = PrecMatrix(MJD_J2000, mjdtt);
							nn = NutMatrix(mjdtt);
							theta = GHAMatrix(mjdut1);
							ee = po*theta*nn*pp;
							//计算卫星和地面目标的地固坐标到惯性坐标（R）
							TargetEci = Inv(ee) * TargetEcef;
							DeltaEci = TargetEci - SatEci[tempf].slice(0, 2);
							DeltaOrbit = ECIOrbitMatrix(Arglat[tempf] * Rad, Inclination[tempf] * Rad, Ascend[tempf] * Rad)*DeltaEci;
							//利用卫星姿态角（Roll, Pitch, Yaw）, 调用BodyOrbitMatrix函数进行转换得到传感器（卫星本体）姿态下矢量差
							DeltaBody = BodyOrbitMatrix(R, P, Y)*DeltaOrbit;
							DeltaBodyXOZ(0) = DeltaBody(0);
							DeltaBodyXOZ(1) = 0;
							DeltaBodyXOZ(2) = DeltaBody(2);
							DeltaBodyYOZ(0) = 0;
							DeltaBodyYOZ(1) = DeltaBody(1);
							DeltaBodyYOZ(2) = DeltaBody(2);
							//根据矢量差分别求xOz和yOz夹角
							AngleXOZ = Ang(DeltaBodyXOZ, PlaneXOZ);
							AngleYOZ = Ang(DeltaBodyYOZ, PlaneYOZ);
							if (AngleXOZ <= theh && AngleYOZ <= thev)
								TimeWindow[j][2] = tempf;
							else
								TimeWindow[j][3] = tempf;
						}
						StartPoint = round((TimeWindow[j][1] + TimeWindow[j][0]) / 2); //开始时刻
						FinishPoint = round((TimeWindow[j][3] + TimeWindow[j][2]) / 2); //结束时刻
						WindowTime[0] = StartPoint;
						WindowTime[1] = FinishPoint;
						CenTime = StartPoint;
						Angle = CentralTime(SatEcef, WindowTime, Target, &CenTime);
						CalDat((Mjd_UTC + CenTime / 86400.0), Cyear, Cmonth, Cday, Chour, Cmin, Csec);
						//cout << "中心成像时刻：" << Cyear << "-" << Cmonth << "-" << Cday << " " << Chour << ":" << Cmin << ":" << Csec << ",侧摆角为：" << Angle << endl;
						CalDat((Mjd_UTC + StartPoint / 86400.0), Syear, Smonth, Sday, Shour, Smin, Ssec);
						CalDat((Mjd_UTC + FinishPoint / 86400.0), Fyear, Fmonth, Fday, Fhour, Fmin, Fsec);
						AccessTime[num][0] = Syear;
						AccessTime[num][1] = Smonth;
						AccessTime[num][2] = Sday;
						AccessTime[num][3] = Shour;
						AccessTime[num][4] = Smin;
						AccessTime[num][5] = Ssec;
						AccessTime[num][6] = Fyear;
						AccessTime[num][7] = Fmonth;
						AccessTime[num][8] = Fday;
						AccessTime[num][9] = Fhour;
						AccessTime[num][10] = Fmin;
						AccessTime[num][11] = Fsec;
						AccessTime[num][12] = StartPoint;
						AccessTime[num][13] = FinishPoint;
						AccessTime[num][14] = CenTime;
						AccessTime[num][15] = Angle;
						cout << "StartTime:";
						cout << AccessTime[num][0] << "-" << AccessTime[num][1] << "-" << AccessTime[num][2] << " " << AccessTime[num][3] << ":" << AccessTime[num][4] << ":" << AccessTime[num][5] << endl;
						cout << "EndTime:";
						cout << AccessTime[num][6] << "-" << AccessTime[num][7] << "-" << AccessTime[num][8] << " " << AccessTime[num][9] << ":" << AccessTime[num][10] << ":" << AccessTime[num][11] << endl;
						num++;
						FLAG = 1;
						break;
					}

				}

			}

		}
		cout << "OrbitCircle:" << j << endl;
		//if (FLAG == 1)
		//	window_test(j, :) = window_final(1, :); //每圈的可见窗口是在一个侧摆角下只有一行，列代表时间窗口信息
		//else
		//window_test(j, 1:2) = [j, R * 180 / pi];
	}
	//delete[] Arglat; //释放动态内存
	//delete[] Inclination;
	//delete[] Ascend;
	return num;
}

int AccessWindowsArea(Vector SatEci[], Vector SatEcef[], Geodetic Target[], Vector Coe[], double Fov[], const int N_Step, const int MaxNumOrbit, double Sway, double Mjd_UTC, double ** AccessTime) //区域目标的参数初始化及其时间窗口计算
{
	double CentralTime(Vector SatEcef[], int WindowTime[], Geodetic Target, int* CenTime);
	double PointsDistance(double lat1, double lon1, double lat2, double lon2); //声明
	int n = 0; //记录轨道圈数,+1
	int num = 0; //记录时间窗口个数
	int(*OrbitCircle)[2] = new int[MaxNumOrbit][2]; //动态分配轨道圈数个数，最大20圈,每圈存储长度为2的开始和结束时刻序号大小
	double thev, theh;
	double mjdtt, mjdut1, mjdutc;
	int N = N_Step + 1; ////总的轨道个数
	//矩形传感器视场角
	thev = Fov[0] * Rad; //对应STK的Vertical，垂直飞行方向，
	theh = Fov[1] * Rad; //对应STK的Horizontal，平行飞行方向
	// 利用读入的coe矩阵生成轨道根数u，i，om，用于Eci2Ev函数
	double *Arglat = new double[N];
	double *Inclination = new double[N];
	double *Ascend = new double[N];
	int Syear, Smonth, Sday, Shour, Smin, Fyear, Fmonth, Fday, Fhour, Fmin;
	double Ssec, Fsec;
	Matrix pp(3, 3), nn(3, 3), theta(3, 3), ee(3, 3), po(3, 3); //eci2ecef的矩阵
	for (int i = 0; i < N; i++)
	{
		Arglat[i] = fmod(Coe[i](4) + Coe[i](5), pi2*Deg); //纬度幅角 = 真近点角 + 近地点幅角，这里单位用了deg
		Inclination[i] = Coe[i](2);
		Ascend[i] = Coe[i](3);
	}
	// 根据纬度幅角划分一天内轨道运行时段,将一天时间T个数据按纬度幅角为90°分为n个时间段
	OrbitCircle[n][0] = 0; //第一圈开始时刻m的值
	for (int i = 0; i < N; i++)
	{
		if (90 - Arglat[i] > 0.0 && 90 - Arglat[i + 1] <= 0.0)
		{
			OrbitCircle[n][1] = i;
			OrbitCircle[n + 1][0] = i + 1;
			n = n + 1;
		}
	}
	OrbitCircle[n][1] = N; //最大圈数标号实际为N-1
	/*for (int test = 0; test <= n; test++)
	{
	cout << OrbitCircle[test][0] << "," << OrbitCircle[test][1] << endl;
	}*/
	//搜索满足可见性的轨道时刻
	double P = 0.0 *Rad; //设置传感器姿态角
	double Y = 0.0 *Rad;
	double R = Sway *Rad;
	Vector PlaneXOZ(0, 0, 1), PlaneYOZ(0, 0, 1);
	int(*TimeWindow)[4] = new int[n + 1][4]; //动态分配时间窗口数量
	double H, MaxAngle, IncludeAngle1, IncludeAngle2, IncludeAngle3, IncludeAngle4;
	int StartPoint, FinishPoint;
	int WindowTime[2];
	int CenTime;
	int Cyear, Cmonth, Cday, Chour, Cmin;
	double Angle, Csec;
	int step = 40;
	//根据矩形视场的半角判断点目标可见性,分成n个时间段进行二分法查找
	//根据区域目标大小确定采样点个数
	double dis;
	dis = PointsDistance(Target[0].lat, Target[0].lon, Target[1].lat, Target[1].lon); //求区域目标上边两点顶点的距离
	//const int m = floor(dis); //1公里1个点
	const int m = 1000;
	double lonstep = abs(Target[0].lon - Target[1].lon) / (m + 0.0);
	double latstep = abs(Target[0].lat - Target[2].lat) / (m + 0.0);
	Geodetic Top[m + 1];
	Geodetic Bot[m + 1];
	Geodetic Lef[m + 1];
	Geodetic Rig[m + 1];
	Vector TopEcef[m + 1];
	Vector BotEcef[m + 1];
	Vector LefEcef[m + 1];
	Vector RigEcef[m + 1];
	Vector TopEci[m + 1];
	Vector BotEci[m + 1];
	Vector LefEci[m + 1];
	Vector RigEci[m + 1];
	Vector DeltaTopEci[m + 1];
	Vector DeltaBotEci[m + 1];
	Vector DeltaLefEci[m + 1];
	Vector DeltaRigEci[m + 1];
	Vector DeltaTopOrbit[m + 1];
	Vector DeltaBotOrbit[m + 1];
	Vector DeltaLefOrbit[m + 1];
	Vector DeltaRigOrbit[m + 1];
	Vector DeltaTopBody[m + 1];
	Vector DeltaBotBody[m + 1];
	Vector DeltaLefBody[m + 1];
	Vector DeltaRigBody[m + 1];
	Vector DeltaTopBodyXOZ(3);
	Vector DeltaTopBodyYOZ(3);
	Vector DeltaBotBodyXOZ(3);
	Vector DeltaBotBodyYOZ(3);
	Vector DeltaLefBodyXOZ(3);
	Vector DeltaLefBodyYOZ(3);
	Vector DeltaRigBodyXOZ(3);
	Vector DeltaRigBodyYOZ(3);
	double AngleXOZTop[m + 1];
	double AngleYOZTop[m + 1];
	double AngleXOZBot[m + 1];
	double AngleYOZBot[m + 1];
	double AngleXOZLef[m + 1];
	double AngleYOZLef[m + 1];
	double AngleXOZRig[m + 1];
	double AngleYOZRig[m + 1];
	Geodetic CenTarget;
	CenTarget.lat = (Target[0].lat + Target[2].lat) / 2;
	CenTarget.lon = (Target[0].lon + Target[2].lon) / 2;
	CenTarget.h = (Target[0].h + Target[2].h) / 2;
	//矩形区域目标每个边采样得到上下左右四个边的点集合(经纬度和地固坐标)
	for (int k = 0; k <= m; k++)
	{
		Top[k].lon = Target[0].lon + k*lonstep;
		Top[k].lat = Target[0].lat;
		Top[k].h = 0;
		TopEcef[k] = Top[k].Position(R_Earth, f_Earth);
		TopEcef[k] = operator / (TopEcef[k], 1000);

		Bot[k].lon = Target[3].lon + k*lonstep;
		Bot[k].lat = Target[3].lat;
		Bot[k].h = 0;
		BotEcef[k] = Bot[k].Position(R_Earth, f_Earth);
		BotEcef[k] = operator / (BotEcef[k], 1000);

		Lef[k].lon = Target[3].lon;
		Lef[k].lat = Target[3].lat + k*latstep;
		Lef[k].h = 0;
		LefEcef[k] = Lef[k].Position(R_Earth, f_Earth);
		LefEcef[k] = operator / (LefEcef[k], 1000);

		Rig[k].lon = Target[2].lon;
		Rig[k].lat = Target[2].lat + k*latstep;
		Rig[k].h = 0;
		RigEcef[k] = Rig[k].Position(R_Earth, f_Earth);
		RigEcef[k] = operator / (RigEcef[k], 1000);
	}
	for (int j = 0; j < n + 1; j++)
	{
		int FLAG = 0;
		TimeWindow[j][0] = OrbitCircle[j][0]; //开始时刻的范围
		TimeWindow[j][1] = OrbitCircle[j][0] + step;
		TimeWindow[j][2] = OrbitCircle[j][0]; //结束时刻的范围
		TimeWindow[j][3] = OrbitCircle[j][0] + step;
		int flag = 0; //判断是否找到开始时刻
		for (int i = OrbitCircle[j][0] + step; i < OrbitCircle[j][1]; i = i + step) //间隔1min，粗搜索, 默认i = 1时不可见
		{
			//判断卫星和目标是否同侧
			H = Norm(SatEcef[i].slice(0, 2)) * 1000; //此时卫星与地心的距离,m
			MaxAngle = acos(R_Earth / H); //卫星与地球的最大夹角
			IncludeAngle1 = Ang(SatEcef[i].slice(0, 2), TopEcef[0]);
			IncludeAngle2 = Ang(SatEcef[i].slice(0, 2), BotEcef[0]);
			IncludeAngle3 = Ang(SatEcef[i].slice(0, 2), RigEcef[0]);
			IncludeAngle4 = Ang(SatEcef[i].slice(0, 2), RigEcef[m]);
			//如果卫星和目标是否同侧，则根据卫星与目标夹角范围进行时间窗口计算
			if (IncludeAngle1 <= MaxAngle || IncludeAngle2 <= MaxAngle || IncludeAngle3 <= MaxAngle || IncludeAngle4 <= MaxAngle) //卫星与目标在同侧
			{
				int k;
				for (k = 0; k <= m; k++)
				{
					//计算时间和转换矩阵
					mjdtt = Mjd_UTC + i / 86400.0 + IERS::TT_UTC(Mjd_UTC + i / 86400.0) / 86400.0;
					mjdutc = Mjd_UTC + i / 86400.0;
					mjdut1 = Mjd_UTC + i / 86400.0 + IERS::UT1_UTC(Mjd_UTC + i / 86400.0) / 86400.0;
					po = PoleMatrix(mjdutc);
					pp = PrecMatrix(MJD_J2000, mjdtt);
					nn = NutMatrix(mjdtt);
					theta = GHAMatrix(mjdut1);
					ee = po*theta*nn*pp;
					//计算卫星和地面目标的地固坐标到惯性坐标（R）
					TopEci[k] = Inv(ee)*TopEcef[k];
					BotEci[k] = Inv(ee)*BotEcef[k];
					LefEci[k] = Inv(ee)*LefEcef[k];
					RigEci[k] = Inv(ee)*RigEcef[k];
					DeltaTopEci[k] = TopEci[k] - SatEci[i].slice(0, 2); //计算J2000到卫星轨道坐标系的矢量差,km
					DeltaBotEci[k] = BotEci[k] - SatEci[i].slice(0, 2);
					DeltaLefEci[k] = LefEci[k] - SatEci[i].slice(0, 2);
					DeltaRigEci[k] = RigEci[k] - SatEci[i].slice(0, 2);
					//利用轨道根数（Arglat, Inclination, Ascend）, 调用ECIOrbitMatrix函数进行转换得到轨道坐标系系下的下矢量差
					DeltaTopOrbit[k] = ECIOrbitMatrix(Arglat[i] * Rad, Inclination[i] * Rad, Ascend[i] * Rad)*DeltaTopEci[k];
					DeltaBotOrbit[k] = ECIOrbitMatrix(Arglat[i] * Rad, Inclination[i] * Rad, Ascend[i] * Rad)*DeltaBotEci[k];
					DeltaLefOrbit[k] = ECIOrbitMatrix(Arglat[i] * Rad, Inclination[i] * Rad, Ascend[i] * Rad)*DeltaLefEci[k];
					DeltaRigOrbit[k] = ECIOrbitMatrix(Arglat[i] * Rad, Inclination[i] * Rad, Ascend[i] * Rad)*DeltaRigEci[k];
					//利用卫星姿态角（Roll, Pitch, Yaw）, 调用BodyOrbitMatrix函数进行转换得到传感器（卫星本体）姿态下矢量差
					DeltaTopBody[k] = BodyOrbitMatrix(R, P, Y)*DeltaTopOrbit[k];
					DeltaBotBody[k] = BodyOrbitMatrix(R, P, Y)*DeltaBotOrbit[k];
					DeltaLefBody[k] = BodyOrbitMatrix(R, P, Y)*DeltaLefOrbit[k];
					DeltaRigBody[k] = BodyOrbitMatrix(R, P, Y)*DeltaRigOrbit[k];
					//根据矢量差分别求xOz和yOz夹角
					DeltaTopBodyXOZ(0) = DeltaTopBody[k](0);
					DeltaTopBodyXOZ(1) = 0;
					DeltaTopBodyXOZ(2) = DeltaTopBody[k](2);
					DeltaTopBodyYOZ(0) = 0;
					DeltaTopBodyYOZ(1) = DeltaTopBody[k](1);
					DeltaTopBodyYOZ(2) = DeltaTopBody[k](2);

					DeltaBotBodyXOZ(0) = DeltaBotBody[k](0);
					DeltaBotBodyXOZ(1) = 0;
					DeltaBotBodyXOZ(2) = DeltaBotBody[k](2);
					DeltaBotBodyYOZ(0) = 0;
					DeltaBotBodyYOZ(1) = DeltaBotBody[k](1);
					DeltaBotBodyYOZ(2) = DeltaBotBody[k](2);

					DeltaLefBodyXOZ(0) = DeltaLefBody[k](0);
					DeltaLefBodyXOZ(1) = 0;
					DeltaLefBodyXOZ(2) = DeltaLefBody[k](2);
					DeltaLefBodyYOZ(0) = 0;
					DeltaLefBodyYOZ(1) = DeltaLefBody[k](1);
					DeltaLefBodyYOZ(2) = DeltaLefBody[k](2);

					DeltaRigBodyXOZ(0) = DeltaRigBody[k](0);
					DeltaRigBodyXOZ(1) = 0;
					DeltaRigBodyXOZ(2) = DeltaRigBody[k](2);
					DeltaRigBodyYOZ(0) = 0;
					DeltaRigBodyYOZ(1) = DeltaRigBody[k](1);
					DeltaRigBodyYOZ(2) = DeltaRigBody[k](2);

					AngleXOZTop[k] = Ang(DeltaTopBodyXOZ, PlaneXOZ);
					AngleYOZTop[k] = Ang(DeltaTopBodyYOZ, PlaneYOZ);
					AngleXOZBot[k] = Ang(DeltaBotBodyXOZ, PlaneXOZ);
					AngleYOZBot[k] = Ang(DeltaBotBodyYOZ, PlaneYOZ);
					AngleXOZLef[k] = Ang(DeltaLefBodyXOZ, PlaneXOZ);
					AngleYOZLef[k] = Ang(DeltaLefBodyYOZ, PlaneYOZ);
					AngleXOZRig[k] = Ang(DeltaRigBodyXOZ, PlaneXOZ);
					AngleYOZRig[k] = Ang(DeltaRigBodyYOZ, PlaneYOZ);
					if (AngleXOZTop[k] <= theh && AngleYOZTop[k] <= thev || AngleXOZBot[k] <= theh && AngleYOZBot[k] <= thev || AngleXOZLef[k] <= theh && AngleYOZLef[k] <= thev || AngleXOZRig[k] <= theh && AngleYOZRig[k] <= thev)
					{
						TimeWindow[j][0] = TimeWindow[j][0];
						TimeWindow[j][1] = TimeWindow[j][0] + step;
						flag = 1;
						TimeWindow[j][2] = i;
						TimeWindow[j][3] = i + step;
						break;
					}
					else
					{
						if (k < m) //边界点采样数还没有用完，则继续下一次边界循环,k最大值为m
							continue;
						else if (flag != 1 && k == m) //如果flag!= 1,说明第一次循环先不可见，并且边界点采样数用完边界始终不可见，则改变时间，继续下一次外循环
						{
							TimeWindow[j][0] = i;
							TimeWindow[j][1] = i + step;
							TimeWindow[j][2] = i;
							TimeWindow[j][3] = i + step;
						}
						else //如果flag == 1；说明先可见后不可见，找到了开始时刻和结束时刻的范围
						{
							while (TimeWindow[j][1] - TimeWindow[j][0] > 1) //二分法查找具体开始时刻，精度小于1秒停止查找
							{
								int temps = round(TimeWindow[j][0] + (TimeWindow[j][1] - TimeWindow[j][0]) / 2); //中间时刻
								for (int k1 = 0; k1 <= m; k1++)
								{
									//计算时间和转换矩阵
									mjdtt = Mjd_UTC + temps / 86400.0 + IERS::TT_UTC(Mjd_UTC + temps / 86400.0) / 86400.0;
									mjdutc = Mjd_UTC + temps / 86400.0;
									mjdut1 = Mjd_UTC + temps / 86400.0 + IERS::UT1_UTC(Mjd_UTC + temps / 86400.0) / 86400.0;
									po = PoleMatrix(mjdutc);
									pp = PrecMatrix(MJD_J2000, mjdtt);
									nn = NutMatrix(mjdtt);
									theta = GHAMatrix(mjdut1);
									ee = po*theta*nn*pp;
									//计算卫星和地面目标的地固坐标到惯性坐标（R）
									TopEci[k1] = Inv(ee) * TopEcef[k1];
									BotEci[k1] = Inv(ee) * BotEcef[k1];
									LefEci[k1] = Inv(ee) * LefEcef[k1];
									RigEci[k1] = Inv(ee) * RigEcef[k1];
									DeltaTopEci[k1] = TopEci[k1] - SatEci[temps].slice(0, 2);
									DeltaBotEci[k1] = BotEci[k1] - SatEci[temps].slice(0, 2);
									DeltaLefEci[k1] = LefEci[k1] - SatEci[temps].slice(0, 2);
									DeltaRigEci[k1] = RigEci[k1] - SatEci[temps].slice(0, 2);
									DeltaTopOrbit[k1] = ECIOrbitMatrix(Arglat[temps] * Rad, Inclination[temps] * Rad, Ascend[temps] * Rad)*DeltaTopEci[k1];
									DeltaBotOrbit[k1] = ECIOrbitMatrix(Arglat[temps] * Rad, Inclination[temps] * Rad, Ascend[temps] * Rad)*DeltaBotEci[k1];
									DeltaLefOrbit[k1] = ECIOrbitMatrix(Arglat[temps] * Rad, Inclination[temps] * Rad, Ascend[temps] * Rad)*DeltaLefEci[k1];
									DeltaRigOrbit[k1] = ECIOrbitMatrix(Arglat[temps] * Rad, Inclination[temps] * Rad, Ascend[temps] * Rad)*DeltaRigEci[k1];
									//利用卫星姿态角（Roll, Pitch, Yaw）, 调用BodyOrbitMatrix函数进行转换得到传感器（卫星本体）姿态下矢量差
									DeltaTopBody[k1] = BodyOrbitMatrix(R, P, Y)*DeltaTopOrbit[k1];
									DeltaBotBody[k1] = BodyOrbitMatrix(R, P, Y)*DeltaBotOrbit[k1];
									DeltaLefBody[k1] = BodyOrbitMatrix(R, P, Y)*DeltaLefOrbit[k1];
									DeltaRigBody[k1] = BodyOrbitMatrix(R, P, Y)*DeltaRigOrbit[k1];
									//根据矢量差分别求xOz和yOz夹角
									DeltaTopBodyXOZ(0) = DeltaTopBody[k1](0);
									DeltaTopBodyXOZ(1) = 0;
									DeltaTopBodyXOZ(2) = DeltaTopBody[k1](2);
									DeltaTopBodyYOZ(0) = 0;
									DeltaTopBodyYOZ(1) = DeltaTopBody[k1](1);
									DeltaTopBodyYOZ(2) = DeltaTopBody[k1](2);

									DeltaBotBodyXOZ(0) = DeltaBotBody[k1](0);
									DeltaBotBodyXOZ(1) = 0;
									DeltaBotBodyXOZ(2) = DeltaBotBody[k1](2);
									DeltaBotBodyYOZ(0) = 0;
									DeltaBotBodyYOZ(1) = DeltaBotBody[k1](1);
									DeltaBotBodyYOZ(2) = DeltaBotBody[k1](2);

									DeltaLefBodyXOZ(0) = DeltaLefBody[k1](0);
									DeltaLefBodyXOZ(1) = 0;
									DeltaLefBodyXOZ(2) = DeltaLefBody[k1](2);
									DeltaLefBodyYOZ(0) = 0;
									DeltaLefBodyYOZ(1) = DeltaLefBody[k1](1);
									DeltaLefBodyYOZ(2) = DeltaLefBody[k1](2);

									DeltaRigBodyXOZ(0) = DeltaRigBody[k1](0);
									DeltaRigBodyXOZ(1) = 0;
									DeltaRigBodyXOZ(2) = DeltaRigBody[k1](2);
									DeltaRigBodyYOZ(0) = 0;
									DeltaRigBodyYOZ(1) = DeltaRigBody[k1](1);
									DeltaRigBodyYOZ(2) = DeltaRigBody[k1](2);

									AngleXOZTop[k1] = Ang(DeltaTopBodyXOZ, PlaneXOZ);
									AngleYOZTop[k1] = Ang(DeltaTopBodyYOZ, PlaneYOZ);
									AngleXOZBot[k1] = Ang(DeltaBotBodyXOZ, PlaneXOZ);
									AngleYOZBot[k1] = Ang(DeltaBotBodyYOZ, PlaneYOZ);
									AngleXOZLef[k1] = Ang(DeltaLefBodyXOZ, PlaneXOZ);
									AngleYOZLef[k1] = Ang(DeltaLefBodyYOZ, PlaneYOZ);
									AngleXOZRig[k1] = Ang(DeltaRigBodyXOZ, PlaneXOZ);
									AngleYOZRig[k1] = Ang(DeltaRigBodyYOZ, PlaneYOZ);
									if (AngleXOZTop[k1] <= theh && AngleYOZTop[k1] <= thev || AngleXOZBot[k1] <= theh && AngleYOZBot[k1] <= thev || AngleXOZLef[k1] <= theh && AngleYOZLef[k1] <= thev || AngleXOZRig[k1] <= theh && AngleYOZRig[k1] <= thev)
									{
										TimeWindow[j][1] = temps;
										break;
									}
									else if (k1 < m)
										continue;
									else
										TimeWindow[j][0] = temps;
								}
							}
							while (TimeWindow[j][3] - TimeWindow[j][2] > 1) //二分法查找具体结束时刻，精度小于1秒停止查找
							{
								int tempf = round(TimeWindow[j][2] + (TimeWindow[j][3] - TimeWindow[j][2]) / 2); //中间时刻
								for (int k1 = 0; k1 <= m; k1++)
								{
									//计算时间和转换矩阵
									mjdtt = Mjd_UTC + tempf / 86400.0 + IERS::TT_UTC(Mjd_UTC + tempf / 86400.0) / 86400.0;
									mjdutc = Mjd_UTC + tempf / 86400.0;
									mjdut1 = Mjd_UTC + tempf / 86400.0 + IERS::UT1_UTC(Mjd_UTC + tempf / 86400.0) / 86400.0;
									po = PoleMatrix(mjdutc);
									pp = PrecMatrix(MJD_J2000, mjdtt);
									nn = NutMatrix(mjdtt);
									theta = GHAMatrix(mjdut1);
									ee = po*theta*nn*pp;
									//计算卫星和地面目标的地固坐标到惯性坐标（R）
									TopEci[k1] = Inv(ee) * TopEcef[k1];
									BotEci[k1] = Inv(ee) * BotEcef[k1];
									LefEci[k1] = Inv(ee) * LefEcef[k1];
									RigEci[k1] = Inv(ee) * RigEcef[k1];
									DeltaTopEci[k1] = TopEci[k1] - SatEci[tempf].slice(0, 2);
									DeltaBotEci[k1] = BotEci[k1] - SatEci[tempf].slice(0, 2);
									DeltaLefEci[k1] = LefEci[k1] - SatEci[tempf].slice(0, 2);
									DeltaRigEci[k1] = RigEci[k1] - SatEci[tempf].slice(0, 2);
									DeltaTopOrbit[k1] = ECIOrbitMatrix(Arglat[tempf] * Rad, Inclination[tempf] * Rad, Ascend[tempf] * Rad)*DeltaTopEci[k1];
									DeltaBotOrbit[k1] = ECIOrbitMatrix(Arglat[tempf] * Rad, Inclination[tempf] * Rad, Ascend[tempf] * Rad)*DeltaBotEci[k1];
									DeltaLefOrbit[k1] = ECIOrbitMatrix(Arglat[tempf] * Rad, Inclination[tempf] * Rad, Ascend[tempf] * Rad)*DeltaLefEci[k1];
									DeltaRigOrbit[k1] = ECIOrbitMatrix(Arglat[tempf] * Rad, Inclination[tempf] * Rad, Ascend[tempf] * Rad)*DeltaRigEci[k1];
									//利用卫星姿态角（Roll, Pitch, Yaw）, 调用BodyOrbitMatrix函数进行转换得到传感器（卫星本体）姿态下矢量差
									DeltaTopBody[k1] = BodyOrbitMatrix(R, P, Y)*DeltaTopOrbit[k1];
									DeltaBotBody[k1] = BodyOrbitMatrix(R, P, Y)*DeltaBotOrbit[k1];
									DeltaLefBody[k1] = BodyOrbitMatrix(R, P, Y)*DeltaLefOrbit[k1];
									DeltaRigBody[k1] = BodyOrbitMatrix(R, P, Y)*DeltaRigOrbit[k1];
									//根据矢量差分别求xOz和yOz夹角
									DeltaTopBodyXOZ(0) = DeltaTopBody[k1](0);
									DeltaTopBodyXOZ(1) = 0;
									DeltaTopBodyXOZ(2) = DeltaTopBody[k1](2);
									DeltaTopBodyYOZ(0) = 0;
									DeltaTopBodyYOZ(1) = DeltaTopBody[k1](1);
									DeltaTopBodyYOZ(2) = DeltaTopBody[k1](2);
									DeltaBotBodyXOZ(0) = DeltaBotBody[k1](0);
									DeltaBotBodyXOZ(1) = 0;
									DeltaBotBodyXOZ(2) = DeltaBotBody[k1](2);
									DeltaBotBodyYOZ(0) = 0;
									DeltaBotBodyYOZ(1) = DeltaBotBody[k1](1);
									DeltaBotBodyYOZ(2) = DeltaBotBody[k1](2);
									DeltaLefBodyXOZ(0) = DeltaLefBody[k1](0);
									DeltaLefBodyXOZ(1) = 0;
									DeltaLefBodyXOZ(2) = DeltaLefBody[k1](2);
									DeltaLefBodyYOZ(0) = 0;
									DeltaLefBodyYOZ(1) = DeltaLefBody[k1](1);
									DeltaLefBodyYOZ(2) = DeltaLefBody[k1](2);
									DeltaRigBodyXOZ(0) = DeltaRigBody[k1](0);
									DeltaRigBodyXOZ(1) = 0;
									DeltaRigBodyXOZ(2) = DeltaRigBody[k1](2);
									DeltaRigBodyYOZ(0) = 0;
									DeltaRigBodyYOZ(1) = DeltaRigBody[k1](1);
									DeltaRigBodyYOZ(2) = DeltaRigBody[k1](2);
									AngleXOZTop[k1] = Ang(DeltaTopBodyXOZ, PlaneXOZ);
									AngleYOZTop[k1] = Ang(DeltaTopBodyYOZ, PlaneYOZ);
									AngleXOZBot[k1] = Ang(DeltaBotBodyXOZ, PlaneXOZ);
									AngleYOZBot[k1] = Ang(DeltaBotBodyYOZ, PlaneYOZ);
									AngleXOZLef[k1] = Ang(DeltaLefBodyXOZ, PlaneXOZ);
									AngleYOZLef[k1] = Ang(DeltaLefBodyYOZ, PlaneYOZ);
									AngleXOZRig[k1] = Ang(DeltaRigBodyXOZ, PlaneXOZ);
									AngleYOZRig[k1] = Ang(DeltaRigBodyYOZ, PlaneYOZ);
									if (AngleXOZTop[k1] <= theh && AngleYOZTop[k1] <= thev || AngleXOZBot[k1] <= theh && AngleYOZBot[k1] <= thev || AngleXOZLef[k1] <= theh && AngleYOZLef[k1] <= thev || AngleXOZRig[k1] <= theh && AngleYOZRig[k1] <= thev)
									{
										TimeWindow[j][2] = tempf;
										break;
									}
									else if (k1 < m)
										continue;
									else
										TimeWindow[j][3] = tempf;
								}
							}
							StartPoint = round((TimeWindow[j][1] + TimeWindow[j][0]) / 2.0); //开始时刻
							FinishPoint = round((TimeWindow[j][3] + TimeWindow[j][2]) / 2.0); //结束时刻
							WindowTime[0] = StartPoint;
							WindowTime[1] = FinishPoint;
							CenTime = StartPoint;
							Angle = CentralTime(SatEcef, WindowTime, CenTarget, &CenTime);
							CalDat((Mjd_UTC + CenTime / 86400.0), Cyear, Cmonth, Cday, Chour, Cmin, Csec);
							//cout << "中心成像时刻：" << Cyear << "-" << Cmonth << "-" << Cday << " " << Chour << ":" << Cmin << ":" << Csec << ",侧摆角为：" << Angle << endl;
							CalDat((Mjd_UTC + StartPoint / 86400.0), Syear, Smonth, Sday, Shour, Smin, Ssec);
							CalDat((Mjd_UTC + FinishPoint / 86400.0), Fyear, Fmonth, Fday, Fhour, Fmin, Fsec);
							AccessTime[num][0] = Syear;
							AccessTime[num][1] = Smonth;
							AccessTime[num][2] = Sday;
							AccessTime[num][3] = Shour;
							AccessTime[num][4] = Smin;
							AccessTime[num][5] = Ssec;
							AccessTime[num][6] = Fyear;
							AccessTime[num][7] = Fmonth;
							AccessTime[num][8] = Fday;
							AccessTime[num][9] = Fhour;
							AccessTime[num][10] = Fmin;
							AccessTime[num][11] = Fsec;
							AccessTime[num][12] = StartPoint;
							AccessTime[num][13] = FinishPoint;
							AccessTime[num][14] = CenTime;
							AccessTime[num][15] = Angle;
							cout << "StartTime:";
							cout << AccessTime[num][0] << "-" << AccessTime[num][1] << "-" << AccessTime[num][2] << " " << AccessTime[num][3] << ":" << AccessTime[num][4] << ":" << AccessTime[num][5] << endl;
							cout << "EndTime:";
							cout << AccessTime[num][6] << "-" << AccessTime[num][7] << "-" << AccessTime[num][8] << " " << AccessTime[num][9] << ":" << AccessTime[num][10] << ":" << AccessTime[num][11] << endl;
							num++;
							FLAG = 1;
							break;

						}
					}
				}
				if (flag == 1 && k == m)
					break; //跳出外循环，避免重复的 window_final赋值
			}

		}
                cout << "OrbitCircle:" << j << endl;
	}
	return num;
}

//------------------------------------------------------------------------------
//
// 根据视线法求解可见窗口下传感器覆盖，开始结束两个时刻形成的大条带
//
//------------------------------------------------------------------------------
// 
void CalculateCov(Vector SatEci[], Vector Coe[], int CovTime[], double Fov[], double Sway, double Mjd_UTC, Geodetic CovBoundary[])
{
	void SensorCov(Vector SatEci[], Vector Coe[], int CovTime, double Fov[], double Sway, double Mjd_UTC, Geodetic Coverage[]); //声明函数
	Geodetic StartCov[4], FinishCov[4];
	SensorCov(SatEci, Coe, CovTime[0], Fov, Sway, Mjd_UTC, StartCov);
	SensorCov(SatEci, Coe, CovTime[1], Fov, Sway, Mjd_UTC, FinishCov);
	for (int i = 0; i < 2; i++)
	{
		CovBoundary[i].h = StartCov[i + 2].h;
		CovBoundary[i].lat = StartCov[i + 2].lat * Deg;
		CovBoundary[i].lon = StartCov[i + 2].lon * Deg;
		CovBoundary[i + 2].h = FinishCov[3 - i].h;
		CovBoundary[i + 2].lat = FinishCov[3 - i].lat * Deg;
		CovBoundary[i + 2].lon = FinishCov[3 - i].lon * Deg;
	}
}

//------------------------------------------------------------------------------
//
// 传感器覆盖区域边界坐标求解, 通过卫星传感器（四个边界点）视线方程计算与地球的交点坐标确定地面覆盖范围
//
//------------------------------------------------------------------------------
// 
void SensorCov(Vector SatEci[], Vector Coe[], int CovTime, double Fov[], double Sway, double Mjd_UTC, Geodetic Coverage[])
{
	//设置成像点的坐标, 长短轴四个端点的视场范围
	double Arglat, Inclination, Ascend;
	double R, P, Y, theh, thev;
	double a, b, temp; //地球长短半径
	double A, B, C, Delta; //以Z为变量的一元二次方程的系数
	double xs, ys, zs, l1, l2, l3; //椭球和射线方程的系数
	Vector cov1(3), cov2(3), s1(3), s2(3);
	double p[4];
	Vector A1(3), A2(3), A3(3), A4(3), S1(3), S2(3), S3(3), S4(3);
	Vector l[4];

	double mjdtt, mjdut1, mjdutc;
	Matrix pp(3, 3), nn(3, 3), theta(3, 3), ee(3, 3), po(3, 3); //eci2ecef的矩阵
	Arglat = fmod(Coe[CovTime](4) + Coe[CovTime](5), pi2*Deg); //纬度幅角 = 真近点角 + 近地点幅角，这里单位用了deg
	Inclination = Coe[CovTime](2);
	Ascend = Coe[CovTime](3);
	R = Sway * Rad;
	P = 0.0 * Rad;
	Y = 0.0 * Rad;
	thev = Fov[0] * Rad;
	theh = Fov[1] * Rad;
	a = R_Earth;//地球椭球的半长轴和半短轴的长度,m
	b = RS_Earth;
	Vector CovEcef[4];
	Vector CovEci[4];
	if (theh == 0) //圆锥形传感器求其投影椭圆长短半轴端点为边界点
	{
		p[2] = 0; //输出
		p[0] = pi / 2;
		p[3] = pi2 / 2; //输出
		p[1] = pi3 / 2;

		S1(0) = tan(thev)*sin(p[0]); //框幅式成像的四个端点
		S1(1) = tan(thev)*cos(p[0]);
		S1(2) = 1;
		S2(0) = tan(thev)*sin(p[1]);
		S2(1) = tan(thev)*cos(p[1]);
		S2(2) = 1;
		S3(0) = tan(thev)*sin(p[2]);
		S3(1) = tan(thev)*cos(p[2]);
		S3(2) = 1;
		S4(0) = tan(thev)*sin(p[3]);
		S4(1) = tan(thev)*cos(p[3]);
		S4(2) = 1;
		A1 = Inv(BodyOrbitMatrix(R, P, Y)) * S1;
		A2 = Inv(BodyOrbitMatrix(R, P, Y)) * S2;
		A3 = Inv(BodyOrbitMatrix(R, P, Y)) * S3;
		A4 = Inv(BodyOrbitMatrix(R, P, Y)) * S4;
	}
	else   //矩形传感器视线矢量表示
	{
		p[0] = acos(tan(thev) / sqrt(tan(thev)*tan(thev) + tan(theh)*tan(theh))); //矩形视场四个端点的角度，以垂直飞行方向为0度
		p[1] = pi - p[0];
		p[2] = pi + p[0]; //前边两个点
		p[3] = -p[0];
		temp = tan(thev) / sqrt(tan(thev)*tan(thev) + tan(theh)*tan(theh));
		//四个端点的卫星本体坐标到卫星轨道坐标系的转换，默认传感器安装矩阵为3 * 3的全1矩阵
		S1(0) = tan(thev) / temp * sin(p[0]); //矩形视场成像的四个端点
		S1(1) = tan(thev) / temp * cos(p[0]);
		S1(2) = 1;
		S2(0) = tan(thev) / temp * sin(p[1]);
		S2(1) = tan(thev) / temp * cos(p[1]);
		S2(2) = 1;
		S3(0) = tan(thev) / temp * sin(p[2]);
		S3(1) = tan(thev) / temp * cos(p[2]);
		S3(2) = 1;
		S4(0) = tan(thev) / temp * sin(p[3]);
		S4(1) = tan(thev) / temp * cos(p[3]);
		S4(2) = 1;
		A1 = Inv(BodyOrbitMatrix(R, P, Y)) * S1; //利用圆锥视场求矩形视场
		A2 = Inv(BodyOrbitMatrix(R, P, Y)) * S2;
		A3 = Inv(BodyOrbitMatrix(R, P, Y)) * S3;
		A4 = Inv(BodyOrbitMatrix(R, P, Y)) * S4;
	}
	//卫星轨道坐标系转换到地心惯性坐标系
	l[0] = Inv(ECIOrbitMatrix(Arglat * Rad, Inclination * Rad, Ascend * Rad))*A1;
	l[1] = Inv(ECIOrbitMatrix(Arglat * Rad, Inclination * Rad, Ascend * Rad))*A2;
	l[2] = Inv(ECIOrbitMatrix(Arglat * Rad, Inclination * Rad, Ascend * Rad))*A3;
	l[3] = Inv(ECIOrbitMatrix(Arglat * Rad, Inclination * Rad, Ascend * Rad))*A4;
	//此时刻卫星的eci坐标
	xs = SatEci[CovTime](0) * 1000; //[m]
	ys = SatEci[CovTime](1) * 1000;
	zs = SatEci[CovTime](2) * 1000;
	for (int j = 0; j < 4; j++)
	{
		//四个端点分别取对应的转换向量
		l1 = l[j](0);
		l2 = l[j](1);
		l3 = l[j](2);
		A = (b*b)*(l1*l1 + l2*l2) + (a*a)*(l3*l3);
		B = 2 * (b*b)*((l1*xs + l2*ys)*l3 - (l1*l1 + l2*l2)*zs);
		C = (b*b)*((l1*l1 + l2*l2)*(zs*zs) + (xs*xs + ys*ys)*(l3*l3) - 2 * (l1*xs + l2*ys)*l3*zs - (a*a)*(l3*l3));
		Delta = B * B - 4 * A * C;
		cov1(2) = (-B + sqrt(Delta)) / (2 * A);  //z
		cov1(0) = (cov1(2) - zs)*l1 / l3 + xs; //x
		cov1(1) = (cov1(2) - zs)*l2 / l3 + ys; //y
		cov2(2) = (-B - sqrt(Delta)) / (2 * A);
		cov2(0) = (cov2(2) - zs)*l1 / l3 + xs;
		cov2(1) = (cov2(2) - zs)*l2 / l3 + ys;
		SatEci[CovTime] = operator * (1000, SatEci[CovTime]);
		s1 = operator - (cov1, SatEci[CovTime].slice(0, 2));
		s2 = operator - (cov2, SatEci[CovTime].slice(0, 2));
		if (Norm(s1) <= Norm(s2)) //[m]
			CovEci[j] = cov1;
		else
			CovEci[j] = cov2;
		//计算时间和转换矩阵
		mjdtt = Mjd_UTC + CovTime / 86400.0 + IERS::TT_UTC(Mjd_UTC + CovTime / 86400.0) / 86400.0;
		mjdutc = Mjd_UTC + CovTime / 86400.0;
		mjdut1 = Mjd_UTC + CovTime / 86400.0 + IERS::UT1_UTC(Mjd_UTC + CovTime / 86400.0) / 86400.0;
		po = PoleMatrix(mjdutc);
		pp = PrecMatrix(MJD_J2000, mjdtt);
		nn = NutMatrix(mjdtt);
		theta = GHAMatrix(mjdut1);
		ee = po*theta*nn*pp;
		CovEcef[j] = ee * CovEci[j];
		Coverage[j] = Geodetic(CovEcef[j], R_Earth, f_Earth);
	}
}

//------------------------------------------------------------------------------
//
// Maximum computation
//
//------------------------------------------------------------------------------
template<class T> const T& Max(const T& a, const T& b)
{
	return (a<b) ? b : a;
}

//------------------------------------------------------------------------------
//
// 对时间窗口数据进行排序
//
//------------------------------------------------------------------------------
void SortAcess(double ** Access, int WindowNum, int AP)
{
	int row = 6;
        if (AP == 1)
        {
            row = 14;
        }
	//row = sizeof(Access[0]) / sizeof(Access[0][0]); //求解每一行元素个数
	double temp;
	for (int i = 0; i<WindowNum; i++) //对时间进行从大到小的排序（冒泡排序）         
	{
		for (int j = i + 1; j<WindowNum; j++)
		{
			if (Access[j][4]<Access[i][4])
			{
				for (int k = 0; k < row; k++)
				{
					temp = Access[i][k];
					Access[i][k] = Access[j][k];
					Access[j][k] = temp;
				}
			}
		}
	}
}

//------------------------------------------------------------------------------
//
// 球面两点之间距离计算(km), 输入两点的经纬度单位为deg
//
//------------------------------------------------------------------------------

double PointsDistance(double lat1, double lon1, double lat2, double lon2)
{
	lat1 = lat1*Rad;
	lon1 = lon1*Rad;
	lat2 = lat2*Rad;
	lon2 = lon2*Rad;
	return R_Earth * acos(sin(lat1) * sin(lat2) + cos(lat1) * cos(lat2) * cos(lon1 - lon2));
}

//------------------------------------------------------------------------------
//
// 传感器幅宽大小计算（km）,输入角度单位为deg,卫星轨道高度h单位为km
//
//------------------------------------------------------------------------------
// 

double Stripwidth(double Sway, double Fov, double SatH)
{
	double Re, Beta1, Beta2, D1, D2;
	double Sw;
	Re = R_Earth / 1000; //地球半径，km
	Beta2 = (Sway - Fov)*Rad;
	Beta1 = (Sway + Fov)*Rad;
	D2 = asin(sin(Beta2)*(Re + SatH) / Re) - Beta2;
	D1 = asin(sin(Beta1)*(Re + SatH) / Re) - Beta1;
	Sw = Re*(D1 - D2);
	return Sw;
}

//------------------------------------------------------------------------------
//
// 根据给定时刻卫星位置，计算卫星和地面目标夹角（侧摆角），deg
//
//------------------------------------------------------------------------------
// 

double SwayAngle(double lamt, double fait, double lams, double fais, double SatH)
{
	double D, Sway, Angle;
	//lamt = lamt*Rad;
	//lams = lams*Rad;
	//fait = fait*Rad;
	//fais = fais*Rad;
	D = R_Earth * acos(sin(fait) * sin(fais) + cos(fait) * cos(fais) * cos(lamt - lams));
	Angle = D / R_Earth;
	Sway = atan(sin(Angle) / ((R_Earth + SatH) / R_Earth - cos(Angle)))*Deg;
	return Sway;
}

//------------------------------------------------------------------------------
//
// 根据计算得到的目标的时间窗口，通过星地最短距离计算中心成像时刻
//
//------------------------------------------------------------------------------
// 

double CentralTime(Vector SatEcef[], int WindowTime[], Geodetic Target, int* CenTime)
{
	Vector TargetEcef(3);
	Geodetic CenSatLLA;
	int STime, FTime;
	double minDistance, tempDistance, Sway;
	TargetEcef = Target.Position(R_Earth, f_Earth); //计算目标的地固坐标
	STime = WindowTime[0];
	FTime = WindowTime[1];
	*CenTime = STime;
	minDistance = Norm(operator * (1000, SatEcef[STime].slice(0, 2)) - TargetEcef); //将最小距离初始化为开始时刻的星地距离[m]
	for (int i = STime; i <= FTime; i++) //遍历搜索最小星地距离，输出对应时刻
	{
		tempDistance = Norm(operator * (1000, SatEcef[i].slice(0, 2)) - TargetEcef);
		if (tempDistance <= minDistance)
		{
			minDistance = tempDistance;
			*CenTime = i;
		}
	}
	CenSatLLA = Geodetic(operator * (1000, SatEcef[*CenTime].slice(0, 2)), R_Earth, f_Earth); //[m]
	Sway = SwayAngle(CenSatLLA.lon, CenSatLLA.lat, Target.lon, Target.lat, CenSatLLA.h);
	return Sway;
}

//------------------------------------------------------------------------------
//
// The scheduler is optimal in the sense of making maximum usage of the ground task being scheduled.
// The inputs are satellite access times(sat, task, sway angle, start time, stop time, duration)that are ordered by stop time.
// Outputs are the scheduled accesses:sat, task, sway angle, contact start time, contact stop time, contact duration.
//
//------------------------------------------------------------------------------
// 
int GreedySchedule(const int Nsat, const int Ntask, int WindowNum, double ** Access, double ** Sched)
{
	int A = WindowNum; //可见窗口数量
	//初始化参数
	int *rqdCntsPerDay = new int[Nsat]; //定义卫星一天的最大开机次数
	for (int i = 0; i<Nsat; i++)
		rqdCntsPerDay[i] = 100;

	int *rqdCntsPerTask = new int[Ntask];//定义任务每天的需求次数；若为周期重复观测任务, 可以增加一个实际单任务观测次数
	for (int i = 0; i<Ntask; i++)
		rqdCntsPerTask[i] = 1;

	int *actCntsPerDay = new int[Nsat]; //每颗星每天实际的观测目标次数
	for (int i = 0; i < Nsat; i++)
		actCntsPerDay[i] = 0;

	int *actCntsPerTask = new int[Ntask]; //实际每个任务的被观测次数	
	for (int i = 0; i < Ntask; i++)
		actCntsPerTask[i] = 0;

	int *sattime = new int[Nsat]; //Satellite loading time
	for (int i = 0; i < Nsat; i++)
		sattime[i] = 0;

	int *tasktime = new int[Ntask]; //task loading time
	for (int i = 0; i < Ntask; i++)
		tasktime[i] = 0;

	//Initialize tsoplast
	int tstoplast = 0;
	int schedsat = -1; //用于判断是否连续两个任务由同一颗卫星调度
	int *schedtstoplast = new int[Nsat]; //the last scheduled contact stop time for each task
	for (int i = 0; i < Nsat; i++)
		schedtstoplast[i] = 0;
	int numsched = 0; //Total Number of scheduled accesses
	double setuptime = 2 * 60; //卫星设置时间[s]
	double durmin = 0.2 * 60; //最小持续时长[s]
	double durmax = 0.4 * 60; //最大持续时长[s]
	//动态初始化调度结果为最大的窗口个数A
	int *task = new int[A];
	int *sat = new int[A];
	double *sway = new double[A];
	double *tstart = new double[A];
	double *tstop = new double[A];
	double *dur = new double[A];
	int N, M; //cuurent sat and task;
	double setupmargin; //同一卫星连续两次成像时间间隔
	//Loop through all accesses
	for (int i = 0; i < A; i++)
	{
		task[i] = Access[i][0];
		sat[i] = Access[i][1];
		sway[i] = Access[i][2];
		tstart[i] = Access[i][3];
		tstop[i] = Access[i][4];
		dur[i] = Access[i][5];
		N = sat[i];
		M = task[i];
		//If we’re trying to schedule on the same sat then we need the margin of setup time; Otherwise there is no margin
		if (schedsat == sat[i])
			setupmargin = setuptime;
		else
			setupmargin = 0;
		//Get the last contact’s stop time that’s been scheduled on the current task.
		tstoplast = schedtstoplast[N]; //当前卫星最后一个接触的停止时间
		//if we haven’t reached the required contacts per day for the current %sat;
		if (actCntsPerDay[N] < rqdCntsPerDay[N] && actCntsPerTask[M] < rqdCntsPerTask[M]) //判断实际接触次数是否满足需求次数
		{
			//Test for the opportunities start time to be later than the previous stop time plus the setup time;
			//Also the duration must be greater than the minimum duration allowed; duration is in minutes;
			if (dur[i] >= durmin) //断持续时间是否充足
			{
				if (tstart[i] >= tstoplast + setupmargin) //& (dur(i) >= durmin)
				{
					//cout << "可见窗口" << i << "来自任务" << task[i] << "和卫星" << N << "，可以被调度" << endl; //输出被调度的窗口
					//----------Set scheduled accesses----------
					Sched[numsched][0] = task[i];
					Sched[numsched][1] = sat[i];
					Sched[numsched][2] = sway[i];
					Sched[numsched][3] = tstart[i];
					//limit the contact duration to the user set maximim duration，约束调度时长为窗口持续时间和最大时间的小者
					if (dur[i] >= durmax)
						Sched[numsched][5] = durmax;
					else
						Sched[numsched][5] = dur[i];
					//new stop time is now start time + actual contact duration;
					Sched[numsched][4] = tstart[i] + Sched[numsched][5]; //Keep track of the sats and tasks that have now been scheduled
					schedsat = sat[i];
					//Keep track of the end of the last scheduled contact;
					schedtstoplast[N] = tstop[i];
					//----------Increment data counters----------增量数据计数器
					actCntsPerDay[N] = actCntsPerDay[N] + 1; //每颗星的每天实际的观测次数
					actCntsPerTask[M] = actCntsPerTask[M] + 1; //每个任务每天实际的观测次数
					sattime[N] = sattime[N] + Sched[numsched][5]; //Satellite loading time(minutes)，卫星观测总时间
					tasktime[M] = tasktime[M] + Sched[numsched][5]; //task loading time(minutes)，任务被观测总时间
					numsched = numsched + 1; //已经调度任务计数器
				}
				else
				{
					//同一卫星连续两次任务时间窗口冲突,不能被调度
					//cout << "任务" << task[i] << "和卫星" << N <<"的可见窗口" << i << "由于和前一任务时间冲突，故不能被调度" << endl;
				}
			}
			else
			{
				//The duration of the pass is insufficient to meet the miminum required contact duration;
				//cout << "任务" << task[i] << "和卫星" << N << "的可见窗口" << i << "由于时长太短，故不能被调度" << endl;
			}
		}
		else
		{
			//number of required contacts for this sat has been met，接触数量满足.find sat name corresponding to current satellite number;
			//cout << "任务" << task[i] << "和卫星" << N << "的可见窗口" << i << "由于该任务已经被调度过，故该窗口被舍弃，不能被调度" << endl;
		}
	}
	return numsched;
}


//------------------------------------------------------------------------------
//
// Main program
//
//------------------------------------------------------------------------------
int main() {

	cout << "\n      Satellite Mission Plan Based on High Precision Orbit Propagator     \n" << endl;

	// 定义变量
	double    Mjd_UTC;
	Vector    Y0(6), Y(6), ECIr(3), ECIv(3), ECEFr(3), ECEFv(3);
	Vector    Kep(6);
	AuxParam  Aux;     // Auxiliary parameters
	int Year, Month, Day, Hour, Min;
	double Sec;
	// 设置外推轨道时长与步长，结果存储变量
	const int Nsat = 120;
	const int Ntask = 10;
	const int AreaNtask = 1;
	double Step = 1; // [s]
	const int N_Step = 3600*3; // 3h
	static Vector Eph[N_Step + 1]; //
	static Vector EphAll[Nsat][N_Step + 1];
	static Vector SatEcef[Nsat][N_Step + 1];
	static Vector SatEci[Nsat][N_Step + 1];
	static Vector mycoe[Nsat][N_Step + 1];
	// 转换矩阵变量定义
	double mjdtt, mjdut1, mjdutc;
	Matrix pp(3, 3), nn(3, 3), theta(3, 3), dtheta(3, 3), ee(3, 3), dee(3, 3), po(3, 3); //eci2ecef的矩阵
	Matrix ss(3, 3);
	Vector vec0(0, -1, 0);
	Vector vec1(1, 0, 0);
	Vector vec2(0, 0, 0);
	ss.SetCol(0, vec0);
	ss.SetCol(1, vec1);
	ss.SetCol(2, vec2);
	// 打开GGM03S.txt文件读取地球重力场参数
	ifstream inp;
	inp.open("GGM03S.txt");
	int z = 0, n = 360;
	double temp;
	do
	{
		for (int x = 0; x <= z; x++)
		{
			inp >> temp;
			inp >> temp;
			inp >> temp;
			cnm(z, x) = temp;
			inp >> temp;
			snm(z, x) = temp;
			inp >> temp;
			inp >> temp;
		}  z++;
	} while (z <= n);
	inp.close();

	//读取轨道初始状态、卫星参数、历元时间（Envisat）
	int Ephyear, Ephmonth, Ephday, Ephhour, Ephmin;
	double Ephsec, Areadrag, Areasolar, Mass;
	inp.open("InitialState_120Sats.txt");
	FILE *f;
	if ((f = fopen("SatelliteStatesAll.txt", "w+")) == NULL)
	{
		fprintf(stdin, "Can't open \"words\" file.\n");
		exit(1);
	}
	// 循环对120颗卫星预报3小时轨道
	for (int ii = 0; ii < Nsat; ii++)
	{
		inp >> temp; Ephyear = temp;
		inp >> temp; Ephmonth = temp;
		inp >> temp; Ephday = temp;
		inp >> temp; Ephhour = temp;
		inp >> temp; Ephmin = temp;
		inp >> temp; Ephsec = temp;
		inp >> temp; Areadrag = temp;
		inp >> temp; Areasolar = temp;
		inp >> temp; Mass = temp;
		for (int j = 0; j < 6; j++)
		{
			inp >> temp;
			Y0(j) = 1e3*temp;
		}

		// 设置初始历元状态 (Envisat)
		Mjd_UTC = Mjd(Ephyear, Ephmonth, Ephday, Ephhour, Ephmin, Ephsec);
		jd = Mjd_UTC + 2400000.5;
		mfme = 1440.0*(Mjd_UTC - floor(Mjd_UTC));
		initeop(eoparr, jdeopstart);
		initspw(spwarr, jdspwstart);
		findeopparam(jd, mfme, interp, eoparr, jdeopstart, dut1, dat, lod, xp, yp, ddpsi, ddeps, dx, dy, x, y, s, deltapsi, deltaeps);
		// Initialize UT1-UTC and UTC-TAI time difference
		IERS::Set(dut1, -dat, xp, yp);

		// 设置摄动力模型参数
		Aux.Mjd0_TT = Mjd_UTC + IERS::TT_UTC(Mjd_UTC) / 86400.0;
		Aux.Area_drag = Areadrag;  // [m^2]
		Aux.Area_solar = Areasolar;   // [m^2]
		Aux.mass = Mass; // [kg]
		Aux.CR = 1.0;
		Aux.CD = 2.7;
		Aux.n = 36;
		Aux.m = 36;
		Aux.Sun = true;
		Aux.Moon = true;
		Aux.SRad = true;
		Aux.Drag = true;

		// 调用Ephemeris函数计算轨道
		Ephemeris(Y0, N_Step, Step, Aux, EphAll[ii]);
		fprintf(f, "%4d", ii); //输出卫星标号
		fprintf(f, "\n");
		// 将Eph中存储的J2000轨道转换成地固系下的轨道坐标
		for (int i = 0; i <= N_Step; i++)
		{
			Y = 1e-3*EphAll[ii][i]; //j2000坐标，km
			ECIr = Y.slice(0, 2); //将eci位置和速度切片分离
			ECIv = Y.slice(3, 5);
			SatEci[ii][i] = Stack(ECIr, ECIv); //卫星ecef地固坐标
			mycoe[ii][i] = COE(ECIr, ECIv);
			//计算ecef坐标
			mjdtt = Mjd_UTC + (Step*i) / 86400.0 + IERS::TT_UTC(Mjd_UTC + (Step*i) / 86400.0) / 86400.0;
			mjdutc = Mjd_UTC + (Step*i) / 86400.0;
			mjdut1 = Mjd_UTC + (Step*i) / 86400.0 + IERS::UT1_UTC(Mjd_UTC + (Step*i) / 86400.0) / 86400.0;
			po = PoleMatrix(mjdutc);
			pp = PrecMatrix(MJD_J2000, mjdtt);
			nn = NutMatrix(mjdtt);
			theta = GHAMatrix(mjdut1);
			dtheta = (omega_Earth - 0.843994809*1e-9*lod)*ss*theta;
			ee = po*theta*nn*pp;
			ECEFr = ee*ECIr;
			dee = po*dtheta*nn*pp;
			ECEFv = ee*ECIv + dee*ECIr;
			SatEcef[ii][i] = Stack(ECEFr, ECEFv); //卫星ecef地固坐标
			CalDat((Mjd_UTC + (Step*i) / 86400.0), Year, Month, Day, Hour, Min, Sec);

			fprintf(f, "%4d/%02d/%02d-", Year, Month, Day);
			fprintf(f, "%02d:%02d:%06.3f", Hour, Min, Sec);

			for (int j = 0; j < 3; j++)
			{
				fprintf(f, "%15.6f", ECIr(j));
			}
			for (int j = 0; j < 3; j++)
			{
				fprintf(f, "%15.9f", ECIv(j));
			}
			fprintf(f, "\n");
		};
		cout << "Orbit prediction for Satellite " << ii << " is done." << endl;
	}
	inp.close();
	fclose(f);

	// 读取任务信息：地理位置、时间约束、分辨率、成像类型
	double MaxFov[2] = { 60, 5};
	double RealFov[2] = { 5, 5 };
	double Sway = 0.0;
	Geodetic TarArea[4]; //定义区域目标的经纬高
	vector<Point>TargetArea(4);
	Geodetic TarPoint; //定义点目标的经纬高
	const int MaxNumOrbit = 25; //最大轨道圈数
	// 定义输出窗口指针AcessTime
	double ** AccessTime;
	AccessTime = new double*[MaxNumOrbit];
	for (int i = 0; i < MaxNumOrbit; i++)
		AccessTime[i] = new double[16];
	for (int i = 0; i < MaxNumOrbit; i++)
		for (int j = 0; j < 16; j++)
			AccessTime[i][j] = 0;
	int TotalNumAll = 0; //所有卫星对所有任务的时间窗口数量
	int TotalNum[Nsat];
	for (int i =0;i<Nsat;i++)
		TotalNum[i] = 0; //单颗卫星对所有任务的总时间窗口数量
	int SatName, TaskName;
	Geodetic ** CovBoundary; //定义条带的二维数组
        int WindowNum[Ntask];
        for (int i =0;i<Ntask;i++)
		WindowNum[i] = 0; //单颗卫星对单个任务的时间窗口数量，
        inp.open("Task.txt");
	int AP; //判断是区域目标/点目标
	inp >> temp; AP = temp;
        inp.close();
        if (AP == 1)
        {
            // 将区域目标时间窗口和侧摆角度数据存储在AccessAreaAll文件中
	    if ((f = fopen("AccessAreaAll.txt", "w+")) == NULL)
	    {
	        fprintf(stdin, "Can't open \"words\" file.\n");
		exit(1);
	    }
        }
	else
        {
            // 将点目标时间窗口和侧摆角度数据存储在AccessPointAll文件中
	    if ((f = fopen("AccessPointAll.txt", "w+")) == NULL)
	    {
		fprintf(stdin, "Can't open \"words\" file.\n");
		exit(1);
	    }
        }
        // 计算时间窗口
	for (int ii = 0; ii < Nsat; ii++)
	{
		inp.open("Task.txt");
		int species; //判断是区域目标/点目标
		inp >> temp; species = temp;
		if (species == 1)
		{
			for (int i = 0; i < AreaNtask; i++)
			{
				for (int k = 0; k < 4; k++)
				{
					// Task三个数据分别对应的：纬度、经度、高度
					inp >> temp; TarArea[k].lat = temp* Rad; TargetArea[k].y = temp;
					inp >> temp; TarArea[k].lon = temp* Rad; TargetArea[k].x = temp;
					inp >> temp; TarArea[k].h = temp;
				}
				WindowNum[i] = AccessWindowsArea(SatEci[ii], SatEcef[ii], TarArea, mycoe[ii], MaxFov, N_Step, MaxNumOrbit, Sway, Mjd_UTC, AccessTime);  //计算区域目标的时间窗口
				//计算传感器覆盖边界
				int ** sftime; //初始化开始结束时间下标
				sftime = new int*[WindowNum[i]];
				for (int j = 0; j < WindowNum[i]; j++)
					sftime[i] = new int[2];
				for (int j = 0; j < WindowNum[i]; j++)
					for (int k = 0; k < 2; k++)
						sftime[j][k] = static_cast<int>(AccessTime[j][12 + k]); //将AcessTime中时间下标赋值给sftime
				//初始化条带覆盖边界
				CovBoundary = new Geodetic*[WindowNum[i]];
				for (int j = 0; j < WindowNum[i]; j++)
					CovBoundary[j] = new Geodetic[4];
				for (int j = 0; j < WindowNum[i]; j++)
					for (int k = 0; k < 4; k++)
					{
						CovBoundary[j][k].lat = 0;
						CovBoundary[j][k].lon = 0;
						CovBoundary[j][k].h = 0;
					}
				// 调用CalculateCov函数计算条带覆盖边界
				for (int j = 0; j < WindowNum[i]; j++)
				{
					CalculateCov(SatEci[ii], mycoe[ii], sftime[j], RealFov, AccessTime[j][15], Mjd_UTC, CovBoundary[j]);
					//for (int k = 0; k < 4; k++) //输出条带坐标
					//{
						//cout << CovBoundary[j][k].lat << "," << CovBoundary[j][k].lon << "," << CovBoundary[j][k].h << endl;
					//}
				}
				SatName = ii; //卫星标号
				TaskName = i; //任务标号
				for (int j = 0; j < WindowNum[i]; j++)
				{
					fprintf(f, "%4d", TaskName);
					fprintf(f, "%4d", SatName);
					fprintf(f, "%9.3f ", AccessTime[j][15]);
					fprintf(f, "%6.f ", AccessTime[j][12]);
					fprintf(f, "%6.f ", AccessTime[j][13]);
					fprintf(f, "%6.f ", AccessTime[j][13] - AccessTime[j][12]);
					for (int k = 0; k < 4; k++) //循环输出条带四个顶点坐标
					{
						fprintf(f, "%9.3f ", CovBoundary[j][k].lat);
						fprintf(f, "%9.3f ", CovBoundary[j][k].lon);
					}
					
					fprintf(f, "\n");
				}
				TotalNum[ii] += WindowNum[i];
			}
			
		}
		else
		{
			//fprintf(f, "Task Sat Angle Stime Ftime Dtime");
			//fprintf(f, "\n");
			for (int i = 0; i < Ntask; i++)
			{
				inp >> temp; TarPoint.lat = temp * Rad;
				inp >> temp; TarPoint.lon = temp * Rad;
				inp >> temp; TarPoint.h = 0;
				WindowNum[i] = AccessWindowsPoint(SatEci[ii], SatEcef[ii], TarPoint, mycoe[ii], MaxFov, N_Step, MaxNumOrbit, Sway, Mjd_UTC, AccessTime);  //计算点目标的时间窗口
				SatName = ii;
				TaskName = i;
				for (int j = 0; j < WindowNum[i]; j++)
				{
					fprintf(f, "%4d", TaskName);
					fprintf(f, "%4d", SatName);
					fprintf(f, "%9.3f ", AccessTime[j][15]);
					fprintf(f, "%6.f ", AccessTime[j][12]);
					fprintf(f, "%6.f ", AccessTime[j][13]);
					fprintf(f, "%6.f ", AccessTime[j][13] - AccessTime[j][12]);
					fprintf(f, "\n");
				}
				TotalNum[ii] += WindowNum[i];
			}
		}
		inp.close();
		TotalNumAll += TotalNum[ii];
		cout << "TimeWindow calculation for Satellite:" << ii << " is done." << endl;
	}
	fclose(f);
        
        // 定义调度时间窗口变量以及Access变量
        int SchedSYear, SchedSMonth, SchedSDay, SchedSHour, SchedSMin, SchedFYear, SchedFMonth, SchedFDay, SchedFHour, SchedFMin;
	double SchedSSec, SchedFSec;
	int SchedNum = 0;
	double ** Access;
	Access = new double*[TotalNumAll];
        // 根据时间窗口分配任务
        if (AP == Ntask) //多个点目标任务
        {
			//输出点目标规划的结果
			if ((f = fopen("result_spot.txt", "w+")) == NULL)
			{
				fprintf(stdin, "Can't open \"words\" file.\n");
				exit(1);
			}
			fprintf(f, "task  sat  slewangle    starttime      finishtime");
			fprintf(f, "\n");
            cout << "PointTargets plan:" << endl;
	    for (int i = 0; i < TotalNumAll; i++)
		Access[i] = new double[6];
	    inp.open("AccessPointAll.txt"); //读取Access数据
            cout << "All the TimeWindow data:" << endl;
            cout << "TaskID " << "SatID " << "Swayangle " << "StartTime" << "EndTime " << "Duration " << endl;
	    for (int i = 0; i < TotalNumAll; i++)
	    {
		for (int j = 0; j < 6; j++)
		{
			inp >> temp;
			Access[i][j] = temp;
			cout << Access[i][j] << " ";
		}
		cout << endl;
	    }
	    inp.close();
	    SortAcess(Access, TotalNumAll, AP); // 对时间窗口数据按结束时间最小进行由小到大排序
	    double ** Sched;
	    Sched = new double*[TotalNumAll];
	    for (int i = 0; i < TotalNumAll; i++)
		Sched[i] = new double[6];
	    SchedNum = GreedySchedule(Nsat, Ntask, TotalNumAll, Access, Sched);
	    for (int i = 0; i < SchedNum; i++)
	    {
		CalDat((Mjd_UTC + (Step*Sched[i][3]) / 86400.0), SchedSYear, SchedSMonth, SchedSDay, SchedSHour, SchedSMin, SchedSSec);
		CalDat((Mjd_UTC + (Step*Sched[i][4]) / 86400.0), SchedFYear, SchedFMonth, SchedFDay, SchedFHour, SchedFMin, SchedFSec);
		cout << "Task " << Sched[i][0] << " is executed by satellite " << Sched[i][1] << ",and its swayangle is " << Sched[i][2];
		cout << ".StratTime is ";
		printf("%4d/%02d/%02d-", SchedSYear, SchedSMonth, SchedSDay);
		printf("%02d:%02d:%06.3f", SchedSHour, SchedSMin, SchedSSec);
		cout << ",EndTime is ";
		printf("%4d/%02d/%02d-", SchedFYear, SchedFMonth, SchedFDay);
		printf("%02d:%02d:%06.3f", SchedFHour, SchedFMin, SchedFSec);
		cout << "." << endl;
		// 输出最终结果到result
		fprintf(f, "%4d", (int)Sched[i][0]); //输出任务、卫星、角度
                fprintf(f, " ");
                fprintf(f, "%4d", (int)Sched[i][1]);
                fprintf(f, "    ");
                fprintf(f, "%6.3f", Sched[i][2]);
                fprintf(f, "    ");
		fprintf(f, "%4d/%02d/%02d-", SchedSYear, SchedSMonth, SchedSDay);
		fprintf(f, "%02d:%02d:%06.3f", SchedSHour, SchedSMin, SchedSSec);
                fprintf(f, " ");
		fprintf(f, "%4d/%02d/%02d-", SchedFYear, SchedFMonth, SchedFDay);
		fprintf(f, "%02d:%02d:%06.3f", SchedFHour, SchedFMin, SchedFSec);
		fprintf(f, "\n");
	    }
            fclose(f);
        }
		

        if (AP == AreaNtask) //单个区域目标任务
{
	//输出区域目标规划的结果
	if ((f = fopen("result_area.txt", "w+")) == NULL)
	{
		fprintf(stdin, "Can't open \"words\" file.\n");
		exit(1);
	}
	
        cout << "AreaTargets plan:" << endl;
	//计算备选条带总覆盖面积
	int StripNum, StripNumMod;
	int AccessTotalNumAllMod = 0;
	double CovArea, AreaArea, CovRatio;
	AreaArea = CalculateAreaNew(TargetArea); //计算观测区域面积
	vector<vector<Point>> StripSet(TotalNumAll);
	vector<vector<Point>> StripSetIn(TotalNumAll);
	vector<vector<Point>> StripSetMod; //去掉观测面积较小的条带
	vector<vector<Point>> StripSetPar; //部分条带
	for (int i = 0; i < StripSet.size(); i++)
		StripSet[i].resize(4);
	double ** AccessMod; //去掉观测面积较小的时间窗口
	for (int i = 0; i < TotalNumAll; i++)
		Access[i] = new double[14];
	AccessMod = new double*[TotalNumAll];
	for (int i = 0; i < TotalNumAll; i++)
		AccessMod[i] = new double[14];
	inp.open("AccessAreaAll.txt"); //读取Access数据
        cout << "All the TimeWindow data:" << endl;
        cout << "TaskID " << "SatID " << "Swayangle " << "StartTime" << "EndTime " << "Duration " << "CoveargeBoundary" << endl;
	for (int i = 0; i < TotalNumAll; i++)
	{
		for (int j = 0; j < 14; j++)
		{
			inp >> temp;
			Access[i][j] = temp;
			cout << Access[i][j] << " ";
		}
		cout << endl;
	}
	inp.close();
	SortAcess(Access, TotalNumAll, AP); //根据结束时间最早进行排序
	//for (int i = 0; i < TotalNumAll; i++)
	//{
		//cout << Access[i][1] << endl;
	//}
	for (int i = 0; i < TotalNumAll; i++)
		for (int j = 0; j < 8; j=j+2)
		{
			StripSet[i][j/2].x = Access[i][j + 7]; //lon
			StripSet[i][j/2].y = Access[i][j + 6]; //lat
		}
	StripNum = StripSet.size();
	for (int i = 0; i < StripNum; i++)
	{
		vector<Point> tempIn; //在循环里初始化，避免坐标累计
		StripSetIn[i] = PolygonClip(StripSet[i], TargetArea, tempIn); //求每个条带对于观测区域覆盖的面积
		if (CalculateAreaNew(StripSetIn[i]) > 0.1*AreaArea) //条带面积小于观测区域的10%时，放弃该次观测
		{
			StripSetMod.push_back(StripSetIn[i]);
			for (int k = 0; k < 14; k++)
				AccessMod[AccessTotalNumAllMod][k] = Access[i][k];
			AccessTotalNumAllMod++;
				
		}	
	}
	StripNumMod = StripSetMod.size();
	for (int i = 0; i < StripNumMod; i++)
	{
		StripSetPar.push_back(StripSetMod[i]);
		CovArea = StripTotalCoverage(StripSetPar, i+1);
		if (CovArea >= 0.8*AreaArea) //当条带集合覆盖总面积大于观测区域面积的80%时，认为区域被观测完成
		{
			SchedNum = i + 1; //记录使用的时间窗口数
			break;
		}
		else
			SchedNum++;
	}
	if (CovArea > AreaArea)
		CovRatio = AreaArea;
	CovRatio = CovArea / AreaArea;
	cout << endl;
	cout << "Total coverarge area is " << CovArea << "km^2,and the coverage ratio is " << CovRatio << "." << endl;
	fprintf(f, "totalcoverarea  totalcoverratio");
	fprintf(f, "\n");
	fprintf(f, "%10.3f", CovArea);
	fprintf(f, "      ");
	fprintf(f, "%6.5f", CovRatio);
	fprintf(f, "\n");
	fprintf(f, "task  sat  slewangle    starttime      finishtime");
        fprintf(f, "\n");
	for (int i = 0; i < SchedNum; i++)
	{
		CalDat((Mjd_UTC + (Step*AccessMod[i][3]) / 86400.0), SchedSYear, SchedSMonth, SchedSDay, SchedSHour, SchedSMin, SchedSSec);
		CalDat((Mjd_UTC + (Step*AccessMod[i][4]) / 86400.0), SchedFYear, SchedFMonth, SchedFDay, SchedFHour, SchedFMin, SchedFSec);
		cout << "Task " << AccessMod[i][0] << " is executed by satellite " << AccessMod[i][1] << ",and its swayangle is " << AccessMod[i][2];
		cout << ".StratTime is ";
		printf("%4d/%02d/%02d-", SchedSYear, SchedSMonth, SchedSDay);
		printf("%02d:%02d:%06.3f", SchedSHour, SchedSMin, SchedSSec);
		cout << ",EndTime is ";
		printf("%4d/%02d/%02d-", SchedFYear, SchedFMonth, SchedFDay);
		printf("%02d:%02d:%06.3f", SchedFHour, SchedFMin, SchedFSec);
		cout << "." << endl;
		// 输出最终结果到result
		fprintf(f, "%4d", (int)AccessMod[i][0]); //输出任务、卫星、角度
                fprintf(f, " ");
                fprintf(f, "%4d", (int)AccessMod[i][1]);
                fprintf(f, "    ");
                fprintf(f, "%6.3f", AccessMod[i][2]);
                fprintf(f, "   ");
		fprintf(f, "%4d/%02d/%02d-", SchedSYear, SchedSMonth, SchedSDay);
		fprintf(f, "%02d:%02d:%06.3f", SchedSHour, SchedSMin, SchedSSec);
                fprintf(f, " ");
		fprintf(f, "%4d/%02d/%02d-", SchedFYear, SchedFMonth, SchedFDay);
		fprintf(f, "%02d:%02d:%06.3f", SchedFHour, SchedFMin, SchedFSec);
		fprintf(f, "\n");
	}
        fclose(f);
}
        
	printf("\nPlan is finished,press any key to exit \n");
	getchar();
	return 0;
}
