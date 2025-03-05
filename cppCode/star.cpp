////////////////////////////////////
#include "star.h"

#include <iostream>
#include <cmath>

#include <TMatrix.h>
using namespace TRiOLD;

////////////////////////////////////
Star::~Star()
{
}

Star::Star()
    : m_GCC(), m_GCV(), m_addData_ptr(nullptr)
{
}

Star::Star(const Equatorial& EqC, const Equatorial& EqV)
    : Star(_calcGC(EqC), _calcGV(EqC, EqV))
{
}

Star::Star(const Galactic& GC, const Galactic& GV)
    : m_GCC(_calcGCC(GC)), m_GCV(_calcGCV(GC, GV)), m_addData_ptr(nullptr)
{
}

Star::Star(const Cartesian& GCC, const Cartesian& GCV)
    : m_GCC(GCC), m_GCV(GCV), m_addData_ptr(nullptr)
{
}

Star::Cartesian Star::getGCC() const
{
    return m_GCC;
}

Star::Cartesian Star::getGCV() const
{
    return m_GCV;
}

void* Star::getAddData_ptr() const
{
    return m_addData_ptr;
}

void Star::setAddData_ptr(void* addData_ptr)
{
    m_addData_ptr = addData_ptr;
}

//void Star::getcalcEqCandV(Equatorial& EqC, Equatorial& EqV) const
//{
//}

void Star::getcalcGCandV(Galactic& GC, Galactic& GV) const
{
    GC = _calcGC(m_GCC);
    GV = _calcGV(GC, m_GCV);
}

void Star::getcalcGCCCandV(Cylindrical& GCCC, Cylindrical& GCCV, double R_Sun, const Cartesian& GCV_Sun) const
{
    GCCC.R = std::sqrt(std::pow(m_GCC.x-R_Sun, 2.0) + std::pow(m_GCC.y, 2.0));
    if(GCCC.R == 0.0)
        GCCC.theta = M_PI;
    else GCCC.theta = std::atan2(m_GCC.y, m_GCC.x-R_Sun);
    if(GCCC.theta < 0.0)
        GCCC.theta += 2*M_PI;
    GCCC.Z = m_GCC.z;

    Cartesian GCV_abs;
    GCV_abs.x = -(m_GCV.x + GCV_Sun.x);
    GCV_abs.y = -(m_GCV.y + GCV_Sun.y);
    GCV_abs.z =   m_GCV.z + GCV_Sun.z;
    double rotAngleXY =_calcRotAngleXY(R_Sun, m_GCC);
    GCCV.R = GCV_abs.x * std::cos(rotAngleXY) - GCV_abs.y * std::sin(rotAngleXY);
    GCCV.theta = GCV_abs.x * std::sin(rotAngleXY) + GCV_abs.y * std::cos(rotAngleXY);
    GCCV.Z = GCV_abs.z;
}

Star::Cartesian Star::getcalcGCC_local(const Cartesian& GCC0) const
{
    return Cartesian( m_GCC.x-GCC0.x, m_GCC.y-GCC0.y, m_GCC.z-GCC0.z);
}

void Star::getcalcCCandV_localRot(Cartesian& GCC_localRot, Cartesian& GCV_localRot,
                                  const Cartesian& GCC0, double R_Sun) const
{
    Cartesian GCC_local = getcalcGCC_local(GCC0);
    double rotAngleXY = _calcRotAngleXY(R_Sun, GCC_local);
    GCC_localRot = _rotateVector(GCC_local, rotAngleXY);
    GCV_localRot = _rotateVector(m_GCV, rotAngleXY);
}

Star::Galactic Star::_calcGC(const Equatorial& EqC)
{
    //Ag tensor
    Matrix<double> Ag(
        -0.0548755601367195, -0.8734370902532698, -0.4838350155472244,
         0.4941094280132430, -0.4448296298016944,  0.7469822445004389,
        -0.8676661489582886, -0.1980763737065720,  0.4559837761713720);

    Matrix<double> q(3,1);
    q[0][0] = std::cos(EqC.ra) * std::cos(EqC.de);
    q[1][0] = std::sin(EqC.ra) * std::cos(EqC.de);
    q[2][0] = std::sin(EqC.de);

    Matrix<double> g = Ag * q;
    Galactic res;
    res.l = std::atan2(g[1][0], g[0][0]);
    res.b = std::atan2(g[2][0], std::sqrt(g[0][0]*g[0][0] + g[1][0]*g[1][0]));
    res.r = EqC.r;
    return res;
}

Star::Galactic Star::_calcGV(const Equatorial& EqC, const Equatorial& EqV)
{
    const double ra_gp = 3.36603291968;
    const double de_gp = 0.47347728280;
    double c1 = std::sin(de_gp)*std::cos(EqC.de)-std::cos(de_gp)*std::sin(EqC.de)*std::cos(EqC.ra-ra_gp);
    double c2 = std::cos(de_gp) * std::sin(EqC.ra - ra_gp);

    Galactic res;
    res.l = (c1*EqV.ra + c2*EqV.de) / std::sqrt(c1*c1 + c2*c2);
    res.b = (c1*EqV.de - c2*EqV.ra) / std::sqrt(c1*c1 + c2*c2);
    res.r = EqV.r;
    return res;
}

Star::Galactic Star::_calcGC(const Cartesian& GCC)
{
    Galactic res;
    res.r = std::sqrt(GCC.x*GCC.x + GCC.y*GCC.y + GCC.z*GCC.z);
    res.l = std::atan2(GCC.y, GCC.x);
    res.b  = std::asin(GCC.z / res.r);
    return res;
}

Star::Galactic Star::_calcGV(const Galactic& GC, const Cartesian& GCV)
{
    const double k = 4.74057;
    Matrix<double> R(
         std::cos(GC.l)*std::cos(GC.b), std::sin(GC.l)*std::cos(GC.b), std::sin(GC.b),
        -std::sin(GC.l), std::cos(GC.l),  0.0,
        -std::cos(GC.l)*std::sin(GC.b), -std::sin(GC.l)*std::sin(GC.b), std::cos(GC.b));


    Galactic res;
    res.r =  R[0][0]*GCV.x + R[0][1]*GCV.y + R[0][2]*GCV.z;
    res.l = (R[1][0]*GCV.x + R[1][1]*GCV.y + R[1][2]*GCV.z) / k / GC.r;
    res.b = (R[2][0]*GCV.x + R[2][1]*GCV.y + R[2][2]*GCV.z) / k / GC.r;
    return res;
}

Star::Cartesian Star::_calcGCC(const Galactic& GC)
{
    Cartesian res;
    res.x = GC.r * std::cos(GC.l) * std::cos(GC.b);
    res.y = GC.r * std::sin(GC.l) * std::cos(GC.b);
    res.z = GC.r * std::sin(GC.b);
    return res;
}

Star::Cartesian Star::_calcGCV(const Galactic& GC, const Galactic& GV)
{
    const double k = 4.74057;
    Matrix<double> R(std::cos(GC.l) * std::cos(GC.b), -std::sin(GC.l), -std::cos(GC.l) * std::sin(GC.b),
                     std::sin(GC.l) * std::cos(GC.b), std::cos(GC.l), -std::sin(GC.l) * std::sin(GC.b),
                     std::sin(GC.b), 0.0, std::cos(GC.b));

    Cartesian res;
    res.x = R[0][0] * GV.r + R[0][1] * k*GC.r * GV.l + R[0][2] * k*GC.r * GV.b;
    res.y = R[1][0] * GV.r + R[1][1] * k*GC.r * GV.l + R[1][2] * k*GC.r * GV.b;
    res.z = R[2][0] * GV.r + R[2][1] * k*GC.r * GV.l + R[2][2] * k*GC.r * GV.b;
    return res;
}

double Star::_calcRotAngleXY(double R_Sun, const Cartesian& GCC)
{
    return -(M_PI + std::atan2(GCC.y, GCC.x-R_Sun));
}

Star::Cartesian Star::_rotateVector(const Cartesian& vector, double angleXY)
{
    Cartesian res;
    res.x = vector.x*std::cos(angleXY) - vector.y*std::sin(angleXY);
    res.y = vector.x*std::sin(angleXY) + vector.y*std::cos(angleXY);
    res.z = vector.z;
    return res;
}

////////////////////////////////////
