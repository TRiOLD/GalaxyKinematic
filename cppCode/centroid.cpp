////////////////////////////////////
#include "centroid.h"

#include <cmath>

////////////////////////////////////
Centroid::Centroid()
    : Star(), m_starsAmount(0), m_V0(), m_Mm(3), m_Mp(3), m_KPsErr()
{

}

Centroid::Centroid(const Cartesian &GCC, const Cartesian &GCV, std::size_t starsAmount,
                   const Cartesian &V0, const Matrix<double> &M)
    : Star(GCC, GCV), m_starsAmount(starsAmount), m_V0(V0), m_KPsErr()
{
    _expandTensor(m_Mm, m_Mp, M);
}

Centroid::Centroid(const Cartesian &GCC, const Cartesian &GCV, std::size_t starsAmount,
                   const Cartesian &V0, const Matrix<double> &Mm, const Matrix<double> &Mp)
    : Star(GCC, GCV), m_starsAmount(starsAmount), m_V0(V0), m_Mm(Mm), m_Mp(Mp), m_KPsErr()
{

}

Centroid::Centroid(const Cartesian &GCC, const Cartesian &GCV, std::size_t starsAmount,
                   const KinematicParameters &KPs)
    : Star(GCC, GCV), m_starsAmount(starsAmount), m_KPsErr()
{
    _expandKPs(m_V0, m_Mm, m_Mp, KPs);
}

std::size_t Centroid::getStarsAmount() const
{
    return m_starsAmount;
}

Star::Cartesian Centroid::getV0() const
{
    return m_V0;
}

Matrix<double> Centroid::getMm() const
{
    return m_Mm;
}

Matrix<double> Centroid::getMp() const
{
    return m_Mp;
}

Matrix<double> Centroid::getcalcM() const
{
    return m_Mm + m_Mp;
}

Star::Cartesian Centroid::getcalcV0_localRot(double R_Sun) const
{
    return _rotateVector(m_V0, _calcRotAngleXY(R_Sun, m_GCC));
}

Matrix<double> Centroid::getcalcMm_localRot(double R_Sun) const
{
    return _rotateMatrix(m_Mm, _calcRotAngleXY(R_Sun, m_GCC));
}

Matrix<double> Centroid::getcalcMp_localRot(double R_Sun) const
{
    return _rotateMatrix(m_Mp, _calcRotAngleXY(R_Sun, m_GCC));
}

Matrix<double> Centroid::getcalcM_localRot(double R_Sun) const
{
    return _rotateMatrix(m_Mm+m_Mp, _calcRotAngleXY(R_Sun, m_GCC));
}

Centroid::KinematicParameters Centroid::getKPs() const
{
    return _toKPs(m_V0, m_Mm, m_Mp);
}

Centroid::KinematicParameters Centroid::getcalcKPs_localRot(double R_Sun) const
{
    Cartesian V0 = getcalcV0_localRot(R_Sun);
    Matrix<double> Mm = getcalcMm_localRot(R_Sun);
    Matrix<double> Mp = getcalcMp_localRot(R_Sun);
    return _toKPs(V0, Mm, Mp);
}

Centroid::KinematicParameters Centroid::getKPsErr() const
{
    return m_KPsErr;
}

void Centroid::setKPsErr(const KinematicParameters &KPsErr)
{
    m_KPsErr = KPsErr;
}

void Centroid::_expandTensor(Matrix<double> &Mm, Matrix<double> &Mp, const Matrix<double> &M)
{
    Mp = Matrix<double>(
        M[0][0], 0.5*(M[0][1]+M[1][0]), 0.5*(M[0][2]+M[2][0]),
        0.5*(M[0][1]+M[1][0]), M[1][1], 0.5*(M[1][2]+M[2][1]),
        0.5*(M[0][2]+M[2][0]), 0.5*(M[1][2]+M[2][1]), M[2][2]);
    Mm = M - Mp;
}

void Centroid::_expandKPs(Cartesian &V0, Matrix<double> &Mm, Matrix<double> &Mp, const KinematicParameters &KPs)
{
    V0 = Cartesian(KPs.u, KPs.v, KPs.w);
    Mm = Matrix<double>(
        0.0, -KPs.w3, KPs.w2,
        KPs.w3, 0.0, -KPs.w1,
       -KPs.w2, KPs.w1, 0.0);
    Mp = Matrix<double>(
        KPs.mp11, KPs.mp12, KPs.mp13,
        KPs.mp12, KPs.mp22, KPs.mp23,
        KPs.mp13, KPs.mp23, KPs.mp33);
}

Centroid::KinematicParameters Centroid::_toKPs(const Cartesian &V0, const Matrix<double> &Mm, const Matrix<double> &Mp)
{
    return KinematicParameters(
        V0.x, V0.y, V0.z, Mm[2][1], Mm[0][2], Mm[1][0],
        Mp[0][1], Mp[1][2], Mp[0][2], Mp[0][0], Mp[1][1], Mp[2][2]);
}

Matrix<double> Centroid::_rotateMatrix(const Matrix<double> &matrix, double angleXY)
{
    Matrix<double> R(
        std::cos(angleXY), -std::sin(angleXY), 0.0,
        std::sin(angleXY),  std::cos(angleXY), 0.0,
        0.0, 0.0, 1.0);
    return R * matrix * R.getcalcTranspose();
}

////////////////////////////////////
