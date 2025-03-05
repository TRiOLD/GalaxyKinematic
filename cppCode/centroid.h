////////////////////////////////////
#ifndef CENTROID_H
#define CENTROID_H

////////////////////////////////////
#include "star.h"

#include <TMatrix.h>
using namespace TRiOLD;

////////////////////////////////////
class Centroid : public Star
{
    ////////////////////////////////
public:
    struct KinematicParameters
    {
        double u, v, w;
        double w1, w2,w3;
        double mp12, mp23, mp13;
        double mp11, mp22, mp33;

        KinematicParameters()
            :u(0.0), v(0.0), w(0.0), w1(0.0), w2(0.0), w3(0.0),
             mp12(0.0), mp23(0.0), mp13(0.0), mp11(0.0), mp22(0.0), mp33(0.0)
        {}

        KinematicParameters(double u, double v, double w, double w1, double w2, double w3,
        double mp12, double mp23, double mp13, double mp11, double mp22, double mp33)
            :u(u), v(v), w(w), w1(w1), w2(w2),w3(w3),
             mp12(mp12), mp23(mp23), mp13(mp13), mp11(mp11), mp22(mp22), mp33(mp33)
        {}
    };
    ////////////////////////////////

public:
    Centroid();
    Centroid(const Cartesian& GCC, const Cartesian& GCV, std::size_t starsAmount,
             const Cartesian& V0, const Matrix<double>& M);
    Centroid(const Cartesian& GCC, const Cartesian& GCV, std::size_t starsAmount,
             const Cartesian& V0, const Matrix<double>& Mm, const Matrix<double>& Mp);
    Centroid(const Cartesian& GCC, const Cartesian& GCV, std::size_t starsAmount,
             const KinematicParameters& KPs);

protected:
    std::size_t m_starsAmount;
    Cartesian m_V0;
    Matrix<double> m_Mm;
    Matrix<double> m_Mp;

    KinematicParameters m_KPsErr;

public:
    std::size_t getStarsAmount() const;

    Cartesian getV0() const;
    Matrix<double> getMm() const;
    Matrix<double> getMp() const;
    Matrix<double> getcalcM() const;
    Cartesian getcalcV0_localRot(double R_Sun) const;
    Matrix<double> getcalcMm_localRot(double R_Sun) const;
    Matrix<double> getcalcMp_localRot(double R_Sun) const;
    Matrix<double> getcalcM_localRot(double R_Sun) const;

    KinematicParameters getKPs() const;
    KinematicParameters getcalcKPs_localRot(double R_Sun) const;

    KinematicParameters getKPsErr() const;
    void setKPsErr(const KinematicParameters& KPsErr);

private:
    static void _expandTensor(Matrix<double>& Mm, Matrix<double>& Mp, const Matrix<double>& M);
    static void _expandKPs(Cartesian& V0, Matrix<double>& Mm, Matrix<double>& Mp, const KinematicParameters& KPs);
    static KinematicParameters _toKPs(const Cartesian& V0, const Matrix<double>& Mm, const Matrix<double>& Mp);

    static Matrix<double> _rotateMatrix(const Matrix<double>& matrix, double angleXY);

};

////////////////////////////////////
#endif // CENTROID_H
