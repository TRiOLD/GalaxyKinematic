////////////////////////////////////
#ifndef STAR_H
#define STAR_H

////////////////////////////////////
class Star
{
    ////////////////////////////////
public:
    struct Equatorial
    {
        double ra;
        double de;
        double r;

        Equatorial()
            : ra(0.0), de(0.0), r(0.0) {}
        Equatorial(double ra, double de, double r)
            : ra(ra), de(de), r(r) {}
    };

    struct Galactic
    {
        double l;
        double b;
        double r;

        Galactic()
            : l(0.0), b(0.0), r(0.0) {}
        Galactic(double l, double b, double r)
            : l(l), b(b), r(r) {}
    };

    struct Cartesian
    {
        double x;
        double y;
        double z;

        Cartesian()
            : x(0.0), y(0.0), z(0.0) {}
        Cartesian(double x, double y, double z)
            : x(x), y(y), z(z) {}

    };

    struct Cylindrical
    {
        double R;
        double theta;
        double Z;

        Cylindrical()
            : R(0.0), theta(0.0), Z(0.0) {}
        Cylindrical(double R, double theta, double Z)
            : R(R), theta(theta), Z(Z) {}
    };

    ////////////////////////////////

public:
    virtual ~Star();

    Star();
    Star(const Equatorial &EqC, const Equatorial &EqV);
    Star(const Galactic &GC, const Galactic &GV);
    Star(const Cartesian &GCC, const Cartesian &GCV);

protected:
    Cartesian m_GCC;
    Cartesian m_GCV;
    void* m_addData_ptr;

public:
    Cartesian getGCC() const;
    Cartesian getGCV() const;

    void *getAddData_ptr() const;
    void setAddData_ptr(void *addData_ptr);

 // void getcalcEqCandV(Equatorial& EqC, Equatorial& EqV) const;
    void getcalcGCandV(Galactic &GC, Galactic& GV) const;
    void getcalcGCCCandV(Cylindrical &GCCC, Cylindrical &GCCV,
                         double R_Sun, const Cartesian &GCV_Sun) const;

    Cartesian getcalcGCC_local(const Cartesian &GCC0) const;
    void getcalcCCandV_localRot(Cartesian &GCC_localRot, Cartesian &GCV_localRot,
                                const Cartesian &GCC0, double R_Sun) const;

private:
    static Galactic _calcGC(const Equatorial &EqC);
    static Galactic _calcGV(const Equatorial &EqC, const Equatorial &EqV);
    static Galactic _calcGC(const Cartesian &GCC);
    static Galactic _calcGV(const Galactic &GC, const Cartesian &GCV);
    static Cartesian _calcGCC(const Galactic &GC);
    static Cartesian _calcGCV(const Galactic &GC, const Galactic &GV);

protected:
    static double _calcRotAngleXY(double R_Sun, const Cartesian &GCC);
    static Cartesian _rotateVector(const Cartesian &vector, double angleXY);
};

////////////////////////////////////
#endif // STAR_H
