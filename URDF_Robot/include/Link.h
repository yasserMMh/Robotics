#ifndef LINK_H_INCLUDED
#define LINK_H_INCLUDED

#include "MainHeaders.h"
#include "SetupCasADi.h"

using namespace SetupCasADi;
using namespace SetupCasADi::Denavit_Hartenberg;
using namespace SetupCasADi::SpecialMatrix;

/*! Link class */
class Link
{
public:
    //!< Custom default constructor.
    /*!
      default constructor has no functionality yet, it was created for later expansion
    */
    Link();

    //!< Custom non default constructor.
    /*!
      it takes:
                argument name       stands for       data type

        1           rpy           Roll Pitch Yaw       string

        2           xyz        positions at x, y, z    string

        3           mass                               double

        4          radius                              double

        5          length                              double

        the method will also initialize Homogeneous Transformation Matrix
        using identity matrix
    */
    Link(std::string rpy, std::string xyz, double mass, double radius, double length);

    //!< Custom default destructor.
    /*!
      default destructor has no functionality yet, it was created for later expansion
    */
    ~Link();

public:
    Vector3Real rpy; /**< rpy is a 3 x 1 vector holds 3 angles (Roll, Pitch, Yaw) */
    Vector3Real xyz; /**< xyz is a 3 x 1 vector holds 3 points (x, y, z) to describe the position */

    Matrix3Symb MOI; /**< MOI is a 3 x 3 matrix describes the moment of inertia tensor */

    double mass;   /**< mass of the link */
    double radius; /**< radius of the link */
    double length; /**< length of the link */

    Matrix4Symb HTM; /**< HTM is a 4 x 4 matrix describes Homogeneous Transformation Matrix of the link */

public:
    //! A normal method taking 6 double arguments and returning no value.
    /*!
      That method is used to trigger 2 methods
      1. set_MOI()
      2. set_HTM()
    */
    void Initialize(double ixy, double ixz, double iyz, double izz, double iyy, double ixx);

    //! A normal method taking 6 double arguments and returning no value.
    /*!
      The 6 double arguments are inertia of
        1. xy
        2. xz
        3. yz
        4. zz
        5. yy
        6. xx
      used to set the Moment of Inertia Matrix
    */
    void set_MOI(double ixy, double ixz, double iyz, double izz, double iyy, double ixx);

    //! A normal method taking no arguments and returning no value.
    /*!
      set_HTM() method uses 2 members
        1. rpy
        2. xyz
      to set the Homogeneous Transformation Matrix
    */
    void set_HTM();

};
#endif // LINK_H_INCLUDED
