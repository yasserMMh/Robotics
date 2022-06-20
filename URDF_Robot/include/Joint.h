#ifndef JOINT_H
#define JOINT_H

#include "MainHeaders.h"
#include "SetupCasADi.h"
using namespace SetupCasADi;
using namespace SetupCasADi::Denavit_Hartenberg;
using namespace SetupCasADi::SpecialMatrix;

//!< macro definition
/*!
    Three macros used to describe Joint type
    1. fixed
    2. revolut
    3. prismatic
    these macros are not used in the code at this stage,
    it was created for later expansion
*/
#define fixed           100
#define Revolut_Joint   101
#define Prismatic_Joint 102

/*! Joint class */
class Joint
{
public:
    //!< Custom non default constructor.
    /*!
      it takes:

                argument name       stands for       data type

        1           rpy           Roll Pitch Yaw       string

        2           xyz        positions at x, y, z    string

        3           axis         axis of rotation      string

        4         nth_joint                            string

        the constructor will also initialize Homogeneous Transformation Matrix
        using identity matrix
    */
    Joint(std::string rpy, std::string xyz, std::string axis, std::string nth_joint);

    //!< Custom default destructor.
    /*!
      default destructor has no functionality yet, it was created for later expansion
    */
    ~Joint();

public:
    int JointType; /**< int Joint Type, this field has no functionality yet, it was created for later expansion */
    int nth_joint; /**< a field used to distinguish a joint in a list */

public:
    Vector3Real rpy;  /**< rpy is a 3 x 1 vector holds 3 angles (Roll, Pitch, Yaw) */
    Vector3Real xyz;  /**< xyz is a 3 x 1 vector holds 3 points (x, y, z) to describe the position */
    Vector3Real axis; /**< axis is a 3 x 1 vector holds rotation axis */

    SymbReal q; /**< q is the generalized coordinate of the joint. q is special data type that has two fields (symbolic, value). */
    SymbReal dot_q; /**< dot_q is the derivative of generalized coordinate wrt time of the joint. dot_q is special data type that has two fields (symbolic, value). */

    Matrix4Symb HTM; /**< HTM is a 4 x 4 matrix describes Homogeneous Transformation Matrix of the joint */


public:
    //! A normal method taking no arguments and returning no value.
    /*!
      That method is used to trigger 2 methods
      1. set_HTM()
      2. moveJoint()
    */
    void Initialize();

    //! A normal method taking no arguments and returning no value.
    /*!
      set_HTM() method uses 2 members
        1. rpy
        2. xyz
      to set the Homogeneous Transformation Matrix
    */
    void set_HTM();

    /*!
      moveJoint() method uses one member
        1. axis
      to calculate the symbolic version of the Homogeneous Transformation Matrix
    */
    void moveJoint();

};
#endif // JOINT_H
