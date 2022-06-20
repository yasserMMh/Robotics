#ifndef CENTERPOINT_H_INCLUDED
#define CENTERPOINT_H_INCLUDED

#include "Joint.h"
#include "Link.h"

using namespace SetupCasADi;
using namespace SetupCasADi::Denavit_Hartenberg;
using namespace SetupCasADi::SpecialMatrix;

/*! CenterPoint class template, where n is the total number of center points. */
/*! Center Point is a group of a joint and a link. */
template <int n>
class CenterPoint
{
public:

    //!< Custom default constructor.
    /*!
      default constructor has no functionality yet, it was created for later expansion
    */
    CenterPoint()
    {

    }

public:
    Joint* joint; /**< pointer to a field of data type Joint. */
    Link* link;   /**< pointer to a field of data type Link. */
    CenterPoint<n>* next; /**< define the next center point in the linked list. */

    int nth_cp; /**< a field used to distinguish a center point in a list. */

    Matrix4Symb HTM; /**< HTM is a 4 x 4 matrix describes Homogeneous Transformation Matrix of the center point. this matrix is a result of HTM of joint x HTM of link. */

    Eigen::Matrix<Symb, 3, n> LVJac; /**< LVJac is a 3 x n matrix describes Jacobian Matrix of Linear Velocity of the center point. */
    Eigen::Matrix<Symb, 3, n> RVJac; /**< RVJac is a 3 x n matrix describes Jacobian Matrix of Rotational Velocity of the center point. */

    Eigen::Matrix<Matrix3Symb, n, 1> RotJac; /**< RotJac is a 3 x 3 x n matrix describes derivative of Rotation matrix wrt generalized coordinates of the center point. */
    Eigen::Matrix<Matrix3Symb, n, 1> SkewRotJac; /**< SkewRotJac is a 3 x 3 x n skew matrix, it is a result of transpose(Rotation matrix) * RotJac matrix. */

    Eigen::Matrix<Symb, 1, 1> PE; /**< scalar Potential Energy at the center point. */

public:

    //! A normal method taking 4 arguments and returning no value.
    /*!
      Inputs:

                argument name       stands for       data type

        1           rpy           Roll Pitch Yaw       string

        2           xyz        positions at x, y, z    string

        3           axis         axis of rotation      string

        4         nth_joint                              int

        the method create data in the heap and attach its address
        to the pointer joint
    */
    void init_Joint(std::string rpy, std::string xyz, std::string axis, int nth_cp)
    {
        this->nth_cp = nth_cp;
        this->joint = new Joint(rpy, xyz, axis, std::to_string(nth_cp));
    }

    //! A normal method taking 5 arguments and returning no value.
    /*!
      Inputs:

                argument name       stands for       data type

        1           rpy           Roll Pitch Yaw       string

        2           xyz        positions at x, y, z    string

        3           mass                               double

        4          radius                              double

        5          length                              double

        the method create data in the heap and attach its address
        to the pointer link
    */
    void init_Link(std::string rpy, std::string xyz, double mass, double radius, double length)
    {
        this->link = new Link(rpy, xyz, mass, radius, length);
    }

    //! A normal method taking 6 double arguments and returning no value.
    /*!
      That method is used to trigger 2 methods
      1. joint->Initialize()
      2. link->Initialize()
    */
    void Initialize(double ixy, double ixz, double iyz, double izz, double iyy, double ixx)
    {
        joint->Initialize();
        link->Initialize(ixy, ixz, iyz, izz, iyy, ixx);

    }

    //! A normal method taking one argument and returning no value.
    /*!
      Inputs:
            1. GeneralizedCoor: is a list of all generalized coordinates of the entire system
      That method calculates the Jacobian of Linear Velocity
    */
    void Calculate_LVJac(std::vector<Symb>& GeneralizedCoor)
    {
        zeros(LVJac, 3, n);
        for(int i = 0; i < this->nth_cp; i++)
        {
            for(int j = 0; j < 3; j++)
            {
                this->LVJac(j, i) = Symb::jacobian(this->HTM(j, 3), GeneralizedCoor[i]);
            }
        }
    }

    //! A normal method taking one argument and returning no value.
    /*!
      Inputs:
            1. GeneralizedCoor: is a list of all generalized coordinates of the entire system
      That method calculates the derivative of Rotation matrix of size (3 x 3)
      wrt generalized coordinates of size (n x 1)
    */
    void Calculate_RotJac(std::vector<Symb>& GeneralizedCoor)
    {
        int Rotsize = 3; /**< size of rotation matrix is 3 x 3 */
        for(int k = 0; k < this->nth_cp; k++)
        {
            for(int i = 0; i < Rotsize; i++)
            {
                for(int j = 0; j < Rotsize; j++)
                {
                    this->RotJac(k)(i, j) = Symb::jacobian(this->HTM(i, j), GeneralizedCoor[k]);
                }
            }
        }
    }

    //! A normal method taking no argument and returning no value.
    /*!
      That method calculates the skew version of the derivative of rotation matrix
    */
    void Calculate_SkewRotJac()
    {
        Matrix3Symb Rot = this->HTM.topLeftCorner(3, 3);/**< extract rotation matrix from Homogeneous Transformation Matrix. */

        for(int i = 0; i < this->nth_cp; i++)
        {
            SkewRotJac(i) = Rot.transpose() * this->RotJac(i);
        }
    }

    //! A normal method taking no argument and returning no value.
    /*!
      That method calculates the Jacobian of Rotation Velocity
    */
    void Calculate_RVJac()
    {
        zeros(RVJac, 3, n);
        for(int i = 0; i < nth_cp; i++)
        {
            this->RVJac.middleCols(i, 1) = SkewToVector(this->SkewRotJac(i)); /**< SkewToVector() is a simple function that convert skew matrix to vector of 3 x 1. */
        }
    }

    //! A normal method taking no argument and returning no value.
    /*!
      That method calculates Potential energy at each center point
    */
    void Calculate_PE()
    {
        this->PE = this->link->mass * g.transpose() * this->HTM.topRightCorner(3, 1);
    }

};

#endif // CENTERPOINT_H_INCLUDED
