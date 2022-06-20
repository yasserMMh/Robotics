#include "Link.h"

Link::Link(std::string rpy, std::string xyz, double mass, double radius, double length)
{
    StringToVector(this->rpy, rpy);
    StringToVector(this->xyz, xyz);

    this->HTM = SpecialMatrix::Identity<Symb, 4>(); /**< Initialize the HTM matrix */

    this->mass   = mass;
    this->radius = radius;
    this->length = length;

}

Link::~Link()
{

}

void Link::Initialize(double ixy, double ixz, double iyz, double izz, double iyy, double ixx)
{
    set_MOI(ixy, ixz, iyz, izz, iyy, ixx);
    set_HTM();

}

void Link::set_MOI(double ixy, double ixz, double iyz, double izz, double iyy, double ixx)
{
    this->MOI << ixx, ixy, ixz,
                 ixy, iyy, iyz,
                 ixz, iyz, izz;
}

void Link::set_HTM()
{
    this->HTM(0, 3) = this->xyz(0);
    this->HTM(1, 3) = this->xyz(1);
    this->HTM(2, 3) = this->xyz(2);

    //! A normal function takes 4 arguments and returning no value.
    /*!
      RotRPY stands for Rotation Roll Pitch Yaw
      the prototype of the function is found in SetupCasADi.h file
      the implementation of the function is found in SetupCasADi.cpp file
      this function takes 4 arguments
      1. a reference of Homogeneous Transformation Matrix
      2. roll angle
      3. pitch angle
      4. yaw angle
    */
    RotRPY(this->HTM , rpy(0), rpy(1), rpy(2));
}

