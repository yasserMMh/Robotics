#include "Joint.h"

Joint::Joint(std::string rpy, std::string xyz, std::string axis, std::string nth_joint)
{
    StringToVector(this->rpy, rpy);
    StringToVector(this->xyz, xyz);
    StringToVector(this->axis, axis);

    this->HTM = SpecialMatrix::Identity<Symb, 4>(); /**< Initialize the HTM matrix */

    this->q.symbolic = Symb::sym( "q" + nth_joint);
    this->dot_q.symbolic = Symb::sym( "dot_q" + nth_joint);
}

Joint::~Joint()
{

}

void Joint::Initialize()
{
    this->set_HTM();
    this->moveJoint();
}

void Joint::set_HTM()
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
    RotRPY(this->HTM, rpy(0), rpy(1), rpy(2));

}

void Joint::moveJoint()
{

    if(this->axis(0) == 1)
    {
        Roll( this->HTM, this->q.symbolic);
    }
    else if(this->axis(1) == 1)
    {
        Pitch(this->HTM, this->q.symbolic);
    }
    else if(this->axis(2) == 1)
    {
        Yaw(  this->HTM, this->q.symbolic);
    }

}
