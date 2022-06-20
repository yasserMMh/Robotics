#include <iostream>
#include <fstream>
#include <typeinfo>
#include "SetupCasADi.h"

#include "Link.h"
#include "Joint.h"
#include "CenterPoint.h"
#include "Robot.h"

#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>

#include <casadi/casadi.hpp>
#include <casadi/core/sx.hpp>
#include <casadi/core/sx_fwd.hpp>
#include <casadi/core/nlpsol.hpp>
#include <casadi/core/sparsity.hpp>
#include <casadi/core/calculus.hpp>

#include <vector>

using namespace SetupCasADi;
using namespace SetupCasADi::Denavit_Hartenberg;
using namespace SetupCasADi::SpecialMatrix;

int main()
{

    const int n = 7;

    Robot<n>* robot;
    robot->Initialize();

    std::vector<Symb> independentValue = {-1, -0.23, 0.43, 0.32, 0, 0, 0, 0, 1, 1, 1, 2.5, 1, 1, 1.3};

    std::cout << "Homogeneous Transformation Matrix: " << std::endl;
    Eigen::Matrix<double, 4, 4> HTM;
    robot->solve_HTM(HTM, independentValue);
    PrintReal(HTM, 4, 4);
    std::cout << std::endl;

    std::cout << "Mass Matrix: " << std::endl;
    Eigen::Matrix<double, n, n> Mass;
    robot->solve_MassMatrix(Mass, independentValue);
    PrintReal(Mass, n, n);
    std::cout << std::endl;

    std::cout << "Test Mass Matrix: " << std::endl;
    PrintReal(Mass.transpose() - Mass, n, n);
    std::cout << std::endl;

    std::cout << "Coriolis Matrix: " << std::endl;
    Eigen::Matrix<double, n, n> Coriolis;
    robot->solve_CoriolisMatrix(Coriolis, independentValue);
    PrintReal(Coriolis, n, n);
    std::cout << std::endl;

    std::cout << "Test Coriolis Matrix: " << std::endl;
    Eigen::Matrix<double, n, n> divMassMatrixTime;
    robot->solve_divMassMatrixTime(divMassMatrixTime, independentValue);
    Eigen::Matrix<double, n, n> temp = divMassMatrixTime - 2*Coriolis;
    Eigen::Matrix<double, n, n> temp2 = temp + temp.transpose();
    PrintReal(temp2, n, n);
    std::cout << std::endl;

    std::cout << "Gravity Matrix: " << std::endl;
    Eigen::Matrix<double, n, 1> Gravity;
    robot->solve_GravityMatrix(Gravity, independentValue);
    PrintReal(Gravity, n, 1);
    std::cout << std::endl;

    std::cout << "Jacobian of Linear Velocity: " << std::endl;
    Eigen::Matrix<double, 3, n> LVJac;
    robot->solve_LVJac(LVJac, independentValue);
    PrintReal(LVJac, 3, n);
    std::cout << std::endl;

    std::cout << "Jacobian of Rotational Velocity: " << std::endl;
    Eigen::Matrix<double, 3, n> RVJac;
    robot->solve_RVJac(RVJac, independentValue);
    PrintReal(RVJac, 3, n);
    std::cout << std::endl;

}

