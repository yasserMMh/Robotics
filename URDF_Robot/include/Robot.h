#ifndef ROBOT_H
#define ROBOT_H

#include <iostream>
#include <fstream>
#include <rapidxml.hpp>

using namespace rapidxml;

#include "CenterPoint.h"
#include "SetupCasADi.h"

using namespace SetupCasADi;
using namespace SetupCasADi::Denavit_Hartenberg;
using namespace SetupCasADi::SpecialMatrix;

/*! Robot class template, where n is the total number of center points. */
/*! Robot is a group of center points, a center point is described by a joint and a link. */
template<int n>
class Robot
{
public:
    CenterPoint<n> head; /**< head is the base, and it is a field of data type CenterPoint. */
    int Current_Num_joints; /**< Current_Num_joints is an integer number used to track current defined joints. */

    Eigen::Matrix<Symb, 1, 1> PE; /**< PE is n x 1 vector describes the Potential Energy at each center point for the entire system. */

    Eigen::Matrix<Symb, n, n> MassMatrix; /**< MassMatrix is n x n matrix describes Mass of the entire system. */
    Eigen::Matrix<Symb, n, n> CoriolisMatrix; /**< CoriolisMatrix is n x n matrix describes Coriolis of the entire system. */
    Eigen::Matrix<Symb, n, 1> GravityMatrix; /**< GravityMatrix is n x 1 vector describes Gravity of the entire system. */

    Eigen::Matrix<Eigen::Matrix<Symb, n, n>, 1, n> divMassMatrix; /**< divMassMatrix is n x n x n matrix describes the derivative of Mass matrix wrt generalized coordinates. */
    Eigen::Matrix<Symb, n, n> divMassMatrixTime; /**< divMassMatrix is n x n matrix describes the derivative of Mass matrix wrt time. */

    std::vector<Symb> GeneralizedCoor; /**< GeneralizedCoor is a n x 1 vector holds generalized coordinates of the entire system. */

    Joint* eef; /**< eef stands for end effector. */


    xml_document<> doc; /**< doc is a field that holds the URDF file. */
	xml_node<>* root_node;

    xml_node<>* link_node;
	xml_node<>* inertial_node;
	xml_node<>* origin_node;
	xml_node<>* inertia_node;
	xml_node<>* mass_node;
	xml_node<>* visual_node;
	xml_node<>* geometry_node;
	xml_node<>* cylinder_node;

    //!< Custom default constructor.
    /*!
      it takes no inputs, it parse the URDF file.
    */
    Robot()
    {
        std::ifstream URDFfile("./URDF/diva_teleop.urdf"); /**< Path of URDF file. */

        if(URDFfile.is_open())
        {
            std::cout << "URDF file is successfully opened" << std::endl;
        }
        else
        {
            std::cout << "Failed to open URDF file" << std::endl;
        }

        std::string content((std::istreambuf_iterator<char>(URDFfile)), std::istreambuf_iterator<char>());
        char* URDFbuffer = doc.allocate_string(content.c_str());

        doc.parse<0>(URDFbuffer); /**< Parse the buffer using the xml file parsing library into doc. */
        root_node = doc.first_node("robot"); /**< Find the root node. */

        Current_Num_joints = 0; /**< Initialize number of joints to zero, since nothing has been yet appended. */
    }

    //!< Custom default destructor.
    /*!
      default destructor destroy all pointers of the Robot
    */
    ~Robot()
    {

        std::cout << "Destroy Ropot" << std::endl;
        CenterPoint<n>* cur = &head; /**< Start from the base. */
        CenterPoint<n>* prev = nullptr;
        while(cur != nullptr)
        {
            prev = cur;
            cur = cur->next;
            if(prev != &head)
            {
                delete prev;
                prev = nullptr;
            }
        }
    }

    //! A normal method taking no arguments and returning no value.
    /*!
      That method triggers 13 functions
         1. append all information from URDF
         2. calculate Homogeneous Transformation Matrix in symbolic form
         3. calculate Jacobian Matrix of Linear Velocity in symbolic form
         4. calculate derivative of Rotation matrix wrt generalized coordinates in symbolic form
         5. calculate skew matrix version of previous step in symbolic form
         6. calculate Jacobian Matrix of Rotational Velocity in symbolic form
         7. calculate Mass Matrix in symbolic form
         8. calculate Potential Energy in symbolic form
         9. calculate Gravity Matrix in symbolic form
        10. calculate derivative of Mass Matrix wrt generalized coordinates in symbolic form
        11. calculate Coriolis Matrix in symbolic form
        12. Combine all generalized coordinates and their derivative wrt time in one list in symbolic form
        13. calculate derivative of Mass Matrix wrt Time in symbolic form
    */
    void Initialize()
    {
        append();
        Initialize_HTM();
        Initialize_LVJac();
        Initialize_RotJac();
        Initialize_SkewRotJac();
        Initialize_RVJac();
        Calculate_MassMatrix();
        Calculate_PE();
        Calculate_GravityMatrix();
        Calculate_divMassMatrix();
        Calculate_CoriolisMatrix();
        Combine_GC();
        Calculate_divMassMatrixTime();
    }

    //! A normal method taking no arguments and returning no value.
    /*!
      That method extract information of the joints and links of the Robot from URDF file.
    */
    void append()
    {
        std::string NodeName;
        CenterPoint<n>* cur = &head;
        CenterPoint<n>* new_cp;

        for (xml_node<> * node = root_node->first_node("joint"); node; node = node->next_sibling("joint"))
        {

            if(Current_Num_joints >= n) /**< append end effector joint. */
            {
                eef = new Joint(
                                node->first_node("origin")->first_attribute("rpy")->value(),
                                node->first_node("origin")->first_attribute("xyz")->value(),
                                node->first_node("axis")->first_attribute("xyz")->value(),
                                std::to_string(n + 1)
                                );

                GeneralizedCoor.push_back(eef->q.symbolic);
                eef->Initialize();

                break;
            }

            int nth_cp = Current_Num_joints + 1;

            new_cp = new CenterPoint<n>();

            new_cp->init_Joint
            (
            node->first_node("origin")->first_attribute("rpy")->value(),
            node->first_node("origin")->first_attribute("xyz")->value(),
            node->first_node("axis")->first_attribute("xyz")->value(),
            nth_cp
            );

            GeneralizedCoor.push_back(new_cp->joint->q.symbolic);

            link_node       = node->next_sibling("link");
            inertial_node   = link_node->first_node("inertial");
            origin_node     = inertial_node->first_node("origin");
            inertia_node    = inertial_node->first_node("inertia");
            mass_node       = inertial_node->first_node("mass");
            visual_node     = link_node->first_node("visual");
            geometry_node   = visual_node->first_node("geometry");
            cylinder_node   = geometry_node->first_node("cylinder");

            new_cp->init_Link
            (
            origin_node->first_attribute("rpy")->value(),
            origin_node->first_attribute("xyz")->value(),
            std::stod(mass_node->first_attribute("value")->value()),
            std::stod(cylinder_node->first_attribute("radius")->value()),
            std::stod(cylinder_node->first_attribute("length")->value())
            );

            new_cp->Initialize
            (
            std::stod(inertia_node->first_attribute("ixy")->value()),
            std::stod(inertia_node->first_attribute("ixz")->value()),
            std::stod(inertia_node->first_attribute("iyz")->value()),
            std::stod(inertia_node->first_attribute("izz")->value()),
            std::stod(inertia_node->first_attribute("iyy")->value()),
            std::stod(inertia_node->first_attribute("ixx")->value())
            );

            cur->next = new_cp;
            cur = cur->next;

            Current_Num_joints++;
        }
    }

    //! A normal method taking one arguments and returning one value.
    /*!
      That method returns a pointer to the desired center point.
    */
    CenterPoint<n>* get_CenterPoint(int i)
    {
        if(Current_Num_joints < i || i <= 0) /**< check if i is a valid index joint. */
        {
            std::cout << "center point is not available in that position" << std::endl;
            return nullptr;
        }

        CenterPoint<n>* cur = &head; /**< Start from the base. */
        int k = 0;
        while(k < i) /**< check if there is more joints in the list or not. */
        {
            cur = cur->next; /**< move forward. */
            k++;
        }

        return cur;
    }

    //! A normal method taking no arguments and returning no value.
    /*!
      That method calculates Homogeneous Transformation Matrices at each center point for the entire system in symbolic form
    */
    void Initialize_HTM()
    {
        Matrix4Symb HTM_Temp = SpecialMatrix::Identity<Symb, 4>(); /**< Initialize the HTM matrix */

        CenterPoint<n>* cur = nullptr;
        for(int i = 1; i <= Current_Num_joints; i++)
        {
            cur = this->get_CenterPoint(i);
            HTM_Temp *= cur->joint->HTM;
            cur->HTM = HTM_Temp * cur->link->HTM; /**< multiply previous series to the HTM matrix of the current joint at center point */
        }
    }

    //! A normal method taking no arguments and returning no value.
    /*!
      That method calculates Jacobian of Linear Velocity Matrices for the entire system in symbolic form
    */
    void Initialize_LVJac()
    {
        for(int i = 1; i <= Current_Num_joints; i++)
        {
            this->get_CenterPoint(i)->Calculate_LVJac(GeneralizedCoor);
        }
    }

    //! A normal method taking no arguments and returning no value.
    /*!
      That method calculates calculates the derivative of Rotation matrix of size (3 x 3)
      wrt generalized coordinates of size (n x 1) for the entire system in symbolic form
    */
    void Initialize_RotJac()
    {
        for(int i = 1; i <= Current_Num_joints; i++)
        {
            this->get_CenterPoint(i)->Calculate_RotJac(GeneralizedCoor);
        }
    }

    //! A normal method taking no arguments and returning no value.
    /*!
      That method calculates the skew version of the derivative of rotation matrix for the entire system in symbolic form
    */
    void Initialize_SkewRotJac()
    {
        for(int i = 1; i <= Current_Num_joints; i++)
        {
            this->get_CenterPoint(i)->Calculate_SkewRotJac();
        }
    }

    //! A normal method taking no arguments and returning no value.
    /*!
      That method calculates the Jacobian of Rotation Velocity for the entire system in symbolic form
    */
    void Initialize_RVJac()
    {
        for(int i = 1; i <= Current_Num_joints; i++)
        {
            this->get_CenterPoint(i)->Calculate_RVJac();
        }
    }

    //! A normal method taking no arguments and returning no value.
    /*!
      That method calculates the Mass Matrix of the entire system in symbolic form
    */
    void Calculate_MassMatrix()
    {

        zeros(MassMatrix, n, n);
        CenterPoint<n>* cur = nullptr;

        for(int i = 1; i <= this->Current_Num_joints; i++)
        {
            cur = get_CenterPoint(i);

            MassMatrix = MassMatrix + 0.5 * cur->link->mass * cur->LVJac.transpose() * cur->LVJac + 0.5 * cur->RVJac.transpose() * cur->link->MOI * cur->RVJac;

        }

    }

    //! A normal method taking no arguments and returning no value.
    /*!
      That method calculates the Potential Energy of the entire system in symbolic form
    */
    void Calculate_PE()
    {
        this->PE(0) = 0;

        CenterPoint<n>* cur = nullptr;

        for(int i = 1; i <= this->Current_Num_joints; i++)
        {
            cur = get_CenterPoint(i);

            cur->Calculate_PE();

            PE += cur->PE;
        }
    }

    //! A normal method taking no arguments and returning no value.
    /*!
      That method calculates the Gravity Vector of the entire system in symbolic form
    */
    void Calculate_GravityMatrix()
    {
        for(int i = 0; i < this->Current_Num_joints; i++)
        {
            GravityMatrix(i) = casadi::SX::jacobian(this->PE(0), GeneralizedCoor[i]);

        }

    }

    //! A normal method taking no arguments and returning no value.
    /*!
      That method calculates the Gravity Vector of the entire system in symbolic form
    */
    void Calculate_divMassMatrix()
    {
        for(int k = 0; k < this->Current_Num_joints; k++)
        {
            for(int i = 0; i < n; i++)
            {
                for(int j = 0; j < n; j++)
                {

                    divMassMatrix(k)(i, j) = casadi::SX::jacobian(this->MassMatrix(i, j), GeneralizedCoor[k]);
                }
            }
        }
    }

    //! A normal method taking no arguments and returning no value.
    /*!
      That method calculates the Coriolis Matrix of the entire system in symbolic form
    */
    void Calculate_CoriolisMatrix()
    {
        zeros(CoriolisMatrix, n, n);
        CenterPoint<n>* cur = nullptr;

        for(int i = 0; i < n; i++)
        {
            for(int j = 0; j < n; j++)
            {
                for(int k = 0; k < n; k++)
                {
                    cur = get_CenterPoint(k+1);

                    CoriolisMatrix(i, j) += ( divMassMatrix(k)(i, j) + divMassMatrix(j)(i, k) - divMassMatrix(i)(k, j) ) *  (cur->joint->dot_q.symbolic);
                }
                CoriolisMatrix(i, j) *= 0.5;
            }
        }
    }

    //! A normal method taking 2 arguments and returning no value.
    /*!
      That method substitutes real numbers in symbolic form of Homogeneous Transformation Matrix
    */
    void solve_HTM(Eigen::Matrix<double, 4, 4>& sol, std::vector<Symb>& independentValue, bool include_eef = true)
    {
        if(include_eef)
        {
            sol = SolveMatrix<4, 4>(this->get_CenterPoint(n)->HTM * this->eef->HTM, GeneralizedCoor, independentValue);
        }
        else
        {
            sol = SolveMatrix<4, 4>(this->get_CenterPoint(n)->HTM, GeneralizedCoor, independentValue);
        }
    }

    //! A normal method taking 2 arguments and returning no value.
    /*!
      That method substitutes real numbers in symbolic form of Mass Matrix
    */
    void solve_MassMatrix(Eigen::Matrix<double, n, n>& sol, std::vector<Symb>& independentValue)
    {
        sol = SolveMatrix<n, n>(this->MassMatrix, GeneralizedCoor, independentValue);
    }

    void Combine_GC()
    {
        CenterPoint<n>* cur = nullptr;
        for(int i = 0; i < n; i++)
        {
            cur = get_CenterPoint(i+1);

            GeneralizedCoor.push_back(cur->joint->dot_q.symbolic);
        }
    }

    //! A normal method taking 2 arguments and returning no value.
    /*!
      That method substitutes real numbers in symbolic form of Coriolis Matrix
    */
    void solve_CoriolisMatrix(Eigen::Matrix<double, n, n>& sol, std::vector<Symb>& independentValue)
    {
        sol = SolveMatrix<n, n>(this->CoriolisMatrix, GeneralizedCoor, independentValue);
    }

    //! A normal method taking 2 arguments and returning no value.
    /*!
      That method substitutes real numbers in symbolic form of Gravity Matrix
    */
    void solve_GravityMatrix(Eigen::Matrix<double, n, 1>& sol, std::vector<Symb>& independentValue)
    {
        sol = SolveMatrix<n, 1>(this->GravityMatrix, GeneralizedCoor, independentValue);
    }

    //! A normal method taking 2 arguments and returning no value.
    /*!
      That method substitutes real numbers in symbolic form of Jacobian Matrix of Linear Velocity
    */
    void solve_LVJac(Eigen::Matrix<double, 3, n>& sol, std::vector<Symb>& independentValue)
    {
        CenterPoint<n>* cur = nullptr;
        sol = Eigen::MatrixXd::Zero(3, n);
        for(int i = 0; i < Current_Num_joints; i++)
        {
            cur = get_CenterPoint(i+1);
            sol += SolveMatrix<3, n>(cur->LVJac, GeneralizedCoor, independentValue);
        }
    }

    //! A normal method taking 2 arguments and returning no value.
    /*!
      That method substitutes real numbers in symbolic form of Jacobian Matrix of Rotational Velocity
    */
    void solve_RVJac(Eigen::Matrix<double, 3, n>& sol, std::vector<Symb>& independentValue)
    {
        CenterPoint<n>* cur = nullptr;
        sol = Eigen::MatrixXd::Zero(3, n);
        for(int i = 0; i < Current_Num_joints; i++)
        {
            cur = get_CenterPoint(i+1);
            sol += SolveMatrix<3, n>(cur->RVJac, GeneralizedCoor, independentValue);
        }
    }

    //! A normal method taking no arguments and returning no value.
    /*!
      That method calculates derivative of Mass Matrix wrt Time in symbolic form
    */
    void Calculate_divMassMatrixTime()
    {
        zeros(divMassMatrixTime, n, n);
        for(int i = 0; i < Current_Num_joints; i++)
        {
            divMassMatrixTime += divMassMatrix(i) * GeneralizedCoor[i+n+1];
        }
    }

    //! A normal method taking 2 arguments and returning no value.
    /*!
      That method substitutes real numbers in symbolic form of derivative of Mass Matrix wrt Time
    */
    void solve_divMassMatrixTime(Eigen::Matrix<double, n, n>& sol, std::vector<Symb>& independentValue)
    {
        sol = SolveMatrix<n, n>(this->divMassMatrixTime, GeneralizedCoor, independentValue);
    }

};
#endif // ROBOT_H
