#ifndef SETUPCASADI_H_INCLUDED
#define SETUPCASADI_H_INCLUDED

#include <iostream>
#include <string>
#include <eigen3/Eigen/Dense>
# include <cppad/cppad.hpp>
#include <cmath>
#include <math.h>

#include <casadi/casadi.hpp>
#include <casadi/core/sx_fwd.hpp>
#include <casadi/core/sparsity.hpp>
#include <casadi/core/calculus.hpp>

namespace SetupCasADi
{
    /****************************************************************************** *
     * for simplicity, set an alias to                                              *
     *      1. casadi::SX ---------------------> Symb                               *
     *      2. Eigen::Matrix<Symb, n, c> ------> Matrix(n)Symb                      *
     * ******************************************************************************/
    typedef casadi::SX Symb; // called be used to declare scalar or vector
    /* Set alias, Matrix4Symb stands for matrix (4 x 4) size symbolic */
    typedef Eigen::Matrix<Symb, 4, 4> Matrix4Symb;
    typedef Eigen::Matrix<Symb, 3, 3> Matrix3Symb;
    typedef Eigen::Matrix<Symb, 3, 1> Vector3Symb;
    typedef Eigen::MatrixX<Symb>      MatrixSymb ;

    /* following data types are alias that hold actual real numbers in matrix of size (nrow x ncol) */
    typedef Eigen::Matrix<double, 4, 4> Matrix4Real;
    typedef Eigen::Matrix<double, 3, 3> Matrix3Real;
    typedef Eigen::Matrix<double, 3, 1> Vector3Real;

    /* Define 3 orthogonal unit axis x, y, z */
    namespace CartesianAxis
    {
        extern Vector3Symb x;
        extern Vector3Symb y;
        extern Vector3Symb z;
    }

    extern Vector3Symb g;

    /********************************************************* *
     * SpecialMatrix namespace has 2 functions used to create  *
     *      1. (Msize x Msize) Identity matrix                 *
     *      2. (nrow  x ncol ) Zero matrix                     *
     * The functions are applicable for any data type          *
     * *********************************************************/
    namespace SpecialMatrix
    {
        template<typename T, int Msize>
        Eigen::Matrix<T, Msize, Msize> Identity()
        {
            /* define the matrix */
            Eigen::Matrix<T, Msize, Msize> IdentityMat;

            for (int i = 0; i < Msize; i++)
            {
                for (int j = 0; j < Msize; j++)
                {
                    if(i == j)
                    {
                        IdentityMat(i , j) = 1;
                    }
                    else
                    {
                        IdentityMat(i , j) = 0;
                    }

                }
            }
            return IdentityMat;
        }

        template<typename T, int nrow, int ncol>
        Eigen::Matrix<T, nrow, ncol> Zero()
        {
            /* define the matrix */
            Eigen::Matrix<T, nrow, ncol> ZeroMat;
            for (int i = 0; i < nrow; i++)
            {
                for (int j = 0; j < ncol; j++)
                {
                    ZeroMat(i , j) = 0;

                }
            }
            return ZeroMat;
        }

        template<typename T>
        void PrintSymb(T& MatSymb, int rows, int cols)
        {
            std::cout << "{" << std::endl;
            for(int i = 0; i < rows; i++)
            {
                std::cout << "  [";
                for(int j = 0; j < cols; j++)
                {
                    std::setw(5);
                    std::cout << ' ' << std::setw(3) << MatSymb(i, j);
                }
                std::cout << "]\n";
            }
            std::cout << "}" << std::endl;
        }

        template<typename T>
        void PrintReal(T& MatSymb, int rows, int cols)
        {
            std::cout << "{" << std::endl;
            for(int i = 0; i < rows; i++)
            {
                std::cout << "  [";
                for(int j = 0; j < cols; j++)
                {
                    std::cout << std::fixed;
                    std::cout << std::setprecision(6);
                    std::setw(5);
                    std::cout << ' ' << std::setw(10) << MatSymb(i, j);
                }
                std::cout << "]\n";
            }
            std::cout << "}" << std::endl;
        }


        template<typename T>
        void zeros(T& mat, int rows, int cols)
        {
            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < cols; j++)
                {
                    mat(i , j) = 0;

                }
            }
        }

    }

    namespace Denavit_Hartenberg
    {
        /* ********************** *
         * set names for          *
         *      1. offset         *
         *      2. theta          *
         *      3. link length    *
         *      4. alpha          *
         * ********************** */
        extern std::string offsetName;
        extern std::string thetaName;
        extern std::string linkLenName;
        extern std::string alphaName;

        extern std::string massName;
        extern std::string radiusName;
        /* ************************************************************* *
         * Define a structure that holds two variables of type           *
         *      1. Symb: casadi::SX                                      *
         *      2. double                                                *
         * this structure would be used to define fields of DH_Parameter *
         * ************************************************************* */
        struct SymbReal
        {
            Symb symbolic;
            double value;
            SymbReal() : value(0) {}; // initialize value to zero
        };

        /* ************************************************************* *
         * Properties of the link material assuming solid cylinder shape *
         * ************************************************************* */
        struct LinkProp
        {
            SymbReal mass;
            SymbReal radius;
            /* Override the Automatic Constructor */
            LinkProp()
            { mass.symbolic = Symb::sym(massName), radius.symbolic = Symb::sym(radiusName); };
            /* user defined Constructor */
            LinkProp(double mass, double radius, std::string nth_joint);
            /* ***************************************************************************** *
             * Before accepting the inputs from the user check the validity of the inputs    *
             * ***************************************************************************** */
            double set_inputs(double input);
        };

        /* ************************************************************** *
         * Merge the 4 Denavit Hartenberg parameters into one structure   *
         * d    : offset                                                  *
         * theta: angle along z axis                                      *
         * r    : link length                                             *
         * alpha: angle along x axis                                      *
         * (d, theta, r, alpha) are names stored in previously defined    *
         * variables -> the names can be changed in SetupCasADi.cpp       *
         * every field (offset, theta, linkLen, alpha) has symbolic and   *
         * numeric data                                                   *
         * ************************************************************** */
        struct DH_Parameter
        {
            /* The parameters are public and accessible by default */
            SymbReal offset;
            SymbReal theta;
            SymbReal linkLen;
            SymbReal alpha;

            /* Override the Automatic Constructor */
            DH_Parameter()
            { offset.symbolic = Symb::sym(offsetName), theta.symbolic = Symb::sym(thetaName), linkLen.symbolic = Symb::sym(linkLenName), alpha.symbolic = Symb::sym(alphaName); };
            /* user defined Constructor */
            DH_Parameter(double offset, double theta, double linkLen, double alpha, std::string nth_joint, bool is_Rad = false);
            /* ***************************************************************************** *
             * Before accepting the inputs from the user check the validity of the inputs    *
             * that is why there are 2 functions                                             *
             *      1. set_Angle()                                                           *
             *      2. set_Length()                                                          *
             * ***************************************************************************** */
            double set_Angle(double angle, bool is_Rad);
            double set_Length(double length);
        };

        Matrix3Symb Skew(const Vector3Symb& r);
        Vector3Symb SkewToVector(const Matrix3Symb& mat);
        bool CheckSkew(const Matrix3Symb& mat);
        Matrix3Symb OuterProduct(const Vector3Symb& r1, const Vector3Symb& r2);
        Matrix4Symb RotSymb(const Vector3Symb& r, Symb& angleSymb);
        Matrix4Symb TransSymb(const Vector3Symb& r, Symb& distanceSymb, bool is_CPOL = false); // CPOL stands for Center Point of the Link
        Matrix4Symb HomogeneousMatrix(DH_Parameter& para, bool is_CPOL = false);
        /* MOIT: Moment of Inertia Tensor */
        Matrix3Symb MOIT(LinkProp& prop, Symb& linkLength);

        void StringToVector(Vector3Real& vector3D, std::string& str);
        void RotRPY(Matrix4Symb& rotMat, double roll, double pitch, double yaw);
        void Roll( Matrix4Symb& rotMat, Symb& angle);
        void Pitch(Matrix4Symb& rotMat, Symb& angle);
        void Yaw(  Matrix4Symb& rotMat, Symb& angle);

        template<int nrows, int ncols>
        Eigen::Matrix<double, nrows, ncols> SolveMatrix(const Eigen::Matrix<Symb, nrows, ncols>& MatSymb, std::vector<Symb>& independentSymb, std::vector<Symb>& independentValue)
        {
            //int MatSize = 4;
            /* Create a (4 x 4) matrix to hold members of real elements matrix */
            Eigen::Matrix<double, nrows, ncols> MatReal = Eigen::MatrixXd::Zero(nrows, ncols);
            std::vector<Symb> OutputVector(nrows * ncols);

            for(int i = 0; i < nrows * ncols; i++)
            {
                OutputVector[i] = MatSymb(i);
            }

            auto f = casadi::Function("f", independentSymb, OutputVector);
            /* Solve functions given the inputs <independentValue> */
            auto sol = f(independentValue);
            for(int i = 0; i < nrows * ncols; i++)
            {
                MatReal(i) = static_cast<double>(sol[i]);
            }

            return MatReal;
        }

    }




}




#endif // SETUPCASADI_H_INCLUDED
