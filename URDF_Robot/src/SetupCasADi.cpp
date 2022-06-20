#include "SetupCasADi.h"

namespace SetupCasADi
{
    /*!
      Define 3 orthogonal unit axis x, y, z
    */
    namespace CartesianAxis
    {
        Vector3Symb x(1, 0, 0);
        Vector3Symb y(0, 1, 0);
        Vector3Symb z(0, 0, 1);
    }
    Vector3Symb g(0, 0, 9.8);


    namespace Denavit_Hartenberg
    {
        std::string offsetName  = "d";
        std::string thetaName   = "t";
        std::string linkLenName = "l";
        std::string alphaName   = "a";

        std::string massName    = "m";
        std::string radiusName  = "r";

        /*!
          Define constructor of struct LinkProp
        */
        LinkProp::LinkProp(double mass, double radius, std::string nth_joint)
        {
            this->mass.symbolic   = Symb::sym(  massName + nth_joint);
            this->radius.symbolic = Symb::sym(radiusName + nth_joint);

            this->mass.value      = set_inputs(mass);
            this->radius.value    = set_inputs(radius);

        }

        /*!
          Check the validity of the mass and radius, both must be positive
        */
        double LinkProp::set_inputs(double input)
        {
            /*!
              Unfortunately, double numbers may vary with small values such as 1e-10
              Define Epsilon variable to hold small number as toleration
            */
            double Epsilon = -1e-10;
            if (input < Epsilon)
            {
                std::cout << "mass or radius of the link is/are a/ negative number/s" << std::endl;
                std::cout << "Correct your input/s" << std::endl;
                exit (EXIT_FAILURE);
            }
            return input;
        }

        /*!
          define constructor of struct DH_Parameter
        */
        DH_Parameter::DH_Parameter(double offset, double theta, double linkLen, double alpha,std::string nth_joint, bool is_Rad)
        {
            this->offset.symbolic  = Symb::sym( offsetName + nth_joint);
            this->theta.symbolic   = Symb::sym(  thetaName + nth_joint);
            this->linkLen.symbolic = Symb::sym(linkLenName + nth_joint);
            this->alpha.symbolic   = Symb::sym(  alphaName + nth_joint);

            this->offset.value     = set_Length(offset);
            this->theta.value      = set_Angle(theta, is_Rad);
            this->linkLen.value    = set_Length(linkLen);
            this->alpha.value      = set_Angle(alpha, is_Rad);
        };

        /*!
          Check the validity of the angle
        */
        double DH_Parameter::set_Angle(double angle, bool is_Rad)
        {
            if(!is_Rad) /**< Check if angle in rad, if not convert angle from degree to rad. */
            {
                angle = angle * M_PI / 180;
            }

            /*!
              Rotation angle and Twist angle have range[-PI, PI]
              Unfortunately, double numbers may vary with small values such as 1e-10
              Define Epsilon variable to hold small number as toleration
            */
            double Epsilon = 1e-10;
            if (abs(angle) > (M_PI + Epsilon))
            {
                std::cout << "Rotation angle and/or Twist angle have/has range[-PI, PI]" << std::endl;
                std::cout << "Correct your input" << std::endl;
                exit (EXIT_FAILURE);
            }
            else
            {
                return angle;
            }
        }

        /*!
          Check the validity of the length
        */
        double DH_Parameter::set_Length(double length)
        {
             /*!
              link length and joint offset are positive numbers
              Unfortunately, double numbers may vary with small values such as 1e-10
              Define Epsilon variable to hold small number as toleration
            */
            double Epsilon = -1e-10;
            if (length < Epsilon)
            {
                std::cout << "Link length and/or Joint offset are/is positive not negative number/s" << std::endl;
                std::cout << "Correct your input" << std::endl;
                exit (EXIT_FAILURE);
            }
            else
            {
                return length;
            }
        }

        /*!
          Calculate Moment of inertia tensor, assuming solid cylinder
        */
        Matrix3Symb MOIT(LinkProp& prop, Symb& linkLength)
        {
            Matrix3Symb MOI = SpecialMatrix::Zero<Symb, 3, 3>();
            /*!
              solid cylinder has only diagonal elements the rest are zeros
            */
            MOI(0, 0) = (1.0/12) * prop.mass.symbolic * (3*Symb::sq(prop.radius.symbolic) + Symb::sq(linkLength)); /**< sq() returns the square of the input. */
            MOI(1, 1) = (1.0/12) * prop.mass.symbolic * (3*Symb::sq(prop.radius.symbolic) + Symb::sq(linkLength)); /**< sq() attention 1/12 = 0 , 1.0/12 = 0.0833. */
            MOI(2, 2) = (1.0/2)  * prop.mass.symbolic * Symb::sq(prop.radius.symbolic);

            return MOI;
        }

        /*!
          Calculate skew matrix version of a given vector
        */
        Matrix3Symb Skew(const Vector3Symb& r)
        {
             /*!
               Create a (3 x 3) matrix to hold members of skew matrix
               r(0) -> i, r(1) -> j, r(2) -> k
             */
            Matrix3Symb r_skew;
            r_skew << 0   , -r(2),  r(1),
                      r(2),  0   , -r(0),
                     -r(1),  r(0),   0  ;

            return r_skew;
        }

        /*!
          Calculate vector version of a given skew matrix
        */
        Vector3Symb SkewToVector(const Matrix3Symb& mat)
        {
            Vector3Symb r; /**< Instantiate a column vector. */
            /*
            if(!CheckSkew(mat)) /**< check the validity of the given skew matrix.
            {
                std::cout << "The input is not a valid skew matrix" << std::endl;
                r << 0, 0, 0;
                return r;
            }
            */

            r << mat(2, 1), mat(0, 2), mat(1, 0);

            return r;
        }

        /*!
          Check if the matrix is a valid skew matrix
        */
        bool CheckSkew(const Matrix3Symb& mat)
        {
            /*!
              Skew matrix has 2 conditions
                    1. diagonal entries are zeros
                    2. mat(i, j), -1*mat(j, i)
              mathematically speaking -> A matrix is
              skew-symmetric if and only if it is the
              opposite of its transpose.
            */
            mat.transpose();
            int mat_size = 3;

            for(int i = 0; i < mat_size; i++) /**< if condition one failed return false. */
            {
                if(!Symb::is_equal(mat(i, i), 0))
                {
                    return false;
                }
            }
            /**< if condition two failed return false. */
            for(int i = 0; i < mat_size; i++)
            {
                for(int j = 0; j < mat_size; j++)
                {
                    if(!Symb::is_equal(mat(i, j), -1*mat(j, i), 1) && i != j)
                    {
                        return false;
                    }
                }
            }
            /**< return true if both conditions passed. */
            return true;
        }

        /*!
          Calculate outer product
        */
        Matrix3Symb OuterProduct(const Vector3Symb& r1, const Vector3Symb& r2)
        {
            return r1 * r2.transpose();
        }

        /*!
          Calculate rotation matrix
          Given:
                1. axis of rotation
                2. angleSymb (reference of the symbolic theta)
        */
        Matrix4Symb RotSymb(const Vector3Symb& r, Symb& angleSymb)
        {
            /*!
              Create a (4 x 4) matrix to hold members of rotation matrix
            */
            Matrix4Symb rotationM = SpecialMatrix::Identity<Symb, 4>(); /**< return true if both conditions passed. */

            /*!
              Rodrigues' rotation formula
            */
            rotationM.topLeftCorner(3, 3) = cos(angleSymb) * SpecialMatrix::Identity<Symb, 3>() + sin(angleSymb) * Skew(r) + (1 - cos(angleSymb)) * OuterProduct(r, r);

            return rotationM;

        }

        Matrix4Symb TransSymb(const Vector3Symb& r, Symb& distanceSymb, bool is_CPOL)
        {
            /*!
              Create a (4 x 4) matrix to hold members of Transition matrix
            */
            Matrix4Symb TransM = SpecialMatrix::Identity<Symb, 4>(); /**< Initialize the matrix. */

            /*!
              The distance traveled along specific vector
              The column vector of size (3 x 1) shows up at the top right corner of the Trans matrix
            */
            TransM.topRightCorner(3, 1) = distanceSymb * r;
            if(is_CPOL)
            {
                TransM.topRightCorner(3, 1) *= 0.5;
            }

            return TransM;
        }

        Matrix4Symb HomogeneousMatrix(DH_Parameter& para, bool is_CPOL)
        {

            return TransSymb(CartesianAxis::z, para.offset.symbolic) * RotSymb(CartesianAxis::z, para.theta.symbolic)
                     * TransSymb(CartesianAxis::x, para.linkLen.symbolic, is_CPOL) * RotSymb(CartesianAxis::x, para.alpha.symbolic);
        }

        void StringToVector(Vector3Real& vector3D, std::string& str)
        {
            char delimiter = ' ';
            int startIndex = 0;
            int vectorIndex = 0;
            int sz = str.size();

            int i = 0;
            while(i <= sz)
            {
                if(str[i] == delimiter || i == sz)
                {
                    vector3D(vectorIndex) = std::stod(str.substr(startIndex));
                    startIndex = i;
                    vectorIndex++;
                }

                i++;
            }

        }


        /* ********************************************************************** *
         * Rot Function calculates the change in rotation matrix, thus rotation   *
         * matrix should be at least initialized.                                 *
         * RotRPY: Rotation Roll Pitch Yaw                                        *
         * ********************************************************************** */
        void RotRPY(Matrix4Symb& rotMat, double roll, double pitch, double yaw)
        {
            Matrix3Symb temp;

            if(roll != 0)
            {
                temp <<  1          ,  0        ,  0         ,
                         0          ,  cos(roll), -sin(roll) ,
                         0          ,  sin(roll),  cos(roll) ;

                rotMat.topLeftCorner(3, 3) *= temp;
            }
            if(pitch != 0)
            {
                temp <<  cos(pitch) ,  0         ,  sin(pitch) ,
                         0          ,  1         ,  0          ,
                        -sin(pitch) ,  0         ,  cos(pitch) ;

                rotMat.topLeftCorner(3, 3) *= temp;
            }
            if(yaw != 0)
            {
                temp <<  cos(yaw) , -sin(yaw),  0           ,
                         sin(yaw) ,  cos(yaw),  0           ,
                         0        ,  0       ,  1           ;

                rotMat.topLeftCorner(3, 3) *= temp;

            }


        }

        void Roll(Matrix4Symb& rotMat, Symb& angle)
        {
            Matrix3Symb temp;
            temp <<  1          ,  0       ,  0        ,
                     0          ,  cos(angle), -sin(angle) ,
                     0          ,  sin(angle),  cos(angle) ;

            rotMat.topLeftCorner(3, 3) *= temp;

        }

        void Pitch(Matrix4Symb& rotMat, Symb& angle)
        {
            Matrix3Symb temp;
            temp <<  cos(angle) ,  0         ,  sin(angle) ,
                     0          ,  1         ,  0          ,
                    -sin(angle) ,  0         ,  cos(angle) ;

            rotMat.topLeftCorner(3, 3) *= temp;

        }

        void Yaw(Matrix4Symb& rotMat, Symb& angle)
        {
            Matrix3Symb temp;
            temp <<  cos(angle) , -sin(angle),  0           ,
                     sin(angle) ,  cos(angle),  0           ,
                     0          ,  0         ,  1           ;

            rotMat.topLeftCorner(3, 3) *= temp;

        }

    }

}
