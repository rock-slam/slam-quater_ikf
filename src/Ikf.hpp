/**\file Ikf.hpp
 *
 * This class has the primitive methods for an Indirect Kalman Filter implementation
 * for an Attitude and Heading Reference System - AHRS. The filter is Quaternion
 * based using accelerometers, gyroscopes and magnetometers. The filter performs the
 * prediction step based on the gyroscopes and therefore quaternion integration.
 * The measurement is formed by two step. First measurement step uses the accelerometers
 * in order to correct the pitch and roll angles. Second measurement step uses the
 * magnetometers only for the yaw angle. The first one estimates external acceleration
 * and compensate it increasing the measurement noise matrix.
 *
 * This indirect Kalman filter is based on the paper:  Young Soo Suh, Member, IEEE
 * "Orientation estimation using a quaternion-based indirect Kalman filter with adaptive estimation of external acceleration"
 *
 * @author Javier Hidalgo Carrio | DFKI RIC Bremen | javier.hidalgo_carrio@dfki.de
 * @date January 2014.
 * @version 2.0.
 */

#ifndef INDIRECT_KALMAN_FILTER_HPP
#define INDIRECT_KALMAN_FILTER_HPP

#include <Eigen/Geometry> /** Quaternion, Angle-Axis,...*/
#include <Eigen/StdVector> /** STL container with Eigen types */
#include <Eigen/LU> /** Lineal algebra of Eigen */
#include <Eigen/SVD> /** Singular Value Decomposition (SVD) of Eigen */
#include <Eigen/Eigenvalues> /** Eigen eigenvalues and eigenvectors **/


#include <vector> /** std_vector **/

/** Boost **/
#include <boost/shared_ptr.hpp> /** For shared pointers **/

/** Adaptive update **/
#include "AdaptiveAttitudeCov.hpp"

//#define INDIRECT_KALMAN_FILTER_DEBUG_PRINTS 1

namespace filter
{
    using namespace Eigen;

    template <typename _Scalar, bool _Accelerometers, bool _Inclinometers>
    class Ikf
    {

    public:
        /** Convert bool to int is standard in C++ **/
        static const int flagAcc = static_cast<int>(_Accelerometers);
        static const int flagIncl = static_cast<int>(_Inclinometers);

        /** Constants definitions **/
        static const unsigned int QUATERNION_SIZE = 4;
        static const unsigned int IKFSTATEVECTORSIZE = 9 + (flagAcc*3) + (flagIncl*3);

    private:
        /** WGS-84 ellipsoid constants (Nominal Gravity Model and Earth angular velocity) **/
        static const int Re = 6378137; /** Equatorial radius in meters **/
        static const int Rp = 6378137; /** Polar radius in meters **/
        static const double ECC = 0.0818191908426; /** First eccentricity **/
        static const double GRAVITY = 9.79766542; /** Mean value of gravity value in m/s^2 **/
        static const double GWGS0 = 9.7803267714; /** Gravity value at the equator in m/s^2 **/
        static const double GWGS1 = 0.00193185138639; /** Gravity formula constant **/
        static const double EARTHW = 7.292115e-05; /** Earth angular velocity in rad/s **/

        /** Magnetic declination **/
        enum DECLINATION_CONSTS {
            EAST = 1, /** EAST is 1 and means positive magnetic declination **/
            WEST = 2 /** WEST is 2 and means negative magnetic declination **/
        };

        /** Data Types for the Adaptive measurements **/
        typedef AdaptiveAttitudeCov<_Scalar, IKFSTATEVECTORSIZE, 3> AdaptiveAttitudeCovAcc;
        typedef AdaptiveAttitudeCov<_Scalar, IKFSTATEVECTORSIZE, 3> AdaptiveAttitudeCovIncl;

        /** Filter members */
        Eigen::Matrix <_Scalar,IKFSTATEVECTORSIZE,1> x; /** State vector */
        Eigen::Matrix <_Scalar,3,1> gtilde; /** gravitation acceleration */
        Eigen::Matrix <_Scalar,3,1> mtilde; /** Magnetic dip angle */
        Eigen::Quaternion <_Scalar> q4;  /** Attitude quaternion */
        Eigen::Matrix <_Scalar,QUATERNION_SIZE,QUATERNION_SIZE> oldomega4; /** Quaternion integration matrix */
        Eigen::Matrix <_Scalar,IKFSTATEVECTORSIZE,IKFSTATEVECTORSIZE> P; /** Error covariance matrix */
        Eigen::Matrix <_Scalar,IKFSTATEVECTORSIZE,IKFSTATEVECTORSIZE> A; /** System matrix */
        Eigen::Matrix <_Scalar,IKFSTATEVECTORSIZE,IKFSTATEVECTORSIZE> Q; /** Process noise covariance matrix */
        Eigen::Matrix <_Scalar,3,3> Ra; /** Measurement noise covariance matrix for acc */
        Eigen::Matrix <_Scalar,3,3> Qg; /** Measurement noise covariance matrix for gyros */
        Eigen::Matrix <_Scalar,3,3> Rm; /** Measurement noise covariance matrix for mag */
        Eigen::Matrix <_Scalar,3,3> Ri; /** Measurement noise covariance matrix for mag */
        Eigen::Matrix <_Scalar,3,IKFSTATEVECTORSIZE> H1; /** Measurement 1 Observation matrix */
        Eigen::Matrix <_Scalar,3,IKFSTATEVECTORSIZE> H2; /** Measurement 2 Observation matrix */
        Eigen::Matrix <_Scalar,3,IKFSTATEVECTORSIZE> H3; /** Measurement 3 Observation matrix */
        Eigen::Matrix <_Scalar,3,1> bghat; /** Estimated bias instability for gyroscope */
        Eigen::Matrix <_Scalar,3,1> kghat; /** Estimated rate random walk for gyroscope */
        Eigen::Matrix <_Scalar,3,1> kahat; /** Estimated acceleration random walk for accelerometers */
        Eigen::Matrix <_Scalar,3,1> kihat; /** Estimated acceleration random walk for inclinometers */

        /** Object of Class for Adaptive Measurement of Attitude Covariance Matrix (Accelerometers) **/
        boost::shared_ptr<AdaptiveAttitudeCovAcc> adapAttAcc;

        /** Object of Class for Adaptive Measurement of Attitude Covariance Matrix (Inclinometers) **/
        boost::shared_ptr<AdaptiveAttitudeCovIncl> adapAttIncl;

    protected:
        void initAdaptiveAttitude(const unsigned int accM1, const unsigned int accM2, const double accGamma,
                                const unsigned int incM1, const unsigned int incM2, const double incGamma)
        {
            this->adapAttAcc.reset(new AdaptiveAttitudeCovAcc (accM1, accM2, accGamma));
            this->adapAttIncl.reset(new AdaptiveAttitudeCovIncl (incM1, incM2, incGamma));
        }


    public:

        /** Default Destructor
        */
        ~Ikf()
        {
            adapAttAcc.reset();
            adapAttIncl.reset();
        }

        /**
        * @brief This function Initialize the vectors and matrix of the IKF
        *
        * This method receives the measurement noise matrix of the sensors
        * The theoretical gravity value and the Dip angle of the location.
        *
        * @author Javier Hidalgo Carrio.
        *
        * @param[in] P_0 Initial state covariance matrix
        * @param[in] Ra measurement noise matrix of Accelerometers.
        * @param[in] Qg covariance noise matrix of Gyroscopes.
        * @param[in] Rm measurement noise matrix of Magnetometers.
        * @param[in] Qbg covariance noise matrix of the gyroscopes bias
        * @param[in] Qba covariance noise matrix of the accelerometers bias 
        * @param[in] g local gravitational value.
        * @param[in] alpha Dip angle
        *
        * @return void
        *
        */
        void Init(const Eigen::Matrix <_Scalar,Ikf::IKFSTATEVECTORSIZE,Ikf::IKFSTATEVECTORSIZE> &P_0,
                const Eigen::Matrix <_Scalar,3,3> &Ra,
                const Eigen::Matrix <_Scalar, 3, 3> Qg,
                const Eigen::Matrix <_Scalar,3,3> &Rm,
                const Eigen::Matrix <_Scalar,3,3> &Ri,
                const Eigen::Matrix <_Scalar, 3, 3> &Qbg,
                const Eigen::Matrix <_Scalar,3,3> &Qkg,
                const Eigen::Matrix <_Scalar,3,3> &Qka,
                const Eigen::Matrix <_Scalar,3,3> &Qki,
                const Eigen::Matrix <_Scalar,3,1> &taub,
                double g, double alpha,
                unsigned int am1, unsigned int am2, double agamma,
                unsigned int im1, unsigned int im2, double igamma)
        {
            /** Gravitation acceleration **/
            gtilde << 0, 0, g;

            /** Dip angle (alpha is in rad) **/
            mtilde(0) = cos(alpha);
            mtilde(1) = 0;
            mtilde(2) = -sin(alpha);


            /** Kalman filter state, error covariance and process noise covariance **/
            x = Eigen::Matrix <_Scalar,Ikf::IKFSTATEVECTORSIZE,1>::Zero();

            Q = Eigen::Matrix <_Scalar,Ikf::IKFSTATEVECTORSIZE,Ikf::IKFSTATEVECTORSIZE>::Zero();
            Q.template block<3, 3> (0,0) = 0.25 * (Qg);
            Q.template block<3, 3> (3,3) = (Qbg);
            Q.template block<3, 3> (6,6) = (Qkg);
            if (_Accelerometers)
                Q.template block<3, 3> (9,9) = (Qka);
            if (_Inclinometers)
                Q.template block<3, 3> (9+(flagAcc*3),6+(flagAcc*3)) = (Qki);

            /** Initial error covariance **/
            P = (P_0);

            /** Static part of the observation matrices **/
            H1 = Eigen::Matrix<_Scalar, 3,Ikf::IKFSTATEVECTORSIZE>::Zero();
            H2 = Eigen::Matrix<_Scalar, 3,Ikf::IKFSTATEVECTORSIZE>::Zero();
            H3 = Eigen::Matrix<_Scalar, 3,Ikf::IKFSTATEVECTORSIZE>::Zero();

            if (_Accelerometers)
            {
                H1(0,9) = 1; H1(1,10) = 1; H1(2,11) = 1;
            }
            if (_Inclinometers)
            {
                H3(0,9) = 1; H3(1,10) = 1; H3(2,11) = 1;
            }
            if (_Accelerometers && _Inclinometers)
            {
                H3 = Eigen::Matrix<_Scalar, 3,Ikf::IKFSTATEVECTORSIZE>::Zero();
                H3(0,12) = 1; H3(1,13) = 1; H3(2,14) = 1;
            }

            /** System matrix A **/
            A = Eigen::Matrix <_Scalar,Ikf::IKFSTATEVECTORSIZE,Ikf::IKFSTATEVECTORSIZE>::Zero();
            A(0,3) = -0.5;A(1,4) = -0.5;A(2,5) = -0.5;
            A(0,6) = -0.5;A(1,7) = -0.5;A(2,8) = -0.5;
            A(3,3) = -(1.0/taub[0]);A(4,4) = -(1.0/taub[1]);A(5,5) = -(1.0/taub[2]);

            /** Initial bias **/
            bghat = Eigen::Matrix <_Scalar,3,1>::Zero();
            kghat = Eigen::Matrix <_Scalar,3,1>::Zero();
            kahat = Eigen::Matrix <_Scalar,3,1>::Zero();
            kihat = Eigen::Matrix <_Scalar,3,1>::Zero();

            /** Default omega matrix **/
            oldomega4 << 0 , 0 , 0 , 0,
              0 , 0 , 0 , 0,
              0 , 0 , 0 , 0,
              0 , 0 , 0 , 0;

            /** Initial quaternion in Init**/
            q4.w() = 1.00;
            q4.x() = 0.00;
            q4.y() = 0.00;
            q4.z() = 0.00;

            /** Fill matrix Qg, Ra, Rm and Ri **/
            this->Ikf::Ra = Ra;
            this->Ikf::Qg = Qg;
            this->Ikf::Rm = Rm;
            this->Ikf::Ri = Ri;

            /** Initialize adaptive object **/
            initAdaptiveAttitude(am1, am2, agamma, im1, im2, igamma);

            /** Print filter information **/
            #ifdef INDIRECT_KALMAN_FILTER_DEBUG_PRINTS
            std::cout<< "P:\n"<<P<<"\n";
            std::cout<< "Q:\n"<<Q<<"\n";
            std::cout<< "H1:\n"<<H1<<"\n";
            std::cout<< "H2:\n"<<H2<<"\n";
            std::cout<< "H3:\n"<<H3<<"\n";
            std::cout<< "A:\n"<<A<<"\n";
            std::cout<< "mtilde:\n"<<mtilde<<"\n";
            std::cout<< "gtilde:\n"<<gtilde<<"\n";
            std::cout<< "Ra:\n"<<Ra<<"\n";
            std::cout<< "Qg:\n"<<Qg<<"\n";
            std::cout<< "Rm:\n"<<Rm<<"\n";
            std::cout<< "Ri:\n"<<Ri<<"\n";
            #endif

            return;

        }

        /**
        * @brief Performs the prediction step of the filter.
        * 
        * It computes the discrete version of the matrix A to propagate forward
        * the state vector x. It computes the Q and P matrix as well as the 
        * quaternion integration from the input vector u and the delta time.
        *
        * @author Javier Hidalgo Carrio.
        *
        * @param[in] u vector with the angular velocity
        * @param[in] dt delta time between samples
        *
        * @return void
        *
        */
        void predict(Eigen::Matrix <_Scalar,3,1>  &u, double dt)
        {
            Eigen::Matrix <_Scalar,3,3> vec2product; /**< Vector 2 product matrix */
            Eigen::Matrix <_Scalar,3,1> angvelo; /**< Vec 2 product  matrix */
            Eigen::Matrix <_Scalar,QUATERNION_SIZE,QUATERNION_SIZE> omega4; /**< Quaternion integration matrix */
            Eigen::Matrix <_Scalar,QUATERNION_SIZE,1> quat; /**< Quaternion integration matrix */
            Eigen::Matrix <_Scalar,IKFSTATEVECTORSIZE,IKFSTATEVECTORSIZE> dA; /**< Discrete System matrix */
            Eigen::Matrix <_Scalar,IKFSTATEVECTORSIZE,IKFSTATEVECTORSIZE> Qd; /**< Discrete Q matrix */

            /** Compute the vector2product matrix with the angular velocity **/
            angvelo = (u) - bghat - kghat; /** Eliminate the Bias instability and rate random walk**/

            vec2product << 0, -angvelo(2), angvelo(1),
                        angvelo(2), 0, -angvelo(0),
                        -angvelo(1), angvelo(0), 0;

            /** Compute the dA Matrix **/
            A.template block<3, 3> (0,0) = -vec2product;
            dA = Eigen::Matrix<_Scalar,IKFSTATEVECTORSIZE,IKFSTATEVECTORSIZE>::Identity() + A * dt + 0.5 *A * A * pow(dt,2);

            /** Propagate the vector through the system **/
            x = dA * x;
            Qd = Q*dt + 0.5*dt*dt*A*Q + 0.5*dt*dt*Q*A.transpose();
            Qd = 0.5*(Qd + Qd.transpose());//Guarantee symmetry
            P = dA*P*dA.transpose() + Qd;

            omega4 << 0,-angvelo(0), -angvelo(1), -angvelo(2),
                    angvelo(0), 0, angvelo(2), -angvelo(1),
                    angvelo(1), -angvelo(2), 0, angvelo(0),
                    angvelo(2), angvelo(1), -angvelo(0), 0;

            quat(0) = q4.w();
            quat(1) = q4.x();
            quat(2) = q4.y();
            quat(3) = q4.z();

            /** Third-order gyroscopes integration accuracy **/
            quat = (Eigen::Matrix<_Scalar,QUATERNION_SIZE,QUATERNION_SIZE>::Identity() +(0.75 * omega4 *dt)-(0.25 * oldomega4 * dt) -
            ((1.0/6.0) * angvelo.squaredNorm() * pow(dt,2) *  Eigen::Matrix<_Scalar,QUATERNION_SIZE,QUATERNION_SIZE>::Identity()) -
            ((1.0/24.0) * omega4 * oldomega4 * pow(dt,2)) - ((1.0/48.0) * angvelo.squaredNorm() * omega4 * pow(dt,3))) * quat;

            q4.w() = quat(0);
            q4.x() = quat(1);
            q4.y() = quat(2);
            q4.z() = quat(3);
            q4.normalize();

            oldomega4 = omega4;

            return;

        }

        void update(const Eigen::Matrix <_Scalar,3,1>  &acc, bool acc_on)
        {
             if (_Accelerometers)
             {
                 Eigen::Matrix <_Scalar,3,1>  aux; aux.setZero();
                 update(acc, acc_on, aux, false, aux, false);
             }
             else
                throw std::runtime_error("IKF without accelerometers and you are updating with accelerometers");
        }

        void update(const Eigen::Matrix <_Scalar,3,1>  &acc, bool acc_on,
                const Eigen::Matrix< _Scalar,3, 1> &measurement, bool measurement_on)
        {
            Eigen::Matrix <_Scalar,3,1>  aux; aux.setZero();
            if (_Inclinometers)
            {
                update(acc, acc_on, measurement, measurement_on, aux, false);
            }
            else
            {
                update(acc, acc_on, aux, false, measurement, measurement_on);
            }
        }


        void update(const Eigen::Matrix <_Scalar,3,1>  &measurement1,
                const Eigen::Matrix< _Scalar,3, 1> &measurement2)
        {
            Eigen::Matrix <_Scalar,3,1>  aux; aux.setZero();
            if (_Accelerometers && _Inclinometers)
            {
                update(measurement1, true, measurement2, true, aux, false);
            }
            else if (_Accelerometers)
            {
                update(measurement1, true, aux, false, measurement2, true);
            }
            else if (_Inclinometers)
            {
                update(aux, false, measurement1, true, measurement2, true);
            }
            else
                throw std::runtime_error("IKF without accelerometers and inclinometers. Only update with magnetometers allowed");
        }


        /**
        * @brief Performs the measurement and correction steps of the filter.
        * 
        * The IKf is based on two measurement step:\
        *  1. Measurement step to correct Pitch and Roll from accelerometers.
        *  2. Measurement step to correct Pitch and Roll from inclinometers.
        *  3. Measurement step to correct Yaw angle from magnetometers.
        * 
        * The first measurement step is dynamics. The noise covariamce matrix
        * of the update is dynamic depending on external accelerations felt on
        * the accelerometers. That means the variance noise increase or decrease
        * depending on the external acceleration. Thas is the main different between 
        * normal EKF.
        * 
        * The second measurement step only affects the Yaw (heading) angle.
        *
        * @author Javier Hidalgo Carrio.
        *
        * @param[in] acc pointer to vector with accelerations
        * @param[in] magn pointer to vector with magnetometers
        * @param[in] magn_on boolean value to connect or disconnect the magnetometers correction
        * @return void
        *
        */
        void update(const Eigen::Matrix <_Scalar,3,1>  &acc, bool acc_on,
                const Eigen::Matrix< _Scalar,3, 1> &incl, bool incl_on,
                const Eigen::Matrix <_Scalar,3, 1> &mag, bool magn_on)
        {
            Eigen::Matrix <_Scalar,3,3> vec2product; /**< Vector 2 product  matrix */
            Eigen::Matrix <_Scalar,3,3> fooR2; /**<  Measurement noise matrix from accelerometers matrix Ra*/
            Eigen::Matrix <_Scalar,IKFSTATEVECTORSIZE,IKFSTATEVECTORSIZE> P1; /**< Error covariance matrix for measurement 1*/
            Eigen::Matrix <_Scalar,IKFSTATEVECTORSIZE,IKFSTATEVECTORSIZE> P2; /**< Error covariance matrix for measurement 2*/
            Eigen::Matrix <_Scalar,IKFSTATEVECTORSIZE,IKFSTATEVECTORSIZE> P3; /**< Error covariance matrix for measurement 3*/
            Eigen::Matrix <_Scalar,IKFSTATEVECTORSIZE,IKFSTATEVECTORSIZE> auxM; /**< Auxiliar matrix for computing Kalman gain in measurement*/
            Eigen::Matrix <_Scalar,IKFSTATEVECTORSIZE, 3> K1; /**< Kalman Gain matrix for measurement 1*/
            Eigen::Matrix <_Scalar,IKFSTATEVECTORSIZE, 3> K2; /**< Kalman Gain matrix for measurement 2*/
            Eigen::Matrix <_Scalar,IKFSTATEVECTORSIZE, 3> K3; /**< Kalman Gain matrix for measurement 3*/
            Eigen::Matrix <_Scalar,3,3> R1; /**< Acceleration covariance matrix */
            Eigen::Matrix <_Scalar,3,3> R3; /**< Inclinometers covariance matrix */
            Eigen::Quaternion <_Scalar> qe;  /**< Attitude error quaternion */
            Eigen::Matrix <_Scalar,3,1> gtilde_body; /**< Gravitation in the body frame */
            Eigen::Matrix <_Scalar,3,1> mtilde_body; /**< Magnetic field in the body frame */
            Eigen::Matrix <_Scalar,3,1> z1; /**< Measurement vector 1 Acc */
            Eigen::Matrix <_Scalar,3,1> z2; /**< Measurement vector 2 Mag */
            Eigen::Matrix <_Scalar,3,1> z3; /**< Measurement vector 3 Incl */
            Eigen::Matrix <_Scalar,3,1> auxvector; /**< Auxiliar vector variable */


            /**----------------------- **/
            /** Measurement step 1 Acc **/
            /**----------------------- **/

            if (acc_on && _Accelerometers)
            {
                /** First measurement step (Pitch and Roll correction from Acc) **/
                gtilde_body = q4.inverse() * gtilde;
                vec2product << 0, -gtilde_body(2), gtilde_body(1),
                            gtilde_body(2), 0, -gtilde_body(0),
                            -gtilde_body(1), gtilde_body(0), 0;

                H1.template block<3, 3> (0,0) = 2*vec2product;

                /** Measurement **/
                z1 = (acc) - kahat - gtilde_body;

                #ifdef INDIRECT_KALMAN_FILTER_DEBUG_PRINTS
                std::cout<<"acc:\n"<<acc<<"\n";
                std::cout<<"z1:\n"<<z1<<"\n";
                std::cout<<"g_body:\n"<<gtilde_body<<"\n";
                #endif

                /** The adaptive algorithm **/
                R1 = adapAttAcc->matrix (x, P, z1, H1, Ra);

                /** Compute the Kalman Gain Matrix **/
                P1 = P;
                Eigen::Matrix<_Scalar, 3, 3> S1, S1_inverse;
                S1 = H1 * P1 * H1.transpose() + R1;
                S1_inverse = S1.inverse();
                K1 = P1 * H1.transpose() * S1_inverse;
                Eigen::Matrix<_Scalar, 3, 1> innovationAcc = (z1 - H1 * x);

                /** Update the state vector and the covariance matrix **/
                x = x + K1 * innovationAcc;
                P = (Eigen::Matrix<_Scalar,IKFSTATEVECTORSIZE,IKFSTATEVECTORSIZE>::Identity()
                      -K1*H1)*P*(Eigen::Matrix<_Scalar,IKFSTATEVECTORSIZE,IKFSTATEVECTORSIZE>::Identity()
                      -K1*H1).transpose() + K1*R1*K1.transpose();

                #ifdef INDIRECT_KALMAN_FILTER_DEBUG_PRINTS
                std::cout<<"x(k+1|k+1):\n"<<x<<"\n";
                std::cout<<"P(k+1|k+1):\n"<<P<<"\n";
                std::cout<<"innovation:\n"<<innovationAcc<<"\n";
                std::cout<<"K1:\n"<<K1<<"\n";
                std::cout<<"R1:\n"<<R1<<"\n";
                #endif

                /** Update the quaternion with the Indirect approach **/
                /** This is necessary mainly because after(in the 2 measurement) C(q) is computed **/
                qe.w() = 1;
                qe.x() = x(0);
                qe.y() = x(1);
                qe.z() = x(2);
                q4 = q4 * qe;

                /** Normalize quaternion **/
                q4.normalize();

                /** Reset the quaternion part of the state vector **/
                x.template block<3,1>(0,0) = Eigen::Matrix<_Scalar, 3, 1>::Zero();
            }

            /**------------------------- **/
            /** Measurement step 2 Mag   **/
            /** It only updates Yaw angle**/
            /**------------------------- **/

            if (magn_on)
            {

                /** Second measurement step **/
                mtilde_body = q4.inverse() * mtilde;
                vec2product << 0, -mtilde_body(2), mtilde_body(1),
                            mtilde_body(2), 0, -mtilde_body(0),
                            -mtilde_body(1), mtilde_body(0), 0;

                /** Observation matrix **/
                H2.template block<3, 3> (0,0) = 2*vec2product;

                /** Measurement vector **/
                z2 = (mag) - mtilde_body;

                P2 = Eigen::Matrix<_Scalar, IKFSTATEVECTORSIZE, IKFSTATEVECTORSIZE>::Zero();
                P2.template block<3, 3>(0,0) = P.template block<3, 3>(0,0);

                auxvector << 0, 0, 1;
                auxvector = q4.inverse() * auxvector;

                /** Compute Kalman Gain **/
                auxM = Eigen::Matrix<_Scalar, IKFSTATEVECTORSIZE, IKFSTATEVECTORSIZE>::Zero();
                auxM.template block<3, 3>(0,0) = auxvector * auxvector.transpose();
                K2 = auxM * P2 * H2.transpose() * (H2*P2*H2.transpose() + Rm).inverse();

                /** Update the state vector and the covariance matrix **/
                x = x + K2*(z2 - (H2*x));
                P = P - K2 * H2 * P - P * H2.transpose() * K2.transpose() + K2*(H2*P*H2.transpose() + Rm)*K2.transpose();

                /** Update the quaternion with the Indirect approach **/
                qe.w() = 1;
                qe.x() = x(0);
                qe.y() = x(1);
                qe.z() = x(2);
                q4 = q4 * qe;

                /** Normalize quaternion **/
                q4.normalize();

                /** Reset the quaternion part of the state vector **/
                x.template block<3,1>(0,0) = Eigen::Matrix<_Scalar, 3, 1>::Zero();
            }

            if (incl_on && _Inclinometers)
            {
                /**--------------------------------- **/
                /** Measurement step 3 Inclinometers **/
                /**--------------------------------- **/

                /** Measurement step (Pitch and Roll correction from Inclinometers) **/
                gtilde_body = q4.inverse() * gtilde;
                vec2product << 0, -gtilde_body(2), gtilde_body(1),
                            gtilde_body(2), 0, -gtilde_body(0),
                            -gtilde_body(1), gtilde_body(0), 0;

                H3.template block<3, 3> (0,0) = 2*vec2product;

                /** Measurement **/
                z3 = (incl) - kihat - gtilde_body;

                #ifdef INDIRECT_KALMAN_FILTER_DEBUG_PRINTS
                std::cout<<"incl:\n"<<incl<<"\n";
                std::cout<<"z3:\n"<<z3<<"\n";
                std::cout<<"g_body:\n"<<gtilde_body<<"\n";
                #endif

                /** The adaptive algorithm **/
                R3 = adapAttIncl->matrix (x, P, z3, H3, Ri);

                /** Compute the Kalman Gain Matrix **/
                P3 = P;
                Eigen::Matrix<_Scalar, 3, 3> S3, S3_inverse;
                S3 = H3 * P3 * H3.transpose() + R3;
                S3_inverse = S3.inverse();
                K3 = P3 * H3.transpose() * S3_inverse;
                Eigen::Matrix<_Scalar, 3, 1> innovationIncl = (z3 - H3 * x);

                /** Update the state vector and the covariance matrix **/
                x = x + K3 * innovationIncl;
                P = (Eigen::Matrix<_Scalar,IKFSTATEVECTORSIZE,IKFSTATEVECTORSIZE>::Identity()
                      -K3*H3)*P*(Eigen::Matrix<_Scalar,IKFSTATEVECTORSIZE,IKFSTATEVECTORSIZE>::Identity()
                      -K3*H3).transpose() + K3*R3*K3.transpose();

                #ifdef INDIRECT_KALMAN_FILTER_DEBUG_PRINTS
                std::cout<<"x(k+1|k+1):\n"<<x<<"\n";
                std::cout<<"P(k+1|k+1):\n"<<P<<"\n";
                std::cout<<"innovation:\n"<<innovationIncl<<"\n";
                std::cout<<"K3:\n"<<K3<<"\n";
                std::cout<<"R3:\n"<<R3<<"\n";
                #endif

                /** Update the quaternion with the Indirect approach **/
                /** This is necessary mainly because after(in the 2 measurement) C(q) is computed **/
                qe.w() = 1;
                qe.x() = x(0);
                qe.y() = x(1);
                qe.z() = x(2);
                q4 = q4 * qe;

                /** Normalize quaternion **/
                q4.normalize();

                /** Reset the quaternion part of the state vector **/
                x.template block<3,1>(0,0) = Eigen::Matrix<_Scalar, 3, 1>::Zero();
            }

            /**---------------------------- **/
            /** Reset the rest of the state **/
            /**---------------------------- **/
            bghat = bghat + x.template block<3, 1> (3,0);
            x.template block<3, 1> (3,0) = Eigen::Matrix <_Scalar, 3, 1>::Zero();
            kghat = kghat + x.template block<3, 1> (6,0);
            x.template block<3, 1> (6,0) = Eigen::Matrix <_Scalar, 3, 1>::Zero();

            if (_Accelerometers)
            {
                kahat = kahat + x.template block<3, 1> (9,0);
                x.template block<3, 1> (6,0) = Eigen::Matrix <_Scalar, 3, 1>::Zero();
            }

            if (_Inclinometers)
            {
                kihat = kihat + x.template block<3, 1> (9+(flagAcc*3),0);
                x.template block<3, 1> (6+(flagAcc*3),0) = Eigen::Matrix <_Scalar, 3, 1>::Zero();
            }

            /** Guarantee SPD Covariance matrix **/
            P = this->guaranteeSPD< Eigen::Matrix <_Scalar,IKFSTATEVECTORSIZE,IKFSTATEVECTORSIZE> >(P);
            P = 0.5 * (P + P.transpose());//Guarantee symmetry

            #ifdef INDIRECT_KALMAN_FILTER_DEBUG_PRINTS
            std::cout<<"bghat:\n"<<bghat<<"\n";
            std::cout<<"kghat:\n"<<kghat<<"\n";
            std::cout<<"kahat:\n"<<kahat<<"\n";
            std::cout<<"kihat:\n"<<kihat<<"\n";
            #endif

            return;

        }

        /**
        * @brief This function Initialize Attitude
        * 
        * Initial orientation value beforeestart the IKF 
        *
        * @author Javier Hidalgo Carrio.
        *
        * @param[in] *initq pointer to quaternion with the initial orientation
        *
        * @return true if everything all right. false on other cases.
        *
        */
        bool setAttitude (const Eigen::Quaternion <double> &initq)
        {
            if (&initq != NULL)
            {
                /** Initial orientation **/
                q4 = initq;

                return true;
            }
            return false;
        }

        /**
        * @brief This function Initialize the State vector
        * 
        * The state vector is formed by 9 element.
        * (0-2) -> the vector part of a error quaternion
        * (3-5) -> gyroscope bias estimation
        * (6-8) -> accelerometer bias estimation
        *
        * @param[in] *x_0 a initial/desired state vector
        *
        * @return OK is everything all right. ERROR on other cases.
        *
        */
        void setState (const Eigen::Matrix <double,IKFSTATEVECTORSIZE,1> &x_0)
        {
            x = x_0;

            return;
        }

        /**
        * @brief This function set the initial Omega matrix
        * 
        * Initial Omega matrix with angular velocity for 
        * quaternion integration.
        *
        * @author Javier Hidalgo Carrio.
        *
        * @param[in] *u pointer to vector with the angular velocity
        *
        * @return true if everything all right. false on other cases.
        *
        */
        bool setOmega (const Eigen::Matrix <_Scalar, 3,1> &u)
        {
            if (&u != NULL)
            {
                /** Initialization for quaternion integration **/
                oldomega4 << 0,-(u)(0), -(u)(1), -(u)(2),
                    (u)(0), 0, (u)(2), -(u)(1),
                    (u)(1), -(u)(2), 0, (u)(0),
                    (u)(2), (u)(1), -(u)(0), 0;

                return true;
            }
            return false;
        }

        /**
        * @brief On/Off initial bias
        *
        * Initial On/Off sensor bias. Otherwise the default
        * in Init is set them to zero.
        *
        * @param[in] gbias vector with initial gyroscopes bias
        * @param[in] abias vector with initial accelerometers bias
        *
        * @return void.
        *
        */
        void setInitBias (const Eigen::Matrix<_Scalar, 3, 1> &gbiasinstability,
            const Eigen::Matrix<_Scalar, 3, 1> &gratewalk,
            const Eigen::Matrix<_Scalar, 3, 1> &aratewalk,
            const Eigen::Matrix<_Scalar, 3, 1> &iratewalk)
        {
            this->bghat = gbiasinstability;
            this->kghat = gratewalk;
            this->kahat = aratewalk;
            this->kihat = iratewalk;

            return;
        }

        /**
        * @brief Initial gravity
        *
        * Set initial gravity. If after initialization a new
        * theoretical or measured gravity is available use this method.
        * Note: always before start running the filter.
        *
        * @param[in] gravity initial gravity value
        *
        * @return void.
        *
        */
        void setGravity(const double gravity)
        {
            gtilde << 0.00, 0.00, gravity;
            return;
        }

        /**
        * @brief Set the filter covariance
        *
        * @param[in] Pk covariance matrix of appropriate dimension
        *
        * @return void.
        *
        */
        void setCovariance(const Eigen::Matrix< double, Ikf::IKFSTATEVECTORSIZE , Ikf::IKFSTATEVECTORSIZE> &Pk)
        {
            P = Pk;
            return;
        }

        /**
        * @brief Gets the current gyroscopes bias instability
        */
        Eigen::Matrix<double, 3, 1> getGyroBiasInstability()
        {
            return this->bghat;
        }

        /**
        * @brief Gets the current gyroscopes random walk
        */
        Eigen::Matrix<double, 3, 1> getGyroRateRandomWalk()
        {
            return this->kghat;
        }


        /**
        * @brief Gets the current accelerometers bias
        */
        Eigen::Matrix<double, 3, 1> getAccAccelerationRandomWalk()
        {
            return this->kahat;
        }

        /**
        * @brief Gets the current inclinometers bias
        */
        Eigen::Matrix<double, 3, 1> getInclAccelerationRandomWalk()
        {
            return this->kihat;
        }

        /**
        * @brief Gets the current state vector of the filter
        * 
        * @author Javier Hidalgo Carrio.
        *
        * @return State Vector
        *
        */
        Eigen::Matrix <double,IKFSTATEVECTORSIZE,1> getState()
        {
            return x;
        }

        /**
        * @brief Gets gravity in the IMU body frame
        */
        Eigen::Matrix<_Scalar, 3, 1> getGravityinBody()
        {
            return q4.inverse() * gtilde;
        }

        /**
        * @brief Gets gravity in the local Geographic frame
        */
        inline Eigen::Matrix<_Scalar, 3, 1> getGravity() const
        {
            return gtilde;
        }

        /**
        * @brief Gets the current orientation in Quaternion
        * 
        * @author Javier Hidalgo Carrio.
        *
        * @return Quaternion with the current orientation.
        *
        */
        Eigen::Quaternion <double> getAttitude()
        {
            return q4;
        }

        /**
        * @brief Gets Noise covariance matrix
        * 
        * @author Javier Hidalgo Carrio.
        *
        * @return Matrix P of the covariance of the state vector
        *
        */
        Eigen::Matrix <double,IKFSTATEVECTORSIZE,IKFSTATEVECTORSIZE> getCovariance()
        {
            return P;
        }

        /**
        * @brief Conversion Quaternion to DCM (Direct Cosine Matrix) (Alternative to Eigen)
        * 
        * Conversion to a transformation matrix from a quaternion
        * The quaternion is represented in Eigen convention:
        * w+xi+yj+zk, first element the scalar and others three are the vectorial part.
        * 
        * @author Javier Hidalgo Carrio.
        *
        * @param[in] *q pointer to a quaternion vector.
        * @param[out] *C pointer to a matrix. The three by three matrix
        *
        * @return void
        *
        */
        void Quaternion2DCM(Eigen::Quaternion< _Scalar >* q, Eigen::Matrix< _Scalar, 3, 3  >*C)
        {
            double q0, q1, q2, q3;

            if (C != NULL)
            {
                /** Take the parameters of the quaternion */
                q0 = q->w();
                q1 = q->x();
                q2 = q->y();
                q3 = q->z();

                /** Create the DCM matrix from the actual quaternion */
                (*C)(0,0) = 2 * q0 * q0 + 2 * q1 * q1 - 1;
                (*C)(0,1) = 2 * q1 * q2 + 2 * q0 * q3;
                (*C)(0,2) = 2 * q1 * q3 - 2 * q0 * q2;
                (*C)(1,0) = 2 * q1 * q2 - 2 * q0 * q3;
                (*C)(1,1) = 2 * q0 * q0 + 2 * q2 * q2 - 1;
                (*C)(1,2) = 2 * q2 * q3 + 2 * q0 * q1;
                (*C)(2,0) = 2 * q1 * q3 + 2 * q0 * q2;
                (*C)(2,1) = 2 * q2 * q3 - 2 * q0 * q1;
                (*C)(2,2) = 2 * q0 * q0 + 2 * q3 * q3 - 1;	
            }
            return;
        }

        /**
        * @brief Method to guarantee Semi-Positive Definite (SPD) matrix
        */
        template <typename _MatrixType>
        static _MatrixType guaranteeSPD (const _MatrixType &A)
        {
            _MatrixType spdA;
            Eigen::VectorXd s;
            s.resize(A.rows(), 1);

            /**
             * Single Value Decomposition
            */
            Eigen::JacobiSVD <Eigen::MatrixXd > svdOfA (A, Eigen::ComputeThinU | Eigen::ComputeThinV);

            s = svdOfA.singularValues(); //!eigenvalues

            #ifdef INDIRECT_KALMAN_FILTER_DEBUG_PRINTS
            std::cout<<"[SPD-SVD] s: \n"<<s<<"\n";
            std::cout<<"[SPD-SVD] svdOfA.matrixU():\n"<<svdOfA.matrixU()<<"\n";
            std::cout<<"[SPD-SVD] svdOfA.matrixV():\n"<<svdOfA.matrixV()<<"\n";

            Eigen::EigenSolver<_MatrixType> eig(A);
            std::cout << "[SPD-SVD] BEFORE: eigen values: " << eig.eigenvalues().transpose() << std::endl;
            #endif

            for (register int i=0; i<s.size(); ++i)
            {
                #ifdef INDIRECT_KALMAN_FILTER_DEBUG_PRINTS
                std::cout<<"[SPD-SVD] i["<<i<<"]\n";
                #endif

                if (s(i) < 0.00)
                    s(i) = 0.00;
            }

            spdA = svdOfA.matrixU() * s.matrix().asDiagonal() * svdOfA.matrixV().transpose();

            #ifdef INDIRECT_KALMAN_FILTER_DEBUG_PRINTS
            Eigen::EigenSolver<_MatrixType> eigSPD(spdA);
            if (eig.eigenvalues() == eigSPD.eigenvalues())
                std::cout<<"[SPD-SVD] EQUAL!!\n";

            std::cout << "[SPD-SVD] AFTER: eigen values: " << eigSPD.eigenvalues().transpose() << std::endl;
            #endif

            return spdA;
        };
    };
} // end namespace filter

#endif
