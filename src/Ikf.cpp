/**\file ikf.cpp
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
 * "Orientation estimation using a quaternion-based indirect Klaman filter with adaptive estimation of external acceleration"
 * A copy if the manuscript can be found in the /doc folder of the library.
 * 
 * @author Javier Hidalgo Carrio | DFKI RIC Bremen | javier.hidalgo_carrio@dfki.de
 * @date June 2011.
 * @version 1.0.
 */

#include <iostream> /**< IO C++ Standard library */
#include <algorithm> /**< Algorithm C++ Standard library */
#include <Eigen/LU> /**< Lineal algebra of Eigen */
#include <Eigen/SVD> /**< Singular Value Decomposition (SVD) of Eigen */
#include <base/Pose.hpp>
#include "Ikf.hpp" /**< Indirect Kalman Filter */

/** WGS-84 ellipsoid constants (Nominal Gravity Model and Earth angular velocity) **/
#ifndef Re
#define Re	6378137 /**< Equatorial radius in meters **/
#endif
#ifndef Rp
#define Rp	6378137 /**< Polar radius in meters **/
#endif
#ifndef ECC
#define ECC  0.0818191908426 /**< First eccentricity **/
#endif
#ifndef GRAVITY
#define GRAVITY 9.79766542 /**< Mean value of gravity value in m/s^2 **/
#endif
#ifndef GWGS0
#define GWGS0 9.7803267714 /**< Gravity value at the equator in m/s^2 **/
#endif
#ifndef GWGS1
#define GWGS1 0.00193185138639 /**< Gravity formula constant **/
#endif
#ifndef EARTHW
#define EARTHW  7.292115e-05 /**< Earth angular velocity in rad/s **/
#endif


//#define DEBUG_PRINTS 1

namespace filter
{

    /** Name spaces to use **/
    using namespace Eigen;
    using namespace std;

    /** Indirect Kalman Filter methods **/

    /**
    * @brief This function Initialize the vectors and matrix of the IKF
    */
    void ikf::Init(const Eigen::Matrix <double,ikf::IKFSTATEVECTORSIZE,ikf::IKFSTATEVECTORSIZE> &P_0,
                    const Eigen::Matrix <double,ikf::NUMAXIS,ikf::NUMAXIS> &Ra,
                    const Eigen::Matrix <double,NUMAXIS,NUMAXIS> &Rg,
                    const Eigen::Matrix <double,NUMAXIS,NUMAXIS> &Rm,
                    const Eigen::Matrix <double,ikf::NUMAXIS,ikf::NUMAXIS> &Qbg,
                    const Eigen::Matrix <double,ikf::NUMAXIS,ikf::NUMAXIS> &Qba,
                    double g, double alpha, unsigned int m1, unsigned int m2, double gamma)
    {

      /** Gravitation acceleration **/
      gtilde << 0, 0, g;

      /** Dip angle (alpha is in rad) **/
      mtilde(0) = cos(alpha);
      mtilde(1) = 0;
      mtilde(2) = -sin(alpha);


      /** Kalman filter state, error covariance and process noise covariance **/
      x = Matrix <double,ikf::IKFSTATEVECTORSIZE,1>::Zero();

      Q = Matrix <double,ikf::IKFSTATEVECTORSIZE,ikf::IKFSTATEVECTORSIZE>::Zero();
      Q.block <NUMAXIS, NUMAXIS> (0,0) = 0.25 * (Rg);
      Q.block <NUMAXIS, NUMAXIS> (3,3) = (Qbg);
      Q.block <NUMAXIS, NUMAXIS> (6,6) = (Qba);

      /** Initial error covariance **/
      P = (P_0);

      H1 = Matrix <double,NUMAXIS,ikf::IKFSTATEVECTORSIZE>::Zero();
      H2 = Matrix <double,NUMAXIS,ikf::IKFSTATEVECTORSIZE>::Zero();
      H1(0,6) = 1; H1(1,7) = 1; H1(2,8) = 1;

      /** System matrix A **/
      A = Matrix <double,ikf::IKFSTATEVECTORSIZE,ikf::IKFSTATEVECTORSIZE>::Zero();
      A(0,3) = -0.5;A(1,4) = -0.5;A(2,5) = -0.5;

      /** Initial bias **/
      bghat = Matrix <double,NUMAXIS,1>::Zero();
      bahat = Matrix <double,NUMAXIS,1>::Zero();

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

      /** Default initial bias **/
      bghat << 0.00, 0.00, 0.00;
      bahat << 0.00, 0.00, 0.00;

      /** Fill matrix Rq, Ra and Rm **/
      this->ikf::Ra = (Ra);
      this->ikf::Rg = (Rg);
      this->ikf::Rm = (Rm);

      /** Initialize adaptive object **/
      initAdaptiveAttitude(m1, m2, gamma);

      /** Print filter information **/
      #ifdef DEBUG_PRINTS
      std::cout<< "P:\n"<<P<<"\n";
      std::cout<< "Q:\n"<<Q<<"\n";
      std::cout<< "H1:\n"<<H1<<"\n";
      std::cout<< "H2:\n"<<H2<<"\n";
      std::cout<< "A:\n"<<A<<"\n";
      std::cout<< "mtilde:\n"<<mtilde<<"\n";
      std::cout<< "gtilde:\n"<<gtilde<<"\n";
      std::cout<< "Ra:\n"<<(*Ra)<<"\n";
      std::cout<< "Rg:\n"<<(*Rg)<<"\n";
      std::cout<< "Rm:\n"<<(*Rm)<<"\n";
      #endif

      return;
    }

    void ikf::initAdaptiveAttitude(const unsigned int M1, const unsigned int M2,
                        const double gamma)
    {
        this->adapAtt.reset(new filter::AdaptiveAttitudeCov (M1, M2, gamma));
    }

    void ikf::setInitBias(const Eigen::Matrix<double, NUMAXIS, 1> &gbias,
            const Eigen::Matrix<double, NUMAXIS, 1> &abias)
    {
        this->bghat = gbias;
        this->bahat = abias;

        return;
    }

    /**
    * @brief Set the current state vector of the filter
    */
    void ikf::setState(Eigen::Matrix< double, ikf::IKFSTATEVECTORSIZE , 1  > *x_0)
    {
      x = (*x_0);

      return;
    }

    /**
    * @brief This function Initilize Attitude
    */
    bool ikf::setAttitude(Eigen::Quaternion< double > *initq)
    {
      if (initq != NULL)
      {
	/** Initial orientation **/
	q4 = (*initq);
	
	return true;
      }

      return false;
    }

    /**
    * @brief This function set the initial Omega matrix
    */
    bool ikf::setOmega(Eigen::Matrix< double, ikf::NUMAXIS , 1  >* u)
    {
      if (u != NULL)
      {
	/** Initialization for quaternion integration **/
	oldomega4 << 0,-(*u)(0), -(*u)(1), -(*u)(2),
	  (*u)(0), 0, (*u)(2), -(*u)(1),
	  (*u)(1), -(*u)(2), 0, (*u)(0),
	  (*u)(2), (*u)(1), -(*u)(0), 0;

	return true;
      }
      return false;
    }

    /**
    * @brief Gravity
    */
    void ikf::setGravity(double gravity)
    {
        gtilde << 0.00, 0.00, gravity;
        return;
    }

    /**
    * @brief Set Noise covariance matrix
    */
    void ikf::setCovariance(const Eigen::Matrix< double, ikf::IKFSTATEVECTORSIZE , ikf::IKFSTATEVECTORSIZE> &Pk)
    {
        P = Pk;
        return;
    }

    /**
    * @brief Gets the current orientation in Euler angles
    */
    Eigen::Matrix< double, ikf::NUMAXIS , 1  > ikf::getEuler()
    {
      Eigen::Matrix <double, NUMAXIS, 1> euler;

      //std::cout << Eigen::Matrix3d(q4) << std::endl; 
      Vector3d e = base::getEuler(q4);
       euler(0) = e[2]; 
       euler(1) = e[1]; 
       euler(2) = e[0]; 
//       std::cout << "Attitude (getEuler): "<< euler(0)<<" "<<euler(1)<<" "<<euler(2)<<"\n";
       //std::cout << "Attitude in degrees (getEuler): "<< euler(0)*R2D<<" "<<euler(1)*R2D<<" "<<euler(2)*R2D<<"\n";

      return euler;
    }

    /**
    * @brief Gets the current orientation in Quaternion
    */
    Eigen::Quaternion< double > ikf::getAttitude()
    {
      return q4;
    }

    /**
    * @brief Gets the current gyroscopes bias
    */
    Eigen::Matrix<double, ikf::NUMAXIS, 1> ikf::getGyroBias()
    {
        return this->bghat;
    }

    /**
    * @brief Gets the current accelerometers bias
    */
    Eigen::Matrix<double, ikf::NUMAXIS, 1> ikf::getAccBias()
    {
        return this->bahat;
    }


    /**
    * @brief Gets the current state vector of the filter
    */
    Eigen::Matrix< double, ikf::IKFSTATEVECTORSIZE , 1  > ikf::getState()
    {
      return x;

    }

    /**
    * @brief Gets gravity in IMU body frame
    */
    Eigen::Matrix<double, ikf::NUMAXIS, 1> ikf::getGravityinBody()
    {
        return q4.inverse() * gtilde;
    }

    /**
    * @brief Gets Noise covariance matrix
    */
    Eigen::Matrix< double, ikf::IKFSTATEVECTORSIZE , ikf::IKFSTATEVECTORSIZE> ikf::getCovariance()
    {
	return P;
    }


    /**
    * @brief Performs the prediction step of the filter.
    */
    void ikf::predict(Eigen::Matrix< double, ikf::NUMAXIS , 1  >* u, double dt)
    {
      Eigen::Matrix <double,NUMAXIS,NUMAXIS> vec2product; /**< Vec 2 product  matrix */
      Eigen::Matrix <double,NUMAXIS,1> angvelo; /**< Vec 2 product  matrix */
      Eigen::Matrix <double,QUATERSIZE,QUATERSIZE> omega4; /**< Quaternion integration matrix */
      Eigen::Matrix <double,QUATERSIZE,1> quat; /**< Quaternion integration matrix */
      Eigen::Matrix <double,IKFSTATEVECTORSIZE,IKFSTATEVECTORSIZE> dA; /**< Discrete System matrix */
      Eigen::Matrix <double,IKFSTATEVECTORSIZE,IKFSTATEVECTORSIZE> Qd; /**< Discrete Q matrix */

      /** Compute the vector2product matrix with the angular velocity **/
      angvelo = (*u) - bghat; /** Eliminate the Bias **/

      vec2product << 0, -angvelo(2), angvelo(1),
		    angvelo(2), 0, -angvelo(0),
		    -angvelo(1), angvelo(0), 0;

      /** Compute the dA Matrix **/
      A.block<NUMAXIS, NUMAXIS> (0,0) = -vec2product;
      dA = Matrix<double,IKFSTATEVECTORSIZE,IKFSTATEVECTORSIZE>::Identity() + A * dt + 0.5 *A * A * pow(dt,2);

      /** Propagate the vector through the system **/
      x = dA * x;
      Qd = Q*dt + 0.5*dt*A*Q + 0.5*dt*Q*A.transpose();
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
      quat = (Matrix<double,QUATERSIZE,QUATERSIZE>::Identity() +(0.75 * omega4 *dt)-(0.25 * oldomega4 * dt) -
      ((1.0/6.0) * angvelo.squaredNorm() * pow(dt,2) *  Matrix<double,QUATERSIZE,QUATERSIZE>::Identity()) -
      ((1.0/24.0) * omega4 * oldomega4 * pow(dt,2)) - ((1.0/48.0) * angvelo.squaredNorm() * omega4 * pow(dt,3))) * quat;

      q4.w() = quat(0);
      q4.x() = quat(1);
      q4.y() = quat(2);
      q4.z() = quat(3);
      q4.normalize();

      oldomega4 = omega4;

      return;

    }


     /**
    * @brief Performs the measurement and correction steps of the filter.
    */
    void ikf::update(Eigen::Matrix< double, ikf::NUMAXIS , 1  >* acc, Eigen::Matrix< double, ikf::NUMAXIS , 1  >* mag, bool magn_on)
    {
      Eigen::Matrix <double,NUMAXIS,NUMAXIS> vec2product; /**< Vec 2 product  matrix */
      Eigen::Matrix <double,NUMAXIS,NUMAXIS> fooR2; /**<  Measurement noise matrix from accelerometers matrix Ra*/
      Eigen::Matrix <double,IKFSTATEVECTORSIZE,IKFSTATEVECTORSIZE> P1; /**< Error convariance matrix for measurement 1*/
      Eigen::Matrix <double,IKFSTATEVECTORSIZE,IKFSTATEVECTORSIZE> P2; /**< Error convariance matrix for measurement 2*/
      Eigen::Matrix <double,IKFSTATEVECTORSIZE,IKFSTATEVECTORSIZE> auxM; /**< Auxiliar matrix for computing Kalman gain in measurement*/
      Eigen::Matrix <double,IKFSTATEVECTORSIZE, NUMAXIS> K1; /**< Kalman Gain matrix for measurement 1*/
      Eigen::Matrix <double,IKFSTATEVECTORSIZE, NUMAXIS> K2; /**< Kalman Gain matrix for measurement 2*/
      Eigen::Matrix <double,NUMAXIS,NUMAXIS> R1; /**< Acceleration covariance matrix */
      Eigen::Quaternion <double> qe;  /**< Attitude error quaternion */
      Eigen::Matrix <double,NUMAXIS,1> gtilde_body; /**< Gravitation in the body frame */
      Eigen::Matrix <double,NUMAXIS,1> mtilde_body; /**< Magnetic field in the body frame */
      Eigen::Matrix <double,NUMAXIS,1> z1; /**< Measurement vector 1 Acc */
      Eigen::Matrix <double,NUMAXIS,1> z2; /**< Measurement vector 2 Mag */
      Eigen::Matrix <double,NUMAXIS,1> auxvector; /**< Auxiliar vector variable */


      /**----------------------- **/
      /** Measurement step 1 Acc **/
      /**----------------------- **/

      /** First measurement step (Pitch and Roll correction from Acc) **/
      gtilde_body = q4.inverse() * gtilde;
      vec2product << 0, -gtilde_body(2), gtilde_body(1),
		    gtilde_body(2), 0, -gtilde_body(0),
		    -gtilde_body(1), gtilde_body(0), 0;

      H1.block<NUMAXIS, NUMAXIS> (0,0) = 2*vec2product;

      /** Measurement **/
      z1 = (*acc) - bahat - gtilde_body;

      #ifdef DEBUG_PRINTS
      std::cout<<"acc:\n"<<*acc<<"\n";
      std::cout<<"z1:\n"<<z1<<"\n";
      std::cout<<"g_body:\n"<<gtilde_body<<"\n";
      #endif

      /** The adaptive algorithm **/
      R1 = adapAtt->matrix<IKFSTATEVECTORSIZE> (x, P, z1, H1, Ra);

      /** Compute the Kalman Gain Matrix **/
      P1 = P;
      Eigen::Matrix<double, NUMAXIS, NUMAXIS> S1, S1_inverse;
      S1 = H1 * P1 * H1.transpose() + R1;
      S1_inverse = S1.inverse();
      K1 = P1 * H1.transpose() * S1_inverse;
      Eigen::Matrix<double, NUMAXIS, 1> innovation = (z1 - H1 * x);

      /** Update the state vector and the covariance matrix **/
      x = x + K1 * innovation;
      P = (Matrix<double,IKFSTATEVECTORSIZE,IKFSTATEVECTORSIZE>::Identity()
              -K1*H1)*P*(Matrix<double,IKFSTATEVECTORSIZE,IKFSTATEVECTORSIZE>::Identity()
              -K1*H1).transpose() + K1*R1*K1.transpose();
      P = 0.5 * (P + P.transpose());//Guarantee symmetry

      #ifdef DEBUG_PRINTS
      std::cout<<"x(k+1|k+1):\n"<<x<<"\n";
      std::cout<<"P(k+1|k+1):\n"<<P<<"\n";
      std::cout<<"innovation:\n"<<innovation<<"\n";
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
      x.block<NUMAXIS,1>(0,0) = Matrix<double, NUMAXIS, 1>::Zero();


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
	H2.block<NUMAXIS, NUMAXIS> (0,0) = 2*vec2product;

	/** Measurement vector **/
	z2 = (*mag) - mtilde_body;

	P2 = Matrix<double, IKFSTATEVECTORSIZE, IKFSTATEVECTORSIZE>::Zero();
	P2.block<NUMAXIS, NUMAXIS>(0,0) = P.block<NUMAXIS, NUMAXIS>(0,0);

	auxvector << 0, 0, 1;
	auxvector = q4.inverse() * auxvector;

	/** Compute Kalman Gain **/
	auxM = Matrix<double, IKFSTATEVECTORSIZE, IKFSTATEVECTORSIZE>::Zero();
	auxM.block<NUMAXIS, NUMAXIS>(0,0) = auxvector * auxvector.transpose();
	K2 = auxM * P2 * H2.transpose() * (H2*P2*H2.transpose() + Rm).inverse();

	/** Update the state vector and the covariance matrix **/
	x = x + K2*(z2 - (H2*x));
	P = P - K2 * H2 * P - P * H2.transpose() * K2.transpose() + K2*(H2*P*H2.transpose() + Rm)*K2.transpose();
	P = 0.5 * (P + P.transpose());

	/** Update the quaternion with the Indirect approach **/
	qe.w() = 1;
	qe.x() = x(0);
	qe.y() = x(1);
	qe.z() = x(2);
	q4 = q4 * qe;

	/** Normalize quaternion **/
	q4.normalize();

	/** Reset the quaternion part of the state vector **/
	x.block<NUMAXIS,1>(0,0) = Matrix<double, NUMAXIS, 1>::Zero();
      }

      /**---------------------------- **/
      /** Reset the rest of the state **/
      /**---------------------------- **/
      bghat = bghat + x.block<NUMAXIS, 1> (3,0);
      x.block<NUMAXIS, 1> (3,0) = Matrix <double, NUMAXIS, 1>::Zero();

      bahat = bahat + x.block<NUMAXIS, 1> (6,0);
      x.block<NUMAXIS, 1> (6,0) = Matrix <double, NUMAXIS, 1>::Zero();

      #ifdef DEBUG_PRINTS
      std::cout<<"bahat:\n"<<bahat<<"\n";
      std::cout<<"bghat:\n"<<bghat<<"\n";
      #endif


      return;
    }


    /**
    * @brief Conversion Quaternion to DCM (Direct Cosine Matrix) (Alternative to Eigen)
    */
    void ikf::Quaternion2DCM(Eigen::Quaternion< double >* q, Eigen::Matrix< double, ikf::NUMAXIS, ikf::NUMAXIS  >*C)
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

}
