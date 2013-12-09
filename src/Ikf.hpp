/**\file ikf.hpp
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
 *
 * @author Javier Hidalgo Carrio | DFKI RIC Bremen | javier.hidalgo_carrio@dfki.de
 * @date June 2011.
 * @version 1.0.
 */

#include "AdaptiveAttitudeCov.hpp"
#include <Eigen/Geometry> /** Eigen data type for Matrix, Quaternion, etc... */

/** Boost **/
#include <boost/shared_ptr.hpp> /** For shared pointers **/

namespace filter
{
  using namespace Eigen;

  class ikf
  {

    /** IKF constant parameters **/
    public:
      enum CONSTS {
        IKFSTATEVECTORSIZE = 9,
        QUATERSIZE = 4,
        NUMAXIS = 3,
      };

    /**
     * Filter members
     */
    private:
      Eigen::Matrix <double,IKFSTATEVECTORSIZE,1> x; /**< State vector */
      Eigen::Matrix <double,NUMAXIS,1> gtilde; /**< gravitation acceleration */
      Eigen::Matrix <double,NUMAXIS,1> mtilde; /**< Magnetic dip angle */
      Eigen::Quaternion <double> q4;  /**< Attitude quaternion */
      Eigen::Matrix <double,QUATERSIZE,QUATERSIZE> oldomega4; /**< Quaternion integration matrix */
      Eigen::Matrix <double,IKFSTATEVECTORSIZE,IKFSTATEVECTORSIZE> P; /**< Error convariance matrix */
      Eigen::Matrix <double,IKFSTATEVECTORSIZE,IKFSTATEVECTORSIZE> A; /**< System matrix */
      Eigen::Matrix <double,IKFSTATEVECTORSIZE,IKFSTATEVECTORSIZE> Q; /**< Process noise convariance matrix */
      Eigen::Matrix <double,NUMAXIS,NUMAXIS> Ra; /**< Measurement noise convariance matrix for acc */
      Eigen::Matrix <double,NUMAXIS,NUMAXIS> Rg; /**< Measurement noise convariance matrix for gyros */
      Eigen::Matrix <double,NUMAXIS,NUMAXIS> Rm; /**< Measurement noise convariance matrix for mag */
      Eigen::Matrix <double,NUMAXIS,IKFSTATEVECTORSIZE> H1; /**< Measurement 1 Observation matrix */
      Eigen::Matrix <double,NUMAXIS,IKFSTATEVECTORSIZE> H2; /**< Measurement 2 Observation matrix */
      Eigen::Matrix <double,NUMAXIS,1> bghat; /**< Estimated bias for gyroscope */
      Eigen::Matrix <double,NUMAXIS,1> bahat; /**< Estimated bias for accelerometer */

      /** Object of Class for Adaptive Measurement of Attitude Covariance Matrix **/
      boost::shared_ptr<filter::AdaptiveAttitudeCov> adapAtt;


    protected:
      void initAdaptiveAttitude(const unsigned int M1, const unsigned int M2, const double gamma);

    public:

    /**
    * @brief Gets the current gyroscopes bias
    */
    Eigen::Matrix<double, ikf::NUMAXIS, 1> getGyroBias();

    /**
    * @brief Gets the current accelerometers bias
    */
    Eigen::Matrix<double, ikf::NUMAXIS, 1> getAccBias();

     /**
      * @brief Gets the current state vector of the filter
      * 
      * @author Javier Hidalgo Carrio.
      *
      * @return State Vector
      *
      */
      Eigen::Matrix <double,IKFSTATEVECTORSIZE,1> getState();

      /**
      * @brief Gets gravity in IMU body frame
      */
      Eigen::Matrix<double, ikf::NUMAXIS, 1> getGravityinBody();


       /**
      * @brief Gets the current orientation in Quaternion
      * 
      * @author Javier Hidalgo Carrio.
      *
      * @return Quaternion with the current orientation.
      *
      */
      Eigen::Quaternion <double> getAttitude();

      /**
      * @brief Gets the current orientation in Euler angles
      * 
      * @author Javier Hidalgo Carrio.
      *
      * @return Current orientation in Euler angles.
      *
      */
      Eigen::Matrix <double, NUMAXIS, 1> getEuler();

      /**
      * @brief Gets Noise covariance matrix
      * 
      * @author Javier Hidalgo Carrio.
      *
      * @return Matrix P of the covariance of the state vector
      *
      */
      Eigen::Matrix <double,IKFSTATEVECTORSIZE,IKFSTATEVECTORSIZE> getCovariance();

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
      bool setAttitude (Eigen::Quaternion <double> *initq);

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
      void setState (Eigen::Matrix <double,IKFSTATEVECTORSIZE,1> *x_0);

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
      bool setOmega (Eigen::Matrix <double,NUMAXIS,1>  *u);

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
      void setInitBias (const Eigen::Matrix<double, NUMAXIS, 1> &gbias,
            const Eigen::Matrix<double, NUMAXIS, 1> &abias);

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

      void setGravity(double gravity);

      /**
      * @brief Set the filter covariance
      *
      * @param[in] Pk covariance matrix of appropriate dimension
      *
      * @return void.
      *
      */

      void setCovariance(const Eigen::Matrix< double, ikf::IKFSTATEVECTORSIZE , ikf::IKFSTATEVECTORSIZE> &Pk);

      /**
      * @brief This function Initilize the vectors and matrix of the IKF
      * 
      * This method receives the measurement noise matrix of the sensors
      * The theoretical gravity value and the Dip angle of the location.
      *
      * @author Javier Hidalgo Carrio.
      *
      * @param[in] P_0 Initial state covariance matrix
      * @param[in] Ra measurement noise matrix of Accelerometers.
      * @param[in] Rg measurement noise matrix of Gyroscopes.
      * @param[in] Rm measurement noise matrix of Magnetometers.
      * @param[in] Qbg covariance noise matrix of the gyroscopes bias
      * @param[in] Qba covariance noise matrix of the accelerometers bias 
      * @param[in] g local gravitational value.
      * @param[in] alpha Dip angle
      *
      * @return void
      *
      */
      void Init(const Eigen::Matrix <double,ikf::IKFSTATEVECTORSIZE,ikf::IKFSTATEVECTORSIZE> &P_0,
                const Eigen::Matrix <double,ikf::NUMAXIS,ikf::NUMAXIS> &Ra,
                const Eigen::Matrix <double,NUMAXIS,NUMAXIS> &Rg,
                const Eigen::Matrix <double,NUMAXIS,NUMAXIS> &Rm,
                const Eigen::Matrix <double,ikf::NUMAXIS,ikf::NUMAXIS> &Qbg,
                const Eigen::Matrix <double,ikf::NUMAXIS,ikf::NUMAXIS> &Qba,
                double g, double alpha, unsigned int m1, unsigned int m2, double gamma);

      /**
      * @brief Performs the prediction step of the filter.
      * 
      * It computes the discrete version of the matrix A to propagate forward
      * the state vector x. It computes the Q and P matrix as well as the 
      * quaternion integration from the input vector u and the delta time.
      *
      * @author Javier Hidalgo Carrio.
      *
      * @param[in] *u pointer to vector with the angular velocity
      * @param[in] dt delta time between samples
      *
      * @return void
      *
      */
      void predict(Eigen::Matrix <double,NUMAXIS,1>  *u, double dt);

      /**
      * @brief Performs the measurement and correction steps of the filter.
      * 
      * The IKf is based on two measurement step:\
      *  1. Measurement step to correct Pitch and Roll from accelerometers.
      *  2. Measurement step to correct Yaw angle from magnetometers.
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
      * @param[in] *acc pointer to vector with accelerations
      * @param[in] *magn pointer to vector with magnetometers
      * @param[in]  magn_on boolean value to connect or disconnect the magnetometers correction
      * @return void
      *
      */
      void update(Eigen::Matrix <double,NUMAXIS,1>  *acc, Eigen::Matrix <double,NUMAXIS,1>  *mag, bool magn_on);

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
      void Quaternion2DCM(Eigen::Quaternion< double >* q, Eigen::Matrix< double, NUMAXIS, NUMAXIS  >*C);

      /**
      * @brief Conversion Quaternion to DCM (Direct Cosine Matrix) (Alternative to Eigen)
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

            #ifdef DEBUG_PRINTS
            std::cout<<"[SPD-SVD] s: \n"<<s<<"\n";
            std::cout<<"[SPD-SVD] svdOfA.matrixU():\n"<<svdOfA.matrixU()<<"\n";
            std::cout<<"[SPD-SVD] svdOfA.matrixV():\n"<<svdOfA.matrixV()<<"\n";

            Eigen::EigenSolver<_MatrixType> eig(A);
            std::cout << "[SPD-SVD] BEFORE: eigen values: " << eig.eigenvalues().transpose() << std::endl;
            #endif

            for (register int i=0; i<s.size(); ++i)
            {
                #ifdef DEBUG_PRINTS
                std::cout<<"[SPD-SVD] i["<<i<<"]\n";
                #endif

                if (s(i) < 0.00)
                    s(i) = 0.00;
            }

            spdA = svdOfA.matrixU() * s.matrix().asDiagonal() * svdOfA.matrixV();

            #ifdef DEBUG_PRINTS
            Eigen::EigenSolver<_MatrixType> eigSPD(spdA);
            if (eig.eigenvalues() == eigSPD.eigenvalues())
                std::cout<<"[SPD-SVD] EQUAL!!\n";

            std::cout << "[SPD-SVD] AFTER: eigen values: " << eigSPD.eigenvalues().transpose() << std::endl;
            #endif

            return spdA;
        };
  };
} // end namespace filter
