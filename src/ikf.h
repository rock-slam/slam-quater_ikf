/**\file ikf.h
 * Header function file and defines
 */
#include <Eigen/Geometry> /**< Eigen data type for Matrix, Quaternion, etc... */


namespace filter
{
  using namespace Eigen;
  
  /** General defines **/
  #ifndef OK
  #define OK	0  /**< Integer value in order to return when everything is all right. */
  #endif
  #ifndef ERROR
  #define ERROR	-1  /**< Integer value in order to return when an error occured. */
  #endif

  /** IKF constant parameters **/
  #ifndef IKFSTATEVECTORSIZE
  #define IKFSTATEVECTORSIZE 9 /**< Number of variables of the vector state-space representation **/
  #endif
  #ifndef QUATERSIZE
  #define QUATERSIZE 4 /**< Number of parameters of a quaternion **/
  #endif

  #ifndef PI
  #define PI 3.141592653589793238462643383279502884197169399375105820974944592307816406286 /**< Pi Number */
  #endif
  #define M1 5 /**< Parameter for adaptive algorithm */
  #define M2 3 /**< Parameter for adaptive algorithm */
  #define GAMMA 0.1 /**< Parameter for adaptive algorithm */
  #define R2COUNT 100 /**< Parameter for adaptive algorithm */
  
  #define D2R PI/180.00 /**< Convert degree to radian **/
  #define R2D 180.00/PI /**< Convert radian to degree **/
  
  /** Sensors constant parameters **/
  #ifndef NUMAXIS
  #define NUMAXIS 3 /**< Number of axis sensed by the sensors **/
  #endif
  
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
  
  /** Commented due to Eigen3 update (with Eigen3 JacobiSVD new class this is not needed ) **/
  // typedef Eigen::Matrix<double,NUMAXIS, NUMAXIS> MatrixMeasurement; /**< Measurement matrix type definition (necessary for SVD decomposition) */
  
  
  class ikf
  {
    /**
     * Filter members
     */
    private:
      int r1count; /**< Variable used in the adaptive algorithm, to compute the Uk matrix for SVD*/
      double r2count; /**< Variable used in the adaptive algorithm, to compute the final Qstart cov. matrix*/
      Eigen::Matrix <double,IKFSTATEVECTORSIZE,1> x; /**< State vector */
      Eigen::Matrix <double,NUMAXIS,1> gtilde; /**< gravitation acceleration */
      Eigen::Matrix <double,NUMAXIS,1> mtilde; /**< Magnetic dip angle */
      Eigen::Quaternion <double> q4;  /**< Attitude quaternion */
      Eigen::Matrix <double,QUATERSIZE,QUATERSIZE> oldomega4; /**< Quaternion integration matrix */
      Eigen::Matrix <double,IKFSTATEVECTORSIZE,IKFSTATEVECTORSIZE> P; /**< Error convariance matrix */
      Eigen::Matrix <double,IKFSTATEVECTORSIZE,IKFSTATEVECTORSIZE> A; /**< System matrix */
      Eigen::Matrix <double,IKFSTATEVECTORSIZE,IKFSTATEVECTORSIZE> Q; /**< Process noise convariance matrix */
      Eigen::Matrix <double,NUMAXIS,NUMAXIS> R; /**< Measurement noise convariance matrix */
      Eigen::Matrix <double,NUMAXIS,NUMAXIS*M1> RHist; /**< History of M1 measurement noise convariance matrix (for the adaptive algorithm) */
      Eigen::Matrix <double,NUMAXIS,NUMAXIS> Ra; /**< Measurement noise convariance matrix for acc */
      Eigen::Matrix <double,NUMAXIS,NUMAXIS> Rg; /**< Measurement noise convariance matrix for gyros */
      Eigen::Matrix <double,NUMAXIS,NUMAXIS> Rm; /**< Measurement noise convariance matrix for mag */
      Eigen::Matrix <double,NUMAXIS,IKFSTATEVECTORSIZE> H1; /**< Measurement 1 Observation matrix */
      Eigen::Matrix <double,NUMAXIS,IKFSTATEVECTORSIZE> H2; /**< Measurement 2 Observation matrix */
      Eigen::Matrix <double,NUMAXIS,1> bghat; /**< Estimated bias for gyroscope */
      Eigen::Matrix <double,NUMAXIS,1> bahat; /**< Estimated bias for accelerometer */

    protected:
    
    public:
      
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
      * @return OK is everything all right. ERROR on other cases.
      *
      */
      int setAttitude (Eigen::Quaternion <double> *initq);
      
      /**
      * @brief This function Initialize the State vector
      * 
      * The state vector is formed by 9 element.
      * (0-2) -> the vector patr of a error quaternion
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
      * @return OK is everything all right. ERROR on other cases.
      *
      */
      int setOmega (Eigen::Matrix <double,NUMAXIS,1>  *u);
      
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
      void Init(Eigen::Matrix <double,IKFSTATEVECTORSIZE,IKFSTATEVECTORSIZE> *P_0, Eigen::Matrix <double,NUMAXIS,NUMAXIS> *Ra, Eigen::Matrix <double,NUMAXIS,NUMAXIS> *Rg, Eigen::Matrix <double,NUMAXIS,NUMAXIS> *Rm,
		   Eigen::Matrix <double,NUMAXIS,NUMAXIS> *Qbg, Eigen::Matrix <double,NUMAXIS,NUMAXIS> *Qba, double g, double alpha);
      
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
      * @param[in]  magn_on_off boolean value to connect or disconnect the magnetometers correction
      * @return void
      *
      */
      void update(Eigen::Matrix <double,NUMAXIS,1>  *acc, Eigen::Matrix <double,NUMAXIS,1>  *mag, bool magn_on_off);
      
      
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

  };

} // end namespace filter
