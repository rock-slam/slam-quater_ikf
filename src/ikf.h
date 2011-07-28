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
  #define STATEVECTORSIZE 9 /**< Number of variables of the vector state-space representation **/
  #define QUATERSIZE 4 /**< Number of parameters of a quaternion **/

  #ifndef PI
  #define PI 3.141592653589793238462643383279502884197169399375105820974944592307816406286 /**< Pi Number */
  #endif
  #define M1 2 /**< Parameter for adaptive algorithm */
  #define M2 3 /**< Parameter for adaptive algorithm */
  #define GAMMA 0.1 /**< Parameter for adaptive algorithm */
  #define R2COUNT 100 /**< Parameter for adaptive algorithm */
  #define EAST 1 /**< EAST is 1 and means positive magnetic declination **/
  #define WEST 2 /**< WEST is 2 and means negative magnetic declination **/
  
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
      double r2count; /**< Variable used in the adaptive algorithm */
      Eigen::Matrix <double,STATEVECTORSIZE,1> x; /**< State vector */
      Eigen::Matrix <double,NUMAXIS,1> gtilde; /**< gravitation acceleration */
      Eigen::Matrix <double,NUMAXIS,1> mtilde; /**< Magnetic dip angle */
      Eigen::Quaternion <double> q4;  /**< Attitude quaternion */
      Eigen::Matrix <double,QUATERSIZE,QUATERSIZE> oldomega4; /**< Quaternion integration matrix */
      Eigen::Matrix <double,STATEVECTORSIZE,STATEVECTORSIZE> P; /**< Error convariance matrix */
      Eigen::Matrix <double,STATEVECTORSIZE,STATEVECTORSIZE> A; /**< System matrix */
      Eigen::Matrix <double,STATEVECTORSIZE,STATEVECTORSIZE> Q; /**< Process noise convariance matrix */
      Eigen::Matrix <double,NUMAXIS,NUMAXIS> R; /**< Measurement noise convariance matrix */
      Eigen::Matrix <double,NUMAXIS,NUMAXIS> Ra; /**< Measurement noise convariance matrix for acc */
      Eigen::Matrix <double,NUMAXIS,NUMAXIS> Rg; /**< Measurement noise convariance matrix for gyros */
      Eigen::Matrix <double,NUMAXIS,NUMAXIS> Rm; /**< Measurement noise convariance matrix for mag */
      Eigen::Matrix <double,NUMAXIS,STATEVECTORSIZE> H1; /**< Measurement 1 Observation matrix */
      Eigen::Matrix <double,NUMAXIS,STATEVECTORSIZE> H2; /**< Measurement 2 Observation matrix */
      Eigen::Matrix <double,NUMAXIS,1> bghat; /**< Estimated bias for gyroscope */
      Eigen::Matrix <double,NUMAXIS,1> bahat; /**< Estimated bias for accelerometer */

    protected:
    
    public:
      Eigen::Matrix <double,STATEVECTORSIZE,1> getState();
      
      Eigen::Quaternion <double> getAttitude();
      
      Eigen::Matrix <double, NUMAXIS, 1> getEuler();
      
      Eigen::Matrix <double,STATEVECTORSIZE,STATEVECTORSIZE> getCovariance();
      
      int setAttitude (Eigen::Quaternion <double> *initq);
      
      int setOmega (Eigen::Matrix <double,NUMAXIS,1>  *u);
      
      void Init (Eigen::Matrix <double,NUMAXIS,NUMAXIS> *Ra, Eigen::Matrix <double,NUMAXIS,NUMAXIS> *Rg, Eigen::Matrix <double,NUMAXIS,NUMAXIS> *Rm, double g, double alpha);
      
      void predict(Eigen::Matrix <double,NUMAXIS,1>  *u, double dt);
      
      void update(Eigen::Matrix <double,NUMAXIS,1>  *acc, Eigen::Matrix <double,NUMAXIS,1>  *mag);
      
      double GravityModel (double latitude, double altitude);
      
      void SubstractEarthRotation(Eigen::Matrix <double, NUMAXIS, 1> *u, Eigen::Quaternion <double> *qb_g, double latitude);
      
      int CorrectMagneticDeclination (Eigen::Quaternion <double> *quat, double magnetic_declination,  int mode);
      
      void Quaternion2DCM(Eigen::Quaternion< double >* q, Eigen::Matrix< double, NUMAXIS, NUMAXIS  >*C);
      
      void Quaternion2Euler(Eigen::Quaternion< double >* q, Eigen::Matrix< double, NUMAXIS , 1  >* euler);
      
      void Euler2Quaternion(Eigen::Matrix< double, NUMAXIS , 1  >* euler, Eigen::Quaternion< double >* q);
  };

} // end namespace filter
