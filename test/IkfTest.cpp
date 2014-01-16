#define BOOST_TEST_MODULE template_for_test_test
#include <boost/test/included/unit_test.hpp>

#include <quater_ikf/Ikf.hpp> /**< IKF header */
#include <iostream> /**< IO C++ Standard library */

#ifndef D2R
#define D2R M_PI/180.00 /** Convert degree to radian **/
#endif
#ifndef R2D
#define R2D 180.00/M_PI /** Convert radian to degree **/
#endif

#ifndef M1
#define M1 10
#endif

#ifndef M2
#define M2 30
#endif

BOOST_AUTO_TEST_CASE( QUATER_IKF_3_3)
{
    Eigen::Matrix3d Ra; /** Measurement noise covariance matrix for acc */
    Eigen::Matrix3d Rg; /** Measurement noise covariance matrix for gyros */
    Eigen::Matrix3d Rm; /** Measurement noise covariance matrix for mag */
    Eigen::Matrix3d Ri; /** Measurement noise covariance matrix for inclinometers */
    Eigen::Matrix <double, 12, 12> P_0; /** Initial covariance matrix **/
    Eigen::Matrix3d Qbg; /** Noise for the gyros bias instability **/
    Eigen::Matrix3d Qba; /** Noise for the acc bias instability **/
    Eigen::Matrix3d Qbi; /** Noise for the inclinometers bias instability **/

    filter::Ikf<double, true, true> myfilter;

    Ra = Eigen::Matrix3d::Zero();

    Rg = Eigen::Matrix3d::Zero();

    Rm = Eigen::Matrix3d::Zero();

    Ri = Eigen::Matrix3d::Zero();

    /** Noise for error in gyros bias instability **/
    Qbg.setZero();

    /** Noise for error in accelerometers bias instability **/
    Qba.setZero();

    /** Noise for error in inclinometers bias instability **/
    Qbi.setZero();


    /** Initial error covariance **/
    P_0.setZero();
    P_0.block <3, 3> (0,0) = 1.0e-06 * Eigen::Matrix3d::Identity();//Error quaternion
    P_0.block <3, 3> (3,3) = 1.0e-06 * Eigen::Matrix3d::Identity();//Gyros bias
    P_0.block <3, 3> (6,6) = 1.0e-06 * Eigen::Matrix3d::Identity();//Accelerometers bias
    P_0.block <3, 3> (9,9) = 1.0e-06 * Eigen::Matrix3d::Identity();//Inclinometers bias

    /** Theoretical Gravity **/
    double gravity = 9.81;
    double dip_angle = 0.81;

    /** Initialize the filter, including the adaptive part **/
    myfilter.Init(P_0, Ra, Rg, Rm, Ri, Qbg, Qba, Qbi, gravity, dip_angle,
            M1, M2, 0.002,
            M1, M2, 0.04);


    /** Predict **/
    Eigen::Vector3d gyro;
    Eigen::Vector3d acc, mag, incl;
    myfilter.predict(gyro, 0.1);

    /** Update **/
    myfilter.update(acc, true, incl, true, mag, false);
    myfilter.update(acc, true, incl, true);
    myfilter.update(acc, true);
    myfilter.update(acc, incl);

    /** Testing Get **/
    std::cout<<"getGyroBias:\n"<< myfilter.getGyroBias()<<"\n";
    std::cout<<"getAccBias:\n"<< myfilter.getAccBias()<<"\n";
    std::cout<<"getInclBias:\n"<<myfilter.getInclBias()<<"\n";
    std::cout<<"getState:\n"<<myfilter.getState()<<"\n";
    std::cout<<"getGravityinBody:\n"<<myfilter.getGravityinBody()<<"\n";
    std::cout<<"getCovariance:\n"<<myfilter.getCovariance()<<"\n";

    /** Testing Set **/
    myfilter.setAttitude(myfilter.getAttitude());
    myfilter.setState(myfilter.getState());
    myfilter.setInitBias(myfilter.getGyroBias(), myfilter.getAccBias(), myfilter.getInclBias());
    myfilter.setGravity(9.81);
    myfilter.setCovariance(myfilter.getCovariance());
}

BOOST_AUTO_TEST_CASE( QUATER_IKF_3_0)
{
    Eigen::Matrix3d Ra; /** Measurement noise covariance matrix for acc */
    Eigen::Matrix3d Rg; /** Measurement noise covariance matrix for gyros */
    Eigen::Matrix3d Rm; /** Measurement noise covariance matrix for mag */
    Eigen::Matrix3d Ri; /** Measurement noise covariance matrix for inclinometers */
    Eigen::Matrix <double, 9, 9> P_0; /** Initial covariance matrix **/
    Eigen::Matrix3d Qbg; /** Noise for the gyros bias instability **/
    Eigen::Matrix3d Qba; /** Noise for the acc bias instability **/
    Eigen::Matrix3d Qbi; /** Noise for the inclinometers bias instability **/

    filter::Ikf<double, true, false> myfilter;

    Ra = Eigen::Matrix3d::Zero();

    Rg = Eigen::Matrix3d::Zero();

    Rm = Eigen::Matrix3d::Zero();

    Ri = Eigen::Matrix3d::Zero();

    /** Noise for error in gyros bias instability **/
    Qbg.setZero();

    /** Noise for error in accelerometers bias instability **/
    Qba.setZero();

    /** Noise for error in inclinometers bias instability **/
    Qbi.setZero();


    /** Initial error covariance **/
    P_0.setZero();
    P_0.block <3, 3> (0,0) = 1.0e-06 * Eigen::Matrix3d::Identity();//Error quaternion
    P_0.block <3, 3> (3,3) = 1.0e-06 * Eigen::Matrix3d::Identity();//Gyros bias
    P_0.block <3, 3> (6,6) = 1.0e-06 * Eigen::Matrix3d::Identity();//Accelerometers bias

    /** Theoretical Gravity **/
    double gravity = 9.81;
    double dip_angle = 0.81;

    /** Initialize the filter, including the adaptive part **/
    myfilter.Init(P_0, Ra, Rg, Rm, Ri, Qbg, Qba, Qbi, gravity, dip_angle,
            M1, M2, 0.002,
            M1, M2, 0.04);

    /** Predict **/
    Eigen::Vector3d gyro;
    Eigen::Vector3d acc, mag, incl;
    myfilter.predict(gyro, 0.1);

    /** Update **/
    myfilter.update(acc, true, incl, false, mag, false);
    myfilter.update(acc, true);

    /** Testing Get **/
    std::cout<<"getGyroBias:\n"<< myfilter.getGyroBias()<<"\n";
    std::cout<<"getAccBias:\n"<< myfilter.getAccBias()<<"\n";
    std::cout<<"getInclBias:\n"<<myfilter.getInclBias()<<"\n";
    std::cout<<"getState:\n"<<myfilter.getState()<<"\n";
    std::cout<<"getGravityinBody:\n"<<myfilter.getGravityinBody()<<"\n";
    std::cout<<"getCovariance:\n"<<myfilter.getCovariance()<<"\n";

    /** Testing Set **/
    myfilter.setAttitude(myfilter.getAttitude());
    myfilter.setState(myfilter.getState());
    myfilter.setInitBias(myfilter.getGyroBias(), myfilter.getAccBias(), myfilter.getInclBias());
    myfilter.setGravity(9.81);
    myfilter.setCovariance(myfilter.getCovariance());
}

BOOST_AUTO_TEST_CASE( QUATER_IKF_0_3)
{
    Eigen::Matrix3d Ra; /** Measurement noise covariance matrix for acc */
    Eigen::Matrix3d Rg; /** Measurement noise covariance matrix for gyros */
    Eigen::Matrix3d Rm; /** Measurement noise covariance matrix for mag */
    Eigen::Matrix3d Ri; /** Measurement noise covariance matrix for inclinometers */
    Eigen::Matrix <double, 9, 9> P_0; /** Initial covariance matrix **/
    Eigen::Matrix3d Qbg; /** Noise for the gyros bias instability **/
    Eigen::Matrix3d Qba; /** Noise for the acc bias instability **/
    Eigen::Matrix3d Qbi; /** Noise for the inclinometers bias instability **/

    filter::Ikf<double, false, true> myfilter;

    Ra = Eigen::Matrix3d::Zero();

    Rg = Eigen::Matrix3d::Zero();

    Rm = Eigen::Matrix3d::Zero();

    Ri = Eigen::Matrix3d::Zero();

    /** Noise for error in gyros bias instability **/
    Qbg.setZero();

    /** Noise for error in accelerometers bias instability **/
    Qba.setZero();

    /** Noise for error in inclinometers bias instability **/
    Qbi.setZero();


    /** Initial error covariance **/
    P_0.setZero();
    P_0.block <3, 3> (0,0) = 1.0e-06 * Eigen::Matrix3d::Identity();//Error quaternion
    P_0.block <3, 3> (3,3) = 1.0e-06 * Eigen::Matrix3d::Identity();//Gyros bias
    P_0.block <3, 3> (6,6) = 1.0e-06 * Eigen::Matrix3d::Identity();//Accelerometers bias

    /** Theoretical Gravity **/
    double gravity = 9.81;
    double dip_angle = 0.81;

    /** Initialize the filter, including the adaptive part **/
    myfilter.Init(P_0, Ra, Rg, Rm, Ri, Qbg, Qba, Qbi, gravity, dip_angle,
            M1, M2, 0.002,
            M1, M2, 0.04);

    /** Predict **/
    Eigen::Vector3d gyro;
    Eigen::Vector3d acc, mag, incl;
    myfilter.predict(gyro, 0.1);

    /** Update **/
    myfilter.update(acc, false, incl, true, mag, false);
    myfilter.update(acc, false, incl, true);

    /** Testing Get **/
    std::cout<<"getGyroBias:\n"<< myfilter.getGyroBias()<<"\n";
    std::cout<<"getAccBias:\n"<< myfilter.getAccBias()<<"\n";
    std::cout<<"getInclBias:\n"<<myfilter.getInclBias()<<"\n";
    std::cout<<"getState:\n"<<myfilter.getState()<<"\n";
    std::cout<<"getGravityinBody:\n"<<myfilter.getGravityinBody()<<"\n";
    std::cout<<"getCovariance:\n"<<myfilter.getCovariance()<<"\n";

    /** Testing Set **/
    myfilter.setAttitude(myfilter.getAttitude());
    myfilter.setState(myfilter.getState());
    myfilter.setInitBias(myfilter.getGyroBias(), myfilter.getAccBias(), myfilter.getInclBias());
    myfilter.setGravity(9.81);
    myfilter.setCovariance(myfilter.getCovariance());
}

