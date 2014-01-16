#ifndef ADAPTIVE_ATTITUDE_COV_HPP
#define ADAPTIVE_ATTITUDE_COV_HPP

#include <Eigen/Geometry> /** Quaternion, Angle-Axis,...*/
#include <Eigen/StdVector> /** STL container with Eigen types */
#include <Eigen/LU> /** Lineal algebra of Eigen */
#include <Eigen/SVD> /** Singular Value Decomposition (SVD) of Eigen */

#include <vector>

//#define ADAPTIVE_ATTITUDE_COV_DEBUG_PRINTS 1

namespace filter
{
    /**@brief Class for Adaptive measurement matrix for the attitude correction in 3D
    */
    template <typename _Scalar, size_t _DoFState, size_t _MeasurementSize>
    class AdaptiveAttitudeCov
    {

    protected:

        typedef Eigen::Matrix<_Scalar, _MeasurementSize, _MeasurementSize> MeasurementMatrix;

        unsigned int r1count; /** Variable used in the adaptive algorithm, to compute the Uk matrix for SVD*/
        unsigned int m1; /** Parameter for adaptive algorithm (to estimate Uk which is not directly observale) */
        unsigned int m2; /** Parameter for adaptive algorithm (to prevent falser entering in no-external acc mode) */
        double gamma; /** Parameter for adaptive algorithm (only entering when Qstart is greater than RHR'+Ra) */
        unsigned int r2count; /** Parameter for adaptive algorithm */

        /** History of M1 measurement noise covariance matrix (for the adaptive algorithm) */
        std::vector <  MeasurementMatrix, Eigen::aligned_allocator < MeasurementMatrix > > RHist;

    public:

        AdaptiveAttitudeCov(const unsigned int M1, const unsigned int M2,
                        const double GAMMA)
            :m1(M1), m2(M2), gamma(GAMMA)
        {
            r1count = 0;
            r2count = M2;
            RHist.resize(M1);
            for (typename std::vector< MeasurementMatrix, Eigen::aligned_allocator < MeasurementMatrix > >::iterator it = RHist.begin()
                    ; it != RHist.end(); ++it)
                (*it).setZero();

            #ifdef ADAPTIVE_ATTITUDE_COV_DEBUG_PRINTS
            std::cout<<"[INIT ADAPTIVE_ATTITUDE] M1: "<<m1<<"\n";
            std::cout<<"[INIT ADAPTIVE_ATTITUDE] M2: "<<m2<<"\n";
            std::cout<<"[INIT ADAPTIVE_ATTITUDE] GAMMA: "<<gamma<<"\n";
            std::cout<<"[INIT ADAPTIVE_ATTITUDE] r1count: "<<r1count<<"\n";
            std::cout<<"[INIT ADAPTIVE_ATTITUDE] r2count: "<<r2count<<"\n";
            std::cout<<"[INIT ADAPTIVE_ATTITUDE] RHist is of size "<<RHist.size()<<"\n";
            #endif

        }

        ~AdaptiveAttitudeCov(){}

        Eigen::Matrix<_Scalar, _MeasurementSize, _MeasurementSize> matrix(const Eigen::Matrix <_Scalar, _DoFState, 1> xk,
                            const Eigen::Matrix<_Scalar, _DoFState, _DoFState> &Pk,
                            const Eigen::Matrix<_Scalar, _MeasurementSize, 1> &z,
                            const Eigen::Matrix<_Scalar, _MeasurementSize, _DoFState> &H,
                            const Eigen::Matrix<_Scalar, _MeasurementSize, _MeasurementSize>  &R)
        {
            MeasurementMatrix R1a; /** Measurement noise covariance matrix for the adaptive algorithm */
            MeasurementMatrix fooR; /**  Measurement noise matrix from accelerometers matrix Ra */
            MeasurementMatrix Uk; /** Uk measurement noise covariance matrix for the adaptive algorithm */
            MeasurementMatrix Qstar; /** External acceleration covariance matrix */
            MeasurementMatrix u; /** Unitary matrix U for the SVD decomposition */
            Eigen::Matrix<_Scalar, _MeasurementSize, 1> lambda; /** Lambda vector for the adaptive algorithm */
            Eigen::Matrix<_Scalar, _MeasurementSize, 1> mu; /** mu vector for the adaptive algorithm */
            Eigen::Matrix<_Scalar, _MeasurementSize, 1> s; /** Unitary matrix V for the SVD decomposition */

            /** Estimation of R **/
            R1a = (z - H*xk) * (z - H*xk).transpose();

            RHist[r1count] = R1a;

            #ifdef ADAPTIVE_ATTITUDE_COV_DEBUG_PRINTS
            std::cout<<"[ADAPTIVE_ATTITUDE] xk:\n"<<xk<<"\n";
            std::cout<<"[ADAPTIVE_ATTITUDE] Pk:\n"<<Pk<<"\n";
            std::cout<<"[ADAPTIVE_ATTITUDE] z:\n"<<z<<"\n";
            std::cout<<"[ADAPTIVE_ATTITUDE] H:\n"<<H<<"\n";
            std::cout<<"[ADAPTIVE_ATTITUDE] R:\n"<<R<<"\n";
            std::cout<<"[ADAPTIVE_ATTITUDE] r1count:\n"<<r1count<<"\n";
            std::cout<<"[ADAPTIVE_ATTITUDE] R1a:\n"<<R1a<<"\n";
            std::cout<<"[ADAPTIVE_ATTITUDE] z:\n"<<z<<"\n";
            #endif


            /** r1count + 1 modulus the number of history M1 **/
            r1count = (r1count+1)%(m1);

            Uk.setZero();

            /** Starting the Uk is R **/
            for (register int j=0; j<static_cast<int>(m1); ++j)
            {
                Uk += RHist[j];
            }

            Uk = Uk/static_cast<_Scalar>(m1);

            fooR = H*Pk*H.transpose() + R;

            /**
            * Single Value Decomposition
            */
            Eigen::JacobiSVD <Eigen::MatrixXd > svdOfUk(Uk, Eigen::ComputeThinU);

            s = svdOfUk.singularValues(); //!eigenvalues
            u = svdOfUk.matrixU();//!eigenvectors

            for (size_t i=0; i<_MeasurementSize; ++i)
            {
                lambda[i] = s(i);
                mu(i) = u.col(i).transpose() * fooR * u.col(i);
            }

            #ifdef ADAPTIVE_ATTITUDE_COV_DEBUG_PRINTS
            std::cout<<"[ADAPTIVE_ATTITUDE] (lambda - mu) is:\n"<<(lambda - mu)<<"\n";
            #endif

            if ((lambda - mu).maxCoeff() > gamma)
            {

                #ifdef ADAPTIVE_ATTITUDE_COV_DEBUG_PRINTS
                std::cout<<"[ADAPTIVE_ATTITUDE] "<<(lambda - mu).maxCoeff() <<" Bigger than Gamma("<<gamma<<")\n";
                #endif

                r2count = 0;
                Eigen::Matrix<_Scalar, _MeasurementSize, 1> auxvector; /** Auxiliary vector variable */
                for(size_t i=0; i<_MeasurementSize; ++i)
                {
                    auxvector(i) = std::max(lambda(i)-mu(i),static_cast<double>(0.00));
                }

                Qstar = auxvector(0) * u.col(0) * u.col(0).transpose() + auxvector(1) * u.col(1) * u.col(1).transpose() + auxvector(2) * u.col(2) * u.col(2).transpose();
            }
            else
            {
                #ifdef ADAPTIVE_ATTITUDE_COV_DEBUG_PRINTS
                std::cout<<"[ADAPTIVE_ATTITUDE] "<<(lambda - mu).maxCoeff() <<" Lower than Gamma("<<gamma<<") r2count: "<<r2count<<"\n";
                #endif

                r2count ++;
                if (r2count < m2)
                {
                    Eigen::Matrix<_Scalar, _MeasurementSize, 1> auxvector; /** Auxiliary vector variable */
                    for(size_t i=0; i<_MeasurementSize; ++i)
                    {
                        auxvector(i) = std::max(lambda(i)-mu(i),static_cast<double>(0.00));
                    }

                    Qstar = auxvector(0) * u.col(0) * u.col(0).transpose() + auxvector(1) * u.col(1) * u.col(1).transpose() + auxvector(2) * u.col(2) * u.col(2).transpose();
                }
                else
                    Qstar = MeasurementMatrix::Zero();
            }

            #ifdef ADAPTIVE_ATTITUDE_COV_DEBUG_PRINTS
            std::cout<<"[ADAPTIVE_ATTITUDE] Qstar:\n"<<Qstar<<"\n";
            #endif

            return R + Qstar; //! R is the static and Qstar is the external acceleration covariance
        }
    };
}
#endif
