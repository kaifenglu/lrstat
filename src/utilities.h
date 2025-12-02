#include <Rcpp.h>
#include <R_ext/Applic.h>

using namespace Rcpp;

#ifndef __UTILITIES__
#define __UTILITIES__

void set_seed(int seed);

NumericVector stl_sort(const NumericVector& x);

IntegerVector which(const LogicalVector& vector);

IntegerVector findInterval3(NumericVector x, NumericVector v,
                            bool rightmost_closed,
                            bool all_inside,
                            bool left_open);

double brent(const std::function<double(double)>& f,
             double x1, double x2, double tol);

double bisect(const std::function<double(double)>& f,
              double x1, double x2, double tol);

double errorSpentcpp(const double t,
                     const double error,
                     const String sf,
                     const double sfpar);

List exitprobcpp(const NumericVector& b,
                 const NumericVector& a,
                 const NumericVector& theta,
                 const NumericVector& I);

NumericVector dtpwexpcpp(const NumericVector& q,
                         const NumericVector& piecewiseSurvivalTime,
                         const NumericVector& lambda,
                         const double lowerBound,
                         const bool logd);

NumericVector ptpwexpcpp(const NumericVector& q,
                         const NumericVector& piecewiseSurvivalTime,
                         const NumericVector& lambda,
                         const double lowerBound,
                         const bool lowertail,
                         const bool logp);

double qtpwexpcpp1(const double p,
                   const NumericVector& piecewiseSurvivalTime,
                   const NumericVector& lambda,
                   const double lowerBound,
                   const bool lowertail,
                   const bool logp);

NumericVector qtpwexpcpp(const NumericVector& p,
                         const NumericVector& piecewiseSurvivalTime,
                         const NumericVector& lambda,
                         const double lowerBound,
                         const bool lowertail,
                         const bool logp);

List mtpwexpcpp(const NumericVector& piecewiseSurvivalTime,
                const NumericVector& lambda,
                const double lowerBound);

NumericVector rtpwexpcpp(
    const int n,
    const NumericVector& piecewiseSurvivalTime,
    const NumericVector& lambda,
    const double lowerBound);

NumericVector getBoundcpp(
    const int k,
    const NumericVector& informationRates,
    const double alpha,
    const String typeAlphaSpending,
    const double parameterAlphaSpending,
    const NumericVector& userAlphaSpending,
    const NumericVector& spendingTime,
    const LogicalVector& efficacyStopping);

List getPower(const double alpha,
              const int kMax,
              const NumericVector& b,
              const NumericVector& theta,
              const NumericVector& I,
              const std::string bsf,
              const double bsfpar,
              const NumericVector& st,
              const LogicalVector& futilityStopping,
              const NumericVector& w);

double intnorm(const std::function<double(double)>& f,
               double mu, double sigma, double a, double b);

NumericVector mini(const std::function<double(double)>& f,
                   double x1, double x2, double tol);

List quad(integr_fn f, void *ex, double lower, double upper,
          double tol);

List quad2d(const std::function<double(double,double)>& f,
            double ax, double bx,
            double ay, double by,
            double tol);

double pbvnormcpp(NumericVector lower, NumericVector upper, double rho);

List hazard_pdcpp(const NumericVector& piecewiseSurvivalTime,
                  const NumericVector& hazard_pfs,
                  const NumericVector& hazard_os,
                  const double rho_pd_os);

List hazard_subcpp(const NumericVector& piecewiseSurvivalTime,
                   const NumericVector& hazard_itt,
                   const NumericVector& hazard_pos,
                   const double p_pos);

List bmini(NumericVector x0, optimfn fn, optimgr gr, void *ex, double eps);

NumericVector accrual(const NumericVector& time,
                      const NumericVector& accrualTime,
                      const NumericVector& accrualIntensity,
                      const double accrualDuration);

NumericVector getAccrualDurationFromN(
    const NumericVector& nsubjects,
    const NumericVector& accrualTime,
    const NumericVector& accrualIntensity);

NumericVector patrisk(const NumericVector& time,
                      const NumericVector& piecewiseSurvivalTime,
                      const NumericVector& lambda,
                      const NumericVector& gamma);

NumericVector pevent(const NumericVector& time,
                     const NumericVector& piecewiseSurvivalTime,
                     const NumericVector& lambda,
                     const NumericVector& gamma);

NumericMatrix natrisk(const NumericVector& time,
                      const double allocationRatioPlanned,
                      const NumericVector& accrualTime,
                      const NumericVector& accrualIntensity,
                      const NumericVector& piecewiseSurvivalTime,
                      const NumericVector& lambda1,
                      const NumericVector& lambda2,
                      const NumericVector& gamma1,
                      const NumericVector& gamma2,
                      const double accrualDuration,
                      const double minFollowupTime,
                      const double maxFollowupTime);

NumericMatrix nevent(const NumericVector& time,
                     const double allocationRatioPlanned,
                     const NumericVector& accrualTime,
                     const NumericVector& accrualIntensity,
                     const NumericVector& piecewiseSurvivalTime,
                     const NumericVector& lambda1,
                     const NumericVector& lambda2,
                     const NumericVector& gamma1,
                     const NumericVector& gamma2,
                     const double accrualDuration,
                     const double minFollowupTime,
                     const double maxFollowupTime);

NumericMatrix nevent2(const NumericVector& time,
                      const double allocationRatioPlanned,
                      const NumericVector& accrualTime,
                      const NumericVector& accrualIntensity,
                      const NumericVector& piecewiseSurvivalTime,
                      const NumericVector& lambda1,
                      const NumericVector& lambda2,
                      const NumericVector& gamma1,
                      const NumericVector& gamma2,
                      const double accrualDuration,
                      const double minFollowupTime,
                      const double maxFollowupTime);

List getDesign(const double beta,
               const double IMax,
               const double theta,
               const int kMax,
               const NumericVector& informationRates,
               const LogicalVector& efficacyStopping,
               const LogicalVector& futilityStopping,
               const NumericVector& criticalValues,
               const double alpha,
               const std::string typeAlphaSpending,
               const double parameterAlphaSpending,
               const NumericVector& userAlphaSpending,
               const NumericVector& futilityBounds,
               const std::string typeBetaSpending,
               const double parameterBetaSpending,
               const NumericVector& userBetaSpending,
               const NumericVector& spendingTime,
               const double varianceRatio);

List getDesignEquiv(const double beta,
                    const double IMax,
                    const double thetaLower,
                    const double thetaUpper,
                    const double theta,
                    const int kMax,
                    const NumericVector& informationRates,
                    const NumericVector& criticalValues,
                    const double alpha,
                    const std::string typeAlphaSpending,
                    const double parameterAlphaSpending,
                    const NumericVector& userAlphaSpending,
                    const NumericVector& spendingTime);

List adaptDesign(double betaNew,
                 double INew,
                 const int L,
                 const double zL,
                 const double theta,
                 const double IMax,
                 const int kMax,
                 const NumericVector& informationRates,
                 const LogicalVector& efficacyStopping,
                 const LogicalVector& futilityStopping,
                 const NumericVector& criticalValues,
                 const double alpha,
                 const std::string typeAlphaSpending,
                 const double parameterAlphaSpending,
                 const NumericVector& userAlphaSpending,
                 const NumericVector& futilityBounds,
                 const std::string typeBetaSpending,
                 const double parameterBetaSpending,
                 const NumericVector& spendingTime,
                 const bool MullerSchafer,
                 const int kNew,
                 const NumericVector& informationRatesNew,
                 const LogicalVector& efficacyStoppingNew,
                 const LogicalVector& futilityStoppingNew,
                 const std::string typeAlphaSpendingNew,
                 const double parameterAlphaSpendingNew,
                 const std::string typeBetaSpendingNew,
                 const double parameterBetaSpendingNew,
                 const NumericVector& userBetaSpendingNew,
                 const NumericVector& spendingTimeNew,
                 const double varianceRatio);

bool hasVariable(DataFrame df, std::string varName);

double quantilecpp(const NumericVector& x, const double p);

double squantilecpp(const std::function<double(double)>& S, double p);

IntegerVector c_vectors_i(IntegerVector vec1, IntegerVector vec2);
NumericVector c_vectors(NumericVector vec1, NumericVector vec2);
NumericMatrix subset_matrix_by_row(NumericMatrix a, IntegerVector q);
NumericMatrix c_matrices(NumericMatrix a1, NumericMatrix a2);

List bygroup(DataFrame data, const StringVector& variables);

int cholesky2(NumericMatrix matrix, int n, double toler);
void chsolve2(NumericMatrix matrix, int n, NumericVector y);
void chinv2(NumericMatrix matrix, int n);
NumericMatrix invsympd(NumericMatrix matrix, int n, double toler);

DataFrame survsplit(NumericVector tstart,
                    NumericVector tstop,
                    NumericVector cut);

bool is_sorted(NumericVector x);

List qrcpp(const NumericMatrix& X, double tol);

List svdcpp(const NumericMatrix& X, const bool outtransform,
            const bool decreasing);

NumericMatrix rmvnormcpp(int n, NumericVector mean, NumericMatrix sigma);

NumericVector float_to_fraction(const double x, const double tol);

#endif // __UTILITIES__
