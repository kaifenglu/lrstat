#ifndef __ENROLLMENT_EVENT__
#define __ENROLLMENT_EVENT__

struct FlatMatrix;

#include <vector>

double accrual1(
    const double time,
    const std::vector<double>& accrualTime,
    const std::vector<double>& accrualIntensity,
    const double accrualDuration);

std::vector<double> accrual(
    const std::vector<double>& time,
    const std::vector<double>& accrualTime,
    const std::vector<double>& accrualIntensity,
    const double accrualDuration);

double getAccrualDurationFromN1(
    const double nsubjects,
    const std::vector<double>& accrualTime,
    const std::vector<double>& accrualIntensity);

std::vector<double> getAccrualDurationFromN(
    const std::vector<double>& nsubjects,
    const std::vector<double>& accrualTime,
    const std::vector<double>& accrualIntensity);

double patrisk1(
    const double time,
    const std::vector<double>& piecewiseSurvivalTime,
    const std::vector<double>& lambda,
    const std::vector<double>& gamma);

std::vector<double> patrisk(
    const std::vector<double>& time,
    const std::vector<double>& piecewiseSurvivalTime,
    const std::vector<double>& lambda,
    const std::vector<double>& gamma);

double pevent1(
    const double time,
    const std::vector<double>& piecewiseSurvivalTime,
    const std::vector<double>& lambda,
    const std::vector<double>& gamma);

std::vector<double> pevent(
    const std::vector<double>& time,
    const std::vector<double>& piecewiseSurvivalTime,
    const std::vector<double>& lambda,
    const std::vector<double>& gamma);

std::pair<double, double> natrisk1cpp(
    const double time,
    const double allocationRatioPlanned,
    const std::vector<double>& accrualTime,
    const std::vector<double>& accrualIntensity,
    const std::vector<double>& piecewiseSurvivalTime,
    const std::vector<double>& lambda1,
    const std::vector<double>& lambda2,
    const std::vector<double>& gamma1,
    const std::vector<double>& gamma2,
    const double accrualDuration,
    const double minFollowupTime,
    const double maxFollowupTime);

FlatMatrix natriskcpp(
    const std::vector<double>& time,
    const double allocationRatioPlanned,
    const std::vector<double>& accrualTime,
    const std::vector<double>& accrualIntensity,
    const std::vector<double>& piecewiseSurvivalTime,
    const std::vector<double>& lambda1,
    const std::vector<double>& lambda2,
    const std::vector<double>& gamma1,
    const std::vector<double>& gamma2,
    const double accrualDuration,
    const double minFollowupTime,
    const double maxFollowupTime);

std::pair<double, double> nevent1cpp(
    const double time,
    const double allocationRatioPlanned,
    const std::vector<double>& accrualTime,
    const std::vector<double>& accrualIntensity,
    const std::vector<double>& piecewiseSurvivalTime,
    const std::vector<double>& lambda1,
    const std::vector<double>& lambda2,
    const std::vector<double>& gamma1,
    const std::vector<double>& gamma2,
    const double accrualDuration,
    const double minFollowupTime,
    const double maxFollowupTime);

FlatMatrix neventcpp(
    const std::vector<double>& time,
    const double allocationRatioPlanned,
    const std::vector<double>& accrualTime,
    const std::vector<double>& accrualIntensity,
    const std::vector<double>& piecewiseSurvivalTime,
    const std::vector<double>& lambda1,
    const std::vector<double>& lambda2,
    const std::vector<double>& gamma1,
    const std::vector<double>& gamma2,
    const double accrualDuration,
    const double minFollowupTime,
    const double maxFollowupTime);

std::pair<double, double> nevent21cpp(
    const double time,
    const double allocationRatioPlanned,
    const std::vector<double>& accrualTime,
    const std::vector<double>& accrualIntensity,
    const std::vector<double>& piecewiseSurvivalTime,
    const std::vector<double>& lambda1,
    const std::vector<double>& lambda2,
    const std::vector<double>& gamma1,
    const std::vector<double>& gamma2,
    const double accrualDuration,
    const double minFollowupTime,
    const double maxFollowupTime);

FlatMatrix nevent2cpp(const std::vector<double>& time,
                      const double allocationRatioPlanned,
                      const std::vector<double>& accrualTime,
                      const std::vector<double>& accrualIntensity,
                      const std::vector<double>& piecewiseSurvivalTime,
                      const std::vector<double>& lambda1,
                      const std::vector<double>& lambda2,
                      const std::vector<double>& gamma1,
                      const std::vector<double>& gamma2,
                      const double accrualDuration,
                      const double minFollowupTime,
                      const double maxFollowupTime);

#endif // __ENROLLMENT_EVENT__

