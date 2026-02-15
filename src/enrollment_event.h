#ifndef __ENROLLMENT_EVENT__
#define __ENROLLMENT_EVENT__

struct FlatMatrix;

#include <vector>

std::vector<double> accrual(const std::vector<double>& time,
                            const std::vector<double>& accrualTime,
                            const std::vector<double>& accrualIntensity,
                            const double accrualDuration);

std::vector<double> getAccrualDurationFromN(
    const std::vector<double>& nsubjects,
    const std::vector<double>& accrualTime,
    const std::vector<double>& accrualIntensity);

std::vector<double> patrisk(const std::vector<double>& time,
                            const std::vector<double>& piecewiseSurvivalTime,
                            const std::vector<double>& lambda,
                            const std::vector<double>& gamma);

#endif // __ENROLLMENT_EVENT__

