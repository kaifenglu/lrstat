#ifndef __SURVIVAL_ANALYSIS__
#define __SURVIVAL_ANALYSIS__

struct FlatMatrix;
struct DataFrameCpp;
struct ListCpp;

#include <string>
#include <vector>

DataFrameCpp survQuantilecpp(const std::vector<double>& time,
                             const std::vector<int>& event,
                             const double cilevel = 0.95,
                             const std::string& transform = "loglog",
                             const std::vector<double>& probs = {0.25, 0.5, 0.75});

DataFrameCpp kmestcpp(const DataFrameCpp& data,
                      const std::vector<std::string>& stratum,
                      const std::string& time = "time",
                      const std::string& time2 = "",
                      const std::string& event = "event",
                      const std::string& weight = "",
                      const std::string& conftype = "log-log",
                      const double conflev = 0.95,
                      const bool keep_censor = false);

DataFrameCpp kmdiffcpp(const DataFrameCpp& data,
                       const std::vector<std::string>& stratum,
                       const std::string& treat = "treat",
                       const std::string& time = "time",
                       const std::string& time2 = "",
                       const std::string& event = "event",
                       const std::string& weight = "",
                       const double milestone = 0,
                       const double survDiffH0 = 0,
                       const double conflev = 0.95);

DataFrameCpp lrtestcpp(const DataFrameCpp& data,
                       const std::vector<std::string>& stratum,
                       const std::string& treat = "treat",
                       const std::string& time = "time",
                       const std::string& time2 = "",
                       const std::string& event = "event",
                       const std::string& weight = "",
                       const bool weight_readj = false,
                       const double rho1 = 0,
                       const double rho2 = 0);

DataFrameCpp rmestcpp(const DataFrameCpp& data,
                      const std::vector<std::string>& stratum,
                      const std::string& time = "time",
                      const std::string& event = "event",
                      const double milestone = 0,
                      const double conflev = 0.95,
                      const bool biascorrection = false);

DataFrameCpp rmdiffcpp(const DataFrameCpp& data,
                       const std::vector<std::string>& stratum,
                       const std::string& treat = "treat",
                       const std::string& time = "time",
                       const std::string& event = "event",
                       const double milestone = 0,
                       const double rmstDiffH0 = 0,
                       const double conflev = 0.95,
                       const bool biascorrection = false);

ListCpp liferegcpp(const DataFrameCpp& data,
                   const std::vector<std::string>& stratum,
                   const std::string time = "time",
                   const std::string time2 = "",
                   const std::string event = "event",
                   const std::vector<std::string>& covariates = {},
                   const std::string weight = "",
                   const std::string offset = "",
                   const std::string id = "",
                   const std::string dist = "weibull",
                   const std::vector<double>& init = {},
                   const bool robust = false,
                   const bool plci = false,
                   const double alpha = 0.05,
                   const int maxiter = 50,
                   const double eps = 1.0e-9);

FlatMatrix residuals_liferegcpp(const std::vector<double>& beta,
                                const FlatMatrix& vbeta,
                                const DataFrameCpp& data,
                                const std::vector<std::string>& stratum,
                                const std::string& time = "time",
                                const std::string& time2 = "",
                                const std::string& event = "event",
                                const std::vector<std::string>& covariates = {},
                                const std::string& weight = "",
                                const std::string& offset = "",
                                const std::string& id = "",
                                const std::string& dist = "weibull",
                                const std::string& type = "response",
                                bool collapse = false,
                                bool weighted = false);

ListCpp phregcpp(const DataFrameCpp& data,
                 const std::vector<std::string>& stratum,
                 const std::string& time = "time",
                 const std::string& time2 = "",
                 const std::string& event = "event",
                 const std::vector<std::string>& covariates = {},
                 const std::string& weight = "",
                 const std::string& offset = "",
                 const std::string& id = "",
                 const std::string& ties = "efron",
                 const std::vector<double>& init = {},
                 const bool robust = false,
                 const bool est_basehaz = true,
                 const bool est_resid = true,
                 const bool firth = false,
                 const bool plci = false,
                 const double alpha = 0.05,
                 const int maxiter = 50,
                 const double eps = 1.0e-9);

DataFrameCpp survfit_phregcpp(const int p,
                              const std::vector<double>& beta,
                              const FlatMatrix& vbeta,
                              const DataFrameCpp& basehaz,
                              const DataFrameCpp& newdata,
                              const std::vector<std::string>& covariates,
                              const std::vector<std::string>& stratum,
                              const std::string& offset = "",
                              const std::string& id = "",
                              const std::string& tstart = "",
                              const std::string& tstop = "",
                              const bool sefit = true,
                              const std::string& conftype = "log-log",
                              const double conflev = 0.95);

ListCpp residuals_phregcpp(const int p,
                           const std::vector<double>& beta,
                           const FlatMatrix& vbeta,
                           const std::vector<double>& resmart,
                           const DataFrameCpp& data,
                           const std::vector<std::string>& stratum,
                           const std::string& time = "time",
                           const std::string& time2 = "",
                           const std::string& event = "event",
                           const std::vector<std::string>& covariates = {},
                           const std::string& weight = "",
                           const std::string& offset = "",
                           const std::string& id = "",
                           const std::string& ties = "efron",
                           const std::string& type = "schoenfeld",
                           const bool collapse = false,
                           const bool weighted = false);

ListCpp assess_phregcpp(const int p,
                        const std::vector<double>& beta,
                        const FlatMatrix& vbeta,
                        const DataFrameCpp& data,
                        const std::vector<std::string>& stratum ,
                        const std::string& time = "time",
                        const std::string& time2 = "",
                        const std::string& event = "event",
                        const std::vector<std::string>& covariates = {},
                        const std::string& weight = "",
                        const std::string& offset = "",
                        const std::string& ties = "efron",
                        const int resample = 1000,
                        const int seed = 0);

ListCpp zph_phregcpp(int p,
                     const std::vector<double>& beta,
                     const FlatMatrix& vbeta,
                     const std::vector<double>& resmart,
                     const DataFrameCpp& data,
                     const std::vector<std::string>& stratum,
                     const std::string& time = "time",
                     const std::string& time2 = "",
                     const std::string& event = "event",
                     const std::vector<std::string>& covariates = {},
                     const std::string& weight = "",
                     const std::string& offset = "",
                     const std::string& ties = "efron",
                     const std::string& transform = "km");

#endif // __SURVIVAL_ANALYSIS__
