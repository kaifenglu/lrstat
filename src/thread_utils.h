#ifndef THREAD_UTILS_H
#define THREAD_UTILS_H

#include <Rcpp.h>
#include <RcppThread.h>
#include <mutex>
#include <string>
#include <vector>

namespace thread_utils {

// Shared mutex and vector (C++17 inline functions ensure single definition across TUs)
inline std::mutex& thread_warning_mutex() {
  static std::mutex m;
  return m;
}

inline std::vector<std::string>& thread_warning_msgs() {
  static std::vector<std::string> v;
  return v;
}

// Called from any thread (worker or main). Records a warning and prints a
// thread-safe console message.
inline void push_thread_warning(const std::string& msg) {
  {
    std::lock_guard<std::mutex> lock(thread_warning_mutex());
    thread_warning_msgs().push_back(msg);
  }
  RcppThread::Rcerr << msg << '\n';
}

// Console-only message (no recording)
inline void push_thread_message_console(const std::string& msg) {
  RcppThread::Rcerr << msg << '\n';
}

// Must be called on the main (R) thread after parallel work finishes.
// Emits collected messages as real R warnings (so they appear in warnings()).
inline void drain_thread_warnings_to_R() {
  std::vector<std::string> local;
  {
    std::lock_guard<std::mutex> lock(thread_warning_mutex());
    if (thread_warning_msgs().empty()) return;
    local.swap(thread_warning_msgs());
  }
  for (auto &m : local) {
    Rcpp::warning("%s", m.c_str());
  }
}

} // namespace thread_utils

#endif // THREAD_UTILS_H