

inline std::string now_hms() {
  std::time_t t = std::time(nullptr);
  char buf[16];
  std::strftime(buf, sizeof(buf), "%H:%M:%S", std::localtime(&t));
  return std::string(buf);
}

struct Timer {
  std::chrono::steady_clock::time_point start;
  void begin() { start = std::chrono::steady_clock::now(); }
  std::tuple<int,int,int> hms() const {
    auto dur = std::chrono::duration_cast<std::chrono::seconds>(
      std::chrono::steady_clock::now() - start
    ).count();
    int h = static_cast<int>(dur / 3600);
    int m = static_cast<int>((dur - h*3600) / 60);
    int s = static_cast<int>(dur - h*3600 - m*60);
    return {h,m,s};
  }
};

inline void print_completed(const char* prefix, const Timer& tm) {
  auto [h,m,s] = tm.hms();
  Rcpp::Rcout << prefix << " completed in: " << h << "h " << m << "m " << s << "s.\n";
}
