#ifndef INC_GIT_VERSION_hpp
#define INC_GIT_VERSION_hpp

extern const char *git_build_time;
extern const char *git_build_sha;
extern const char *git_build_branch;
extern const char *git_log_stat;

#include <string>

inline std::string FixedLength(const char *s)
{
    std::string tmp(s);
    tmp.resize(41, ' ');
    return tmp;
}

#define Log_Git_Info(basename)                                                                  \
{                                                                                               \
    std::string git_version("##############################################################\n");\
    git_version +=          "#  libpotential's Git versioning:                            #\n"; \
    git_version +=          "#    Branch:       " + FixedLength(git_build_branch)     + " #\n"; \
    git_version +=          "#    Commit id:    " + FixedLength(git_build_sha)        + " #\n"; \
    git_version +=          "#    Build time:   " + FixedLength(git_build_time)       + " #\n"; \
    git_version +=          "##############################################################\n"; \
                                                                                                \
    std_cout << git_version << "\n";                                                            \
                                                                                                \
    git_version += git_log_stat;                                                                \
                                                                                                \
    std::string filename = basename + "/git_version_libpotentials.log";                         \
    std::ofstream git_file;                                                                     \
    git_file.open(filename.c_str());                                                            \
    git_file << git_version << "\n";                                                            \
    git_file.close();                                                                           \
                                                                                                \
    std::string gitdiff = reinterpret_cast<const char*>(src_Git_Diff_patch);                    \
    filename = basename + "/git_diff_libpotentials.patch";                                      \
    git_file.open(filename.c_str(), std::fstream::out);                                         \
    git_file << gitdiff << "\n";                                                                \
    git_file.close();                                                                           \
}

#endif // INC_GIT_VERSION_hpp

// ********** End of file ***************************************
