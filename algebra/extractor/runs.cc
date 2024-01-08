#include <string>   // string
#include "runs.h"


struct PMR
{
    size_t start;
    size_t period;
    size_t count;
    size_t remainder;
}; // PMR


extern "C" {


size_t
find_runs(char const word[], PMR pmrs[])
{
    std::string const cword = "!" + std::string(word) + "!";

    auto const runs = linear_time_runs::compute_all_runs(cword.data(), cword.size());

    size_t idx = 0;
    for (auto run : runs)
    {
        pmrs[idx].start = run.start - 1;
        pmrs[idx].period = run.period;
        pmrs[idx].count = run.length / run.period;
        pmrs[idx].remainder = run.length % run.period;
        idx += 1;
    } // for

    return runs.size();
} // find_runs


} // extern "C"
