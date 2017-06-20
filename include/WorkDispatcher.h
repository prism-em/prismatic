// Copyright Alan (AJ) Pryor, Jr. 2017
// Transcribed from MATLAB code by Colin Ophus
// PRISM is distributed under the GNU General Public License (GPL)
// If you use PRISM, we ask that you cite the following papers:

#ifndef PRISM_WORKDISPATCHER_H
#define PRISM_WORKDISPATCHER_H
#include "params.h"
#include "configure.h"
#include <mutex>
namespace Prismatic {
    class WorkDispatcher {
    public:
        WorkDispatcher(size_t _current,
                       size_t _stop);

        bool getWork(size_t& job_start, size_t& job_stop, size_t num_requested=1, size_t cpu_early_stop=SIZE_MAX);
    private:
        std::mutex lock;
        size_t current, stop;
    };
}
#endif //PRISM_WORKDISPATCHER_H
