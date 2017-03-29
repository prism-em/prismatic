//
// Created by AJ Pryor on 3/13/17.
//

#ifndef PRISM_WORKDISPATCHER_H
#define PRISM_WORKDISPATCHER_H
#include "params.h"
#include "configure.h"
#include <mutex>
namespace PRISM {
    class WorkDispatcher {
    public:
        WorkDispatcher(size_t _current,
                       size_t _stop,
                       size_t _num_per_call);

        bool getWork(size_t& job_start, size_t& job_stop);
    private:
        std::mutex lock;
        size_t current, stop, num_per_call;
    };
}
#endif //PRISM_WORKDISPATCHER_H
