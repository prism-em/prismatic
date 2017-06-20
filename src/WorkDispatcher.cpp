// Copyright Alan (AJ) Pryor, Jr. 2017
// Transcribed from MATLAB code by Colin Ophus
// PRISM is distributed under the GNU General Public License (GPL)
// If you use PRISM, we ask that you cite the following papers:

#include "WorkDispatcher.h"
#include <mutex>
// helper function for dispatching work

namespace Prismatic {
        WorkDispatcher::WorkDispatcher(size_t _current,
					   size_t _stop) :
					   current(_current),
					   stop(_stop){};

		bool WorkDispatcher::getWork(size_t& job_start, size_t& job_stop, size_t num_requested, size_t early_cpu_stop){
			std::lock_guard<std::mutex> gatekeeper(lock);
			if (job_start >= stop | current>=early_cpu_stop) return false; // all jobs done, terminate
			job_start = current;
			job_stop = std::min(stop, current + num_requested);
			current = job_stop;
			return true;
		}
}
