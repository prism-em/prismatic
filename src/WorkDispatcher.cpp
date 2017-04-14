//
// Created by AJ Pryor on 3/13/17.
//

#include "WorkDispatcher.h"
#include <mutex>
// helper function for dispatching work to various CPU/GPU worker threads

namespace PRISM {
        WorkDispatcher::WorkDispatcher(size_t _current,
					   size_t _stop,
					   size_t _num_per_call) :
					   current(_current),
					   stop(_stop),
					   num_per_call(_num_per_call) {};

		bool WorkDispatcher::getWork(size_t& job_start, size_t& job_stop){
			std::lock_guard<std::mutex> gatekeeper(lock);
			if (job_start >= stop) return false; // all jobs done, terminate
			job_start = current;
			job_stop = std::min(stop, current + num_per_call);
			current = job_stop;
			return true;
		}
}
