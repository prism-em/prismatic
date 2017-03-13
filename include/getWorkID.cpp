//
// Created by AJ Pryor on 3/13/17.
//

#include "getWorkID.h"
#include "params.h"
#include "configure.h"
#include <mutex>
// helper function for dispatching work

static size_t N_current;
static size_t N_total;
void setWorkStartStop(const size_t& new_N_current, const size_t& new_N_total){
	N_current = new_N_current;
	N_total   = new_N_total;
}
bool getWorkID_probePos(const PRISM::Parameters<PRISM_FLOAT_PRECISION>& pars, size_t& Nstart, size_t& Nstop){
	static std::mutex lock; // mutex to synchronize reading/incrementing job ID
//	static size_t N_curddrent = 0; // number of next job
//	static const size_t N_total = pars.xp.size() * pars.yp.size(); // total number of jobs
	constexpr size_t NUM_JOBS_PER_CALL = 5; // number of pieces of work to assign to each calling thread
	std::lock_guard<std::mutex> supervisor(lock);
	if (N_current >= N_total) return false; // all jobs done, terminate
	Nstart = N_current;
	Nstop = std::min(N_total, Nstart + NUM_JOBS_PER_CALL);
	N_current = Nstop;
//	std::cout << " N_current = " << N_current<< std::endl;
	return true;
}