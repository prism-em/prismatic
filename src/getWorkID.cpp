//
// Created by AJ Pryor on 3/13/17.
//

#include "getWorkID.h"
#include <mutex>
// helper function for dispatching work

static size_t N_current;
static size_t N_total;
static size_t N_per_call;

void setWorkStartStop(const size_t& new_N_current, const size_t& new_N_total, const size_t& new_N_per_call){
	N_current  = new_N_current;
	N_total    = new_N_total;
	N_per_call = new_N_per_call;
}
bool getWorkID(const PRISM::Parameters<PRISM_FLOAT_PRECISION>& pars, size_t& Nstart, size_t& Nstop){
	static std::mutex lock; // mutex to synchronize reading/incrementing job ID
	std::lock_guard<std::mutex> supervisor(lock);
	if (N_current >= N_total) return false; // all jobs done, terminate
	Nstart = N_current;
	Nstop = std::min(N_total, Nstart + N_per_call);
	N_current = Nstop;
	return true;
}