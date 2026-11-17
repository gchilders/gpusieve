// BOINC runtime glue for gpusieve, see boinc_glue.h.
#include "boinc_glue.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <chrono>
#include <string>
#include <deque>

// a checkpoint every standalone_s seconds, for runs without a client
static bool standalone_checkpoint_due(double standalone_s)
{
	static auto last = std::chrono::steady_clock::now();
	const auto now = std::chrono::steady_clock::now();
	if (std::chrono::duration<double>(now - last).count() < standalone_s) return false;
	last = now;
	return true;
}

#ifdef HAVE_BOINC
#include "boinc_api.h"

static bool g_up = false;              // boinc_init_options succeeded
static double g_last_fraction = 0;

int gs_boinc_init()
{
	BOINC_OPTIONS opt;
	boinc_options_defaults(opt);
	opt.normal_thread_priority = 1;    // GPU app: keep the host thread that feeds the GPU at normal priority
	int rc = boinc_init_options(&opt);
	if (rc) { fprintf(stderr, "BOINC: boinc_init_options failed (%d)\n", rc); return rc; }
	g_up = true;
	boinc_fraction_done(0);
	fprintf(stderr, "BOINC: API initialised (%s)\n", boinc_is_standalone() ? "standalone" : "under a client");
	return 0;
}

bool gs_boinc_managed() { return g_up && !boinc_is_standalone(); }

int gs_boinc_gpu_device()
{
	if (!g_up) return -1;
	// The assignment is in init_data.xml. The client also appends "--device N", but only for app versions with an
	// api_version below 7.5, so gpusieve takes init_data.xml first and the command line only as a fallback.
	APP_INIT_DATA aid;
	boinc_get_init_data(aid);
	if (aid.gpu_device_num < 0) return -1;
	if (aid.gpu_type[0] && strncmp(aid.gpu_type, "NVIDIA", 6) != 0) {
		fprintf(stderr, "BOINC: this is a CUDA application, but the client assigned a '%s' device (%d); ignored\n", aid.gpu_type,
			aid.gpu_device_num);
		return -1;
	}
	return aid.gpu_device_num;
}

double gs_boinc_ncpus()
{
	if (!g_up) return 0;
	APP_INIT_DATA aid;
	boinc_get_init_data(aid);
	return aid.ncpus > 0 ? aid.ncpus : 0;
}

const char* gs_boinc_resolve(const char* logical)
{
	static std::deque<std::string> keep;     // the returned pointers must stay valid
	if (!g_up || !logical) return logical;
	std::string phys;
	if (boinc_resolve_filename_s(logical, phys) || phys.empty()) {
		fprintf(stderr, "BOINC: cannot resolve the file name '%s'\n", logical);
		gs_boinc_finish(1);
	}
	keep.push_back(phys);
	return keep.back().c_str();
}

void gs_boinc_fraction_done(double f)
{
	if (!g_up || !(f >= 0)) return;
	if (f > 1) f = 1;
	if (f <= g_last_fraction) return;      // the client wants nondecreasing values
	g_last_fraction = f;
	boinc_fraction_done(f);
}

bool gs_boinc_time_to_checkpoint(double standalone_s)
{
	if (gs_boinc_managed()) return boinc_time_to_checkpoint() != 0;
	return standalone_checkpoint_due(standalone_s);
}

void gs_boinc_checkpoint_completed() { if (gs_boinc_managed()) boinc_checkpoint_completed(); }   // pairs with boinc_time_to_checkpoint()

void gs_boinc_temporary_exit(int delay_s, const char* reason)
{
	if (!gs_boinc_managed()) return;
	fprintf(stderr, "BOINC: temporary exit for %d s: %s\n", delay_s, reason);
	fflush(stderr); fflush(stdout);
	boinc_temporary_exit(delay_s, reason, false);   // does not return
}

void gs_boinc_finish(int status)
{
	fflush(stderr); fflush(stdout);
	if (g_up) {
		if (status == 0) boinc_fraction_done(1.0);
		boinc_finish(status);              // does not return
	}
	exit(status);
}

#else   // no BOINC: standalone behaviour

int gs_boinc_init() { return 0; }
bool gs_boinc_managed() { return false; }
int gs_boinc_gpu_device() { return -1; }
double gs_boinc_ncpus() { return 0; }
const char* gs_boinc_resolve(const char* logical) { return logical; }
void gs_boinc_fraction_done(double) {}
bool gs_boinc_time_to_checkpoint(double standalone_s) { return standalone_checkpoint_due(standalone_s); }
void gs_boinc_checkpoint_completed() {}
void gs_boinc_temporary_exit(int, const char*) {}
void gs_boinc_finish(int status) { fflush(stderr); fflush(stdout); exit(status); }

#endif
