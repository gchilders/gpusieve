// BOINC runtime glue for gpusieve. Compiled into every build; without HAVE_BOINC the functions are no-ops (or report "not
// managed"), so gpusieve.cu needs no #ifdef around its calls. With HAVE_BOINC they wrap the BOINC API: initialisation, the
// client's GPU assignment, logical file names, fraction done, checkpoint timing and the exits.
#pragma once

// boinc_init_options with normal host-thread priority (GPU apps must not be starved by the CPU tasks); 0 on success.
int gs_boinc_init();
// true when the BOINC API is up and a client manages this run (false standalone and in non-BOINC builds)
bool gs_boinc_managed();
// the CUDA device the client assigned (init_data.xml gpu_device_num), or -1 for none / not managed
int gs_boinc_gpu_device();
// ncpus the client gives this task (init_data.xml), or 0 if unknown
double gs_boinc_ncpus();
// physical path of a logical file name (the name itself when not managed); the string stays valid for the whole run
const char* gs_boinc_resolve(const char* logical);
// monotone fraction done in [0, 1]
void gs_boinc_fraction_done(double f);
// true when a checkpoint is due: the client's checkpoint interval when managed, else every `standalone_s` seconds
bool gs_boinc_time_to_checkpoint(double standalone_s);
void gs_boinc_checkpoint_completed();
// transient failure (GPU busy, out of memory): ask the client to restart the task after `delay_s` seconds. Returns only
// when not managed; the caller then exits with an error itself.
void gs_boinc_temporary_exit(int delay_s, const char* reason);
// end of the run: boinc_finish(status) when the API is up, else exit(status). Does not return.
[[noreturn]] void gs_boinc_finish(int status);
