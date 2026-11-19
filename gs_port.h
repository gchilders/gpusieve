// Portability of gpusieve's host code between GCC/Clang (Linux) and MSVC (nvcc's host compiler on Windows). MSVC has
// no __int128, and Windows has no ftruncate, fsync or fseeko (__builtin_clzll: gsc::clz64 in gsc.h). The 128-bit
// arithmetic is needed only for moduli of 2^32 and more and for the special-q lattice reduction, so the portable
// version below is plain code rather than intrinsics; build with -DGS_NO_INT128 to use it with GCC as well (that is how
// it is tested against __int128). Replacing a file over an existing one is gs_replace_file() in boinc_glue.cpp, which
// keeps <windows.h> out of gpusieve.cu.
#pragma once
#include <cstdint>
#include <cstdio>
#include <sys/types.h>
#include <sys/stat.h>
#ifdef _WIN32
#include <io.h>         // _chsize_s, _commit, _fileno
#else
#include <unistd.h>     // ftruncate, fsync
#endif
#if defined(_MSC_VER) || defined(GS_NO_INT128)
#define GS_PORTABLE_INT128 1
#endif

#ifndef GS_PORTABLE_INT128

typedef __int128 gs_i128;

// (a * b + c) mod m, for any 64-bit a, b, c and m > 0
inline uint64_t gs_muladdmod(uint64_t a, uint64_t b, uint64_t c, uint64_t m)
{
	return (uint64_t)(((unsigned __int128)a * b + c) % m);
}

#else

// the 128-bit product of a and b
inline void gs_mul64(uint64_t a, uint64_t b, uint64_t& hi, uint64_t& lo)
{
	const uint64_t a0 = (uint32_t)a, a1 = a >> 32, b0 = (uint32_t)b, b1 = b >> 32;
	const uint64_t p00 = a0 * b0, p01 = a0 * b1, p10 = a1 * b0, p11 = a1 * b1;
	const uint64_t mid = (p00 >> 32) + (uint32_t)p01 + (uint32_t)p10;
	lo = mid << 32 | (uint32_t)p00;
	hi = p11 + (p01 >> 32) + (p10 >> 32) + (mid >> 32);
}

// (hi * 2^64 + lo) mod m, one bit at a time (only used for moduli of 2^32 and more, which are rare)
inline uint64_t gs_mod128(uint64_t hi, uint64_t lo, uint64_t m)
{
	uint64_t r = hi % m;
	for (int i = 63; i >= 0; i--) {
		const bool carry = r >> 63;
		r = r << 1 | (lo >> i & 1);
		if (carry || r >= m) r -= m;      // 2r + bit < 2m: one subtraction (modulo 2^64 when it carried out)
	}
	return r;
}

inline uint64_t gs_muladdmod(uint64_t a, uint64_t b, uint64_t c, uint64_t m)
{
	uint64_t hi, lo;
	gs_mul64(a, b, hi, lo);
	lo += c;
	hi += lo < c;                         // a * b + c < 2^128
	return gs_mod128(hi, lo, m);
}

// a signed 128-bit integer (two's complement, wrapping) with the operations reduce_q() uses
struct gs_i128 {
	uint64_t lo = 0, hi = 0;
	gs_i128() {}
	gs_i128(int v) : lo((uint64_t)(int64_t)v), hi(v < 0 ? ~0ull : 0) {}
	gs_i128(long v) : lo((uint64_t)(int64_t)v), hi(v < 0 ? ~0ull : 0) {}
	gs_i128(long long v) : lo((uint64_t)v), hi(v < 0 ? ~0ull : 0) {}
	gs_i128(unsigned long v) : lo((uint64_t)v), hi(0) {}
	gs_i128(unsigned long long v) : lo((uint64_t)v), hi(0) {}
	explicit operator long long() const { return (long long)lo; }
	explicit operator long() const { return (long)(int64_t)lo; }
	bool negative() const { return (int64_t)hi < 0; }

	friend gs_i128 operator+(gs_i128 a, gs_i128 b) { gs_i128 r; r.lo = a.lo + b.lo; r.hi = a.hi + b.hi + (r.lo < a.lo); return r; }
	friend gs_i128 operator-(gs_i128 a) { gs_i128 r; r.lo = ~a.lo + 1; r.hi = ~a.hi + (r.lo == 0); return r; }
	friend gs_i128 operator-(gs_i128 a, gs_i128 b) { return a + -b; }
	friend gs_i128 operator*(gs_i128 a, gs_i128 b)
	{
		gs_i128 r;
		gs_mul64(a.lo, b.lo, r.hi, r.lo);
		r.hi += a.lo * b.hi + a.hi * b.lo;
		return r;
	}
	friend gs_i128 operator<<(gs_i128 a, int s)     // 0 <= s < 128
	{
		gs_i128 r;
		if (s == 0) return a;
		if (s >= 64) { r.hi = a.lo << (s - 64); r.lo = 0; }
		else { r.hi = a.hi << s | a.lo >> (64 - s); r.lo = a.lo << s; }
		return r;
	}
	static bool ult(gs_i128 a, gs_i128 b) { return a.hi != b.hi ? a.hi < b.hi : a.lo < b.lo; }
	// truncating division, as for the built-in integers; b != 0
	friend gs_i128 operator/(gs_i128 a, gs_i128 b)
	{
		const bool neg = a.negative() != b.negative();
		if (a.negative()) a = -a;
		if (b.negative()) b = -b;
		gs_i128 q, r;
		for (int i = 127; i >= 0; i--) {
			r = r << 1;
			r.lo |= (i >= 64 ? a.hi >> (i - 64) : a.lo >> i) & 1;
			if (!ult(r, b)) {
				r = r - b;
				if (i >= 64) q.hi |= 1ull << (i - 64); else q.lo |= 1ull << i;
			}
		}
		return neg ? -q : q;
	}
	friend bool operator<(gs_i128 a, gs_i128 b) { return a.hi != b.hi ? (int64_t)a.hi < (int64_t)b.hi : a.lo < b.lo; }
	friend bool operator>(gs_i128 a, gs_i128 b) { return b < a; }
	friend bool operator==(gs_i128 a, gs_i128 b) { return a.lo == b.lo && a.hi == b.hi; }
	friend bool operator!=(gs_i128 a, gs_i128 b) { return !(a == b); }
	gs_i128& operator+=(gs_i128 b) { return *this = *this + b; }
	gs_i128& operator-=(gs_i128 b) { return *this = *this - b; }
	gs_i128& operator<<=(int s) { return *this = *this << s; }
};

#endif

// file operations of the checkpointing (64-bit offsets on both systems)
inline int64_t gs_ftell(FILE* f)
{
#ifdef _WIN32
	return _ftelli64(f);
#else
	return (int64_t)ftello(f);
#endif
}

inline int gs_fseek(FILE* f, int64_t off)
{
#ifdef _WIN32
	return _fseeki64(f, off, SEEK_SET);
#else
	return fseeko(f, (off_t)off, SEEK_SET);
#endif
}

// truncate the open file to off bytes; 0 on success
inline int gs_ftruncate(FILE* f, int64_t off)
{
#ifdef _WIN32
	return _chsize_s(_fileno(f), off) == 0 ? 0 : -1;
#else
	return ftruncate(fileno(f), (off_t)off);
#endif
}

// flush the file's data to the disk (call fflush first); 0 on success, else -1 with errno set. On Windows a device or
// pipe (NUL, a console) has nothing to sync and gives 0, as fsync's EINVAL does on Linux, which the callers accept.
inline int gs_fsync(FILE* f)
{
#ifdef _WIN32
	struct _stat64 st;
	if (_fstat64(_fileno(f), &st) == 0 && !(st.st_mode & _S_IFREG)) return 0;
	return _commit(_fileno(f));
#else
	return fsync(fileno(f));
#endif
}

// the size of a file; false if it cannot be read
inline bool gs_file_size(const char* path, uint64_t& size)
{
#ifdef _WIN32
	struct _stat64 st;
	if (_stat64(path, &st) != 0) return false;
#else
	struct stat st;
	if (stat(path, &st) != 0) return false;
#endif
	size = (uint64_t)st.st_size;
	return true;
}
