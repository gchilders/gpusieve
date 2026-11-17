// gsc2ab: turns compact survivor files of gpusieve -compact (.gsc, see gsc.h) into an .ab file, or checks them.
//
//   gsc2ab [-o out.ab] in.gsc ...      the pairs of all inputs, in their order, as one .ab file (-o -: standard output)
//   gsc2ab -info in.gsc ...            per file: special-q blocks, first and last special-q, pairs, bits per pair
//
// The inputs must have the same side and I. Every block's CRC is checked; a damaged or truncated file stops gsc2ab with
// exit status 2 (so a BOINC validator can reject the result), naming the file and the special-q after which it broke.
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include "gsc.h"

static bool read_file(const char* path, std::vector<uint8_t>& d)
{
	FILE* f = strcmp(path, "-") ? fopen(path, "rb") : stdin;
	if (!f) { perror(path); return false; }
	d.clear();
	uint8_t buf[1 << 16]; size_t r;
	while ((r = fread(buf, 1, sizeof buf, f)) > 0) d.insert(d.end(), buf, buf + r);
	if (f != stdin) fclose(f);
	return true;
}

int main(int argc, char** argv)
{
	const char* outf = nullptr; bool info = false;
	std::vector<const char*> in;
	for (int k = 1; k < argc; k++) {
		if (!strcmp(argv[k], "-o") && k + 1 < argc) outf = argv[++k];
		else if (!strcmp(argv[k], "-info")) info = true;
		else in.push_back(argv[k]);
	}
	if (in.empty() || (!info && !outf)) {
		fprintf(stderr, "usage: gsc2ab [-o out.ab | -o -] in.gsc ...   or   gsc2ab -info in.gsc ...\n");
		return 1;
	}
	FILE* out = nullptr;
	if (!info) {
		out = strcmp(outf, "-") ? fopen(outf, "wb") : stdout;
		if (!out) { perror(outf); return 1; }
	}
	uint32_t side0 = ~0u, logI0 = 0;
	unsigned long long total = 0;
	std::vector<uint8_t> d; std::vector<uint8_t> wbuf;
	for (const char* path : in) {
		if (!read_file(path, d)) return 1;
		uint32_t h[4];
		if (d.size() < 16 || memcmp(d.data(), "GSC1", 4)) { fprintf(stderr, "%s: not a .gsc file\n", path); return 2; }
		memcpy(h, d.data(), 16);
		const uint32_t version = h[1], side = h[2], logI = h[3];
		if (version != 1 || side > 2 || logI < 10 || logI > 19) { fprintf(stderr, "%s: unknown version or header\n", path); return 2; }
		if (side0 == ~0u) {
			side0 = side; logI0 = logI;
			if (out) { uint32_t ah[4] = {0, 1, side, 0}; memcpy(ah, "NFAB", 4); fwrite(ah, sizeof ah, 1, out); }
		} else if (side != side0 || logI != logI0) { fprintf(stderr, "%s: side or I differs from the first input\n", path); return 1; }
		const uint64_t I = 1ull << logI, half = I >> 1;
		const uint8_t* p = d.data() + 16; const uint8_t* end = d.data() + d.size();
		gsc::Block B;
		unsigned long long nblocks = 0, npairs = 0, q_first = 0, q_last = 0;
		for (;;) {
			const int rc = gsc::decode_block(p, end, B);
			if (rc == 1) break;
			if (rc < 0) {
				fprintf(stderr, "%s: damaged or truncated after %llu blocks (last good special-q %llu)\n", path, nblocks, q_last);
				return 2;
			}
			if (!nblocks) q_first = B.q;
			q_last = B.q; nblocks++; npairs += B.keys.size();
			if (out) {
				wbuf.resize(B.keys.size() * 12); uint8_t* w = wbuf.data();
				for (uint64_t key : B.keys) {
					const int64_t j = (int64_t)(key >> logI), i = (int64_t)(key & (I - 1)) - (int64_t)half;
					int64_t a = i * B.a0 + j * B.a1, b = i * B.b0 + j * B.b1;
					if (b < 0) { a = -a; b = -b; }
					const uint32_t b32 = (uint32_t)b; memcpy(w, &a, 8); memcpy(w + 8, &b32, 4); w += 12;
				}
				fwrite(wbuf.data(), 1, w - wbuf.data(), out);
			}
		}
		total += npairs;
		if (info) printf("%s: %llu special-q (%llu .. %llu), %llu pairs, %.1f bits per pair (side %u, I = 2^%u)\n", path, nblocks,
			q_first, q_last, npairs, npairs ? 8.0 * (d.size() - 16) / npairs : 0.0, side, logI);
	}
	if (out && (fflush(out) != 0 || (out != stdout && fclose(out) != 0))) { perror(outf); return 1; }
	if (!info) fprintf(stderr, "%llu pairs written\n", total);
	return 0;
}
