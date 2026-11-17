// gsc2ab: turns compact survivor files of gpusieve -compact (.gsc, see gsc.h) into an .ab file, or checks them.
//
//   gsc2ab [-o out.ab] in.gsc ...      the pairs of all inputs, in their order, as one .ab file (-o -: standard output)
//   gsc2ab -info in.gsc ...            per file: special-q blocks, first and last special-q, pairs, bits per pair
//
// The inputs must have the same side and I. They are read as streams (a block at a time), and every block's CRC is
// checked; a damaged or truncated input stops gsc2ab with exit status 2 (so a BOINC validator can reject the result),
// naming the file and the special-q after which it broke. The output is written to out.ab.part and renamed to out.ab only
// when every input was read completely, so a failed conversion leaves no output behind.
#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include "gsc.h"

// a file read in chunks: buf[pos..] holds the bytes not yet decoded
struct Reader {
	FILE* f = nullptr; std::vector<uint8_t> buf; size_t pos = 0; bool eof = false;
	bool open(const char* path) { f = strcmp(path, "-") ? fopen(path, "rb") : stdin; return f != nullptr; }
	void close() { if (f && f != stdin) fclose(f); f = nullptr; }
	bool more() {                            // read another chunk; false at the end of the file
		if (eof) return false;
		if (pos > buf.size() / 2) { buf.erase(buf.begin(), buf.begin() + pos); pos = 0; }
		const size_t old = buf.size(), want = std::max<size_t>(1 << 20, old);
		buf.resize(old + want);
		const size_t r = fread(buf.data() + old, 1, want, f);
		buf.resize(old + r);
		if (r == 0) eof = true;
		return r > 0;
	}
	size_t avail() const { return buf.size() - pos; }
};

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
	std::string part;
	if (outf) {
		if (strcmp(outf, "-")) { part = std::string(outf) + ".part"; out = fopen(part.c_str(), "wb"); }
		else out = stdout;
		if (!out) { perror(part.c_str()); return 1; }
	}
	auto fail = [&](int status) {          // no output from a failed conversion
		if (out && out != stdout) { fclose(out); remove(part.c_str()); }
		return status;
	};
	uint32_t side0 = ~0u, logI0 = 0;
	unsigned long long total = 0;
	std::vector<uint8_t> wbuf;
	for (const char* path : in) {
		Reader R;
		if (!R.open(path)) { perror(path); return fail(1); }
		while (R.avail() < 16 && R.more()) {}
		uint32_t h[4];
		if (R.avail() < 16 || memcmp(R.buf.data(), "GSC1", 4)) { fprintf(stderr, "%s: not a .gsc file\n", path); R.close(); return fail(2); }
		memcpy(h, R.buf.data(), 16); R.pos = 16;
		const uint32_t version = h[1], side = h[2], logI = h[3];
		if (version != 1 || side > 2 || logI < 10 || logI > 19) { fprintf(stderr, "%s: unknown version or header\n", path); R.close(); return fail(2); }
		if (side0 == ~0u) {
			side0 = side; logI0 = logI;
			if (out) { uint32_t ah[4] = {0, 1, side, 0}; memcpy(ah, "NFAB", 4); fwrite(ah, sizeof ah, 1, out); }
		} else if (side != side0 || logI != logI0) {
			fprintf(stderr, "%s: side or I differs from the first input\n", path); R.close(); return fail(1);
		}
		const uint64_t I = 1ull << logI, half = I >> 1;
		gsc::Block B;
		unsigned long long nblocks = 0, npairs = 0, nbytes = 16, q_first = 0, q_last = 0;
		for (;;) {
			const uint8_t* p = R.buf.data() + R.pos;
			const int rc = gsc::decode_block(p, R.buf.data() + R.buf.size(), B);
			if (rc == 1 || rc == 2) {                // the end of what was read: read on, unless the file has ended
				if (R.more()) continue;
				if (rc == 1) break;
			}
			if (rc != 0) {
				fprintf(stderr, "%s: %s after %llu blocks (last good special-q %llu)\n", path, rc == 2 ? "truncated" : "damaged",
					nblocks, q_last);
				R.close(); return fail(2);
			}
			nbytes += p - (R.buf.data() + R.pos);
			R.pos = p - R.buf.data();
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
		R.close();
		total += npairs;
		if (info) printf("%s: %llu special-q (%llu .. %llu), %llu pairs, %.1f bits per pair (side %u, I = 2^%u)\n", path, nblocks,
			q_first, q_last, npairs, npairs ? 8.0 * (nbytes - 16) / npairs : 0.0, side, logI);
	}
	if (out) {
		if (fflush(out) != 0) { perror(outf); return fail(1); }
		if (out != stdout) {
			if (fclose(out) != 0) { out = nullptr; remove(part.c_str()); perror(outf); return 1; }
			out = nullptr;
			if (rename(part.c_str(), outf) != 0) { perror(outf); remove(part.c_str()); return 1; }
		}
		fprintf(stderr, "%llu pairs written\n", total);
	}
	return 0;
}
