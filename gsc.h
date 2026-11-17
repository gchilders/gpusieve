// Compact survivor files (.gsc), written by gpusieve -compact and turned into .ab files by gsc2ab.
//
// A 16-byte header: "GSC1", then version (1), the side of the .ab header to write (the side NFS_factory rebuilds) and logI,
// as little-endian 32-bit integers. Then one block per special-q, in the order sieved:
//
//   varint q, varint r, zigzag varints a0 b0 a1 b1       the special-q, its root and its reduced lattice basis
//   varint n, byte k, varint nbytes, nbytes of bits       n survivors as Golomb-Rice coded gaps (parameter k)
//   4 bytes CRC-32 of the block's bytes before it
//
// A survivor is the lattice point (i,j), (a,b) = i*(a0,b0) + j*(a1,b1) with the sign making b > 0. Its key is
// j*2^logI + (i + 2^(logI-1)); the keys are sorted, and the gaps key[0] and key[t] - key[t-1] - 1 are written as Rice
// codes: the quotient gap >> k in unary (ones, then a zero; 32 ones are followed by an Elias-gamma code of quotient - 31),
// then the low k bits. Survivors are about uniform in the sieve region, so a gap costs about log2(region / n) + 1.5 bits:
// 20 bits at I = 16 with 8000 survivors per special-q, against 96 bits in an .ab file. Blocks are byte-aligned, so a
// checkpoint can truncate the file after any of them.
#pragma once
#include <cstdint>
#include <cstring>
#include <vector>
#include <algorithm>
#include <cmath>

namespace gsc {

inline uint32_t crc32(const uint8_t* d, size_t n)
{
	static uint32_t tab[256];
	static bool init = false;
	if (!init) {
		for (uint32_t i = 0; i < 256; i++) { uint32_t c = i; for (int k = 0; k < 8; k++) c = c & 1 ? 0xedb88320u ^ (c >> 1) : c >> 1; tab[i] = c; }
		init = true;
	}
	uint32_t c = 0xffffffffu;
	for (size_t i = 0; i < n; i++) c = tab[(c ^ d[i]) & 255] ^ (c >> 8);
	return c ^ 0xffffffffu;
}

inline void put_varint(std::vector<uint8_t>& o, uint64_t v) { while (v >= 128) { o.push_back((uint8_t)(v | 128)); v >>= 7; } o.push_back((uint8_t)v); }
inline uint64_t zz(int64_t v) { return ((uint64_t)v << 1) ^ (uint64_t)(v >> 63); }
inline int64_t unzz(uint64_t v) { return (int64_t)(v >> 1) ^ -(int64_t)(v & 1); }
// 1 ok, 0 at the end of the buffer (incomplete), -1 overlong
inline int get_varint(const uint8_t*& p, const uint8_t* end, uint64_t& v)
{
	v = 0;
	for (int s = 0; s < 64; s += 7) { if (p >= end) return 0; const uint8_t c = *p++; v |= (uint64_t)(c & 127) << s; if (!(c & 128)) return 1; }
	return -1;
}

struct BitWriter {
	std::vector<uint8_t>& o; uint64_t acc = 0; int nb = 0;
	explicit BitWriter(std::vector<uint8_t>& out) : o(out) {}
	void put(uint64_t v, int n) {            // the n low bits of v, n <= 56
		acc |= (v & ((n == 64) ? ~0ull : ((1ull << n) - 1))) << nb; nb += n;
		while (nb >= 8) { o.push_back((uint8_t)acc); acc >>= 8; nb -= 8; }
	}
	void ones(uint64_t n) { while (n >= 32) { put(0xffffffffu, 32); n -= 32; } put((1ull << n) - 1, (int)n); }
	void flush() { if (nb) { o.push_back((uint8_t)acc); acc = 0; nb = 0; } }
};
struct BitReader {
	const uint8_t* p; const uint8_t* end; uint64_t acc = 0; int nb = 0; bool bad = false;
	BitReader(const uint8_t* b, const uint8_t* e) : p(b), end(e) {}
	void fill() { while (nb <= 56) { if (p < end) acc |= (uint64_t)*p++ << nb; else if (nb == 0) { bad = true; return; } else return; nb += 8; } }
	uint64_t get(int n) {                    // n <= 56
		if (n == 0) return 0;
		if (nb < n) fill();
		if (nb < n) { bad = true; return 0; }
		const uint64_t v = acc & ((1ull << n) - 1); acc >>= n; nb -= n; return v;
	}
	int bit() { return (int)get(1); }
};

inline void put_gamma(BitWriter& w, uint64_t v)      // v >= 1
{
	int L = 63 - __builtin_clzll(v);
	for (int t = 0; t < L; t++) w.put(0, 1);
	for (int t = L; t >= 0; t--) w.put(v >> t & 1, 1);
}
inline uint64_t get_gamma(BitReader& r)
{
	int L = 0;
	while (!r.bad && r.bit() == 0) if (++L > 63) { r.bad = true; return 0; }
	uint64_t v = 1;
	for (int t = 0; t < L && !r.bad; t++) v = v << 1 | (uint64_t)r.bit();
	return v;
}

// the Rice parameter with the fewest bits for these gaps
inline int best_k(const std::vector<uint64_t>& gaps)
{
	if (gaps.empty()) return 0;
	unsigned __int128 sum = 0; for (uint64_t g : gaps) sum += g;
	const double mean = (double)sum / gaps.size();
	int k0 = mean > 1 ? (int)(log2(mean)) : 0;
	int best = 0; uint64_t bestbits = ~0ull;
	for (int k = std::max(0, k0 - 2); k <= std::min(56, k0 + 2); k++) {
		uint64_t bits = 0;
		for (uint64_t g : gaps) { const uint64_t q = g >> k; bits += k + (q < 32 ? q + 1 : 32 + 2 * (64 - __builtin_clzll(q - 31)) - 1); }
		if (bits < bestbits) { bestbits = bits; best = k; }
	}
	return best;
}

// sort keys ascending: LSD radix sort in passes of 13 bits for large blocks (the survivors of a special-q can number
// millions, and std::sort then took several times as long as the GPU), std::sort for small ones
inline void sort_keys(std::vector<uint64_t>& keys, std::vector<uint64_t>& tmp)
{
	const size_t n = keys.size();
	if (n < 65536) { std::sort(keys.begin(), keys.end()); return; }
	uint64_t mx = 0; for (uint64_t k : keys) mx |= k;
	const int bits = mx ? 64 - __builtin_clzll(mx) : 1;
	tmp.resize(n);
	std::vector<size_t> cnt(1 << 13);
	for (int sh = 0; sh < bits; sh += 13) {
		std::fill(cnt.begin(), cnt.end(), 0);
		for (uint64_t k : keys) cnt[(k >> sh) & 8191]++;
		size_t sum = 0; for (size_t& c : cnt) { const size_t t = c; c = sum; sum += t; }
		for (uint64_t k : keys) tmp[cnt[(k >> sh) & 8191]++] = k;
		keys.swap(tmp);
	}
}

// one block; keys sorted ascending and distinct
inline void encode_block(std::vector<uint8_t>& o, uint64_t q, uint64_t r, int64_t a0, int64_t b0, int64_t a1, int64_t b1,
                         const std::vector<uint64_t>& keys)
{
	const size_t start = o.size();
	put_varint(o, q); put_varint(o, r);
	put_varint(o, zz(a0)); put_varint(o, zz(b0)); put_varint(o, zz(a1)); put_varint(o, zz(b1));
	put_varint(o, keys.size());
	std::vector<uint64_t> gaps(keys.size());
	for (size_t t = 0; t < keys.size(); t++) gaps[t] = t ? keys[t] - keys[t - 1] - 1 : keys[t];
	const int k = best_k(gaps);
	o.push_back((uint8_t)k);
	std::vector<uint8_t> bits;
	bits.reserve(keys.size() * (k + 2) / 8 + 64);
	{
		BitWriter w(bits);
		for (uint64_t g : gaps) {
			const uint64_t qq = g >> k;
			if (qq < 32) { w.ones(qq); w.put(0, 1); } else { w.ones(32); put_gamma(w, qq - 31); }
			w.put(g, k);
		}
		w.flush();
	}
	put_varint(o, bits.size());
	o.insert(o.end(), bits.begin(), bits.end());
	const uint32_t c = crc32(o.data() + start, o.size() - start);
	for (int t = 0; t < 4; t++) o.push_back((uint8_t)(c >> (8 * t)));
}

struct Block { uint64_t q, r; int64_t a0, b0, a1, b1; std::vector<uint64_t> keys; };

// decode the block at p; 0 ok (p advanced past it), 1 clean end of data, 2 the data ends inside the block (a stream
// reader reads more and retries; at the end of a file it is truncated), -1 damaged. p moves only on success.
inline int decode_block(const uint8_t*& p0, const uint8_t* end, Block& B)
{
	const uint8_t* p = p0;
	if (p == end) return 1;
	const uint8_t* start = p;
	uint64_t v[7];
	for (int t = 0; t < 7; t++) { const int rc = get_varint(p, end, v[t]); if (rc <= 0) return rc == 0 ? 2 : -1; }
	B.q = v[0]; B.r = v[1]; B.a0 = unzz(v[2]); B.b0 = unzz(v[3]); B.a1 = unzz(v[4]); B.b1 = unzz(v[5]);
	const uint64_t n = v[6];
	if (p >= end) return 2;
	const int k = *p++;
	if (k > 56) return -1;
	uint64_t nbytes;
	{ const int rc = get_varint(p, end, nbytes); if (rc <= 0) return rc == 0 ? 2 : -1; }
	if (nbytes > (1ull << 30) || n > nbytes * 8) return -1;     // no block is that large; each code has at least one bit
	if (nbytes > (uint64_t)(end - p) || (uint64_t)(end - p) - nbytes < 4) return 2;
	const uint8_t* bits = p; p += nbytes;
	uint32_t c = 0; for (int t = 0; t < 4; t++) c |= (uint32_t)p[t] << (8 * t);
	if (crc32(start, p - start) != c) return -1;
	p += 4;
	B.keys.resize(n);
	BitReader rd(bits, bits + nbytes);
	uint64_t key = 0;
	for (uint64_t t = 0; t < n; t++) {
		uint64_t qq = 0;
		while (qq < 32 && rd.bit() == 1 && !rd.bad) qq++;
		if (qq == 32) qq = 31 + get_gamma(rd);
		const uint64_t g = qq << k | rd.get(k);
		if (rd.bad) return -1;
		key = t ? key + g + 1 : g;
		B.keys[t] = key;
	}
	p0 = p;
	return 0;
}

}  // namespace gsc
