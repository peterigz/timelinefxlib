//The constants and hashes the GPU modules share, each the twin of the CPU function it names

//tfxPI2
static const float tfx_two_pi = 6.283185307;

//tfx__seedgen_u32, a murmur3 finalizer. tfx__wide_seedgen_base(base, h) is tfx_seedgen(base ^ h)
uint tfx_seedgen(uint hash) {
	hash ^= hash >> 15;
	hash *= 0x85ebca6bu;
	hash ^= hash >> 13;
	hash *= 0xc2b2ae35u;
	hash ^= hash >> 16;
	return hash;
}

//tfx__white_unit and tfx__wide_white_unit: the signed reinterpret is what centres it on 0
float tfx_white_unit(uint seed, uint axis) {
	return (float(int(tfx_seedgen(seed ^ axis))) / 4294967295.0) * 2.0;
}

//tfx__hash_unit, exact on both sides because 24 bits convert to a float without rounding
float tfx_hash_unit(uint seed, uint axis) {
	return float(tfx_seedgen(seed ^ axis) >> 8) * (1.0 / 16777216.0);
}
