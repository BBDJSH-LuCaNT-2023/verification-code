static long kronecker(long D,unsigned long p) {
#define mulredc(x,y) (t.d=(__uint128_t)x*y,t.d+=(t.l[0]*pbar)*(__uint128_t)p,\
t.l[1]-=p,t.l[1]+(p&((long)t.l[1]>>63)))
	static long chi8[] = {0,1,0,-1,0,-1,0,1};
	unsigned long q,x,y,pbar;
	int k;
	union { __uint128_t d; unsigned long l[2]; } t;

	if (p == 2) return chi8[D&7];
	if ( !(x=D%(long)p) ) return 0;
	if ((long)x < 0) x += p;

	// compute pbar = -p^{-1} mod 2^64 by Newton's method
	k = ((p+1)>>2<<3)-p; // correct mod 2^4
	k *= 2+k*(int)p; k *= 2+k*(int)p; k *= 2+k*(int)p;
	pbar = (unsigned long)k*(2+(unsigned long)k*p);

	q = p>>1;
	while (!(q & 1)) {
		x = mulredc(x,x);
		q>>=1;
	}

	y=x, q>>=1;
	while (q) {
		x = mulredc(x,x);
		if (q & 1) y = mulredc(y,x);
		q>>=1;
	}

	return (-y) % p ? -1 : 1;
}
