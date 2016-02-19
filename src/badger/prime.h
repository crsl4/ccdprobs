#ifndef PRIMEHDR
#define PRIMEHDR

using namespace std;

class FindPrime {
public:
  static int nextPrime(int s) {
    // Returns the first (likely) prime number >= s.
    if(s%2==0)
      s++;
    while(!isPrime(s))
      s += 2;
    return s;
  }
private:
  static unsigned long long int witness(const unsigned long long int &a, const unsigned long long int &i, 
					const unsigned long long int &n) {
    // Returns 0 if a is definitely a composite number. Otherwise returns a^i % n.
    if(i == 0)
      return 1;
    unsigned long long int x = witness(a,i/2,n);
    if(x == 0)
      return 0;
    unsigned long long int y = (x * x) % n;
    if(y == 1 && x != 1 && x != n-1)
      return 0;
    if(i % 2 != 0)
      y = (a * y) % n;
    return y;
  }

  static int isPrime(unsigned long long int p) {
    // Returns 0 if p is definitely a composite number. Returns 1 if p is probably prime.
    for(int i=0;i<50;i++) {
      unsigned long long int x = rand()%(p-3)+1;
      if(witness(x,p-1,p) != 1) 
	return 0;
    }
    return 1;
  }
};

#endif
