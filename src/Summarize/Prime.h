#ifndef PRIMEHDR
#define PRIMEHDR

using namespace std;

class Prime {
public:
  static int nthprime(int n);

  static int nextPrime(int s);

private:
  static unsigned long long int witness(const unsigned long long int &a, const unsigned long long int &i, 
					const unsigned long long int &n);

  static int isPrime(unsigned long long int p);
};

#endif
