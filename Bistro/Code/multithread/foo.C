#include <iostream>
#include <iomanip>
#include <thread>
#include <vector>
#include <random>

using namespace std;

void thread_function(mt19937_64& rng,int id,int n,vector<int>& v)
{
  cout << "Beginning thread " << id << ", also known as " << this_thread::get_id() << endl;
  for ( int i=0; i<n; ++i )
    {
      uniform_int_distribution<> rint(0,99);
      v[n*id + i] = rint(rng);
    }
  cout << "Ending thread " << id << endl;
 }
      
  int main()
{
  unsigned int cores = thread::hardware_concurrency();
  cout << cores << " cores." << endl;

  random_device rd;
  unsigned int initial_seed = rd();
//  unsigned int initial_seed = 747705983;
  cout << "Seed = " << initial_seed << endl;
  minstd_rand seed_rng(initial_seed);
  uniform_int_distribution<> rint(0,4294967295);
  vector<unsigned int> seeds(cores);
  vector<mt19937_64*> rng;
  for ( vector<unsigned int>::iterator p=seeds.begin(); p!=seeds.end(); ++p )
  {
    *p = rint(seed_rng);
    rng.push_back( new mt19937_64(*p) );
  }
  for ( vector<unsigned int>::iterator p=seeds.begin(); p!=seeds.end(); ++p )
    cout << setw(15) << *p << endl;

  int n = 10;
  vector<int> v(cores*n);
  vector<thread> threads;
  for ( int i=0; i<cores; ++i )
  {
    threads.push_back(thread(thread_function,ref(*(rng[i])),i,n,ref(v)));
  }
  for ( int i=0; i<cores; ++i )
    threads.at(i).join();
  
  cout << "v:";
  for ( vector<int>::iterator p=v.begin(); p!=v.end(); ++p )
    cout << " " << dec << *p;
  cout << endl;
  return 0;
}


