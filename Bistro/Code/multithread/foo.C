#include <iostream>
#include <iomanip>
#include <thread>
#include <vector>

using namespace std;

void thread_function(int id,int n,vector<int>& v)
{
  cout << "Beginning thread " << id << ", also known as " << this_thread::get_id() << endl;
  for ( int i=0; i<n; ++i )
  {
    v[n*id + i] = id;
  }
  cout << "Ending thread " << id << endl;
}

int main()
{
  unsigned int cores = thread::hardware_concurrency();
  cout << cores << " cores." << endl;

  int n = 10;
  vector<int> v(cores*n);
  vector<thread> threads;
  for ( int i=0; i<cores; ++i )
  {
    threads.push_back(thread(thread_function,i,n,ref(v)));
  }
  for ( int i=0; i<cores; ++i )
    threads.at(i).join();
  
  cout << "v:";
  for ( vector<int>::iterator p=v.begin(); p!=v.end(); ++p )
    cout << " " << dec << *p;
  cout << endl;
  return 0;
}


