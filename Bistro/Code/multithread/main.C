#include <iostream>
#include <thread>

using namespace std;

// creating a thread using a Function pointer
// void thread_function()
// {
//   //cout << "entering thread function" << endl;
//   for ( int i=0; i < 10; i++)
//     cout << "thread function executing" << endl;
// }

// int main()
// {
//   //cout << "entering main" << endl;
//   thread threadObj(thread_function);
//   //cout << "in main, after calling thread function" << endl;
//   for( int i=0; i < 10; i++)
//     cout << "Display from main thread" << endl;
//   threadObj.join();
//   cout << "exit from main function" << endl;
//   return 0;
// }

// creating a thread using function objects
// class DisplayThread
// {
// public:
//   void operator()()
//   {
//     for( int i = 0; i < 100; i++)
//       cout << "Display thread executing " << i << endl;
//   }
// };

// int main()
// {
//   cout << "entering main" << endl;
//   thread threadObj( (DisplayThread()) );
//   for(int i=0; i < 10; i++)
//     cout << "Display from the main thread " << i << endl;
//   cout << "waiting for thread to complete" << endl;
//   threadObj.join();
//   cout << "Exiting from main thread" << endl;
//   return 0;
// }

// creating a thread using lambda functions
// int main()
// {
//   //int x = 9;
//   //cout << "entering main" << endl;
//   thread threadObj([]{
//       for( int i=0; i<10; i++)
// 	cout << "Display thread executing " << i << endl;
//     });
  
//   for(int i=0; i < 10; i++)
//     cout << "Display from main thread " << i << endl;

//   threadObj.join();
//   cout << "Exiting from main" << endl;
//   return 0;
// }


// thread id
// void thread_function()
// {
//   cout << "inside thread :: ID = " << this_thread::get_id() << endl;
// }

// int main()
// {
//   thread threadObj1(thread_function);
//   thread threadObj2(thread_function);

//   if(threadObj1.get_id() != threadObj2.get_id())
//     cout << "Both threads have different IDs" << endl;

//   cout << "from main :: ID of thread 1: " << threadObj1.get_id() << endl;
//   cout << "from main :: ID of thread 2: " << threadObj2.get_id() << endl;

//   threadObj1.join();
//   threadObj2.join();
//   return 0;
// }


// passing arguments to a thread
void threadCallback(int x, string str)
{
  cout << "entering thread callback" << endl;
  cout << "Passed number = " << x << endl;
  cout << "Passed string = " << str << endl;
}

int main()
{
  //may return 0 when not able to detect
  unsigned concurrentThreadsSupported = std::thread::hardware_concurrency();
  cout << "hardware concurrency: " << concurrentThreadsSupported << endl;
  int x = 10;
  string str = "claudia";
  cout << "main before calling thread" << endl;
  thread threadObj(threadCallback,x,str);
  cout << "main after calling thread" << endl;
  threadObj.join();
  return 0;
}
