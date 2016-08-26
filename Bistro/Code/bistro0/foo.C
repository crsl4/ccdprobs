#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include <string>

using namespace std;

int main()
{
  multimap<string,double> m;
  m.insert(pair<string,double>("A",0.1));
  m.insert(pair<string,double>("A",0.2));
  m.insert(pair<string,double>("B",0.3));

  cout << "Before:" << endl;
  for ( multimap<string,double>::iterator p=m.begin(); p!=m.end(); ++p )
  {
    cout << p->first << " --> " << p->second << endl;
  }
  cout << endl;
  for ( multimap<string,double>::iterator p=m.begin(); p!=m.end(); ++p )
  {
    p->second -= 0.1;
  }
  cout << "After:" << endl;
  for ( multimap<string,double>::iterator p=m.begin(); p!=m.end(); ++p )
  {
    cout << p->first << " --> " << p->second << endl;
  }
  
  return 0;
}
