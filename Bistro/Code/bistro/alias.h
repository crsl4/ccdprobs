#ifndef ALIAS_H_
#define ALIAS_H_

#include <vector>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <map>
#include <random>

using namespace std;

// Program uses the BOOST dynamic bitset library in order to
// store clades more efficiently and without a limit on the size of the tree
#include <boost/dynamic_bitset.hpp>

using namespace boost;

template <typename T>
struct Cell {
    T alias;
    T index;
    double cutoff;
};

// class alias uses the alias algorithm to randomly select a random dynamic_bitset from a given multinomial distribution
// class alias is the same, but the objects selected are integers
template <typename T>
class Alias {
private:
  vector<Cell<T> > cells;
  void computeTwoPointDistributions();
public:
  Alias() {}
  Alias(vector<double>&, vector<T>&);
  virtual ~Alias() {}
  void initiate(vector<double>&, vector<T>&);
  T pick(mt19937_64&);
  void printTwoPointDistributions();
};

template <typename T>
Alias<T>::Alias(vector<double>& ps, vector<T>& is)
{
  initiate(ps,is);
}

template <typename T>
void Alias<T>::initiate(vector<double>& ps, vector<T>& is)
{
  int index = 0;
  int nCells = ps.size();
  cells.resize(nCells);
  for(vector<double>::iterator pit = ps.begin(); pit != ps.end(); ++pit, ++index) {
    cells[index].cutoff = *pit * nCells;
    cells[index].alias = is[index];
    cells[index].index = is[index];
  }
  computeTwoPointDistributions();
}

template <typename T>
T Alias<T>::pick(mt19937_64& rng)
{
  uniform_int_distribution<> rint(0,cells.size()-1);
  int k = rint(rng);
  uniform_real_distribution<> runif(0,1);
  return ( runif(rng) < cells[k].cutoff ? cells[k].index : cells[k].alias );
}

template <typename T>
void Alias<T>::computeTwoPointDistributions()
{
  typename vector< Cell<T> >::iterator litr = cells.end() - 1;
  typename vector< Cell<T> >::iterator ritr = cells.begin();
   
    // sort the vector of multiplied probabilities such that all cells
    // having probability less than 1 are at the end, cells with probabilities
    // greater than 1 are at the beginning
    // e.g., if initial cells are 0.8, 0.9, 1.1, 1.2
    // after this while loop, cells become 1.1, 1.2, 0.8, 0.9
    while (true) {
        while (ritr != cells.end() && ritr->cutoff >= 1.00000) {
            ritr++;
        }
        if (ritr == cells.end()) {
            break;
        }

        while (litr != cells.begin() && litr->cutoff <= 1.00000) {
            litr--;
        }

        if (litr == cells.begin() && litr->cutoff <= 1.00000) {
            break;
        }

        if (litr < ritr)
            break;

        Cell<T> temp = *ritr;
        *ritr = *litr;
        *litr = temp;
        litr--; ritr++;
    }

    litr = cells.end() - 1;
    ritr = cells.end() - 1;
    while (true) {
        //find a cell with probability initial cutoff greater than 1.0
        while (ritr != cells.begin() && ritr->cutoff <= 1.00000) {
            ritr--;
        }
        if (ritr == cells.begin() && ritr->cutoff <= 1.00000) {
            break;
        }

        //find a cell with probability initial cutoff less than 1.0
        while (litr != cells.begin() && litr->cutoff >= 1.00000) {
            litr--;
        }

        if (litr == cells.begin() && litr->cutoff >= 1.00000) {
            break;
        }

        litr->alias = ritr->index;
        ritr->cutoff -= (1 - litr->cutoff);
        litr--;
    }
}

template <typename T>
void Alias<T>::printTwoPointDistributions()
{
  cerr << "Two point distribution\n";
  cerr << "Cutoff:\t";
  for (int i = 0; i < cells.size(); i++) {
    cerr << cells[i].cutoff << "\t";
  }
  cerr << endl;
  cerr << "Index:\t";
  for (int i = 0; i < cells.size(); i++) {
    cerr << cells[i].index << "\t";
  }
  cerr << endl;
  cerr << "Alias:\t";
  for (int i = 0; i < cells.size(); i++) {
    cerr << cells[i].alias << "\t";
  }
  cerr << endl;
}

#endif /* ALIAS_H_ */
