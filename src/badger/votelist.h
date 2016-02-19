#ifndef VOTELISTHDR
#define VOTELISTHDR

#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <string>
#include "prime.h"

using namespace std;

template <class K, class D> class Database {
  // A simple database implemented by hashing with linear probing and deletion compaction.
  // K must have a hash() method.

public:

  Database(int s, const K& e) : size(s), empty(e), keys(s,e), data(s) {}

  int find(const K &key, D &datum, int &m) {
    // Find key and set datum to its data. Returns 1 if found, 0 otherwise.
    // m is set to the position where it is, or the position where it should go.
    for(m=key.hash()%size;keys[m]!=empty;m=(m+1)%size) {
      if(keys[m]==key) {
	datum = data[m];
	return 1;
      }
    }
    return 0;
  }
  
  void add(const K &key, const D &datum, int m) {
    // Add (key,data). key is not already in the database. m is its position.
    keys[m] = key;
    data[m] = datum;
  }

  void remove(const K &key) {
    // Remove (key,-). key must be in the database.
    int orig1 = key.hash()%size, pos1;
    for(pos1=orig1;keys[pos1]!=key;pos1=(pos1+1)%size)
      {}
    remove2(orig1,pos1);
  }

private:
  int size;		// Size of the hash table.
  K empty;		// Key used to indicate empty slot.
  vector<K> keys;	// The keys.
  vector<D> data;	// The data.

  void remove2(int orig1, int pos1) {
    // replace keys[pos1] by a key that would be hashed there, or empty if no such keys.
    int orig2,pos2;
    for(pos2=(pos1+1)%size;keys[pos2]!=empty;pos2=(pos2+1)%size) {
      orig2 = keys[pos2].hash()%size;
      if(orig1<pos2 ? orig2<=pos1 : (pos2<pos1 ? (pos2<orig2 && orig2<=pos1) : (orig2<pos1 || pos2<orig2))) {
	keys[pos1] = keys[pos2];
	data[pos1] = data[pos2];
	remove2(orig2,pos2);
	return;
      }
    }
    keys[pos1]=empty;
  }
};

template <class K> class Contender {
 public:
  Contender() {}
  Contender(K i, int v, int n) : item(i), votes(v), id(n) {}
  K item;
  int votes;
  int id;
};

// Based on MJRTY - A fast majorty vote algorithm by Robert S. Boyer and J Strother Moore. 

template <class K> class VoteList {
  // Find top m vote-getters in an election.

 private:
  class Tally;
  class Candidate;

  typedef typename VoteList<K>::Tally TallyK;
  typedef typename VoteList<K>::Candidate CandidateK;
  typedef typename vector<TallyK>::iterator TallyPtr;
  typedef typename vector<CandidateK>::iterator CandidatePtr;

 public:

  VoteList(int m, const K& empty) :
    // Initialize everything. 
    // m is the number of vote-getters desired and should be < 20000.
    // empty is a special instance of K that is not valid.
    // It is needed to indicate an empty slot in the hash table.
    // K must have a method hash().
    
    numCandidates(0),
    maxCandidates(m),
    tallyList(m),
    candidateList(m),
    firstTally(tallyList.end()),
    endTally(tallyList.end()),
    tallyFirstFree(tallyList.end()),
    talliesUsed(0),
    candidateFirstFree(candidateList.end()),
    endCand(candidateList.end()),
    beginCand(candidateList.begin()),
    candidatesUsed(0),
    db(FindPrime::nextPrime(3*m+1),empty)
    {}

  void voteForCandidate(K vote, int n) {
    // Called on first pass to find candidates.
    // Add n to the number of votes for the candidate represented by vote.
    // May add candidate to the tally list or remove lower vote-getting candidates from the tally list.
    //checkTallies();
    CandidatePtr cand;
    int m;
    if(db.find(vote,cand,m))
      addVotes(cand,n);
    else {
      if(numCandidates>=maxCandidates)
	subtractVotes(n);
      if(n>0) {
	cand = newCandidate(vote);
	db.add(vote,cand,m);
	addCandidate(cand,n);
      }
    }
  }

  void checkCandidate(K vote, int n) {
    // Called on second pass to verify winners.
    CandidatePtr cand;
    int m;
    if(db.find(vote,cand,m))
      cand->contender.votes += n;
  }

  int getContender(K vote, Contender<K>& contender) {
    // Return the id of a candidate.
    CandidatePtr cand;
    int m;
    if(db.find(vote,cand,m)) {
      contender = cand->contender;
      return 1;
    }
    else 
      return 0;
  }

  class iterator {
  public:
    iterator();

    iterator(TallyPtr t, CandidatePtr c, TallyPtr et, CandidatePtr ec, CandidatePtr bc) :
      tally(t), cand(c), endTally(et), endCand(ec), beginCand(bc) {}

    iterator(const iterator& x) : 
      tally(x.tally), cand(x.cand), endTally(x.endTally), endCand(x.endCand), beginCand(x.beginCand) {}

    iterator& operator=(const iterator& x) { 
      tally=x.tally; cand=x.cand; endTally=x.endTally; endCand=x.endCand; beginCand=x.beginCand; 
      return *this; 
    }

    bool operator==(const iterator& x) { 
      return (tally==x.tally && cand==x.cand && endTally==x.endTally && 
	      endCand==x.endCand && beginCand==x.beginCand);
    }

    bool operator!=(const iterator& x) { 
      return (tally!=x.tally || cand!=x.cand || endTally!=x.endTally || 
	      endCand!=x.endCand || beginCand!=x.beginCand);
    }

    const Contender<K>& operator*() const { return cand->contender; }

    const Contender<K>* operator->() const { return &(cand->contender); }

    iterator& operator++() { 
      if(cand!=endCand)
	if(cand->next!=endCand)
	  cand = cand->next;
	else {
	  tally = tally->next;
	  if(tally!=endTally) 
	    cand = tally->candidates;
	  else
	    cand = endCand;
	}
      return *this;
    }

    iterator operator++(int) {
      iterator temp(*this);
      operator++();
      return temp;
    }
  private:
    TallyPtr tally,endTally;
    CandidatePtr cand,endCand,beginCand;
  };

  iterator begin() const { 
    TallyPtr tally = firstTally;
    CandidatePtr cand = (firstTally==endTally ? endCand : firstTally->candidates);
    return iterator(tally,cand,endTally,endCand,beginCand);
  }

  iterator end() const { return iterator(endTally,endCand,endTally,endCand,beginCand); }

 private:
  class Candidate { 
    // The candidate record. All candidates with the same number of votes are stored
    // as an unordered doubly-linked list attached to the tally record for that number of votes.
  public:
    Contender<K> contender; // The candidate information: the individual (K), votes, and an id.
    TallyPtr tally;	    // Pointer to the tally record.
    CandidatePtr prev;	    // Previous candidate record in the tally's list of candidates.
    CandidatePtr next;	    // Next candidate record in the tally's list of candidates.
  };

  class Tally {
    // List of all candidates with the same number of votes. Tallies are stored in a doubly
    // linked list ordered by number of votes.
  public:
    int count;			// Number of votes - number of votes in previous tally record 
    int numCands;		// Number of candidates in tally record's list.
    CandidatePtr candidates;	// Head of list of candidates.
    TallyPtr prev;		// Previous tally record.
    TallyPtr next;		// Next tally record.
  };

  int numCandidates;		// Number of current candidates.
  int maxCandidates;		// Maximum number of candidates.
  vector<Tally> tallyList;	// Freelist of tally records.
  vector<Candidate> candidateList;// Freelist of candidate records.
  TallyPtr firstTally;		// Head of tally list.
  TallyPtr endTally;		// End of the tally list.
  TallyPtr tallyFirstFree;	// Head of tally freelist.
  int talliesUsed;		// Number of tally freelist elements in use.
  CandidatePtr candidateFirstFree;// Head of candidate freelist.
  CandidatePtr endCand;		// End of the candidate list.
  CandidatePtr beginCand;	// Beginning of the candidate list.
  int candidatesUsed;		// Number of candidate freelist elements in use.
  Database<K,CandidatePtr> db;	// Database used to find candidate records.

  CandidatePtr newCandidate(K n) {
    // Get a new candidate record from the freelist.
    CandidatePtr cand = candidateFirstFree;
    if(cand==endCand) 
      cand = candidateList.begin() + candidatesUsed++;
    else
      candidateFirstFree = cand->next;
    cand->contender.item = n;
    cand->contender.votes = 0;
    cand->contender.id = cand - beginCand;
    return cand;
  }

  void deleteCandidate(CandidatePtr cand) {
    // Return a candidate record to the freelist.
    cand->next = candidateFirstFree;
    candidateFirstFree = cand;
  }

  TallyPtr newTally(int c, CandidatePtr cand, TallyPtr p, TallyPtr nx) {
    // Get a new tally record from the freelist.
    TallyPtr tally = tallyFirstFree;
    if(tally==endTally) 
      tally = tallyList.begin() + talliesUsed++;
    else
      tallyFirstFree = tally->next;
    cand->prev = cand->next = endCand;
    cand->tally = tally;

    tally->count = c;
    tally->numCands = 1;
    tally->candidates = cand;
    tally->prev = p;
    tally->next = nx;
    return tally;
  }

  void deleteTally(TallyPtr tally) {
    // Return a tally record to the freelist.
    tally->next = tallyFirstFree;
    tallyFirstFree = tally;
  }

  void addCandidateToTally(CandidatePtr cand, TallyPtr tally) {
    cand->prev = endCand;
    cand->next = tally->candidates;
    cand->tally = tally;
    tally->candidates->prev = cand;
    tally->candidates = cand;
    tally->numCands++;
  }

  TallyPtr findTally(TallyPtr tally, TallyPtr& prev, int& n) {
    // Finds the first tally with >= n cumulative votes. tally must not be endTally.
    // n will be decremented by the intervening tallys' votes.
    // prev should be the tally's prev.
    TallyPtr p;
    for(p=tally;p!=endTally && n>p->count;prev=p,p=p->next)
      n -= p->count;
    return p;
  }

  void addBefore(TallyPtr tally, TallyPtr prev, TallyPtr next) {
    // Insert tally between prev and next in the doubly-linked list.
    tally->next = next;
    tally->prev = prev;
    addOtherLinks(tally,prev,next);
  }

  void addOtherLinks(TallyPtr tally, TallyPtr prev, TallyPtr next) {
    // Insert tally between prev and next. Assumes tally's links have been set.
    if(prev==endTally)
      firstTally = tally;
    else
      prev->next = tally;
    if(next!=endTally)
      next->prev = tally;
  }

  void addVotes(CandidatePtr cand, int n) {
    // Add n votes to candidate cand.
    TallyPtr tally=cand->tally;
    TallyPtr prev;
    if(tally->numCands==1) { // Only candidate in the tally record.
      if(tally->next==endTally)
	tally->count += n;
      else if(tally->next->count>n) {
	tally->count += n;
	tally->next->count -= n;
      }
      else {
	tally->next->count += tally->count;

	// Delink tally.
	if(tally->prev!=endTally)
	  tally->prev->next = tally->next;
	else
	  firstTally =  tally->next;
	tally->next->prev = tally->prev;
	prev = tally->prev;
	TallyPtr tally2 = findTally(tally->next,prev,n);
	if(tally2!=endTally && tally2->count==n) { // Add candidate to tally2's list
	  addCandidateToTally(cand,tally2);
	  deleteTally(tally);
	}
	else { // Link tally before tally2
	  tally->count = n;
	  if(tally2!=endTally)
	    tally2->count -= n;
	  addBefore(tally,prev,tally2);
	}
      }
    }
    else { // More than one candidate in the tally record
      // Delink cand from list.
      if(cand->prev!=endCand)
	cand->prev->next = cand->next;
      else
	tally->candidates = cand->next;
      if(cand->next!=endCand)
	cand->next->prev = cand->prev;
      tally->numCands--;
      // Either create new tally record or add cand to appropriate tally.
      prev = tally;
      TallyPtr tally2 = findTally(tally->next,prev,n);
      if(tally2!=endTally && tally2->count==n)
	addCandidateToTally(cand,tally2);
      else {
	TallyPtr tally3 = newTally(n,cand,prev,tally2);
	if(tally2!=endTally)
	  tally2->count -= n;
	addOtherLinks(tally3,prev,tally2);
      }
    }
  }

  void addCandidate(CandidatePtr cand, int n) {
    // Add new candidate.
    numCandidates++;
    if(firstTally == endTally)
      firstTally = newTally(n,cand,endTally,endTally);
    else {
      TallyPtr prev = endTally;
      TallyPtr tally = findTally(firstTally,prev,n);
      if(tally!=endTally && tally->count==n)
	addCandidateToTally(cand,tally);
      else {
	TallyPtr tally2 = newTally(n,cand,prev,tally);
	if(tally!=endTally)
	  tally->count -= n;
	addOtherLinks(tally2,prev,tally);
      }
    }
  }

  void subtractVotes(int& n) {
    // Subtract n votes from first tally record, removing it if vote count becomes 0.
    // firstTally should not be endTally; n is decremented by the number of deleted votes.
    if(firstTally!=endTally)
      if(firstTally->count>n) {
	firstTally->count -= n;
	n = 0;
      }
      else {
	n -= firstTally->count;
	numCandidates -= firstTally->numCands;
	// Remove the candidate records first.
	for(CandidatePtr c=firstTally->candidates,m; c!=endCand; c=m) {
	  m = c->next;
	  db.remove(c->contender.item);
	  deleteCandidate(c);
	}
	// Remove tally record.
	TallyPtr temp=firstTally;
	firstTally = firstTally->next;
	if(firstTally!=endTally)
	  firstTally->prev = endTally;
	deleteTally(temp);
	subtractVotes(n);
      }
  }

  void errorMsg(string msg, int x) { cerr << msg << " (" << x << ")" << endl; exit(1); }

  void printTallies(ostream& c) {
    for(TallyPtr p=firstTally;p!=endTally;p=p->next) {
      c << "(" << p-tallyList.begin() << "," << p->count << ") ";
      for(CandidatePtr q=p->candidates;q!=endCand;q=q->next) {
	q->contender.item.print(c);
	c << " ";
      }
    }
    c << endl;
  }

  void checkTallies() {
    for(TallyPtr p=firstTally;p!=endTally;p=p->next) {
      if(!(tallyList.begin()<=p && p<tallyList.end())) 
	errorMsg("Error in tally list, bad linked list element", p-tallyList.begin());
      int n=0;
      for(CandidatePtr q=p->candidates;q!=endCand;q=q->next) {
	n++;
	if(!(candidateList.begin()<=q && q<candidateList.end())) 
	  errorMsg("Error in candidate list, bad linked list element", p-tallyList.begin());
	if(q->next!=endCand && q->next->prev!=q)
	  errorMsg("Error in candidate list, bad prev pointer", q-candidateList.begin());
	if(!(q==p->candidates?q->prev==endCand:q->prev->next==q))
	  errorMsg("Error in candidate list, bad next pointer", q-candidateList.begin());
	if(q->tally!=p)
	  errorMsg("Error in candidate list, bad tally pointer", q-candidateList.begin());
	if(q->prev==q->next && q->prev!=endCand)
	  errorMsg("Error in candidate list, prev==next", q-candidateList.begin());
      }
      if(p->next!=endTally && p->next->prev!=p)
	errorMsg("Error in tally list, bad prev pointer", p-tallyList.begin());
      if(!(p==firstTally?p->prev==endTally:p->prev->next==p))
	errorMsg("Error in tally list, bad next pointer", p-tallyList.begin());
      if(p->count<=0)
	errorMsg("Error in tally list, bad count", p-tallyList.begin());
      if(n!=p->numCands)
	errorMsg("Error in tally list, bad numCands", p-tallyList.begin());
    }
  }
};

#endif
