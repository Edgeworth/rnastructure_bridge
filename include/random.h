#if !defined(RANDOM_H)
#define RANDOM_H
#define TABSIZE 32
#include <limits>
#include <vector>

using namespace std;

typedef unsigned long long UINT64; // replace with uint64_t  when c++11 is used

class randomnumber {
 private:
  long IMM1,DIVISOR,/*RNMX,*/idum2,idum,iy,iv[TABSIZE];
  double IM1INV;

 public:
  /*	make a random number generator	*/
  randomnumber();
  /*	seed the random number generator	*/
  void seed(long seeddouble);
  /*	generate a random number between 0.0 and 1.0	*/
  double roll();
  /* generate a random integer in a range
     Seems to work correctly BUT I make no
     guarantees about rigorously correct behavior
     MS 11/2014 */
  int roll_int(int min, int max);
};

// A Pseudo-Random Number Generator that is extremely fast on modern computer architectures
// and has very good random-number qualities (distribution, period etc)
// This is based on xorshift* generator as suggested by Marsaglia.
// It is a 64-bit generator with 64 bits of state has a maximal period of 2^64−1 and passes all the "Diehard tests" and fails only the MatrixRank test of BigCrush.
// See: https://en.wikipedia.org/wiki/Xorshift
// Author: Richard M. Watson 2017
class rand64 {
 private:
  UINT64 state;
  static UINT64 calc_next(UINT64 &state);
 public:
  //!	make a random number generator with a seed value. Specify 0 for the seed to seed with time(0)
  rand64(UINT64 seed);
  //!	seed the random number generator.  Specify 0 for the seed to seed with time(0)
  void seed(UINT64 seed);

  //!	generate a random number in the range [0, 1.0)
  double nextDouble();

  //!	generate a random int in the range [-INT_MAX, INT_MAX)
  int nextInt();
  //!	generate a random int in the range [0, max)
  int nextInt(const int max);
  //!	generate a random int in the range [min, max)  i.e.  min <= value < max
  int nextInt(const int min, const int max);

  //! return the next random 64-bit unsigned integer
  UINT64 next();
  //! return the current random 64-bit unsigned integer (i.e. the last one generated by a call to next*)
  UINT64 current() const;

  // function operator returns the same as next()
  UINT64 operator()();
  // function operator returns the same as nextInt(int max)
  int operator()(const int max);
  // function operator returns a random element from the vector
  template<class T> T &operator()(vector<T> &list) { return list[(size_t)(next() % list.size())]; }
  template<class T> const T &operator()(const vector<T> &list) { return list[(size_t)(next() % list.size())]; }
};


class RandomOutput {
 public:

  //! A template used to generate a double functor() from a PRNG like rand64.
  //! example:
  //!
  //! rand64 r(1234);
  //! randomDouble<rand64> d(r);
  //! cout << d() << endl; // outputs a random double
  template<typename T> class Double {
   private:
    T _generator;
   public:
    inline Double(const T generator) : _generator(generator) { }
    inline double operator()() { return _generator.nextDouble(); }
  };

  //! A template used to generate a int functor() from a PRNG like rand64.
  //! The functor then outputs a random int within a specified range.
  //! example:
  //!
  //! rand64 r(1234);
  //! randomInt<rand64> i(r, 10, 20);
  //! cout << i() << endl; // outputs a random int in the range [10, 20)
  template<typename T> class Int {
   private:
    T rnd;
    int minVal, maxVal;
   public:
#undef min
#undef max
    inline Int(const T generator, const int min=numeric_limits<int>::min(), const int max=numeric_limits<int>::max()) :
                                                                                                                            rnd(generator),
                                                                                                                            minVal(min),
                                                                                                                            maxVal(max) { }
    inline int operator()() { return rnd.nextInt(minVal, maxVal); }
  };

  template<typename T> inline static Double<T> asDouble(const T generator) { return Double<T>(generator); }
  template<typename T> inline static Int<T> asInt(const T generator) { return Int<T>(generator); }
};

#endif
