/* 
   LightMutex.h
   Copied from the Open Shogi Lib (OSL).
   See http://gps.tanaka.ecc.u-tokyo.ac.jp/gpsshogi/pukiwiki.php?OpenShogiLib
   See their license at http://gps.tanaka.ecc.u-tokyo.ac.jp/gpsshogi/index.php?osl%2F%A5%E9%A5%A4%A5%BB%A5%F3%A5%B9
*/

#ifndef LIGHT_MUTEX_H
#define LIGHT_MUTEX_H

#include <boost/thread/thread.hpp>
#include <boost/utility.hpp>

class LightMutex : boost::noncopyable {
  volatile int lock;
public:
  class scoped_lock;
  LightMutex() :lock(0) {}
  bool tryLock(){
    int dummy;
    asm __volatile__(" movl $1,%0" "\n\t"
                     " xchgl (%1),%0" "\n\t"
                     : "=&q"(dummy)
                     : "q"(&lock)
                     : "cc");
    return dummy==0;
  }
  bool waitLock(int counter){
    for(int i=0;i<counter;i++){
      asm __volatile__(" pause" "\n\t");
      if(lock==0)
        return true;
    }
    return false;
  }
  void unLock(){
    lock=0;
  }
};

class LightMutex::scoped_lock : boost::noncopyable {
  LightMutex& m;
public:
  scoped_lock(LightMutex& m) :m(m){
    while(!m.tryLock()){
      for(int i=0;i<2;i++){
        if(!m.waitLock(100)) break;
        if(m.tryLock()) return;
      }
      boost::thread::yield();
    }
  }
  ~scoped_lock(){
    m.unLock();
  }
};

#endif /* LIGHT_MUTEX_H */
// ;;; Local Variables:
// ;;; mode:c++
// ;;; c-basic-offset:2
// ;;; End:

