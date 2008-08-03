//---------------------------------------------------------------------------
// a2d.h
// Author: Pary WU
// Defines a2d map multidimensional array class
//---------------------------------------------------------------------------
#ifndef __A2D_HPP_040422163308
#define __A2D_HPP_040422163308

#ifndef __cplusplus
#error "For C++ only"
#endif

#include <map>

#define __A2D_DUMP__
#ifdef __A2D_DUMP__
#include <iostream>
#endif
#include "oberror.h"
using namespace mesmer;
using namespace std;

template <class V, class K=int>
class a2d_t{
public:
  typedef std::map< K, V > leaf_t; //the leaf (container) of the whole structure
  typedef std::map< K, void* > proxy_map_t; //the branch

  template <class T>
  struct proxy_t: public proxy_map_t{ //this proxy_t is based on the proxy_map_t
    typedef typename proxy_map_t::iterator iterator; //it has iterator the same type as the proxy_map_t
    proxy_t(): proxy_map_t(){} //copy this proxy_t is just the same as copy a proxy_map_t
    void remove(){	//in order remove this proxy_t one has to follow the iterator
      for(iterator i = this->begin(); i != this->end(); i++)
        delete (T*)(i->second);    
    }
    ~proxy_t(){
      remove();
    }

    T& operator[](const K i){
      proxy_map_t& m = *this;
      if(m.find(i) == m.end())
        m[i] = new T();
      return *((T*)m[i]);
    }
    T* find_cast(const K key){
      iterator i = find(key);
      if(i != this->end())
        return (T*)(i->second);
      return 0;
    }
  };

  typedef proxy_t<leaf_t> proxy2_t;
  proxy2_t* root;



public:

  void print(const K rows, const K columns){
    ctest << "{" << endl;
    for(int i(0); i<rows; ++i){
      for(int j(0); j<columns; ++j){
        ctest << setw(20) << (*this)[i][j];
      }
      ctest << endl;
    }
    ctest << "}" << endl;
  }

  a2d_t(const int x_size, const int y_size):root(new proxy2_t()){
    /* here x_size, y_size are dummy */
  }
  a2d_t():root(new proxy2_t()){}
  ~a2d_t(){
    delete root;
  }
  leaf_t& operator[](const K i){
    return (*root)[i]; 
  }
  bool exists(const K x, const K y){
    leaf_t* leaf = root->find_cast(x);
    if(leaf)
      if(leaf->find(y) != leaf->end())
        return true;
    return false;
  }
  bool remove(const K x, const K y){
    leaf_t* leaf = root->find_cast(x);
    if(leaf){
      typename leaf_t::iterator j = leaf->find(y);
      if(j != leaf->end()){
        leaf->erase(j);
        if(leaf->size() == 0)
          root->erase(x);
        return true;
      }
    }
    return false;
  }
#ifdef __A2D_DUMP__
  void dump(){
    cout << "{" << endl;
    for(typename proxy2_t::iterator i = root->begin(); i != root->end(); i++){
      leaf_t* l = (leaf_t*)(i->second);
      for(typename leaf_t::iterator j = l->begin(); j!=  l->end(); j++){
        cout << i->first << "," << j->first << ": " << j->second << endl;
      }

    }
    cout << "}" << endl;
  }
#else
  void dump(){}
#endif
};

#endif
