// ref_vector.h

#ifndef REF_VECTOR_H
#define REF_VECTOR_H

#include <TRefArray.h>

namespace g2esoft{

  template<typename T> class ref_vector_iterator;
/*  
  template<typename T> class ref_vector : public TObject {
  public:
    ref_vector(){}
    ~ref_vector(){}
    
    typedef ref_vector_iterator<T> iterator;
    typedef ref_vector_iterator<const T> const_iterator;

    // iterator action
    iterator begin(){return iterator(a.begin());}
    const_iterator begin()const{return const_iterator(a.begin());}
    
    iterator end(){return iterator(a.end());}
    const_iterator end()const{return const_iterator(a.end());}
    
    // area functions
    size_t size()const{return a.GetSize();}
    bool empty()const{return a.IsEmpty();}
    void resize(size_t sz){a.Expand(sz);}
    
    // access to elements
    T& operator[](size_t n){return static_cast<T&>(a[n]);}
    const T& operator[](size_t n)const{return static_cast<const T&>(a[n]);}
    
    T& at(size_t n){return static_cast<T&>(a.At(n));}
    const T& at(size_t n)const {return static_cast<const T&>(a.At(n));}
    
    T& front(){return static_cast<T&>(a.First());}
    const T& front()const{return static_cast<const T&>(a.First());}
    
    T& back(){return static_cast<T&>(a.Last());}
    const T& back()const{return static_cast<const T&>(a.Last());}
    
    // modification
    void push_back(T &x){a.Add(x);}
    void pop_back(){a.RemoveAt(a.GetLast());}
    // insertion/removal is not supported in TRefArray!
    void clear(){a.Clear();}
    
  private:
    TRefArray a;
    
    ClassDef(ref_vector, 1);
  };  
*/  

  template<typename T> class const_ref_vector : public TObject {
  public:
    const_ref_vector(){}
    ~const_ref_vector(){}
    
    typedef ref_vector_iterator<const T *> iterator;

    // iterator action
    iterator begin()const{return iterator(a.begin());}
    iterator end()const{return iterator(a.end());}
    
    // area functions
    //size_t size()const{return a.GetSize();} // allocated slots size?
    size_t size()const{return a.GetLast()+1;} // size of non-empty slots
    bool empty()const{return a.IsEmpty();}
    void resize(size_t sz){a.Expand(sz);}
    
    // access to elements
    const T * operator[](size_t n)const{return static_cast<const T *>(a[n]);}
    const T * at(size_t n)const{return static_cast<const T *>(a.At(n));}
    
    const T * front()const{return static_cast<const T *>(a.First());}
    const T * back()const{return static_cast<const T *>(a.Last());}
    
    // modification
    void push_back(const T *x){a.Add(const_cast<T *>(x));}
    void pop_back(){a.RemoveAt(a.GetLast());}
    // insertion/removal is not supported in TRefArray!
    void clear(){a.Clear();}
    
  private:
    TRefArray a;
    
    ClassDefNV(const_ref_vector, 1);
  };  
  
  template<typename T> class ref_vector_iterator : public std::iterator<std::forward_iterator_tag, T> {
  public:
    ref_vector_iterator(const TIter &i) : it(i) {}
    ref_vector_iterator(const ref_vector_iterator<T> &i) : it(i.it) {}
    ~ref_vector_iterator();
    
    ref_vector_iterator<T> &operator=(const TIter &i){it = i; return *this;}
    ref_vector_iterator<T> &operator=(const ref_vector_iterator<T> &i){it = i.it; return *this;}
    
    T& operator*(){return static_cast<T>(*it);}
    T& operator->(){return static_cast<T>(*it);}
    
    ref_vector_iterator<T> &operator++(){it.Next(); return *this;}
    
    bool operator==(ref_vector_iterator<T> &ref){return (it == ref.it);}
    bool operator!=(ref_vector_iterator<T> &ref){return (it != ref.it);}

  private:
    TIter it;    
  };
  
}


#endif
