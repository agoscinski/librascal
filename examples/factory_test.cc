#include "json_io.hh"
// #include "type_list.hpp"
// #include "meta_counter.hpp"
// #include "meta_list.hpp"

#include <iostream>
#include <memory>
#include <string>
#include <unordered_map>

#include <cstdlib>
#include <cxxabi.h>
#include <regex>
#include <map>
#include <type_traits>

std::string demangle(const char *name) {

  int status = -4; // some arbitrary value to eliminate the compiler warning

  std::unique_ptr<char, void (*)(void *)> res{
      abi::__cxa_demangle(name, NULL, NULL, &status), std::free};

  return (status == 0) ? res.get() : name;
}


template<class Impl>
class StrMngr {
    public:
    StrMngr() {std::cout<<"StrMngr"<<std::endl;}
};

class Impl1 : public StrMngr<Impl1> {
    public:
    Impl1() {std::cout<<"Impl1"<<std::endl;}

    void update(std::vector<double>){
      std::cout<<"update Impl1"<<std::endl;
    }

    decltype(auto) get_name(){
        return demangle(typeid(Impl1).name());
    }

};

template<class Impl>
class AdImpl1 : public StrMngr<AdImpl1<Impl>> {
    public:
    using ImplPtr = std::shared_ptr<Impl>;
    using Manager_t = AdImpl1<Impl>;
    using ManagerPtr = std::shared_ptr<Manager_t>;

    AdImpl1(ImplPtr& m, double):manager{m}
    {}

    AdImpl1(ImplPtr& m, std::tuple< double> a)
    :AdImpl1(m, std::get<0>(a) )
    {}

    ImplPtr manager;

    void update() {
      std::cout<<"update "<< this->get_name() <<std::endl;
    }

    //! Updates the underlying manager as well as the adaptor
    template <class... Args>
    void update(Args &&... arguments) {
      this->manager->update(arguments...);
      this->update();
    }


    decltype(auto) get_name(){
        return demangle(typeid(Manager_t).name());
    }
};

template<class Impl>
class AdImpl2 : public StrMngr<AdImpl2<Impl>> {
    public:
    using ImplPtr = std::shared_ptr<Impl>;
    using Manager_t = AdImpl2<Impl>;
    using ManagerPtr = std::shared_ptr<Manager_t>;
    AdImpl2(ImplPtr& m, double, int)
    :manager{m}
    {}

    AdImpl2(ImplPtr& m, std::tuple< double, int> a)
    :AdImpl2(m, std::get<0>(a),std::get<1>(a))
    {}

    ImplPtr manager;

    void update() {
      std::cout<<"update "<< this->get_name() <<std::endl;
    }

    //! Updates the underlying manager as well as the adaptor
    template <class... Args>
    void update(Args &&... arguments) {
      this->manager->update(arguments...);
      this->update();
    }

    decltype(auto) get_name(){
        return demangle(typeid(Manager_t).name());
    }
};


template< typename Impl,  template<class> class AdImpl,
                            template<class> class ... AdImplr>
struct Recursion {
    using Manager_t = AdImpl<Impl>;
    using ImplPtr = std::shared_ptr<Impl>;
    using ManagerPtr = std::shared_ptr<Manager_t>;
    using type = Recursion<Manager_t,AdImplr...>;

    template<typename Arg, typename ...Args>
    Recursion(ImplPtr& m, Arg& arg, Args& ...args)
    :manager{std::make_shared<Manager_t>(m,arg)},
    next_stack{manager, args...}
    {std::cout<<this->manager->get_name()<<std::endl;}

    template<typename Arg1, typename Arg2>
    Recursion(ImplPtr& m, Arg1& arg1, Arg2& arg2)
    :manager{std::make_shared<Manager_t>(m,arg1)},
    next_stack{manager, arg2}
    {std::cout<<this->manager->get_name()<<std::endl;}

    ManagerPtr manager;
    type next_stack;


    decltype(auto) get_manager() {
        return this->next_stack.get_manager();
    }
};

template< typename Impl,  template<class> class AdImpl>
struct Recursion<Impl,AdImpl> {
    using Manager_t = AdImpl<Impl>;
    using ImplPtr = std::shared_ptr<Impl>;
    using ManagerPtr = std::shared_ptr<Manager_t>;
    using type = Manager_t;

    template<typename Arg>
    Recursion(ImplPtr& m, Arg& arg)
    :manager{std::make_shared<Manager_t>(m,arg)}
    {std::cout<<this->manager->get_name()<<std::endl;}

    ManagerPtr manager;

    decltype(auto) get_manager() {
        return this->manager;
    }
};

template <typename Impl, template<class> class ... AdImplr, typename Arg, typename ...Args >
decltype(auto) structure_manager_factory(Arg arg, Args ...args) {
  auto manager_base = std::make_shared<Impl>(Impl{});

  auto test = Recursion<Impl1, AdImplr... >(manager_base, args...);

  auto manager = test.get_manager();

  std::cout  <<std::endl;
  std::cout << manager->get_name() <<std::endl;

  manager->update(arg);

  return manager;

}


template<typename T>
auto get_info() {
  return typeid(T).name();
}


int main() {
  std::vector<double> vv{};
  vv.push_back(14);
  // auto manager_base = std::make_shared<Impl1>(Impl1{});

//   double d{10},e{10};
  auto d = std::make_tuple<double>(1.5);
  auto tp = std::make_tuple<double,int>(1.5,10);
  auto tp1 = std::make_tuple<double,int>(1.5,10);


  auto manager = structure_manager_factory<Impl1,AdImpl1,AdImpl2,AdImpl2>(vv,d,tp,tp1);




  return(0);
}
