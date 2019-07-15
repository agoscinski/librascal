#ifndef SRC_MATH_INTERPOLATOR_HH_
#define SRC_MATH_INTERPOLATOR_HH_

#include <functional>
#include <forward_list>
#include <iostream>
#include <limits>
#include "math_utils.hh"

namespace rascal {
  namespace math {
    using Vector_Ref = typename Eigen::Ref<const Vector_t>;

    enum class GridType_t {Uniform};
    enum class RefinementMethod_t {HeapBased, Uniform, Adaptive};

    // TODO(alex) make plots of hyp1f1 normalized
    // TODO(alex) look at the graphs again, and make a grid type which is similar
    // to the shape
    // is similar to the graph

    // TODO(alex) currently the grid rational could be static
    template <GridType_t Type, RefinementMethod_t Method>
    struct GridRational {
      constexpr static GridType_t GridType { Type };
      constexpr static RefinementMethod_t RefinementMethod { Method };
    };

    //TODO(alex) adaptive: tree structure blocks{start,end,error,blocks}
    // refine where k largest error exist
    // get grid from (tree -> grid)
    template <>
    struct GridRational<GridType_t::Uniform, RefinementMethod_t::Adaptive> {
      GridRational<GridType_t::Uniform, RefinementMethod_t::Adaptive>() : grid_meshes{},grid_meshes_error{}, grid_size{0}, max_it{} {}

      //}
      // TODO for this grid rational it makes sense to save the test grid,
      // because it is one step ahead and we can add more points 
      // TODO make hyperparameter k highest error
      // TODO make grid rational init compute, increase_finness() and then remove x1,x2,fineness
      Vector_t compute_grid(double x1, double x2, int fineness) {
        if (fineness == 0) {
          this->grid_meshes = {x1, x2};
          this->grid_size = 2;
          return this->grid_from_meshes();
        }
        this->refine();
        return this->grid_from_meshes();
      }

      void refine() {
        auto next_it{this->max_it};
        next_it++;
        double cur{*this->max_it};
        double mid_point{cur+(*next_it-cur)/2};
        this->grid_meshes.emplace_after(this->max_it, mid_point);
        this->grid_size++;
        //if (this->grid_size % 1000 == 0) {
        //  std::cout << this->grid_size << std::endl;
        //}

      }

      Vector_t grid_from_meshes() {
        Vector_t grid = Vector_t::Zero(this->grid_size);
        int i{0};
        for (auto it=this->grid_meshes.begin(); it!=this->grid_meshes.end(); ++it) {
          grid(i) = (*it);
          i++;
        }
        return grid;
      }
      
      Vector_t compute_test_grid(double, double, int) {
        Vector_t test_grid = Vector_t::Zero(this->grid_size-1);
        int i{0};
        auto next_it{this->grid_meshes.begin()};
        next_it++;
        // if auto does not work use: std::forward_list<double>::iterator  
        for (auto it=this->grid_meshes.begin(); next_it!=this->grid_meshes.end(); ++it) {
          double cur{*it};
          double mid_point{cur+(*next_it-cur)/2};
          test_grid(i) = mid_point;

          next_it++;
          i++;
        }
        return test_grid;
      }

      void update_errors(Vector_Ref error_grid) {
        if (error_grid.size() != this->grid_size -1) {
          std::runtime_error("Gridsize does not match with error grid");
        }
        double max{0};
        auto it{this->grid_meshes.begin()};
        for (int i{0}; i < error_grid.size(); i++) {
          if (std::abs(error_grid(i)) > max) {
            max = std::abs(error_grid(i));          
            this->max_it = it;
          }
          it++;
        }
      }

      // two sortings
      // for refinement key = error
      // for constructing grid x1
      std::forward_list<double> grid_meshes;
      // TODO priority queue does not work well because splitting meshes does.
      // One could make Mesh only owning the x1 and error but the. For the
      // current test grid we need some clever binary tree structure where
      // leafs have errors and their parents add up, Node{x1,error}
      //std::priority_queue<Mesh> grid_meshes_error;
      std::forward_list<double> grid_meshes_error;
      int grid_size;
      std::forward_list<double>::iterator max_it;
    };

    // TODO(alex) I think they all can be static
    template <>
    struct GridRational<GridType_t::Uniform, RefinementMethod_t::HeapBased> {
      Vector_t compute_grid(double x1, double x2, int fineness) {
        double nb_grid_points = 2 << fineness;
        return Vector_t::LinSpaced(nb_grid_points, x1, x2);
      }
      Vector_t compute_test_grid(double x1, double x2, int fineness) {
        double nb_grid_points = 2 << fineness;
        // half of a step size in the grid to get the inner grid points
        // step size = (x2-x1)/(nb_grid_points-1)
        double offset{(x2-x1)/(2*(nb_grid_points-1))};
        return Vector_t::LinSpaced(nb_grid_points-1, x1+offset, x2-offset);
      }

      void update_errors(Vector_Ref) {} 
      int grid_size{0};
    };

    enum class InterpolationMethod_t {CubicSpline};

    template <InterpolationMethod_t Type>
    struct InterpolationMethod{};

    template <>
    class InterpolationMethod<InterpolationMethod_t::CubicSpline> {
     public:
      InterpolationMethod<InterpolationMethod_t::CubicSpline>(){}

      void initialize(const Vector_Ref & grid,
          const Vector_Ref & evaluated_grid){
        this->compute_second_derivatives_on_grid(grid, evaluated_grid);
      }

      // TODO(felix) the numerical recipes gives the option to set the first
      // derivative's starting and end point,
      // for now I did not include this option, I do not see now where
      // we would use it
      // TODO(alex) reference numerical recipes
      void compute_second_derivatives_on_grid(
          const Vector_Ref & grid, const Vector_Ref & evaluated_grid) {
        this->second_derivatives = this->sety2(grid, evaluated_grid);
      }

      double interpolate(const Vector_Ref & grid,
          const Vector_Ref & evaluated_grid,
          double x, size_t nearest_grid_index_to_x) {
        return this->rawinterp(grid, evaluated_grid,
            nearest_grid_index_to_x, x);
      }
     private:
      // This is done to be close to the numerical recipes implementation in
      // naming while making it more readable.
      // TODO(alex) reference numerical recipes
      Vector_t sety2(const Vector_Ref & xv, const Vector_Ref & yv) {
        int n{static_cast<int>(xv.size())};
        Vector_t y2 = Vector_t::Zero(n);
        Vector_t u = Vector_t::Zero(n);
        size_t sig;
        double p;
        y2(0) = 0.0;
        u(0) = 0.0;
        for (int i{1}; i<n-1; i++) {
          sig=(xv(i)-xv(i-1))/(xv(i+1)-xv(i-1));
          p=sig*y2(i-1)+2.0;
          y2(i)=(sig-1.0)/p;
          u(i)=(yv(i+1)-yv(i))/(xv(i+1)-xv(i)) - (yv(i)-yv(i-1))/(xv(i)-xv(i-1));
          u(i)=(6.0*u(i)/(xv(i+1)-xv(i-1))-sig*u(i-1))/p;
        }
        u(n-1) = 0.0;
        p=0.0;
        y2(n-1)=(u(n-1)-p*u(n-2))/(p*y2(n-2)+1.0);
        for (int k{n-2};k>0;k--) {
          y2(k)=y2(k)*y2(k+1)+u(k);
        }
        y2(0)=y2(0)*y2(1)+u(0);
        return y2;
      }

      double rawinterp(const Vector_Ref & xx, const Vector_Ref & yy,
          size_t j1, double x){
        size_t klo{j1}, khi{j1+1};
        const Vector_Ref y2 = Vector_Ref(this->second_derivatives);
        double h{xx(khi)-xx(klo)};
        if (h == 0.0) { throw ("Bad xa input to routine splint");}
        double a{(xx(khi)-x)/h};
        double b{(x-xx(klo))/h};
        return a*yy(klo)+b*yy(khi)+((a*a*a-a)*y2(klo)
          +(b*b*b-b)*y2(khi))*(h*h)/6.0;
      }

      Vector_t second_derivatives{};
    };

    enum class SearchMethod_t {Hunt};

    template <SearchMethod_t Type>
    struct SearchMethod{};

    template <>
    struct SearchMethod<SearchMethod_t::Hunt> {

      SearchMethod<SearchMethod_t::Hunt>() : correlated{false},
          nb_support_points{2}, last_accessed_index{0} {} 

      // If the requests to locate seem correlated, then the heuristic is used
      size_t search(double x, const Vector_Ref & grid) {
        return this->correlated ? this->hunt(x, grid) : this->locate(x, grid);
      } 

      // TODO(alex) move this to a Base class if we want to implement more
      // search methods
      // TODO(alex) ref numerical recipes
      size_t locate(double x, const Vector_Ref & xx) {
        int n{static_cast<int>(xx.size())};
        int mm{static_cast<int>(nb_support_points)};
        int jsav{static_cast<int>(this->last_accessed_index)};

        // TODO(alex) is this faster than pow(n, 0.25) ?
        int dj = std::min(1, 
            static_cast<int>(std::round(std::sqrt(std::sqrt(n)))));
        int ju,jm,jl;
        if (n < 2 || mm < 2 || mm > n) throw("locate size error");
        bool ascnd=(xx[n-1] >= xx[0]);
        jl=0;
        ju=n-1;
        while (ju-jl > 1) {
          jm = (ju+jl) >> 1;
          if ((x >= xx[jm]) == ascnd)
            jl=jm;
          else
            ju=jm;
        }
        this->correlated = abs(jl-jsav) > dj ? 0 : 1;
        jsav = jl;

        this->last_accessed_index = jsav;
        return std::max(0,std::min(n-mm,jl-((mm-2)>>1)));
      }

      // TODO(alex) change xx to const Ref & because it is not modified 
      size_t hunt(double x, const Vector_Ref & xx){
        int n{static_cast<int>(xx.size())};
        int mm{static_cast<int>(nb_support_points)};
        int dj = std::min(1, 
            static_cast<int>(std::round(std::sqrt(std::sqrt(n)))));
        int jsav{static_cast<int>(this->last_accessed_index)};

        int jl=jsav, jm, ju, inc=1;
        if (n < 2 || mm < 2 || mm > n) throw("hunt size error");
        bool ascnd=(xx[n-1] >= xx[0]);
        if (jl < 0 || jl > n-1) {
          jl=0;
          ju=n-1;
        } else {
          if ((x >= xx[jl]) == ascnd) {
            for (;;) {
              ju = jl + inc;
              if (ju >= n-1) { ju = n-1; break;}
              else if ((x < xx[ju]) == ascnd) break;
              else {
                jl = ju;
                inc += inc;
              }
            }
          } else {
            ju = jl;
            for (;;) {
              jl = jl - inc;
              if (jl <= 0) { jl = 0; break;}
              else if ((x >= xx[jl]) == ascnd) break;
              else {
                ju = jl;
                inc += inc;
              }
            }
          }
        }
        while (ju-jl > 1) {
          jm = (ju+jl) >> 1;
          if ((x >= xx[jm]) == ascnd)
            jl=jm;
          else
            ju=jm;
        }
        this->correlated = abs(jl-jsav) > dj ? 0 : 1;
        jsav = jl;
        return std::max(0,std::min(n-mm,jl-((mm-2)>>1)));
      }

      bool correlated;
      size_t nb_support_points;
      size_t last_accessed_index;
    };

    template<class InterpolationMethod, class GridRational, class SearchMethod>
    class Interpolator {
     public: 
      Interpolator() : mean_error{0}, max_error{0}, intp_method{InterpolationMethod()}, grid_rational{GridRational()}, search_method{SearchMethod()} {}

      void initalize(std::function<double(double)> function, double x1, double x2, double precision) {
        if (x2<x1) {
          throw std::runtime_error("x2 must be greater x1");
        }
        this->function = function;
        this->x1 = x1;
        this->x2 = x2;
        this->precision = precision;

        this->initialize_interpolator();
      }

      void initialize_interpolator() {
        // Fineness starts with zero and is incremently increased
        // this definition is arbitrary but make computation more readable
        this->fineness = 0;
        double error{this->compute_grid_error()};
        // TODO(alex) add some procedure to not get locked if precision is too
        // high
        while (error > this->precision) {
          this->fineness++;
          error = this->compute_grid_error();
        }
      }

      // TODO(alex) if I use temporary variables instead of this, does the
      // compiler optimize this?
      double compute_grid_error() {
        this->grid = 
            this->grid_rational.compute_grid(this->x1,this->x2, this->fineness);
        this->evaluated_grid = this->eval(this->grid);

        this->intp_method.initialize(this->grid, this->evaluated_grid);

        Vector_t test_grid{this->grid_rational.compute_test_grid(this->x1,this->x2,this->fineness)};
        Vector_t test_grid_interpolated{this->interpolate(test_grid)};
        Vector_t test_grid_evaluated{this->eval(test_grid)}; 
        // computes the relative error
        Vector_t error_grid{2*((test_grid_interpolated - test_grid_evaluated).array()/
          (std::numeric_limits< double >::min()+test_grid_interpolated.array().abs() + test_grid_evaluated.array().abs())).abs()};
        grid_rational.update_errors(Vector_Ref(error_grid));
        this->max_error = error_grid.maxCoeff();
        this->mean_error = error_grid.mean();
        //if (this->grid.size() % 1000==0) {
        //  std::cout << "grid_size=" << this->grid.size() << std::endl;
        //  std::cout << "mean_error=" << this->mean_error << std::endl;
        //}
        //if (grid_rational.grid_size % 50 == 0) {
        //  std::cout << "fineness=" << this->fineness << std::endl;
        //  std::cout << "mean error=" << this->mean_error << std::endl;
        //  std::cout << "max error=" << this->max_error << std::endl;
        //}
        return this->mean_error;
      }

      double eval(double x) {return this->function(x);}

      // We use evaluate when the function is used and interpolate when the
      // interpolation method is used
      Vector_t eval(const Vector_Ref & grid) {
        Vector_t evaluated_grid = Vector_t::Zero(grid.size());
        for (int i{0}; i<evaluated_grid.size(); i++) {
          evaluated_grid(i) = this->eval(grid(i));
        }
        return evaluated_grid;
      }

      // TODO(alex) test this assumption:
      // this should save one copy operation instead of the above used for 
      // this->evaluated_grid 
      //void eval() {
      //  for (size_t i{0}; i<this->evaluated_grid.size(); i++) {
      //    this->evaluated_grid(i) = this->function(grid(i));
      //  }
      //}

      double interpolate(double x) {
        // TODO(alex) throw runtime error, what is diff?
        if (x<this->x1) { throw std::runtime_error ("x is outside of range, below x1"); }
        if (x>this->x2) { throw std::runtime_error ("x is outside of range, above x2"); }
        size_t nearest_grid_index_to_x{this->search_method.search(x, this->grid)};
        return intp_method.interpolate(
            this->grid, this->evaluated_grid,
            x, nearest_grid_index_to_x);
      }

      Vector_t interpolate(const Vector_Ref & points) {
        Vector_t interpolated_points = Vector_t::Zero(points.size());
        for (int i{0}; i<points.size(); i++) {
          interpolated_points(i) = this->interpolate(points(i));
        }
        return interpolated_points;
      }
     
      std::function<double(double)> function{};
      double x1{0};
      double x2{1};
      double precision{1e-5};
      double mean_error;
      double max_error;
      int fineness{0};
      Vector_t grid{};
      Vector_t evaluated_grid{};

      InterpolationMethod intp_method;
      GridRational grid_rational;
      SearchMethod search_method;
    };


  // TODO(alex) make a CRTP calculator and check if this can be merged with the GradientCalutaro stuff

  }  // namespace math
}  // namespace rascal
#endif  // SRC_MATH_INTERPOLATOR_HH_