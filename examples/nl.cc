#include <../src/neighbourhood_managers/neighbourhood_manager_cell.hh>
#include <../src/neighbourhood_managers/neighbourhood_manager_lammps.hh>
#include <iostream>
#include <../src/basic_types.h>
using namespace std;

using Manager_t = proteus::NeighbourhoodManagerCell;
constexpr static int Natom{8};
constexpr static int dim{3};
using ptr_t = double**;

int main()
{

  cout << Natom << " "<< dim << " test " << endl;
  Eigen::MatrixXd pos(3,8);
  pos << 0.00,2.68,1.79,0.00,1.79,0.89,0.89,2.68,
            0.00,2.68,0.00,1.79,1.79,2.68,0.89,0.89,
            1.79,2.68,0.00,0.00,1.79,0.89,2.68,0.89;
  cout << "Now the array a is:" << endl << pos << endl;
  Eigen::MatrixXd cell(3,3);
  cell << 3.57,0.00,0.00,
            0.00,3.57,0.00,
            0.00,0.00,3.57;
  Eigen::VectorXi num(8);
  num << 6, 6, 6, 6, 6, 6, 6, 6;
  std::array<bool,3> pbc = {1,1,1};
    //int nb_pairs;
  Manager_t manager;
  
  double rc_max{4};
  manager.build(pos,cell,pbc,rc_max);
  //manager.set_positions(pos);
  
  /*
  for (auto center:manager){
      cout << "Center id: " << center.get_index() << endl;
      cout << "Neighbour ids: " ;
      
      for (auto neigh : center){
          cout << neigh.get_index() << ", "<< endl;
      }
      cout <<  endl;
  }
  */
  /*
  for (auto atom: manager) {
    for (auto pair: atom) {
        cout << " distance between atom " << atom.get_index() << " and atom " << pair.get_index() << "  " << (atom.get_position() - pair.get_position()).norm()  << endl;
    }
  }
  
  cout << manager.get_nb_clusters(1) << " test " << endl;
  */
  //cout <<  nb_pairs << " test " << endl;


    return(0);
}