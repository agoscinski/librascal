#include <../src/neighbourhood_managers/neighbourhood_manager_cell.hh>
#include <../src/neighbourhood_managers/neighbourhood_manager_lammps.hh>
#include <iostream>
#include <../src/basic_types.h>
#include <Eigen/StdVector>
#include <cmath>
using namespace std;

using Manager_t = rascal::NeighbourhoodManagerCell;
constexpr static int Natom{8};
constexpr static int dim{3};
using ptr_t = double**;
using vVector3d = std::vector<Eigen::Vector3d,Eigen::aligned_allocator<Eigen::Vector3d> >;
int main()
{   
    Manager_t manager;
    Eigen::MatrixXd cell(3,3);
    cell << 6.19,2.41,0.21,
            0.00,6.15,1.02,
            0.00,0.00,7.31;
    Eigen::MatrixXd positions(3,22); // 3,22
    positions << 3.689540159937393,5.123016813620886,1.994119731169116,6.818437242389163,2.630056617829216,6.182500355729062,2.114977334498767,6.697579639059512,1.392155450018263,7.420401523540017,2.432242071439904,6.380314902118375,1.112656394115962,7.699900579442317,3.569715877854675,5.242841095703604,3.122826344932127,5.689730628626151,3.248684682453303,5.563872291104976,2.608353462112637,6.204203511445642,
                5.035681855581504,2.134827911489532,0.946910011088814,6.223599755982222,4.168634519120968,3.001875247950068,1.980327734683430,5.190182032387606,2.943861424421339,4.226648342649697,5.457161501166098,1.713348265904937,1.501663178733906,5.668846588337130,5.208365510425203,1.962144256645833,2.728127406527150,4.442382360543885,2.839975217222644,4.330534549848392,0.744216089807768,6.426293677263268,
                4.643695520786083,2.662204050783991,1.250682335857938,6.055217235712136,0.860905287815103,6.444994283754972,4.536108843695142,2.769790727874932,5.609177455068640,1.696722116501434,6.703053268421970,0.602846303148105,3.487609972580834,3.818289598989240,1.436734374347541,5.869165197222533,1.054504320562138,6.251395251007936,3.998423858825871,3.307475712744203,5.323662899811682,1.982236671758393;
    std::array<bool,3> pbc{{true,true,true}};
    double cutoff_max{3};  
    std::vector<int> center_ids;
    for (int ii{0}; ii < 22; ++ii) center_ids.push_back(ii);
    manager.build(positions,center_ids,cell,pbc,cutoff_max);
    
    Eigen::MatrixXd positions_test(3,22); // 3,22
    positions_test << 3.689540159937393,5.123016813620886,1.994119731169116,6.818437242389163,2.630056617829216,6.182500355729062,2.114977334498767,6.697579639059512,1.392155450018263,7.420401523540017,2.432242071439904,6.380314902118375,1.112656394115962,7.699900579442317,3.569715877854675,5.242841095703604,3.122826344932127,5.689730628626151,3.248684682453303,5.563872291104976,2.608353462112637,6.204203511445642,
                5.035681855581504,2.134827911489532,0.946910011088814,6.223599755982222,4.168634519120968,3.001875247950068,1.980327734683430,5.190182032387606,2.943861424421339,4.226648342649697,5.457161501166098,1.713348265904937,1.501663178733906,5.668846588337130,5.208365510425203,1.962144256645833,2.728127406527150,4.442382360543885,2.839975217222644,4.330534549848392,0.744216089807768,6.426293677263268,
                4.643695520786083,2.662204050783991,1.250682335857938,6.055217235712136,0.860905287815103,6.444994283754972,4.536108843695142,2.769790727874932,5.609177455068640,1.696722116501434,6.703053268421970,0.602846303148105,3.487609972580834,3.818289598989240,1.436734374347541,5.869165197222533,1.054504320562138,6.251395251007936,3.998423858825871,3.307475712744203,5.323662899811682,1.982236671758393;
        
    rascal::VecXi neighlist;
    
    int atom_counter{};
    int pair_counter{};
    constexpr bool verbose{false};
    for (auto center: manager) {
        if (atom_counter - center.get_index() != 0 ){
            cout << "index " << center.get_index() << endl;
        }
        ++atom_counter;
        cout << "Center atom index: " << center.get_index() << endl;
        for (int ii{3};ii<3;++ii){
            if (positions_test(ii,center.get_index()) - center.get_position()[ii] != 0 ){
                cout << "position " << center.get_index() << endl;
            }
        }
        cout << "Neighbour indices: ";
        for (auto neigh : center){
            cout  << neigh.get_index() << ", ";
        }
        cout << endl;
    }

    return(0);
}