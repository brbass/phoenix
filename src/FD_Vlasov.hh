#ifndef FD_Vlasov_hh
#define FD_Vlasov_hh

#include <memory>
#include <string>
#include <vector>

#include <Amesos.h>
#include <AztecOO.h>
#include <AztecOO_Version.h>
// #ifdef EPETRA_MPI
// #  include <mpi.h>
// #  include <Epetra_MpiComm.h>
// #else
# include <Epetra_SerialComm.h>
// #endif
#include <Epetra_Map.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Vector.h>
#include <Epetra_LinearProblem.h>

#include "Array.hh"

using std::string;
using std::unique_ptr;
using std::vector;

class FD_Vlasov_hh
{
private:

    int const number_of_cells_;
    int const number_of_groups_;
    int const number_of_angles_;
    int const number_of_elements_;
    int const number_of_time_steps_;
    int const dump_number_;
    
    double const qm_ = 1;

    Array<double> force_;
    Array<double> electric_field_;
    Array<double> magnetic_field_;

    Array<double> angle_;
    Array<double> position_;
    Array<double> velocity_;
    Array<double> angle_length_;
    Array<double> position_length_;
    Array<double> velocity_length_;
    
    Array<double> density_;
    Array<double> density_old_;

    const int index_base_ = 0;

    vector<int> number_of_entries_per_row_;
    vector<int> my_global_elements_; // global indeces of elements (rows)
    vector<int> row_map_; // position of elements from diagonal
    
    unique_ptr<Epetra_SerialComm> comm_;
    unique_ptr<Epectra_Map> map_;
    unique_ptr<Epetra_CrsMatrix> matrix_;
    unique_ptr<Epetra_Vector> lhs_;
    unique_ptr<Epetra_Vector> rhs_;
    unique_ptr<Epetra_LinearProblem> problem_;
    
    string solver_type_ = "Klu";
    Amesos factory_;
    unique_ptr<Amesos_BaseSolver> solver_;

    void fill_matrix();
    void fill_lhs();
    void fill_rhs();
    
    void calculate_charge_density();
    void calculate_electric_field();
    void calculate_force();
    
    void step();
    void solve();

public:
    
    FD_Vlasov(int number_of_cells,
              int number_of_groups,
              int number_of_angles,
              int number_of_time_steps);

    void solve();
    
}

#endif
