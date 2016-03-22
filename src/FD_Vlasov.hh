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

using std::string;
using std::unique_ptr;
using std::vector;

class FD_Vlasov_hh
{
private:

    int number_of_cells_;
    int number_of_groups_;
    int number_of_ordinates_;
    int number_of_elements_;
    int number_of_time_steps_;
    
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
    void calculate_force();
    void step();

public:
    
    FD_Vlasov(int number_of_cells,
              int number_of_groups,
              int number_of_ordinates,
              int number_of_time_steps);

    void solve();
    
}

#endif
