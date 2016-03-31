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

    string input_path_;
    string output_path_;
    pugi::xml_document output_file_;

    int number_of_points_;
    int number_of_velocities_;
    int number_of_angles_;
    int number_of_elements_;
    int number_of_time_steps_;
    int dump_number_;
    
    double point_distance_;
    double velocity_distance_;
    double angle_distance_;
    double time_step_;
    double const q_ = 1;
    double const qm_ = 1;

    Array<double> position_;
    Array<double> velocity_;
    Array<double> angle_;
    
    Array<double> electric_field_;
    Array<double> magnetic_field_;
    Array<double> charge_density_;
    
    Array<double> density_;
    Array<double> mean_density_;

    const int index_base_ = 0;

    vector<int> number_of_entries_per_row_;
    vector<int> small_number_of_entries_per_row_;
    
    string solver_type_ = "Klu";
    Amesos factory_;

    unique_ptr<Epetra_SerialComm> comm_;
    unique_ptr<Epectra_Map> map_;
    unique_ptr<Epetra_CrsMatrix> matrix_;
    unique_ptr<Epetra_Vector> lhs_;
    unique_ptr<Epetra_Vector> rhs_;
    unique_ptr<Epetra_LinearProblem> problem_;
    unique_ptr<Amesos_BaseSolver> solver_;

    unique_ptr<Epetra_SerialComm> small_comm_;
    unique_ptr<Epectra_Map> small_map_;
    unique_ptr<Epetra_CrsMatrix> small_matrix_;
    unique_ptr<Epetra_Vector> small_lhs_;
    unique_ptr<Epetra_Vector> small_rhs_;
    unique_ptr<Epetra_LinearProblem> small_problem_;
    
    unique_ptr<Amesos_BaseSolver> small_solver_;

    void fill_matrix();
    void fill_lhs();
    void fill_rhs();
    
    void calculate_charge_density();
    void calculate_electric_field();
    
    void step();
    void solve();

    int check_point(int p);
    int check_velocity(int g);
    int check_angle(int o);
    
    double spatial_constant(int g, int o);
    double velocity_constant(int p, int g, int o);
    double angle_constant(int p, int g, int o);
    
    int subscript_to_index(int p, int g, int o)
    {
        return o + number_of_angles * (g + number_of_velocities * p);
    }
    
public:
    
    FD_Vlasov(int number_of_cells,
              int number_of_velocities,
              int number_of_angles,
              int number_of_time_steps);

    void solve();
    
}

#endif
