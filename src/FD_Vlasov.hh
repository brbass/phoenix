#ifndef FD_Vlasov_hh
#define FD_Vlasov_hh

#include <cmath>
#include <memory>
#include <numeric>
#include <string>
#include <vector>

// #include <Amesos.h>
#include <AztecOO.h>
#ifdef EPETRA_MPI
#  include <mpi.h>
#  include <Epetra_MpiComm.h>
#else
# include <Epetra_SerialComm.h>
#endif
#include <Epetra_Map.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_Vector.h>
#include <Epetra_LinearProblem.h>

#include "pugixml.hpp"

#include "Check.hh"

using std::floor;
using std::string;
using std::unique_ptr;
using std::vector;

class FD_Vlasov
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

    vector<double> position_;
    vector<double> velocity_;
    vector<double> angle_;
    
    vector<double> electric_field_;
    vector<double> magnetic_field_;
    vector<double> charge_density_;
    
    vector<double> density_;
    vector<double> mean_density_;

    const int index_base_ = 0;

    vector<int> number_of_entries_per_row_;
    vector<int> charge_number_of_entries_per_row_;
    
    int max_iterations_ = 5000;
    double tolerance_ = 1e-6;
    int solver_print_;
    int kspace_;
    int poly_ord_;

#ifdef EPETRA_MPI
    unique_ptr<Epetra_MpiComm> comm_;
#else
    unique_ptr<Epetra_SerialComm> comm_;
#endif
    unique_ptr<Epetra_Map> map_;
    unique_ptr<Epetra_CrsMatrix> matrix_;
    unique_ptr<Epetra_Vector> lhs_;
    unique_ptr<Epetra_Vector> rhs_;
    unique_ptr<Epetra_LinearProblem> problem_;
    unique_ptr<AztecOO> solver_;

#ifdef EPETRA_MPI
    unique_ptr<Epetra_MpiComm> charge_comm_;
#else
    unique_ptr<Epetra_SerialComm> charge_comm_;
#endif
    unique_ptr<Epetra_Map> charge_map_;
    unique_ptr<Epetra_CrsMatrix> charge_matrix_;
    unique_ptr<Epetra_Vector> charge_lhs_;
    unique_ptr<Epetra_Vector> charge_rhs_;
    unique_ptr<Epetra_LinearProblem> charge_problem_;
    unique_ptr<AztecOO> charge_solver_;

    void parse_xml();
    void initialize_trilinos();
    void fill_matrix();
    void fill_charge_matrix();
    
    void calculate_charge_density();
    void calculate_electric_field();
    void calculate_density();

    int check_point(int p);
    int check_velocity(int g);
    int check_angle(int o);
    
    double spatial_constant(int g, int o);
    double velocity_constant(int p, int g, int o);
    double angle_constant(int p, int g, int o);

    void initial_xml();
    void dump_xml(int i, double t);
    void save_xml();

    int subscript_to_index(int p, int g, int o)
    {
        return o + number_of_angles_ * (g + number_of_velocities_ * p);
    }
    
    void index_to_subscript(int n, vector<int> &subscript)
    {
        Check(subscript.size() == 3);

        int product = number_of_velocities_*number_of_angles_;
        int sum = n;
        
        subscript[0] = floor(static_cast<double>(sum) / (number_of_velocities_ * number_of_angles_));
        sum -= number_of_velocities_ * number_of_angles_ * subscript[0];
        subscript[1] = floor(static_cast<double>(sum) / number_of_angles_);
        sum -= number_of_angles_ * subscript[1];
        subscript[2] = sum;
    }

    
public:
    
    FD_Vlasov(string input_path);
    
    void check_class_invariants();
    void solve();
};

#endif
