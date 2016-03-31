#include "FD_Vlasov.hh"

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "pugixml.hpp"

#include "Array.hh"
#include "Check.hh"

using namespace std;

FD_Vlasov::
FD_Vlasov(string filename):
    input_path_(filename)
{
    parse_xml();
    initialize_trilinos();
    check_class_invariants();
}

void FD_Vlasov::
check_class_invariants()
{
    Check(number_of_elements_ == number_of_points_ * number_of_velocities_ * number_of_angles_, "number_of_elements");

    Check(electric_field_.size() == number_of_points_, "electric_field size");
    Check(magnetic_field_.size() == number_of_points_, "magnetic_field size");
    Check(charge_density_.size() == number_of_points_, "charge_density size");

    Check(angle_.size() == number_of_ordinates_, "angle size");
    Check(position_.size() == number_of_points_, "position size");
    Check(velocity_.size() == number_of_velocities_, "velocity size");

    Check(density_.size() == number_of_elements_, "density size");
    Check(mean_density_.size() == number_of_points_ * number_of_velocities_, "mean_density size");

    Check(number_of_entries_per_row_.size() == number_of_elements_, "number_of_entries_per_row size");
    Check(small_number_of_entries_per_row_.size() == number_of_points_, "small_number_of_entries_per_row size");
}

void FD_Vlasov::
parse_xml(string filename)
{
    output_path_ = input_path_.substr(0, input_path_.find_last_of(".")) + "-out.xml";
    output_path_.append_child("output");
    
    pugi::xml_document doc;
    
    if (!doc.load_file(filename.c_str()))
    {
        Check(false, "could not open document");
    }
    
    pugi::xml_node input = doc.child("input");
    
    pugi::xml_node point_node = input.child("point_discretization");
    pugi::xml_node velocity_node = input.child("velocity_discretization");
    pugi::xml_node angle_node = input.child("angle_discretization");
    pugi::xml_node time_node = input.child("time_discretization");
    
    number_of_points_ = child_value<int>(point_node, "number_of_points");
    number_of_velocities_ = child_value<int>(velocity_node, "number_of_velocities");
    number_of_angles_ = child_value<int>(angle_node, "number_of_angles");
    number_of_elements_ = number_of_points_ * number_of_velocities_ * number_of_angles_;
    number_of_time_steps_ = child_value<int>(time_node, "number_of_time_steps");
    dump_number_ = child_value<int>(time_node, "dump_number");
    
    point_distance_ = child_value<double>(point_node, "point_distance");
    velocity_distance_ = child_value<double>(velocity_node, "velocity_distance");
    angle_distance_ = 2.0 * M_PI / static_cast<double>(number_of_angles_);
    time_step_ = child_value<double>(time_node, "time_step");
    
    position_.resize(number_of_points_);
    for (unsigned i = 0; i < number_of_points_; ++i)
    {
        position_(i) = i * point_distance_;
    }
    
    velocity_.resize(number_of_velocities_);
    for (unsigned i = 0; i < number_of_velocities_; ++i)
    {
        velocity_(i) = i * velocity_distance_;
    }

    angle_.resize(number_of_angles_);
    for (unsigned i = 0; i < number_of_angles_; ++i)
    {
        angle_(i) = i * angle_distance_;
    }
    
    // start again at "electric field"
}

Array<double> FD_Vlasov::
get_length(int number,
           Array<double> &array)
{
    Array<double> length(number - 1);
    
    for (unsigned i = 0; i < number - 1; ++i)
    {
        length(i) = array(i + 1) - array(i);
    }
    
    return length;
}

void FD_Vlasov::
initialize_trilinos()
{
    // main matrix

    number_of_entries_per_row_.resize(number_of_elements_, 7);

    // initialize communications classes

    comm_ = unique_ptr<Epetra_SerialComm> (new Epetra_SerialComm);
    map_ = unique_ptr<Epetra_Map> ( new Epetra_Map(number_of_elements_, index_base_, *comm_));

    // initialize and fill matrix

    matrix_ = unique_ptr<Epetra_CRSMatrix> (new Epetra_CrsMatrix(Copy, *map_, &number_of_entries_per_row_[0], true));
    lhs_ = unique_ptr<Epetra_Vector> (new Epetra_Vector(*map_));
    rhs_ = unique_ptr<Epetra_Vector> (new Epetra_Vector(*map_));
    fill_matrix();
    
    problem_ = unique_ptr<Epetra_LinearProblem> (new Epetra_LinearProblem(matrix_.get(), lhs_.get(), rhs_.get()));
    solver_ = unique_ptr<Amesos_BaseSolver> (factory_.Create(solver_type_, *problem_));

    solver_->SymbolicFactorization();
    
    // charge density matrix

    small_number_of_entries_per_row_.resize(number_of_points_, 3);

    // initialize communications classes

    small_comm_ = unique_ptr<Epetra_SerialComm> (new Epetra_SerialComm);
    small_map_ = unique_ptr<Epetra_Map> ( new Epetra_Map(number_of_elements(), index_base, *comm));

    // initialize and fill matrix
    
    small_matrix_ = unique_ptr<Epetra_CRSMatrix> (new Epetra_CrsMatrix(Copy, *small_map_, &small_number_of_entries_per_row_[0], true));
    small_lhs_ = unique_ptr<Epetra_Vector> (new Epetra_Vector(*small_map_));
    small_rhs_ = unique_ptr<Epetra_Vector> (new Epetra_Vector(*small_map_));
    fill_small_matrix();
    
    small_problem_ = unique_ptr<Epetra_LinearProblem> (new Epetra_LinearProblem(small_matrix_.get(), small_lhs_.get(), small_rhs_.get()));
    small_solver_ = unique_ptr<Amesos_BaseSolver> (factory_.Create(solver_type_, *small_problem_));

    small_solver_->SymbolicFactorization();
    
}

void FD_Vlasov::
fill_matrix()
{
    for (int i = 0; i < number_of_points_; ++i)
    {
        int im1 = check_point(i - 1);
        int ip1 = check_point(i + 1);
        
        for (int j = 0; j < number_of_velocities_; ++j)
        {
            int jm1 = check_group(j - 1);
            int jp1 = check_group(j + 1);
            
            for (int k = 0; k < number_of_angles_; ++k)
            {
                int km1 = check_angle(k - 1);
                int kp1 = check_angle(k + 1);
                
                int n = subscript_to_index(i, j, k);
                int nim = subscript_to_index(im1, j, k);
                int nip = subscript_to_index(ip1, j, k);
                int njm = subscript_to_index(i, jm1, k);
                int njp = subscript_to_index(i, jp1, k);
                int nkm = subscript_to_index(i, j, km1);
                int nkp = subscript_to_index(i, j, kp1);
                
                vector<int> column_indices(number_of_entries_per_row_[n], 0);
                vector<double> fill_vector(number_of_entries_per_row_[n], 0);
                
                column_indices[0] = nim;
                column_indices[1] = njm;
                column_indices[2] = nkm;
                column_indices[3] = n;
                column_indices[4] = nkp;
                column_indices[5] = njp;
                column_indices[6] = nip;
                
                fill_vector[0] = -spatial_constant(j, k) / (2 * point_distance_);
                fill_vector[1] = -velocity_constant(i, jm1, k) / (2 * velocity_distance_);
                fill_vector[2] = -angle_constant(i, j, km1) / (2 * angle_distance_);
                fill_vector[3] = 2 / time_step_;
                fill_vector[4] = angle_constant(i, j, kp1) / (2 * angle_distance_);
                fill_vector[5] = velocity_constant(i, jp1, k) / (2 * velocity_distance_);
                fill_vector[6] = spatial_constant(j, k) / (2 * point_distance_);

                // insert values into the matrix for the lhs
                matrix_->InsertGlobalValues(n, number_of_entries_per_row_[n], &fill_vector[0], &column_indices[0]);
                
                // insert values into the rhs
                double sum = 0;
                sum += fill_vector[0] * density_(im1, j, k);
                sum += fill_vector[1] * density_(i, jm1, k);
                sum += fill_vector[2] * density_(i, j, km1);
                sum += fill_vector[3] * density_(i, j, k);
                sum += fill_vector[4] * density_(i, j, kp1);
                sum += fill_vector[5] * density_(i, jp1, k);
                sum += fill_vector[6] * density_(ip1, j, k);
                sum += 2 * density_(i, j, k) / time_step_;
                
                (*rhs_)[n] = sum;
            }
        }
    }
    
    matrix_->FillComplete();
    lhs_->PutScalar(1.0);
}

void FD_Vlasov::
fill_small_matrix()
{
    for (unsigned i = 0; i < number_of_cells_; ++i)
    {
        vector<int> column_indices(small_number_of_entries_per_row_[i], 0);
        vector<double> fill_vector(small_number_of_entries_per_row_[i], 0);
        
        int im1 = check_point(i - 1);
        int ip1 = check_point(i + 1);

        column_indices[0] = im1;
        column_indices[1] = i;
        column_indices[2] = ip1;
        
        fill_vector[0] = -electric_field_(im1) / (2 * point_distance_);
        fill_vector[1] = 0;
        fill_vector[2] = electric_field_(ip1) / (2 * point_distance_);
        
        small_matrix_->InsertGlobalValues(i, small_number_of_entries_per_row_[i], &fill_vector[0], &column_indices[0]);
        
        (*small_rhs_)[i] = charge_density_(i);
    }
    
    small_matrix_->FillComplete();
    small_lhs_->PutScalar(1.0);
}

void FD_Vlasov::
calculate_charge_density()
{
    // integrate over angle
    
    for (unsigned i = 0; i < number_of_cells_; ++i)
    {
        for (unsigned j = 0; j < number_of_velocities_; ++j)
        {
            double sum = 0;
            
            for (unsigned k = 0; k < number_of_angles_ - 1; ++k)
            {
                sum += (angle_(k + 1) - angle_(k)) * (density_(k + 1) + density_(k)) / 2;
            }
            
            mean_density_(i, j) = sum;
        }
    }
    
    // integrate over velocity

    for (unsigned i = 0; i < number_of_cells_; ++i)
    {
        double sum = 0;
        
        for (unsigned j = 0; j < number_of_velocities_ - 1; ++j)
        {
            double vj = velocity(j);
            double vjp1 = velocity(j + 1);
            
            sum += (vjp1 - vj) * (vj * vj * mean_density_(i, j) + vjp1 * vjp1 * mean_density_(i, j + 1)) / 2;
        }
        
        charge_density_(i) = q_ * sum;
    }
}

void FD_Vlasov::
calculate_electric_field()
{
    fill_small_matrix();
    
    small_solver_->NumericFactorization();
    small_solver_->Solve();
    
    for (unsigned i = 0; i < number_of_points_; ++i)
    {
        electric_field_(i) = (*small_lhs)[i];
    }
}

void FD_Vlasov::
calculate_density()
{
    fill_matrix();
    
    solver_->NumericFactorization();
    solver_->Solve();
    
    for (int n = 0; n < number_of_elements_; ++n)
    {
        density_(n) = (*lhs_)[n];
    }
}

void FD_Vlasov::
solve()
{
    for (int i = 0; i < number_of_time_steps(); ++i)
    {
        if (i % dump_number_ == 0)
        {
            double t = time_step() * i;
            
            dump_xml(i, t);
        }
        
        calculate_charge_density();
        calculate_electric_field();
        calculate_density();
    }
}

int FD_Vlasov::
check_point(int p)
{
    return (p + number_of_points_) % number_of_points_;
}

int FD_Vlasov::
check_velocity(int g)
{
    if (g >= 0 && g < number_of_velocities_)
    {
        return g;
    }
    else if (g == -1)
    {
        return 1;
    }
    else if (g == number_of_velocities_)
    {
        return number_of_velocities_ - 2;
    }
    else
    {
        Check(false, "group not found");
        
        return -1;
    }
}

int FD_Vlasov::
check_angle(int o)
{
    return (o + number_of_angles_) % number_of_angles_;
}

double FD_Vlasov::
spatial_constant(int g, int o)
{
    return velocity_(g) * cos(angle_(o));
}

double FD_Vlasov::
velocity_constant(int p, int g, int o)
{
    return electric_field_(p) * cos(angle_(o));
}

double FD_Vlasov::
angle_constant(int p, int g, int o)
{
    return 1 / velocity_(g) * (-electric_field_(p) - velocity_(g) * magnetic_field_(p) * cos(2 * angle_(o)));
}

void FD_Vlasov::
dump_xml(int i, double t)
{
}
