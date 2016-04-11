#include "FD_Vlasov.hh"

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "pugixml.hpp"

#include "Check.hh"
#include "XML_Child_Value.hh"

using namespace std;

/*
  Solves the Vlasov equation in 1D
  input_path: filepath for xml input file
*/
FD_Vlasov::
FD_Vlasov(string input_path):
    input_path_(input_path)
{
    parse_xml();
    initialize_trilinos();
    check_class_invariants();
}

/* 
   Check sizes of objects
*/
void FD_Vlasov::
check_class_invariants()
{
    Check(number_of_elements_ == number_of_points_ * number_of_velocities_ * number_of_angles_, "number_of_elements");

    Check(electric_field_.size() == number_of_points_, "electric_field size");
    Check(magnetic_field_.size() == number_of_points_, "magnetic_field size");
    Check(charge_density_.size() == number_of_points_, "charge_density size");

    Check(position_.size() == number_of_points_, "position size");
    Check(velocity_.size() == number_of_velocities_, "velocity size");
    Check(angle_.size() == number_of_angles_, "angle size");

    Check(density_.size() == number_of_elements_, "density size");
    Check(mean_density_.size() == number_of_points_ * number_of_velocities_, "mean_density size");

    Check(number_of_entries_per_row_.size() == number_of_elements_, "number_of_entries_per_row size");
    Check(charge_number_of_entries_per_row_.size() == number_of_points_, "charge_number_of_entries_per_row size");
}

/*
  Get problem data from xml file
*/
void FD_Vlasov::
parse_xml()
{
    output_path_ = input_path_.substr(0, input_path_.find_last_of(".")) + "-out.xml";
    output_file_.append_child("output");
    
    pugi::xml_document doc;
    
    if (!doc.load_file(input_path_.c_str()))
    {
        Check(false, "could not open document");
    }
    
    pugi::xml_node input = doc.child("input");
    
    // discretization information

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
    for (int i = 0; i < number_of_points_; ++i)
    {
        position_[i] = i * point_distance_;
    }
    
    velocity_.resize(number_of_velocities_);
    for (int i = 0; i < number_of_velocities_; ++i)
    {
        velocity_[i] = (i + 1) * velocity_distance_;
    }
    
    angle_.resize(number_of_angles_);
    for (int i = 0; i < number_of_angles_; ++i)
    {
        angle_[i] = i * angle_distance_;
    }
    
    // problem data

    pugi::xml_node data = input.child("data");
    pugi::xml_node magnetic_field = data.child("magnetic_field");
    pugi::xml_node initial_density = data.child("initial_density");
    
    electric_field_.resize(number_of_points_, 0);

    string magnetic_field_type = child_value<string>(magnetic_field, "type");
    if (magnetic_field_type == "constant")
    {
        double constant_magnetic_field = child_value<double>(magnetic_field, "value");
        magnetic_field_.resize(number_of_points_, constant_magnetic_field);
    }
    else
    {
        Check(false, "magnetic_field type not found");
    }

    charge_density_.resize(number_of_points_, 0);

    string initial_density_type = child_value<string>(initial_density, "type");
    if (initial_density_type == "constant")
    {
        double constant_density = child_value<double>(initial_density, "value");
        density_.resize(number_of_elements_, constant_density);
    }
    else
    {
        Check(false, "initial_density type not found");
    }
    
    mean_density_.resize(number_of_points_ * number_of_velocities_, 0);
}

/*
  Initialize Epetra and Aztec data
*/
void FD_Vlasov::
initialize_trilinos()
{
    //
    // initialize matrix for density solves
    //

    // periodic boundary conditions for angle and space mean that the only variation in
    // the number of entries per row is due to the velocity dirichlet boundary conditions
    // f_j = 0 for j = -1, number_of_velocities
    number_of_entries_per_row_.resize(number_of_elements_, 7);
    for (int i = 0; i < number_of_points_; ++i)
    {
        for (int k = 0; k < number_of_angles_; ++k)
        {
            int njm = subscript_to_index(i, 0, k);
            int njp = subscript_to_index(i, number_of_velocities_ - 1, k);

            number_of_entries_per_row_[njm] -= 1;
            number_of_entries_per_row_[njp] -= 1;
        }
    }
    
    // initialize communications classes

    comm_ = unique_ptr<Epetra_SerialComm> (new Epetra_SerialComm);
    map_ = unique_ptr<Epetra_Map> ( new Epetra_Map(number_of_elements_, index_base_, *comm_));

    // initialize and fill matrix

    matrix_ = unique_ptr<Epetra_CrsMatrix> (new Epetra_CrsMatrix(Copy, *map_, &number_of_entries_per_row_[0], true));
    lhs_ = unique_ptr<Epetra_Vector> (new Epetra_Vector(*map_));
    rhs_ = unique_ptr<Epetra_Vector> (new Epetra_Vector(*map_));
    fill_matrix();
    
    problem_ = unique_ptr<Epetra_LinearProblem> (new Epetra_LinearProblem(matrix_.get(), lhs_.get(), rhs_.get()));
    solver_ = unique_ptr<AztecOO> (new AztecOO(*problem_));

    solver_->SetAztecOption(AZ_precond, AZ_Jacobi);
    solver_->SetAztecOption(AZ_poly_ord, 3);
    solver_->SetAztecOption(AZ_solver, AZ_gmres);
    solver_->SetAztecOption(AZ_kspace, 1000);

    //
    // initialize matrix for charge density solves
    //
    
    charge_number_of_entries_per_row_.resize(number_of_points_, 3);

    // initialize communications classes

    charge_comm_ = unique_ptr<Epetra_SerialComm> (new Epetra_SerialComm);
    charge_map_ = unique_ptr<Epetra_Map> ( new Epetra_Map(number_of_elements_, index_base_, *comm_));

    // initialize and fill matrix
    
    charge_matrix_ = unique_ptr<Epetra_CrsMatrix> (new Epetra_CrsMatrix(Copy, *charge_map_, &charge_number_of_entries_per_row_[0], true));
    charge_lhs_ = unique_ptr<Epetra_Vector> (new Epetra_Vector(*charge_map_));
    charge_rhs_ = unique_ptr<Epetra_Vector> (new Epetra_Vector(*charge_map_));
    fill_charge_matrix();
    
    charge_problem_ = unique_ptr<Epetra_LinearProblem> (new Epetra_LinearProblem(charge_matrix_.get(), charge_lhs_.get(), charge_rhs_.get()));
    charge_solver_ = unique_ptr<AztecOO> (new AztecOO(*charge_problem_));

    charge_solver_->SetAztecOption(AZ_precond, AZ_Jacobi);
    charge_solver_->SetAztecOption(AZ_poly_ord, 3);
    charge_solver_->SetAztecOption(AZ_solver, AZ_gmres);
    charge_solver_->SetAztecOption(AZ_kspace, 1000);
}

/*
  Fill density matrix
*/
void FD_Vlasov::
fill_matrix()
{
    vector<int> subscript(3);
    
    for (int l = 0; l < number_of_elements_; ++l) // local index
    {
        int n = l; // global index
        
        index_to_subscript(l, subscript);
        
        int i = subscript[0];
        int j = subscript[1];
        int k = subscript[2];
        
        int im1 = check_point(i - 1);
        int ip1 = check_point(i + 1);
        int jm1 = check_velocity(j - 1);
        int jp1 = check_velocity(j + 1);
        int km1 = check_angle(k - 1);
        int kp1 = check_angle(k + 1);

        int nim = subscript_to_index(im1, j, k);
        int nip = subscript_to_index(ip1, j, k);
        int njm = subscript_to_index(i, jm1, k);
        int njp = subscript_to_index(i, jp1, k);
        int nkm = subscript_to_index(i, j, km1);
        int nkp = subscript_to_index(i, j, kp1);
        
        vector<int> column_indices;
        vector<double> fill_vector;
        double rhs_sum = 0;
        double value;

        column_indices.push_back(nim);
        fill_vector.push_back(-spatial_constant(j, k) / (2 * point_distance_));
        rhs_sum += fill_vector.back() * density_[nim];
        
        if (j != 0)
        {
            column_indices.push_back(njm);
            fill_vector.push_back(-velocity_constant(i, jm1, k) / (2 * velocity_distance_));
            rhs_sum -= fill_vector.back() * density_[njm];
        }
        
        column_indices.push_back(nkm);
        fill_vector.push_back(-angle_constant(i, j, km1) / (2 * angle_distance_));
        rhs_sum -= fill_vector.back() * density_[nkm];
        
        column_indices.push_back(n);
        fill_vector.push_back(2 / time_step_);
        rhs_sum -= fill_vector.back() * density_[n];
        
        column_indices.push_back(nkp);
        fill_vector.push_back(angle_constant(i, j, kp1) / (2 * angle_distance_));
        rhs_sum -= fill_vector.back() * density_[nkp];
        
        if (j != number_of_velocities_ - 1)
        {
            column_indices.push_back(njp);
            fill_vector.push_back(velocity_constant(i, jp1, k) / (2 * velocity_distance_));
            rhs_sum -= fill_vector.back() * density_[njp];
        }
        
        column_indices.push_back(nip);
        fill_vector.push_back(spatial_constant(j, k) / (2 * point_distance_));
        rhs_sum -= fill_vector.back() * density_[nip];

        Check(column_indices.size() == number_of_entries_per_row_[l], "column_indices size");
        Check(fill_vector.size() == number_of_entries_per_row_[l], "fill_vector size");
        
        // insert values into the matrix for the lhs
        if (matrix_->Filled())
        {
            matrix_->ReplaceGlobalValues(n, number_of_entries_per_row_[l], &fill_vector[0], &column_indices[0]);
        }
        else
        {
            matrix_->InsertGlobalValues(n, number_of_entries_per_row_[l], &fill_vector[0], &column_indices[0]);
        }
        
        // insert values into the rhs
        int num_entries = 1;
        vector<int> global_index = {n};
        vector<double> fill_value = {sum};
        rhs_->ReplaceGlobalValues(num_entries, &fill_value[0], &global_index[0]);
    }

    if (! matrix_->Filled())
    {
        matrix_->FillComplete();
    }
}

/* 
   Fill charge density matrix
*/
void FD_Vlasov::
fill_charge_matrix()
{
    for (int l = 0; l < number_of_points_; ++l)
    {
        int i = l;
        
        vector<int> column_indices(charge_number_of_entries_per_row_[l], 0);
        vector<double> fill_vector(charge_number_of_entries_per_row_[l], 0);
        
        int im1 = check_point(i - 1);
        int ip1 = check_point(i + 1);

        column_indices[0] = im1;
        column_indices[1] = i;
        column_indices[2] = ip1;
        
        fill_vector[0] = -electric_field_[im1] / (2 * point_distance_);
        fill_vector[1] = 0;
        fill_vector[2] = electric_field_[ip1] / (2 * point_distance_);
        
        if (matrix_->Filled())
        {
            charge_matrix_->ReplaceGlobalValues(i, charge_number_of_entries_per_row_[l], &fill_vector[0], &column_indices[0]);
        }
        else
        {
            charge_matrix_->InsertGlobalValues(i, charge_number_of_entries_per_row_[l], &fill_vector[0], &column_indices[0]);
        }
        
        int num_entries = 1;
        vector<int> global_index = {i};
        vector<double> fill_value = {charge_density_[i]};
        
        charge_rhs_->ReplaceGlobalValues(num_entries, &fill_value[0], &global_index[0]);
    }
    
    if (! charge_matrix_->Filled())
    {
        charge_matrix_->FillComplete();
    }
}

/* 
   Calculate charge density from particle density
*/
void FD_Vlasov::
calculate_charge_density()
{
    // integrate over angle
    
    for (int i = 0; i < number_of_points_; ++i)
    {
        for (int j = 0; j < number_of_velocities_; ++j)
        {
            double sum = 0;
            
            for (int k = 0; k < number_of_angles_ - 1; ++k)
            {
                int n = subscript_to_index(i, j, k);
                int nkp = subscript_to_index(i, j, k + 1);
                
                sum += (angle_[k + 1] - angle_[k]) * (density_[n] + density_[nkp]) / 2;
            }
     
            int n = j + number_of_points_ * i;
            
            mean_density_[n] = sum;
        }
    }
    
    // integrate over velocity

    for (int i = 0; i < number_of_points_; ++i)
    {
        double sum = 0;
        
        for (int j = 0; j < number_of_velocities_ - 1; ++j)
        {
            int n = j + number_of_points_ * i;
            int njp = j + 1 + number_of_points_ * i;

            double vj = velocity_[j];
            double vjp1 = velocity_[j + 1];

            sum += (vjp1 - vj) * (vj * vj * mean_density_[n] + vjp1 * vjp1 * mean_density_[njp]) / 2;
        }
        
        charge_density_[i] = q_ * sum;
    }
}

/*
  Calculate electric field from charge density
*/
void FD_Vlasov::
calculate_electric_field()
{
    fill_charge_matrix();
    
    charge_solver_->Iterate(max_iterations_, tolerance_);

    charge_lhs_->ExtractCopy(&electric_field_[0]);
}

/*
  Calculate particle density after time step
*/
void FD_Vlasov::
calculate_density()
{
    fill_matrix();
    
    solver_->Iterate(max_iterations_, tolerance_);
    
    lhs_->ExtractCopy(&density_[0]);
}

/* 
   Run problem for all time values
*/
void FD_Vlasov::
solve()
{
    for (int i = 0; i < number_of_time_steps_; ++i)
    {
        if (i % dump_number_ == 0)
        {
            double t = time_step_ * i;
            
            dump_xml(i, t);
        }
        
        calculate_charge_density();
        calculate_electric_field();
        calculate_density();
    }
}

/*
  Ensure that point index is inside problem
  Use periodic boundaries
*/
int FD_Vlasov::
check_point(int p)
{
    return (p + number_of_points_) % number_of_points_;
}

/* 
   Check that velocity index is inside problem
   Use reflecting boundaries
*/
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
        Check(false, "velocity not found");
        
        return -1;
    }
}

/* 
   Check that angle index is inside problem
   Use periodic boundaries
*/
int FD_Vlasov::
check_angle(int o)
{
    return (o + number_of_angles_) % number_of_angles_;
}

/* 
   Constant in front of spatial derivative terms
*/
double FD_Vlasov::
spatial_constant(int g, int o)
{
    return velocity_[g] * cos(angle_[o]);
}

/* 
   Constant in front of velocity derivative terms
*/
double FD_Vlasov::
velocity_constant(int p, int g, int o)
{
    return electric_field_[p] * cos(angle_[o]);
}

/*
  Constant in front of angular derivative terms
*/
double FD_Vlasov::
angle_constant(int p, int g, int o)
{
    return 1 / velocity_[g] * (-electric_field_[p] - velocity_[g] * magnetic_field_[p] * cos(2 * angle_[o]));
}

/* 
   Write current data to XML output file
*/
void FD_Vlasov::
dump_xml(int i, double t)
{
}
