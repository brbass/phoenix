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
    solve();
    save_xml();
}

/* 
   Check sizes of objects
*/
void FD_Vlasov::
check_class_invariants()
{
    Assert(number_of_elements_ == number_of_points_ * number_of_velocities_ * number_of_angles_);

    Assert(position_.size() == number_of_points_);
    Assert(velocity_.size() == number_of_velocities_);
    Assert(angle_.size() == number_of_angles_);

    Assert(electric_field_x_.size() == number_of_points_);
    Assert(electric_field_y_.size() == number_of_points_);
    Assert(magnetic_field_z_.size() == number_of_points_);
    Assert(charge_density_.size() == number_of_points_);

    Assert(density_.size() == number_of_elements_);
    Assert(mean_density_.size() == number_of_points_ * number_of_velocities_);

    Assert(number_of_entries_per_row_.size() == number_of_elements_);
    Assert(charge_number_of_entries_per_row_.size() == number_of_points_);
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
        AssertMsg(false, "could not open document");
    }
    
    pugi::xml_node input = doc.child("input");
    
    // discretization information

    pugi::xml_node solver_node = input.child("solver");
    pugi::xml_node point_node = input.child("point_discretization");
    pugi::xml_node velocity_node = input.child("velocity_discretization");
    pugi::xml_node angle_node = input.child("angle_discretization");
    pugi::xml_node time_node = input.child("time_discretization");
    
    max_iterations_ = child_value<int>(solver_node, "max_iterations");
    tolerance_ = child_value<double>(solver_node, "tolerance");
    solver_print_ = child_value<int>(solver_node, "print");
    kspace_ = child_value<int>(solver_node, "kspace");
    poly_ord_ = child_value<int>(solver_node, "poly_ord");

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
        angle_[i] = static_cast<double>(i) * angle_distance_;
    }
    
    // problem data

    pugi::xml_node data = input.child("data");
    pugi::xml_node magnetic_field = data.child("magnetic_field");
    pugi::xml_node electric_field = data.child("electric_field");
    pugi::xml_node initial_density = data.child("initial_density");
    
    electric_field_x_.assign(number_of_points_, 0);

    electric_field_y_.resize(number_of_points_, 0);
    magnetic_field_z_.resize(number_of_points_, 0);

    string electric_field_type = child_value<string>(electric_field, "type");
    if (electric_field_type == "constant")
    {
        double constant_electric_field = child_value<double>(electric_field, "value");
        electric_field_y_.assign(number_of_points_, constant_electric_field);
    }
    else
    {
        AssertMsg(false, "electric_field type not found");
    }

    string magnetic_field_type = child_value<string>(magnetic_field, "type");
    if (magnetic_field_type == "constant")
    {
        double constant_magnetic_field = child_value<double>(magnetic_field, "value");
        magnetic_field_z_.assign(number_of_points_, constant_magnetic_field);
    }
    else
    {
        AssertMsg(false, "magnetic_field type not found");
    }

    charge_density_.resize(number_of_points_, 0);

    string initial_density_type = child_value<string>(initial_density, "type");
    if (initial_density_type == "constant")
    {
        double constant_density = child_value<double>(initial_density, "value");
        density_.resize(number_of_elements_, constant_density);
    }
    else if (initial_density_type == "monodirectional")
    {
        double monodirectional_density = child_value<double>(initial_density, "value");
        int k = child_value<int>(initial_density, "direction");
        
        density_.assign(number_of_elements_, 0);
        for (unsigned i = 0; i < number_of_points_; ++i)
        {
            for (unsigned j = 0; j < number_of_velocities_; ++j)
            {
                int n = subscript_to_index(i, j, k);
                density_[n] = monodirectional_density;
            }
        }
    }
    else if (initial_density_type == "monoenergetic")
    {
        double monoenergetic_density = child_value<double>(initial_density, "value");
        int j = child_value<int>(initial_density, "velocity");
        
        density_.assign(number_of_elements_, 0);
        for (unsigned i = 0; i < number_of_points_; ++i)
        {
            for (unsigned k = 0; k < number_of_angles_; ++k)
            {
                int n = subscript_to_index(i, j, k);
                density_[n] = monoenergetic_density;
            }
        }
    }
    else
    {
        AssertMsg(false, "initial_density type not found");
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
    // the number of entries per row is due to the velocity reflecting boundary conditions
    // f_j = 0 for j = -1, number_of_velocities
    number_of_entries_per_row_.resize(number_of_elements_, 7);
    for (int i = 0; i < number_of_points_; ++i)
    {
        for (int k = 0; k < number_of_angles_; ++k)
        {
            int njm = subscript_to_index(i, 0, k);
            int njp = subscript_to_index(i, number_of_velocities_ - 1, k);

            number_of_entries_per_row_[njm] = 5;
            number_of_entries_per_row_[njp] = 5;
        }
    }
    
    // initialize communications classes

#ifdef EPETRA_MPI
    comm_ = unique_ptr<Epetra_MpiComm> (new Epetra_MpiComm(MPI_COMM_WORLD));
#else
    comm_ = unique_ptr<Epetra_SerialComm> (new Epetra_SerialComm);
#endif
    map_ = unique_ptr<Epetra_Map> ( new Epetra_Map(number_of_elements_, index_base_, *comm_));

    // initialize and fill matrix

    matrix_ = unique_ptr<Epetra_CrsMatrix> (new Epetra_CrsMatrix(Copy, *map_, &number_of_entries_per_row_[0], true));
    lhs_ = unique_ptr<Epetra_Vector> (new Epetra_Vector(*map_));
    rhs_ = unique_ptr<Epetra_Vector> (new Epetra_Vector(*map_));
    fill_matrix();
    
    problem_ = unique_ptr<Epetra_LinearProblem> (new Epetra_LinearProblem(matrix_.get(), lhs_.get(), rhs_.get()));
    solver_ = unique_ptr<AztecOO> (new AztecOO(*problem_));

    solver_->SetAztecOption(AZ_precond, AZ_Jacobi);
    solver_->SetAztecOption(AZ_poly_ord, poly_ord_);
    solver_->SetAztecOption(AZ_solver, AZ_gmres);
    solver_->SetAztecOption(AZ_kspace, kspace_);
    if (solver_print_)
    {
        solver_->SetAztecOption(AZ_output, AZ_all);
    }
    else
    {
        solver_->SetAztecOption(AZ_output, AZ_none);
    }

    //
    // initialize matrix for charge density solves
    //
    
    charge_number_of_entries_per_row_.resize(number_of_points_, 3);

    // initialize communications classes

#ifdef EPETRA_MPI
    charge_comm_ = unique_ptr<Epetra_MpiComm> (new Epetra_MpiComm(MPI_COMM_WORLD));
#else
    charge_comm_ = unique_ptr<Epetra_SerialComm> (new Epetra_SerialComm);
#endif
    charge_map_ = unique_ptr<Epetra_Map> ( new Epetra_Map(number_of_points_, index_base_, *charge_comm_));

    // initialize and fill matrix
    
    charge_matrix_ = unique_ptr<Epetra_CrsMatrix> (new Epetra_CrsMatrix(Copy, *charge_map_, &charge_number_of_entries_per_row_[0], true));
    charge_lhs_ = unique_ptr<Epetra_Vector> (new Epetra_Vector(*charge_map_));
    charge_rhs_ = unique_ptr<Epetra_Vector> (new Epetra_Vector(*charge_map_));
    fill_charge_matrix();
    
    charge_problem_ = unique_ptr<Epetra_LinearProblem> (new Epetra_LinearProblem(charge_matrix_.get(), charge_lhs_.get(), charge_rhs_.get()));
    charge_solver_ = unique_ptr<AztecOO> (new AztecOO(*charge_problem_));

    charge_solver_->SetAztecOption(AZ_precond, AZ_Jacobi);
    charge_solver_->SetAztecOption(AZ_poly_ord, poly_ord_);
    charge_solver_->SetAztecOption(AZ_solver, AZ_gmres);
    charge_solver_->SetAztecOption(AZ_kspace, kspace_);
    if (solver_print_)
    {
        charge_solver_->SetAztecOption(AZ_output, AZ_all);
    }
    else
    {
        charge_solver_->SetAztecOption(AZ_output, AZ_none);
    }
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

        // time derivative
        
        column_indices.push_back(n);
        value = 2 / time_step_;
        fill_vector.push_back(value);
        rhs_sum += value * density_[n];

        // spatial derivative
        
        column_indices.push_back(nip);
        value = spatial_constant(j, k) / (2 * point_distance_);
        fill_vector.push_back(value);
        rhs_sum -= value * density_[nip];
            
        column_indices.push_back(nim);
        value = -spatial_constant(j, k) / (2 * point_distance_);
        fill_vector.push_back(value);
        rhs_sum -= value * density_[nim];

        // velocity derivative

        if (j == 0 || j == number_of_velocities_ - 1)
        {
            
        }
        else 
        {
            column_indices.push_back(njp);
            value = velocity_constant(i, jp1, k) / (2 * velocity_distance_);
            fill_vector.push_back(value);
            rhs_sum -= value * density_[njp];
        
            column_indices.push_back(njm);
            value = -velocity_constant(i, jm1, k) / (2 * velocity_distance_);
            fill_vector.push_back(value);
            rhs_sum -= value * density_[njm];
        }
        // angular derivative
        
        column_indices.push_back(nkm);
        value = -angle_constant(i, j, km1) / (2 * angle_distance_);
        fill_vector.push_back(value);
        rhs_sum -= value * density_[nkm];
        
        column_indices.push_back(nkp);
        value = angle_constant(i, j, kp1) / (2 * angle_distance_);
        fill_vector.push_back(value);
        rhs_sum -= value * density_[nkp];
        
        Check(column_indices.size() == number_of_entries_per_row_[l]);
        Check(fill_vector.size() == number_of_entries_per_row_[l]);
        
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
        vector<double> fill_value = {rhs_sum};
        rhs_->ReplaceGlobalValues(num_entries, &fill_value[0], &global_index[0]);
    }

    if (! matrix_->Filled())
    {
        matrix_->FillComplete();
    }
    
    // cout << *matrix_ << endl;
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
        
        fill_vector[0] = -electric_field_x_[im1] / (2 * point_distance_);
        fill_vector[1] = 0;
        fill_vector[2] = electric_field_x_[ip1] / (2 * point_distance_);
        
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

    // cout << *charge_matrix_ << endl;
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
            
            for (int k = 0; k < number_of_angles_; ++k)
            {
                int k1 = check_angle(k+1);
                int n = subscript_to_index(i, j, k);
                int nkp = subscript_to_index(i, j, k1);
                
                sum += angle_distance_ * (density_[n] + density_[nkp]) / 2;
            }
            
            int n = j + number_of_velocities_ * i;
            
            mean_density_[n] = sum;
        }
    }
    
    // integrate over velocity

    for (int i = 0; i < number_of_points_; ++i)
    {
        double sum = 0;
        
        for (int j = 0; j < number_of_velocities_ - 1; ++j)
        {
            int n = j + number_of_velocities_ * i;
            int njp = j + 1 + number_of_velocities_ * i;

            sum += velocity_distance_ * (velocity_[j] * mean_density_[n] + velocity_[j + 1] * mean_density_[njp]) / 2;
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

    charge_lhs_->ExtractCopy(&electric_field_x_[0]);
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
    calculate_charge_density();
    calculate_electric_field();

    initial_xml();

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
    if (g < 0)
    {
        g = -g;
    }
    else if (g >= number_of_velocities_)
    {
        g = number_of_velocities_ - 1 - (g - number_of_velocities_ - 1);
    }
    else
    {
        return g;
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
    return qm_ * (electric_field_x_[p] * cos(angle_[o]) + electric_field_y_[p] * sin(angle_[o]));
}

/*
  Constant in front of angular derivative terms: MAY HAVE ERROR
*/
double FD_Vlasov::
angle_constant(int p, int g, int o)
{
    return qm_ / velocity_[g] * (-electric_field_x_[p] * sin(angle_[o]) + electric_field_y_[p] * cos(angle_[o])) - qm_ * magnetic_field_z_[p];
}

/* 
   Write current data to XML output file
*/
void FD_Vlasov::
initial_xml()
{
    pugi::xml_node output = output_file_.child("output");
    string input_file = input_path_.substr(input_path_.find_last_of("/\\") + 1);
    append_child(output, input_file, "input_file");

    pugi::xml_node data = output.append_child("data");
    
    append_child(data, number_of_points_, "number_of_points");
    append_child(data, number_of_velocities_, "number_of_velocities");
    append_child(data, number_of_angles_, "number_of_angles");
    append_child(data, number_of_elements_, "number_of_elements");
    append_child(data, number_of_time_steps_, "number_of_time_steps");
    append_child(data, dump_number_, "dump_number");
    append_child(data, point_distance_, "point_distance");
    append_child(data, velocity_distance_, "velocity_distance");
    append_child(data, angle_distance_, "angle_distance");
    append_child(data, time_step_, "time_step");
    append_child(data, q_, "q");
    append_child(data, qm_, "qm");
    append_child(data, position_, "position");
    append_child(data, velocity_, "velocity");
    append_child(data, angle_, "angle");
    append_child(data, number_of_entries_per_row_, "number_of_entries_per_row");
    append_child(data, charge_number_of_entries_per_row_, "charge_number_of_entries_per_row");
    append_child(data, max_iterations_, "max_iterations");
    append_child(data, tolerance_, "tolerance");
}


/* 
   Write current data to XML output file
*/
void FD_Vlasov::
dump_xml(int i, double t)
{
    pugi::xml_node output = output_file_.child("output");
    pugi::xml_node dump = output.append_child("dump");

    pugi::xml_attribute it = dump.append_attribute("time_step");
    it.set_value(to_string(i).c_str());
    pugi::xml_attribute ti = dump.append_attribute("time");
    ti.set_value(to_string(t).c_str());

    append_child(dump, electric_field_x_, "electric_field_x");
    append_child(dump, electric_field_y_, "electric_field_y");
    append_child(dump, magnetic_field_z_, "magnetic_field_z");
    append_child(dump, charge_density_, "charge_density");
    append_child(dump, density_, "density");
    append_child(dump, mean_density_, "mean_density");
}

void FD_Vlasov::
save_xml()
{
    // pugi::xml_node output = output_.child("output");
    // pugi::xml_node final = output.append_child("final");

    // append_child(final, time_values_, "time_values");
    // append_child(final, kinetic_energy_, "kinetic_energy");
    
    output_file_.save_file(output_path_.c_str());
}
