#include "Array.hh"
#include "FD_Vlasov.hh"

FD_Vlasov::
FD_Vlasov(int number_of_points,
          int number_of_groups,
          int number_of_ordinates,
          int number_of_time_steps,
          int dump_number,
          double time_step,
          Array<double> magnetic_field,
          Array<double> angle,
          Array<double> position,
          Array<double> velocity,
          Array<double> initial_density):
    number_of_points_(number_of_points),
    number_of_groups_(number_of_groups),
    number_of_ordinates_(number_of_ordinates),
    number_of_elements_(number_of_points * number_of_groups * number_of_ordinates),
    number_of_time_steps_(number_of_time_steps),
    dump_number_(dump_number),
    magnetic_field_(magnetic_field),
    angle_(angle),
    position_(position),
    velocity_(velocity),
    angle_length_(get_length(number_of_angles,
                             angle)),
    position_length_(get_length(number_of_points, 
                                position)),
    velocity_length_(get_length(number_of_groups,
                                velocity))
    density_(initial_density)
{
    number_of_entries_per_row_.resize(number_of_elements(), 7);

    initialize_trilinos();
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
    // initialize communications classes

    comm_ = unique_ptr<Epetra_SerialComm> (new Epetra_SerialComm);
    map_ = unique_ptr<Epetra_Map> ( new Epetra_Map(number_of_elements(), index_base, *comm));

    // initialize and fill matrix

    matrix_ = unique_ptr<Epetra_CRSMatrix> (new Epetra_CrsMatrix(Copy, *map_, &number_of_entries_per_row_[0], true));
    fill_matrix();

    lhs_ = unique_ptr<Epetra_Vector> (new Epetra_Vector(*map_));
    fill_lhs();
    
    rhs_ = unique_ptr<Epetra_Vector> (new Epetra_Vector(*map_));
    fill_rhs();

    problem_ = unique_ptr<Epetra_LinearProblem> (new Epetra_LinearProblem(matrix_.get(), lhs_.get(), rhs_.get()));
    solver_ = unique_ptr<Amesos_BaseSolver> (factory_.Create(solver_type_, *problem_));

    solver_->SymbolicFactorization();
}

void FD_Vlasov::
fill_matrix()
{
    for (int n = 0; n < number_of_elements(); ++n)
    {
        int i = local_index(n, 0);
        int j = local_index(n, 1);
        int k = local_index(n, 2);

        vector<double> initial_vector(number_of_entries_per_row_[i], 0);
        
        int sum = 0;
        for (int i1 = i - 1; i < i + 2; ++i)
        {
            for (int j1 = j - 1; j1 < j + 2; ++j)
            {
                for (int k1 = k - 1; k1 < k + 2; ++k)
                {
                    column_indices[sum] = check_ordinate_index(k1) + number_of_ordinates() * (check_group_index(j1) + number_of_groups() * check_cell_index(i1));
                }
            }
        }
        
        vector<double> initial_vector(number_of_entries_per_row_[i], 0);
        matrix.InsertGlobalValues(i, number_of_entries_per_row_[i], &initial_vector[0], &column_indices[0]);
    }
}

void FD_Vlasov::
fill_lhs()
{
    
}

void FD_Vlasov::
fill_rhs();
{
    
}

void FD_Vlasov::
calculate_charge_density()
{
    
}

void FD_Vlasov::
calculate_electric_field()
{
    
}

void FD_Vlasov::
calculate_force()
{
    calculate_charge_density();
    calculate_eletric_field();

    
}

void FD_Vlasov::
step()
{
    
}

void FD_Vlasov::
solve()
{
    for (int i = 0; i < number_of_time_steps(); ++i)
    {
        calculate_force();
        
        step();
    }
}
