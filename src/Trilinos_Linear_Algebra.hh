#ifndef Trilinos_Linear_Algebra_hh
#define Trilinos_Linear_Algebra_hh

#include <vector>

using std::vector;

class Trilinos_Linear_Algebra
{
private:
    
public:

    Trilinos_Linear_Algebra();
    
    void epetra_solve(vector<double> &a_data,
                      vector<double> &b_data,
                      vector<double> &x_data,
                      unsigned number_of_elements);
  
    void amesos_dense_solve(vector<double> &a_data,
                            vector<double> &b_data,
                            vector<double> &x_data,
                            unsigned number_of_elements);
  
    void aztec_dense_solve(vector<double> &a_data,
                           vector<double> &b_data,
                           vector<double> &x_data,
                           unsigned number_of_elements);

};

#endif
