#include <iostream>
#include <string>

#ifdef EPETRA_MPI
#  include "mpi.h"
#endif

#include "FD_Vlasov.hh"

using namespace std;

int main(int argc, char** argv)
{
#ifdef EPETRA_MPI
    MPI_Init(&argc, &argv);
#endif

    if (argc != 2)
    {
        cerr << "usage: phoenix [input.xml]" << endl;
        return 1;
    }

    string filename = argv[1];
    
    FD_Vlasov fd_vlasov(filename);

#ifdef EPETRA_MPI
    MPI_Finalize();
#endif
}

