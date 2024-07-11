/* Esempio: 
 * il modello del problema viene fornito da un file in input e 
 * viene aggiunta una nuova colonna al problema tramite la funzione 
 * add_column. */

#include <ilcplex/cplex.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>

void free_and_null(char** ptr);
int add_column(CPXENVptr env, CPXLPptr lp);

int main(int argc, char* argv[]) {
    int solstat; // Solution status.
    double objval; // Objective value.
    double* x = NULL; // Variabiles value.

    CPXENVptr env = NULL;
    CPXLPptr lp = NULL;
    int status;
    int i;
    int cur_numcols;

    // Check command line arguments.
    if (argc != 2) {
        fprintf(stderr, "Missing input file.\nUsage: %s <file>\n", argv[0]);
        goto TERMINATE;
    }

    // Inizialize the CPLEX environment.
    env = CPXopenCPLEX(&status);
    if (env == NULL) {
        char errrmsg[CPXMESSAGEBUFSIZE];
        fprintf(stderr, "Could not open CPLEX environment.\n");
        CPXgeterrorstring(env, status, errrmsg);
        fprintf(stderr, "%s", errrmsg);
        goto TERMINATE;
    }

    // Turn on output to the screen.
    //status = CPXsetintparam(env, CPXPARAM_ScreenOutput, CPX_ON);
    //if (status) {
    //    fprintf(stderr, "Failure to turn on screen indicator, error %d\n", status);
    //    goto TERMINATE;
    //}

    // Create the problem, using the filnemane as the problem name.
    lp = CPXcreateprob(env, &status, argv[1]);
    if (lp == NULL) {
        fprintf(stderr, "Failedto create LP.\n");
        goto TERMINATE;
    }

    // Read file and copy the data into the created lp.
    status = CPXreadcopyprob(env, lp, argv[1], NULL);
    if (status) {
        fprintf(stderr, "Failed to read and copy the problem data.\n");
        goto TERMINATE;
    }

    // Add new col to the problem.
    status = add_column(env, lp);
    if (status) {
        fprintf(stderr, "Failed to add new column.\n");
        goto TERMINATE;
    }

    // Optimize the problem and obtain solution.
    status = CPXmipopt(env, lp);
    if (status) {
        fprintf(stderr, "Failed to optimize MIP.\n");
        goto TERMINATE;
    }

    // Get solution status.
    solstat = CPXgetstat(env, lp);
    printf("Solution status %d.\n", status);

    // Print objective value.
    status = CPXgetobjval(env, lp, &objval);
    if (status) {
        fprintf(stderr, "Failed to obtain objective value.\n");
        goto TERMINATE;
    }
    printf("Objective value %.10g\n", objval);

    // Print variabile values.
    cur_numcols = CPXgetnumcols(env, lp); // Get number of variabiles.
    x = (double*) malloc(cur_numcols * sizeof(double));
    if (x == NULL) {
        fprintf(stderr, "No memory for solution values.\n");
        goto TERMINATE;
    }

    status = CPXgetx(env, lp, x, 0, cur_numcols - 1);
    if (status) {
        fprintf(stderr, "Failed to obtain solution.\n");
        goto TERMINATE;
    }

    // Print x[i] values.
    for (i = 0; i < cur_numcols; i++) {
        printf("x[%d] = %17.10g\n", i, x[i]);
    }

    // Write a copy of the problem to a file.
    status = CPXwriteprob(env, lp, "addCol.lp", NULL);
    if (status) {
        fprintf (stderr, "Failed to write LP to disk.\n");
    }


TERMINATE:

    // Free up solution.
    free_and_null((char**) &x);

    // Free up the problem allocated by CPXcreateprob, if necessary.
    if (lp != NULL) {
        status = CPXfreeprob(env, &lp);
        if (status) {
            fprintf(stderr, "CPXfreeprob failed, error code %d.\n", status);
        }
    }

    // Free up the CPLEX environment, if necessary.
    if (env != NULL) {
        status = CPXcloseCPLEX(&env);
        if (status) {
            char errmsg[CPXMESSAGEBUFSIZE];
            fprintf(stderr, "Could not close CPLEX environment.\n");
            CPXgeterrorstring(env, status, errmsg);
            fprintf(stderr, "%s", errmsg);
        }
    }

    return status;
}

// Add a new integer column to the problem.
// In questo caso vado ad aggiungere una nuova colonna (variabile x3) alla
// matrice dei vincoli.
int add_column(CPXENVptr env, CPXLPptr lp) {
    int status = 0;

    int ccnt  = 1; // Number of new columns being added to the constraint 
                   // matrix. 
    int nzcnt = 2; // Number of nonzero constraint coefficients to be added to 
                   // the constraint matrix.

    double lb = 0.0; // Limite inferiore variabile.
    double ub = 10.0; // Limite superiore variabile.

    int matbeg = 0; // Inizio dei coefficienti nella matrice 
                    // dei vincoli.
    int matind[] = {0, 1}; // Indici delle righe in cui la variabile ha un coeff
                            // diverso da zero.
    double matval[] = {9, 10}; // Coefficienti nei vincoli.

    status = CPXaddcols(env, lp, ccnt, nzcnt, NULL, &matbeg, matind, matval, &lb, &ub, NULL);
    if (status) {
        fprintf(stderr, "Failed to add a new column to the problem.\n");
        return status;
    }

    // Specifico il tipo di variabile aggiunta.
    char ctype = CPX_INTEGER;
    int var_index = CPXgetnumcols(env, lp) - 1; // Indice della variabile
                                                // aggiunta (x3 nel nostro
                                                // caso).
    status = CPXchgctype(env, lp, 1, &var_index, &ctype);
    if (status) {
        fprintf(stderr, "Failed to change ctype.\n");
        return status;
    }

    return status;
}

// Frees up pointer *ptr and sets it to NULL.
void free_and_null(char** ptr) {
    if (*ptr != NULL) {
        free(*ptr);
        *ptr = NULL;
    }
}
