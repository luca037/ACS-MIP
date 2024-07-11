/* Esempio: 
 * vengono forniti due programmi in input e vengono risolti entrambi.
 * Un esempio di output si trova in mip1_mip3_output.txt. */

#include <ilcplex/cplex.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>

// Frees up pointer *ptr and sets it to NULL.
void free_and_null(char** ptr) {
    if (*ptr != NULL) {
        free(*ptr);
        *ptr = NULL;
    }
}

int optimze_and_print_results(CPXENVptr env, CPXLPptr lp, double *x) {
    int status = 0;
    int solstat;

    double objval;
    
    int cur_numcols;

    // Optimize lp and obtain solution.
    status = CPXmipopt(env, lp);
    if (status) {
        fprintf(stderr, "Failed to optimize lp.\n");
        return status;
    }

    // Get solution status.
    solstat = CPXgetstat(env, lp);
    printf("Solution status %d.\n", status);

    // Print objective value.
    status = CPXgetobjval(env, lp, &objval);
    if (status) {
        fprintf(stderr, "Failed to obtain objective value (lp).\n");
        return status;
    }
    printf("Objective value: %.10g\n", objval);

    // Print variabile values.
    cur_numcols = CPXgetnumcols(env, lp); // Get number of variabiles.
    x = (double*) malloc(cur_numcols * sizeof(double));
    if (x == NULL) {
        fprintf(stderr, "No memory for solution values for lp.\n");
        return status;
    }

    status = CPXgetx(env, lp, x, 0, cur_numcols - 1);
    if (status) {
        fprintf(stderr, "Failed to obtain solution.\n");
        return status;
    }

    // Print x[i] values.
    for (int i = 0; i < cur_numcols; i++) {
        printf("x[%d] = %17.10g\n", i, x[i]);
    }
    
    return status;
};

int main(int argc, char* argv[]) {
    int solstat; // Solution status.
    //double objval; // Objective value.
    double *x_mip1 = NULL, *x_mip2 = NULL; // Variabiles value.

    CPXENVptr env = NULL;
    CPXLPptr mip1 = NULL, mip2 = NULL;
    int status;
    int i;
    int cur_numcols;

    // Check command line arguments.
    if (argc != 3) {
        fprintf(stderr, "Missing input file.\nUsage: %s <file1> <file2>\n", argv[0]);
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
    status = CPXsetintparam(env, CPXPARAM_ScreenOutput, CPX_ON);
    if (status) {
        fprintf(stderr, "Failure to turn on screen indicator, error %d\n", status);
        goto TERMINATE;
    }

// ### MIP1 ###
    printf("### Optimize %s ###\n", argv[1]);
    // Create the first problem, using the filnemane as the problem name.
    mip1 = CPXcreateprob(env, &status, argv[1]);
    if (mip1 == NULL) {
        fprintf(stderr, "Failedto create MIP1.\n");
        goto TERMINATE;
    }

    // Read file and copy the data into the created mip1.
    status = CPXreadcopyprob(env, mip1, argv[1], NULL);
    if (status) {
        fprintf(stderr, "Failed to read and copy the problem data (mip1).\n");
        goto TERMINATE;
    }

    // ### Optimize mip1 and print results ###
    status = optimze_and_print_results(env, mip1, x_mip1);
    if (status) {
        fprintf(stderr, "Failed to optimize mip1.\n");
        goto TERMINATE;
    }

// ### MIP1 ###
    printf("\n### Optimize %s ###\n", argv[2]);
    // Create the second problem, using the filnemane as the problem name.
    mip2 = CPXcreateprob(env, &status, argv[2]);
    if (mip1 == NULL) {
        fprintf(stderr, "Failedto create MIP2.\n");
        goto TERMINATE;
    }

    // Read file and copy the data into the created mip2.
    status = CPXreadcopyprob(env, mip2, argv[2], NULL);
    if (status) {
        fprintf(stderr, "Failed to read and copy the problem data (mip2).\n");
        goto TERMINATE;
    }

    // ### Optimize mip2 and print results ###
    status = optimze_and_print_results(env, mip2, x_mip2);
    if (status) {
        fprintf(stderr, "Failed to optimize mip2.\n");
        goto TERMINATE;
    }


TERMINATE:

    // Free up solution.
    free_and_null((char**) &x_mip1);
    free_and_null((char**) &x_mip2);

    // Free up the problem allocated by CPXcreateprob, if necessary.
    if (mip1 != NULL) {
        status = CPXfreeprob(env, &mip1);
        if (status) {
            fprintf(stderr, "CPXfreeprob failed, error code %d.\n", status);
        }
    }

    if (mip2 != NULL) {
        status = CPXfreeprob(env, &mip2);
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
