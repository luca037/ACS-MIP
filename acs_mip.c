#include <ilcplex/cplex.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>

/* Maximum lenght of the slack variables names (considering '\0').
 * The slack vectors Delta+ and Delta- have the same length.
 * The variables ara named as follow:
 *     Delta+ = [dp1, dp2, dp3, ..., dpi, ...,]
 *     Delta- = [dn1, dn2, dn3, ..., dni, ...,]
 * Note that 2 characters are taken by 'd' and 'p' or 'n' and 1 character
 * is taken by '\0'. */
#define MAX_SLACK_NAMES_LEN 9

/* Maximum column's name length. */
#define MAX_COLNAME_LEN 9

/* Maximum number of attempts to solve FMIP/OMIP that has no solution due to
 * infeasability. */
#define MAX_ATTEMPTS 2000

/* Maximum number of nodes explored during optimization. */
#define NODE_LIMIT 5000

/* Optimizer time limit (seconds). */
#define TIME_LIMIT 180.0

/* Maximum number of run. A run is complete when the solver has found an
 * optimal solution for FMIP and then an optimal solution for OMIP. */
#define MAX_RUN 50

/* It's used for initializing the starting vector. The value of the variable
 * is raondomly choosen from the range
 *      [max(lb, BOUND_CONSTANT), min(ub, BOUND_CONSTANT)]
 *  where lb is the lower bound of the variable and ub is its upper bound. */
#define BOUND_CONSTANT 10E6


/**
 * Print help message.
 *
 * progname: The executable's name.
 */
void print_usage(char *progname) {
   fprintf (stderr, "Usage: %s <intput> <output>\n"
                    "   input: is a file with extension \n"
                    "      MPS, SAV, or LP (lower case is allowed)\n"
                    "   output: is a file with extension csv\n"
                    "  This program uses the CPLEX MIP optimizer.\n"
                    " Exiting...\n", progname
   );
}


/**
 * Frees up pointer and set it to NULL.
 *
 * ptr: The pointer to free up.
 */
void free_and_null(char** ptr) {
    if (*ptr != NULL) {
        free(*ptr);
        *ptr = NULL;
    }
}


/**
 * Sums elements from beg to end.
 *
 * x: array with at least (end - beg + 1) elements.
 * beg: Position where the sum begins.
 * end: Position where the sum ends.
 *
 * return: The sum of the elements.
 */
double sum(double *x, int beg, int end) {
    if (beg > end) {
        return 0;
    }
    double s = x[beg];
    for (int i = beg + 1; i < end; i++) {
        s += x[i];
    }
    return s;
}


/**
 * Returns the maximum between a and b.
 */
double max(double a, double b) {
    return (a >= b)? a : b;
}

/**
 * Returns the minimum between a and b.
 */
double min(double a, double b) {
    return (a <= b)? a : b;
}

/**
 * Accesses column's name.
 *
 * env: A pointer to the CPLEX environment.
 * lp: A pointer to a CPLEX problem object.
 * index: Index of the column.
 * colname: An array where the specified column names are to be returned.
 *
 * return: Returns 0 if successful and nonzero if an error occurs.
 */
int get_colname(CPXENVptr env, CPXLPptr lp, int index, char *colname) {
    int surplus, status;
    char namestore[MAX_COLNAME_LEN];

    // Ricavo il nome della colonna.
    status = CPXgetcolname(
        env,
        lp,
        (char**) &namestore,
        namestore,
        MAX_COLNAME_LEN,
        &surplus,
        index,
        index
    );
    if (status) {
        fprintf(stderr, "Failed get column name.\n");
        return status;
    }

    strncpy(colname, namestore, MAX_COLNAME_LEN);

    return 0;
}


/**
 * Inizializes lower bounds, upper bounds and the array with indices of the
 * integer variables.
 *
 * env: A pointer to the CPLEX environment.
 * mip: A pointer to a CPLEX problem object.
 * lp: An array where the lower bounds of mip are to be returned.
 *     This array must be of length at least numcols.
 * up: An array where the upper bounds of mip are to be returned.
 *     This array must be of length at least numcols.
 * numcols: The length of the arrays lp and up.
 * int_indices: An array where the indices of integer variables are to be
 *              returned.
 *              This array must be of length at least num_int_vars.
 * num_int_vars: The length of the array int_indices.
 *
 * return: Returns 0 if successful and nonzero if an error occurs.
 */
int init_mip_bds_and_indices(
    CPXENVptr env,
    CPXLPptr mip,
    double *lb,
    double *ub,
    int numcols,
    int *int_indices,
    int num_int_vars
) {
    int i, j, status;

    // Salvo gli indici delle variabili intere di MIP.
    for (i = 0, j = 0; i < numcols; i++) {
        char type;
        status = CPXgetctype(env, mip, &type, i, i);
        if (status) {
            fprintf(stderr, "Failed to get variable %d type.", i);
            return status;
        }

        if (type == CPX_BINARY || type == CPX_INTEGER) {
            int_indices[j] = i; // Salvo l'indice della variabile.
            j += 1;
        }
    }

    // Salvo lb e up di MIP.
    status = CPXgetub(env, mip, ub, 0, numcols - 1);
    if (status) {
        fprintf(stderr, "Failed to copy upper bounds coefficients of SRC.\n");
        return status;
    }

    status = CPXgetlb(env, mip, lb, 0, numcols - 1);
    if (status) {
        fprintf(stderr, "Failed to copy upper bounds coefficients of SRC.\n");
        return status;
    }

    return 0;
}


/**
 * Restores upper and lower bounds of lp with the specified values only for the
 * variables that was fixed by variable_fixing.
 *
 * env: A pointer to the CPLEX environment.
 * lp: A pointer to a CPLEX problem object.
 * lb: An array where are stored all the lower bounds of lp.
 *     This array must be of length at least numcols.
 * ub: An array where are stored all the upper bounds of lp.
 *     This array must be of length at least numcols.
 * numcols: The length of the arrays lb and up.
 * int_indices: An array where are stored the indices of the integer variables
 *              of lp.
 *              This array must be of length at least num_int_vars.
 * fixed_indices: An array that tells if a variable was fixed.
 *                A variable with index int_indices[i] was fixed if and only if 
 *                fixed_indices[i] is a nonzero value.
 *                This array must be of length at least num_int_vars.
 * num_int_vars: The length of the arrays int_indices and fixed_indices.
 *
 * return: Returns 0 if successful and nonzero if an error occurs.
 */
int restore_bounds(
    CPXENVptr env,
    CPXLPptr lp,
    double *lb,
    double *ub,
    int numcols,
    int *int_indices,
    int *fixed_indices,
    int num_int_vars
) {
    int i, status;
    char lower = 'L', upper = 'U';

    for (i = 0; i < num_int_vars; i++) {
        if (fixed_indices[i]) { // Se è stata fissata.
            // Reset lower bound.
            status = CPXchgbds(
                env,
                lp,
                1,
                &int_indices[i],
                &lower,
                &lb[int_indices[i]]
            );
            if (status) {
                fprintf(stderr, "Failed to reset lower bound "
                                "of variable %d.\n", int_indices[i]);
                return status;
            }
            // Reset upper bound.
            status = CPXchgbds(
                env,
                lp,
                1,
                &int_indices[i],
                &upper,
                &ub[int_indices[i]]
            );
            if (status) {
                fprintf(stderr, "Failed to reset upper bound "
                                "of variable %d.\n", int_indices[i]);
                return status;
            }
        }
    }

    return 0;
}


/**
 * Generate the starting point, thus the initial vector for the first
 * FMIP problem to solve.
 *
 * initial_vector: An array where the generate random values are to be returned.
 * int_indices: An array where are stored the indices of the integer variables
 *              of lp.
 *              This array must be of length at least num_int_vars.
 * lb: An array where are stored all the lower bounds.
 * ub: An array where are stored all the upper bounds.
 * cb: The fixed bound constant.
 */
void init_initial_vector(
    double *initial_vector,
    int *int_indices,
    int num_int_vars,
    double *lb,
    double *ub,
    double cb
) {
    int i, upper, lower, val;
    
    printf("\n");
    for (i = 0; i < num_int_vars; i++) {
        upper = (int) min(ub[int_indices[i]], cb);
        lower = (int) max(lb[int_indices[i]], -cb);
        val = (int) (rand() % (upper - lower + 1)) + lower;
        printf("Col: %d, upper: %d, lower: %d, val: %d\n", int_indices[i], upper, lower, val);
        // (rand() % (upper – lower + 1)) + lower 
        initial_vector[int_indices[i]] = val;
    }
}


/**
 * Fixes some of the problem's variables to a value specified by the array
 * named x. The index of the variable that will be fixed is randomly choose
 * from the indices specifed in the array named int_indices.
 * The ammount of variables that is fixed is set by the percentage parameter.
 *
 * For example:
 *      int_indices = {1, 4, 20}, x = {100, 49, 66}, fixed_indices = {0, 0, 0}
 *      
 *      If the random generetor generates the value 2 then the variable
 *      with index int_indices[2] is fixed to the value x[int_indices[2]].
 *      Then in fixed_indices[2] will be stored 1.
 *
 *      In other words variable with index 20 will be fixed to 66 and 
 *      fixed_indices = {0, 0, 1}.
 *
 * env: A pointer to the CPLEX environment.
 * lp: A pointer to a CPLEX problem object.
 * int_indices: An array where are stored the indices of the integer variables
 *              of lp. This array must be of length at least num_int_vars.
 * fixed_indices: An array that tells if a variable was fixed.
 *                A variable with index int_indices[i] was fixed if and only if 
 *                fixed_indices[i] is a nonzero value.
 *                This array must be of length at least num_int_vars.
 * num_int_vars: The length of the arrays int_indices and fixed_indices.
 * x: An array where are stored the values to use for fixing the variables.
 *    This array must be of length at least the ammount of columns of the
 *    starting MIP problem.
 * percent: Percentage of integer variable to fix. This value must be in the
 *          range [0, 100].
 *
 * return: Returns 0 if successful and nonzero if an error occurs.
 */
int variable_fixing(
    CPXENVptr env,
    CPXLPptr lp,
    int *int_indices,
    int *fixed_indices,
    int num_int_vars,
    double *x,
    int percentage
) {
    int i, tmp, cnt, rnd, status;
    double val; // Valore a cui fisso la var scelta.
    char lu = 'B';

    // Numero di variabili da fissare.
    cnt = num_int_vars * percentage / 100;
    printf("Variables to fix: %d\n", cnt);

    // Genero una posizione random di 'int_indices'.
    rnd = rand() % num_int_vars;

    // Fisso cnt variabili intere a partire dall'indice estratto.
    for (i = 0; i < cnt; i++) {
        // Prelevo il valore.
        tmp = (rnd + i) % num_int_vars;
        val = x[int_indices[tmp]];

        // Fisso il valore della variabile.
        status = CPXchgbds(env, lp, 1, &int_indices[tmp], &lu, &val);
        if (status) {
            fprintf(stderr, 
                    "Failed to fix variable with index %d.\n", int_indices[rnd]);
            return status;
        }
        fixed_indices[tmp] = 1;
    }

    // Stampa variabili fissate.
    //printf("Variabili fissate:\n");
    //for (int j = 0; j < num_int_vars; j++) {
    //    char name[MAX_COLNAME_LEN];
    //    if (fixed_indices[j]) {
    //        get_colname(env, lp, int_indices[j], name);
    //        printf("Var fissata -> %s, valore -> %.2f\n", name, x[int_indices[j]]);
    //    }
    //}
    //printf("\n");


    //for (i = 0; i < cnt; ) {
    //    // Genero una posizione random di 'index'.
    //    rnd = rand() % num_int_vars;
    //    // Se è già stata fissata la variabile genero un altro valore.
    //    if (fixed_indices[rnd] == 1) {
    //        continue;
    //    }

    //    val = x[int_indices[rnd]];

    //    // Fisso il valore della variabile.
    //    status = CPXchgbds(env, lp, 1, &int_indices[rnd], &lu, &val);
    //    if (status) {
    //        fprintf(stderr, 
    //                "Failed to fix variable with index %d.\n", int_indices[rnd]);
    //        return status;
    //    }
    //    fixed_indices[rnd] = 1;

    //    i += 1;
    //}

    return 0;
}


/**
 * Optimizes the problem pointed by lp.
 *
 * env: A pointer to the CPLEX environment.
 * lp: A pointer to a CPLEX problem object.
 * objval: A pointer where objective optimal value is to be returned.
 * solstat: A pointer where the status code is to be returned. 
 * x: An array where the values of the variables are to be returned. 
 *    This array must be of length at least (end - beg + 1).
 * beg: Specifies the beginning of the range of the variables values to be
 *      returned.
 * end: Specifies the end of the range of the variables values to be returned.
 * verbose: Manage output verbosity. If it's set to a nonzero value then
 *          objval and solstat are printed to standard output.
 *
 * return: Returns 0 if successful and nonzero if an error occurs.
 */
int optimize_prob(
    CPXENVptr env,
    CPXLPptr lp,
    double *objval,
    int *solstat,
    double *x,
    int beg,
    int end,
    int verbose
) {
    int status;
    int numcols = CPXgetnumcols(env, lp);

    // Optimize lp and obtain solution.
    status = CPXmipopt(env, lp);
    if (status) {
        fprintf(stderr, "Failed to optimize lp.\n");
        return status;
    }

    // Get solution status.
    *solstat = CPXgetstat(env, lp);

    // Get objective value.
    status = CPXgetobjval(env, lp, objval);
    if (status) {
        fprintf(stderr, "Failed to obtain objective value.\n");
        if (verbose) {
            printf("Solution status: %d\n", *solstat);
        }
        return status;
    }

    // Get variables values.
    status = CPXgetx(env, lp, x, beg, end);
    if (status) {
        fprintf(stderr, "Failed to obtain solution.\n");
        if (verbose) {
            printf("Solution status: %d\n", *solstat);
        }
        return status;
    }

    if (verbose) {
        printf("Solution status: %d\n", *solstat);
        printf("Objective value: %.2f\n", *objval);
        //printf("Variables values:\n");
        //char colname[MAX_COLNAME_LEN];
        //for (int i = 0; i < numcols; i++) { 
        //    if (get_colname(env, lp, i, colname) == 0) {
        //        printf("Column = %d,\tValue = %.2f,\tVar name = %s\n", i, x[i], colname);
        //    } else {
        //        printf("Column = %d,\tValue = %.2f,\tVar name = *error*\n", i, x[i]);
        //    }
        //}
    }
    
    return 0;
}


/**
 * Adds slack variables to the problem pointed by lp. The total number
 * of variables added is given by the 2 * number of columns of the problem.
 *
 * env: A pointer to the CPLEX environment.
 * lp: A pointer to a CPLEX problem object.
 *
 * return: Returns 0 if successful and nonzero if an error occurs.
 */
int add_slack_cols(CPXENVptr env, CPXLPptr lp) {
    int i, tmp, status;
    char ctype;

    int numrows = CPXgetnumrows(env, lp);
    int numcols = CPXgetnumcols(env, lp);

    int ccnt = 2 * numrows; // Numero di colonne (slack) da aggiungere.
                            // I vettori slack hanno dimensione pari al
                            // numero di righe della matrice dei vincoli.
    int nzcnt = 2 * numrows * numrows; // Non zero coefficients counter.
                                       // Tutti i coefficienti sono diversi 
                                       // da zero.

    printf("Slack variables to add: %d\n", ccnt);

    double *matval = NULL;
    int *matbeg = NULL, *matind = NULL;
    char **colnames = NULL;

    matbeg = (int*) malloc(ccnt * sizeof(int));
    matind = (int*) malloc(nzcnt * sizeof(int));
    matval = (double*) malloc(nzcnt * sizeof(double));

    colnames = (char**) malloc(ccnt * sizeof(char*));

    if (matbeg == NULL ||
        matind == NULL ||
        matval == NULL || 
        colnames == NULL) {
        fprintf(stderr, "No memory for adding slack variabiles.\n");
        status = 1;
        goto TERMINATE;
    }

    // Alloco spazio per i nomi.
    for (i = 0; i < ccnt; i++) {
        colnames[i] = (char*) malloc(MAX_SLACK_NAMES_LEN * sizeof(char));
        if (colnames[i] == NULL) {
            fprintf(stderr, "No memory for slack variabiles names.\n");
            status = 1;
            goto TERMINATE;
        }
    }

    // Assegno i nomi alle variabili slack.
    for (i = 0; i < numrows; i++) {
        snprintf(colnames[i], MAX_SLACK_NAMES_LEN, "dp%d", i + 1);
        snprintf(colnames[i + numrows], MAX_SLACK_NAMES_LEN, "dn%d", i + 1);
    }

    // Inizializzo: matbeg.
    for (i = 0; i < ccnt; i++) {
        matbeg[i] = i * numrows;
    }

    // Inizializzo: matind; matval.
    for (i = 0; i < nzcnt; i++) {
        // matind = [0,...,m-1,0,...,m-1,...] dove m = numrows.
        matind[i] = i % numrows; 
        // matval ha i primi nzcnt/2 valori a 1 e i restanti nzcnt/2 a -1.
        matval[i] = (i < nzcnt / 2)? 1 : -1;
    }

    // Ora posso aggiungere le colonne di slack.
    status = CPXaddcols(
        env,
        lp,
        ccnt,
        nzcnt,
        NULL,
        matbeg,
        matind,
        matval,
        NULL,
        NULL,
        colnames
    );
    if (status) {
        fprintf(stderr, "Failed to add new columns to the problem.\n");
        goto TERMINATE;
    }
    // Nota che non aggiorno il numero di colonne. Quindi numcols indica
    // il numero di colonne senza considerare le variabili slack aggiunte.

    // Setto il tipo delle slack.
    // Gli indici delle variabili slack nella matrice dei vincoli, appartengono
    // al seguente range:
    //      [numcols, ..., numcols + ccnt]
    ctype = CPX_CONTINUOUS;
    tmp = numcols + ccnt;
    for (i = numcols; i < tmp; i++) {
        status = CPXchgctype(env, lp, 1, &i, &ctype);
        if (status) {
            fprintf(stderr, "Failed to set slack variable type.\n");
        }
    }

TERMINATE:

    free_and_null((char**) &matbeg);
    free_and_null((char**) &matind);
    free_and_null((char**) &matval);
    for (i = 0; i < ccnt; i++) {
        free_and_null(&colnames[i]);
    }
    free_and_null(colnames);

    return 0;
}


/**
 * Copies the problem pointed by src to the problem pointed by dst.
 *
 * env: A pointer to the CPLEX environment.
 * src: A pointer to a CPLEX problem object. This is the source problem.
 * dst: A pointer to a CPLEX problem object. The copy of src.
 *
 * return: Returns 0 if successful and nonzero if an error occurs.
 */
int copy_prob(CPXENVptr env, CPXLPptr src, CPXLPptr dst) {
    int i, status;

    int numcols = CPXgetnumcols(env, src);
    int numrows = CPXgetnumrows(env, src);
    int numnz = CPXgetnumnz(env, src);

    double *obj = NULL, *rhs = NULL, *matval = NULL, *ub = NULL, *lb = NULL;
    char *sense = NULL, *xctype = NULL;
    int *matbeg = NULL, *matind = NULL, *matcnt = NULL;

    char **colnames = NULL, **rownames = NULL;
    char *colnamestore = NULL, *rownamestore = NULL;

    int nzcnt, surplus;

    // Salvo i nomi delle variabili di SRC.
    status = CPXgetcolname(env, src, NULL, NULL, 0, &surplus, 0, numcols - 1);
    if (status != CPXERR_NEGATIVE_SURPLUS) {
        fprintf(stderr, "Failed to get columns names size.\n");
        goto TERMINATE;
    }

    colnamestore = (char*) malloc(-surplus * sizeof(char));
    colnames = (char**) malloc(numcols * sizeof(char*));
    if (colnames == NULL || colnamestore == NULL) {
        fprintf(stderr, "No memory for saving variables names of SRC.\n");
        status = 1;
        goto TERMINATE;
    }

    status = CPXgetcolname(
        env,
        src,
        colnames,
        colnamestore,
        -surplus,
        &surplus,
        0,
        numcols - 1
    );
    if (status) {
        fprintf(stderr, "Failed to get variables names of SRC.\n");
        goto TERMINATE;
    }

    // Salvo i nomi dei vincoli di SRC.
    status = CPXgetrowname(env, src, NULL, NULL, 0, &surplus, 0, numrows - 1);
    if (status != CPXERR_NEGATIVE_SURPLUS) {
        fprintf(stderr, "Failed to get constraints names size.\n");
        goto TERMINATE;
    }

    rownamestore = (char*) malloc(-surplus * sizeof(char));
    rownames = (char**) malloc(numrows * sizeof(char*));
    if (rownames == NULL || rownamestore == NULL) {
        fprintf(stderr, "No memory for saving constraints names of SRC.\n");
        status = 1;
        goto TERMINATE;
    }

    status = CPXgetrowname(
        env,
        src,
        rownames,
        rownamestore,
        -surplus,
        &surplus,
        0,
        numrows - 1
    );
    if (status) {
        fprintf(stderr, "Failed to get constraints names of SRC.\n");
        goto TERMINATE;
    }

    // Salvo i coefficienti della funzione obiettivo del problema SRC.
    obj = (double*) malloc(numcols * sizeof(double));
    if (obj == NULL) {
        fprintf(stderr, "No memory for saving obj of SRC.\n");
        status = 1;
        goto TERMINATE;
    }

    status = CPXgetobj(env, src, obj, 0, numcols - 1);
    if (status) {
        fprintf(stderr, "Failed to copy obj coefficients of SRC into array.\n");
        goto TERMINATE;
    }

    // Salvo i coefficienti rhs dei vincoli del problema SRC.
    rhs = (double*) malloc(numrows * sizeof(double));
    if (rhs == NULL) {
        fprintf(stderr, "No memory for saving rhs of SRC.\n");
        status = 1;
        goto TERMINATE;
    }

    status = CPXgetrhs(env, src, rhs, 0, numrows - 1);
    if (status) {
        fprintf(stderr, "Failed to copy rhs coefficients of SRC into array.\n");
        goto TERMINATE;
    }

    // Salvo il senso (>=, <=, =) dei vincoli del problema SRC.
    sense = (char*) malloc(numrows * sizeof(double));
    if (sense == NULL) {
        fprintf(stderr, "No memory for saving sense of SRC.\n");
        status = 1;
        goto TERMINATE;
    }

    status = CPXgetsense(env, src, sense, 0, numrows - 1);
    if (status) {
        fprintf(stderr, "Failed to copy sense coefficients of SRC into array.\n");
        goto TERMINATE;
    }

    // Salvo i dati della matrice dei vincoli di SRC.
    matbeg = (int*) malloc(numcols * sizeof(int));
    matind = (int*) malloc(numnz * sizeof(int));
    matcnt = (int*) malloc(numcols * sizeof(int));
    matval = (double*) malloc(numnz * sizeof(double));

    if (matbeg == NULL || matind == NULL || matval == NULL || matcnt == NULL) {
        fprintf(stderr, "No memory for saving matbeg "
                        "or matind or matval or matcnt of SRC.\n");
        status = 1;
        goto TERMINATE;
    }

    status = CPXgetcols(env, src, &nzcnt, matbeg, matind, matval, numnz, &surplus, 0, numcols - 1);
    if (status) {
        fprintf(stderr, "Failed to copy matbeg""or matind "
                        "or matval coefficients of SRC.\n");
        goto TERMINATE;
    }

    // Calcolo delle componenti di matcnt.
    for (i = 0; i < numcols - 1; i++) {
        matcnt[i] = matbeg[i+1] - matbeg[i];
    }
    matcnt[numcols - 1] = nzcnt - matbeg[numcols - 1];

    // Salvo i lower bound e upper bound delle variabili.
    ub = (double*) malloc(numcols * sizeof(double));
    lb = (double*) malloc(numcols * sizeof(double));
    if (ub == NULL || lb == NULL) {
        fprintf(stderr, "No memory for saving lower bounds "
                        "or upper bounds of SRC.\n");
        status = 1;
        goto TERMINATE;
    }

    status = CPXgetub(env, src, ub, 0, numcols - 1);
    if (status) {
        fprintf(stderr, "Failed to copy upper bounds coefficients of SRC.\n");
        goto TERMINATE;
    }

    status = CPXgetlb(env, src, lb, 0, numcols - 1);
    if (status) {
        fprintf(stderr, "Failed to copy upper bounds coefficients of SRC.\n");
        goto TERMINATE;
    }

    // Salvo la tipologia delle variabili del problema SRC.
    xctype = (char*) malloc(numcols * sizeof(char));
    if (xctype == NULL) {
        fprintf(stderr, "No memory for saving variables types of SRC.\n");
        status = 1;
        goto TERMINATE;
    }

    status = CPXgetctype(env, src, xctype, 0, numcols - 1);
    if (status) {
        fprintf(stderr, "Failed to copy upper bounds coefficients of SRC.\n");
        goto TERMINATE;
    }

    // Ora posso crere la copia del problema SRC in DST.
    // TODO: non gestisco rngval (settato a NULL).
    status = CPXcopylpwnames(
        env,
        dst,
        numcols,
        numrows,
        CPXgetobjsen(env, src),
        obj,
        rhs,
        sense,
        matbeg,
        matcnt,
        matind,
        matval,
        lb,
        ub,
        NULL,
        colnames,
        rownames
    );
    if (status) {
        fprintf(stderr, "Failed to populate DST from SRC.\n");
        goto TERMINATE;
    }

    // Setto la tipologia delle variabili del problema DST.
    status = CPXcopyctype(env, dst, xctype);
    if (status) {
        fprintf(stderr, "Failed to set variables types of DST.\n");
        goto TERMINATE;
    }

TERMINATE:

    free_and_null((char**) &obj);
    free_and_null((char**) &rhs);
    free_and_null(&sense);
    free_and_null((char**) &matbeg);
    free_and_null((char**) &matind);
    free_and_null((char**) &matval);
    free_and_null((char**) &ub);
    free_and_null((char**) &lb);
    free_and_null(&xctype);
    free_and_null((char**) &colnames);
    free_and_null(&colnamestore);
    free_and_null((char**) &rownames);
    free_and_null(&rownamestore);

    return 0;
}


/**
 * Creates FMIP starting from MIP. Note that this function doesn't call
 * variable_fixing.
 *
 * env: A pointer to the CPLEX environment.
 * mip: A pointer to a CPLEX problem object.
 * fmip: A pointer to a CPXLPptr that points to a CPLEX problem. 
 *
 * return: Returns 0 if successful and nonzero if an error occurs.
 */
int create_fmip(CPXENVptr env, CPXLPptr mip, CPXLPptr *fmip) {
    int i, cnt, status;
    double tmp;
    
    int numcols = CPXgetnumcols(env, mip);

    *fmip = CPXcreateprob(env, &status, "FMIP");
    if (*fmip == NULL) {
        fprintf(stderr, "Failed to create FMIP.\n");
        return 1;
    }

    status = copy_prob(env, mip, *fmip);
    if (status) {
        fprintf(stderr, "Failed populate FMIP with MIP data.\n");
        return status;
    }

    // Introduco le variabili slack.
    status = add_slack_cols(env, *fmip);
    if (status) {
        fprintf(stderr, "Failed to add slack columns to FMIP.\n");
        return status;
    }
    // Nota che il numero di colonne di FMIP non viene aggiornato dopo
    // aver aggiunto le variabili slack. Ciò significa che numcols indica
    // il numero di colonne senza considerare le variabili di slack.

    // Modifico la funzone obiettivo.
    // Prima setto i coeffcienti delle variabili xi a zero.
    for (i = 0, tmp = 0; i < numcols; i++) {
        status = CPXchgobj(env, *fmip, 1, &i, &tmp);
        if (status) {
            fprintf(stderr, "Failed change variable index %d "
                            "coefficient in obj.\n", i);
            return status;
        }
    }

    // Ora setto i coefficienti delle variabili slack a 1.
    // Gli indici delle variabili slack sono: [numcols, ..., cnt - 1]
    cnt = CPXgetnumcols(env, *fmip);
    for (i = numcols, tmp = 1; i < cnt; i++) {
        status = CPXchgobj(env, *fmip, 1, &i, &tmp);
        if (status) {
            fprintf(stderr, "Failed change variable index %d (slack) "
                            "coefficient in obj.\n", i);
            return status;
        }
    }

    // Cambio il senso del problema a minimizzazione.
    if (CPXgetobjsen(env, *fmip) != CPX_MIN) {
        status = CPXchgobjsen(env, *fmip, CPX_MIN);
        if (status) {
            fprintf(stderr, "Failed setting FMIP obj sense.\n");
            return status;
        }
    }

    return 0;
}


/**
 * Creates OMIP starting from MIP. Note that this function doesn't call
 * variable_fixing.
 *
 * env: A pointer to the CPLEX environment.
 * mip: A pointer to a CPLEX problem object.
 * omip: A pointer to a CPXLPptr that points to a CPLEX problem. 
 *
 * return: Returns 0 if successful and nonzero if an error occurs.
 */
int create_omip(CPXENVptr env, CPXLPptr mip, CPXLPptr *omip, double rhs_slack) {
    int i, cnt, status;
    double tmp;
    
    int numcols = CPXgetnumcols(env, mip);

    char sense;
    char *rowname = "AD";

    int slack_cnt, rmatbeg;
    int *rmatind = NULL;
    double *rmatval = NULL;

    *omip = CPXcreateprob(env, &status, "OMIP");
    if (*omip == NULL) {
        fprintf(stderr, "Failed to create OMIP.\n");
        status = 1;
        goto TERMINATE;
    }

    status = copy_prob(env, mip, *omip);
    if (status) {
        fprintf(stderr, "Failed populate OMIP with MIP data.\n");
        goto TERMINATE;
    }

    // Introduco le variabili slack.
    status = add_slack_cols(env, *omip);
    if (status) {
        fprintf(stderr, "Failed to add slack columns to OMIP.\n");
        goto TERMINATE;
    }
    // Nota che il numero di colonne di OMIP non viene aggiornato dopo
    // aver aggiunto le variabili slack. Ciò significa che numcols indica
    // il numero di colonne senza considerare le variabili di slack.

    // Aggiungo il vincolo che limita il grado della non ammissibilità.
    sense = 'L';
    slack_cnt = CPXgetnumcols(env, *omip) - numcols; // Slack da aggiungere.

    rmatbeg = 0;
    rmatind = (int*) malloc(slack_cnt * sizeof(int));
    rmatval = (double*) malloc(slack_cnt * sizeof(double));
    
    if (rmatind == NULL || rmatval == NULL) {
        fprintf(stderr, "No memory for adding row to OMIP.\n");
        status = 1;
        goto TERMINATE;
    }

    // Inizializzo i valori di matind e matval.
    for (i = 0; i < slack_cnt; i++) {
        rmatind[i] = i + numcols; // Gli indici delle slack: 
                                  //    [numcols, ..., numcols + slack_cnt - 1]
        rmatval[i] = 1;
    }

    // Ora posso aggiungere la nuova riga della matrice dei vincoli.
    status = CPXaddrows(
        env,
        *omip,
        0,
        1,
        slack_cnt,
        &rhs_slack,
        &sense,
        &rmatbeg,
        rmatind,
        rmatval,
        NULL,
        &rowname
    );
    if (status) {
        fprintf(stderr, "Failed to add row to OMIP.\n");
        goto TERMINATE;
    }

TERMINATE:

    free_and_null((char**) &rmatind);
    free_and_null((char**) &rmatval);

    return 0;
}


int main(int argc, char* argv[]) {
    FILE *out_csv;

    CPXENVptr env = NULL;
    CPXLPptr mip = NULL, fmip = NULL, omip = NULL;

    // Variabili per le accedere alle soluzioni dei problemi.
    int numcols_mip, numcols_submip, solstat_fmip, solstat_omip;
    double objval_fmip, objval_omip; // Objective value.
    double *x_fmip = NULL, *x_omip = NULL; // Variabiles value.

    // Variabili per il variable fixing.
    int *int_indices = NULL; // Indici delle variabili intere.
    int *fixed_indices = NULL; // Indici delle variabili che sono state fissate.
    int num_int_vars; // Dimensione dei due array sopra.
    double *initial_vector = NULL; // Vettore di valori iniziale.
    double *lb_mip = NULL, *ub_mip = NULL; // Necessari per resettare i lb e ub
                                           // dopo un variable fixing.
    double slack_sum; // Somma variabili slack (per OMIP).

    int i, j, tmp, cnt, status;
    int sol_loop_detect = 0;
    double old_objv_fmip, old_objv_omip;

    // Setto il seme.
    srand(5);

    // Check command line arguments.
    if (argc != 3) {
        print_usage(argv[0]);
        goto TERMINATE;
    }

    out_csv = fopen(argv[2], "w+");
    if (out_csv == NULL) {
        fprintf(stderr, "The output file %s could not be opened.\n", argv[2]);
        goto TERMINATE;
    }

    printf("INPUT FILE: %s\n", argv[1]);

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

    // Set limit to explorable nodes.
    status = CPXsetintparam(env, CPX_PARAM_NODELIM, NODE_LIMIT);
    if (status) {
        fprintf(stderr, "Failure to set node limit, error %d\n", status);
        goto TERMINATE;
    }

    // Set time limit.
    status = CPXsetdblparam(env, CPXPARAM_TimeLimit, TIME_LIMIT);
    if (status) {
        fprintf(stderr, "Failure to set time limit, error %d\n", status);
        goto TERMINATE;
    }

    // ### Creazione del MIP, FMIP e OMIP ###
    printf("\nCreating MIP.\n");
    mip = CPXcreateprob(env, &status, argv[1]);
    if (mip == NULL) {
        fprintf(stderr, "Failed to create MIP.\n");
        goto TERMINATE;
    }

    // Read file and copy the data into the created mip1.
    status = CPXreadcopyprob(env, mip, argv[1], NULL);
    if (status) {
        fprintf(stderr, "Failed to read and copy the problem data (MIP).\n");
        goto TERMINATE;
    }

    printf("\nCreating FMIP.\n");
    status = create_fmip(env, mip, &fmip);
    if (status) {
        fprintf(stderr, "Failed to create FMIP.\n");
        goto TERMINATE;
    }

    // Salvo il problema FMIP in un file.
    //status = CPXwriteprob(env, fmip, "fmip.lp", NULL);
    //if (status) {
    //    fprintf (stderr, "Failed to write FMIP to disk.\n");
    //    goto TERMINATE;
    //}

    printf("\nCreating OMIP.\n");
    status = create_omip(env, mip, &omip, 200);
    if (status) {
        fprintf(stderr, "Failed to create OMIP.\n");
        goto TERMINATE;
    }

    // ### Inizializzazione array con informazioni problema MIP ###
    // Trovo quante sono le variabili intere di MIP.
    numcols_mip = CPXgetnumcols(env, mip);
    for (i = 0, num_int_vars = 0; i < numcols_mip; i++) {
        char type;

        status = CPXgetctype(env, mip, &type, i, i);
        if (status) {
            fprintf(stderr, "Failed to get variable %d type.", i);
            goto TERMINATE;
        }

        // TODO: Da capire cosa sono le semi-integer variables e se devo
        //       considerarle.
        if (type == CPX_BINARY || type == CPX_INTEGER) {
            num_int_vars += 1;
        }
    }

    // Inizializzo i vettori contenenti le informazioni del problema MIP:
    // lower bound, upper bound e indici delle variabili intere.
    int_indices = (int*) malloc(num_int_vars * sizeof(int));
    fixed_indices = (int*) malloc(num_int_vars * sizeof(int));

    lb_mip = (double*) malloc(numcols_mip * sizeof(double));
    ub_mip = (double*) malloc(numcols_mip * sizeof(double));

    if (int_indices == NULL ||
        fixed_indices == NULL ||
        lb_mip == NULL ||
        ub_mip == NULL) {
        fprintf(stderr, "No memory for int_indices, fixed_indices,"
                        "lb_mip and ub_mip.\n");
        goto TERMINATE;
    }

    status = init_mip_bds_and_indices(
        env,
        mip,
        lb_mip,
        ub_mip,
        numcols_mip,
        int_indices,
        num_int_vars
    );
    if (status) {
        fprintf(stderr, "Failed to initialize arrays with MIP information.\n");
        goto TERMINATE;
    }

    // Inizializione vettore di valori iniziale.
    initial_vector = (double*) malloc(numcols_mip * sizeof(double));
    if (initial_vector == NULL) {
        fprintf(stderr, "No memory for initial_vector.\n");
        goto TERMINATE;
    }

    init_initial_vector(
        initial_vector,
        int_indices,
        num_int_vars,
        lb_mip,
        ub_mip,
        BOUND_CONSTANT
    );

    // Alloco lo spazio per le soluzioni di FMIP.
    numcols_submip = CPXgetnumcols(env, fmip);
    x_fmip = (double*) malloc(numcols_submip * sizeof(double));
    if (x_fmip == NULL) {
        fprintf(stderr, "No memory for solution values for FMIP.\n");
        goto TERMINATE;
    }

    // Alloco lo spazio per le soluzioni di OMIP.
    x_omip = (double*) malloc(numcols_submip * sizeof(double));
    if (x_omip == NULL) {
        fprintf(stderr, "No memory for solution values for OMIP.\n");
        goto TERMINATE;
    }

    // ### Inizio ciclo di risoluzione ###
    for (cnt = 0; cnt < MAX_RUN; cnt++) {
        // Risolvo l'FMIP: massimo 'MAX_ATTEMPTS' tentativi.
        for (i = 0; i < MAX_ATTEMPTS; i++) {
            // Fisso le variabili di FMIP.
            printf("\n### Variable fixing on FMIP "
                   "- Attempt %d - Run %d ###\n", i, cnt);
            bzero(fixed_indices, num_int_vars * sizeof(int)); // Inizializzo a
                                                              // zero.
            if (cnt == 0) { // Uso l'initial vector al primo ciclo.
                status = variable_fixing(
                    env,
                    fmip,
                    int_indices,
                    fixed_indices,
                    num_int_vars,
                    initial_vector,
                    30
                );
            } else { // Altrimenti uso i valori di OMIP.
                if (sol_loop_detect == 1) { // Gestione loop.
                    status = variable_fixing(
                        env,
                        fmip,
                        int_indices,
                        fixed_indices,
                        num_int_vars,
                        x_omip,
                        30
                    );
                    sol_loop_detect = 0;
                } else {
                    status = variable_fixing(
                        env,
                        fmip,
                        int_indices,
                        fixed_indices,
                        num_int_vars,
                        x_omip,
                        50
                    );
                }
            }
            if (status) {
                fprintf(stderr, "Failed to fix variables of FMIP.\n");
                goto TERMINATE;
            }

            // Ottimizzo FMIP.
            status = optimize_prob(
                env,
                fmip,
                &objval_fmip,
                &solstat_fmip,
                x_fmip,
                0,
                numcols_submip - 1,
                1
            );
            if (status) {
                if (solstat_fmip == CPXMIP_INFEASIBLE ||
                    solstat_fmip == CPXMIP_NODE_LIM_INFEAS ||
                    solstat_fmip == CPXMIP_TIME_LIM_INFEAS) {
                    printf("FMIP is infeasible.\n");
                } else {
                    fprintf(stderr, "Failed to optimize FMIP.\n"
                                    "Solution status: %d\n", solstat_fmip);
                    goto TERMINATE;
                }
            }

            // Restore dei bounds di FMIP.
            status = restore_bounds(
                env,
                fmip,
                lb_mip,
                ub_mip,
                numcols_mip,
                int_indices,
                fixed_indices,
                num_int_vars
            );
            if (status) {
                fprintf (stderr, "Failed to restore FMIP bounds.\n");
                goto TERMINATE;
            }

            // Se ho trovato una soluzione ottima.
            if (solstat_fmip == CPXMIP_OPTIMAL ||
                solstat_fmip == CPXMIP_OPTIMAL_TOL ||
                solstat_fmip == CPXMIP_NODE_LIM_FEAS ||
                solstat_fmip == CPXMIP_TIME_LIM_FEAS) {
                printf("Found an optimal solution for FMIP.\n");
                fprintf(out_csv, "fmip,%.2f\n", objval_fmip); // Salvo in output.
                break;
            } 
        }

        // Se nessun FMIP è stato risolto esco.
        if (solstat_fmip == CPXMIP_INFEASIBLE ||
            solstat_fmip == CPXMIP_NODE_LIM_INFEAS ||
            solstat_fmip == CPXMIP_TIME_LIM_INFEAS) {
            printf("All FMIPs generated was infeasibile.\n");
            break;
        }
        // TODO: gestire il caso in cui objval_fmip è zero: cosa faccio?

        // Aggiorno il valore rhs del vincolo sulle slack di OMIP.
        tmp = CPXgetnumrows(env, omip) - 1;
        status = CPXchgrhs(env, omip, 1, &tmp, &objval_fmip);
        if (status) {
            fprintf(stderr, "Failed to update rhs value of OMIP.\n");
            goto TERMINATE;
        }

        // Risolvo l'OMIP: massimo 'MAX_ATTEMPTS' tentativi.
        for (i = 0; i < MAX_ATTEMPTS; i++) {
            // Fisso le variabili di OMIP.
            printf("\n### Variable fixing on OMIP "
                   "- Attempt %d - Run %d ###\n", i, cnt);
            bzero(fixed_indices, num_int_vars * sizeof(int)); // Inizializzo a 
                                                              // zero.
            status = variable_fixing(
                env,
                omip,
                int_indices,
                fixed_indices,
                num_int_vars,
                x_fmip,
                50
            );
            if (status) {
                fprintf(stderr, "Failed to fix variables of FMIP.\n");
                goto TERMINATE;
            }

            // Ottimizzo OMIP.
            status = optimize_prob(
                env,
                omip,
                &objval_omip,
                &solstat_omip,
                x_omip,
                0,
                numcols_submip - 1,
                1
            );
            if (status) {
                if (solstat_omip == CPXMIP_INFEASIBLE ||
                    solstat_omip == CPXMIP_NODE_LIM_INFEAS ||
                    solstat_omip == CPXMIP_TIME_LIM_INFEAS) {
                    printf("OMIP is infeasible.\n");
                } else {
                    fprintf(stderr, "Failed to optimize OMIP."
                                    "\n Solution status: %d", solstat_omip);
                    goto TERMINATE;
                }
            }

            // Restore dei bounds di OMIP.
            status = restore_bounds(
                env,
                omip,
                lb_mip,
                ub_mip,
                numcols_mip,
                int_indices,
                fixed_indices,
                num_int_vars
            );
            if (status) {
                fprintf (stderr, "Failed to restore OMIP bounds.\n");
                goto TERMINATE;
            }

            // Se ho trovato una soluzione ottima.
            if (solstat_omip == CPXMIP_OPTIMAL ||
                solstat_omip == CPXMIP_OPTIMAL_TOL ||
                solstat_omip == CPXMIP_NODE_LIM_FEAS ||
                solstat_omip == CPXMIP_TIME_LIM_FEAS) {
                slack_sum = sum(x_omip, numcols_mip, numcols_submip - 1);
                printf("Slack sum: %.2f\n", slack_sum);
                printf("Found an optimal solution for OMIP.\n");
                fprintf(out_csv, "omip,%.2f\n", objval_omip); // Salvo in output.
                if (slack_sum == 0) {
                    printf("Found a feasibile solution for MIP.\n");
                    goto TERMINATE;
                } else {
                    break;
                }
            }
            // TODO: devono essere gestiti anche gli altri status che forniscono
            //       una soluzione ottima.
        }

        // Se nessun OMIP è stato risolto esco.
        if (solstat_omip == CPXMIP_INFEASIBLE ||
            solstat_omip == CPXMIP_NODE_LIM_INFEAS ||
            solstat_omip == CPXMIP_TIME_LIM_INFEAS) {
            printf("All OMIPs generated was infeasibile.\n");
            break;
        }

        // Controllo la presenza di loop.
        if (cnt == 0) {
            old_objv_fmip = objval_fmip;
            old_objv_omip = objval_omip;
        } else {
            if (old_objv_fmip == objval_fmip && old_objv_omip == objval_omip) {
                sol_loop_detect = 1;
                printf("LOOP DETECT\n");
            } else {
                old_objv_fmip = objval_fmip;
                old_objv_omip = objval_omip;
            }
        }
    }

TERMINATE:

    // Free up solution.
    free_and_null((char**) &x_fmip);
    free_and_null((char**) &x_omip);
    free_and_null((char**) &int_indices);
    free_and_null((char**) &fixed_indices);
    free_and_null((char**) &lb_mip);
    free_and_null((char**) &ub_mip);

    // Free up the problem allocated by CPXcreateprob, if necessary.
    if (mip != NULL) {
        status = CPXfreeprob(env, &mip);
        if (status) {
            fprintf(stderr, "CPXfreeprob failed, error code %d.\n", status);
        }
    }

    if (fmip != NULL) {
        status = CPXfreeprob(env, &fmip);
        if (status) {
            fprintf(stderr, "CPXfreeprob failed, error code %d.\n", status);
        }
    }

    if (omip != NULL) {
        status = CPXfreeprob(env, &fmip);
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
