#include <ilcplex/cplex.h>
#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <unistd.h>

/* Maximum lenght of the slack variables names (considering '\0').
 * The slack vectors Delta+ and Delta- have the same length.
 * The variables ara named as follow:
 *     Delta+ = [dp1, dp2, dp3, ..., dpi, ...,]
 *     Delta- = [dn1, dn2, dn3, ..., dni, ...,]
 * Note that 2 characters are taken by 'd' and 'p' or 'n' and 1 character
 * is taken by '\0'. */
#define MAX_SLACK_NAME_LEN 9

/* Maximum column's name length. */
#define MAX_COLNAME_LEN 9

/* Maximum number of attempts to solve OMIP. */
#define MAX_ATTEMPTS 1

/* Maximum number of nodes explored during optimization. */
#define NODE_LIMIT 1000 /* Not set */

/* Optimizer time limit (seconds). */
#define TIME_LIMIT 60 /* Not set */

/* Maximum number of iteration. An iteration is complete when the solver has 
 * found a feasibile solution for FMIP and than a feasibile solution for OMIP. */
#define MAX_ITR 10

/* It's used for initializing the starting vector. The value of the variable
 * is raondomly choosen from the range
 *      [max(lb, BOUND_CONSTANT), min(ub, BOUND_CONSTANT)]
 *  where lb is the lower bound of the variable and ub is its upper bound. */
#define BOUND_CONSTANT 10E5


/**
 * Print help message.
 *
 * progname: The executable's name.
 */
void print_usage(char *progname) {
   fprintf (stderr, "Usage: %s -i <intput> -o <output> -s <seed> -p <%%varfixing>\n"
                    "   input: is a file with extension \n"
                    "      MPS, SAV, or LP (lower case is allowed)\n"
                    "   output: is a file with extension csv\n"
                    "   seed: any integer number\n"
                    "   %%varfixing: is the percentage of fixed variables\n"
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
 *              returned. This array must be of length at least num_int_vars.
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

    // Save the indices of the integer variables of MIP.
    for (i = 0, j = 0; i < numcols && j < num_int_vars; i++) {
        char type;
        status = CPXgetctype(env, mip, &type, i, i);
        if (status) {
            fprintf(stderr, "Failed to get variable %d type.\n", i);
            return status;
        }

        if (type == CPX_BINARY || type == CPX_INTEGER) {
            int_indices[j] = i; // Save index.
            j += 1;
        }
    }

    // Save lower bounds and upper bounds of MIP.
    status = CPXgetub(env, mip, ub, 0, numcols - 1);
    if (status) {
        fprintf(stderr, "Failed to copy upper bounds coefficients of MIP.\n");
        return status;
    }

    status = CPXgetlb(env, mip, lb, 0, numcols - 1);
    if (status) {
        fprintf(stderr, "Failed to copy lower bounds coefficients of MIP.\n");
        return status;
    }

    return 0;
}


/**
 * Restores upper and lower bounds of lp with the specified values only for the
 * variables that was fixed by variable_fixing.
 *
 * env: A pointer to the CPLEX environment.
 * mip: A pointer to a CPLEX problem object.
 * lb: An array where are stored all the lower bounds of lp.
 *     This array must be of length at least numcols.
 * ub: An array where are stored all the upper bounds of lp.
 *     This array must be of length at least numcols.
 * numcols: The length of the arrays lb and up.
 * int_indices: An array where are stored the indices of the integer variables
 *              of lp. This array must be of length at least num_int_vars.
 * is_fixed: An array that tells if a variable was fixed.
 *           A variable with index int_indices[i] was fixed if and only if 
 *           fixed_indices[i] is a nonzero value.
 *           This array must be of length at least num_int_vars.
 * num_int_vars: The length of the arrays int_indices and fixed_indices.
 *
 * return: Returns 0 if successful and nonzero if an error occurs.
 */
int restore_bounds(
    CPXENVptr env,
    CPXLPptr mip,
    double *lb,
    double *ub,
    int *int_indices,
    int *is_fixed,
    int num_int_vars
) {
    int i, status;
    char lower = 'L', upper = 'U';

    for (i = 0; i < num_int_vars; i++) {
        if (is_fixed[i]) {
            // Reset lower bounds.
            status = CPXchgbds(
                env,
                mip,
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
            // Reset upper bounds.
            status = CPXchgbds(
                env,
                mip,
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
 * Fixes some of the problem variables to a value specified by the array
 * named x. An index is randomly choose from the indices stored in the array 
 * named int_indices. The variables are fixed starting from the choosen index.
 * The ammount of variables that is fixed is set by the percentage parameter.
 *
 * For example:
 *      int_indices = {0, 2, 3}, x = {100, ?, 66, 4},
 *      fixed_indices = {0, 0, 0}, percentage = 50
 *      
 *      If the random generetor generates the value 2 then the variable
 *      with index int_indices[2] is fixed to the value x[int_indices[2]].
 *      Then in fixed_indices[2] will be stored 1.
 *
 *      In other words variable with index 3 will be fixed to 4 and 
 *      fixed_indices = {0, 0, 1}.
 *
 * env: A pointer to the CPLEX environment.
 * lp: A pointer to a CPLEX problem object.
 * int_indices: An array where are stored the indices of the integer variables
 *              of lp. This array must be of length at least num_int_vars.
 * is_fixed: An array that tells if a variable was fixed.
 *           A variable with index int_indices[i] was fixed if and only if 
 *           fixed_indices[i] is a nonzero value. This array must be of length
 *           at least num_int_vars.
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
    int *is_fixed,
    int num_int_vars,
    double *x,
    int rho
) {
    int i, tmp, cnt, rnd, status;
    double val;
    char lu = 'B';

    // Number of variables to fix.
    cnt = num_int_vars * rho / 100;
    printf("Variables to fix: %d\n", cnt);

    // Choose a random position of 'int_indices', 
    // aka choose an integer variable's index.
    rnd = rand() % num_int_vars;

    // Fix 'cnt' integer variables starting form the previously choosen index.
    for (i = 0; i < cnt; i++) {
        // Prelevo il valore.
        tmp = (rnd + i) % num_int_vars;
        val = x[int_indices[tmp]];

        // Fix variable's value.
        status = CPXchgbds(env, lp, 1, &int_indices[tmp], &lu, &val);
        if (status) {
            fprintf(stderr, 
                    "Failed to fix variable with index %d.\n", int_indices[rnd]);
            return status;
        }
        is_fixed[tmp] = 1;
    }

    return 0;
}


/**
 * Optimizes the problem pointed by lp.
 *
 * env: A pointer to the CPLEX environment.
 * lp: A pointer to a CPLEX problem object. The problem pointed by lb can be
 *     both a mixed-integer or a linear problem.
 * objval: A pointer where objective optimal value is to be returned.
 * solstat: A pointer where the status code is to be returned. 
 * x: An array where the values of the variables are to be returned. 
 *    This array must be of length at least (end - beg + 1).
 * beg: Specifies the beginning of the range of the variables values to be
 *      returned.
 * end: Specifies the end of the range of the variables values to be returned.
 * verbose: Manage output verbosity. If it's set to 1 objval and solstat 
 *          are printed to standard output. If the value is greater than 1
 *          objval, solstat and all the variables values are printed to
 *          standard output.
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
    int numcols = CPXgetnumcols(env, lp), prob_t = CPXgetprobtype(env, lp);

    // Optimize lp.
    switch (prob_t) {
        case CPXPROB_LP:
            status = CPXlpopt(env, lp);
            break;
        case CPXPROB_MILP:
            status = CPXmipopt(env, lp);
            break;
        default:
            fprintf(stderr, "Type of problem (%d) not supported.\n", prob_t);
            return 1;
    }
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

    // If verbose is set print solution info.
    if (verbose) {
        printf("Solution status: %d\n", *solstat);
        printf("Objective value: %f\n", *objval);
        if (verbose > 1) {
            printf("Variables values:\n");
            char colname[MAX_COLNAME_LEN];
            for (int i = 0; i < numcols; i++) { 
                if (get_colname(env, lp, i, colname) == 0) {
                    printf("Column = %d,\tValue = %f,\tVar name = %s\n",
                            i, x[i], colname);
                } else {
                    printf("Column = %d,\tValue = %f,\tVar name = *error*\n",
                            i, x[i]);
                }
            }
        }
    }
    
    return 0;
}


/**
 * Adds slack variables to the problem pointed by lp. The total number
 * of variables added is given by the 2 * number of columns of the problem
 * pointed by lp.
 *
 * env: A pointer to the CPLEX environment.
 * lp: A pointer to a CPLEX problem object.
 *
 * return: Returns 0 if successful and nonzero if an error occurs.
 */
int add_slack_cols(CPXENVptr env, CPXLPptr lp) {
    int i, tmp, status;
    char ctype, sense;

    int numrows = CPXgetnumrows(env, lp);
    int numcols = CPXgetnumcols(env, lp);

    int ccnt = 2 * numrows; // Tot slack vars to add.
    int nzcnt = 2 * numrows; // Non zero coefficients counter.

    double *matval = NULL;
    int *matbeg = NULL, *matind = NULL;
    char **colnames = NULL;

    printf("Slack variables to add: %d\n", ccnt);

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

    // Space for the names.
    for (i = 0; i < ccnt; i++) {
        colnames[i] = (char*) malloc(MAX_SLACK_NAME_LEN * sizeof(char));
        if (colnames[i] == NULL) {
            fprintf(stderr, "No memory for slack variabiles names.\n");
            status = 1;
            goto TERMINATE;
        }
    }

    // Assign slack variables names.
    for (i = 0; i < numrows; i++) {
        snprintf(colnames[i], MAX_SLACK_NAME_LEN, "dp%d", i + 1);
        snprintf(colnames[i + numrows], MAX_SLACK_NAME_LEN, "dn%d", i + 1);
    }

    // Init: matbeg.
    for (i = 0; i < ccnt; i++) {
        matbeg[i] = i;
    }

    // Init: matind; matval.
    for (i = 0; i < numrows; i++) {
        // matind = [0,...,m-1,0,...,m-1,...] where m = numrows.
        matind[i] = matind[i + numrows] = i; 

        // matval = first numrows values set to 1, the remaning numrows to -1.
        status = CPXgetsense(env, lp, &sense, i, i);
        if (status) {
            fprintf(stderr, "Failed to get constraint %d sense.\n", i);
            goto TERMINATE;
        }

        // Scheme: (s := slack var)
        //      if   ax <= b   then   ax - s       <= b
        //      if   ax >= b   then   ax + s       >= b
        //      if   ax  = b   then   ax + s1 - s2  = b
        switch (sense) {
            case 'L':
                matval[i] = 0;
                matval[i + numrows] = -1;
                break;
            case 'G':
                matval[i] = 1;
                matval[i + numrows] = 0;
                break;
            case 'E':
                matval[i] = 1;
                matval[i + numrows] = -1;
                break;
        }
        // TODO: Non viene gestito il ranged constraint.
    }

    // Adding slack variables (columns).
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
    // Note that 'numcols' is not updated, thus it contains the number of
    // columns without considering the slack variables that we added before.

    // Set the type of the slack variables.
    // Their indices belong to the following range:
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
    free_and_null((char**) &colnames);

    return status;
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

    // Save the variables names of SRC.
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

    // Save the constraints names of SRC.
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

    // Save the objective function's coefficients of SRC.
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

    // Save the right hand side coefficients of SRC.
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

    // Save the righ hand side sense (thus >=, <=, =) of SRC.
    sense = (char*) malloc(numrows * sizeof(char));
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

    // Save the constraints matrix's information of SRC.
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

    status = CPXgetcols(
        env,
        src,
        &nzcnt,
        matbeg,
        matind,
        matval,
        numnz,
        &surplus,
        0,
        numcols - 1
    );
    if (status) {
        fprintf(stderr, "Failed to copy matbeg""or matind "
                        "or matval coefficients of SRC.\n");
        goto TERMINATE;
    }

    // Calculate matcnt's values.
    for (i = 0; i < numcols - 1; i++) {
        matcnt[i] = matbeg[i+1] - matbeg[i];
    }
    matcnt[numcols - 1] = nzcnt - matbeg[numcols - 1];

    // Save lower bounds and upperbounds of SRC.
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
        fprintf(stderr, "Failed to copy lower bounds coefficients of SRC.\n");
        goto TERMINATE;
    }

    // Save variables types of SRC.
    xctype = (char*) malloc(numcols * sizeof(char));
    if (xctype == NULL) {
        fprintf(stderr, "No memory for saving variables types of SRC.\n");
        status = 1;
        goto TERMINATE;
    }

    status = CPXgetctype(env, src, xctype, 0, numcols - 1);
    if (status) {
        fprintf(stderr, "Failed to copy variables types of SRC.\n");
        goto TERMINATE;
    }

    // Create DST, the copy of SRC.
    // TODO: manage rngval (is the NULL parameter).
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

    // Set variables types of DST.
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
    free_and_null((char**) &matcnt);
    free_and_null((char**) &ub);
    free_and_null((char**) &lb);
    free_and_null(&xctype);
    free_and_null((char**) &colnames);
    free_and_null(&colnamestore);
    free_and_null((char**) &rownames);
    free_and_null(&rownamestore);

    return status;
}


/**
 * Generate the starting point, thus the initial vector for the first
 * FMIP problem to solve.
 *
 * env: A pointer to the CPLEX environment.
 * fmip: A pointer to a CPLEX problem object.
 * initial_vector: An array where the generate values are to be returned.
 *                 This array must be of length at least the same as the number
 *                 of columns of the starting MIP problem.
 * int_indices: An array where are stored the indices of the integer variables
 *              of fmip. This array must be of length at least num_int_vars.
 * num_int_vars: The length of int_indices.
 * lb: An array where are stored all the lower bounds. This array must be of 
 *     length at least as the number of columns of the starting MIP problem.
 * ub: An array where are stored all the upper bounds. This array must be of 
 *     length at least as the number of columns of the starting MIP problem.
 * cb: The fixed bound constant.
 * theta: Minimum percentage of integer variables fixed at every iteration.
 */
int generate_starting_vector(
    CPXENVptr env,
    CPXLPptr fmip,
    double *initial_vector,
    int *int_indices,
    int num_int_vars,
    double *lb,
    double *ub,
    double cb,
    double theta
) {
    int i, j, beg, prev_beg, cnt, tmp, tot_fixed, upper, lower;
    double rnd, objval, *x_relax;

    int status, solstat, numcols = CPXgetnumcols(env, fmip);;

    char *var_type = NULL;
    int *is_fixed = NULL;
    char lu = 'B', l = 'L', u = 'U', c_type = CPX_CONTINUOUS;

    int tot_R, tot_O, tot_0;

    // Space for variables types of FMIP.
    var_type = (char*) malloc(num_int_vars * sizeof(char));
    if (var_type == NULL) {
        fprintf(stderr, "No memory for var_type.\n");
        goto  TERMINATE;
    }

    // Relax FMIP: set all integer variables to continous.
    for (i = 0; i < num_int_vars; i++) {
        // Save current type.
        status = CPXgetctype(
            env,
            fmip,
            &var_type[i],
            int_indices[i],
            int_indices[i]
        );
        if (status) {
            fprintf(stderr, "Failed to get variable index %d type.\n",
                            int_indices[i]);
            goto TERMINATE;
        }
        // Change to continous.
        status = CPXchgctype(env, fmip, 1, &int_indices[i], &c_type);
        if (status) {
            fprintf(stderr, "Failed to change variable index %d type.\n", 
                             int_indices[i]);
            goto TERMINATE;
        }
    }

    // Set FMIP to linear problem.
    status = CPXchgprobtype(env, fmip, CPXPROB_LP);
    if (status) {
        fprintf(stderr, "Failed to change FMIP problem type.\n");
        goto TERMINATE;
    }

    // Space for the FMIP relax solutions.
    x_relax = (double*) malloc(numcols * sizeof(double));
    is_fixed = (int*) malloc(num_int_vars * sizeof(int));
    if (x_relax == NULL || is_fixed == NULL) {
        fprintf(stderr, "No memory for x_relax and is_fixed (FMIP relax).\n");
        goto  TERMINATE;
    }

    // Init is_fixed with zeros (no value fixed).
    bzero(is_fixed, num_int_vars * sizeof(int));

    // Generate starting vector algorithm.
    prev_beg = beg = tot_fixed = 0;
    tmp = ceil(num_int_vars * theta / 100);
    for (j = 0; tot_fixed != num_int_vars; j++) {
        // Fix integer variables to random integer within bounds.
        printf("Generating random values . . .\n");
        prev_beg = beg;
        for (cnt = 0; beg < num_int_vars && cnt < tmp; beg++){
            // If value was already fixed, go ahead.
            if (is_fixed[beg]) {
                continue;
            }

            // Generate random value.
            upper = (int) min(ub[int_indices[beg]], cb);
            lower = (int) max(lb[int_indices[beg]], -cb);
            rnd = (int) (rand() % (upper - lower + 1)) + lower;

            // Save the value.
            initial_vector[int_indices[beg]] = rnd;
            is_fixed[beg] = 'R'; // Random fix.

            // Fix the variable to the generated value in FMIP relax.
            status = CPXchgbds(env, fmip, 1, &int_indices[beg], &lu, &rnd);
            if (status) {
                fprintf(stderr, 
                        "Failed to fix variable %d.\n", int_indices[beg]);
                goto TERMINATE;
            }

            cnt += 1; // A variable was fixed.
        }
        tot_fixed += cnt; // Update total fixed counter.

        // If all values has been fixed don't resolve the relaxation.
        if (tot_fixed == num_int_vars) {
            for (tot_R = 0, tot_O = 0, tot_0 = 0, i = 0; i < num_int_vars; i++) {
                if (is_fixed[i] == 'R') tot_R += 1;
                else if (is_fixed[i] == 'O') tot_O += 1;
                else if (!is_fixed[i]) tot_0 += 1;
            }
            printf("Fixed values situation: Random: %d, Optimize: %d, "
                   "Not fixed: %d\n", tot_R, tot_O, tot_0);
            break;
        }

        // Optimize FMIP relax.
        printf("Optimize FMIP relaxation . . . - Iteration %d\n", j);
        status = optimize_prob(
            env,
            fmip,
            &objval,
            &solstat,
            x_relax,
            0,
            numcols - 1,
            1
        );
        if (status) {
            fprintf(stderr, "Error during FMIP relax optimization.\n");
        } 

        // Check solution status.
        switch (solstat) {
            case CPX_STAT_OPTIMAL:
            case CPX_STAT_FEASIBLE:
                printf("Found a feasibile solution for FMIP relax.\n");
                // Fix extra values: take them from relaxation solution.
                for (i = beg; i < num_int_vars; i++) {
                    if (!is_fixed[i] &&                       // If is int.
                        round(x_relax[int_indices[i]]) == 
                        x_relax[int_indices[i]]) {
                        initial_vector[int_indices[i]] = x_relax[int_indices[i]];
                        is_fixed[i] = 'O'; // Optimize fix.
                        tot_fixed += 1; // Update total fixed counter.

                        // Fix the variable in FMIP relax.
                        status = CPXchgbds(
                            env,
                            fmip,
                            1,
                            &int_indices[i],
                            &lu,
                            &x_relax[int_indices[i]]
                        );
                        if (status) {
                            fprintf(stderr, 
                                    "Failed to fix variable "
                                    "with index %d.\n", int_indices[i]);
                            goto TERMINATE;
                        }
                    }
                }
                break;
            case CPX_STAT_ABORT_DETTIME_LIM:
                tmp *= 50;
                printf("Increased fixed random variables "
                       "at each iteration due to abort time error.\n");
            default:
                fprintf(stderr, "Failed to optimize FMIP relax.\n");
        }

        for (tot_R = 0, tot_O = 0, tot_0 = 0, i = 0; i < num_int_vars; i++) {
            if (is_fixed[i] == 'R') tot_R += 1;
            else if (is_fixed[i] == 'O') tot_O += 1;
            else if (!is_fixed[i]) tot_0 += 1;
        }
        printf("Fixed values situation: Random: %d, Optimize: %d, "
               "Not fixed: %d\n", tot_R, tot_O, tot_0);
    }

    // Set FMIP back to mixed-integer.
    status = CPXchgprobtype(env, fmip, CPXPROB_MILP);
    if (status) {
        fprintf(stderr, "Failed to change FMIP problem type.\n");
        goto TERMINATE;
    }

    // Restore FMIP variables types.
    for (i = 0; i < num_int_vars; i++) {
        status = CPXchgctype(env, fmip, 1, &int_indices[i], &var_type[i]);
        if (status) {
            fprintf(stderr, "Failed to change variable %d type.\n", 
                            int_indices[i]);
            goto TERMINATE;
        }
    }
    
    // Restore FMIP bounds.
    status = restore_bounds(
        env,
        fmip,
        lb,
        ub,
        int_indices,
        is_fixed,
        num_int_vars
    );
    if (status) {
        fprintf(stderr, "Failed to restore FMIP bounds.\n");
    }

TERMINATE:

    free_and_null((char**) &x_relax);
    free_and_null((char**) &is_fixed);
    free_and_null(&var_type);

    return status;
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

    // Add slack variables.
    status = add_slack_cols(env, *fmip);
    if (status) {
        fprintf(stderr, "Failed to add slack columns to FMIP.\n");
        return status;
    }
    // Note that 'numcols' is not updated, thus it contains the number of
    // columns without considering the slack variables that we added before.

    // Modify objective function.
    // First set all variables coefficients to zero.
    for (i = 0, tmp = 0; i < numcols; i++) {
        status = CPXchgobj(env, *fmip, 1, &i, &tmp);
        if (status) {
            fprintf(stderr, "Failed change variable index %d "
                            "coefficient in obj.\n", i);
            return status;
        }
    }

    // Then set slack coefficients to 1. Slack's indices goes from:
    //      [numcols, ..., cnt - 1]
    cnt = CPXgetnumcols(env, *fmip);
    for (i = numcols, tmp = 1; i < cnt; i++) {
        status = CPXchgobj(env, *fmip, 1, &i, &tmp);
        if (status) {
            fprintf(stderr, "Failed change variable index %d (slack) "
                            "coefficient in obj.\n", i);
            return status;
        }
    }

    // Change objective sense to minimization.
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

    // Add slack variables.
    status = add_slack_cols(env, *omip);
    if (status) {
        fprintf(stderr, "Failed to add slack columns to OMIP.\n");
        goto TERMINATE;
    }
    // Note that 'numcols' is not updated, thus it contains the number of
    // columns without considering the slack variables that we added before.

    // Add constraint (row) that limits the sum of the slack variables.
    sense = 'L';
    slack_cnt = CPXgetnumcols(env, *omip) - numcols; // Num of slack to add.

    rmatbeg = 0;
    rmatind = (int*) malloc(slack_cnt * sizeof(int));
    rmatval = (double*) malloc(slack_cnt * sizeof(double));
    
    if (rmatind == NULL || rmatval == NULL) {
        fprintf(stderr, "No memory for adding row to OMIP.\n");
        status = 1;
        goto TERMINATE;
    }

    // Init: matind; matval.
    for (i = 0; i < slack_cnt; i++) {
        rmatind[i] = i + numcols; // Slack's indices goes from: 
                                  //   [numcols, ..., numcols + slack_cnt - 1].
        rmatval[i] = 1; // All coefficients are set to 1.
    }

    // Add the new constraint to the matrix.
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

    return status;
}


int main(int argc, char *argv[]) {
    CPXENVptr env = NULL;
    CPXLPptr mip = NULL, fmip = NULL, omip = NULL;

    int numcols_mip, numcols_submip; // Number of columns.
    int solstat_fmip, solstat_omip;  // Solution status.
    double objval_fmip, objval_omip; // Objective values.
    int numnz_mip;                   // Number of nonzero elements.

    double *x_fmip = NULL, *x_omip = NULL; // Variables values.

    int *int_indices = NULL; // Indices of the integer variables of MIP.
    int *is_fixed = NULL;    // Manages fixed variables.
    int num_int_vars;        // Ammount of integer variables of MIP.

    int *varindices_submip = NULL; // Indices of the variables in the subMIPs.

    double *starting_vector = NULL; // Used for the first variable fixing.

    double *lb_mip = NULL; // Lower bounds of MIP.
    double *ub_mip = NULL; // Upper bounds of MIP.

    double slack_sum; // For OMIP.

    char *in_fn = NULL, *out_fn = NULL; // Input and output file names.
    FILE *out_csv = NULL;               // Csv output file.
    int seed = 0, rho = -1;

    int i, j, tmp, cnt, opt, status;

    int quality_loop = 0;
    double dettime_lim; // Deterministic time parameter.

    int effortlevel = CPX_MIPSTART_REPAIR; // Effort level for the mip starts.

    // Check command line options.
    while ((opt = getopt(argc, argv, "i:o:s:p:")) != -1) {
        switch (opt) {
            case 'i':
                in_fn = optarg;
                break;
            case 'o':
                out_fn = optarg;
                break;
            case 's':
                seed = atoi(optarg);
                break;
            case 'p':
                rho = atoi(optarg);
                if (rho > 100) {
                    rho = 100;
                }
                break;
            default:
                print_usage(argv[0]);
                return 1;
        }
    }

    if (in_fn == NULL || out_fn == NULL || seed == 0 || rho < 0) {
        print_usage(argv[0]);
        return 1;
    }

    // Set seed.
    srand(seed);

    out_csv = fopen(out_fn, "w+");
    if (out_csv == NULL) {
        fprintf(stderr, "The output file %s could not be opened.\n", out_fn);
        goto TERMINATE;
    }

    printf("INPUT FILE: %s\n", in_fn);

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

    // Set node limit.
    //status = CPXsetintparam(env, CPX_PARAM_NODELIM, NODE_LIMIT);
    //if (status) {
    //    fprintf(stderr, "Failure to set node limit, error %d\n", status);
    //    goto TERMINATE;
    //}

    // Set time limit.
    //status = CPXsetdblparam(env, CPXPARAM_TimeLimit, TIME_LIMIT);
    //if (status) {
    //    fprintf(stderr, "Failure to set time limit, error %d\n", status);
    //    goto TERMINATE;
    //}

    // Create MIP from the input file.
    printf("\nCreating MIP.\n");
    mip = CPXcreateprob(env, &status, in_fn);
    if (mip == NULL) {
        fprintf(stderr, "Failed to create MIP.\n");
        goto TERMINATE;
    }

    // Read file and copy the data into the created MIP.
    status = CPXreadcopyprob(env, mip, in_fn, NULL);
    if (status) {
        fprintf(stderr, "Failed to read and copy the problem data (MIP).\n");
        goto TERMINATE;
    }

    // Get number of nonzero elements in MIP.
    numnz_mip = CPXgetnumnz(env, mip);
    if (!numnz_mip) {
        fprintf(stderr, "Failed to read number of nonzero elements of MIP.\n");
        goto TERMINATE;
    }

    // Calculate deterministic time limit and set the parameter.
    dettime_lim = (double) numnz_mip / 100.0;
    if (dettime_lim < 20000) {
        dettime_lim = 20000;
    } else if (dettime_lim > 100000) {
        dettime_lim = 100000;
    }

    status = CPXsetdblparam(env, CPXPARAM_DetTimeLimit, dettime_lim);
    if (status) {
        fprintf(stderr, 
                "Failure to set deterministic time limit, error %d\n", status);
        goto TERMINATE;
    }

    // Create FMIP.
    printf("\nCreating FMIP.\n");
    status = create_fmip(env, mip, &fmip);
    if (status) {
        fprintf(stderr, "Failed to create FMIP.\n");
        goto TERMINATE;
    }

    // Init all the arrays that store info about MIP.
    // Calculate the total ammount of integer variables of MIP.
    numcols_mip = CPXgetnumcols(env, mip);
    for (i = 0, num_int_vars = 0; i < numcols_mip; i++) {
        char type;

        status = CPXgetctype(env, mip, &type, i, i);
        if (status) {
            fprintf(stderr, "Failed to get variable %d type.\n", i);
            goto TERMINATE;
        }

        if (type == CPX_BINARY || type == CPX_INTEGER) {
            num_int_vars += 1;
        }
    }

    // Init: int_indices; fixed_indices; lb_mip; ub_mip.
    int_indices = (int*) malloc(num_int_vars * sizeof(int));
    is_fixed = (int*) malloc(num_int_vars * sizeof(int));

    lb_mip = (double*) malloc(numcols_mip * sizeof(double));
    ub_mip = (double*) malloc(numcols_mip * sizeof(double));

    if (int_indices == NULL ||
        is_fixed == NULL ||
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
        fprintf(stderr, "Failed to initialize arrays with MIP info.\n");
        goto TERMINATE;
    }

    // Space for the initial vector (the starting point).
    starting_vector = (double*) malloc(numcols_mip * sizeof(double));
    if (starting_vector == NULL) {
        fprintf(stderr, "No memory for initial_vector.\n");
        goto TERMINATE;
    }

    // Space for the solution of FMIP.
    numcols_submip = CPXgetnumcols(env, fmip);
    x_fmip = (double*) malloc(numcols_submip * sizeof(double));
    if (x_fmip == NULL) {
        fprintf(stderr, "No memory for solution values for FMIP.\n");
        goto TERMINATE;
    }

    // Space for the solution of OMIP.
    x_omip = (double*) malloc(numcols_submip * sizeof(double));
    if (x_omip == NULL) {
        fprintf(stderr, "No memory for solution values for OMIP.\n");
        goto TERMINATE;
    }

    // Init indices of the variables in the subMIPs.
    varindices_submip = (int*) malloc(numcols_submip * sizeof(int));
    if (varindices_submip == NULL) {
        fprintf(stderr, "No memory for sumMip var indices.\n");
        goto TERMINATE;
    }
    
    for (i = 0; i < numcols_submip; i++) {
        varindices_submip[i] = i;
    }

    // Generate the starting point.
    printf("\n### Generate starting point ###\n");
    status = generate_starting_vector(
        env,
        fmip,
        starting_vector,
        int_indices,
        num_int_vars,
        lb_mip,
        ub_mip,
        BOUND_CONSTANT,
        1
    );
    if (status) {
        fprintf(stderr, "Failed to generate initial vector.\n");
        goto TERMINATE;
    }

    // ACS algorithm.
    for (cnt = 0; cnt < MAX_ITR; cnt++) {
        // Variable fixing on FMIP.
        printf("\n### Variable fixing on FMIP - Iteration %d ###\n", cnt);
        bzero(is_fixed, num_int_vars * sizeof(int));
        if (cnt == 0) { // Use starting vector only in the first iteration.
            status = variable_fixing(
                env,
                fmip,
                int_indices,
                is_fixed,
                num_int_vars,
                starting_vector,
                rho
            );
        } else { // Otherwise use OMIP's solution.
            status = variable_fixing(
                env,
                fmip,
                int_indices,
                is_fixed,
                num_int_vars,
                x_omip,
                rho 
            );
        }
        if (status) {
            fprintf(stderr, "Failed to fix variables of FMIP.\n");
            goto TERMINATE;
        }

        // Optimize FMIP.
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
            // If FMIP is infeasible.
            switch (solstat_fmip) {
                case CPXMIP_NODE_LIM_INFEAS:
                    printf("FMIP is infeasible (Node limit).\n");
                    goto TERMINATE;
                case CPXMIP_TIME_LIM_INFEAS:
                    printf("FMIP is infeasible (Time limit).\n");
                    goto TERMINATE;
                case CPXMIP_DETTIME_LIM_INFEAS:
                    printf("FMIP is infeasible (DetTime limit).\n");
                    goto TERMINATE;
                default:
                    fprintf(stderr, "Failed to optimize FMIP "
                                    "(Unknown status)\n"
                                    "Solution status: %d\n", solstat_fmip);
                    goto TERMINATE;
            }
        }

        // Restore bounds of FMIP.
        status = restore_bounds(
            env,
            fmip,
            lb_mip,
            ub_mip,
            int_indices,
            is_fixed,
            num_int_vars
        );
        if (status) {
            fprintf (stderr, "Failed to restore FMIP bounds.\n");
            goto TERMINATE;
        }

        // If the FMIP's solution is feasibile.
        switch (solstat_fmip) {
            case CPXMIP_OPTIMAL:
                printf("Found a feasibile solution for FMIP (Optimal).\n");
                break;
            case CPXMIP_OPTIMAL_TOL:
                printf("Found a feasibile solution for FMIP (Optimal tollerance).\n");
                break;
            case CPXMIP_NODE_LIM_FEAS:
                printf("Found a feasibile solution for FMIP (Node limit).\n");
                break;
            case CPXMIP_TIME_LIM_FEAS:
                printf("Found a feasibile solution for FMIP (Time limit).\n");
                break;
            case CPXMIP_DETTIME_LIM_FEAS:
                printf("Found a feasibile solution for FMIP (DetTime limit).\n");
                break;
        }

        // Create OMIP in the first iteration.
        if (cnt == 0) {
            printf("\nCreating OMIP.\n");
            status = create_omip(env, mip, &omip, objval_fmip);
            if (status) {
                fprintf(stderr, "Failed to create OMIP.\n");
                goto TERMINATE;
            }
        } else { // Update the slack constraint of OMIP.
            tmp = CPXgetnumrows(env, omip) - 1;
            status = CPXchgrhs(env, omip, 1, &tmp, &objval_fmip);
            if (status) {
                fprintf(stderr, "Failed to update rhs value (slack) of OMIP.\n");
                goto TERMINATE;
            }
        }

        if (cnt != 0) {
            // Remove all mip starts from OMIP.
            status = CPXdelmipstarts(
                env,
                omip,
                0,
                CPXgetnummipstarts(env, omip) - 1
            );
            if (status) {
                fprintf(stderr, "Error deleting MIP starts.\n");
                goto TERMINATE;
            }
        }

        // Add the FMIP solution as a mip start for OMIP.
        tmp = 0;
        status = CPXaddmipstarts(
            env,
            omip,
            1,
            numcols_submip,
            &tmp,
            varindices_submip,
            x_fmip,
            &effortlevel,
            NULL
        );
        if (status) {
            fprintf(stderr, "Failed to add mip start for OMIP.\n");
            goto TERMINATE;
        }
        

        // TODO: Remove the for loop below. Using mip starts avoid the
        //       infeasability caused by the deterministic time limit.
        // Try to solve OMIP.
        for (i = 0; i < MAX_ATTEMPTS; i++) {
            // Variable fixing on OMIP.
            printf("\n### Variable fixing on OMIP "
                   "- Attempt %d - Iteration %d ###\n", i, cnt);
            bzero(is_fixed, num_int_vars * sizeof(int));
            status = variable_fixing(
                env,
                omip,
                int_indices,
                is_fixed,
                num_int_vars,
                x_fmip,
                rho
            );
            if (status) {
                fprintf(stderr, "Failed to fix variables of OMIP.\n");
                goto TERMINATE;
            }

            // Optimize OMIP.
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
                // If OMIP is infeasible.
                switch (solstat_omip) {
                    case CPXMIP_NODE_LIM_INFEAS:
                        printf("OMIP is infeasible (Node limit).\n");
                        break;
                    case CPXMIP_TIME_LIM_INFEAS:
                        printf("OMIP is infeasible (Time limit).\n");
                        break;
                    case CPXMIP_DETTIME_LIM_INFEAS:
                        printf("OMIP is infeasible (DetTime limit).\n");
                        break;
                    default:
                        fprintf(stderr, "Failed to optimize OMIP (Unknown)\n"
                                        "Solution status: %d\n", solstat_fmip);
                        goto TERMINATE;
                }
            }

            // Restore bounds of OMIP.
            status = restore_bounds(
                env,
                omip,
                lb_mip,
                ub_mip,
                int_indices,
                is_fixed,
                num_int_vars
            );
            if (status) {
                fprintf (stderr, "Failed to restore OMIP bounds.\n");
                goto TERMINATE;
            }

            // If the OMIP's solution is feasibile.
            switch (solstat_omip) {
                case CPXMIP_OPTIMAL:
                    printf("Found a feasibile solution for OMIP (Optimal).\n");
                    break;
                case CPXMIP_OPTIMAL_TOL:
                    printf("Found a feasibile solution for OMIP (Optimal tollerance).\n");
                    break;
                case CPXMIP_NODE_LIM_FEAS:
                    printf("Found a feasibile solution for OMIP (Node limit).\n");
                    break;
                case CPXMIP_TIME_LIM_FEAS:
                    printf("Found a feasibile solution for OMIP (Time limit).\n");
                    break;
                case CPXMIP_DETTIME_LIM_FEAS:
                    printf("Found a feasibile solution for OMIP (DetTime limit).\n");
                    break;
                default:
                    continue; // Try to solve another OMIP.
            }

            // Check sum of slack variables.
            slack_sum = sum(x_omip, numcols_mip, numcols_submip - 1);
            printf("Slack sum: %f\n", slack_sum);
            if (slack_sum == 0) {
                printf("Found a feasibile solution for MIP.\n");
                quality_loop = 1; // Quit ACS.

                // Save results to the csv output file.
                fprintf(out_csv, "fmip,%f\nomip,%f\n", slack_sum, objval_omip);
            } else {
                // Save results to the csv output file.
                fprintf(out_csv, "fmip,%f\nomip,%f\n", objval_fmip, objval_omip);
            }

            break; // Only if OMIP's solution is feasibile.
        }

        if (MAX_ATTEMPTS == i) { // If no OMIP has been resolved.
            printf("All OMIPs were infeasibile.\n");
            break;
        } else if (quality_loop) { // If a (MIP) feasibile solution was found.
            break;
        }

        // Remove all mip starts from FMIP.
        status = CPXdelmipstarts(
            env,
            fmip,
            0,
            CPXgetnummipstarts(env, fmip) - 1
        );
        if (status) {
            fprintf(stderr, "Error deleting MIP starts.\n");
            goto TERMINATE;
        }

        // Add OMIP solution as a mip start for FMIP.
        tmp = 0;
        status = CPXaddmipstarts(
            env,
            fmip,
            1,
            numcols_submip,
            &tmp,
            varindices_submip,
            x_omip,
            &effortlevel,
            NULL
        );
        if (status) {
            fprintf(stderr, "Failed to add mip start for FMIP.\n");
            goto TERMINATE;
        }
    }

    // Improve the quality of the solution.
    if (quality_loop) {
        // Update the slack constraint of OMIP, if necessary.
        // (slack_sum could be 0 while objval_fmip could be grater than 0)
        if (slack_sum != objval_fmip) {
            tmp = CPXgetnumrows(env, omip) - 1;
            status = CPXchgrhs(env, omip, 1, &tmp, &slack_sum);
            if (status) {
                fprintf(stderr, "Failed to update rhs value (slack) of OMIP.\n");
                goto TERMINATE;
            }
        }

        printf("\n### Improve feasibile solution quality ###");
        int obj_sen = CPXgetobjsen(env, omip);
        double prev_objval = objval_omip;
        for (i = 0; ; i++) {
            printf("\nIteration: %d\n", i);

            bzero(is_fixed, num_int_vars * sizeof(int));

            status = variable_fixing(
                env,
                omip,
                int_indices,
                is_fixed,
                num_int_vars,
                x_omip,
                rho
            );
            if (status) {
                fprintf(stderr, "Failed to fix variables of OMIP.\n");
                goto TERMINATE;
            }

            // Optimize OMIP.
            optimize_prob(
                env,
                omip,
                &objval_omip,
                &solstat_omip,
                x_omip,
                0,
                numcols_submip - 1,
                1
            );

            // Restore bounds of OMIP.
            status = restore_bounds(
                env,
                omip,
                lb_mip,
                ub_mip,
                int_indices,
                is_fixed,
                num_int_vars
            );
            if (status) {
                fprintf (stderr, "Failed to restore OMIP bounds.\n");
                goto TERMINATE;
            }

            // If the OMIP's solution is feasibile.
            switch (solstat_omip) {
                case CPXMIP_OPTIMAL:
                    printf("Found a feasibile solution for OMIP (Optimal).\n");
                    break;
                case CPXMIP_OPTIMAL_TOL:
                    printf("Found a feasibile solution for OMIP (Optimal tollerance).\n");
                    break;
                case CPXMIP_NODE_LIM_FEAS:
                    printf("Found a feasibile solution for OMIP (Node limit).\n");
                    break;
                case CPXMIP_TIME_LIM_FEAS:
                    printf("Found a feasibile solution for OMIP (Time limit).\n");
                    break;
                case CPXMIP_DETTIME_LIM_FEAS:
                    printf("Found a feasibile solution for OMIP (DetTime limit).\n");
                    break;
                default:
                    goto TERMINATE; // Quit.
            }

            // Check if the new solution is better.
            if ((obj_sen == CPX_MIN && objval_omip < prev_objval) ||
                (obj_sen == CPX_MAX && objval_omip > prev_objval)) {
                prev_objval = objval_omip;
            } else { // Otherwise save the solution and quit.
                objval_omip = prev_objval;
                fprintf(out_csv, "quality,%f\n", objval_omip);
                break;
            }
        }
    }

TERMINATE:
    
    // Close output file, if necessary.
    if (out_csv != NULL) {
        fclose(out_csv);
    }

    // Free up all the arrays.
    free_and_null((char**) &x_fmip);
    free_and_null((char**) &x_omip);
    free_and_null((char**) &int_indices);
    free_and_null((char**) &is_fixed);
    free_and_null((char**) &starting_vector);
    free_and_null((char**) &lb_mip);
    free_and_null((char**) &ub_mip);
    free_and_null((char**) &varindices_submip);

    // Free up the problems allocated by CPXcreateprob, if necessary.
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
        status = CPXfreeprob(env, &omip);
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
