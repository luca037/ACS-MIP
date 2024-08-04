#include <ilcplex/cplex.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>

// Lunghezza massima del nome delle variabili slack (compreso il terminatore
// di stringa).
// Delta+ e Delta- sono due vettori di lunghezza uguale al numero delle 
// righe della matrice dei vincoli.
// Le loro componenti hanno i segunti nomi:
//      Delta+ = [dp1, dp2, dp3, ..., dpi, ...,]
//      Delta- = [dn1, dn2, dn3, ..., dni, ...,]
// Se 9 è la lunghezza massima allora 9-3 è il numero di cifre massime
// che può avere il nome delle variabili (-3 per la presenza del carattere
// 'd', del carattere 'p' e di '\0').
#define MAX_SLACK_NAMES_LEN 9

// Lunghezza massima nome variabili e vincoli.
#define MAX_COLNAME_LEN 9

// Numero di tentativi che vengono effettuati per risolvere un problema che
// risolta non risolvibile.
#define MAX_ATTEMPTS 2000

// Massimo numero nodi esplorabili in un problema.
#define NODE_LIMIT 5000

// Massimo numero di run.
#define MAX_RUN 50


/**
 * Print help message.
 *
 * progname: The executable's name.
 */
void print_usage(char *progname) {
   fprintf (stderr, "Usage: %s <filename>\n"
                    "   where filename is a file with extension \n"
                    "      MPS, SAV, or LP (lower case is allowed)\n"
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


// Ritorna la somma degli elementi dell'array passato.
// La lunghezza di 'x' deve essere almeno 'end' - 'beg' + 1.
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
 * Fixes some of the problem's variables to a value specified by the array
 * named x. The index of the variable that will be fixed is randomly choose
 * from the indices specifed in the array named int_indices.
 *
 * For example:
 *      int_indices = {1, 4, 20}, x = {100, 49, 66}, fixed_indices = {0, 0, 0}
 *      
 *      If the random generetor generates the value 2 then the variable
 *      with index int_indices[2] is fixed to the value x[int_indices[2]].
 *      Then in fixed_indices[2] will be stored 1.
 *
 *      In this our example the variable with index 20 will be fixed to 66
 *      and fixed_indices = {0, 0, 1}.
 *
 * env: A pointer to the CPLEX environment.
 * lp: A pointer to a CPLEX problem object.
 * int_indices: An array where are stored the indices of the integer variables
 *              of lp.
 *              This array must be of length at least num_int_vars.
 * fixed_indices: An array that tells if a variable was fixed.
 *                A variable with index int_indices[i] was fixed if and only if 
 *                fixed_indices[i] is a nonzero value.
 *                This array must be of length at least num_int_vars.
 * num_int_vars: The length of the arrays int_indices and fixed_indices.
 * x: An array where are stored the values to use for fixing the variables.
 *    This array must be of length at least num_int_vars.
 *
 * return: Returns 0 if successful and nonzero if an error occurs.
 */
int variable_fixing(
    CPXENVptr env,
    CPXLPptr lp,
    int *int_indices,
    int *fixed_indices,
    int num_int_vars,
    double *x
) {
    int i, cnt, rnd, status;
    double val; // Valore a cui fisso la var scelta.
    char lu = 'B';

    cnt = num_int_vars / 5; // Numero di variabili da fissare.
    printf("Variables to fix: %d\n", cnt);
    for (i = 0; i < cnt; ) {
        // Genero una posizione random di 'index'.
        rnd = rand() % num_int_vars;
        // Se è già stata fissata la variabile genero un altro valore.
        if (fixed_indices[rnd] == 1) {
            continue;
        }

        val = x[int_indices[rnd]];

        // Fisso il valore della variabile.
        status = CPXchgbds(env, lp, 1, &int_indices[rnd], &lu, &val);
        if (status) {
            fprintf(stderr, 
                    "Failed to fix variable with index %d.\n", int_indices[rnd]);
            return status;
        }
        fixed_indices[rnd] = 1;

        i += 1;
    }

    return 0;
}


// Ottimizza il problema passato e scrive i risulati nelle variabili passate.
int optimize_prob(
    CPXENVptr env,
    CPXLPptr lp,
    double *objval,
    int *solstat,
    double *x,
    int begin,
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
        return status;
    }

    // Get variables values.
    status = CPXgetx(env, lp, x, begin, end);
    if (status) {
        fprintf(stderr, "Failed to obtain solution.\n");
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


// Aggiunge le variabili slack alla matrice dei vincoli di lp.
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


// Copia i dati del problema src nel problema dst.
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


// Permette di creare il problema FMIP partendo dal problema MIP fornito.
// Non effettua il variable fixing.
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


// Permette di creare l'OMIP partendo dal problema MIP fornito.
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
    CPXENVptr env = NULL;
    CPXLPptr mip = NULL, fmip = NULL, omip = NULL;

    // Variabili per le accedere alle soluzioni dei problemi.
    int numcols_mip, numcols_submip, solstat_fmip, solstat_omip;
    double objval_fmip, objval_omip; // Objective value.
    double *x_fmip = NULL, *x_omip = NULL; // Variabiles value.

    // Variabili per il variable fixing.
    int *int_indices = NULL; // Indici delle variabili intere.
    int *fixed_indices = NULL; // Indici delle variabili che sono state fissate.
    int num_int_vars; // Dimensione dei due array.
    double *lb_mip = NULL, *ub_mip = NULL; // Necessari per resettare i lb e ub
                                           // dopo un variable fixing.
    double slack_sum; // Somma variabili slack (per OMIP).

    int i, j, tmp, cnt, status;

    // Setto il seme.
    srand(5);

    // Check command line arguments.
    if (argc != 2) {
        print_usage(argv[0]);
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

    // ### Creazione del MIP, FMIP e OMIP ###
    printf("\nCreating MIP.\n");
    mip = CPXcreateprob(env, &status, argv[1]);
    if (mip == NULL) {
        fprintf(stderr, "Failedto create MIP.\n");
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
            if (cnt == 0) { // Genero l'initial vector al primo ciclo.
                status = variable_fixing(
                    env,
                    fmip,
                    int_indices,
                    fixed_indices,
                    num_int_vars, 
                    ub_mip
                );
            } else { // Altrimenti uso i valori di OMIP.
                status = variable_fixing(
                    env,
                    fmip,
                    int_indices,
                    fixed_indices,
                    num_int_vars,
                    x_omip
                );
            }
            if (status) {
                fprintf(stderr, "Failed to fix variables of FMIP.\n");
                goto TERMINATE;
            }

            // Stampa variabili fissate.
            //printf("Variabili fissate:\n");
            //for (j = 0; j < num_int_vars; j++) {
            //    char name[MAX_COLNAME_LEN];
            //    if (fixed_indices[j]) {
            //        get_colname(env, fmip, int_indices[j], name);
            //        printf("Var fissata -> %s\n", name);
            //    }
            //}
            //printf("\n");

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
                if (solstat_fmip == CPXMIP_INFEASIBLE) {
                    printf("FMIP is infeasible.\n");
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
                    continue;
                } else {
                    fprintf(stderr, "Failed to optimize FMIP.\n"
                                    "Solution status: %d", solstat_fmip);
                    goto TERMINATE;
                }
            }

            // Se ho trovato una soluzione ottima.
            if (solstat_fmip == CPXMIP_OPTIMAL ||
                solstat_fmip == CPXMIP_OPTIMAL_TOL ||
                solstat_fmip == CPXMIP_NODE_LIM_FEAS) {
                printf("Found an optimal solution for FMIP.\n");
                break;
            } 
        }

        // Se nessun FMIP è stato risolto esco.
        if (solstat_fmip == CPXMIP_INFEASIBLE) {
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
                x_fmip
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
                if (solstat_omip == CPXMIP_INFEASIBLE) {
                    printf("OMIP is infeasible.\n");
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
                    continue;
                } else {
                    fprintf(stderr, "Failed to optimize OMIP."
                                    "\n Solution status: %d", solstat_omip);
                    goto TERMINATE;
                }
            }

            // Se ho trovato una soluzione ottima.
            slack_sum = sum(x_omip, numcols_mip, numcols_submip - 1);
            printf("Slack sum: %.2f\n", slack_sum);
            if (solstat_omip == CPXMIP_OPTIMAL ||
                solstat_omip == CPXMIP_OPTIMAL_TOL ||
                solstat_omip == CPXMIP_NODE_LIM_FEAS) {
                printf("Found an optimal solution for OMIP.\n");
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
        if (solstat_fmip == CPXMIP_INFEASIBLE) {
            printf("All OMIPs generated was infeasibile.\n");
            break;
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
