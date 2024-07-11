/*
 * FMIP si ottiene partendo dal MIP:
 *      - aggiungendo due colonne (variabili slack) alla matrice dei vincoli
 *      - fissando il valore di alcune variabili xi (aggiungendo righe)
 *      - aggiungendo i bounud per gli slack (entrambi >= 0)
 *      - modificando la funzione obiettivo (minizzando la somma delle var
 *        slack).
 *
 * OMIP si ottiene partendo dal MIP:
 *      - aggiungendo due colonne (variabili slack) alla matrice dei vincoli
 *        come nell'FMIP
 *      - fissando il valore di alcune variabili xi (aggiungendo righe)
 *        (come nell'FMIP ma non le stesse var)
 *      - aggiungendo i bounud per gli slack (entrambi >= 0)
 *      - aggiungendo il vincolo per cui la somma delle variabili
 *        slack sia minore della somma delle variabili di slack
 *        fornita in input
 */

#include <ilcplex/cplex.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Lunghezza del nome delle variabili slack.
// Delta+ e Delta- sono due vettori di lunghezza uguale al numero delle 
// righe della matrice dei vincoli.
// Le loro componenti hanno i segunti nomi:
//      Delta+ = [dp1, dp2, dp3, ..., dpi, ...,]
//      Delta- = [dn1, dn2, dn3, ..., dni, ...,]
// Tutti i nomi delle componenti sono lunghi 3 + 1 (il +1 per il terminatore
// di stringa).
#define SLACK_NAMES_LEN 4

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

// Aggiunge le variabili slack alla matrice dei vincoli di lp.
int add_slack_cols(CPXENVptr env, CPXLPptr lp) {
    int i, status = 0;

    int numrows = CPXgetnumrows(env, lp);
    int numcols = CPXgetnumcols(env, lp);

    int ccnt = 2 * numrows; // Numero di colonne (slack) da aggiungere.
                            // I vettori slack hanno dimensione pari al
                            // numero di righe della matrice dei vincoli.
    int nzcnt = 2 * numrows * numrows; // Non zero coefficients counter.
                                       // Tutti i coefficienti sono diversi 
                                       // da zero.

    double *obj = NULL, *lb = NULL, *ub = NULL, *matval = NULL;
    int *matbeg = NULL, *matind = NULL;
    char **colnames = NULL;

    obj = (double*) malloc(ccnt * sizeof(double));
    lb = (double*) malloc(ccnt * sizeof(double));
    ub = (double*) malloc(ccnt * sizeof(double));

    matbeg = (int*) malloc(ccnt * sizeof(int));
    matind = (int*) malloc(nzcnt * sizeof(int));
    matval = (double*) malloc(nzcnt * sizeof(double));

    colnames = (char**) malloc(ccnt * sizeof(char*));

    if (obj == NULL || lb == NULL || ub == NULL || 
        matbeg == NULL || matind == NULL || matval == NULL || colnames == NULL) {
        fprintf(stderr, "No memory for adding slack variabiles.\n");
        status = 1;
        goto TERMINATE;
    }

    // Alloco spazio per i nomi.
    for (i = 0; i < ccnt; i++) {
        colnames[i] = (char*) malloc(SLACK_NAMES_LEN * sizeof(char));
        if (colnames[i] == NULL) {
            fprintf(stderr, "No memory for slack variabiles names.\n");
            status = 1;
            goto TERMINATE;
        }
    }

    // Assegno i nomi alle variabili slack.
    for (i = 0; i < numrows; i++) {
        snprintf(colnames[i], SLACK_NAMES_LEN, "dp%d", i + 1);
        snprintf(colnames[i + numrows], SLACK_NAMES_LEN, "dn%d", i + 1);
    }

    // Inizializzo: coefficienti nella funzione obiettivo; lowerbound; 
    // upperbound; matbeg.
    for (i = 0; i < ccnt; i++) {
        obj[i] = 0;
        lb[i] = 0;
        ub[i] = CPX_INFBOUND;
        matbeg[i] = i * numrows;
    }

    // Inizializzo: matind; matval.
    for (i = 0; i < nzcnt; i++) {
        // matind = [0,...,m-1,0,...,m-1,...] dove m = numrows.
        matind[i] = i % numrows; 
        // matval ha i primi nzcnt/2 valori a 1 e i restanti nzcnt/2 a -1.
        matval[i] = (i < nzcnt / 2)? 1 : -1;
    }

    // Ora posso aggiungere le due colonne di slack.
    status = CPXaddcols(env, lp, ccnt, nzcnt, NULL, matbeg, matind, matval, lb, ub, colnames);
    if (status) {
        fprintf(stderr, "Failed to add new columns to the problem.\n");
        goto TERMINATE;
    }

    // Aggiorno il numero di colonne della matrice dei vincoli.
    numcols = CPXgetnumcols(env, lp);

    // Setto il tipo delle slack.
    char ctype[] = {CPX_CONTINUOUS, CPX_CONTINUOUS};
    int indices[] = {numcols - 2, numcols - 1};

    status = CPXchgctype(env, lp, 2, indices, ctype);
    if (status) {
        fprintf(stderr, "Failed to change ctype.\n");
        goto TERMINATE;
    }

TERMINATE:

    free_and_null((char**) &obj);
    free_and_null((char**) &lb);
    free_and_null((char**) &ub);
    free_and_null((char**) &matbeg);
    free_and_null((char**) &matind);
    free_and_null((char**) &matval);
    for (i = 0; i < ccnt; i++) {
        free_and_null(&colnames[i]);
    }
    free_and_null(colnames);

    return status;
}

// Permette di creare il problema FMIP partendo dal suo MIP.
int create_fmip(CPXENVptr env, CPXLPptr mip, CPXLPptr *fmip) {
    int i, cnt, status = 0;
    double tmp;
    
    int numcols = CPXgetnumcols(env, mip);
    int numrows = CPXgetnumrows(env, mip);
    int numnz = CPXgetnumnz(env, mip);

    double *obj = NULL, *rhs = NULL, *matval = NULL, *ub = NULL, *lb = NULL;
    char *sense = NULL, *xctype = NULL;
    int *matbeg = NULL, *matind = NULL, *matcnt = NULL;

    int nzcnt, surplus;

    *fmip = CPXcreateprob(env, &status, "FMIP");
    if (*fmip == NULL) {
        fprintf(stderr, "Failed to create FMIP.\n");
        status = 1;
        goto TERMINATE;
    }

    // Salvo i coefficienti della funzione obiettivo del problema MIP.
    obj = (double*) malloc(numcols * sizeof(double));
    if (obj == NULL) {
        fprintf(stderr, "No memory for saving obj of MIP.\n");
        goto TERMINATE;
    }

    status = CPXgetobj(env, mip, obj, 0, numcols - 1);
    if (status) {
        fprintf(stderr, "Failed to copy obj coefficients of MIP into array.\n");
        goto TERMINATE;
    }

    //printf("Coefficienti del problema %s:\n", argv[1]);
    //for (i = 0; i < cur_numcols; i++) {
    //    printf("x[%d] = %0.3g\n", i, obj[i]);
    //}

    // Buffer in cui salvo i coefficienti rhs dei vincoli del problema MIP.
    rhs = (double*) malloc(numrows * sizeof(double));
    if (rhs == NULL) {
        fprintf(stderr, "No memory for saving rhs of MIP.\n");
        status = 1;
        goto TERMINATE;
    }

    status = CPXgetrhs(env, mip, rhs, 0, numrows - 1);
    if (status) {
        fprintf(stderr, "Failed to copy rhs coefficients of MIP into array.\n");
        goto TERMINATE;
    }

    //printf("Coefficienti rhs vincoli del problema %s:\n", argv[1]);
    //for (i = 0; i < cur_numrows; i++) {
    //    printf("%0.2f\n", rhs[i]);
    //}

    // Salvo in un buffer il senso dei vincoli del problema MIP.
    sense = (char*) malloc(numrows * sizeof(double));
    if (sense == NULL) {
        fprintf(stderr, "No memory for saving sense of MIP.\n");
        status = 1;
        goto TERMINATE;
    }

    status = CPXgetsense(env, mip, sense, 0, numrows - 1);
    if (status) {
        fprintf(stderr, "Failed to copy sense coefficients of MIP into array.\n");
        goto TERMINATE;
    }

    //printf("Coefficienti sense vincoli del problema %s:\n", argv[1]);
    //for (i = 0; i < cur_numrows; i++) {
    //    printf("%c\n", sense[i]);
    //}

    // Salvo i dati della matrice dei vincoli di MIP.
    matbeg = (int*) malloc(numcols * sizeof(int));
    matind = (int*) malloc(numnz * sizeof(int));
    matcnt = (int*) malloc(numcols * sizeof(int));
    matval = (double*) malloc(numnz * sizeof(double));

    if (matbeg == NULL || matind == NULL || matval == NULL || matcnt == NULL) {
        fprintf(stderr, "No memory for saving matbeg or matind or matval or matcnt of MIP.\n");
        status = 1;
        goto TERMINATE;
    }

    status = CPXgetcols(env, mip, &nzcnt, matbeg, matind, matval, numnz, &surplus, 0, numcols - 1);
    if (status) {
        fprintf(stderr, "Failed to copy matbeg or matind or matval coefficients of MIP.\n");
        goto TERMINATE;
    }

    // Calcolo delle componenti di matcnt.
    for (i = 0; i < numcols - 1; i++) {
        matcnt[i] = matbeg[i+1] - matbeg[i];
    }
    matcnt[numcols - 1] = nzcnt - matbeg[numcols - 1];

    //printf("matbeg:\n");
    //for (i = 0; i < cur_numcols; i++) {
    //    printf("%d, ", matbeg[i]);
    //}
    //printf("\n");

    //printf("matind:\n");
    //for (i = 0; i < cur_numnz; i++) {
    //    printf("%d, ", matind[i]);
    //}
    //printf("\n");

    //printf("matval:\n");
    //for (i = 0; i < cur_numnz; i++) {
    //    printf("%0.2f, ", matval[i]);
    //}
    //printf("\n");

    //printf("matcnt:\n");
    //for (i = 0; i < cur_numcols; i++) {
    //    printf("%d, ", matcnt[i]);
    //}
    //printf("\n");

    // Salvo i lower bound e upper bound delle variabili.
    ub = (double*) malloc(numcols * sizeof(double));
    lb = (double*) malloc(numcols * sizeof(double));
    if (ub == NULL || lb == NULL) {
        fprintf(stderr, "No memory for saving lower bounds or upper bounds of MIP.\n");
        status = 1;
        goto TERMINATE;
    }

    status = CPXgetub(env, mip, ub, 0, numcols - 1);
    if (status) {
        fprintf(stderr, "Failed to copy upper bounds coefficients of MIP.\n");
        goto TERMINATE;
    }

    status = CPXgetlb(env, mip, lb, 0, numcols - 1);
    if (status) {
        fprintf(stderr, "Failed to copy upper bounds coefficients of MIP.\n");
        goto TERMINATE;
    }

    //printf("Coefficienti upper bounds del problema %s\n", argv[1]);
    //for (i = 0; i < cur_numcols; i++) {
    //    printf("%.2f\n", ub[i]);
    //}
    //printf("Coefficienti lower bounds del problema %s\n", argv[1]);
    //for (i = 0; i < cur_numcols; i++) {
    //    printf("%.2f\n", lb[i]);
    //}
    
    // Salvo la tipologia delle variabili del problema MIP.
    xctype = (char*) malloc(numcols * sizeof(char));
    if (xctype == NULL) {
        fprintf(stderr, "No memory for saving variables types of MIP.\n");
        status = 1;
        goto TERMINATE;
    }

    status = CPXgetctype(env, mip, xctype, 0, numcols - 1);
    if (status) {
        fprintf(stderr, "Failed to copy upper bounds coefficients of MIP.\n");
        goto TERMINATE;
    }

    // Ora posso crere la copia del problema MIP in FMIP.
    // TODO: non gestisco rngval (settato a NULL).
    status = CPXcopylp(
            env, *fmip, numcols, numrows, CPX_MIN, 
            obj, rhs, sense, matbeg, matcnt, matind, matval, lb, ub, NULL
    );
    if (status) {
        fprintf(stderr, "Failed to populate FMIP from MIP.\n");
        goto TERMINATE;
    }

    // Setto la tipologia delle variabili del problema FMIP.
    status = CPXcopyctype(env, *fmip, xctype);
    if (status) {
        fprintf(stderr, "Failed to set variables types of FMIP.\n");
        goto TERMINATE;
    }

    // Introduco le variabili slack.
    status = add_slack_cols(env, *fmip);
    if (status) {
        fprintf(stderr, "Failed to add slack columns to FMIP.\n");
        goto TERMINATE;
    }
    // Nota che il numero di colonne di FMIP non viene aggiornato dopo
    // aver aggiunto le variabili slack. Ciò significa che numcols indica
    // il numero di colonne senza considerare le variabili di slack.

    // Modifico la funzone obiettivo.
    // Prima setto i coeffcienti delle variabili xi a zero.
    for (i = 0, tmp = 0; i < numcols; i++) {
        status = CPXchgobj(env, *fmip, 1, &i, &tmp);
        if (status) {
            fprintf(stderr, "Failed change variable index %d coefficient in obj.\n", i);
            goto TERMINATE;
        }
    }
    // Ora setto i coefficienti delle variabili slack a 1.
    // Le variabili slack in totale sono 2 * numrows.
    // Il range di indici delle variabili slack è:
    //      [numcols, numcols + 2 * numrows]
    cnt = numcols + 2 * numrows;
    for (i = numcols, tmp = 1; i < cnt; i++) {
        status = CPXchgobj(env, *fmip, 1, &i, &tmp);
        if (status) {
            fprintf(stderr, "Failed change variable index %d (slack) coefficient in obj.\n", i);
            goto TERMINATE;
        }
    }

    // Fisso il valore di alcune variabili xi.
    // Esempio semplice: ne fisso metà partendo dall'inizio (suppongo che
    // queste var siano tutte intere). Fisso il valore delle variabili al 
    // loro upperbound.
// ---------------------------------------------------------------------------
    // TODO: scrivere l'algoritmo che prende la decisione sulle variabili
    // da fissare. Quello qui sotto è solo un esempio per testare.
    cnt = numcols / 2;
    int *indices = (int*) malloc(cnt * sizeof(int));
    char *lu= (char*) malloc(cnt * sizeof(char));
    double *bd= (double*) malloc(cnt * sizeof(double));

    for (i = 0; i < cnt; i++) {
        indices[i] = i;
        lu[i] = 'B';
        status = CPXgetub(env, *fmip, &bd[i], i, i);
        if (status) {
            fprintf(stderr, "Failed get ub of variabile index %d.\n", i);
            goto TERMINATE;
        }
    }

    status = CPXchgbds(env, *fmip, cnt, indices, lu, bd);
    if (status) {
        fprintf(stderr, "Failed to fix variables of FMIP.\n");
        goto TERMINATE;
    }
// ---------------------------------------------------------------------------


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

    free_and_null((char**) &indices);
    free_and_null(&lu);
    free_and_null((char**) &bd);

    return status;
}

int main(int argc, char* argv[]) {
    int solstat; // Solution status.
    //double objval; // Objective value.
    double *x_mip1 = NULL, *x_mip2 = NULL; // Variabiles value.

    CPXENVptr env = NULL;
    CPXLPptr mip = NULL, fmip = NULL, omip = NULL;
    int status;
    int i;
    int cur_numcols, cur_numrows, cur_numnz;

    // Check command line arguments.
    if (argc != 2) {
        fprintf(stderr, "Missing input file.\nUsage: %s <file1>\n", argv[0]);
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

// ### Creazione del MIP originale (file di input) ###
    printf("Creazione del mip originale.\n");
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

    // Optimize MIP and print results
    //status = optimze_and_print_results(env, mip, x_mip1);
    //if (status) {
    //    fprintf(stderr, "Failed to optimize MIP.\n");
    //    goto TERMINATE;
    //}

// ### Creazione dell'FMIP ###
    printf("\nCreazione dell'FMIP.\n");
    status = create_fmip(env, mip, &fmip);
    if (status) {
        fprintf(stderr, "Failed to create FMIP.\n");
        goto TERMINATE;
    }

    // Optimize fmip and print results.
    //status = optimze_and_print_results(env, fmip, x_mip2);
    //if (status) {
    //    fprintf(stderr, "Failed to optimize fmip.\n");
    //    goto TERMINATE;
    //}

    // Salvo il problema FMIP in un file.
    status = CPXwriteprob(env, fmip, "fmip.lp", NULL);
    if (status) {
        fprintf (stderr, "Failed to write LP to disk.\n");
    }

TERMINATE:

    // Free up solution.
    free_and_null((char**) &x_mip1);
    free_and_null((char**) &x_mip2);

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
