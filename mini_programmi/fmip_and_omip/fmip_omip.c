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

// Lunghezza massima nome variabili.
#define MAX_COLNAME_LEN 9

// Frees up pointer *ptr and sets it to NULL.
void free_and_null(char** ptr) {
    if (*ptr != NULL) {
        free(*ptr);
        *ptr = NULL;
    }
}

// Ritorna la somma degli elementi dell'array passato.
double sum(double *x, int len) {
    if (len == 0) return 0;
    double s = x[0];
    for (int i = 1; i < len; i++) {
        s += x[i];
    }

    return s;
}

// Permette di ottenere il nome di una variabile, dato il suo indice.
// Salva il nome della variabile in colname.
int get_colname(CPXENVptr env, CPXLPptr lp, int index, char *colname) {
    int surplus, status = 0;
    char namestore[MAX_COLNAME_LEN];

    // Ricavo il nome della colonna.
    status = CPXgetcolname(env, lp, (char**) &namestore, namestore, MAX_COLNAME_LEN, &surplus, index, index);
    if (status) {
        fprintf(stderr, "Failed get column name.\n");
        return status;
    }

    strncpy(colname, namestore, MAX_COLNAME_LEN);

    return status;
}


// Ottimizza il problema passato e scrive i risulati nelle variabili passate.
int optimize_prob(CPXENVptr env, CPXLPptr lp, double *objval, int *solstat, double *x, int begin, int end) {
    int status = 0;
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
        fprintf(stderr, "Failed to obtain objective value (lp).\n");
        return status;
    }

    // Print variabile values.
    status = CPXgetx(env, lp, x, begin, end);
    if (status) {
        fprintf(stderr, "Failed to obtain solution.\n");
        return status;
    }

    // Stampo il valore delle variabili.
    char colname[MAX_COLNAME_LEN];
    for (int i = 0; i < numcols; i++) { 
        if (get_colname(env, lp, i, colname) == 0) {
            printf("Column = %d,\tValue = %.2f,\tVar name = %s\n", i, x[i], colname);
        } else {
            printf("Column = %d,\tValue = %.2f,\tVar name = *error*\n", i, x[i]);
        }
    }
    
    return status;
}


// Aggiunge le variabili slack alla matrice dei vincoli di lp.
int add_slack_cols(CPXENVptr env, CPXLPptr lp) {
    int i, tmp, status = 0;
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

    if (matbeg == NULL || matind == NULL || matval == NULL || colnames == NULL) {
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
    status = CPXaddcols(env, lp, ccnt, nzcnt, NULL, matbeg, matind, matval, NULL, NULL, colnames);
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

    return status;
}


// Copia i dati del problema src nel problema dst.
int copy_prob(CPXENVptr env, CPXLPptr src, CPXLPptr dst) {
    int i, status = 0;

    int numcols = CPXgetnumcols(env, src);
    int numrows = CPXgetnumrows(env, src);
    int numnz = CPXgetnumnz(env, src);

    double *obj = NULL, *rhs = NULL, *matval = NULL, *ub = NULL, *lb = NULL;
    char *sense = NULL, *xctype = NULL;
    int *matbeg = NULL, *matind = NULL, *matcnt = NULL;

    int nzcnt, surplus;

    // Salvo i coefficienti della funzione obiettivo del problema SRC.
    obj = (double*) malloc(numcols * sizeof(double));
    if (obj == NULL) {
        fprintf(stderr, "No memory for saving obj of SRC.\n");
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
        fprintf(stderr, "No memory for saving matbeg or matind or matval or matcnt of SRC.\n");
        status = 1;
        goto TERMINATE;
    }

    status = CPXgetcols(env, src, &nzcnt, matbeg, matind, matval, numnz, &surplus, 0, numcols - 1);
    if (status) {
        fprintf(stderr, "Failed to copy matbeg or matind or matval coefficients of SRC.\n");
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
        fprintf(stderr, "No memory for saving lower bounds or upper bounds of SRC.\n");
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
    // TODO: non copio il vero nome delle variabili. Per ora non è necessario.
    //       Ciò implica che prob originale e copia hanno nomi di var diversi.
    status = CPXcopylp(
            env, dst, numcols, numrows, CPXgetobjsen(env, src), 
            obj, rhs, sense, matbeg, matcnt, matind, matval, lb, ub, NULL
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

    return status;
}


// Permette di creare il problema FMIP partendo dal problema MIP fornito.
// Non effettua il variable fixing.
int create_fmip(CPXENVptr env, CPXLPptr mip, CPXLPptr *fmip) {
    int i, cnt, status = 0;
    double tmp;
    
    int numcols = CPXgetnumcols(env, mip);
    int numrows = CPXgetnumrows(env, mip);

    *fmip = CPXcreateprob(env, &status, "FMIP");
    if (*fmip == NULL) {
        fprintf(stderr, "Failed to create FMIP.\n");
        status = 1;
        goto TERMINATE;
    }

    status = copy_prob(env, mip, *fmip);
    if (status) {
        fprintf(stderr, "Failed populate FMIP with MIP data.\n");
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
    // Gli indici delle variabili slack sono: [numcols, ..., cnt - 1]
    cnt = CPXgetnumcols(env, *fmip);
    for (i = numcols, tmp = 1; i < cnt; i++) {
        status = CPXchgobj(env, *fmip, 1, &i, &tmp);
        if (status) {
            fprintf(stderr, "Failed change variable index %d (slack) coefficient in obj.\n", i);
            goto TERMINATE;
        }
    }

    // Cambio il senso del problema a minimizzazione.
    if (CPXgetobjsen(env, *fmip) != CPX_MIN) {
        status = CPXchgobjsen(env, *fmip, CPX_MIN);
        if (status) {
            fprintf(stderr, "Failed setting FMIP obj sense.\n");
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

    //status = CPXchgbds(env, *fmip, cnt, indices, lu, bd);
    //if (status) {
    //    fprintf(stderr, "Failed to fix variables of FMIP.\n");
    //    goto TERMINATE;
    //}
// ---------------------------------------------------------------------------


TERMINATE:

    free_and_null((char**) &indices);
    free_and_null(&lu);
    free_and_null((char**) &bd);

    return status;
}


// Permette di creare l'OMIP partendo dal problema MIP fornito.
int create_omip(CPXENVptr env, CPXLPptr mip, CPXLPptr *omip, double rhs_slack) {
    int i, cnt, status = 0;
    double tmp;
    
    int numcols = CPXgetnumcols(env, mip);
    int numrows = CPXgetnumrows(env, mip);

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


    // Fisso il valore di alcune variabili xi.
// ---------------------------------------------------------------------------
    // TODO: scrivere l'algoritmo che prende la decisione sulle variabili
    // da fissare. Quello qui sotto è solo un esempio per testare.
// ---------------------------------------------------------------------------

    // Aggiungo il vincolo che limita il grado della non ammissibilità.
    sense = 'L';
    slack_cnt = CPXgetnumcols(env, *omip) - numcols; // Numero di variabili slack aggiunte.

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
    status = CPXaddrows(env, *omip, 0, 1, slack_cnt, &rhs_slack, &sense, &rmatbeg, rmatind, rmatval, NULL, &rowname);
    if (status) {
        fprintf(stderr, "Failed to add row to OMIP.\n");
        goto TERMINATE;
    }

TERMINATE:

    free_and_null((char**) &rmatind);
    free_and_null((char**) &rmatval);

    return status;
}


int main(int argc, char* argv[]) {
    int numcols, solstat; // Solution status.
    double objval; // Objective value.
    double *x_mip = NULL, *x_fmip = NULL, *x_omip = NULL; // Variabiles value.

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
    printf("### Creating MIP ###\n");
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

    // Optimize MIP.
    //numcols = CPXgetnumcols(env, mip); // Get number of variabiles.
    //x_mip = (double*) malloc(numcols * sizeof(double));
    //if (x_mip == NULL) {
    //    fprintf(stderr, "No memory for solution values for MIP.\n");
    //    goto TERMINATE;
    //}

    //status = optimize_prob(env, mip, &objval, &solstat, x_mip, 0, numcols - 1);
    //if (status) {
    //    fprintf(stderr, "Failed to optimize MIP.\n");
    //    goto TERMINATE;
    //}

    // Salvo il problema MIP in un file.
    //status = CPXwriteprob(env, mip, "original_MIP.lp", NULL);
    //if (status) {
    //    fprintf (stderr, "Failed to write MIP to disk.\n");
    //    goto TERMINATE;
    //}

// ### Creazione dell'FMIP ###
    printf("\n### Creating FMIP ###\n");
    status = create_fmip(env, mip, &fmip);
    if (status) {
        fprintf(stderr, "Failed to create FMIP.\n");
        goto TERMINATE;
    }

    // Optimize fmip.
    //numcols = CPXgetnumcols(env, fmip);
    //x_fmip= (double*) malloc(numcols * sizeof(double));
    //if (x_fmip == NULL) {
    //    fprintf(stderr, "No memory for solution values for FMIP.\n");
    //    goto TERMINATE;
    //}

    //status = optimize_prob(env, fmip, &objval, &solstat, x_fmip, 0, numcols - 1);
    //if (status) {
    //    fprintf(stderr, "Failed to optimize FMIP.\n");
    //    goto TERMINATE;
    //}

    // Salvo il problema FMIP in un file.
    status = CPXwriteprob(env, fmip, "fmip.lp", NULL);
    if (status) {
        fprintf (stderr, "Failed to write FMIP to disk.\n");
        goto TERMINATE;
    }

// ### Creazione dell'OMIP ###
    printf("\n### Creazione OMIP ###\n");
    status = create_omip(env, mip, &omip, 200);
    if (status) {
        fprintf(stderr, "Failed to create OMIP.\n");
        goto TERMINATE;
    }

    // Optimize omip.
    //numcols = CPXgetnumcols(env, omip);
    //x_omip= (double*) malloc(numcols * sizeof(double));
    //if (x_omip == NULL) {
    //    fprintf(stderr, "No memory for solution values for OMIP.\n");
    //    goto TERMINATE;
    //}

    //status = optimize_prob(env, omip, &objval, &solstat, x_omip, 0, numcols - 1);
    //if (status) {
    //    fprintf(stderr, "Failed to optimize OMIP.\n");
    //    goto TERMINATE;
    //}

    // Salvo il problema OMIP in un file.
    status = CPXwriteprob(env, omip, "omip.lp", NULL);
    if (status) {
        fprintf (stderr, "Failed to write OMIP to disk.\n");
        goto TERMINATE;
    }

TERMINATE:

    // Free up solution.
    free_and_null((char**) &x_mip);
    free_and_null((char**) &x_fmip);
    free_and_null((char**) &x_omip);

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
