#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>


const double beta = 0.1;
const double alpha = 0.2;
double hittingMatrix[328][328][7];
int katzMatrix[328][328][7];
char jac[] = "Jaccard.txt";
char kat[] = "Katz.txt";
char hit[] = "HittingTime.txt";
char page[] = "PageRank.txt";
char inputPath[] = "contact-high-school-proj-graph.txt";



struct edge {
    int v;
    int w;
    struct edge *next;
};

struct coeff {
    double val;
    int u;
    int v;
};

int getUE(int n, int e){
    return n*(n-1)/2 - e;
}

struct node {
    struct edge *start;
    struct edge *end;
};

void addEdge(int u, int v, int w, struct node * graphAdj[]){
    struct edge *e = malloc(sizeof(struct edge));
    e->v = v;
    e->w = w;
    if(graphAdj[u] == NULL){
        graphAdj[u] = malloc(sizeof(struct node));
        graphAdj[u]->start = e;
        graphAdj[u]->end = e;
    }else{
        graphAdj[u]->end->next = e;
        graphAdj[u]->end = e;
    }
}

void addBiEdge(int u, int v, int w, struct node * graphAdj[]){
    addEdge(u, v, w, graphAdj);
    addEdge(v, u, w, graphAdj);
}

void printEdges(int u, struct node * graphAdj[]){
    struct node *node = graphAdj[u];
    printf("Node %d : ", u);

    if(node != NULL) {
        struct edge * next = node->start;
        while (next != NULL) {
            printf("%d -> ", next->v);
            next = next->next;
        }
    }
    printf("\n");
}

void printGraphAdj(int n, struct node * graphAdj[]){
    printf("graph (adjacent list)\n");
    for(int u=1; u<=n; u++){
        printEdges(u, graphAdj);
    }
}

//jaccards
int commonNeighbors(int u, int v, struct node * graphAdj[]){
    int commonNeighbors = 0;
    struct edge * uNext = graphAdj[u]->start;
    struct edge * vNext = graphAdj[v]->start;
    while (uNext != NULL && vNext != NULL){
        if(uNext->v == vNext->v){
            commonNeighbors++;
            uNext = uNext->next;
            vNext = vNext->next;
        }else if(uNext->v < vNext->v){
            uNext = uNext->next;
        }else{
            vNext = vNext->next;
        }
    }
    return commonNeighbors;
}

int unionOfNeighbors(int u, int v, int n, struct node * graphAdj[]){
    int unionOfNeighbors[n+1];
    for(int i=1; i<=n; i++){
        unionOfNeighbors[i] = 0;
    }
    struct edge * uNext = graphAdj[u]->start;
    struct edge * vNext = graphAdj[v]->start;
    while(uNext != NULL){
        unionOfNeighbors[uNext->v] = 1;
        uNext = uNext->next;
    }
    while(vNext != NULL){
        unionOfNeighbors[vNext->v] = 1;
        vNext = vNext->next;
    }

    int count = 0;
    for(int i=1; i <= n; i++){
        if(unionOfNeighbors[i] == 1){
            count++;
        }
    }
    return count;
}

double jaccardsCoefficient(int u, int v, int n, struct node * graphAdj[]){
    int c = commonNeighbors(u, v, graphAdj);
    if(c == 0){
        return 0.0;
    }
    int b = unionOfNeighbors(u, v, n, graphAdj);
    return 1.0 * c / b;
}

int cmp(const void *x, const void *y)
{
    const struct coeff *xx = x, *yy = y;

    double p = xx->val, q = yy->val;
    if (p < q) return 1;
    if (p > q) return -1;
    return 0;
}

void writeToFile(char fileName[], struct coeff coeffs[], int K, int size){
    qsort(coeffs, size, sizeof(struct coeff), cmp);
    FILE* filePointer = fopen(fileName, "w");
    for(int k=0; k<K; k++){
        fprintf(filePointer, "%f %d %d\n", coeffs[k].val, coeffs[k].u, coeffs[k].v);
    }
    fclose(filePointer);
}

void generateJaccards(int n, int K, struct node * graphAdj[], struct coeff jaccards[]){
    int i = 0;

    for(int u=1; u<=n; u++){
        struct edge * uNext = graphAdj[u]->start;
        for(int v=u+1; v<=n; v++){
            if(uNext != NULL && uNext->v < v){
                uNext = uNext->next;
                v--;
            }else if(uNext != NULL && uNext->v == v){
                uNext = uNext->next;
            }else{
                jaccards[i].val = jaccardsCoefficient(u, v, n, graphAdj);
                jaccards[i].u = u;
                jaccards[i].v = v;
                i++;
            }
        }
    }
    writeToFile(jac, jaccards, K, i);
}

void clear(struct coeff coeffs[], int n){
    for(int i=0; i<n; i++){
        coeffs[i].val = 0.0;
        coeffs[i].v = 0;
        coeffs[i].u = 0;
    }
}

//katz

void clearMatrix(int n, int m[n+1][n+1]){
    for(int i=1; i<=n; i++){
        for(int j=1; j<=n; j++){
            m[i][j] = 0;
        }
    }
}

void multiply(int n, int a[n+1][n+1], int b[n+1][n+1], int c[n+1][n+1]){
    for(int i=1; i<=n; i++){
        for (int j=1; j<=n; j++){
            c[i][j] = 0;
        }
    }

    for(int i=1; i<=n; i++){
        for(int j=1; j<=n; j++){
            for(int k=1; k<=n; k++){
                c[i][j] += a[i][k] * b[k][j];
            }
        }
    }
}

void fillMatrix(int n, struct node * graphAdj[], int a[n+1][n+1]){
    for(int i=1; i<=n; i++){
        if(graphAdj[i] != NULL) {
            struct edge *uNext = graphAdj[i]->start;
            while (uNext != NULL){
                a[i][uNext->v] = 1;
                uNext = uNext->next;
            }
        }
    }
}

void copyMatrix(int n, int a[n+1][n+1], int b[n+1][n+1]){
    for(int i=1; i<=n; i++){
        for(int j=1; j<=n; j++){
            b[i][j] = a[i][j];
        }
    }
}

double getKatzCoefficient(int n, int u, int v){
    double score = 0.0;
    for(int l=2; l<=6; l++){
        score += pow(beta,l) * katzMatrix[u][v][l];
    }
    return score;
}

void fillKatzCoefficient(int n, int l, int mul [n+1][n+1]){
    for(int i=1; i<=n; i++){
        for(int j=1; j<=n; j++){
            katzMatrix[i][j][l] += mul[i][j];
        }
    }
}

void generateKatz(int n, int K, struct node * graphAdj[], struct coeff katz[]){
    int matrix [n+1][n+1];
    clearMatrix(n, matrix);
    fillMatrix(n, graphAdj, matrix);


    for(int i=1; i<=n; i++){
        for(int j=1; j<=n; j++){
            for(int k=2; k<=6; k++){
                katzMatrix[i][j][k] = 0;
            }
        }
    }

    int mul [n+1][n+1];
    int temp [n+1][n+1];

    copyMatrix(n, matrix, temp);
    for(int l=2; l<=6; l++){
        clearMatrix(n, mul);
        multiply(n, matrix, temp, mul);
        fillKatzCoefficient(n, l, mul);
        copyMatrix(n, mul, temp);
    }

    int start = 0;
    for (int i=1; i<=n; i++){
        for(int j=i+1; j<=n; j++){
            if(matrix[i][j] != 1){
                katz[start].val = getKatzCoefficient(n, i, j);
                katz[start].u = i;
                katz[start].v = j;
                start++;
             }
        }
    }

    writeToFile(kat, katz, K, start);
}


//hitting time
void clearMatrixDouble(int n, double m[n+1][n+1]){
    for(int i=1; i<=n; i++){
        for(int j=1; j<=n; j++){
            m[i][j] = 0.0;
        }
    }
}

void multiplyDouble(int n, double a[n+1][n+1], double b[n+1][n+1], double c[n+1][n+1]){
    for(int i=1; i<=n; i++){
        for (int j=1; j<=n; j++){
            c[i][j] = 0;
        }
    }

    for(int i=1; i<=n; i++){
        for(int j=1; j<=n; j++){
            for(int k=1; k<=n; k++){
                c[i][j] += a[i][k] * b[k][j];
            }
        }
    }
}

int getSize(struct node * graphAdj[], int i){
    int size = 0;
    struct edge *uNext = graphAdj[i]->start;
    while (uNext != NULL){
        uNext = uNext->next;
        size++;
    }
    return size;
}
void fillMatrixDouble(int n, struct node * graphAdj[], double a[n+1][n+1]){
    for(int i=1; i<=n; i++){
        if(graphAdj[i] != NULL) {
            int size = getSize(graphAdj, i);
            struct edge *uNext = graphAdj[i]->start;
            while (uNext != NULL){
                a[i][uNext->v] = 1.0/size;
                uNext = uNext->next;
            }
        }
    }
}

void fillMatrixDoubleForPagerank(int n, struct node * graphAdj[], double a[n+1][n+1]){
    for(int i=1; i<=n; i++){
        a[i][i] = alpha;
        if(graphAdj[i] != NULL) {
            int size = getSize(graphAdj, i);
            struct edge *uNext = graphAdj[i]->start;
            while (uNext != NULL){
                a[i][uNext->v] = (1.0-alpha)/size;
                uNext = uNext->next;
            }
        }
    }
}

void copyMatrixDouble(int n, double a[n+1][n+1], double b[n+1][n+1]){
    for(int i=1; i<=n; i++){
        for(int j=1; j<=n; j++){
            b[i][j] = a[i][j];
        }
    }
}

double getHittingCoefficientDouble(int n, int u, int v){
    double score = 0.0;
    for(int l=2; l<=6; l++){
        score += l * hittingMatrix[u][v][l];
    }
    for(int l=2; l<=6; l++){
        score += l * hittingMatrix[v][u][l];
    }
    return score;
}

void fillHittingCoefficient(int n, int l, double mul [n+1][n+1]){
    for(int i=1; i<=n; i++){
        for(int j=1; j<=n; j++){
            hittingMatrix[i][j][l] += mul[i][j];
        }
    }
}

void generateHitting(int n, int K, struct node * graphAdj[], struct coeff hitting[]){
    double matrix [n+1][n+1];
    clearMatrixDouble(n, matrix);
    fillMatrixDouble(n, graphAdj, matrix);
    for(int i=1; i<=n; i++){
        for(int j=1; j<=n; j++){
            for(int k=0; k<=6; k++){
                hittingMatrix[i][j][k] = 0;
            }
        }
    }

    double mul [n+1][n+1];
    double temp [n+1][n+1];

    copyMatrixDouble(n, matrix, temp);
    for(int l=2; l<=6; l++){
        clearMatrixDouble(n, mul);
        multiplyDouble(n, matrix, temp, mul);
        fillHittingCoefficient(n, l, mul);
        copyMatrixDouble(n, mul, temp);
    }

    int start = 0;
    for (int i=1; i<=n; i++){
        for(int j=i+1; j<=n; j++){
            if(matrix[i][j] == 0.0){
                double score = -getHittingCoefficientDouble(n, i, j);
                hitting[start].val = score;
                hitting[start].u = i;
                hitting[start].v = j;
                start++;
            }
        }
    }

    writeToFile(hit, hitting, K, start);
}

void generatePagerank(int n, int K, struct node * graphAdj[], struct coeff pagerank[]){
    double matrix [n+1][n+1];
    clearMatrixDouble(n, matrix);
    fillMatrixDoubleForPagerank(n, graphAdj, matrix);
    for(int i=1; i<=n; i++){
        for(int j=1; j<=n; j++){
            for(int k=0; k<=6; k++){
                hittingMatrix[i][j][k] = 0;
            }
        }
    }

    double mul [n+1][n+1];
    double temp [n+1][n+1];

    copyMatrixDouble(n, matrix, temp);
    for(int l=2; l<=6; l++){
        clearMatrixDouble(n, mul);
        multiplyDouble(n, matrix, temp, mul);
        fillHittingCoefficient(n, l, mul);
        copyMatrixDouble(n, mul, temp);
    }

    int start = 0;
    for (int i=1; i<=n; i++){
        for(int j=i+1; j<=n; j++){
            if(matrix[i][j] == 0.0){
                double score = -getHittingCoefficientDouble(n, i, j);
                pagerank[start].val = score;
                pagerank[start].u = i;
                pagerank[start].v = j;
                start++;
            }
        }
    }

    writeToFile(page, pagerank, K, start);
}

int main(){
    FILE* filePointer;
    char buffer[255];
    printf("Please make sure - file 'contact-high-school-proj-graph.txt' is present in current directory\n");
    int N = 0, K = 0;
    printf("Enter no. of vertex : ");
    scanf("%d", &N);
    printf("Enter K : ");
    scanf("%d", &K);
    printf("\n");
    struct node * graphAdj[N+1];
    for(int i=1; i<=N; i++){
        graphAdj[i] = NULL;
    }



    filePointer = fopen(inputPath, "r+");


    int E = 0;
    while(fgets(buffer, 255, filePointer)) {
        int u, v, w;
        sscanf(buffer, "%d %d %d", &u, &v, &w);
        if(u > N || u < 1 || v > N || v < 1){
            printf("vertex range is (1, %d). invalid edge found %d %d %d\n", N, u, v, w);
            return 1;
        }
        addBiEdge(u, v, w, graphAdj);
        E++;
    }

    int ue = getUE(N, E);
    struct coeff jaccards[ue];
    struct coeff katz[ue] ;
    struct coeff hitting[ue] ;
    struct coeff pagerank[ue] ;
    clear(jaccards, ue);
    clear(katz, ue);
    clear(hitting, ue);
    clear(pagerank, ue);


    fclose(filePointer);


    //printGraphAdj(N, graphAdj);
    printf("Generating Jaccard.txt\n");
    generateJaccards(N, K, graphAdj, jaccards);
    printf("Generating Katz.txt\n");
    generateKatz(N, K, graphAdj, katz);
    printf("Generating HittingTime.txt\n");
    generateHitting(N, K, graphAdj, hitting);
    printf("Generating PageRank.txt\n");
    generatePagerank(N,K,graphAdj, pagerank);
}