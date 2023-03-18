#define _CRT_SECURE_NO_WARNINGS
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wstrict-prototypes"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <Accelerate/Accelerate.h>
#include <malloc/malloc.h>
#include <f2c.h>
#include "matx.h" /* マトリックスのかけ算 */


#define FL 500 /* t@C? */
#define DM 2  /*  */
#define NP 7000 /* ?_?? */
#define NM 2  /* ?? */
#define NE 3000  /* vf? */
#define NB 200 /* S?_? */
#define NN DM*NP
#define E8 8 /* 1vf??_ */
#define E16 16 /* EK}gbNX?p16 */
#define M3 3 /* D}gbNX?p3  */
#define GP 3 /* KEX|Cg?  */

#define NVF 50 /* VF? *///
#define SHV 0.001
#define DE 1.0
#define DNU 0.01
#define DMODULUS (0.0001)//0.005
#define DN 0.01

#define SHIKII (0.01)

#define PI (3.14159265359)
#define M (5) /* ??????gp */
#define N (3) /* ??????gp */
#define C1 (8.86)
#define C2 (101.6)
#define TREF (303.0)
#define NUM_OF_TERM (20)
#define NUM_OF_DATA (2000)//各ひずみ成分と板厚方向ひずみの変数サイズの設定に使われてる
#define NUM_OF_MODEL (50)//20
#define POISSON (0.468)
#define TEMP (273.0)//温度変更項
#define NUM_OF_UNKNOWN (50)

int num_of_data;
int num_of_model;
int num_of_vf;
int num_of_term;
int np, nm, ne, nb, nn;
int nip, nnip;//np : number of integration point
double temp;

int readmodulus(char *file, double *time, double *gcoef, double *kcoef, double *ge, double *ke)
{
    FILE *fp;
    int i;
    i = 0;
    if (NULL == (fp = fopen(file, "r"))) {
        printf("\7\n Cannot open modulusfile..\n\n");
        exit(1);
    }
    fscanf(fp, "%lf", ge);
    fscanf(fp, "%lf", ke);
    while (3 == fscanf(fp, "%lf %lf %lf", time, gcoef, kcoef)) {
        //printf("%lf\t%lf\t%lf\n", *time, *gcoef, *kcoef);
        time++;
        gcoef++;
        kcoef++;
        i++;
    }
    fclose(fp);
    //num_of_term = i;
    return i;
}

void readmodel(char *filen, double *pt, double *pf, double *mat, int *el, int *ib, int *ibn, double *disp)
{
    int i;
    int j;
    FILE *fp;
    if (NULL == (fp = fopen(filen, "r"))) {
        printf(" Cannot open model file\n");
        exit(-1);
    }
    fscanf(fp, "%d %d %d %d", &np, &nm, &ne, &nb);
    nn = np*DM;
    nip = ne*(E8 + 1);
    nnip = 3 * nip;
    
    for (i = 0; i < np; i++, pt += 2, pf += 2){
        fscanf(fp, "%lf %lf %lf %lf", pt, pt + 1, pf, pf + 1);
        //printf("%lf %lf %lf %lf\n", *pt, *(pt + 1), *pf, *(pf + 1));
    }
    for (i = 0; i < nm; i++, mat += 3){
        fscanf(fp, "%lf %lf %lf", mat, mat + 1, mat + 2);
    }
    for (i = 0; i < ne; i++) {
        for (j = 0; j < E8; j++) {
            fscanf(fp, "%d", el++);
        }
        fscanf(fp, "%d", el++);
    }
    for (i = 0; i < nb; i++, ibn += 1, ib += 2, disp += 2) {
        fscanf(fp, "%d %d %d %lf %lf", ibn, ib, ib + 1, disp, disp + 1);
    }
    fclose(fp);
}

void read_displacement(char *filedata, double *uuord)
{
    int i;
    FILE *fp;
    double a, b;
    
    i = 0;
    if (NULL == (fp = fopen(filedata, "r"))) {
        printf(" Cannot open displacement file\n");
        exit(1);
    }
    for (i = 0; i < np; i++, uuord += 2) {
        fscanf(fp, "%lf %lf %lf %lf", &a, &b, uuord, uuord + 1);
    }
    fclose(fp);
}

void make_gp(double *g_po, double *g_we)
{
    *g_po++ = -0.7746;
    *g_po++ = 0.0;
    *g_po = 0.7746;
    *g_we++ = 0.5556;
    *g_we++ = 0.8889;
    *g_we = 0.5556;
}

void make_pt(double *elpt, int *ell, double *pt)
{
    int l;
    for(l = 0; l < E8; l++, ell++){
        *elpt++ = *(pt + *ell*DM);
        *elpt++ = *(pt + *ell*DM + 1);
    }
}

void make_sp(double *shape, double *deriv, double s, double t)
{
    double s2 = s*2;
    double t2 = t*2;
    double ss = s*s;
    double tt = t*t;
    double st = s*t;
    double sst = s*st;
    double stt = st*t;
    double st2 = st*2;

    *shape++ = (-1 + st + ss + tt - sst - stt) / 4;
    *shape++ = (1 - t - ss + sst) / 2;
    *shape++ = (-1 - st + ss + tt - sst + stt) / 4;
    *shape++ = (1 + s - tt - stt) / 2;
    *shape++ = (-1 + st + ss + tt + sst + stt) / 4;
    *shape++ = (1 + t - ss - sst) / 2;
    *shape++ = (-1 - st + ss + tt + sst - stt) / 4;
    *shape = (1 - s - tt + stt) / 2;
    
    *deriv++ = (t + s2 - st2 - tt) / 4;
    *deriv++ = -s + st;
    *deriv++ = (-t + s2 - st2 + tt) / 4;
    *deriv++ = (1 - tt) / 2;
    *deriv++ = (t + s2 + st2 + tt) / 4;
    *deriv++ = -s - st;
    *deriv++ = (-t + s2 + st2 - tt) / 4;
    *deriv++ = (-1 + tt) / 2;
    *deriv++ = (s + t2 - ss - st2) / 4;
    *deriv++ = (-1 + ss) / 2;
    *deriv++ = (-s + t2 - ss + st2) / 4;
    *deriv++ = -t - st;
    *deriv++ = (s + t2 + ss + st2) / 4;
    *deriv++ = (1 - ss) / 2;
    *deriv++ = (-s + t2 + ss - st2) / 4;
    *deriv = -t + st;
}

void jacob(double *cart, double *djacb, double *elpt, double *deriv)
{
    double *jacm;
    double *jaci;
    
    jacm=(double *)malloc(sizeof(double)*DM*DM);
    jaci=(double *)malloc(sizeof(double)*DM*DM);
    
    matx(jacm, deriv, elpt, DM, E8, DM);
    *djacb = (*(jacm + 0))*(*(jacm + 3)) - (*(jacm + 1))*(*(jacm + 2));
    if(*djacb < 0.0){
        printf("Determinant < zero\n");
        printf("%f\t%f\n", (*(elpt + 0)), (*(elpt + 1)));
        printf("%f\t%f\n", (*(elpt + 2)), (*(elpt + 3)));
        printf("%f\t%f\n", (*(elpt + 4)), (*(elpt + 5)));
        printf("%f\t%f\n", (*(elpt + 6)), (*(elpt + 7)));
        printf("%f\t%f\n", (*(elpt + 8)), (*(elpt + 9)));
        printf("%f\t%f\n", (*(elpt + 10)), (*(elpt + 11)));
        printf("%f\t%f\n", (*(elpt + 12)), (*(elpt + 13)));
        printf("%f\t%f\n", (*(elpt + 14)), (*(elpt + 15)));
        exit(1);
    }
    (*(jaci + 0)) = (*(jacm + 3)) / (*djacb);
    (*(jaci + 3)) = (*(jacm + 0)) / (*djacb);
    (*(jaci + 1)) = -(*(jacm + 1)) / (*djacb);
    (*(jaci + 2)) = -(*(jacm + 2)) / (*djacb);
    matx(cart, jaci, deriv, DM, DM, E8);
    
    free(jacm);
    free(jaci);
}

void make_b(double *b, double *cart)
{
    double *b1,*b2,*cart1;
    int i;

    for (i = 0, b1 = b + E16, b2 = b1 + E16, cart1 = cart + E8; i < E8; i++) {
        *b++ = *cart;
        *b++ = 0.0;
        *b1++ = 0.0;
        *b1++ = *cart1;
        *b2++ = *cart1++;
        *b2++ = *cart++;
    }
}

void strain(double *uuord, int *el, double *pt, double *points, double *epspoints)
{
    double *elpt, *disp, *bmat, *eps;
    double *shape, *deriv, *cart, *g_po, *g_we;
    double *stcoord;
    double px, py, det, xcoord, ycoord, dx, dy;
    int i, l, k, j, ig, jg;
    g_po = (double *)malloc(sizeof(double)*GP);
    g_we = (double *)malloc(sizeof(double)*GP);
    cart = (double *)malloc(sizeof(double)*DM*E8);
    shape = (double *)malloc(sizeof(double)*E8);
    deriv = (double *)malloc(sizeof(double)*DM*E8);
    stcoord = (double *)malloc(sizeof(double)*DM*E8);
    elpt = (double *)malloc(sizeof(double)*DM*E8);
    disp = (double *)malloc(sizeof(double)*DM*E8);
    bmat = (double *)malloc(sizeof(double) * 16 * 3);
    eps = (double *)malloc(sizeof(double) * 3 * 10);

    make_gp(g_po, g_we);
    for (j = 0, l = 0; l < ne; l++, el += 9) {
        make_pt(elpt, el, pt);
        for (i = 0; i < E8; i++) {
            k = *(el + i);
            *(disp + 2 * i) = (*(uuord + 2 * k));
            *(disp + 2 * i + 1) = (*(uuord + 2 * k + 1));
        }
        for (ig = 0; ig < GP; ig++) {
            for (jg = 0; jg < GP; jg++) {
                px = (*(g_po + ig));
                py = (*(g_po + jg));
                make_sp(shape, deriv, px, py);
                jacob(cart, &det, elpt, deriv);
                make_b(bmat, cart);
                matx(eps, bmat, disp, 3, 16, 1);
                xcoord = (*(shape))*(*(elpt)) + (*(shape + 1))*(*(elpt + 2)) + (*(shape + 2))*(*(elpt + 4)) + (*(shape + 3))*(*(elpt + 6)) + (*(shape + 4))*(*(elpt + 8)) + (*(shape + 5))*(*(elpt + 10)) + (*(shape + 6))*(*(elpt + 12)) + (*(shape + 7))*(*(elpt + 14));
                ycoord = (*(shape))*(*(elpt + 1)) + (*(shape + 1))*(*(elpt + 3)) + (*(shape + 2))*(*(elpt + 5)) + (*(shape + 3))*(*(elpt + 7)) + (*(shape + 4))*(*(elpt + 9)) + (*(shape + 5))*(*(elpt + 11)) + (*(shape + 6))*(*(elpt + 13)) + (*(shape + 7))*(*(elpt + 15));
                dx = (*(shape))*(*(disp)) + (*(shape + 1))*(*(disp + 2)) + (*(shape + 2))*(*(disp + 4)) + (*(shape + 3))*(*(disp + 6)) + (*(shape + 4))*(*(disp + 8)) + (*(shape + 5))*(*(disp + 10)) + (*(shape + 6))*(*(disp + 12)) + (*(shape + 7))*(*(disp + 14));
                dy = (*(shape))*(*(disp + 1)) + (*(shape + 1))*(*(disp + 3)) + (*(shape + 2))*(*(disp + 5)) + (*(shape + 3))*(*(disp + 7)) + (*(shape + 4))*(*(disp + 9)) + (*(shape + 5))*(*(disp + 11)) + (*(shape + 6))*(*(disp + 13)) + (*(shape + 7))*(*(disp + 15));
                (*(points + 2 * j)) = xcoord;
                (*(points + 2 * j + 1)) = ycoord;
                (*(epspoints + 3 * j)) = (*(eps + 0));
                (*(epspoints + 3 * j + 1)) = (*(eps + 1));
                (*(epspoints + 3 * j + 2)) = (*(eps + 2));
                j++;
            }
        }
    }
    free(g_po);
    free(g_we);
    free(cart);
    free(shape);
    free(deriv);
    free(stcoord);
    free(elpt);
    free(disp);
    free(bmat);
    free(eps);
}

void differentiatelapack(double *x, double *y, double *dydx)
{
    int i;
    double a, b, c;
    double *mat_a, *mat_f;
    int nrhs = 1, lwork, info, lda, ldb, m, n;
    double *work;
    char trans = 'N';
    mat_a = (double *)malloc(sizeof(double)*M*N);
    mat_f = (double *)malloc(sizeof(double)*M);
    work = (double *)malloc(sizeof(double)*(M + N));
    
    for (i = 2; i < num_of_data - 2; i++) {
            (*(mat_a + 0)) = (*(x + i - 2))*(*(x + i - 2));
            (*(mat_a + 1)) = (*(x + i - 1))*(*(x + i - 1));
            (*(mat_a + 2)) = (*(x + i))*(*(x + i));
            (*(mat_a + 3)) = (*(x + i + 1))*(*(x + i + 1));
            (*(mat_a + 4)) = (*(x + i + 2))*(*(x + i + 2));
            (*(mat_a + 5)) = (*(x + i - 2));
            (*(mat_a + 6)) = (*(x + i - 1));
            (*(mat_a + 7)) = (*(x + i));
            (*(mat_a + 8)) = (*(x + i + 1));
            (*(mat_a + 9)) = (*(x + i + 2));
            (*(mat_a + 10)) = 1;
            (*(mat_a + 11)) = 1;
            (*(mat_a + 12)) = 1;
            (*(mat_a + 13)) = 1;
            (*(mat_a + 14)) = 1;
            (*(mat_f + 0)) = (*(y + i - 2));
            (*(mat_f + 1)) = (*(y + i - 1));
            (*(mat_f + 2)) = (*(y + i));
            (*(mat_f + 3)) = (*(y + i + 1));
            (*(mat_f + 4)) = (*(y + i + 2));
            n = N;
            m = M;
            lwork = m + n;
            lda = m;
            ldb = m;
            dgels_(&trans, &m, &n, &nrhs, mat_a, &lda, mat_f, &ldb, work, &lwork, &info);
            a = (*(mat_f));
            b = (*(mat_f + 1));
            c = (*(mat_f + 2));
            (*(dydx + i)) = (2 * a*(*(x + i)) + b);
            if (i == 2) {
                (*(dydx + 0)) = (2 * a*(*(x + 0)) + b);
                (*(dydx + 1)) = (2 * a*(*(x + 1)) + b);
            }
            if (i == num_of_data - 3) {
                (*(dydx + num_of_data - 2)) = (2 * a*(*(x + num_of_data - 2)) + b);
                (*(dydx + num_of_data - 1)) = (2 * a*(*(x + num_of_data - 1)) + b);
            }
        }
    free(mat_a);
    free(mat_f);
    free(work);
}

double epssearchkai(double *x, double *y, double t, int num_of_data, double t0, double t1, double t2, int i)
{
    double eps = 0;
 
    if (t1 == *(x + 0) && t2 == *(x + 1)) {
        eps = (*(y + 1) - *(y + 0)) / (*(x + 1) - *(x + 0))*(t - *(x + 0)) + *(y + 0);
    }
    else {
        eps = *(y + i - 1)*((t - t1))*(t - t2)/((t0 - t1)*(t0 - t2)) +
            *(y + i)*(t - t0)*(t - t2)/((t1 - t0)*(t1 - t2)) +
            *(y + i + 1)*(t - t0)*(t - t1)/((t2 - t0)*(t2 - t1));
    }
    return eps;
}

double epssearch(double *x, double *y, double t, int num_of_data)
{
    int i;
    double eps = 0;
    double x1 = 0, x2 = 0, y1 = 0, y2 = 0;
    
    if (t == (*(x + 0)))
        eps = (*(y + 0));
    for (i = 1; i < num_of_data; i++) {
        if (t == (*(x + i))) {
            eps = (*(y + i));
            break;
        }
        else {
            if ((t >(*(x + i - 1)) && t < (*(x + i)))) {
                x1 = (*(x + i - 1));
                x2 = (*(x + i));
                y1 = (*(y + i - 1));
                y2 = (*(y + i));
                eps = (y2 - y1) / (x2 - x1)*(t - x1) + y1;
                break;
            }
        }
    }
    return eps;
}

double relaxationmodulus(double *tau, double *modulus, double constant, int num_of_term, double t, double at)
{
    int i;
    double et;

    for (et = constant, i = 0; i < num_of_term; i++) {
        et += (*(modulus + i))*exp(-t / ((*(tau + i))*at));
    }
    return et;
}

double integrate3(double *x, double *y, double *dydx, int num_of_data, double *tau, double *modulus, double constant, int num_of_term, double t, double temp)/* 台形公式による積分の関数 材料システム2.90式 */
{
    int i, n;
    double s, et1, et2;
    double t0, t1, t2, t3, t4, deltat;
    double at;
    double eps1, eps2;
    
    at = pow(10.0, -C1*(temp - TREF) / (C2 + (temp - TREF)));
    deltat = *(x + num_of_data - 1)/(num_of_data-1);//角周波数間隔
    for (s = 0, i = 0; i < num_of_data - 1; i++) {
        t1 = (*(x + i));
        t2 = (*(x + i + 1));
        if (i == 0 && t2 > t) {
            break;
        }
        if (t2 == t) {
            for (n = 0; n <100; n++) {
                t0 = *(x + i - 1);
                t3 = *(x + i) + deltat*n /100;
                t4 = *(x + i) + deltat*(n + 1) / 100;
                eps1 = epssearchkai(x, dydx, t3, num_of_data, t0, t1, t2, i);
                eps2 = epssearchkai(x, dydx, t4, num_of_data, t0, t1, t2, i);
                et1 = relaxationmodulus(tau, modulus, constant, num_of_term, t - t3, at);
                et2 = relaxationmodulus(tau, modulus, constant, num_of_term, t - t4, at);
                s += (et1*eps1 + et2*eps2)*(t4 - t3) / 2;
                if (t4 == t2) {
                    break;
                }
            }
        }
        else {
            eps1 = epssearch(x, dydx, t1, num_of_data);
            eps2 = epssearch(x, dydx, t2, num_of_data);
            et1 = relaxationmodulus(tau, modulus, constant, num_of_term, t - t1, at);
            et2 = relaxationmodulus(tau, modulus, constant, num_of_term, t - t2, at);
            s += (et1*eps1 + et2*eps2)*(t2 - t1) / 2;
        }
        if (t2 >= t)
            break;
    }
    return s;
}

double through_thickness_strain(double *time, double *timez,double * tau, double *gmodulus, double *kmodulus, double *ex, double *ey, double *ez, double ge, double ke, int *number, int num_of_term)
{
    double *Re11, *Im11, *Re22, *Im22, *omega, coefgamma = 5.0, gamma, seikou = 4545, *AnsRe33, *AnsIm33, *ResGs, *ImsGs, *ResKs, *ImsKs;
    int  n = 0;
    int k = 0;
    int i = 0;
    double at;
    double *Reet, *Imet, *Reez, *Imez;
    double  pr1, pr2, c, cc, s, ss;
    Re11 = (double *)malloc(sizeof(double)*NUM_OF_DATA);
    Im11 = (double *)malloc(sizeof(double)*NUM_OF_DATA);
    Re22 = (double *)malloc(sizeof(double)*NUM_OF_DATA);
    Im22 = (double *)malloc(sizeof(double)*NUM_OF_DATA);
    AnsRe33 = (double *)malloc(sizeof(double)*NUM_OF_DATA);
    AnsIm33 = (double *)malloc(sizeof(double)*NUM_OF_DATA);
    omega = (double *)malloc(sizeof(double)*NUM_OF_DATA);
    Reet = (double *)malloc(sizeof(double)*NUM_OF_DATA);
    Imet = (double *)malloc(sizeof(double)*NUM_OF_DATA);
    ResGs = (double *)malloc(sizeof(double)*NUM_OF_DATA);
    ImsGs = (double *)malloc(sizeof(double)*NUM_OF_DATA);
    ResKs = (double *)malloc(sizeof(double)*NUM_OF_DATA);
    ImsKs = (double *)malloc(sizeof(double)*NUM_OF_DATA);
    Reez = (double *)malloc(sizeof(double)*NUM_OF_DATA);
    Imez = (double *)malloc(sizeof(double)*NUM_OF_DATA);
    /*double deltat = 0.1;*/
    double deltat = (*(time + 1)) - (*(time + 0));
    double deltaomega = 2.0* PI / (*(number + 0)) / deltat;
    
    gamma = coefgamma / ((*(time + (*(number + 0)) - 1)) - deltat);
    for (n = 0; n < (*(number + 0)) / 2 + 1; n++) {//ε11をFFT
        (*(Re11 + n)) = 0;
        (*(Im11 + n)) = 0;
        for (k = 0; k < (*(number + 0)); k++) {
            (*(Re11 + n)) += (*(ex + k)) * exp(-1.0 * gamma * k*deltat)*deltat*cos(2.0 * PI*k*n / (*(number + 0)));//実部を計算
            (*(Im11 + n)) += -(*(ex + k)) * exp(-1.0 * gamma * k*deltat)*deltat*sin(2.0 * PI*k*n / (*(number + 0)));//虚部を計算
        }
        (*(omega + n)) = deltaomega*n;
    }
    n = 0;
    k = 0;
    for (n = 0; n < (*(number + 0)) / 2 + 1; n++) {//ε22をFFT
        (*(Re22 + n)) = 0;
        (*(Im22 + n)) = 0;
        for (k = 0; k < (*(number + 0)); k++) {
            (*(Re22 + n)) += (*(ey + k)) * exp(-1.0 * gamma * k*deltat)*deltat*cos(2.0 * PI*k*n / (*(number + 0)));//実部を計算
            (*(Im22 + n)) += -(*(ey + k)) * exp(-1.0 * gamma * k*deltat)*deltat*sin(2.0 * PI*k*n / (*(number + 0)));//虚部を計算
        }
    }
    *(Im11 + 0) = 0;
    *(Im22 + 0) = 0;
    at = pow(10.0, -C1*(temp - TREF) / (C2 + (temp - TREF)));
    for (n = 0; n < (*(number + 0)) / 2 + 1; n++) {//sGsの計算
        for ((*(ResGs + n)) = ge, i = 0; i < num_of_term; i++) {
            (*(ResGs + n)) += (gamma*gamma*(*(tau + i))*(*(tau + i)) * at * at*(*(gmodulus + i)) + gamma*(*(tau + i)) * at*(*(gmodulus + i)) + (*(omega + n))* (*(omega + n))*(*(tau + i))*(*(tau + i)) * at * at*(*(gmodulus + i))) / (((*(tau + i)) * at*gamma + 1.00000000)* ((*(tau + i)) * at*gamma + 1.0000000) + ((*(omega + n))*(*(tau + i)) * at)*((*(omega + n))*(*(tau + i)) * at));
        }
        for ((*(ImsGs + n)) = 0, i = 0; i < num_of_term; i++) {
            (*(ImsGs + n)) += ((*(omega + n))*(*(gmodulus + i))*(*(tau + i)) * at) / (((*(tau + i)) * at*gamma + 1.0)* ((*(tau + i)) * at*gamma + 1.0) + ((*(omega + n))*(*(tau + i)) * at)*((*(omega + n))*(*(tau + i)) * at));
        }
    }
    for (n = 0; n < (*(number + 0)) / 2 + 1; n++) {//sKsの計算
        for ((*(ResKs + n)) = ke, i = 0; i < num_of_term; i++) {
            (*(ResKs + n)) += (gamma*gamma*(*(tau + i))*(*(tau + i)) * at * at*(*(kmodulus + i)) + gamma*(*(tau + i)) * at*(*(kmodulus + i)) + (*(omega + n))* (*(omega + n))*(*(tau + i))*(*(tau + i)) * at * at*(*(kmodulus + i))) / (((*(tau + i)) * at*gamma + 1.00000000)* ((*(tau + i)) * at*gamma + 1.0000000) + ((*(omega + n))*(*(tau + i)) * at)*((*(omega + n))*(*(tau + i)) * at));
        }
        for ((*(ImsKs + n)) = 0, i = 0; i < num_of_term; i++) {
            (*(ImsKs + n)) += ((*(omega + n))*(*(kmodulus + i))*(*(tau + i)) * at) / (((*(tau + i)) * at*gamma + 1.0)* ((*(tau + i)) * at*gamma + 1.0) + ((*(omega + n))*(*(tau + i)) * at)*((*(omega + n))*(*(tau + i)) * at));
        }
    }
    for (n = 0, (*(Reez + n)) = 0, (*(Imez + n)) = 0; n < (*(number + 0)) / 2 + 1; n++) {//εz(s)の計算(ポアソン比一定の場合)
        (*(Reet + n)) = POISSON;//一定
        (*(Imet + n)) = 0;//一定
        (*(Reez + n)) = -(((*(Re11 + n)) + (*(Re22 + n)))*((*(Reet + n)) - (*(Reet + n))*(*(Reet + n)) - (*(Imet + n))*(*(Imet + n))) - (*(Imet + n))*((*(Im11 + n)) + (*(Im22 + n)))) / ((1.00000000 - (*(Reet + n)))*(1.000000000 - (*(Reet + n))) + (*(Imet + n))*(*(Imet + n)));
        (*(Imez + n)) = -(((*(Im11 + n)) + (*(Im22 + n)))*((*(Reet + n)) - (*(Reet + n))*(*(Reet + n)) - (*(Imet + n))*(*(Imet + n))) + (*(Imet + n))*((*(Re11 + n)) + (*(Re22 + n)))) / ((1.0000000000 - (*(Reet + n)))*(1.000000000 - (*(Reet + n))) + (*(Imet + n))*(*(Imet + n)));
    }
    for (n = 0, (*(Reez + n)) = 0, (*(Imez + n)) = 0; n < (*(number + 0)) / 2 + 1; n++) {//εz(s)の計算(ポアソン比一定の場合)
        (*(Reet + n)) = POISSON;;//一定
        (*(Imet + n)) = 0;//一定
        (*(Reez + n)) = -(((*(Re11 + n)) + (*(Re22 + n)))*((*(Reet + n)) - (*(Reet + n))*(*(Reet + n)) - (*(Imet + n))*(*(Imet + n))) - (*(Imet + n))*((*(Im11 + n)) + (*(Im22 + n)))) / ((1.00000000 - (*(Reet + n)))*(1.000000000 - (*(Reet + n))) + (*(Imet + n))*(*(Imet + n)));
        (*(Imez + n)) = -(((*(Im11 + n)) + (*(Im22 + n)))*((*(Reet + n)) - (*(Reet + n))*(*(Reet + n)) - (*(Imet + n))*(*(Imet + n))) + (*(Imet + n))*((*(Re11 + n)) + (*(Re22 + n)))) / ((1.0000000000 - (*(Reet + n)))*(1.000000000 - (*(Reet + n))) + (*(Imet + n))*(*(Imet + n)));
    }
    for (n = 0, (*(Reez + n)) = 0, (*(Imez + n)) = 0; n < (*(number + 0)) / 2 + 1; n++) {//εz(s)の計算
        (*(Reez + n)) = -((3 * (*(ResKs + n)) + 4 * (*(ResGs + n)))*(3 * ((*(ResKs + n))*((*(Re11 + n)) + (*(Re22 + n))) - (*(ImsKs + n))*((*(Im11 + n)) + (*(Im22 + n)))) - 2 * ((*(ResGs + n))*((*(Re11 + n)) + (*(Re22 + n))) - (*(ImsGs + n))*((*(Im11 + n)) + (*(Im22 + n))))) + (3 * (*(ImsKs + n)) + 4 * (*(ImsGs + n)))*(3 * ((*(ResKs + n))*((*(Im11 + n)) + (*(Im22 + n))) + (*(ImsKs + n))*((*(Re11 + n)) + (*(Re22 + n)))) - 2 * ((*(ResGs + n))*((*(Im11 + n)) + (*(Im22 + n))) + (*(ImsGs + n))*((*(Re11 + n)) + (*(Re22 + n)))))) / ((3 * (*(ResKs + n)) + 4 * (*(ResGs + n)))*(3 * (*(ResKs + n)) + 4 * (*(ResGs + n))) + (3 * (*(ImsKs + n)) + 4 * (*(ImsGs + n)))*(3 * (*(ImsKs + n)) + 4 * (*(ImsGs + n))));
        (*(Imez + n)) = -((3 * (*(ResKs + n)) + 4 * (*(ResGs + n)))*(3 * ((*(ResKs + n))*((*(Im11 + n)) + (*(Im22 + n))) + (*(ImsKs + n))*((*(Re11 + n)) + (*(Re22 + n)))) - 2 * ((*(ResGs + n))*((*(Im11 + n)) + (*(Im22 + n))) + (*(ImsGs + n))*((*(Re11 + n)) + (*(Re22 + n))))) - (3 * (*(ImsKs + n)) + 4 * (*(ImsGs + n)))*(3 * ((*(ResKs + n))*((*(Re11 + n)) + (*(Re22 + n))) - (*(ImsKs + n))*((*(Im11 + n)) + (*(Im22 + n)))) - 2 * ((*(ResGs + n))*((*(Re11 + n)) + (*(Re22 + n))) - (*(ImsGs + n))*((*(Im11 + n)) + (*(Im22 + n)))))) / ((3 * (*(ResKs + n)) + 4 * (*(ResGs + n)))*(3 * (*(ResKs + n)) + 4 * (*(ResGs + n))) + (3 * (*(ImsKs + n)) + 4 * (*(ImsGs + n)))*(3 * (*(ImsKs + n)) + 4 * (*(ImsGs + n))));
    }
    for (n = 0, (*(Reez + n)) = 0, (*(Imez + n)) = 0; n < (*(number + 0)) / 2 + 1; n++) {//εz(s)の計算
        (*(Reez + n)) = -(((((3 * (*(ResKs + n)) - 2 * (*(ResGs + n)))*(3 * (*(ResKs + n)) + 4 * (*(ResGs + n)))) + ((3 * (*(ImsKs + n)) - 2 * (*(ImsGs + n)))*(3 * (*(ImsKs + n)) + 4 * (*(ImsGs + n)))))*((*(Re11 + n)) + (*(Re22 + n)))) - ((((3 * (*(ImsKs + n)) - 2 * (*(ImsGs + n)))*(3 * (*(ResKs + n)) + 4 * (*(ResGs + n)))) - ((3 * (*(ResKs + n)) - 2 * (*(ResGs + n)))*(3 * (*(ImsKs + n)) + 4 * (*(ImsGs + n)))))*((*(Im11 + n)) + (*(Im22 + n))))) / ((3 * (*(ResKs + n)) + 4 * (*(ResGs + n)))*(3 * (*(ResKs + n)) + 4 * (*(ResGs + n))) + (3 * (*(ImsKs + n)) + 4 * (*(ImsGs + n)))*(3 * (*(ImsKs + n)) + 4 * (*(ImsGs + n))));
        (*(Imez + n)) = -(((((3 * (*(ResKs + n)) - 2 * (*(ResGs + n)))*(3 * (*(ResKs + n)) + 4 * (*(ResGs + n)))) + ((3 * (*(ImsKs + n)) - 2 * (*(ImsGs + n)))*(3 * (*(ImsKs + n)) + 4 * (*(ImsGs + n)))))*((*(Im11 + n)) + (*(Im22 + n)))) + ((((3 * (*(ImsKs + n)) - 2 * (*(ImsGs + n)))*(3 * (*(ResKs + n)) + 4 * (*(ResGs + n)))) - ((3 * (*(ResKs + n)) - 2 * (*(ResGs + n)))*(3 * (*(ImsKs + n)) + 4 * (*(ImsGs + n)))))*((*(Re11 + n)) + (*(Re22 + n))))) / ((3 * (*(ResKs + n)) + 4 * (*(ResGs + n)))*(3 * (*(ResKs + n)) + 4 * (*(ResGs + n))) + (3 * (*(ImsKs + n)) + 4 * (*(ImsGs + n)))*(3 * (*(ImsKs + n)) + 4 * (*(ImsGs + n))));
    }
    for (n = 0; n < (*(number + 0)) / 2 - 1; n++) {
        (*(Reez + n + (*(number + 0)) / 2 + 1)) = (*(Reez + (*(number + 0)) / 2 - n - 1));
        (*(Imez + n + (*(number + 0)) / 2 + 1)) = (*(Imez + (*(number + 0)) / 2 - n - 1))*(-1);
    }
    deltat = 2.0* PI / (*(number + 0)) / deltaomega;//時間間隔
    for (n = 0; n < (*(number + 0)); n++) {//33成分を逆変換
        (*(AnsRe33 + n)) = 0;
        (*(AnsIm33 + n)) = 0;
        for (k = 0; k < (*(number + 0)); k++) {
            c = cos((2.0 * PI*k*n) / (*(number + 0)));
            cc = fabs(c);
            if (cc <= 0.000001) {
                c = 0;
            }
            s = sin((2.0 * PI*k*n) / (*(number + 0)));
            ss = fabs(s);
            if (ss <= 0.000001) {
                s = 0;
            }
            (*(AnsRe33 + n)) += ((*(Reez + k)) * c - (*(Imez + k))* s) * exp(gamma * n*deltat) * deltaomega / 2.0 / PI;//実部を計算
            pr1 = ((*(Reez + k)) * c - (*(Imez + k)) * s) * exp(gamma * n*deltat) * deltaomega / 2.0 / PI;
            (*(AnsIm33 + n)) += ((*(Imez + k)) * c + (*(Reez + k)) * s) * exp(gamma * n*deltat)  * deltaomega / 2.0 / PI;//虚部を計算
            pr2 = ((*(Imez + k))* c + (*(Reez + k)) * s) * exp(gamma * n*deltat)  * deltaomega / 2.0 / PI;
        }
    }
    for (n = 0; n < (*(number + 0)); n++) {
        (*(ez + n)) = (*(AnsRe33 + n));
    }
    seikou = (*(time + 1));
    
    free(Re11);
    free(Im11);
    free(Re22);
    free(Im22);
    free(AnsRe33);
    free(AnsIm33);
    free(omega);
    free(Reet);
    free(Imet);
    free(ResGs);
    free(ImsGs);
    free(ResKs);
    free(ImsKs);
    free(Reez);
    free(Imez);
    
    return seikou;
}

void strain2stresskai(double *time, double t, double *epsall, int intpoint, double *tau, double *mater, double *sig, double *epsz)
{
    int i, j, n, *number;
    double *ex, *ey, *ez, *exy, *ekk, *dexdt, *deydt, *dexydt, *dekkdt, *timez;
    double sx, sy, sxy, skk, seikou;
    double ke, ge, *kmodulus, *gmodulus, *nmodulus, *nyumodulus;
    
    ex = (double *)malloc(sizeof(double)*NUM_OF_DATA);
    ey = (double *)malloc(sizeof(double)*NUM_OF_DATA);
    ez = (double *)malloc(sizeof(double)*NUM_OF_DATA);
    exy = (double *)malloc(sizeof(double)*NUM_OF_DATA);
    ekk = (double *)malloc(sizeof(double)*NUM_OF_DATA);
    dexdt = (double *)malloc(sizeof(double)*NUM_OF_DATA);
    deydt = (double *)malloc(sizeof(double)*NUM_OF_DATA);
    dexydt = (double *)malloc(sizeof(double)*NUM_OF_DATA);
    dekkdt = (double *)malloc(sizeof(double)*NUM_OF_DATA);
    kmodulus = (double *)malloc(sizeof(double)*NUM_OF_TERM);
    gmodulus = (double *)malloc(sizeof(double)*NUM_OF_TERM);
    nyumodulus = (double *)malloc(sizeof(double)*NUM_OF_TERM);
    nmodulus = (double *)malloc(sizeof(double)*NUM_OF_TERM);
    timez = (double *)malloc(sizeof(double)*NUM_OF_DATA);
    number = (int *)malloc(sizeof(int)*NUM_OF_TERM);
        
    ge = (*(mater + 0));
    ke = (*(mater + 1));
    for (i = 0; i < num_of_term; i++) {
        (*(gmodulus + i)) = (*(mater + 2 + 2 * i));
        (*(kmodulus + i)) = (*(mater + 2 + 2 * i + 1));
    }
    j = intpoint;
    for (i = 0; i < num_of_data; i++) {
        (*(ex + i)) = (*(epsall + 3 * j + nnip*i));
        (*(ey + i)) = (*(epsall + 3 * j + nnip*i + 1));
    }
    (*(number + 0)) = num_of_data;
    seikou = through_thickness_strain(time, timez, tau, gmodulus, kmodulus, ex, ey, ez, ge, ke, number, num_of_term);//板厚計算
    for (n = 0; n < (*(number + 0)); n++) {//元の時間間隔を保存
        (*(timez + n)) = (*(time + n));
    }
    for (i = 0; i < num_of_data; i++) {
        (*(ekk + i)) = (*(ex + i)) + (*(ey + i)) + (*(ez + i));//1108隠した
        (*(epsz + i)) = (*(ez + i));
    }
    for (i = 0; i < num_of_data; i++) {
        (*(ex + i)) = (*(ex + i)) - (*(ekk + i)) / 3.000000000000;
        (*(ey + i)) = (*(ey + i)) - (*(ekk + i)) / 3.000000000000;
        (*(exy + i)) = (*(epsall + 3 * j + nnip*i + 2)) / 2.000000000000;
    }
    differentiatelapack(time, ex, dexdt);
    differentiatelapack(time, ey, deydt);
    differentiatelapack(time, exy, dexydt);
    differentiatelapack(time, ekk, dekkdt);
    sx = 2 * integrate3(time, ex, dexdt, num_of_data, tau, gmodulus, ge, num_of_term, t, temp);
    sy = 2 * integrate3(time, ey, deydt, num_of_data, tau, gmodulus, ge, num_of_term, t, temp);
    sxy = 2 * integrate3(time, exy, dexydt, num_of_data, tau, gmodulus, ge, num_of_term, t, temp);
    skk = 3 * integrate3(time, ekk, dekkdt, num_of_data, tau, kmodulus, ke, num_of_term, t, temp);//shotaro
    (*(sig + 0)) = sx + skk / 3;
    (*(sig + 1)) = sy + skk / 3;
    (*(sig + 2)) = sxy;
    (*(sig + 3)) = skk;
    
    free(ex);
    free(ey);
    free(ez);
    free(exy);
    free(ekk);
    free(dexdt);
    free(deydt);
    free(dexydt);
    free(dekkdt);
    free(kmodulus);
    free(gmodulus);
    free(nyumodulus);
    free(nmodulus);
    free(timez);
    free(number);
}

void strain2stresskai2(double *time, double t, double *epsall, int intpoint, double *tau, double *mater, double *sig, double *epsz)
{
    int i, j, n, *number;
    double *ex, *ey, *ez, *exy, *ekk, *dexdt, *deydt, *dexydt, *dekkdt, *timez;
    double sx, sy, sxy, skk, seikou;
    double ke, ge, *kmodulus, *gmodulus, *nmodulus, *nyumodulus;

    ex = (double *)malloc(sizeof(double)*NUM_OF_DATA);
    ey = (double *)malloc(sizeof(double)*NUM_OF_DATA);
    ez = (double *)malloc(sizeof(double)*NUM_OF_DATA);
    exy = (double *)malloc(sizeof(double)*NUM_OF_DATA);
    ekk = (double *)malloc(sizeof(double)*NUM_OF_DATA);
    dexdt = (double *)malloc(sizeof(double)*NUM_OF_DATA);
    deydt = (double *)malloc(sizeof(double)*NUM_OF_DATA);
    dexydt = (double *)malloc(sizeof(double)*NUM_OF_DATA);
    dekkdt = (double *)malloc(sizeof(double)*NUM_OF_DATA);
    kmodulus = (double *)malloc(sizeof(double)*NUM_OF_TERM);
    gmodulus = (double *)malloc(sizeof(double)*NUM_OF_TERM);
    nyumodulus = (double *)malloc(sizeof(double)*NUM_OF_TERM);
    nmodulus = (double *)malloc(sizeof(double)*NUM_OF_TERM);
    timez = (double *)malloc(sizeof(double)*NUM_OF_DATA);
    number = (int *)malloc(sizeof(int)*NUM_OF_TERM);
    
    ge = (*(mater + 0));
    ke = (*(mater + 1));
    for (i = 0; i < num_of_term; i++) {
        (*(gmodulus + i)) = (*(mater + 2 + 2 * i));
        (*(kmodulus + i)) = (*(mater + 2 + 2 * i + 1));
    }
    j = intpoint;
    for (i = 0; i < num_of_data; i++) {
        (*(ex + i)) = (*(epsall + 3 * j + nnip*i));
        (*(ey + i)) = (*(epsall + 3 * j + nnip*i + 1));
    }
    (*(number + 0)) = num_of_data;
    seikou = through_thickness_strain(time, timez, tau, gmodulus, kmodulus, ex, ey, ez, ge, ke, number, num_of_term);
    for (n = 0; n < (*(number + 0)); n++) {//元の時間間隔を保存
        (*(timez + n)) = (*(time + n));
    }

    for (i = 0; i < num_of_data; i++) {
        (*(ekk + i)) = (*(ex + i)) + (*(ey + i)) + (*(ez + i));
    }

    for (i = 0; i < num_of_data; i++) {
        (*(ex + i)) = (*(ex + i)) - (*(ekk + i))/3.000000000000;
        (*(ey + i)) = (*(ey + i)) - (*(ekk + i))/3.000000000000;
        (*(exy + i)) = (*(epsall + 3 * j + nnip*i + 2))/2.000000000000;
    }
    differentiatelapack(time, ex, dexdt);
    differentiatelapack(time, ey, deydt);
    differentiatelapack(time, exy, dexydt);
    differentiatelapack(time, ekk, dekkdt);
    sx = 2 * integrate3(time, ex, dexdt, num_of_data, tau, gmodulus, ge, num_of_term, t, temp);
    sy = 2 * integrate3(time, ey, deydt, num_of_data, tau, gmodulus, ge, num_of_term, t, temp);
    sxy = 2 * integrate3(time, exy, dexydt, num_of_data, tau, gmodulus, ge, num_of_term, t, temp);
    skk = 3 * integrate3(time, ekk, dekkdt, num_of_data, tau, kmodulus, ke, num_of_term, t, temp);
    (*(sig + 0)) = sx + skk/3.000000000000;
    (*(sig + 1)) = sy + skk/3.000000000000;
    (*(sig + 2)) = sxy;
    (*(sig + 3)) = skk;

    free(ex);
    free(ey);
    free(ez);
    free(exy);
    free(ekk);
    free(dexdt);
    free(deydt);
    free(dexydt);
    free(dekkdt);
    free(nyumodulus);
    free(nmodulus);
    free(kmodulus);
    free(gmodulus);
    free(timez);
}

void fullfieldstrain2stresskai(double *time, double *epsall, double *tau, double *gmodulus, double *kmodulus, double ge, double ke, double *sigmaall, double *zeps)
{
    int i, j, n;
    int *number;
    double t, seikou;
    double *ex, *ey, *ez, *exy, *ekk, *dexdt, *deydt, *dexydt, *dekkdt;
    double *sigmax, *sigmay, *tauxy, *timez;
    double sx, sy, sxy, skk;

    ex = (double *)malloc(sizeof(double)*NUM_OF_DATA * 10);
    ey = (double *)malloc(sizeof(double)*NUM_OF_DATA * 10);
    ez = (double *)malloc(sizeof(double)*NUM_OF_DATA * 10);
    exy = (double *)malloc(sizeof(double)*NUM_OF_DATA * 10);
    ekk = (double *)malloc(sizeof(double)*NUM_OF_DATA * 10);
    dexdt = (double *)malloc(sizeof(double)*NUM_OF_DATA * 10);
    deydt = (double *)malloc(sizeof(double)*NUM_OF_DATA * 10);
    dexydt = (double *)malloc(sizeof(double)*NUM_OF_DATA * 10);
    dekkdt = (double *)malloc(sizeof(double)*NUM_OF_DATA * 10);
    sigmax = (double *)malloc(sizeof(double)*NUM_OF_DATA * 10);
    sigmay = (double *)malloc(sizeof(double)*NUM_OF_DATA * 10);
    tauxy = (double *)malloc(sizeof(double)*NUM_OF_DATA * 10);
    timez = (double *)malloc(sizeof(double)*NUM_OF_DATA);
    number = (int *)malloc(sizeof(int)*NUM_OF_TERM);
    
    for (j = 0; j < nip; j++) {
        for (i = 0; i < num_of_data; i++) {
            (*(ex + i)) = (*(epsall + 3 * j + nnip*i));
            (*(ey + i)) = (*(epsall + 3 * j + nnip*i + 1));
        }
        (*(number + 0)) = num_of_data;
        seikou = through_thickness_strain(time, timez, tau, gmodulus, kmodulus, ex, ey, ez, ge, ke, number, num_of_term);//板厚計算
        for (n = 0; n < (*(number + 0)); n++) {
            (*(timez + n)) = (*(time + n));
        }
        for (i = 0; i < num_of_data; i++) {
            (*(ekk + i)) = (*(ex + i)) + (*(ey + i)) + (*(ez + i));
        }
        for (i = 0; i < num_of_data; i++) {
            t = (*(time + i));
            (*(ex + i)) = (*(ex + i)) - ((*(ekk + i))/3.000000000000);
            (*(ey + i)) = (*(ey + i)) - ((*(ekk + i))/3.000000000000);
            (*(exy + i)) = (*(epsall + 3 * j + nnip*i + 2))/2.000000000000;
        }
        differentiatelapack(time, ex, dexdt);
        differentiatelapack(time, ey, deydt);
        differentiatelapack(time, exy, dexydt);
        differentiatelapack(time, ekk, dekkdt);
        t = 0;
        for (i = 0; i < num_of_data; i++) {
            t = (*(time + i));
            sx = 2 * integrate3(time, ex, dexdt, num_of_data, tau, gmodulus, ge, num_of_term, t, temp);
            sy = 2 * integrate3(time, ey, deydt, num_of_data, tau, gmodulus, ge, num_of_term, t, temp);
            sxy = 2 * integrate3(time, exy, dexydt, num_of_data, tau, gmodulus, ge, num_of_term, t, temp);
            skk = 3 * integrate3(time, ekk, dekkdt, num_of_data, tau, kmodulus, ke, num_of_term, t, temp);
            (*(zeps + j + nip*i)) = (*(ez + i));
            (*(sigmaall + 3 * j + nnip*i)) = sx + skk/3.000000000000;
            (*(sigmaall + 3 * j + nnip*i + 1)) = sy + skk/3.000000000000;
            (*(sigmaall + 3 * j + nnip*i + 2)) = sxy;
        }
    }
    free(ex);
    free(ey);
    free(ez);
    free(exy);
    free(ekk);
    free(dexdt);
    free(deydt);
    free(dexydt);
    free(dekkdt);
    free(sigmax);
    free(sigmay);
    free(tauxy);
    free(timez);
    free(number);
}

void zero(double *a, int n)
{
    double *al;
    for (al = a + n; a < al;){
        *(a++) = 0.0;
    }
}

int mkmatrix1(double *tau, double *mater, double *pt, int *el, double *time, double *uuord, double *epsall, double *vuuord, double *vftime, int index, int indexmodel, int indextime, double *mat_a, double *mat_f, double *sigmaallkai, double *pointskai)
{
    double *elpt, *eldisp, *bmat, *eps, *sig, *ez, *materkai;
    double *eldisp1, *eldisp2, *eldisp3, *eldisp4;
    double *eps1, *eps2, *eps3, *eps4;
    double *shape, *deriv, *cart, *g_po, *g_we;
    double xcoord, ycoord;
    double px, py, det;
    int i, l, k, j, ig, jg;
    double epsx1, epsy1, gammaxy1;
    double sx, sy, txy;
    double dsx, dsy, dtxy, dsx1, dsy1, dtxy1, dsx2, dsy2, dtxy2, dskk1, dskk2, dskk;
    double t;
    g_po = (double *)malloc(sizeof(double)*GP);
    g_we = (double *)malloc(sizeof(double)*GP);
    cart = (double *)malloc(sizeof(double)*DM*E8);
    shape = (double *)malloc(sizeof(double)*E8);
    deriv = (double *)malloc(sizeof(double)*DM*E8);
    elpt = (double *)malloc(sizeof(double)*DM*E8);
    eldisp = (double *)malloc(sizeof(double)*DM*E8);
    eldisp1 = (double *)malloc(sizeof(double)*DM*E8);
    eldisp2 = (double *)malloc(sizeof(double)*DM*E8);
    eldisp3 = (double *)malloc(sizeof(double)*DM*E8);
    eldisp4 = (double *)malloc(sizeof(double)*DM*E8);
    bmat = (double *)malloc(sizeof(double) * 16 * 3);
    eps = (double *)malloc(sizeof(double) * 3);
    sig = (double *)malloc(sizeof(double) * 4);
    eps1 = (double *)malloc(sizeof(double) * 3);
    eps2 = (double *)malloc(sizeof(double) * 3);
    eps3 = (double *)malloc(sizeof(double) * 3);
    eps4 = (double *)malloc(sizeof(double) * 3);
    ez = (double *)malloc(sizeof(double) * NUM_OF_DATA * 10);
    materkai = (double *)malloc(sizeof(double)*(NUM_OF_TERM * 3 + 2));

    make_gp(g_po, g_we);
    for (j = 0, l = 0; l < ne; l++, el += 9) {
        make_pt(elpt, el, pt);
        /* elpt: Coordinate of elements */
        for (i = 0; i < E8; i++) {
            k = *(el + i);
            *(eldisp1 + 2 * i) = (*(vuuord + 2 * k + index*nn));
            *(eldisp1 + 2 * i + 1) = (*(vuuord + 2 * k + 1 + index*nn));
        }
        /*printf("mkmatrix");*/
        for (ig = 0; ig < GP; ig++) {
            for (jg = 0; jg < GP; jg++) {
                px = (*(g_po + ig));
                py = (*(g_po + jg));
                make_sp(shape, deriv, px, py);
                jacob(cart, &det, elpt, deriv); /* cart 2x8 */
                make_b(bmat, cart);
                /*matx(eps,bmat,eldisp,3,16,1);*///ここをあけてる
                matx(eps1, bmat, eldisp1, 3, 16, 1);
                xcoord = (*(shape))*(*(elpt)) + (*(shape + 1))*(*(elpt + 2)) + (*(shape + 2))*(*(elpt + 4)) + (*(shape + 3))*(*(elpt + 6)) + (*(shape + 4))*(*(elpt + 8)) + (*(shape + 5))*(*(elpt + 10)) + (*(shape + 6))*(*(elpt + 12)) + (*(shape + 7))*(*(elpt + 14));
                ycoord = (*(shape))*(*(elpt + 1)) + (*(shape + 1))*(*(elpt + 3)) + (*(shape + 2))*(*(elpt + 5)) + (*(shape + 3))*(*(elpt + 7)) + (*(shape + 4))*(*(elpt + 9)) + (*(shape + 5))*(*(elpt + 11)) + (*(shape + 6))*(*(elpt + 13)) + (*(shape + 7))*(*(elpt + 15));
                *(pointskai + 2 * j + 0) = xcoord;
                *(pointskai + 2 * j + 1) = ycoord;
                t = (*(vftime + indexmodel));
                strain2stresskai(time, t, epsall, j, tau, mater, sig, ez);
                sx = (*(sig + 0));
                /*sx = 0;*/
                sy = (*(sig + 1));
                txy = (*(sig + 2));
                *(sigmaallkai + 3 * j + 0) = sx;
                *(sigmaallkai + 3 * j + 1) = sy;
                *(sigmaallkai + 3 * j + 2) = txy;
                epsx1 = (*(eps1 + 0));
                epsy1 = (*(eps1 + 1));
                gammaxy1 = (*(eps1 + 2));
                (*(mat_f + index + num_of_vf*indexmodel)) -= (sx*epsx1 + sy*epsy1 + txy*gammaxy1)*det*(*(g_we + ig))*(*(g_we + jg));
                for (i = 0; i < num_of_term * 2 + 2; i++) {
                    (*(materkai + i)) = (*(mater + i));
                }
                for (i = 0; i < num_of_term * 2 + 2; i++) {
                    (*(mater + i)) = (*(mater + i)) - DMODULUS/2;
                    strain2stresskai2(time, t, epsall, j, tau, mater, sig, ez);
                    dsx1 = (*(sig + 0));
                    dsy1 = (*(sig + 1));
                    dtxy1 = (*(sig + 2));
                    dskk1 = (*(sig + 3));//shotaro
                    (*(mater + i)) = (*(mater + i)) + DMODULUS;
                    strain2stresskai2(time, t, epsall, j, tau, mater, sig, ez);
                    dsx2 = (*(sig + 0));
                    dsy2 = (*(sig + 1));
                    dtxy2 = (*(sig + 2));
                    dskk2 = (*(sig + 3));//shotaro
                    dsx = (dsx2 - dsx1) / DMODULUS;
                    dsy = (dsy2 - dsy1) / DMODULUS;
                    dtxy = (dtxy2 - dtxy1) / DMODULUS;
                    dskk = (dskk2 - dskk1) / DMODULUS;
                    (*(mater + i)) = (*(mater + i)) - DMODULUS/2.000000000000;
                    (*(mat_a + (index + num_of_vf*indexmodel) + (num_of_vf*num_of_model)*i)) += (dsx*epsx1 + dsy*epsy1 + dtxy*gammaxy1)*det*(*(g_we + ig))*(*(g_we + jg));
                }
                j++;
            }
        }
    }
    printf("df/dge\t   df/dke\t   df/dg1\t   df/dk1\t   df/dg2\t   df/dk2\n");
    for (i = 0; i < num_of_term * 2 + 2; i++) {
        printf("%lf\t", (*(mat_a + (index + num_of_vf*indexmodel) + (num_of_vf*num_of_model) * i)));
    }
    printf("\n");

    free(g_po);
    free(g_we);
    free(cart);
    free(shape);
    free(deriv);
    free(elpt);
    free(eldisp);
    free(eldisp1);
    free(eldisp2);
    free(eldisp3);
    free(eldisp4);
    free(bmat);
    free(eps);
    free(sig);
    free(eps1);
    free(eps2);
    free(eps3);
    free(eps4);
    free(ez);
    free(materkai);
    
    return 0;
}

int mkmatrix2(double *pt, double *pf, double *vuuord, int index, int indexmodel, double *mat_b, double *mat_vf)
{
    printf("u%d ∫σ＊ε　∫T＊u\n", index + 1);
    printf("%lf\t%lf\n", (*(mat_b + index + indexmodel*num_of_vf)), (*(mat_vf + index + indexmodel*num_of_vf)));
    printf("\n");
    (*(mat_b + index + indexmodel*num_of_vf)) -= (*(mat_vf + index + indexmodel*num_of_vf));
    /*(*(mat_b + index + indexmodel*num_of_vf)) -= (*(mat_b + index + indexmodel*num_of_vf));*/
    return 0;
}


void mkvf(double *pt, double *pf, double *vuuord, int index, int indexmodel, double *mat_vf)
{
    int k;
    double x, y, u, v;
    double vf1;

    for (vf1 = 0, k = 0; k < np; k++) {
        if ((*(pf + 2 * k)) != 0 || (*(pf + 2 * k + 1)) != 0) {
            x = (*(pt + 2 * k));
            y = (*(pt + 2 * k + 1));
            u = (*(vuuord + 2 * k + index*nn));
            v = (*(vuuord + 2 * k + 1 + index*nn));
            vf1 -= u*(*(pf + 2 * k + indexmodel*nn)) + v*(*(pf + 2 * k + 1 + indexmodel*nn));
        }
    }
    (*(mat_vf + index + indexmodel*num_of_vf)) = vf1;
    //return 0;
}

void solvelms(double *mat_a, double *mat_b)
{
    int i;
    int info, nrhs = 1, lwork;
    int n, m, lda, ldb;
    double *work;
    char trans = 'N';
    
    work = (double *)malloc(sizeof(double)*NVF*NUM_OF_MODEL);
    m = (int)num_of_vf*(int)num_of_model;
    n = (int)num_of_term * 2 + 2;
    lda = m;
    ldb = m;
    lwork = m + n;
    
    dgels_(&trans, &m, &n, &nrhs, mat_a, &lda, mat_b, &ldb, work, &lwork, &info);
    for (i = 0; i < ldb; i++) {
        printf("mat_b %lf\n", (*(mat_b + i)));
    }
    free(work);
}

int main(int argc, char *argv[]){
    double *pt, *pftmp, *pf, *mat, *mater;
    int *el, *ib;
    int *ibn;
    double *disp;
    double *uuord, *vuuord;
    double *xzahyo, *yzahyo;
    int i, j, k, l, a, z;
    double *mat_b, *fk, *mat_a, *mat_vf;
    char filen[FL];
    char file[FL];
    char filemodel[FL];
    char filedata1[FL];
    char filedata5[FL];
    char filedata6[FL];
    //char filedata7[FL];
    FILE *fp;
    FILE *fp1;
    FILE *fp2;
    FILE *fp3;
    //FILE *fp4;
    double *dispall, *eps, *epsall, *points, *pointskai, *sigmaall, *sigmaallkai, *stressall, *vdispall, *zeps;
    double *tau, *gcoef, *kcoef, *ncoef, ke, ge;
    double t, *time, *vftime, sum, *vftimekai;

    pt = (double *)malloc(sizeof(double)*NN);
    pftmp = (double *)malloc(sizeof(double)*NN);
    pf = (double *)malloc(sizeof(double)*NN*NUM_OF_MODEL);
    mat = (double *)malloc(sizeof(double) * 3);
    mater = (double *)malloc(sizeof(double)*(NUM_OF_TERM * 3 + 2));
    el = (int *)malloc(sizeof(int)*NE*(E8 + 1));
    ib = (int *)malloc(sizeof(int)*NB * 3);
    ibn = (int *)malloc(sizeof(int)*NB);
    disp = (double *)malloc(sizeof(double)*NN);
    uuord = (double *)malloc(sizeof(double)*NN);
    xzahyo = (double *)malloc(sizeof(double)*NN);
    yzahyo = (double *)malloc(sizeof(double)*NN);
    uuord = (double *)malloc(sizeof(double)*NN);
    vuuord = (double *)malloc(sizeof(double)*NN*NVF);
    fk = (double *)malloc(sizeof(double)*NN);
    dispall = (double *)malloc(sizeof(double)*NN*NUM_OF_DATA);
    vdispall = (double *)malloc(sizeof(double)*NN*NUM_OF_DATA);
    epsall = (double *)malloc(sizeof(double)*NE*(E8 + 1) * 3 * NUM_OF_DATA * 10);
    points = (double *)malloc(sizeof(double)*NE*(E8 + 1) * 2 * NUM_OF_DATA + 100);
    pointskai = (double *)malloc(sizeof(double)*NE*(E8 + 1) * 2 * NUM_OF_DATA);
    eps = (double *)malloc(sizeof(double)*NE*(E8 + 1) * 3 * 10);
    zeps = (double *)malloc(sizeof(double)*NE*(E8 + 1) * 3 * NUM_OF_DATA * 10);
    time = (double *)malloc(sizeof(double)*NUM_OF_DATA);
    gcoef = (double *)malloc(sizeof(double)*NUM_OF_TERM);
    kcoef = (double *)malloc(sizeof(double)*NUM_OF_TERM);
    ncoef = (double *)malloc(sizeof(double)*NUM_OF_TERM);
    tau = (double *)malloc(sizeof(double)*NUM_OF_TERM);
    sigmaall = (double *)malloc(sizeof(double)*NE*(E8 + 1) * 3 * NUM_OF_DATA * 10);
    stressall = (double *)malloc(sizeof(double)*NE*(E8 + 1) * 3 * NUM_OF_DATA * 10);
    sigmaallkai = (double *)malloc(sizeof(double)*NE*(E8 + 1) * 3 * NUM_OF_DATA * 10);
    mat_vf = (double *)malloc(sizeof(double)*NVF*NUM_OF_MODEL);
    vftime = (double *)malloc(sizeof(double)*NVF);
    vftimekai = (double *)malloc(sizeof(double)*NVF);
    mat_a = (double *)malloc(sizeof(double)*NUM_OF_UNKNOWN*NVF*NUM_OF_MODEL);
    mat_b = (double *)malloc(sizeof(double)*NVF*NUM_OF_MODEL);

    if (argc < 2) {
        printf("File-in: ");
        scanf("%s", filen);
    }
    else {
        sprintf(filen, "%s", argv[1]);
    }
    if (NULL == (fp = fopen(filen, "r"))) {
        printf(" Cannot open file\n");
        exit(0);
    }
    fscanf(fp, "%s", file);
    printf("%s\n", file);
    num_of_term = readmodulus(file, tau, gcoef, kcoef, &ge, &ke);
    printf("%d\n", num_of_term);
    //printf("%lf\t%lf\t%lf\n%lf\t%lf\t%lf\n", *tau, *gcoef, *kcoef, (*(tau+1)), (*(gcoef+1)), (*(kcoef+1)));
    fscanf(fp, "%d", &num_of_model);
    for (i = 0; i < num_of_model; i++) {
        fscanf(fp, "%lf %s", &t, filemodel);
        printf("%s\n", filemodel);
        (*(vftime + i)) = t;
        readmodel(filemodel, pt, pftmp, mat, el, ib, ibn, disp);
        for (j = 0; j < nn; j++) {
            (*(pf + j + i*nn)) = (*(pftmp + j));
        }
    }
    fscanf(fp, "%lf", &temp);
    printf("%lg\n", temp);
    fscanf(fp, "%d", &num_of_vf);
    printf("num_of_vf: %d\n", num_of_vf);
    for (i = 0; i < num_of_vf; i++) {
        fscanf(fp, "%s", filedata1);
        printf("vcoord: %s\n", filedata1);
        read_displacement(filedata1, vuuord);
        for (k = 0; k < nn; k++) {
            (*(vdispall + k + i*nn)) = (*(vuuord + k));
            (*(vuuord + k)) = (*(vuuord + k));//* 33.15仮想変位をpixelに変換
        }
    }
    for (j = 0; j < num_of_model; j++) {
        for (i = 0; i < num_of_vf; i++) {
            mkvf(pt, pf, vdispall, i, j, mat_vf);
        }
    }
    fscanf(fp, "%d", &num_of_data);
    printf("num_of_data %d\n", num_of_data);
    for (i = 0; i < num_of_data; i++) {
        fscanf(fp, "%lg %s", &t, filedata1);
        (*(time + i)) = t;
        read_displacement(filedata1, uuord);
        strain(uuord, el, pt, points, eps);
        for (j = 0; j < nn; j++) {
            (*(dispall + j + i*nn)) = (*(uuord + j));
        }
        for (j = 0; j < nnip; j++) {
            (*(epsall + j + i*nnip)) = (*(eps + j));
        }
        /*sprintf(filedata7, "strain%d.txt", i);
        fp4 = fopen(filedata7, "w");
        if (NULL == (fp4 = fopen(filedata7, "w"))) {
            printf("\7\n Cannot open strain file..\n\n");
            exit(1);
        }
        for (j = 0; j < nip + num_of_data; j++) {
            fprintf(fp4, "%lg\t%lg\t%lg\t%lg\t%lg\n", (*(points + 2 * j)), (*(points + 2 * j + 1)), (*(epsall + 3 * j + nnip*i)), (*(epsall + 3 * j + nnip*i + 1)), (*(epsall + 3 * j + nnip*i + 2)));
        }
        fclose(fp4);*/
    }
    fclose(fp);
    (*(mater + 0)) = ge;
    (*(mater + 1)) = ke;
    //printf("%lf\t%lf ", ge, ke);
    printf("ge: %lg\t", (*(mater + 0)));
    printf("ke: %lg\n", (*(mater + 1)));
    for (i = 0; i < num_of_term; i++) {
        (*(mater + 2 + 2 * i)) = (*(gcoef + i));
        (*(mater + 2 + 2 * i + 1)) = (*(kcoef + i));
        printf("g%d:\t%lg\tk%d:\t%lg\n", i + 1, (*(mater + 2 + 2 * i)), i + 1, (*(mater + 2 + 2 * i + 1)));
    }
    fullfieldstrain2stresskai(time, epsall, tau, gcoef, kcoef, ge, ke, sigmaall, zeps);
    printf("要素数%d\n", nip / 9);
    for (i = 0; i < num_of_data; i++) {
        sprintf(filedata5, "stress%d.txt", i);
        fp1 = fopen(filedata5, "w");
        if (NULL == (fp1 = fopen(filedata5, "w"))) {
            printf("\7\n Cannot open stress file..\n\n");
            exit(1);
        }
        for (j = 0; j < nip + num_of_data; j++) {
            fprintf(fp1, "%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", (*(points + 2 * j)), (*(points + 2 * j + 1)), (*(zeps + j + nip*i)), (*(sigmaall + 3 * j + nnip*i)), (*(sigmaall + 3 * j + nnip*i + 1)), (*(sigmaall + 3 * j + nnip*i + 2)));
        }
        fclose(fp1);
    }
    printf("fullfieldstress計算終了 データ数：%d\n", i);
    printf("\n");
    for (a = 0; a < 500; a++) {
        zero(mat_a, NUM_OF_UNKNOWN*NVF*NUM_OF_MODEL);
        zero(mat_b, NVF*NUM_OF_MODEL);
        for (j = 0; j < num_of_model; j++) {
            for (i = 0; i < num_of_vf; i++) {
                for (k = 0, l = 0; l < num_of_data; l++) {
                    if ((*(time + l)) == (*(vftime + j))) {
                        k = l;
                        printf("t=%lg\n", (*(time + l)));
                        break;
                    }
                }
                mkmatrix1(tau, mater, pt, el, time, uuord, epsall, vdispall, vftime, i, j, k, mat_a, mat_b, sigmaallkai, pointskai);
                mkmatrix2(pt, pf, vdispall, i, j, mat_b, mat_vf);
            }
            sprintf(filedata6, "stresssekibunntenn.txt");
            if (NULL == (fp2 = fopen(filedata6, "w"))) {
                printf("\7\n Cannot open stresssekibunntenn file..\n\n");
                exit(1);
            }
            for (z = 0; z < nip * 2; z++) {
                fprintf(fp2, "%lf\t%lf\t%lf\t%lf\t%lf\n", (*(pointskai + 2 * z)), (*(pointskai + 2 * z + 1)), (*(sigmaallkai + 3 * z)), (*(sigmaallkai + 3 * z + 1)), (*(sigmaallkai + 3 * z + 2)));
            }
            fclose(fp2);
        }
        solvelms(mat_a, mat_b);
        printf("mat_b:修正項\n");
        printf("\n");
        printf("繰り返し数：%d\n", a);
        printf("\n");
        if (NULL == (fp3 = fopen("output1.txt", "a+"))) {
                printf("\7\n No search for outputfile\n\n");
                exit(1);
        }
        for (sum = 0, i = 0; i < 2 * num_of_term + 2; i++) {
            (*(mater + i)) += (*(mat_b + i));
            printf("mater%d: %g\n", i, (*(mater + i)));
            sum += (*(mat_b + i))*(*(mat_b + i));
        }
        printf("--------------------------------------");
        printf("\n");
    if (sqrt(sum) < SHIKII) {
            printf("Converged! %d\n", a);
            printf("\n");
            break;
        }
        fprintf(fp3, "%s\t%s\t%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", filen, file, a, (*(mater)), (*(mater + 1)), (*(mater + 2)),(*(mater + 3)),(*(mater + 4)),(*(mater + 5)));
        fclose(fp3);
    }
    printf("ge\tke\n");
    for (i = 0; i < 2 * num_of_term - 2; i++) {
        printf("tau%d\tg%d\tk%d\n", i + 1, i + 1, i + 1);
    }
        printf("%lf\t%lf\n", (*(mater + 0)), (*(mater + 1)));
    for (i = 1; i < 2 * num_of_term - 1; i++) {
        printf("%lf\t%lf\t%lf\n", (*(tau + i - 1)), (*(mater + 2 * i)), (*(mater + 2 * i + 1)));//20200826shotaro隠した
    }
    
    free(pt);
    free(pftmp);
    free(pf);
    free(mat);
    free(mater);
    free(el);
    free(ib);
    free(ibn);
    free(disp);
    free(uuord);
    free(vuuord);
    free(fk);
    free(eps);
    free(dispall);
    free(vdispall);
    free(points);
    free(epsall);
    free(zeps);
    free(tau);
    free(kcoef);
    free(gcoef);
    free(sigmaall);
    free(sigmaallkai);
    free(stressall);
    free(mat_vf);
    free(vftime);
    free(vftimekai);
    free(mat_a);
    free(mat_b);
        
    return 0;
}

