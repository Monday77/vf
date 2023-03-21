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
#include "mattx.h" /* マトリックスの転置のかけ算 */

#define FL 500
#define DM 2  /* 次元 */
#define NP 7000 /* 節点個数の上限 */
#define NM 2  /* 材料数の上限 */
#define NE 3000  /* 要素数の上限 */
#define NB 200 /* 拘束節点数の上限 */
#define NN DM*NP
#define BW 30 /* バンド幅の上限 */
#define E8 8 /* 1要素の節点数 */
#define E16 16 /* EKマトリックスに用いる16 */
#define M3 3 /* Dマトリックスに用いる3  */
#define GP 3 /* ガウスポイントの数  */

#define Ncon 1 /* 材料定数の未知数 */
#define NVF 50
#define ND 100 /* 変位データの数 */
#define SHV 0.001
#define DE 1.0
#define DNU 0.01
#define DMODULUS (0.0001)//0.005
#define DN 0.01

#define SHIKII (0.0001)

#define PI (3.14159265359)
#define M (5) /* ??????gp */
#define N (3) /* ??????gp */
#define C1 (8.86)
#define C2 (101.6)
#define TREF (303.0)
#define NUM_OF_TERM (20)
#define NUM_OF_DATA (2000)//各ひずみ成分と板厚方向ひずみの変数サイズの設定に使われてる
#define NUM_OF_MODEL (50)//20
#define POISSON (0.304276)
//#define TEMP (263.0)//温度変更項
#define NUM_OF_UNKNOWN (50)
#define NUM_OF_DATA1 (500)
#define NC (1)

int _stack=12000;

int num_of_data;
int num_of_model;
int num_of_vf;
int num_of_term;
int np, nm, ne, nb, nn, nbc;
int nip, nnip;//np : number of integration point
double temp;
double e,nu;

double readmodel(char *file1, double *pt, double *pf, double *mat, int *el, int *ib, int *ibn, double *disp, int *ibc, int *ibnc)
{
    FILE *fp;
    int i, j;
    
    if( NULL == (fp = fopen(file1, "r"))) {
        printf(" Cannot open models--file \n");
        exit(1);
    }

  /*  printf("%s\n",filen);*/
    fscanf(fp, "%d %d %d %d", &np, &nm, &ne, &nb);
    //printf("%d\t %d\t %d\t %d\n", np, nm, ne, nb);
    //printf("\n");
    
    nn = np * DM;
    nip = ne * (E8+1);
    nnip = 3 * nip;

  
    for(i = 0; i < np; i++, pt += 2, pf += 2){
        fscanf(fp, "%lf %lf %lf %lf", pt, pt + 1, pf, pf + 1);
        //printf("%lf\t %lf\t %lf\t %lf\n", *pt, *(pt + 1), *pf, *(pf + 1));
      /* 節点座標と節点力 */
    }
    //printf("\n");


    for(i = 0; i < nm; i++, mat += 3){
        fscanf(fp, "%lf %lf %lf", mat, mat + 1, mat + 2);
       // printf("%lf\t %lf\t %lf \n", *mat, *(mat + 1), *(mat + 2));
      /* ヤング率，ポアソン比，板厚 */
    }
    printf("\n");
    
    for(i = 0; i < ne; i++){
        for(j = 0; j < E8; j++){
            fscanf(fp, "%d", el++);
        }
        fscanf(fp, "%d", el++);
    }
   // printf("\n");

  /* 要素の構成 */

    for(i = 0; i < nb; i++, ibn += 1, ib += 2, disp += 2){
        fscanf(fp, "%d %d %d %lf %lf", ibn, ib, ib + 1, disp, disp + 1);
        // printf("%d\t %d\t %d\t %lf\t %lf\n", *ibn, *ib, *(ib + 1), *disp, *(disp + 1));
    }
   /* 固定条件 */
    fscanf(fp, "%d", &nbc);
   // printf("%d\n", nbc);

    for(i = 0; i < nbc; i++, ibnc += 1, ibc += 2){
        fscanf(fp, "%d %d %d", ibnc, ibc, ibc + 1);
     // printf("%d\t %d\t %d\n", *ibnc, *ibc, *(ibc + 1));
    }
    fclose(fp);
    return i;
}

double readmodulus(char *file2, double *time, double *gcoef, double *kcoef, double *mater)
//read K, G
{
 
    FILE* fp;
    int i = 0;
    
    if (NULL == (fp = fopen(file2, "r"))) {
        printf("Cannot open modulus file..\n\n");
        exit(1);
    }
    fscanf(fp, "%lf", mater);
    printf("Ge = %lg\t  ", *mater);
    fscanf(fp, "%lf", (mater + 1));
    printf("Ke = %lg\n", *(mater + 1));
    
    while (3 == fscanf(fp, "%lf %lf %lf", time, gcoef, kcoef)) {
        printf(" [%d] \t: %lg\t%lg\t%lg\n", i + 1, *time, *gcoef, *kcoef);
        //without ++, data saved in the same address
        time++;
        gcoef++;
        kcoef++;
        i++;
    }
    fclose(fp);
    num_of_term = i;
    return i;
}

double read_dispdata(char *file3, double *uuord)
{
    FILE *fp;
    int i;
    double a, b;

    if (NULL == (fp = fopen(file3, "r"))) {
        printf(" Cannot open displacement-file\n");
        exit(1);
    }
    for (i = 0; i < np; i++, uuord += 2) {
        fscanf(fp, "%lf %lf %lf %lf", &a, &b, uuord, uuord + 1);
    }
    fclose(fp);
    return i;
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
}

void make_gp(double *g_po, double *g_we)
{
    *g_po ++= -0.7746;
    *g_po ++= 0.0;
    *g_po = 0.7746;
    *g_we ++= 0.5556;
    *g_we ++= 0.8889;
    *g_we = 0.5556;
}

void make_pt(double *elpt, int *ell, double *pt)
{
    int l;
    for(l = 0; l < E8; l++, ell++){
        *elpt ++= *(pt+ *ell*DM);
        *elpt ++= *(pt+ *ell*DM+1);
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

    *shape++=(-1+st+ss+tt-sst-stt)/4;
    *shape++=(1-t-ss+sst)/2;
    *shape++=(-1-st+ss+tt-sst+stt)/4;
    *shape++=(1+s-tt-stt)/2;
    *shape++=(-1+st+ss+tt+sst+stt)/4;
    *shape++=(1+t-ss-sst)/2;
    *shape++=(-1-st+ss+tt+sst-stt)/4;
    *shape=(1-s-tt+stt)/2;

    *deriv++=(t+s2-st2-tt)/4;
    *deriv++=-s+st;
    *deriv++=(-t+s2-st2+tt)/4;
    *deriv++=(1-tt)/2;
    *deriv++=(t+s2+st2+tt)/4;
    *deriv++=-s-st;
    *deriv++=(-t+s2+st2-tt)/4;
    *deriv++=(-1+tt)/2;
    *deriv++=(s+t2-ss-st2)/4;
    *deriv++=(-1+ss)/2;
    *deriv++=(-s+t2-ss+st2)/4;
    *deriv++=-t-st;
    *deriv++=(s+t2+ss+st2)/4;
    *deriv++=(1-ss)/2;
    *deriv++=(-s+t2+ss-st2)/4;
    *deriv=-t+st;
}
void jacob(double *cart, double *djacb, double *elpt, double *deriv)
{

    double *jacm,*jaci;

    jacm=(double *)malloc(sizeof(double)*DM*DM);
    jaci=(double *)malloc(sizeof(double)*DM*DM);

    matx(jacm,deriv,elpt,DM,E8,DM);

    *djacb=(*(jacm+0))*(*(jacm+3))-(*(jacm+1))*(*(jacm+2));


    if(*djacb<0.0){
        printf("Determinant < zero\n");
        exit(1);
    }
    
    (*(jaci+0))=(*(jacm+3))/(*djacb);
    (*(jaci+3))=(*(jacm+0))/(*djacb);
    (*(jaci+1))=-(*(jacm+1))/(*djacb);
    (*(jaci+2))=-(*(jacm+2))/(*djacb);

    matx(cart,jaci,deriv,DM,DM,E8);
    
    free(jacm);
    free(jaci);
}

void make_b(double *b, double *cart)
{
    double *b1,*b2,*cart1;
    int i;

    for(i=0,b1=b+E16,b2=b1+E16,cart1=cart+E8;i<E8;i++){
        *b++=*cart;
        *b++=0.0;
        *b1++=0.0;
        *b1++=*cart1;
        *b2++=*cart1++;
        *b2++=*cart++;
    }
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

void make_b_displacementgradient(double *bdispgrad, double *cart)
{
  double *b1,*b2,*b3,*cart1;
  int i;

  for(i=0,b1=bdispgrad+E16,b2=b1+E16,b3=b2+E16,cart1=cart+E8;i<E8;i++){
    *bdispgrad++=*cart;
    *bdispgrad++=0.0;
    *b1++=0.0;
    *b1++=*cart1;
    *b2++=*cart1++;
    *b2++=0.0;
    *b3++=0.0;
    *b3++=*cart++;
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
    free(elpt);
    free(disp);
    free(bmat);
    free(eps);
}

double through_thickness_strain(double *time, double *timez,double * tau, double *gmodulus, double *kmodulus, double *ex, double *ey, double *ez, double ge, double ke, int *number, int num_of_term)
{
    double *Re11, *Im11, *Re22, *Im22, *omega, coefgamma = 5.0, gamma, seikou = 4545, *AnsRe33, *AnsIm33, *ResGs, *ImsGs, *ResKs, *ImsKs;
    int n = 0;
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
    gamma = coefgamma / ((*(time + (*(number + 0)) - 1)) - deltat);
    double deltaomega = 2.0* PI / (*(number + 0)) / deltat;
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

    for (n = 0, (*(Reez + n)) = 0, (*(Imez + n)) = 0; n < (*(number + 0)) / 2 + 1; n++) {//εz(s)の計算
                                                                                         //計算１
        (*(Reez + n)) = -((3 * (*(ResKs + n)) + 4 * (*(ResGs + n)))*(3 * ((*(ResKs + n))*((*(Re11 + n)) + (*(Re22 + n))) - (*(ImsKs + n))*((*(Im11 + n)) + (*(Im22 + n)))) - 2 * ((*(ResGs + n))*((*(Re11 + n)) + (*(Re22 + n))) - (*(ImsGs + n))*((*(Im11 + n)) + (*(Im22 + n))))) + (3 * (*(ImsKs + n)) + 4 * (*(ImsGs + n)))*(3 * ((*(ResKs + n))*((*(Im11 + n)) + (*(Im22 + n))) + (*(ImsKs + n))*((*(Re11 + n)) + (*(Re22 + n)))) - 2 * ((*(ResGs + n))*((*(Im11 + n)) + (*(Im22 + n))) + (*(ImsGs + n))*((*(Re11 + n)) + (*(Re22 + n)))))) / ((3 * (*(ResKs + n)) + 4 * (*(ResGs + n)))*(3 * (*(ResKs + n)) + 4 * (*(ResGs + n))) + (3 * (*(ImsKs + n)) + 4 * (*(ImsGs + n)))*(3 * (*(ImsKs + n)) + 4 * (*(ImsGs + n))));

        (*(Imez + n)) = -((3 * (*(ResKs + n)) + 4 * (*(ResGs + n)))*(3 * ((*(ResKs + n))*((*(Im11 + n)) + (*(Im22 + n))) + (*(ImsKs + n))*((*(Re11 + n)) + (*(Re22 + n)))) - 2 * ((*(ResGs + n))*((*(Im11 + n)) + (*(Im22 + n))) + (*(ImsGs + n))*((*(Re11 + n)) + (*(Re22 + n))))) - (3 * (*(ImsKs + n)) + 4 * (*(ImsGs + n)))*(3 * ((*(ResKs + n))*((*(Re11 + n)) + (*(Re22 + n))) - (*(ImsKs + n))*((*(Im11 + n)) + (*(Im22 + n)))) - 2 * ((*(ResGs + n))*((*(Re11 + n)) + (*(Re22 + n))) - (*(ImsGs + n))*((*(Im11 + n)) + (*(Im22 + n)))))) / ((3 * (*(ResKs + n)) + 4 * (*(ResGs + n)))*(3 * (*(ResKs + n)) + 4 * (*(ResGs + n))) + (3 * (*(ImsKs + n)) + 4 * (*(ImsGs + n)))*(3 * (*(ImsKs + n)) + 4 * (*(ImsGs + n))));

    }

    for (n = 0, (*(Reez + n)) = 0, (*(Imez + n)) = 0; n < (*(number + 0)) / 2 + 1; n++) {//εz(s)の計算//改定//計算１
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
            (*(AnsRe33 + n)) += ((*(Reez + k)) * c - (*(Imez + k))* s) * exp(gamma * n*deltat) * deltaomega / 2.000000000 / PI;//実部を計算
            pr1 = ((*(Reez + k)) * c - (*(Imez + k)) * s) * exp(gamma * n*deltat) * deltaomega / 2.0000000000000 / PI;
     
            (*(AnsIm33 + n)) += ((*(Imez + k)) * c + (*(Reez + k)) * s) * exp(gamma * n*deltat)  * deltaomega / 2.00000000000 / PI;//虚部を計算
            pr2 = ((*(Imez + k))* c + (*(Reez + k)) * s) * exp(gamma * n*deltat)  * deltaomega / 2.000000000000 / PI;
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
    free(Reez);
    free(Imez);
    free(ResGs);
    free(ImsGs);
    free(ResKs);
    free(ImsKs);
        
    return seikou;
}

void differentiatelapack(double *x, double *y, double *dydx)
{
    int i;
    double a, b, c;
    double x0;

    double *mat_a, *mat_f;
    int nrhs = 1, lwork, info, lda, ldb, m, n;
    double *work;
    char trans = 'N';
    
    mat_a = (double *)malloc(sizeof(double)*M*N);
    mat_f = (double *)malloc(sizeof(double)*M);
    work = (double *)malloc(sizeof(double)*(M + N));
    
    for (i = 2; i < num_of_data - 2; i++) {
        x0 = (*(x + i));
        (*(mat_a + 0)) = ((*(x + i - 2)) - x0)*((*(x + i - 2)) - x0);
        (*(mat_a + 1)) = ((*(x + i - 1)) - x0)*((*(x + i - 1)) - x0);
        (*(mat_a + 2)) = ((*(x + i)) - x0)*((*(x + i)) - 0);
        (*(mat_a + 3)) = ((*(x + i + 1)) - x0)*((*(x + i + 1)) - x0);
        (*(mat_a + 4)) = ((*(x + i + 2)) - x0)*((*(x + i + 2)) - x0);
        (*(mat_a + 5)) = ((*(x + i - 2)) - x0);
        (*(mat_a + 6)) = ((*(x + i - 1)) - x0);
        (*(mat_a + 7)) = ((*(x + i)) - x0);
        (*(mat_a + 8)) = ((*(x + i + 1)) - x0);
        (*(mat_a + 9)) = ((*(x + i + 2)) - x0);
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
        (*(dydx + i)) = b;

        if (i == 2) {
            (*(dydx + 0)) = 2 * a*((*(x + 0)) - x0) + b;
            (*(dydx + 1)) = 2 * a*((*(x + 1)) - x0) + b;
        }
        if (i == num_of_data - 3) {
            (*(dydx + num_of_data - 2)) = 2 * a*((*(x + num_of_data - 2)) - x0) + b;
            (*(dydx + num_of_data - 1)) = 2 * a*((*(x + num_of_data - 1)) - x0) + b;
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
        eps = ((*(y + i - 1)*((t - t1))*(t - t2)) / ((t0 - t1)*(t0 - t2))) + ((*(y + i)*(t - t0)*(t - t2)) / ((t1 - t0)*(t1 - t2))) + ((*(y + i + 1)*(t - t0)*(t - t1)) / ((t2 - t0)*(t2 - t1)));
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
        if (t2 >= t){
                break;
        }
        }

    return s;
}

double fullfieldstrain2stresskai(double *time, double *epsall, double *tau, double *gmodulus, double *kmodulus, double ge, double ke, double *sigmaall, double *zeps)
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
            (*(ex + i)) = (*(ex + i)) - ((*(ekk + i)) / 3);
            (*(ey + i)) = (*(ey + i)) - ((*(ekk + i)) / 3);
            (*(exy + i)) = (*(epsall + 3 * j + nnip*i + 2)) / 2;
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
            (*(sigmaall + 3 * j + nnip*i)) = sx + skk/3;
            (*(sigmaall + 3 * j + nnip*i + 1)) = sy + skk/3;
            (*(sigmaall + 3 * j + nnip*i + 2)) = sxy;
        }
    }
    return 0;
}

void disp2stress(double *time, double *epsall, double *tau, double *gcoef, double *kcoef, double *mater, double *sigmaall, double *zeps){
    
    int i;
    double ge;
    double ke;

    ge = (*(mater + 0));
    ke = (*(mater + 1));
    (*(mater + 0)) = ge;
    (*(mater + 1)) = ke;
    //printf("ge: %lg\t  ", (*(mater)));
    //printf("ke: %lg\n", (*(mater + 1)));
    for (i = 0; i < num_of_term; i++) {
        (*(mater + 2 + 2 * i)) = (*(gcoef + i));
        (*(mater + 2 + 2 * i + 1)) = (*(kcoef + i));
        //printf("%lg\t%lg\n", (*(mater + 2 + 2 * i)), (*(mater + 2 + 2 * i + 1)));
    }
    fullfieldstrain2stresskai(time, epsall, tau, gcoef, kcoef, ge, ke, sigmaall, zeps);
}

int mkmatrix(double *pt, int *el, double *mat_gb)
{
    double *elpt, *eldisp, *bmat, *eps;
    double *shape, *deriv, *cart, *g_po, *g_we;
    double px, py, det;
    int i, l, k, j, ig, jg;
    int ix;

    g_po=(double *)malloc(sizeof(double)*GP);
    g_we=(double *)malloc(sizeof(double)*GP);
    cart=(double *)malloc(sizeof(double)*DM*E8);
    shape=(double *)malloc(sizeof(double)*E8);
    deriv=(double *)malloc(sizeof(double)*DM*E8);
    elpt=(double *)malloc(sizeof(double)*DM*E8);
    eldisp=(double *)malloc(sizeof(double)*DM*E8);
    bmat=(double *)malloc(sizeof(double)*16*3);
    eps=(double *)malloc(sizeof(double)*3);

    make_gp(g_po,g_we);
    for(j = 0; j < NN; j++){
        for(i = 0; i < NN; i++){
            (*(mat_gb + i + j * NN)) = 0;
        }
    }
    for(j = 0, l = 0; l < ne; l++, el += 9){
        make_pt(elpt, el, pt);
        for(i = 0; i < E8; i++){
            k = *(el + i);
            //printf("k: %d %d %d\n", k, *(el + i), i);
        }
    for(ig = 0; ig < GP; ig++){
        for(jg = 0; jg < GP; jg++){
            px = (*(g_po + ig));
            py = (*(g_po + jg));
            make_sp(shape, deriv, px, py);
            jacob(cart, &det, elpt, deriv);
            make_b(bmat, cart);
            for(ix = 0; ix < 8; ix++){
                k = *(el + ix);
                (*(mat_gb + 2 * k + 3 * j * (nn))) = (*(bmat + 2 * ix + 0 * 16));
                (*(mat_gb + 2 * k + 1 + 3 * j * (nn))) = (*(bmat + 2 * ix + 1 + 0 * 16));
                (*(mat_gb + 2 * k + (3 * j + 1) * (nn))) = (*(bmat + 2 * ix + 1 * 16));
                (*(mat_gb + 2 * k + 1 + (3 * j + 1) * (nn))) = (*(bmat + 2 * ix + 1 + 1 * 16));
                (*(mat_gb + 2 * k + (3 * j + 2) * (nn))) = (*(bmat + 2 * ix + 2 * 16));
                (*(mat_gb + 2 * k + 1 + (3 * j + 2) * (nn))) = (*(bmat + 2 * ix + 1 + 2 * 16));
                }
                j++;
            }
        }
    }
    //printf("nnip %d %d\n", nnip, j);

    free(g_po);
    free(g_we);
    free(cart);
    free(shape);
    free(deriv);
    free(elpt);
    free(eldisp);
    free(bmat);
    free(eps);
    
    return 0;
}

int matrixsolverxy(double *pt, int *ibn, int *ib, double *sigma, int *ibnc, int *ibc, double *uuord, double *mat_m)
{
    
    long l, i, j, k, kk, a, ak, b;
    double *tktmp;
    int *ibtmp;
    int *ind;
    double *u;
    
    long *ipiv;
    int info, nrhs = 1;
    double *kuuvec;

    char trans = 'N';
    int m, n, lda, ldb, lwork;
    double *work;

    int *ibtmpc, *indc;
    double *maa, *mab, *mac, *mad, *da;
    double *mba, *mbb, *mbc, *db, *ub;
    int kkb;
    double *macu, *damacu;
    int num_of_cvdisp;

    tktmp = (double *)malloc(sizeof(double)*NE*(E8 + 1) * 3 * NN);

    ibtmp = (int *)malloc(sizeof(int)*NE*(E8 + 1) * 3);
    ind = (int *)malloc(sizeof(int)*NE*(E8 + 1) * 3);
    ipiv = (long *)malloc(sizeof(long)*NE*(E8 + 1) * 3);
    kuuvec = (double *)malloc(sizeof(double)*NE*(E8 + 1) * 3 * NN);

    u = (double *)malloc(sizeof(double)*NE*(E8 + 1) * 3);

    ibtmpc = (int *)malloc(sizeof(int)*NE*(E8 + 1) * 3);
    indc = (int *)malloc(sizeof(int)*NE*(E8 + 1) * 3);
    maa = (double *)malloc(sizeof(double)*NE*(E8 + 1) * 3 * NN);
    mab = (double *)malloc(sizeof(double)*NE*(E8 + 1) * 3 * NN);
    mac = (double *)malloc(sizeof(double)*NE*(E8 + 1) * 3 * NN);
    mad = (double *)malloc(sizeof(double)*NE*(E8 + 1) * 3 * NN);
    da = (double *)malloc(sizeof(double)*NE*(E8 + 1) * 3);
    mba = (double *)malloc(sizeof(double)*NE*(E8 + 1) * 3 * NN);
    mbb = (double *)malloc(sizeof(double)*NE*(E8 + 1) * 3 * NN);
    mbc = (double *)malloc(sizeof(double)*NE*(E8 + 1) * 3 * NN);
    db = (double *)malloc(sizeof(double)*NE*(E8 + 1) * 3);
    ub = (double *)malloc(sizeof(double)*NE*(E8 + 1) * 3);
    macu = (double *)malloc(sizeof(double)*NE*(E8 + 1) * 3);
    damacu = (double *)malloc(sizeof(double)*NE*(E8 + 1) * 3);
    work = (double *)malloc(sizeof(double)*NE*(E8 + 1) * 3 * NN);
    
    for (i = 0; i < nnip; i++) {
        (*(ibtmp + i)) = 0;
        (*(ind + i)) = 0;
        (*(ibtmpc + i)) = 0;
        (*(indc + i)) = 0;
    }


    for (j = 0; j < nnip; j++) {
        for (i = 0; i < nn; i++) {
            a = nn*j + i;
            *(tktmp + a) = 0;
            *(maa + a) = 0;
            *(mab + a) = 0;
            *(mac + a) = 0;
            *(mad + a) = 0;
        }
    }

    for (j = 0, i = 0; i < nb; i++, j++) {
        k = 2 * (*(ibn + i));
        *(ibtmp + k) = (*(ib + 2 * i));
        *(ibtmp + k + 1) = (*(ib + 2 * i + 1));
    }

    num_of_cvdisp = 0;
    for (j = 0, i = 0; i < nbc; i++, j++) {
        kk = 2 * (*(ibnc + i));
        *(ibtmpc + kk) = (*(ibc + 2 * i));
        *(ibtmpc + kk + 1) = (*(ibc + 2 * i + 1));
        if (num_of_cvdisp<(*(ibc + 2 * i + 1)) || num_of_cvdisp<(*(ibc + 2 * i)))
            num_of_cvdisp += 1;
    }
    //printf("num_of_const_vdisp %d\n", num_of_cvdisp);
  
    k = 0;
    kk = 0;
    for (i = 0; i < nnip; i++) {
        if (*(ibtmp + i) == 1) {
            (*(ind + i)) = -1;
        }
        if (*(ibtmpc + i) != 0) {
            (*(ind + i)) = -1 - (*(ibtmpc + i));
        }
    }
    for (i = 0, j = 0; i < nnip; i++) {
        if ((*(ind + i)) == 0) {
            (*(ind + i)) = (int)j;
            j++;
        }
    }

    for (kk = 0, i = 0; i < nn; i++) {
        if (*(ind + i) >= 0) {
            for (j = 0; j < nnip; j++) {
                a = nn * j + i;
                ak = nn * j + kk;
                (*(maa + ak)) = (*(mat_m + a));
            }
            kk++;
        }
    }

    for (l = 0; l < num_of_cvdisp; l++) {
        for (kkb = 0, i = 0; i < nn; i++) {
            if ((*(ind + i)) == -l - 2) {
                for (j = 0; j < nnip; j++) {
                    a = nn * j + i;
                    (*(maa + j * nn + kk + l)) += (*(mat_m + a));
                }
            }
        }
    }
    //printf("kk+num_of_cvdisp %d\n", (int)kk + num_of_cvdisp);
    for (j = 0; j < nnip; j++) {
        for (i = 0; i < kk + num_of_cvdisp; i++) {
            a = nn * j + i;
            b = nnip * i + j;
            (*(kuuvec + b)) = (*(maa + a));
        }
    }
    
    m = nnip;
    n = (int)kk + num_of_cvdisp;
    lda = m;
    ldb = m;
    lwork = n + m;
    
    //printf("%d %d\n", (int)n, (int)m);
    
    dgels_(&trans, &m, &n, &nrhs, kuuvec, &lda, sigma, &ldb, work, &lwork, &info);
    
    //printf("info %d\n", (int)info);

    for (l = 0; l < num_of_cvdisp; l++) {
        for (j = 0; j < nn; j++) {
            if ((*(ind + j)) == -l - 2) {
                (*(uuord + j)) = (*(sigma + kk + l));
                //printf("uuord %ld %ld %lf\n", j, kk + l, (*(uuord + j)));
            }
        }
    }

    for (k = 0, j = 0; j < nn; j++) {
        if (*(ind + j) >= 0) {
            *(uuord + j) = (*(sigma + k));
            k++;
        }
        else
            if ((*(ind + j)) == -1)
                *(uuord + j) = 0;
    }
    
    free(tktmp);
    free(ibtmp);
    free(ind);
    free(ipiv);
    free(kuuvec);
    free(u);
    free(ibtmpc);
    free(indc);
    free(maa);
    free(mab);
    free(mac);
    free(mad);
    free(da);
    free(mba);
    free(mbb);
    free(mbc);
    free(db);
    free(ub);
    free(macu);
    free(damacu);
    free(work);
    
    return 0;
}

int main(int argc, char **argv){
    
    FILE *fp1;
    FILE *fp2;
    FILE *fp3;
    FILE *fp4;
    FILE *fp5;
    FILE *fp6;
    FILE *fp7;
    FILE *fp8;
    
    char filen[FL];
    char file1[FL];
    char file2[FL];
    char file3[FL];
    char file4[FL];
    char file5[FL];
    char file6[FL];
    char newfile2[FL];
    char newfile5[FL];
    char newfile6[FL];
    char uuordfilex[FL];
    char uuordfiley[FL];
    char senstivitystressxx[FL];
    char senstivitystressyy[FL];
    char senstivitystressxy[FL];
    
    double *pt;
    int *el;
    double *pf;
    double *mat;
    int *ib;
    int *ibn;
    double *disp;
    int *ibc;
    int *ibnc;
    double *points;
    
    double *tau1;
    double *gcoef1;
    double *kcoef1;
    double *mater1;
    
    double *tau2;
    double *gcoef2;
    double *kcoef2;
    double *mater2;
    
    double t;
    double *time;
    double *uuord;
    double *eps;
    double *dispall;
    double *epsall;
    
    double *x;
    double *y;
    
    double *sigmaall1;
    double *zeps1;
    double *sigmaxall1;
    double *sigmayall1;
    double *gammaxyall1;
    
    double *sigmaall2;
    double *zeps2;
    double *sigmaxall2;
    double *sigmayall2;
    double *gammaxyall2;
    
    double *sigmaall3;
    double *zeps3;
    double *sigmaxall3;
    double *sigmayall3;
    double *gammaxyall3;
    
    double *sigmaall4;
    double *zeps4;
    double *sigmaxall4;
    double *sigmayall4;
    double *gammaxyall4;
    
    int modulusnumber;
    int t1;
    int t2;
    double dt;
    
    double *dsigmax;
    double *dsigmay;
    double *dgammaxy;
    
    double *mat_gb;
    double *vuuord;
    double *dsigmall;
    
    int i;
    int j;
    
    pt=(double *)malloc(sizeof(double)*NN);
    pf=(double *)malloc(sizeof(double)*NN);
    mat=(double *)malloc(sizeof(double)*NM*3);
    el=(int *)malloc(sizeof(int)*NE*(E8+1));
    ib=(int *)malloc(sizeof(int)*NB*3);
    ibn=(int *)malloc(sizeof(int)*NB);
    ibc=(int *)malloc(sizeof(int)*NB*3);
    ibnc=(int *)malloc(sizeof(int)*NB);
    disp=(double *)malloc(sizeof(double)*NN);
    points=(double *)malloc(sizeof(double)*NE*(E8+1)*2*NUM_OF_DATA1);
    
    mater1 = (double *)malloc(sizeof(double)*(NUM_OF_TERM * 3 + 2));
    tau1 = (double *)malloc(sizeof(double)*NUM_OF_TERM);
    gcoef1 = (double *)malloc(sizeof(double)*NUM_OF_TERM);
    kcoef1 = (double *)malloc(sizeof(double)*NUM_OF_TERM);
    
    mater2 = (double *)malloc(sizeof(double)*(NUM_OF_TERM * 3 + 2));
    tau2 = (double *)malloc(sizeof(double)*NUM_OF_TERM);
    gcoef2 = (double *)malloc(sizeof(double)*NUM_OF_TERM);
    kcoef2 = (double *)malloc(sizeof(double)*NUM_OF_TERM);
    
    time = (double *)malloc(sizeof(double)*NUM_OF_DATA);
    uuord=(double *)malloc(sizeof(double)*NN);
    eps = (double *)malloc(sizeof(double)*NE*(E8 + 1) * 3 * 10);
    dispall = (double *)malloc(sizeof(double)*NN*NUM_OF_DATA);
    epsall = (double *)malloc(sizeof(double)*NE*(E8 + 1) * 3 * NUM_OF_DATA * 10);
    
    x = (double *)malloc(sizeof(double)*NE*(E8 + 1) * 3 * NUM_OF_DATA * 10);
    y = (double *)malloc(sizeof(double)*NE*(E8 + 1) * 3 * NUM_OF_DATA * 10);
    
    sigmaall1 = (double *)malloc(sizeof(double)*NE*(E8 + 1) * 3 * NUM_OF_DATA * 10);
    zeps1 = (double *)malloc(sizeof(double)*NE*(E8 + 1) * 3 * NUM_OF_DATA * 10);
    sigmaxall1 = (double *)malloc(sizeof(double)*NE*(E8 + 1) * 3 * NUM_OF_DATA * 10);
    sigmayall1 = (double *)malloc(sizeof(double)*NE*(E8 + 1) * 3 * NUM_OF_DATA * 10);
    gammaxyall1 = (double *)malloc(sizeof(double)*NE*(E8 + 1) * 3 * NUM_OF_DATA * 10);

    sigmaall2 = (double *)malloc(sizeof(double)*NE*(E8 + 1) * 3 * NUM_OF_DATA * 10);
    zeps2 = (double *)malloc(sizeof(double)*NE*(E8 + 1) * 3 * NUM_OF_DATA * 10);
    sigmaxall2 = (double *)malloc(sizeof(double)*NE*(E8 + 1) * 3 * NUM_OF_DATA * 10);
    sigmayall2 = (double *)malloc(sizeof(double)*NE*(E8 + 1) * 3 * NUM_OF_DATA * 10);
    gammaxyall2 = (double *)malloc(sizeof(double)*NE*(E8 + 1) * 3 * NUM_OF_DATA * 10);
    
    sigmaall3 = (double *)malloc(sizeof(double)*NE*(E8 + 1) * 3 * NUM_OF_DATA * 10);
    zeps3 = (double *)malloc(sizeof(double)*NE*(E8 + 1) * 3 * NUM_OF_DATA * 10);
    sigmaxall3 = (double *)malloc(sizeof(double)*NE*(E8 + 1) * 3 * NUM_OF_DATA * 10);
    sigmayall3 = (double *)malloc(sizeof(double)*NE*(E8 + 1) * 3 * NUM_OF_DATA * 10);
    gammaxyall3 = (double *)malloc(sizeof(double)*NE*(E8 + 1) * 3 * NUM_OF_DATA * 10);
    
    sigmaall4 = (double *)malloc(sizeof(double)*NE*(E8 + 1) * 3 * NUM_OF_DATA * 10);
    zeps4 = (double *)malloc(sizeof(double)*NE*(E8 + 1) * 3 * NUM_OF_DATA * 10);
    sigmaxall4 = (double *)malloc(sizeof(double)*NE*(E8 + 1) * 3 * NUM_OF_DATA * 10);
    sigmayall4 = (double *)malloc(sizeof(double)*NE*(E8 + 1) * 3 * NUM_OF_DATA * 10);
    gammaxyall4 = (double *)malloc(sizeof(double)*NE*(E8 + 1) * 3 * NUM_OF_DATA * 10);
    
    dsigmax = (double *)malloc(sizeof(double)*NE*(E8 + 1) * 3 * NUM_OF_DATA * 10);
    dsigmay = (double *)malloc(sizeof(double)*NE*(E8 + 1) * 3 * NUM_OF_DATA * 10);
    dgammaxy = (double *)malloc(sizeof(double)*NE*(E8 + 1) * 3 * NUM_OF_DATA * 10);
    
    mat_gb = (double *)malloc(sizeof(double)*NN*NE*(E8 + 1) * 3);
    vuuord = (double *)malloc(sizeof(double)*NN);
    dsigmall = (double *)malloc(sizeof(double)*NE*(E8+1)*3*NUM_OF_DATA1);
    
    if (argc < 2) {
        printf("Sensitivity-Based Virtual Fields\nInput-File-in: ");
        scanf("%s", filen);
    }
    else {
        sprintf(filen, "%s", argv[1]);
    }
    //入力ファイルを読み込む
    if (NULL == (fp1 = fopen(filen, "r"))) {
        printf(" Cannot open file\n");
        exit(0);
    }
    //モデルファイル
    fscanf(fp1, "%s", file1);
    printf("%s", file1);
    readmodel(file1, pt, pf, mat, el, ib, ibn, disp, ibc, ibnc);
    //　材料特性係数ファイル
    fscanf(fp1, "%s", file2);
    //　x＋dx材料特性ファイル名
    printf("EnterModulusNUM: ");
    scanf("%d", &modulusnumber);
    sprintf(newfile2,"%s%d.txt",file2,modulusnumber);
    printf("%s\n", newfile2);
    readmodulus(newfile2, tau2, gcoef2, kcoef2, mater2);
    fscanf(fp1, "%s", file3);
    printf("%s\n", file3);
    readmodulus(file3, tau1, gcoef1, kcoef1, mater1);
    //　温度
    fscanf(fp1, "%lf", &temp);
    printf("\n");
    printf("T = %lg K\n", temp);
    //　面内変位ファイルの数
    fscanf(fp1, "%d", &num_of_data);
    printf("num_of_data %d", num_of_data);
    //　時刻
    printf("\n");
    printf("Enter the t1:");
    scanf("%d", &t1);
    printf("Enter the t2:");
    scanf("%d", &t2);
    dt = t2 - t1;
    printf("time1 = %d\t  time2 = %d\t  dt = %lg \n", t1, t2, dt);
    //　応力変化値ファイル名と仮想変位場ファイル名
    fscanf(fp1, "%s", file5);
    fscanf(fp1, "%s", file6);
    sprintf(newfile5,"%d%s%d%d.txt",modulusnumber,file5,t1,t2);
    sprintf(newfile6,"%d%s%d%d.txt",modulusnumber,file6,t1,t2);
    printf("Stress sensitivity filename : %s \n", newfile5);
    printf("Virtual fields filename : %s \n", newfile6);
    sprintf(senstivitystressxx,"%d%s%d%d-xx.csv",modulusnumber, file5, t1, t2);
    sprintf(senstivitystressyy,"%d%s%d%d-yy.csv",modulusnumber, file5, t1, t2);
    sprintf(senstivitystressxy,"%d%s%d%d-xy.csv",modulusnumber, file5, t1, t2);
    printf("CsvFile for GraphR : %s \n \t\t\t\t\t %s \n \t\t\t\t\t %s \n", senstivitystressxx, senstivitystressyy, senstivitystressxy);
    //　x方向とy方向仮想変位の.csvファイル名
    sprintf(uuordfilex,"%d%s%d%d-x.csv",modulusnumber, file6, t1, t2);
    sprintf(uuordfiley,"%d%s%d%d-y.csv",modulusnumber, file6, t1, t2);
    printf("CsvFile for GraphR : %s \n \t\t\t\t\t %s \n", uuordfilex, uuordfiley);
    for (i = 0; i < num_of_data; i++) {
        fscanf(fp1, "%lf %s", &t, file4);//面内変位を読み込む
        (*(time + i)) = t;
        read_dispdata(file4, uuord);
        strain(uuord, el, pt, points, eps);
        for (j = 0; j < nn; j++) {
            (*(dispall + j + i*nn)) = (*(uuord + j));
        }
        for (j = 0; j < nnip; j++) {
            (*(epsall + j + i*nnip)) = (*(eps + j));
        }
    }
    fclose(fp1);
    disp2stress(time, epsall, tau1, gcoef1, kcoef1, mater1, sigmaall1, zeps1);
    disp2stress(time, epsall, tau2, gcoef2, kcoef2, mater2, sigmaall2, zeps2);
    for (i = 0; i < num_of_data; i++) {
        if (NULL == (fp2 = fopen(newfile5, "w"))) {
            printf("Cannot open dsigma file..\n\n");
            exit(1);
        }
        if( (fp6 = fopen(senstivitystressxx, "w+" )) == NULL)
        {
            printf("Can not open senstivitystressxx \n");
            exit(1);
        }
        if( (fp7 = fopen(senstivitystressyy, "w+" )) == NULL)
        {
            printf("Can not open senstivitystressyy \n");
            exit(1);
        }
        if( (fp8 = fopen(senstivitystressxy, "w+" )) == NULL)
        {
            printf("Can not open senstivitystressxy \n");
            exit(1);
        }
        fprintf(fp6, "%s\t%d\n%s\n%s\t%s\t%s\n", "DataFormat", 2, "memo", "X", "Y", "Data");
        fprintf(fp7, "%s\t%d\n%s\n%s\t%s\t%s\n", "DataFormat", 2, "memo", "X", "Y", "Data");
        fprintf(fp8, "%s\t%d\n%s\n%s\t%s\t%s\n", "DataFormat", 2, "memo", "X", "Y", "Data");
        
        for (j = 0; j < nip; j++) {
            (*(sigmaxall1 + j)) = (*(sigmaall1 + 3 * j + nnip*t1));
            (*(sigmayall1 + j)) = (*(sigmaall1 + 3 * j + nnip*t1 + 1));
            (*(gammaxyall1 + j)) = (*(sigmaall1 + 3 * j + nnip*t1 + 2));
            
            (*(sigmaxall2 + j)) = (*(sigmaall2 + 3 * j + nnip*t1));
            (*(sigmayall2 + j)) = (*(sigmaall2 + 3 * j + nnip*t1 + 1));
            (*(gammaxyall2 + j)) = (*(sigmaall2 + 3 * j + nnip*t1 + 2));
            
            
            (*(sigmaxall3 + j)) = (*(sigmaall1 + 3 * j + nnip*t2));
            (*(sigmayall3 + j)) = (*(sigmaall1 + 3 * j + nnip*t2 + 1));
            (*(gammaxyall3 + j)) = (*(sigmaall1 + 3 * j + nnip*t2 + 2));
            
            (*(sigmaxall4 + j)) = (*(sigmaall2 + 3 * j + nnip*t2));
            (*(sigmayall4 + j)) = (*(sigmaall2 + 3 * j + nnip*t2 + 1));
            (*(gammaxyall4 + j)) = (*(sigmaall2 + 3 * j + nnip*t2 + 2));

            (*(dsigmax + j)) = ( ((*(sigmaxall3 + j)) - (*(sigmaxall4 + j))) - ((*(sigmaxall1 + j)) - (*(sigmaxall2 + j))) );
            (*(dsigmay + j)) = ( ((*(sigmayall3 + j)) - (*(sigmayall4 + j))) - ((*(sigmayall1 + j)) - (*(sigmayall2 + j))) );
            (*(dgammaxy + j)) = ( ((*(gammaxyall3 + j)) - (*(gammaxyall4 + j))) - ((*(gammaxyall1 + j)) - (*(gammaxyall2 + j))) );
            
            (*(dsigmall + 3 * j )) = (*(dsigmax + j));
            (*(dsigmall + 3 * j + 1)) = (*(dsigmay + j));
            (*(dsigmall + 3 * j + 2)) = (*(dgammaxy + j));
            
            (*(x + j)) = (*(points + 2 * j));
            (*(y + j)) = (*(points + 2 * j + 1));
            
            fprintf(fp2, "%lg\t%lg\t%lg\t%lg\t%lg\n", (*(x + j)), (*(y + j)), (*(dsigmall + 3 * j )), (*(dsigmall + 3 * j + 1)), (*(dsigmall + 3 * j + 2)));
            fprintf(fp6, "%lg\t%lg\t%lg\n", (*(x + j)), (*(y + j)), (*(dsigmall + 3 * j )));
            fprintf(fp7, "%lg\t%lg\t%lg\n", (*(x + j)), (*(y + j)), (*(dsigmall + 3 * j + 1)));
            fprintf(fp8, "%lg\t%lg\t%lg\n", (*(x + j)), (*(y + j)), (*(dsigmall + 3 * j + 2)));
        }
        fclose(fp2);
        fclose(fp6);
        fclose(fp7);
        fclose(fp8);
        //inprintf("\n\n\n");
        mkmatrix(pt, el, mat_gb);
        matrixsolverxy(pt, ibn, ib, dsigmall, ibnc, ibc, vuuord, mat_gb);
        if( (fp3 = fopen(newfile6, "w+" )) == NULL)
        {
            printf("Can not open virtual-fields-file \n");
            exit(1);
        }
        if( (fp4 = fopen(uuordfilex, "w+" )) == NULL)
        {
            printf("Can not open uuordfilex \n");
            exit(1);
        }
        if( (fp5 = fopen(uuordfiley, "w+" )) == NULL)
        {
            printf("Can not open uuordfiley \n");
            exit(1);
        }
        fprintf(fp4, "%s\t%d\n%s\n%s\t%s\t%s\n", "DataFormat", 2, "memo", "X", "Y", "Data");
        fprintf(fp5, "%s\t%d\n%s\n%s\t%s\t%s\n", "DataFormat", 2, "memo", "X", "Y", "Data");
        for (i = 0; i < np; i++) {
            fprintf(fp3, "%lg\t%lg\t%lg\t%lg\n", (*(pt + 2 * i)), (*(pt + 2 * i + 1)), (*(vuuord + 2 * i)), (*(vuuord + 2 * i + 1)));
            fprintf(fp4, "%lg\t%lg\t%lg\n", (*(pt + 2 * i)), (*(pt + 2 * i + 1)), (*(vuuord + 2 * i)));
            fprintf(fp5, "%lg\t%lg\t%lg\n", (*(pt + 2 * i)), (*(pt + 2 * i + 1)), (*(vuuord + 2 * i + 1)));
        }
        fclose(fp3);
        fclose(fp4);
        fclose(fp5);
    }
    
    free(pt);
    free(pf);
    free(mat);
    free(el);
    free(ib);
    free(ibn);
    free(ibc);
    free(ibnc);
    free(disp);
    free(points);
    
    free(tau1);
    free(kcoef1);
    free(gcoef1);
    free(mater1);
    
    free(tau2);
    free(kcoef2);
    free(gcoef2);
    free(mater2);
    
    free(time);
    free(uuord);
    free(eps);
    free(dispall);
    free(epsall);
    
    free(x);
    free(y);
    
    free(sigmaall1);
    free(zeps1);
    free(sigmaxall1);
    free(sigmayall1);
    free(gammaxyall1);
    
    free(sigmaall2);
    free(zeps2);
    free(sigmaxall2);
    free(sigmayall2);
    free(gammaxyall2);
    
    free(sigmaall3);
    free(zeps3);
    free(sigmaxall3);
    free(sigmayall3);
    free(gammaxyall3);
    
    free(sigmaall4);
    free(zeps4);
    free(sigmaxall4);
    free(sigmayall4);
    free(gammaxyall4);
    
    free(dsigmax);
    free(dsigmay);
    free(dgammaxy);
    
    free(mat_gb);
    free(vuuord);
    free(dsigmall);
    
    return 0;
}
