#define import_spp
#include <iraf.h>
#include <math.h>
#include <stdlib.h>

#define abs(x) ((x) < 0 ? -(x) : (x))
#define POW2(x) ((x)*(x))

extern int dgetrf_ (XINT* m, XINT* n, double* a, XINT* lda, XINT* ipiv, XINT* info);
extern int dgetri_ (XINT* n, double* a, XINT* lda, XINT* ipiv, double* work, XINT* lwork, XINT* info);
int
pfregs_ (
    float* x, float* y, float* sigmay,
    XINT* npts, XINT* nterms, XINT* m, XINT* mode,
    float* yfit, float* a0, float* a, float* sigma0, float* sigmaa,
    float* r, float* rmul, float* chisqr, float* ftest,
    float (*fctn)(float*, XINT*, XINT*, XINT*)
)
{
    double xmean[10], sigmax[10];
    double sum = 0., ymean = 0., sigma = 0., chisq;
    float weight[1000];

    double* array = (double*)malloc((*nterms) * (*nterms) * sizeof(double));

    chisq = 0.;
    *rmul = 0.;

    for (int i = 0; i < *npts; i++) {
        yfit[i] = 0.;
    }

    for (int i = 0; i < *nterms; i++) {
        xmean[i] = 0.;
        sigmax[i] = 0.;
        r[i] = 0.;
        a[i] = 0.;
        sigmaa[i] = 0.;
        for (int j =0; j < *nterms; j++) {
            array[i+j*(*nterms)] = 0.;
        }
    }

    for (int i = 0; i < *npts; i++) {

        if (*mode == 0 || y[i] == 0) {
            weight[i] = 1.;
        } else if (*mode < 0) {
            weight[i] = 1. / abs(y[i]);
        } else {
            weight[i] = 1. / POW2(sigmay[i]);
        }

        sum += weight[i];
        ymean += weight[i] * y[i];

        for (int j = 0; j < *nterms; j++) {
            XINT I = i+1, J = j+1;
            xmean[j] += weight[i] * fctn(x,&I,&J,m);
        }
    }

    ymean /= sum;

    for (int i = 0; i < *nterms; i++) {
        xmean[i] /= sum;
    }

    float fnpts = *npts;
    float wmean = sum / fnpts;

    for (int i =0; i < *npts; i++) {
        weight[i] /= wmean;
    }

    for (int i = 0; i < *npts; i++) {
        sigma += weight[i] * (y[i] - ymean) * (y[i] - ymean);
        for (int j = 0; j < *nterms; j++) {
            XINT I = i+1, J = j+1;
            sigmax[j] += weight[i] * POW2(fctn(x,&I,&J,m) - xmean[j]);
            r[j] += weight[i] * (fctn(x,&I,&J,m) - xmean[j]) * (y[i] - ymean);
            for (int k = 0; k <= j; k++) {
                XINT K = k+1;
                array[j+k*(*nterms)] += weight[i] * (fctn(x,&I,&J,m) - xmean[j])
                    * (fctn(x, &I,&K,m) - xmean[k]);
            }
        }
    }

    float free1 = *npts - 1;
    sigma = sqrt(sigma / free1);

    for (int i = 0; i < *nterms; i++) {
        sigmax[i] = sqrt(sigmax[i] / free1);
        r[i] /= free1 * sigmax[i] * sigma;
        for (int j = 0; j <= i; j++) {
            array[i+j*(*nterms)] /= free1 * sigmax[i] * sigmax[j];
            array[j+i*(*nterms)] = array[i+j*(*nterms)];
        }
    }

    XINT info;
    XINT* ipiv = (XINT*)malloc((*nterms)*sizeof(XINT));
    dgetrf_(nterms, nterms, array, nterms, ipiv, &info);
    if (info > 0) {
        free(ipiv);
        free(array);
        *a0 = 0.;
        *sigma0 = 0.;
        *rmul = 0.;
        *chisqr = 0.;
        *ftest = 0.;
        return 0;
    }
    XINT lwork = POW2(*nterms);
    double* work = (double*)malloc(lwork * sizeof(double));
    dgetri_(nterms, array, nterms, ipiv, work, &lwork, &info);
    free(work);
    free(ipiv);


    *a0 = ymean;
    for (int i = 0; i < *nterms; i++) {
        for (int j = 0; j <= i; j++) {
            a[i] += r[j] * array[i+j*(*nterms)];
        }
        a[i] *= sigma / sigmax[i];
        *a0 -= a[i] * xmean[i];
        for (int j = 0; j < *npts; j++) {
            XINT I = i+1, J = j+1;
            yfit[j] += a[i] * fctn(x,&J,&I,m);
        }
    }
    for (int i = 0; i < *npts; i++) {
        yfit[i] += *a0;
        chisq += weight[i] * POW2(y[i] - yfit[i]);
    }
    float freen = *npts - *nterms - 1;
    *chisqr = chisq * wmean / freen;

    float varnce;
    if (*mode == 0) {
        varnce = *chisqr;
    } else {
        varnce = 1. / wmean;
    }

    float freej;
    for (int i = 0; i < *nterms; i++) {
        sigmaa[i] = varnce * array[i+i*(*nterms)] / (free1 * POW2(sigmax[i]));
        if (sigmaa[i] > 0) {
            sigmaa[i] = sqrt(sigmaa[i]);
        } else {
            sigmaa[i] = 0.;
        }
        *rmul += a[i] * r[i] * sigmax[i] / sigma;
    }
    freej = *nterms;

    if (*rmul < 0) {
        *ftest = -99999.;
        *rmul = -99999.;
    } else {
        if (abs(*rmul) < 1.) {
            *ftest = (*rmul / freej) / ((1. - *rmul) / freen);
            *rmul = sqrt(*rmul);
        } else if (abs(*rmul) == 1.) {
            *ftest = -99999;
            *rmul = 1.;
        } else {
            *ftest = -99999;
            *rmul = 99999.;
        }
    }

    *sigma0 = varnce / fnpts;

    for (int i = 0; i < *nterms; i++) {
        for (int j = 0; j < *nterms; j++) {
            *sigma0 += varnce * xmean[i] * xmean[j] * array[i + j*(*nterms)]
                / (free1 * sigmax[i] * sigmax[j]);
        }
    }

    *sigma0 = sqrt(*sigma0);
    free(array);
    return 0;
}