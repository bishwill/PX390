#include <stdlib.h>
#include <stdio.h>
#include <math.h>


void printArray(int *arr, int arr_length);
int *reverse_order(int *arr, int arr_length);
int  *cumsum(int *arr, int arr_length);
int arr_dot_product(int *arr1, int *arr2, int len);
void varr_swap(char *arr, int a, int b);
void en_pot(double *posx, double *posy, double *posz,  long ncharges, double *res);


int main(void)
{
   	
    
    double arr1[] = {-1.0, 2.0, -3.0};
    double arr2[] = {4.0, -5.0, 6.0};
    double arr3[] = {-7.0, 8.0, -9.0};


    double output;
    double *res = &output;

    en_pot(arr1, arr2, arr3, 3, res);
    
    printf("%lf\n", *res);
    
   
    
    return(0);
}

int *reverse_order(int *arr, int arr_length) {
    int *newArray = malloc(sizeof(int) * arr_length);
    if (!newArray) {
        return NULL;
    }
    for (int i = arr_length - 1; i > -1; i--) {
        newArray[arr_length - i - 1] = arr[i];
    }
    return newArray;
}

int *cumsum(int *arr, int arr_length) {
    int *newArray = malloc(sizeof(int) * (arr_length + 1));
    newArray[0] = 0;
    for (int i = 0; i < arr_length; i++) {
        newArray[i+1] = newArray[i] + arr[i];
    }
    return newArray;
}

int arr_dot_product(int *arr1, int *arr2, int len) {
    int sum = 0;
    for (int i = 0; i < len; i ++) {
        sum += arr1[i] * arr2[i];
    }
    return sum;
}

void varr_swap(char *arr, int a, int b) {
    char tmp = arr[a];
    arr[a] = arr[b];
    arr[b] = tmp;
}

void en_pot(double *posx, double *posy, double *posz,  long ncharges, double *res) {
    *res = 0;

    for (long i = 0; i < ncharges; i++) {
        for (long j = 0; j < ncharges; j++) {
            if (i != j) {
                double x2 = pow(posx[i] - posx[j], 2);
                double y2 = pow(posy[i] - posy[j], 2);
                double z2 = pow(posz[i] - posz[j], 2);
                double invNorm = pow(x2 + y2 + z2, -0.5);
                *res += invNorm;
            }
        }
    }
    *res /= 2;
}