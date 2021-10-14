#include <stdio.h>
#include <stdlib.h>
#include <math.h>


int main(void)
{
    long int A, B;
    scanf("%ld, %ld", &A, &B);
    printf("%ld\n", (long int) powl(A, B));
}