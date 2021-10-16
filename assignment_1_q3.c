#include <stdio.h>
#include <stdlib.h>
#include <math.h>


int main(void)
{
    long A, B;
    scanf("%ld %ld", &A, &B);
    printf("%ld\n", (long) powl(A, B));
}