#include <stdio.h>
#include <stdlib.h>
#include <math.h>


int main(void)
{
    long int a, b;
    scanf("%ld %ld", &a, &b);
    printf("%ld\n", labs(a-b));
}