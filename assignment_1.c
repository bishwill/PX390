#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


int isprime(long k);

int main(void)
{
    for (int k = 1; k <= 20; k++)
    {
        printf("%d %d\n", k, isprime(k));
    }
    return(0);
}

/*
A fumction to check whether a given number is prime or not
*/
int isprime(long k)
{
    long int i;
    if (k == 1 || k == 2)
    {
        return(1);
    }
    else
    {
        for (i = 2; i < k; i++)
        {
            if (k % i == 0)
            {
                return(0);
            }
        }
        return(1);
    }
}