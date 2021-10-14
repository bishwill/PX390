#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


int isprime(long k);
void happy_meter(int size);

int main(void)
{
    happy_meter(0);
    happy_meter(1);
    happy_meter(5);
    happy_meter(10);
}

/*
A function to check whether a given number is prime or not
*/
int isprime(long k)
{
    if (k == 1 || k == 2)
    {
        return(1);
    }
    else
    {
        for (long int i = 2; i < k; i++)
        {
            if (k % i == 0)
            {
                return(0);
            }
        }
        return(1);
    }
}

void happy_meter(int size)
{
    printf("%s", "I am happy about this assignment ");
    for (int i = 0; i < size; i++)
    {
        printf("%s", ":-) ");
    }
    printf("\n");
}
