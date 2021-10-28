#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


int isprime(long k);
void happy_meter(int size);
int num_conseq_digits(long k);

int main(void)
{
    return(0);
}

/*
A function to check whether a given number is prime or not
*/
int isprime(long k)
{
    if (k == 1)
    {
        return(0);
    }
    if (k == 2)
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

/*
A function to print smiley faces :-)
*/
void happy_meter(int size)
{
    printf("I am happy about this assignment ");
    for (int i = 0; i < size; i++)
    {
        printf("%s", ":-) ");
    }
    printf("\n");
}


int num_conseq_digits(long k)
{
    int numIter = log10(k) + 1;
    long x = k;
    int lastNum = -1;
    int count = 1;
    int tempCount = 1;
    for (int i = 0; i < numIter; i++)
    {
        int currentNum = x % 10;
        if (currentNum == lastNum)
        {
            tempCount += 1;
        }
        else
        {
            if (tempCount > count)
            {
                count = tempCount;
            }
            tempCount = 1;
            lastNum = currentNum;
        }
        x /= 10;
    }
    if (tempCount > count)
    {
        count = tempCount;
    }
    return(count);
}