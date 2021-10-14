#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define TRUE 1;
#define FALSE 0;

int Print_Border(char Symbol, int NumSymbols);

int main(void)
{
    int check_flag, msg_length;
    /* Prefix 'const' indicates that this should not be chnaged
    in the program. */
    const char* message = "Welcome to Scientific Computing 2020/21";

    /* Find out how many characters are in the string */
    msg_length = strlen(message);
    /* Print a border */
    check_flag = Print_Border('*', msg_length + 4);
    /* Print the message */
    printf("\n* %s *\n", message);
    /* Print another border */
    check_flag = Print_Border('*', msg_length + 4);
    printf("\n");

    return(0);
}

int Print_Border(char Symbol, int NumSymbols)
{
    int i;
    /* Print 'NumSymbols' times a character 'Symbol' */
    for(i = 0; i < NumSymbols; i++)
    {
     	/* Note the format for a charcter: '%c' */
        printf("%c", Symbol);
    }

    return(i);
}