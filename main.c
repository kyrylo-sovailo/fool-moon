#include "foolmoon.h"
#include <stdio.h>

int main()
{
    if (is_mercury_retrograde()) printf("Mercury is retrograde\n");
    else printf("Mercury is prograde\n");
    return 0;
}