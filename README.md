# Fool Moon

Do you want to make sure that Mercury is prograde before executing your code?

Or maybe you don't want your demo program to be used in production?

Either way, **Fool Moon** is the library for you!

# Example of usage
```
#include <foolmoon.h>

int main()
{
    if (is_mercury_retrograde())
    {
        printf("Mercury is retrograde\n");
        printf("Critical failure!\n");
        return 1;
    }
    ///Your high quality code here:
    do_some_mercury_dependent_job();
}
```