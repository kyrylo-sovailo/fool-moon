#include "foolmoon.h"
#include <stdio.h>

int main()
{
    printf("Mercury is %s\n", is_mercury_retrograde() ? "retrograde" : "prograde");
    printf("Venus   is %s\n", is_venus_retrograde()   ? "retrograde" : "prograde");
    printf("Mars    is %s\n", is_mars_retrograde()    ? "retrograde" : "prograde");
    printf("Jupyter is %s\n", is_jupyter_retrograde() ? "retrograde" : "prograde");
    printf("Saturn  is %s\n", is_saturn_retrograde()  ? "retrograde" : "prograde");
    printf("Uranus  is %s\n", is_uranus_retrograde()  ? "retrograde" : "prograde");
    printf("Neptune is %s\n", is_neptune_retrograde() ? "retrograde" : "prograde");
    printf("Pluto   is %s\n", is_pluto_retrograde()   ? "retrograde" : "prograde");
    return 0;
}