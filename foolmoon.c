#include "foolmoon.h"
#include <math.h>
#include <time.h>

#define DEG (M_PI / 180.0)
#define AU (1.495978707e+11)
#define DAY (24.0 * 60.0 * 60.0)
#define KM (10e+3)
#define sqr(A) ((A) * (A))
#define cot(A) (-tan((A) - M_PI / 2))
#define sign(A) (((A) >= 0) ? (1.0) : (-1.0))

enum Body
{
    Body_mercury,
    Body_venus,
    Body_earth,
    Body_mars,
    Body_jupyter,
    Body_saturn,
    Body_uranus,
    Body_neptune,
    Body_pluto
};

struct Differentiable
{
    double value;
    double derivative;
};

struct Orbit
{
    //https://nssdc.gsfc.nasa.gov/planetary/factsheet/fact_notes.html
    double period;
    double semimajor;
    double eccentricity;
    double inclination;
    double longitude_of_ascending;
    double longitude_of_perihelion;
    double longitude_j2000;
};

struct Orbit orbits[] = {
    {
        //https://nssdc.gsfc.nasa.gov/planetary/factsheet/mercuryfact.html
        .period = 87.969 * DAY,
        .semimajor = 0.38709893 * AU,
        .eccentricity = 0.20563069,
        .inclination = 7.00487 * DEG,
        .longitude_of_ascending = 48.33167  * DEG,
        .longitude_of_perihelion = 77.45645 * DEG,
        .longitude_j2000 = 252.25084 * DEG
    },
    {
        //https://nssdc.gsfc.nasa.gov/planetary/factsheet/venusfact.html
        .period = 224.701 * DAY,
        .semimajor = 0.72333199 * AU,
        .eccentricity = 0.00677323,
        .inclination = 3.39471 * DEG,
        .longitude_of_ascending = 76.68069 * DEG,
        .longitude_of_perihelion = 131.53298 * DEG,
        .longitude_j2000 = 181.97973 * DEG
    },
    {
        //https://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html
        .period = 365.256 * DAY,
        .semimajor = 1.00000011 * AU,
        .eccentricity = 0.01671022,
        .inclination = 0.00005 * DEG,
        .longitude_of_ascending = -11.26064 * DEG,
        .longitude_of_perihelion = 102.94719 * DEG,
        .longitude_j2000 = 100.46435 * DEG
    },
    {
        //https://nssdc.gsfc.nasa.gov/planetary/factsheet/marsfact.html
        .period = 686.980 * DAY,
        .semimajor = 1.52366231 * AU,
        .eccentricity = 0.09341233,
        .inclination = 1.85061 * DEG,
        .longitude_of_ascending = 49.57854 * DEG,
        .longitude_of_perihelion = 336.04084 * DEG,
        .longitude_j2000 = 355.45332 * DEG
    },
    {
        //https://nssdc.gsfc.nasa.gov/planetary/factsheet/jupiterfact.html
        .period = 4332.589 * DAY,
        .semimajor = 5.20336301 * AU,
        .eccentricity = 0.04839266,
        .inclination = 1.30530 * DEG,
        .longitude_of_ascending = 100.55615 * DEG,
        .longitude_of_perihelion = 14.75385 * DEG,
        .longitude_j2000 = 34.40438 * DEG
    },
    {
        //https://nssdc.gsfc.nasa.gov/planetary/factsheet/saturnfact.html
        .period = 10759.22 * DAY,
        .semimajor = 9.53707032 * AU,
        .eccentricity = 0.05415060,
        .inclination = 2.48446 * DEG,
        .longitude_of_ascending = 113.71504 * DEG,
        .longitude_of_perihelion = 92.43194 * DEG,
        .longitude_j2000 = 49.94432 * DEG
    },
    {
        //https://nssdc.gsfc.nasa.gov/planetary/factsheet/uranusfact.html
        .period = 30685.4 * DAY,
        .semimajor = 19.19126393 * AU,
        .eccentricity = 0.04716771,
        .inclination = 0.76986 * DEG,
        .longitude_of_ascending = 74.22988 * DEG,
        .longitude_of_perihelion = 170.96424 * DEG,
        .longitude_j2000 = 313.23218 * DEG
    },
    {
        //https://nssdc.gsfc.nasa.gov/planetary/factsheet/neptunefact.html
        .period = 60189.0 * DAY,
        .semimajor = 30.06896348 * AU,
        .eccentricity = 0.00858587,
        .inclination = 1.76917 * DEG,
        .longitude_of_ascending = 131.72169 * DEG,
        .longitude_of_perihelion = 44.97135 * DEG,
        .longitude_j2000 = 304.88003 * DEG
    },
    {
        //https://nssdc.gsfc.nasa.gov/planetary/factsheet/plutofact.html
        .period = 90560.0 * DAY,
        .semimajor = 39.48168677 * AU,
        .eccentricity = 0.24880766,
        .inclination = 17.14175 * DEG,
        .longitude_of_ascending = 110.30347 * DEG,
        .longitude_of_perihelion = 224.06676 * DEG,
        .longitude_j2000 = 238.92881 * DEG
    }
};

struct Differentiable add(struct Differentiable a, struct Differentiable b)
{
    struct Differentiable c = { a.value + b.value, a.derivative + b.derivative };
    return c;
}

struct Differentiable multiply(double a, struct Differentiable b)
{
    struct Differentiable c = { a * b.value, a * b.derivative };
    return c;
}

struct Differentiable sind(struct Differentiable a)
{
    struct Differentiable c = { sin(a.value), cos(a.value) * a.derivative };
    return c;
}

struct Differentiable cosd(struct Differentiable a)
{
    struct Differentiable c = { cos(a.value), -sin(a.value) * a.derivative };
    return c;
}

struct Differentiable atan2d(struct Differentiable y, struct Differentiable x)
{
    struct Differentiable c;
    c.value = atan2(y.value, x.value);
    c.derivative = (x.value * y.derivative - y.value * x.derivative) / (sqr(y.value) + sqr(x.value));
    return c;
}

void rotate_x(double angle, struct Differentiable space[3])
{
    struct Differentiable new_y = add(multiply(cos(angle), space[1]), multiply(-sin(angle), space[2]));
    struct Differentiable new_z = add(multiply(sin(angle), space[1]), multiply( cos(angle), space[2]));
    space[1] = new_y;
    space[2] = new_z;
}

void rotate_y(double angle, struct Differentiable space[3])
{
    struct Differentiable new_x = add(multiply( cos(angle), space[0]), multiply(sin(angle), space[2]));
    struct Differentiable new_z = add(multiply(-sin(angle), space[0]), multiply(cos(angle), space[2]));
    space[0] = new_x;
    space[2] = new_z;
}

void rotate_z(double angle, struct Differentiable space[3])
{
    struct Differentiable new_x = add(multiply(cos(angle), space[0]), multiply(-sin(angle), space[1]));
    struct Differentiable new_y = add(multiply(sin(angle), space[0]), multiply( cos(angle), space[1]));
    space[0] = new_x;
    space[1] = new_y;
}

void true_anomaly_to_plane(struct Orbit *orbit, double true_anomaly, double plane[2])
{
    double semimajor = orbit->semimajor;
    double semidistance = orbit->semimajor * orbit->eccentricity;
    double semimanor = sqrt(sqr(semimajor) - sqr(semidistance));

    double c = sqr(semidistance) / sqr(semimajor) - 1;
    if (fabs(sin(true_anomaly)) < fabs(cos(true_anomaly)))
    {
        double a = 1 / sqr(semimajor) + sqr(tan(true_anomaly)) / sqr(semimanor);
        double b = 2 * semidistance / sqr(semimajor);    
        plane[0] = (-b + sign(cos(true_anomaly)) * sqrt(sqr(b) - 4 * a * c)) / (2 * a);
        plane[1] = plane[0] * tan(true_anomaly);
    }
    else
    {
        double a = sqr(cot(true_anomaly)) / sqr(semimajor) + 1 / sqr(semimanor);
        double b = 2 * cot(true_anomaly) * semidistance / sqr(semimajor);
        plane[1] = (-b + sign(sin(true_anomaly)) * sqrt(sqr(b) - 4 * a * c)) / (2 * a);
        plane[0] = plane[1] * cot(true_anomaly);
    }
}

double plane_to_mean_anomaly(struct Orbit *orbit, double plane[2])
{
    double semimajor = orbit->semimajor;
    double semidistance = orbit->semimajor * orbit->eccentricity;
    double semimanor = sqrt(sqr(semimajor) - sqr(semidistance));

    double eccentric_anomaly = atan2(semimajor * plane[1] / semimajor, plane[0] + semidistance);
    return eccentric_anomaly - orbit->eccentricity * sin(eccentric_anomaly);
}

void mean_anomaly_to_plane(struct Orbit *orbit, struct Differentiable mean_anomaly, struct Differentiable plane[2])
{
    double semimajor = orbit->semimajor;
    double semidistance = orbit->semimajor * orbit->eccentricity;
    double semimanor = sqrt(sqr(semimajor) - sqr(semidistance));

    struct Differentiable eccentric_anomaly = { 0, 0 };
    for (int i = 0; i < 100; i++)
    {
        struct Differentiable new_eccentric_anomaly = add(mean_anomaly, multiply(orbit->eccentricity, sind(eccentric_anomaly)));
        if (new_eccentric_anomaly.value == eccentric_anomaly.value) break;
        else eccentric_anomaly = new_eccentric_anomaly;
    }

    plane[0] = multiply(semimajor, cosd(eccentric_anomaly)); plane[0].value -= semidistance;
    plane[1] = multiply(semimajor, sind(eccentric_anomaly));
}

double longitude_to_argument(struct Orbit *orbit, double longitude)
{
    double argument = longitude - orbit->longitude_of_ascending;
    argument = atan2(sqrt(1 + sqr(tan(orbit->inclination))) * sin(argument), cos(argument));
    return argument;
}

void plane_to_space(struct Orbit *orbit, struct Differentiable plane[2], struct Differentiable space[3])
{
    space[0] = plane[0];
    space[1] = plane[1];
    space[2].value = 0; space[2].derivative = 0;
    rotate_z(longitude_to_argument(orbit, orbit->longitude_of_perihelion), space);
    rotate_x(orbit->inclination, space);
    rotate_z(orbit->longitude_of_ascending, space);
}

struct Differentiable time_to_mean_anomaly(struct Orbit *orbit, double time)
{
    double true_anomaly_j2000 = longitude_to_argument(orbit, orbit->longitude_j2000) - longitude_to_argument(orbit, orbit->longitude_of_perihelion);
    double plane_j2000[2];
    true_anomaly_to_plane(orbit, true_anomaly_j2000, plane_j2000);
    double mean_anomaly_j2000 = plane_to_mean_anomaly(orbit, plane_j2000);
    double mean_anomaly_2pi = mean_anomaly_j2000 / (2 * M_PI) + time / orbit->period;
    mean_anomaly_2pi -= floor(mean_anomaly_2pi);
    struct Differentiable mean_anomaly = { 2 * M_PI * mean_anomaly_2pi, 2 * M_PI / orbit->period };
    return mean_anomaly;
}

void time_to_space(struct Orbit *orbit, double time, struct Differentiable space[3])
{
    struct Differentiable mean_anomaly = time_to_mean_anomaly(orbit, time);
    struct Differentiable plane[2];
    mean_anomaly_to_plane(orbit, mean_anomaly, plane);
    plane_to_space(orbit, plane, space);
}

double get_time()
{
    struct tm j2000_calender;
    j2000_calender.tm_sec = 0;
    j2000_calender.tm_min = 0;
    j2000_calender.tm_hour = 12;
    j2000_calender.tm_mday = 1;
    j2000_calender.tm_mon = 0;
    j2000_calender.tm_year = 100;
    j2000_calender.tm_isdst = 0;
    time_t j2000 = mktime(&j2000_calender);
    time_t now = time(NULL);
    return now - j2000;
}

bool is_planet_retrograde(enum Body body)
{
    double now = get_time();
    struct Differentiable earth[3];
    time_to_space(orbits + Body_earth, now, earth);
    struct Differentiable planet[3];
    time_to_space(orbits + body, now, planet);
    for (unsigned int i = 0; i < 3; i++) planet[i] = add(planet[i], multiply(-1, earth[i]));
    struct Differentiable angle = atan2d(planet[1], planet[0]);
    return angle.derivative < 0;
}

bool is_mercury_retrograde(void) { return is_planet_retrograde(Body_mercury); }
bool is_venus_retrograde(void) { return is_planet_retrograde(Body_venus); }
bool is_mars_retrograde(void) { return is_planet_retrograde(Body_mars); }
bool is_jupyter_retrograde(void) { return is_planet_retrograde(Body_jupyter); }
bool is_saturn_retrograde(void) { return is_planet_retrograde(Body_saturn); }
bool is_uranus_retrograde(void) { return is_planet_retrograde(Body_uranus); }
bool is_neptune_retrograde(void) { return is_planet_retrograde(Body_neptune); }
bool is_pluto_retrograde(void) { return is_planet_retrograde(Body_pluto); }