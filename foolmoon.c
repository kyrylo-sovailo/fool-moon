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
    double mean_longitude_j2000;
};

struct Orbit orbits[2] = {
    {
        //https://nssdc.gsfc.nasa.gov/planetary/factsheet/mercuryfact.html
        .period = 87.969 * DAY,
        .semimajor = 0.38709893 * AU,
        .eccentricity = 0.20563069,
        .inclination = 7.00487 * DEG,
        .longitude_of_ascending = 48.33167  * DEG,
        .longitude_of_perihelion = 77.45645 * DEG,
        .mean_longitude_j2000 = 252.25084 * DEG
    },
    {
        //https://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html
        .period = 365.256 * DAY,
        .semimajor = 1.00000011 * AU,
        .eccentricity = 0.01671022,
        .inclination = 0.00005 * DEG,
        .longitude_of_ascending = -11.26064 * DEG,
        .longitude_of_perihelion = 102.94719 * DEG,
        .mean_longitude_j2000 = 100.46435 * DEG
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

void mean_longitude_to_plane(struct Orbit *orbit, double mean_longitude, double plane[2])
{
    double semimajor = orbit->semimajor;
    double semidistance = orbit->semimajor * orbit->eccentricity;
    double semimanor = sqrt(sqr(semimajor) - sqr(semidistance));

    double c = sqr(semidistance) / sqr(semimajor) - 1;
    if (fabs(sin(mean_longitude)) < fabs(cos(mean_longitude)))
    {
        double a = 1 / sqr(semimajor) + sqr(tan(mean_longitude)) / sqr(semimanor);
        double b = 2 * semidistance / sqr(semimajor);    
        plane[0] = (-b + sign(cos(mean_longitude)) * sqrt(sqr(b) - 4 * a * c)) / (2 * a);
        plane[1] = plane[0] * tan(mean_longitude);
    }
    else
    {
        double a = sqr(cot(mean_longitude)) / sqr(semimajor) + 1 / sqr(semimanor);
        double b = 2 * cot(mean_longitude) * semidistance / sqr(semimajor);
        plane[1] = (-b + sign(sin(mean_longitude)) * sqrt(sqr(b) - 4 * a * c)) / (2 * a);
        plane[0] = plane[1] * cot(mean_longitude);
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
    while (true)
    {
        struct Differentiable new_eccentric_anomaly = add(mean_anomaly, multiply(orbit->eccentricity, sind(eccentric_anomaly)));
        if (new_eccentric_anomaly.value == eccentric_anomaly.value) break;
        else eccentric_anomaly = new_eccentric_anomaly;
    }

    plane[0] = multiply(semimajor, cosd(eccentric_anomaly)); plane[0].value -= semidistance;
    plane[1] = multiply(semimajor, sind(eccentric_anomaly));
}

void plane_to_space(struct Orbit *orbit, struct Differentiable plane[2], struct Differentiable space[3])
{
    space[0] = plane[0];
    space[1] = plane[1];
    space[2].value = 0; space[2].derivative = 0;
    rotate_z(orbit->longitude_of_perihelion, space);
    rotate_x(orbit->inclination, space);
    rotate_z(orbit->longitude_of_ascending, space);
}

struct Differentiable time_to_mean_anomaly(struct Orbit *orbit, double time)
{
    double j2000_plane[2];
    mean_longitude_to_plane(orbit, orbit->mean_longitude_j2000, j2000_plane);
    double mean_anomaly_2pi = plane_to_mean_anomaly(orbit, j2000_plane) / (2 * M_PI) + time / orbit->period;
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

enum body
{
    body_mercury,
    body_earth
};

/*
Formulae:
mean_anomaly_j2000 = 
*/

bool is_mercury_retrograde(void)
{
    double now = get_time();
    struct Differentiable earth[3];
    time_to_space(orbits + body_earth, now, earth);
    struct Differentiable mercury[3];
    time_to_space(orbits + body_mercury, now, mercury);
    for (unsigned int i = 0; i < 3; i++) mercury[i] = add(mercury[i], multiply(-1, earth[i]));
    struct Differentiable angle = atan2d(mercury[1], mercury[0]);
    return angle.derivative < 0;
}