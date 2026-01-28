#include "roots.hpp"

#include <cmath>
#include <functional>

/*
 * Common numerical parameters
 */
static const double TOL = 1e-6;
static const int MAX_ITER = 1000;

/*
 * Bisection Method
 */
bool bisection(std::function<double(double)> f,
               double a, double b,
               double *root)
{
    double fa = f(a);
    double fb = f(b);

    // Must bracket a root
    if (fa * fb > 0.0)
        return false;

    for (int i = 0; i < MAX_ITER; ++i)
    {
        double c = 0.5 * (a + b);
        double fc = f(c);

        // Convergence check
        if (std::abs(fc) < TOL || 0.5 * (b - a) < TOL)
        {
            *root = c;
            return true;
        }

        // Shrink bracket
        if (fa * fc < 0.0)
        {
            b = c;
            fb = fc;
        }
        else
        {
            a = c;
            fa = fc;
        }
    }

    return false; // did not converge
}

/*
 * Regula Falsi (False Position) Method
 */
bool regula_falsi(std::function<double(double)> f,
                  double a, double b,
                  double *root)
{
    double fa = f(a);
    double fb = f(b);

    // Must bracket a root
    if (fa * fb > 0.0)
        return false;

    for (int i = 0; i < MAX_ITER; ++i)
    {
        // Linear interpolation
        double c = (a * fb - b * fa) / (fb - fa);
        double fc = f(c);

        // Convergence check
        if (std::abs(fc) < TOL)
        {
            *root = c;
            return true;
        }

        // Update bracket
        if (fa * fc < 0.0)
        {
            b = c;
            fb = fc;
        }
        else
        {
            a = c;
            fa = fc;
        }
    }

    return false; // did not converge
}

/*
 * Newton–Raphson Method
 */
bool newton_raphson(std::function<double(double)> f,
                    std::function<double(double)> g,
                    double a, double b, double c,
                    double *root)
{
    double x = c;

    // Initial guess must lie in interval
    if (x < a || x > b)
        return false;

    for (int i = 0; i < MAX_ITER; ++i)
    {
        double fx = f(x);
        double gx = g(x);

        // Derivative must not be zero
        if (std::abs(gx) < 1e-12)
            return false;

        double x_new = x - fx / gx;

        // Iteration must stay inside bracket
        if (x_new < a || x_new > b)
            return false;

        // Convergence check
        if (std::abs(x_new - x) < TOL || std::abs(f(x_new)) < TOL)
        {
            *root = x_new;
            return true;
        }

        x = x_new;
    }

    return false; // did not converge
}

/*
 * Secant Method
 */
bool secant(std::function<double(double)> f,
            double a, double b, double c,
            double *root)
{
    // Choose second starting point slightly offset from c
    double x0 = c;
    double x1 = c + 1e-3;

    // Ensure initial points lie in interval
    if (x0 < a || x0 > b || x1 < a || x1 > b)
        return false;

    double f0 = f(x0);
    double f1 = f(x1);

    for (int i = 0; i < MAX_ITER; ++i)
    {
        // Prevent division by zero
        if (std::abs(f1 - f0) < 1e-12)
            return false;

        double x2 = x1 - f1 * (x1 - x0) / (f1 - f0);

        // Must stay within bracket
        if (x2 < a || x2 > b)
            return false;

        // Convergence check
        if (std::abs(x2 - x1) < TOL || std::abs(f(x2)) < TOL)
        {
            *root = x2;
            return true;
        }

        // Shift points
        x0 = x1;
        f0 = f1;
        x1 = x2;
        f1 = f(x1);
    }

    return false; // did not converge
}
