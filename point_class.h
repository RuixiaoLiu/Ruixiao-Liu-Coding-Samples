#include <cmath>

//Point class
// an empty constructor of point
template <class myT>
point<myT>::point()
{
    x = 1.;
    y = 1.;
    value = 1.;
}

// an constructor or point class that only consists of usefu info about x and y value
template <class myT>
point<myT>::point(myT xp, myT yp)
{

    x = xp;
    y = yp;
    value = 1.;
    ;
}

// normal constructor of point, with x, y value and the "value" value obtained
// from the equation to optimize
template <class myT>
point<myT>::point(myT xp, myT yp, myT no_of_equ)
{
    // input the independent x and y values
    x = xp;
    y = yp;

    // for the value, need to compute the function value at the given x and y
    // if statements are set to get the function user specifies
    if (no_of_equ == 1)
        value = f1(x, y);
    else if (no_of_equ == 2)
        value = f2(x, y);
    else if (no_of_equ == 3)
        value = f3(x, y);
}

// the three functions for user to choose to evaluate
// 2*x^-2 + y^-2 + 3*y^2 - 2*x^-1 + 4*x^2
template <class myT>
myT point<myT>::f1(myT x, myT y)
{
    myT val;
    val = pow(x, -2.) * 2. + pow(y, -2.) + pow(y, 2.) * 3 - pow(x, -1) * 2. + 4. * pow(x, 2.);      // function to optimize, need to revise for further usage
    return val;
}

// 0.3*(x^2+y^2) - 0.48*x*y, the classic Matyas function
template <class myT>
myT point<myT>::f2(myT x, myT y)
{
    myT val;
    val = 0.3 * (pow(x, 2.) + pow(y, 2.)) - 0.48 * x * y;
    return val;
}

// sin(xy)
template <class myT>
myT point<myT>::f3(myT x, myT y)
{
    myT val;
    val = sin(x * y);
    return val;
}

// shifting of minimum point by the step residual
// the direction and magnitude of shiting is provided by the point rhs
// this is called by another function from grad class, called step
template <class myT>
point<myT> point<myT>::diff(point<myT> rhs, myT no_of_equ)
{
    return point(x - rhs.x, y - rhs.y, no_of_equ);
}
