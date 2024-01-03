#include <cmath>


// grad class
// empty grad constructor
template <class myT>
grad<myT>::grad()
{
    dzdx = myT();
    dzdy = myT();
}

// normal grad constructor
// calculation of gradient using central difference method
// the step size is not divided due to later calculation (function matmul in hess class)
// multiplies the same step size, therefore redundant
template <class myT>
grad<myT>::grad(point<myT> tar, myT step_size, myT no_of_equ)
{
    if (no_of_equ == 1)
    {
        dzdx = (tar.f1(tar.x + step_size, tar.y) - tar.f1(tar.x - step_size, tar.y)) / 2.;
        dzdy = (tar.f1(tar.x, tar.y + step_size) - tar.f1(tar.x, tar.y - step_size)) / 2.;
    }
    else if (no_of_equ == 2)
    {
        dzdx = (tar.f2(tar.x + step_size, tar.y) - tar.f2(tar.x - step_size, tar.y)) / 2.;
        dzdy = (tar.f2(tar.x, tar.y + step_size) - tar.f2(tar.x, tar.y - step_size)) / 2.;
    }
    else if (no_of_equ == 3)
    {
        dzdx = (tar.f3(tar.x + step_size, tar.y) - tar.f3(tar.x - step_size, tar.y)) / 2.;
        dzdy = (tar.f3(tar.x, tar.y + step_size) - tar.f3(tar.x, tar.y - step_size)) / 2.;
    }
}

// operator of gen_mul function
template <class myT>
grad<myT> grad<myT>::operator*(myT& rhs)
{
    return gen_mul(rhs);
}

// the gen_mul function returns the dot product of a vector (grad) and a scaler (step size)
template <class myT>
grad<myT> grad<myT>::gen_mul(myT rhs)
{
    grad temp;
    temp.dzdx = dzdx * rhs;
    temp.dzdy = dzdy * rhs;
    return temp;
}

// the step function computes the next stationary point which has lower value of "value"
// using the residual step calculated, it will call the function diff from the point 
// class to complete the shift of point
template <class myT>
point<myT> grad<myT>::step(point<myT> rhs, myT no_of_equ)
{
    // because using the diff function needs the amount of shifting in the form of point,
    // a 'temporary' overload of point constructor is used
    rhs = rhs.diff(point<myT>(dzdx, dzdy), no_of_equ);
    return rhs;
}

// residual function calculate the norm of the grad vector
// this serves as a measurment of absolute magnitude of shifting of stationary point after 
// an iteration
template <class myT>
myT grad<myT>::residual()
{
    return pow(pow(dzdx, 2.) + pow(dzdy, 2.), 0.5);
}
