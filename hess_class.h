#include <cmath>

//hess class
// an empty constructor of hess
template <class myT>
hess<myT>::hess()
{
    r0c0 = myT();
    r0c1 = myT();
    r1c0 = myT();
    r1c1 = myT();
}


// the normal constructor of hess, which computes the inverse of hessian of function
// in order to compute the residual step. see documentation for more information
template <class myT>
hess<myT>::hess(point<myT> rhs, myT step_size, myT no_of_equ)
{
    myT R0C0, R0C1, R1C0, R1C1, factor;

    // the secondary derivatives are computed based on extension of central difference method and finite element analysis
    if (no_of_equ == 1)
    {
        R0C0 = (rhs.f1(rhs.x + step_size, rhs.y) - 2 * rhs.f1(rhs.x, rhs.y) + rhs.f1(rhs.x - step_size, rhs.y)) / pow(0.05, 2.);
        R0C1 = (rhs.f1(rhs.x + step_size, rhs.y + step_size) - rhs.f1(rhs.x + step_size, rhs.y) - rhs.f1(rhs.x, rhs.y + step_size) + 2. * rhs.f1(rhs.x, rhs.y) - rhs.f1(rhs.x, rhs.y - step_size) - rhs.f1(rhs.x - step_size, rhs.y) + rhs.f1(rhs.x - step_size, rhs.y - step_size)) / (2 * step_size * step_size);
        R1C0 = (rhs.f1(rhs.x + step_size, rhs.y + step_size) - rhs.f1(rhs.x + step_size, rhs.y) - rhs.f1(rhs.x, rhs.y + step_size) + 2. * rhs.f1(rhs.x, rhs.y) - rhs.f1(rhs.x, rhs.y - step_size) - rhs.f1(rhs.x - step_size, rhs.y) + rhs.f1(rhs.x - step_size, rhs.y - step_size)) / (2 * step_size * step_size);
        R1C1 = (rhs.f1(rhs.x, rhs.y + step_size) - 2 * rhs.f1(rhs.x, rhs.y) + rhs.f1(rhs.x, rhs.y - step_size)) / pow(0.05, 2.);
    }
    else if (no_of_equ == 2)
    {
        R0C0 = (rhs.f2(rhs.x + step_size, rhs.y) - 2 * rhs.f2(rhs.x, rhs.y) + rhs.f2(rhs.x - step_size, rhs.y)) / pow(0.05, 2.);
        R0C1 = (rhs.f2(rhs.x + step_size, rhs.y + step_size) - rhs.f2(rhs.x + step_size, rhs.y) - rhs.f2(rhs.x, rhs.y + step_size) + 2. * rhs.f2(rhs.x, rhs.y) - rhs.f2(rhs.x, rhs.y - step_size) - rhs.f2(rhs.x - step_size, rhs.y) + rhs.f2(rhs.x - step_size, rhs.y - step_size)) / (2 * step_size * step_size);
        R1C0 = (rhs.f2(rhs.x + step_size, rhs.y + step_size) - rhs.f2(rhs.x + step_size, rhs.y) - rhs.f2(rhs.x, rhs.y + step_size) + 2. * rhs.f2(rhs.x, rhs.y) - rhs.f2(rhs.x, rhs.y - step_size) - rhs.f2(rhs.x - step_size, rhs.y) + rhs.f2(rhs.x - step_size, rhs.y - step_size)) / (2 * step_size * step_size);
        R1C1 = (rhs.f2(rhs.x, rhs.y + step_size) - 2 * rhs.f2(rhs.x, rhs.y) + rhs.f2(rhs.x, rhs.y - step_size)) / pow(0.05, 2.);
    }
    else if (no_of_equ == 3)
    {
        R0C0 = (rhs.f3(rhs.x + step_size, rhs.y) - 2 * rhs.f3(rhs.x, rhs.y) + rhs.f3(rhs.x - step_size, rhs.y)) / pow(0.05, 2.);
        R0C1 = (rhs.f3(rhs.x + step_size, rhs.y + step_size) - rhs.f3(rhs.x + step_size, rhs.y) - rhs.f3(rhs.x, rhs.y + step_size) + 2. * rhs.f3(rhs.x, rhs.y) - rhs.f3(rhs.x, rhs.y - step_size) - rhs.f3(rhs.x - step_size, rhs.y) + rhs.f3(rhs.x - step_size, rhs.y - step_size)) / (2 * step_size * step_size);
        R1C0 = (rhs.f3(rhs.x + step_size, rhs.y + step_size) - rhs.f3(rhs.x + step_size, rhs.y) - rhs.f3(rhs.x, rhs.y + step_size) + 2. * rhs.f3(rhs.x, rhs.y) - rhs.f3(rhs.x, rhs.y - step_size) - rhs.f3(rhs.x - step_size, rhs.y) + rhs.f3(rhs.x - step_size, rhs.y - step_size)) / (2 * step_size * step_size);
        R1C1 = (rhs.f3(rhs.x, rhs.y + step_size) - 2 * rhs.f3(rhs.x, rhs.y) + rhs.f3(rhs.x, rhs.y - step_size)) / pow(0.05, 2.);
    }

    // computing the inverse of hessian
    // this is the determinant of the normal hessian
    factor = 1 / ((R0C0 * R1C1) - (R1C0 * R0C1));

    r0c0 = factor * R1C1;
    r0c1 = -factor * R1C0;
    r1c0 = r0c1;
    r1c1 = factor * R0C0;
}

// determinant function computes the determinant of the inverse of hessian
template <class myT>
myT hess<myT>::determinant()
{
    return 1 / ((r0c0 * r1c1) - (r1c0 * r0c1));
}

// the matmul function computes the residual step of quasi-newton method
// mathematically it conducts a mathematical multiplication
template <class myT>
grad<myT> hess<myT>::matmul(grad<myT> rhs)
{
    grad<myT> temp;
    temp.dzdx = rhs.dzdx * r0c0 + rhs.dzdy * r0c1;
    temp.dzdy = rhs.dzdx * r1c0 + rhs.dzdy * r1c1;
    return temp;
}

// operator for matmul function
template <class myT>
grad<myT> hess<myT>::operator*(grad<myT>& rhs)        // very important
{
    return matmul(rhs);
}

// the alpha function returns the largest value of eigenvalues of the inverse of hess (matrix)
template <class myT>
myT hess<myT>::alpha()
{
    myT sig1, sig2;
    // computation of eigenvalues
    sig1 = 0.5 * (r0c0 + r1c1) + pow(pow(r0c1, 2) + 0.25 * (r0c0 - r1c1), 1 / 2);
    sig2 = 0.5 * (r0c0 + r1c1) - pow(pow(r0c1, 2) + 0.25 * (r0c0 - r1c1), 1 / 2);
    if (sig1 >= sig2)
    {
        return sig1;
    }
    else
    {
        return sig2;
    }
}
