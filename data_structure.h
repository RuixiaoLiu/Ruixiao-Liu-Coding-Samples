#ifndef DATA_STRUCTURE_H
#define DATA_STRUCTURE_H


/*
In the header file presents the data types that are used in the main programme
There are four data types: point, fir_dir, hess, and iteration_hist, all decorated with template myT

The point class stores the x, y, and value, which are the two independent and the dependent variables.
Since the value is obtained by functions of x and y, its constructor requires the x and y coordiante
of the point of interest, and the number noting the equation. Also, the functions of value with x and y
are member functions of the point class. Lastly, the "diff" function shifts the new minimum point by
the step residual

The grad class mainly stores the first derivative (gradient), [d(value)/dx, d(value)/dy], of function at the local point.
It also stores the step residual of the minimization which will be applied to the point. The "get_mul" function returns the
dot product of the vector and a scalar.

The hess class mainly stores the inverse of hessian, [[d2(value)/dx2, d2(value)/dxdy],[d2(value)/dydx, d2(value/dy2]]^-1, of
the function at local point. The "matmul" function returns the dot product of a vector and a matrix. The "alpha"
function returns the maximum eigenvalue of the hessian. The "determinant" function returns the determinant of the matrix.

The iteration_hist record the iteration steps and final results and relevent informations of each steps, including
the step of iteration, residual and function value, into a doubly linked list.

*/


// the point class, which stores the values of x, y and value at stationary points
template <class myT>
class point
{
public:
    // an empty constructor of point, returns x = 1., y = 1., value = 1.
    point();

    // a 'temporary' constructor of point when only data of interest are the x and y value
    // will be used in the step function in grad class
    point(myT xp, myT yp);

    // normal constructor of point, input x value and y value, and "value" is implicitly computed
    // using the no. of equation provided.
    point(myT xp, myT yp, myT no_of_equ);

    // the three values stored in point class
    myT x, y, value;

    // shifting of minimum point by the step residual, which is the point rhs
    point<myT> diff(point<myT> rhs, myT no_of_equ);

    // the objective functions
    myT f1(myT x, myT y);
    myT f2(myT x, myT y);
    myT f3(myT x, myT y);
};


// class iteration_hist, a doubly linked list that records the points
// during optimization
template <class myT>
class iteration_hist :public point<myT>
{
public:
    // the iteration_hist is a doubly linked list which stores the
    // coordinate of the point data
    point<myT> data;

    // the current iteration time
    myT time;

    // and the size of residual step size
    myT residual;

    //here are the pointers to the linked previous and next data point
    iteration_hist<myT>* next;
    iteration_hist<myT>* prev;

    // the destructor of the iteration_hist
    ~iteration_hist();
};



// the grad class is desisgned to store and process the gradient of the stationary points
// [d(value)/dx, d(value)/dy]
template <class myT>
class grad :public point<myT>
{
public:
    // the empty constructor of grad
    grad();

    // the normal constuctor of grad needs the x and y values of the stationary point and the equation
    // for computation
    grad(point<myT> tar, myT step_size, myT no_of_equ);

    // the two values stored are the first order derivative in both directions
    myT dzdx, dzdy;

    // get_mul returns the dot product of an vector and a scalar, so a gradient and a step size
    // the central equation of gradient descend method
    grad<myT> gen_mul(myT rhs);

    // step function revises the stationary point coordinates and values by the first order derivatives stored
    // in the grad
    point<myT> step(point<myT> rhs, myT no_of_eq);

    // residual function calculates the norm of the grad, useful in deciding the magnitude, and therefore the
    // significance, of shifting stationary point further
    myT residual();

    // "*" is the operator for the gen_mul function
    grad operator*(myT& rhs);
};


// hess class stores the inverse of secondary derivatives of the function, which is also refered as a Hessian
template <class myT>
class hess : public grad<myT>
{
public:
    // the empty constructor of the hess class
    hess();

    // the normal constructor of hess, which computes the inverse of secondary derivative at the stationary point
    hess(point<myT> rhs, myT step_size, myT no_of_equ);

    myT r0c0, r0c1, r1c0, r1c1;

    // the mathmatical multiplcation of a vector and a matric, here a grad and a hess
    // the central equation of of quasi-newton method
    grad<myT> matmul(grad<myT> rhs);

    //the operator of matmul
    grad<myT> operator*(grad<myT>& rhs);

    // alpha is the function that returns the maximum eigenvalue of the hessian hess
    myT alpha();

    // this returns the determinant of the matrix of hess, which decides if the hess (matrix)
    // is ill conditioned
    myT determinant();

};


#include "point_class.h"
#include "grad_class.h"
#include "hess_class.h"
#include "iter_hist_class.h"
#endif