#include <iostream>
#include <fstream>
#include <string.h>
#include <cmath>
#include "data_structure.h"


using namespace std;
/*
A quick manual

Before using the programme, the user needs to have the input.txt file stored in the same folder as the programme. There are 4
parameters in the file separaeted by spaces. From left to right: step size for derivative computation, tolerence of the determinant
of hessian to be ill-conditioned, tolerence of the residual step size (convergence), and the maximum number of iteration.

The user is asked to select from the three pre-set equations to minimize, and enter the x and y value of the initial guess, the
optimization method, and the number of iteration before each recording of data. After entering these information, the programme will
run and return the result. After recording the data in the result.txt file, it will delete the history recursively

Although there is no repeating of input to avoid invalid input, the programme will immediately terminates ones it detects an invalid
input from the user. A warning of termination of the programme will also show on the result.txt file to record the incidence.

please refer to the documentation for more elaborate descriptions of the code

*/





// all variables used in the main function
string method, y_n, start;
double step_size, norm_tol, converge_tol, start_x, start_y;
double result_x, result_y, result_value, result_residual;
int max_iter, no_of_equ, record_seq, j;


/*
read certain variables form the input.txt file
the step size of computing of derivative, the tolerance of norm of hessian matrix
the tolerance of covergence, and the maximum number of iteration
it returns a string, which is a flag of wherether an input file is read successfully
if no input file detected, then the codes terminates
*/
string read_file(const char* fname)
{
    fstream file;
    string flag;
    flag = "valid";

    // the flag that indicates whether file read successfully
    file.open(fname, fstream::in);
    if (file.fail())
    {
        cout << "Unable to read file" << endl;
        flag = "void";
    }

    // read in the variables
    file >> step_size;
    file >> norm_tol;
    file >> converge_tol;
    file >> max_iter;
    cout << "=======================================================================" << endl;
    cout << "Here are parameters collected from file:" << endl;
    cout << "step size is " << step_size << endl;
    cout << "tolerence to norm of Hessian is " << norm_tol << endl;
    cout << "tolerence to convergence is " << converge_tol << endl;
    cout << "maximum number of iteration is " << max_iter << endl;
    cout << "=======================================================================" << endl;
    file.flush();
    file.close();
    return flag;
}



/*
the main optimization algorithm for the inversion project
this function does only one operation of minimizing search, the iteration loop is not here
ther are two optimization algorithms: quasi-newton method and gradient descend method

when the method selected is quasi-newton, it monitors the determinant of the inverse of hessian, since
it can easily get ill-conditioned (i.e. det=0). to avoid optimization faliure, the method will be switched
to gradient descend method automatically
*/
double quasi_newton_method(string& method, int i, double step_size, double norm_tol, int no_of_equ, double start_x, double start_y, double& result_x, double& result_y, double& result_residual)
{
    // creation of temporary parameters and data for the operation
    double* alp, * norm;
    point<double>* ini;
    grad<double>* dir, * res;
    hess<double>* sec_dir;

    // getting the value, the gradient and hessian at the stationary point
    ini = new point<double>(start_x, start_y, no_of_equ);
    dir = new grad<double>(*ini, step_size, no_of_equ);
    sec_dir = new hess<double>(*ini, step_size, no_of_equ);


    norm = new double;
    *norm = sec_dir->determinant();
    res = new grad<double>;
    alp = new double;

    // the checking of determinant of hessian
    // if ill-conditioned, the method is switched to gradient descend method there after
    if (*norm < norm_tol && method == "quasi_newton")
    {
        cout << "The Hessian is ill-conditioned at " << i + 1 << "th iteration, switch to gradient descend method thereafter" << endl;
        cout << "=======================================================================" << endl;
        method = "grad_des";
    }

    // the gradient descend method for optimization
    if (method == "grad_des")
    {

        *alp = sec_dir->alpha();
        *res = dir->gen_mul(*alp);
        delete alp;
    }

    // the quasi-newton method for optimization
    if (method == "quasi_newton")
    {
        *res = sec_dir->matmul(*dir);
        delete alp;
    }

    // obtain the new stationary point by shifting from the old one
    *ini = res->step(*ini, no_of_equ);


    // the resulting x, y and residuals are the starting point of the next iteration
    // they arer changed by reference
    result_x = ini->x;
    result_y = ini->y;
    result_residual = res->residual();

    // delete tempoaray data
    delete ini, dir, norm, alp, sec_dir, res;
    return 0.;
}


// the del_data function delets all the nodes of a doubly linked list, which is the recorded iterations recursively
// this is used after a single optimization operation and recording of data
void del_data(iteration_hist<double>* first)
{
    if (first == nullptr)
        return;
    else
    {
        // recursively delete data by reaching the head node first
        del_data(first->prev);
        free(first);
    }
    return;
}



// the main function
int main(void)
{
    // start of programme
    cout << "Optimization code initializing" << endl;

    // read file function to obtain basic parameters, returns a flag
    start = read_file("input.txt");

    // if the flag returned suggests no input file is detected, the programme will terminate instantly
    if (start == "void")
    {
        cout << "there is no input.txt found in the folder" << endl;
        cout << "please set up the input.txt with five numbers splitting with one space from each other" << endl;
        cout << "'step size of differenciation' 'tolerence of norm of hessian' 'tolerance of convergence' 'maximum number of iteration'" << endl;
    }

    // if read file is successful
    else
    {
        cout << "are these the right parameters? (y/n)" << endl;

        // y_n is the flag for whether the programme will run
        cin >> y_n;

        // open the file that stores the results
        fstream results;
        results.open("results.txt", fstream::out);

        // error check
        if (results.fail())
        {
            cout << "the results.txt is still open, please close it" << endl;
            exit(0);
        }

        // writing the title in the result file
        results << "Optimization Results" << endl;


        // iteration of optimization process
        // continues until user wants to stop
        int all;
        for (all = 0; all < max_iter; all++)
        {
            // constructing the doubly linked list that records the iteration steps
            // this is the "head node"
            iteration_hist<double>* first;
            first = new iteration_hist<double>();
            first->prev = nullptr;

            // if the user want to start an optimization attempt
            if (y_n == "y")
            {
                // selection of target function
                cout << "=======================================================================" << endl;
                cout << "parameter setting finished" << endl;
                cout << "please select the no. of function to minimize" << endl;
                cout << "there are 3 functions to select:" << endl;

                cout << "1, 2*x^-2 + y^-2 + 3*y^2 - 2*x^-1 + 4*x^2" << endl;
                cout << "2, 0.3*(x^2+y^2) - 0.48*x*y" << endl;
                cout << "3, sin(x*y)" << endl;
                cout << "enter the number in front of the equation" << endl;
                cin >> no_of_equ;

                // error check for equation input
                if (no_of_equ != 1 && no_of_equ != 2 && no_of_equ != 3)
                {
                    cout << "=======================================================================" << endl;
                    cout << "invalid input for minimizing equation and could cause a crash" << endl;
                    results << "=======================================================================" << endl;
                    results << "***unexpected shut down due to wrong input***" << endl;
                    break;
                }


                // initial guesses of x and y values
                // for recommended initial guess for the three functions please refer to the documentation
                cout << "please enter the x value of starting point" << endl;
                cin >> start_x;
                cout << "please enter the y value of starting point" << endl;
                cin >> start_y;

                // the selection of opt method
                cout << "please choose the method of optimization, 1 for quasi-newton method, 2 for gradient descend method" << endl;
                cin >> j;
                if (j == 1)
                {
                    method = "quasi_newton";
                }
                else if (j == 2)
                {
                    method = "grad_des";
                }
                else
                {
                    cout << "input method is wrong, please try again" << endl;
                    results << "=======================================================================" << endl;
                    results << "***unexpected shut down due to wrong input***" << endl;
                    break;
                }

                // this allows user to record the information of stationary points during the iteration
                cout << "please enter the number of iterations before each data recording" << endl;
                cin >> record_seq;
                cout << "=======================================================================" << endl;


                // the main optimization 
                int i;
                for (i = 0; i < max_iter; i++)
                {
                    // conduct point shift and obtain new values of result_x, result_y, result_residual
                    quasi_newton_method(method, i, step_size, norm_tol, no_of_equ, start_x, start_y, result_x, result_y, result_residual);

                    // if the convergence tolerence is achieved, the minimum point is found and the process stops
                    if (result_residual < converge_tol)
                    {

                        cout << "convergence successful" << endl;
                        cout << "the result is" << " x = " << result_x << ", y = " << result_y << ", residual = " << result_residual << endl;
                        cout << "=======================================================================" << endl;

                        // record the last iteration information to the iteration_hist (doubly linked list)
                        point<double>* temp;
                        temp = new point<double>(result_x, result_y, no_of_equ);
                        first->data = *temp;
                        first->time = i;
                        first->residual = result_residual;
                        first->next = nullptr;

                        break;
                    }

                    // if the maximum num of iteration is reached before the minimum point is found
                    // i.e. convergence fails
                    else if (i == max_iter - 1)
                    {
                        cout << "convergence failed after maximimum iteration, the last results is saved" << endl;
                        cout << "the last result is" << " x = " << result_x << ", y = " << result_y << ", residual = " << result_residual << endl;
                        cout << "=======================================================================" << endl;
                        point<double>* temp;
                        temp = new point<double>(result_x, result_y, no_of_equ);
                        first->next = nullptr;
                        first->data = *temp;
                        first->time = i;
                        first->residual = result_residual;

                        break;
                    }

                    // recording of stationary points during iteration into the doubly linked list
                    else if (i % record_seq == 0)
                    {
                        point<double>* temp;
                        iteration_hist<double>* current;
                        temp = new point<double>(result_x, result_y, no_of_equ);
                        current = new iteration_hist<double>();

                        // data recording
                        first->data = *temp;
                        first->time = i;
                        first->residual = result_residual;
                        first->next = current;
                        current->prev = first;
                        first = current;

                    }

                    // update the next initial point if loop is not break
                    start_x = result_x;
                    start_y = result_y;


                }

                // jumping to the firstly constructed node of the list
                // i.e. the earliest iteration step
                while (first->prev != nullptr)
                {
                    first = first->prev;
                }

                // writing result file
                // writing of header
                double x, y, val, ti, res;
                point<double>* po;
                po = new point<double>();
                results << "=======================================================================" << endl;
                results << all + 1 << "th optimization attempt" << endl;
                results << "=======================================================================" << endl;
                if (no_of_equ == 1)
                {
                    results << "minimizing function is : 2*x^-2 + y^-2 + 3*y^2 - 2*x^-1 + 4*x^2" << endl;
                }
                else if (no_of_equ == 2)
                {
                    results << "minimizing function is : 0.3*(x^2+y^2) - 0.48*x*y" << endl;
                }
                else if (no_of_equ == 3)
                {
                    results << "minimizing function is: sin(x*y)" << endl;
                }
                results << "=======================================================================" << endl;

                // writing of each recorded iteration steps
                while (first->next != nullptr)
                {
                    // get the data
                    *po = first->data;
                    x = po->x;
                    y = po->y;
                    val = po->value;
                    ti = first->time;
                    res = first->residual;

                    // writing data in the file
                    results << "the " << ti + 1 << "th iteration, x = " << x << ", y = " << y << ", with residual of " << res << endl;
                    first = first->next;

                    // write the final iteration time
                    if (first->next == nullptr)
                    {
                        *po = first->data;
                        x = po->x;
                        y = po->y;
                        val = po->value;
                        res = first->residual;
                        ti = first->time;
                        results << "the " << ti + 1 << "th (final) iteration, x = " << x << ", y = " << y << ", with residual of " << res << ", the value of equation is " << val << endl;
                    }

                }

                // delete the temporary point for writing file
                delete po;

                // recursively delet the iteration history
                del_data(first);

                // ask if the user wants to continue the optimization or terminate the code
                cout << "do you want to run for other equations or initial values? (y/n)" << endl;
                cin >> y_n;
            }

            // if the user wants to terminate the code
            else if (y_n == "n")
            {
                // write in the file that the process is terminated
                cout << "=======================================================================" << endl;
                cout << "session terminated, the optimization history is recorded in the result.txt" << endl;
                cout << "=======================================================================" << endl;
                cout << "please modify your parameter input in the input.txt file if you wish" << endl;
                results << "=======================================================================" << endl;
                results << "function terminated successfully" << endl;

                // closing down the file
                results.flush();
                results.close();
                break;
            }
            else
            {
                // writing accidental determination in the result file
                cout << "=======================================================================" << endl;
                cout << "session terminated due to wrong input, please try again";
                results << "=======================================================================" << endl;
                results << "***unexpected shut down due to wrong input***" << endl;

                // closing down the file
                results.flush();
                results.close();
                break;
            }
        }
    }
}

