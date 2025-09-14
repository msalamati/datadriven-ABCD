#include <array>
#include <iostream>
#include <stdio.h>      /* printf */
#include <stdlib.h>     /* system, NULL, EXIT_FAILURE */
#include <unistd.h>     /* fork */
#include <fstream>     /* write to file*/
#include <sstream>
#include <string>
#include <cmath>

//#include "MatlabDataArray.hpp"
//#include "MatlabEngine.hpp"

#include "cuddObj.hh"
#include "SymbolicSet.hh"
#include "SymbolicModel.hh"
#include "SymbolicModelGrowthBound.hh"

#include "FixedPoint.hh"
#include "RungeKutta4.hh"
#include "TicToc.hh"

/* state space dim */
#define sDIM 3
#define iDIM 1
#define dDIM 2

using namespace std;

/*Data type for ODE solver*/
typedef std::array<double,3> state_type;
typedef std::array<double,1> input_type;
typedef std::array<double,2> dist_type;

double A[3][3] = {
    {0.00027563,  0,       0},
    {0         , -0.3951,  0.687},
    {0         , -0.6869, -0.016}
};

double B[3][1] = {
    {0.00031166},
    {0.1359},
    {0.0230}
};

double E[3][2] = {
    {0.00033103, 0.00031244},
    {0.1309,     0.1308},
    {0.0250,     0.0233}
};

double C[1][3] = {
    {-0.0115, -0.2296, 0.0412}
};

double d1 = 0.3; //disturbance 1
double d2 = 0.2; //disturbance 2

/* state space interval*/
const double samples = 3290;
const double bias = 0.0019976;
const double h_eta = 0.0015;

/*sampling time*/
const double tau = 0.4;  /*CHANGE TO 0.2 for final simulation!*/
const int nint = 5; /*number of intermediate steps*/
// nint * tau = time intervals


OdeSolver ode_solver(sDIM, nint, tau);

auto freqpost = [](state_type &x, input_type &u) -> void {
    auto system_ode = [](state_type &dxdt, const state_type &x, const input_type &u) -> void {
        dxdt[0] = A[0][0] * x[0] + A[0][1] * x[1] + A[0][2] * x[2] + B[0][0] * -u[0];
        dxdt[1] = A[1][0] * x[0] + A[1][1] * x[1] + A[1][2] * x[2] + B[1][0] * -u[0];
        dxdt[2] = A[2][0] * x[0] + A[2][1] * x[1] + A[2][2] * x[2] + B[2][0] * -u[0];
    };
    ode_solver(system_ode, x, u);
};

/* computation of the growth bound (the result is stored in r)  */
auto radpost = [](state_type &r, input_type &u) -> void {
    auto growth_bound_ode = [](state_type &drdt,  const state_type &r, const input_type &u) {
    double H = 1.5715; //Lipschitz calculated by RWD
    //L(u) = 1/tau ln(H*identity matrix);
    double L = (1/tau)*log(H);
    drdt[0] = L*r[0] + E[0][0]*d1 + E[0][1]*d2;
    drdt[1] = L*r[1] + E[1][0]*d1 + E[1][1]*d2;
    drdt[2] = L*r[2] + E[2][0]*d1 + E[2][1]*d2;
    };
    ode_solver(growth_bound_ode,r,u);
};

scots::SymbolicSet createInputSpace(Cudd &mgr) {
    
    double ilb[iDIM] = {0}; /* participation between 0% and 100%*/
    double iub[iDIM] = {0.5};
    double ieta[iDIM] = {0.025};  /*1% increments*/
    
    scots::SymbolicSet is(mgr,iDIM,ilb,iub,ieta);
    
    is.addGridPoints();
    std::cout << is.getSize() << std::endl;
    return is;
};
/****************************************************************************/
/* MAIN COMPUTATION */
/****************************************************************************/
scots::SymbolicSet createInputSpace(Cudd &mgr);
int main(){
    TicToc tt;
    Cudd mgr;															/*bdd manager to organise the variables*/
    
    /****************************************************************************/
    /* construct SymbolicSet for the state space */
    /****************************************************************************/
    
    double a[sDIM] = {-0.02, -0.05, -0.12};									/*lower bounds of the hyper rectangle*/
    double b[sDIM] = {0.02, 0.05, 0.12};								/*upper bounds of the hyper rectangle*/
    double eta[sDIM] = {h_eta, h_eta, h_eta};						/*grid node distance diameter (0.00625 for 2) */
    scots::SymbolicSet ss(mgr,sDIM,a,b,eta);			/*instantiate a symbolic set*/
    ss.addGridPoints();									/*add grid points to the SymbolicSet ss*/
    std::cout << std::endl << "Added State grid points " << std::endl;
    //ss.writeToFile("freqset.bdd");						/*writes the bdd representing the grid points in ss to freqset.bdd */
    std::cout << ss.getSize() << std::endl;
    
    /****************************************************************************/
    /* the target set */
    /****************************************************************************/
    scots::SymbolicSet target(ss); /*subset of state space*/
    target.clear(); //remove current states in target
    
    auto f=[](double* x) -> bool {
        return(
               (C[0][0] * x[0] + C[0][1]*x[1] + C[0][2]*x[2]) > -0.008
               &&
               (C[0][0] * x[0] + C[0][1]*x[1] + C[0][2]*x[2]) < 0.008);
        //return (x[0] >= -0.01 && x[1] >= -0.01 && x[2] >= -0.01);
    };
    
    
    
    target.addByFunction(f);
    std::cout << std::endl << target.getSize() << std::endl;
    /*double H[2*sDIM] = {1,1,1,
     -1,-1,-1};
     double h[3] = {0.005, 0.01};
     
     target.addPolytope(2,H,h,scots::INNER);*/
    //target.complement();
    //target.writeToFile("targetset.bdd");
    
    /****************************************************************************/
    /* the avoid set */
    /****************************************************************************/
    scots::SymbolicSet avoid(ss); /*subset of state space*/
    avoid.clear(); //remove current states in set
    
    auto g=[](double* x) -> bool {
        return((C[0][0] * x[0] + C[0][1]*x[1] + C[0][2]*x[2]) < -0.015);
    };
    
    avoid.addByFunction(g);
    
    std::cout << std::endl << avoid.getSize() << std::endl;
    
    //avoid.writeToFile("avoidset.bdd");
    /****************************************************************************/
    /* the input set - PARTICIPATION OF EVs*/
    /****************************************************************************/
    scots::SymbolicSet is=createInputSpace(mgr);
    //is.writeToFile("participationset.bdd");				/*writes the bdd representing the grid points in ss to freqset.bdd */
    
    /****************************************************************************/
    /* setup class for symbolic model computation */
    /****************************************************************************/
    scots::SymbolicSet sspost(ss,1); /*creates state space for post variables*/
    std::cout << std::endl << "Made post state " << std::endl;
    scots::SymbolicModelGrowthBound<state_type,input_type> abstraction(&ss, &is, &sspost);
    std::cout << std::endl << "SymbolicModelGrowthBound done" << std::endl;
    /*compute transition relation*/
    tt.tic();
    abstraction.computeTransitionRelation(freqpost, radpost);
    std::cout << std::endl;
    tt.toc();
    /*Get number of elements in the transition relation*/
    std::cout << std::endl << "Number of elements in the transition relation: " << abstraction.getSize() << std::endl;
    
    /****************************************************************************/
    /* controller synthesis */
    /****************************************************************************/
    scots::FixedPoint fp(&abstraction);
    std::cout << std::endl << "Building controller " << std::endl;
    BDD T = target.getSymbolicSet();
    BDD A = avoid.getSymbolicSet();
    tt.tic();
    BDD C = fp.reachAvoid(T,A,1);
    tt.toc();
    std::cout << "Controller built" << std::endl;
    /****************************************************************************/
    /* last we store the controller as a SymbolicSet */
    /****************************************************************************/
    scots::SymbolicSet controller(ss,is);
    controller.setSymbolicSet(C);
    std::cout << "Domain size: " << controller.getSize() << std::endl;
    controller.writeToFile("L_d.bdd");
    
    return 1;
    
};




// g++ freq.cc -I/opt/local/include -I//home/campus.ncl.ac.uk/b5013322/Documents/MATLAB/scots-master/bdd -I//home/campus.ncl.ac.uk/b5013322/Documents/MATLAB/scots-master/utils -L/opt/local/lib -lcudd


//./a.out
