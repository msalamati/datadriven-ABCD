#include <array>
#include <iostream>
#include <stdio.h>      /* printf */
#include <stdlib.h>     /* system, NULL, EXIT_FAILURE */
#include <unistd.h>     /* fork */
#include <fstream>     /* write to file*/
#include <sstream>
#include <string>
#include <cmath>
#include <list>
#include "gurobi_c++.h"

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

/*Adaptive Method*/
typedef std::array<double,17> mat;
list<mat> map;
int global_flag = 0;
double h_eta = 0.012;

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
const double samples = 4460;
const double bias = 0.0029321;

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

auto samplepost = [](state_type &x, input_type &u, dist_type &w) -> void {
    auto system_ode = [](state_type &dxdt, const state_type &x, const input_type &u, dist_type &w) -> void {
        dxdt[0] = A[0][0] * x[0] + A[0][1] * x[1] + A[0][2] * x[2] + B[0][0] * -u[0] + E[0][0]*w[0] + E[0][1]*w[1];
        dxdt[1] = A[1][0] * x[0] + A[1][1] * x[1] + A[1][2] * x[2] + B[1][0] * -u[0] + E[1][0]*w[0] + E[1][1]*w[1];
        dxdt[2] = A[2][0] * x[0] + A[2][1] * x[1] + A[2][2] * x[2] + B[2][0] * -u[0] + E[2][0]*w[0] + E[2][1]*w[1];
    };
    ode_solver(system_ode, x, u, w);
};

/* computation of the growth bound (the result is stored in r)  */
auto radpost = [](state_type &r, input_type &u) -> void {
    if(global_flag == 0){
        try {
            double DIST_MAX1 = d1;
            double DIST_MAX2 = d2;
            
            // Create an environment
            GRBEnv env = GRBEnv(true);
            //remove all logs
            env.set(GRB_IntParam_OutputFlag, 0);
            env.set(GRB_IntParam_LogToConsole,0);
            //start environment
            env.start();
            
            // Create an empty model
            GRBModel model = GRBModel(env);
            // Create variables
            
            /*
             M =             U =
             [m11 m12 m13;       [-  u2  u3;
             m21 m22 m23;        u4  -  u6;
             m31 m32 m33];       u7  u8  -];
             */
            
            // diagonal values (can be anything)
            GRBVar m11 = model.addVar(-100, 100, 0, GRB_CONTINUOUS, "m11");
            GRBVar m22 = model.addVar(-100, 100, 0, GRB_CONTINUOUS, "m22");
            GRBVar m33 = model.addVar(-100, 100, 0, GRB_CONTINUOUS, "m33");
            
            // off-diagonal values (only positive)
            GRBVar m12 = model.addVar(0, 100, 0, GRB_CONTINUOUS, "m12");
            GRBVar m13 = model.addVar(0, 100, 0, GRB_CONTINUOUS, "m13");
            GRBVar m21 = model.addVar(0, 100, 0, GRB_CONTINUOUS, "m21");
            GRBVar m23 = model.addVar(0, 100, 0, GRB_CONTINUOUS, "m23");
            GRBVar m31 = model.addVar(0, 100, 0, GRB_CONTINUOUS, "m31");
            GRBVar m32 = model.addVar(0, 100, 0, GRB_CONTINUOUS, "m32");
            
            // boundary values for diagonals (for optimisation)
            GRBVar u11 = model.addVar(0, 100, 0, GRB_CONTINUOUS, "u11");
            GRBVar u22 = model.addVar(0, 100, 0, GRB_CONTINUOUS, "u22");
            GRBVar u33 = model.addVar(0, 100, 0, GRB_CONTINUOUS, "u33");
            
            //disturbance
            GRBVar w1 = model.addVar(0, 100, 0, GRB_CONTINUOUS, "w1");
            GRBVar w2 = model.addVar(0, 100, 0, GRB_CONTINUOUS, "w2");
            GRBVar w3 = model.addVar(0, 100, 0, GRB_CONTINUOUS, "w3");
            
            // Set objective: maximize m11 + m22 + m33 + u2 + u3 + u4 + u6 + u7 + u8
            model.setObjective(u11 + u22 + u33 + m12 + m13 + m21 + m23 + m31 + m32 + w1 + w2 + w3, GRB_MINIMIZE);
            
            // Transitions modelled as from pre to post
            // rather than from first to next
            
            
            //setup other values
            state_type init_pre;
            state_type init_post;
            state_type other_pre;
            state_type other_post;
            state_type lhs;
            state_type rhs;
            double random;
            dist_type randI;
            dist_type randO;// Add constraints by for loop
            for (int i = 0; i < samples; i++)
            {
                for (int j = 0; j < sDIM; j++)
                {
                    random = ((double) rand()) / (double) RAND_MAX;
                    //cout << random << endl;
                    random = 2*random - 1;
                    random = random * h_eta/2;
                    //random value +/- h of r for each dimension
                    other_pre[j] = r[j] + random;
                }
                
                for (int j = 0; j < sDIM; j++)
                {
                    random = ((double) rand()) / (double) RAND_MAX;
                    //cout << random << endl;
                    random = 2*random - 1;
                    random = random * h_eta/2;
                    //random value +/- h of r for each dimension
                    init_pre[j] = r[j] + random;
                }
                
                init_post = init_pre;
                other_post = other_pre;
                
                random =  ((double) rand()) / (double) RAND_MAX;
                random = 2*random - 1;
                random = random * DIST_MAX1;
                randI[0] = random;
                
                random =  ((double) rand()) / (double) RAND_MAX;
                random = 2*random - 1;
                random = random * DIST_MAX2;
                randI[1] = random;
                
                random =  ((double) rand()) / (double) RAND_MAX;
                random = 2*random - 1;
                random = random * DIST_MAX1;
                randO[0] = random;
                
                random =  ((double) rand()) / (double) RAND_MAX;
                random = 2*random - 1;
                random = random * DIST_MAX2;
                randO[1] = random;
                
                
                samplepost(init_post, u, randI);
                samplepost(other_post, u, randO); //updates value other_post with actual post value
                
                for (int j = 0; j < sDIM; j++)
                {
                    rhs[j] = abs(other_post[j] - init_post[j]);
                    lhs[j] = abs(other_pre[j] - init_pre[j]);
                }
                
                model.addConstr(-(lhs[0]*m11 + lhs[1]*m12 + lhs[2]*m13 + w1) <= -rhs[0]);
                model.addConstr(-(lhs[0]*m21 + lhs[1]*m22 + lhs[2]*m23 + w2) <= -rhs[1]);
                model.addConstr(-(lhs[0]*m31 + lhs[1]*m32 + lhs[2]*m33+ w3) <= -rhs[2]);
            }
            
            //absolute values for optimising diagonals of matrix
            model.addConstr( m11 - u11 <= 0);
            model.addConstr(-m11 - u11 <= 0);
            model.addConstr( m22 - u22 <= 0);
            model.addConstr(-m22 - u22 <= 0);
            model.addConstr( m33 - u33 <= 0);
            model.addConstr(-m33 - u33 <= 0);
            
            // Optimize model
            model.optimize();
            
            
            map.push_front({h_eta,r[0],r[1],r[2],u[0],m11.get(GRB_DoubleAttr_X), m12.get(GRB_DoubleAttr_X),m13.get(GRB_DoubleAttr_X),m21.get(GRB_DoubleAttr_X), m22.get(GRB_DoubleAttr_X),m23.get(GRB_DoubleAttr_X),m31.get(GRB_DoubleAttr_X),m32.get(GRB_DoubleAttr_X),m33.get(GRB_DoubleAttr_X),w1.get(GRB_DoubleAttr_X),w2.get(GRB_DoubleAttr_X),w3.get(GRB_DoubleAttr_X)});
            
            //int optimstatus = model.get(GRB_IntAttr_Status);
            //if (optimstatus == GRB_OPTIMAL) {
            //return r values based on optimised growth bound
            r[0] = m11.get(GRB_DoubleAttr_X)*h_eta/2   + m12.get(GRB_DoubleAttr_X)*h_eta/2   + m13.get(GRB_DoubleAttr_X)*h_eta/2 + w1.get(GRB_DoubleAttr_X) + bias;
            r[1] = m21.get(GRB_DoubleAttr_X)*h_eta/2   + m22.get(GRB_DoubleAttr_X)*h_eta/2   + m23.get(GRB_DoubleAttr_X)*h_eta/2 + w2.get(GRB_DoubleAttr_X) + bias;
            r[2] = m31.get(GRB_DoubleAttr_X)*h_eta/2   + m32.get(GRB_DoubleAttr_X)*h_eta/2  + m33.get(GRB_DoubleAttr_X)*h_eta/2 + w3.get(GRB_DoubleAttr_X) + bias;
            //}
        }
        catch(GRBException e) {
            cout << "Error code = " << e.getErrorCode() << endl;
            cout << e.getMessage() << endl;
        } catch(...) {
            cout << "Exception during optimization" << endl;
        }
    }else{
        for (auto const& i : map) {
            
            if ((r[0] >= i[1] - i[0]/2) && (r[0] <= i[1] + i[0]/2) && (r[1] >= i[2] - i[0]/2) && (r[1] <= i[2] + i[0]/2) && (r[2] >= i[3] - i[0]/2) && (r[2] <= i[3] + i[0]/2) && (i[4] == u[0]) && (i[5] == u[1])){
                r[0] = i[6]*h_eta/2    + i[7]*h_eta/2  +i[8]*h_eta/2 + i[15]+bias;
                r[1] = i[9]*h_eta/2    + i[10]*h_eta/2 +i[11]*h_eta/2 +i[16] +bias;
                r[2] = i[12]*h_eta/2   + i[13]*h_eta/2 +i[14]*h_eta/2 +i[17]+bias;
            }
        }
    }
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

int run(Cudd &mgr){
    /* to measure time */
    TicToc tt;
    /****************************************************************************/
    /* construct SymbolicSet for the state space */
    /****************************************************************************/
    
    double a[sDIM] = {-0.02, -0.05, -0.12};                                    /*lower bounds of the hyper rectangle*/
    double b[sDIM] = {0.02, 0.07, 0.12};                                /*upper bounds of the hyper rectangle*/
    double eta[sDIM] = {h_eta, h_eta, h_eta};                        /*grid node distance diameter (0.00625 for 2) */
    scots::SymbolicSet ss(mgr,sDIM,a,b,eta);            /*instantiate a symbolic set*/
    ss.addGridPoints();                                    /*add grid points to the SymbolicSet ss*/
    std::cout << std::endl << "Added State grid points " << std::endl;
    //ss.writeToFile("freqset.bdd");                        /*writes the bdd representing the grid points in ss to freqset.bdd */
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
    //is.writeToFile("participationset.bdd");                /*writes the bdd representing the grid points in ss to freqset.bdd */
    
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
    if(controller.getSize() > 0){
        controller.writeToFile("04sec_t008a015.bdd");
        return 1;
    }else{
        return 0;
    }
}

int main() {
    Cudd mgr;
    run(mgr);
    h_eta = h_eta/2;
    global_flag = 1;
    Cudd mgr2;
    run(mgr2);
    h_eta = h_eta/2;
    Cudd mgr3;
    run(mgr3);
    h_eta = h_eta/2;
    Cudd mgr4;
    run(mgr4);
    //h_eta = h_eta/2;
    //Cudd mgr5;
    //run(mgr5);
    return 1;
}


// g++ freq.cc -I/opt/local/include -I//home/campus.ncl.ac.uk/b5013322/Documents/MATLAB/scots-master/bdd -I//home/campus.ncl.ac.uk/b5013322/Documents/MATLAB/scots-master/utils -L/opt/local/lib -lcudd


//./a.out
