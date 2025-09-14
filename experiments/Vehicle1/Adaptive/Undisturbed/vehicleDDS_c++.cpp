/*
 * vehicle.cc
 *
 *  created on: 30.09.2015
 *      author: rungger
 */

/*
 * information about this example is given in the readme file
 *
 */




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
#include "SymbolicModelGrowthBound.hh"

#include "TicToc.hh"
#include "RungeKutta4.hh"
#include "FixedPoint.hh"

#ifndef M_PI
#define M_PI 3.14159265359
#endif

/* state space dim */
#define sDIM 3
#define iDIM 2

using namespace std;

/* data types for the ode solver */
typedef std::array<double,3> state_type;
typedef std::array<double,2> input_type;

/*Adaptive Method*/
typedef std::array<double,15> mat;
int global_flag = 0;
list<mat> map;

const double samples = 3127;
const double bias = 0.067767;

/* state space interval*/
double h_eta = 1.6;
/* sampling time */
const double tau = 0.3;
/* number of intermediate steps in the ode solver */
const int nint=5;
OdeSolver ode_solver(sDIM,nint,tau);

/* we integrate the vehicle ode by 0.3 sec (the result is stored in x)  */
auto  vehicle_post = [](state_type &x, input_type &u) -> void {
    
    /* the ode describing the vehicle */
    auto rhs =[](state_type& xx,  const state_type &x, input_type &u) {
        double alpha=std::atan(std::tan(u[1])/2.0);
        xx[0] = u[0]*std::cos(alpha+x[2])/std::cos(alpha);
        xx[1] = u[0]*std::sin(alpha+x[2])/std::cos(alpha);
        xx[2] = u[0]*std::tan(u[1]);
    };
    ode_solver(rhs,x,u);
};

/* computation of the growth bound (the result is stored in r)  */
auto radius_post = [](state_type &r, input_type &u) {
    if(global_flag == 0){
        try {
            
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
            
            // Set objective: maximize m11 + m22 + m33 + u2 + u3 + u4 + u6 + u7 + u8
            model.setObjective(u11 + u22 + u33 + m12 + m13 + m21 + m23 + m31 + m32, GRB_MINIMIZE);
            
            // Transitions modelled as from pre to post
            // rather than from first to next
            
            //set initial value to state
            state_type init_pre = r;
            //set post to state value
            state_type init_post = r;
            vehicle_post(init_post, u); //updates value init_post with actual post value
            //setup other values
            state_type other_pre;
            state_type other_post;
            state_type lhs;
            state_type rhs;
            double random;
            // Add constraints by for loop
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
                other_post = other_pre;
                vehicle_post(other_post,u); //updates value other_post with actual post value
                
                for (int j = 0; j < sDIM; j++)
                {
                    rhs[j] = abs(other_post[j] - init_post[j]);
                    lhs[j] = abs(other_pre[j] - init_pre[j]);
                }
                
                model.addConstr(-(lhs[0]*m11 + lhs[1]*m12 + lhs[2]*m13) <= -rhs[0]);
                model.addConstr(-(lhs[0]*m21 + lhs[1]*m22 + lhs[2]*m23) <= -rhs[1]);
                model.addConstr(-(lhs[0]*m31 + lhs[1]*m32 + lhs[2]*m33) <= -rhs[2]);
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
            
            map.push_front({h_eta,r[0],r[1],r[2],u[0],u[1],m11.get(GRB_DoubleAttr_X), m12.get(GRB_DoubleAttr_X),m13.get(GRB_DoubleAttr_X),m21.get(GRB_DoubleAttr_X), m22.get(GRB_DoubleAttr_X),m23.get(GRB_DoubleAttr_X),m31.get(GRB_DoubleAttr_X),m32.get(GRB_DoubleAttr_X),m33.get(GRB_DoubleAttr_X)});
            
            
            //compute infeasibility constraints (if any)
            //model.computeIIS();
            
            //int optimstatus = model.get(GRB_IntAttr_Status);
            //if (optimstatus == GRB_OPTIMAL) {
            //return r values based on optimised growth bound
            r[0] = m11.get(GRB_DoubleAttr_X)*h_eta/2   + m12.get(GRB_DoubleAttr_X)*h_eta/2   + m13.get(GRB_DoubleAttr_X)*h_eta/2 + bias;
            r[1] = m21.get(GRB_DoubleAttr_X)*h_eta/2   + m22.get(GRB_DoubleAttr_X)*h_eta/2   + m23.get(GRB_DoubleAttr_X)*h_eta/2 + bias;
            r[2] = m31.get(GRB_DoubleAttr_X)*h_eta/2   + m32.get(GRB_DoubleAttr_X)*h_eta/2  + m33.get(GRB_DoubleAttr_X)*h_eta/2 + bias;
            
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
                r[0] = i[6]*h_eta/2    + i[7]*h_eta/2  +i[8]*h_eta/2 + bias;
                r[1] = i[9]*h_eta/2    + i[10]*h_eta/2 +i[11]*h_eta/2 + bias;
                r[2] = i[12]*h_eta/2   + i[13]*h_eta/2 +i[14]*h_eta/2 + bias;
            }
        }
    }
    
};

/* forward declaration of the functions to setup the state space
 * input space and obstacles of the vehicle example */
scots::SymbolicSet vehicleCreateStateSpace(Cudd &mgr);
scots::SymbolicSet vehicleCreateInputSpace(Cudd &mgr);

void vehicleCreateObstacles(scots::SymbolicSet &obs);


int run(Cudd &mgr){
    /* to measure time */
    TicToc tt;
    /****************************************************************************/
    /* construct SymbolicSet for the state space */
    /****************************************************************************/
    scots::SymbolicSet ss=vehicleCreateStateSpace(mgr);
    ss.writeToFile("Target_Obst_SS/vehicle_ss_DDS.bdd");
    
    /****************************************************************************/
    /* construct SymbolicSet for the obstacles */
    /****************************************************************************/
    /* first make a copy of the state space so that we obtain the grid
     * information in the new symbolic set */
    scots::SymbolicSet obs(ss);
    vehicleCreateObstacles(obs);
    obs.writeToFile("Target_Obst_SS/vehicle_obst_DDS.bdd");
    
    /****************************************************************************/
    /* we define the target set */
    /****************************************************************************/
    /* first make a copy of the state space so that we obtain the grid
     * information in the new symbolic set */
    scots::SymbolicSet ts(ss);
    /* define the target set as a symbolic set */
    double H[4*sDIM]={-1, 0, 0,
        1, 0, 0,
        0,-1, 0,
        0, 1, 0};
    /* compute inner approximation of P={ x | H x<= h1 }  */
    double h[4] = {-9,9.51,-0, .51};
    ts.addPolytope(4,H,h, scots::INNER);
    ts.writeToFile("Target_Obst_SS/vehicle_target_DDS.bdd");
    
    /****************************************************************************/
    /* construct SymbolicSet for the input space */
    /****************************************************************************/
    scots::SymbolicSet is=vehicleCreateInputSpace(mgr);
    
    /****************************************************************************/
    /* setup class for symbolic model computation */
    /****************************************************************************/
    /* first create SymbolicSet of post variables
     * by copying the SymbolicSet of the state space and assigning new BDD IDs */
    scots::SymbolicSet sspost(ss,1);
    /* instantiate the SymbolicModel */
    scots::SymbolicModelGrowthBound<state_type,input_type> abstraction(&ss, &is, &sspost);
    /* compute the transition relation */
    tt.tic();
    abstraction.computeTransitionRelation(vehicle_post, radius_post);
    std::cout << std::endl;
    tt.toc();
    /* get the number of elements in the transition relation */
    std::cout << std::endl << "Number of elements in the transition relation: " << abstraction.getSize() << std::endl;
    
    /****************************************************************************/
    /* we continue with the controller synthesis */
    /****************************************************************************/
    /* we setup a fixed point object to compute reachabilty controller */
    scots::FixedPoint fp(&abstraction);
    /* the fixed point algorithm operates on the BDD directly */
    BDD T = ts.getSymbolicSet();
    BDD O = obs.getSymbolicSet();
    tt.tic();
    /* compute controller */
    BDD C=fp.reachAvoid(T,O,1);
    tt.toc();
    /****************************************************************************/
    /* last we store the controller as a SymbolicSet
     * the underlying uniform grid is given by the Cartesian product of
     * the uniform gird of the space and uniform gird of the input space */
    /****************************************************************************/
    scots::SymbolicSet controller(ss,is);
    controller.setSymbolicSet(C);
    cout << "Controller size: " << controller.getSize() << endl;
    if(controller.getSize() > 0){
    controller.writeToFile("Controllers/vehicle_controller_DDS.bdd");
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
    h_eta = h_eta/2;
    Cudd mgr5;
    run(mgr5);
    return 1;
}

scots::SymbolicSet vehicleCreateStateSpace(Cudd &mgr) {
    
    /* setup the workspace of the synthesis problem and the uniform grid */
    /* lower bounds of the hyper rectangle */
    double lb[sDIM]={0,0,-M_PI-0.4};
    /* upper bounds of the hyper rectangle */
    double ub[sDIM]={10,10,M_PI+0.4};
    /* grid node distance diameter */
    double eta[sDIM]={h_eta,h_eta,h_eta};
    
    
    scots::SymbolicSet ss(mgr,sDIM,lb,ub,eta);
    
    /* add the grid points to the SymbolicSet ss */
    ss.addGridPoints();
    
    return ss;
}

void vehicleCreateObstacles(scots::SymbolicSet &obs) {
    
    /* add the obstacles to the symbolic set */
    /* the obstacles are defined as polytopes */
    /* define H* x <= h */
    double H[4*sDIM]={-1, 0, 0,
        1, 0, 0,
        0,-1, 0,
        0, 1, 0};
    /* add outer approximation of P={ x | H x<= h1 } form state space */
    double h1[4] = {-1,1.2,-0, 9};
    obs.addPolytope(4,H,h1, scots::OUTER);
    /* add outer approximation of P={ x | H x<= h2 } form state space */
    double h2[4] = {-2.2,2.4,-0,5};
    obs.addPolytope(4,H,h2, scots::OUTER);
    /* add outer approximation of P={ x | H x<= h3 } form state space */
    double h3[4] = {-2.2,2.4,-6,10};
    obs.addPolytope(4,H,h3, scots::OUTER);
    /* add outer approximation of P={ x | H x<= h4 } form state space */
    double h4[4] = {-3.4,3.6,-0,9};
    obs.addPolytope(4,H,h4, scots::OUTER);
    /* add outer approximation of P={ x | H x<= h5 } form state space */
    double h5[4] = {-4.6 ,4.8,-1,10};
    obs.addPolytope(4,H,h5, scots::OUTER);
    /* add outer approximation of P={ x | H x<= h6 } form state space */
    double h6[4] = {-5.8,6,-0,6};
    obs.addPolytope(4,H,h6, scots::OUTER);
    /* add outer approximation of P={ x | H x<= h7 } form state space */
    double h7[4] = {-5.8,6,-7,10};
    obs.addPolytope(4,H,h7, scots::OUTER);
    /* add outer approximation of P={ x | H x<= h8 } form state space */
    double h8[4] = {-7,7.2,-1,10};
    obs.addPolytope(4,H,h8, scots::OUTER);
    /* add outer approximation of P={ x | H x<= h9 } form state space */
    double h9[4] = {-8.2,8.4,-0,8.5};
    obs.addPolytope(4,H,h9, scots::OUTER);
    /* add outer approximation of P={ x | H x<= h10 } form state space */
    double h10[4] = {-8.4,9.3,-8.3,8.5};
    obs.addPolytope(4,H,h10, scots::OUTER);
    /* add outer approximation of P={ x | H x<= h11 } form state space */
    double h11[4] = {-9.3,10,-7.1,7.3};
    obs.addPolytope(4,H,h11, scots::OUTER);
    /* add outer approximation of P={ x | H x<= h12 } form state space */
    double h12[4] = {-8.4,9.3,-5.9,6.1};
    obs.addPolytope(4,H,h12, scots::OUTER);
    /* add outer approximation of P={ x | H x<= h13 } form state space */
    double h13[4] = {-9.3,10 ,-4.7,4.9};
    obs.addPolytope(4,H,h13, scots::OUTER);
    /* add outer approximation of P={ x | H x<= h14 } form state space */
    double h14[4] = {-8.4,9.3,-3.5,3.7};
    obs.addPolytope(4,H,h14, scots::OUTER);
    /* add outer approximation of P={ x | H x<= h15 } form state space */
    double h15[4] = {-9.3,10 ,-2.3,2.5};
    obs.addPolytope(4,H,h15, scots::OUTER);
    
}

scots::SymbolicSet vehicleCreateInputSpace(Cudd &mgr) {
    
    /* lower bounds of the hyper rectangle */
    double lb[sDIM]={-1,-1};
    /* upper bounds of the hyper rectangle */
    double ub[sDIM]={1,1};
    /* grid node distance diameter */
    double eta[sDIM]={.3,.3};
    
    scots::SymbolicSet is(mgr,iDIM,lb,ub,eta);
    is.addGridPoints();
    
    return is;
}
