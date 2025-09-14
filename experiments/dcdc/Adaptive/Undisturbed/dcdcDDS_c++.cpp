/*
 * dcdc.cc
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
#include <cmath>
#include <list>
#include "gurobi_c++.h"

#include "cuddObj.hh"

#include "SymbolicSet.hh"
#include "SymbolicModelGrowthBound.hh"

#include "TicToc.hh"
#include "RungeKutta4.hh"
#include "FixedPoint.hh"

/* state space dim */
#define sDIM 2
#define iDIM 1

using namespace std;

/* data types for the ode solver */
typedef std::array<double,2> state_type;
typedef std::array<double,1> input_type;
typedef std::array<double,8> mat;

/* some model constants */
double h_eta = 0.01;
double const bias = 9.6687e-5;
double const samples = 2285;
const double r0 = 1.0;
const double vs = 1.0;
const double rl = 0.05;
const double rc = rl / 10;
const double xl = 3.0;
const double xc = 70.0;
const double b[2]={vs/xl, 0};

int global_flag = 0;
list<mat> map;

/* sampling time */
const double tau = 0.5;
/* number of intermediate steps in the ode solver */
const int nint=5;
OdeSolver ode_solver(sDIM,nint,tau);

/* we integrate the dcdc ode by 0.5 sec (the result is stored in x)  */
auto  dcdc_post = [](state_type &x, input_type &u) -> void {
    /* the ode describing the system */
    auto  system_ode = [](state_type &dxdt, const state_type &x, const input_type &u) -> void {
        double a[2][2];
        if(u[0]==1) {
            a[0][0] = -rl / xl;
            a[0][1] = 0;
            a[1][0] = 0;
            a[1][1] = (-1 / xc) * (1 / (r0 + rc));
        } else {
            a[0][0] = (-1 / xl) * (rl + ((r0 * rc) / (r0 + rc))) ;
            a[0][1] =  ((-1 / xl) * (r0 / (r0 + rc))) / 5 ;
            a[1][0] = 5 * (r0 / (r0 + rc)) * (1 / xc);
            a[1][1] =(-1 / xc) * (1 / (r0 + rc)) ;
        }
        dxdt[0] = a[0][0]*x[0]+a[0][1]*x[1] + b[0];
        dxdt[1] = a[1][0]*x[0]+a[1][1]*x[1] + b[1];
    };
    ode_solver(system_ode,x,u);
};

/* computation of the growth bound (the result is stored in r)  */
auto radius_post = [](state_type &r, input_type &u) -> void {
    if (global_flag == 0){
        try {
            //int samples = 2285;
            
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
             [m11 m12       [u1  -;
             m21 m22];       -  u4];
             */
            
            // diagonal values (can be anything)
            GRBVar m11 = model.addVar(-100, 100, 0, GRB_CONTINUOUS, "m11");
            GRBVar m22 = model.addVar(-100, 100, 0, GRB_CONTINUOUS, "m22");
            
            // off-diagonal values (only positive)
            GRBVar m12 = model.addVar(0, 100, 0, GRB_CONTINUOUS, "m12");
            GRBVar m21 = model.addVar(0, 100, 0, GRB_CONTINUOUS, "m21");
            
            // boundary values for diagonals (for optimisation)
            GRBVar u11 = model.addVar(0, 100, 0, GRB_CONTINUOUS, "u11");
            GRBVar u22 = model.addVar(0, 100, 0, GRB_CONTINUOUS, "u22");
            // Set objective: maximize m11 + m22 + m33 + u2 + u3 + u4 + u6 + u7 + u8
            model.setObjective(u11 + u22 + m12 + m21, GRB_MINIMIZE);
            
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
                    init_pre[j] = r[j] + random;
                }
                
                for (int j = 0; j < sDIM; j++)
                {
                    random = ((double) rand()) / (double) RAND_MAX;
                    //cout << random << endl;
                    random = 2*random - 1;
                    random = random * h_eta/2;
                    //random value +/- h of r for each dimension
                    other_pre[j] = r[j] + random;
                }
                init_post = init_pre;
                other_post = other_pre;
                dcdc_post(init_post,u);  // updates value init_post with actual post value
                dcdc_post(other_post,u); // updates value other_post with actual post value
                
                for (int j = 0; j < sDIM; j++)
                {
                    rhs[j] = abs(other_post[j] - init_post[j]);
                    lhs[j] = abs(other_pre[j] - init_pre[j]);
                }
                
                model.addConstr(-(lhs[0]*m11 + lhs[1]*m12) <= -rhs[0]);
                model.addConstr(-(lhs[0]*m21 + lhs[1]*m22) <= -rhs[1]);
            }
            
            //absolute values for optimising diagonals of matrix
            model.addConstr( m11 - u11 <= 0);
            model.addConstr(-m11 - u11 <= 0);
            model.addConstr( m22 - u22 <= 0);
            model.addConstr(-m22 - u22 <= 0);
            
            // Optimize model
            model.optimize();
            
            map.push_front({h_eta,r[0],r[1],u[0],m11.get(GRB_DoubleAttr_X), m12.get(GRB_DoubleAttr_X), m21.get(GRB_DoubleAttr_X), m22.get(GRB_DoubleAttr_X)});
            
            //compute infeasibility constraints (if any)
            //model.computeIIS();
            
            //int optimstatus = model.get(GRB_IntAttr_Status);
            //if (optimstatus == GRB_OPTIMAL) {
            //return r values based on optimised growth bound
            r[0] = m11.get(GRB_DoubleAttr_X)*h_eta/2   + m12.get(GRB_DoubleAttr_X)*h_eta/2+ bias;
            r[1] = m21.get(GRB_DoubleAttr_X)*h_eta/2   + m22.get(GRB_DoubleAttr_X)*h_eta/2+ bias;
            //}
            
            
        }
        catch(GRBException e) {
            cout << "Error code = " << e.getErrorCode() << endl;
            cout << e.getMessage() << endl;
        } catch(...) {
            cout << "Exception during optimization" << endl;
        }}
    else{
        for (auto const& i : map) {
            
            if ((r[0] >= i[1] - i[0]/2) && (r[0] <= i[1] + i[0]/2) && (r[1] >= i[2] - i[0]/2) && (r[1] <= i[2] + i[0]/2) && (i[3] == u[0])){
                r[0] = i[4]*h_eta/2   + i[5]*h_eta/2 + bias;
                r[1] = i[6]*h_eta/2   + i[7]*h_eta/2 + bias;
            }
        }
    }
};

int run(Cudd &mgr){
    /* to measure time */
    TicToc tt;
    
    /****************************************************************************/
    /* construct SymbolicSet for the state space */
    /****************************************************************************/
    /* setup the workspace of the synthesis problem and the uniform grid */
    /* lower bounds of the hyper rectangle */
    double lb[sDIM]={0.65,4.95};
    /* upper bounds of the hyper rectangle */
    double ub[sDIM]={1.65,5.95};
    /* grid node distance diameter */
    double eta[sDIM]={h_eta,h_eta};
    scots::SymbolicSet ss(mgr,sDIM,lb,ub,eta);
    ss.addGridPoints();
    /****************************************************************************/
    /* construct SymbolicSet for the input space */
    /****************************************************************************/
    double ilb[iDIM]={1};
    double iub[iDIM]={2};
    double ieta[iDIM]={1};
    scots::SymbolicSet is(mgr,iDIM,ilb,iub,ieta);
    is.addGridPoints();
    /****************************************************************************/
    /* setup class for symbolic model computation */
    /****************************************************************************/
    scots::SymbolicSet sspost(ss,1); /* create state space for post variables */
    scots::SymbolicModelGrowthBound<state_type,input_type> abstraction(&ss, &is, &sspost);
    /* compute the transition relation */
    tt.tic();
    abstraction.computeTransitionRelation(dcdc_post, radius_post);
    std::cout << std::endl;
    tt.toc();
    /* get the number of elements in the transition relation */
    std::cout << std::endl << "Number of elements in the transition relation: " << abstraction.getSize() << std::endl;
    
    /****************************************************************************/
    /* we continue with the controller synthesis for FG (target) */
    /****************************************************************************/
    /* construct SymbolicSet for target (it is a subset of the state space)  */
    scots::SymbolicSet target(ss);
    /* add inner approximation of P={ x | H x<= h } form state space */
    double H[4*sDIM]={-1, 0,
        1, 0,
        0,-1,
        0, 1};
    double h[4] = {-1.1,1.6,-5.4, 5.9};
    target.addPolytope(4,H,h, scots::INNER);
    std::cout << "Target set details:" << std::endl;
    target.writeToFile("Target_Obst_SS/dcdcDDS_target.bdd");
    
    /* we setup a fixed point object to compute reach and stay controller */
    scots::FixedPoint fp(&abstraction);
    /* the fixed point algorithm operates on the BDD directly */
    BDD T = target.getSymbolicSet();
    tt.tic();
    BDD C=fp.reachStay(T,1);
    tt.toc();
    /****************************************************************************/
    /* last we store the controller as a SymbolicSet
     * the underlying uniform grid is given by the Cartesian product of
     * the uniform gird of the space and uniform gird of the input space */
    /****************************************************************************/
    scots::SymbolicSet controller(ss,is);
    controller.setSymbolicSet(C);
    cout << "Controller Size: " << controller.getSize() << endl;
    if(controller.getSize() > 0){
        controller.writeToFile("Controllers/dcdcDDS_controller.bdd");
        return 1;
    }
    else{
        return 0;
    }
}


/****************************************************************************/
/* main computation */
/****************************************************************************/
int main() {
    Cudd mgr;
    run(mgr);
    global_flag = 1;
    h_eta = h_eta/2;
    Cudd mgr2;
    run(mgr2);
    h_eta = h_eta/2;
    Cudd mgr3;
    run(mgr3);
    h_eta = h_eta/2;
    //Cudd mgr4;
    //run(mgr4);
    //h_eta = h_eta/2;
    //Cudd mgr5;
    //run(mgr5);
    return 1;
}

