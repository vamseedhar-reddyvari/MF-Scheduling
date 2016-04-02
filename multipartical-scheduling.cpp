//
// Created by Vamseedhar on 3/26/16.
/*
 * - Create User class with state as 'q'
 * - Assume a policy and calculate emperical rho
 * - Do value iteration and find the optimal policy
 * - From optimal policy conduct auction and do the state change
 * - Calculate state distribution and using this calculate gamma
 * - If distance between rho and gamma is small terminate else do loop again
 *
 */
#include"multipartical-scheduling.h"
string OUTPUT_FILE_PATH = "/Users/Vamsee/Work/GitHub/MFE-Scheduling/";
//Classes
class User{
public:
    User();
    static uint user_id;
    uint id;
    double queue_state;
};
User::User(){
//        id = User::user_id;
//        User::user_id = User::user_id +1;
    queue_state = rand()%MAXQUEUE;
}

void write_func_queue(darray &function,string file_name);
void write_func_bid(darray &function,string file_name);
double update_rho_function(darray &Pi, darray &optimal_policy, darray &rho);
double update_pi(vector<User *> &user_array,darray &PI);
void value_iteration(darray &rho, darray &V, darray const &Phi,darray &optimal_policy);
void state_change(vector<User *> user_array,darray &optimal_policy);
double cost( double st);
void Initialize(darray  &V);
double avg_over_arrivals(double st, darray & V, darray const &Phi);
double integratesuccprob(darray & ro, double uplimit);
double find_norm_diff(darray& V, darray& Vn, int tillS);
double min(double x,double y){
    if (x<=y) return x;
    else return y;
}
double max(double x,double y){
    if (x>=y) return x;
    else return y;
}

int main(int argc, const char * argv[]) {
//    User::user_id = 0;
    // Create an array of user pointers
    vector<User *> user_array(Total_Users);
    for (auto &user_ptr : user_array) {
        user_ptr = new User();
    }

    darray Phi(MAXARRIVAL+1);
    for (unsigned int q = 0; q < MAXARRIVAL; q++){
        Phi[q] = 1.0 / (MAXARRIVAL);         // Initialize Phi = U [0, MAXARRIVAL*PRECISION]
    }

    darray rho(MAXBID);
    rho[0] = 0;
    for(int i=1;i<MAXBID;i++){
        rho[i] = min(rho[i-1] + 0.001,1);
    }

    darray optimal_policy(MAXQUEUE-MAXARRIVAL);
    darray V(MAXQUEUE);
    darray PI(MAXQUEUE);

    double distance = 10;
    while(distance >pow(10,-2)) {
        // Value iteration step. Find optimal policy
        value_iteration(rho, V, Phi, optimal_policy);
        write_func_queue(V,"value-func.txt");
        // State change step by conducting auction
        state_change(user_array, optimal_policy);
        write_func_bid(optimal_policy,"optimal-policy.txt");
        // PI Update
        update_pi(user_array, PI);
        write_func_queue(PI,"pi-func.txt");
        // rho update
        distance = update_rho_function(PI, optimal_policy, rho);
        write_func_bid(rho,"rho.txt");

        cout << "MFE Iteration Distance: " << distance << endl;
    }
    return 0;
}


/********************************************************
 *         Find steady state and gamma                  *
 *******************************************************/
/*
 * - Find the PDF PI with q
 * - Then for each bid find q and then probability from PI to get gamma
 */
double update_pi(vector<User *> &user_array,darray &PI){
    for(int i=0;i<MAXQUEUE;i++){
        PI[i] = 0;
    }
    for (auto user:user_array){
        double q = user->queue_state;
        PI[q] = PI[q]+1;
    }
    for(int i=0;i<MAXQUEUE;i++){
        PI[i] = (double) PI[i]/Total_Users;
    }
    cout<<"PI update done"<<endl;
}

/**************************************************************************************/
/*                     compute new rho                                                */
/**************************************************************************************/

double update_rho_function(darray &Pi, darray &optimal_policy, darray &rho)
{
    int x=0,q=0;
    double currbid;
    darray newrho(MAXBID);

    // compute Pi[{ q: bid[q] < x}] for all x
    for( int x=0; x< MAXBID; x++)
    {

        currbid = (double) x*PRECISION_BID;
        if(x>0) newrho[x] = newrho[x-1];
        while((optimal_policy[q] <= currbid) && (q < QFINAL))
        {
            newrho[x] += Pi[q];
            q++;
            if( q >= QFINAL) break;
        }

    }
    // compute difference between new and the old one
    double diff = find_norm_diff(newrho, rho, QFINAL);
    for(x=0;x<MAXBID;x++)
        rho[x] = newrho[x];
    return(diff);
}

/********************************************************
 *         State Change and Auction Step                *
 *******************************************************/
void state_change(vector<User *> user_array,darray &optimal_policy){
// we have an optimal policy. Based on their state we take bid from all servers and appropriately allocated the job and change the states

    uarray indexes(Total_Users);
    for(int i=0;i<Total_Users;i++){
        indexes[i] = i;
    }
    random_shuffle(indexes.begin(), indexes.end() );

    for(int cell_id =0;cell_id < No_Cells;cell_id++) {

        int index_starting_server = (cell_id * No_Users_Per_Cell);
        User *user_allocated= user_array[ indexes[index_starting_server] ];
        double winner_bid = 0;
        for (int j = 0; j < No_Users_Per_Cell; j++) {
            User *user= user_array[ indexes[index_starting_server + j] ];
            // Find the maximum bid and corresponding server
            // Instead of this we should conduct auction in each cell
            uint queue_state = user->queue_state;
            double user_bid = optimal_policy[queue_state];
            if (user_bid > winner_bid) {
                winner_bid = user_bid;
                user_allocated = user;
            }
        }
        // User found. Apply queue dynamics now.

        // Remove packet from the server who won
        user_allocated->queue_state = round(max((user_allocated->queue_state)*PRECISION_STATE -1,0)/PRECISION_STATE);
        // Add the arrival packets to each user
        for (auto user: user_array){
            double arrival_pkt_size = (rand()%100)*PRECISION_STATE;
            user->queue_state  = min(round( ((user_allocated->queue_state)*PRECISION_STATE + arrival_pkt_size)/PRECISION_STATE), MAXQUEUE);
            //Regenration
            if(rand()%100<BETA*100){
                user->queue_state  = (rand()%MAXREGEN);
            }
        }
    }

    cout<<" Job service update Done"<<endl;
}


/********************************************************
 *         Value Iteration Step                         *
 *******************************************************/
void value_iteration(darray &rho, darray &V, darray const &Phi, darray &Optimal_Policy)
{
    double curstate, nextstate_on_succ;
    double service =SERVICE*PRECISION_STATE,diff;
    double beta = BETA;
    double Vsucc=0, Vfail =0,DeltaV;
    int nextq,count=0;
    darray Vnew(MAXQUEUE);
    int  LQueue = MAXQUEUE-MAXARRIVAL;  /* V vector is updated from [0, LQueue]*/

    // Initialize Value function
    Initialize(V);
    // Value iteration starts
    while(count <1000)
    {
        for(int q=0; q< LQueue;q++)
        {
            curstate    	    =  ((double)q)* PRECISION_STATE;
            nextstate_on_succ   =  (curstate -service <0)? 0 : (curstate-service);
            Vsucc    	    =  avg_over_arrivals(nextstate_on_succ, V, Phi); // E_A [ V [ max(q-service,0) + A ]
            Vfail    	    =  avg_over_arrivals(curstate, V, Phi);        // E_A [ V [ q + A ]
            DeltaV   	    =  Vfail - Vsucc; //  E_A [ V [ q + A ] - E_A [ V [ max(q-service,0) + A ]
            //printf("%f\t%f\n", curstate,nextstate_on_succ);
            Vnew[q]     	    =  cost(curstate) + beta*Vfail - integratesuccprob(rho, beta*DeltaV); //  V = C(q) + E_A [ V [ q + A ] - int_{0}^{beta Delta V} p_{\rho}(x)
            Optimal_Policy[q] =  beta*DeltaV;
        }
        /* find difference between Vnew and V */
        diff = find_norm_diff(V, Vnew, QFINAL);

        /* copy to original */
        for(int q=0;q < LQueue;q++)
        {
            V[q] = Vnew[q];
        }
        /* break if converge */
        if(diff < 1) break;
        count = count+1;
        printf("%d\t%f\n", count, diff);
    }
    cout<<"Value Iteration Converged"<<endl;
}

/**************************************************************************************/
/*       Initialize V such that T_{\rho} V_0 < V_0                                    */
/*       This guarantees monotonic decrase of value function in each step
         of value iteration ( Puterman 6.10 )                                         */
/**************************************************************************************/

void Initialize(darray  &V)
{
    double state, V1,beta = BETA;
    int count;
    for(unsigned int q=0; q< MAXQUEUE; q++)
    {
        count =1;
        state = (double)(q*PRECISION_STATE);
        // V[q] = sum_i \beta^i ( cost(state+ i* MAXARRIVAl))
        V[q]  = cost(state);
        while(count < 1000)
        {
            V1 = pow(beta,(double)count)*cost(state + count*MAXARRIVAL*PRECISION_STATE);
            //printf("%d\t%f\n", count, V1);
            V[q] += V1;
            if( V1 < PRECISION_STATE) break;
            count++;
        }
    }

}
/**************************************************************************************/
/*                     Computes holding cost function at each state                   */
/**************************************************************************************/

double cost( double st)
{
    // C(x) = x*x
    return (st*st);
}
/**************************************************************************************/
/*                     Computes expectation over arrival process                      */
/**************************************************************************************/

double avg_over_arrivals(double st, darray & V, darray const &Phi)
{
    unsigned int q  = (unsigned int)(st/PRECISION_STATE);
    double avg      = 0;
    int maxarrivals = MAXARRIVAL;

    /* averaging over the arrival process   */
    for(unsigned int j=0; j < maxarrivals; j++)
        avg += (V[q+j]*Phi[j]);
    avg = avg;
    return(avg);
}
/**************************************************************************************/
/*                     Computes the integral of success probability                   */
/**************************************************************************************/

double integratesuccprob(darray & ro, double uplimit)
{
    double integral =0;
    int ul = (int)(uplimit/PRECISION_BID);

    //printf("ul:%lu\n", ul);
    for(int j=0; j< ul; j++)
    {
        if( j < MAXBID)
            integral += (pow(ro[j],(double)(M-1))); // succ prob, rho(x)^ (M-1)
        else
            integral += 1.0;

    }
    integral = integral * PRECISION_BID;
    return(integral);
}
/**************************************************************************************/
/*       find distance |V - Vn| till "tillS"                                          */
/**************************************************************************************/

double find_norm_diff(darray& V, darray& Vn, int tillS)
{
    //sup norm
    double wnorm;
    double d=0;
    for (int q=0; q< tillS; q++)
    {
        wnorm = 1.0; // max( cost(q*PRECISION),1);
        double d1 = fabs( V[q] - Vn[q] );
        d1 = d1/wnorm;
        if( d1 > d) d= d1;
    }
    return(d);
}
darray generate_linspace(double start,double end, int No_Values){
    darray list(No_Values);
    if(No_Values ==1){
        list[0] = start;
    }
    else {
        double precision = (double) (end - start) / (No_Values - 1);
        list[0] = start;
        for (int i = 0; i < No_Values; i++) {
            list[i] = start + i * precision;
        }
    }
    return list;
}
iarray generate_intspace(double start,double end, int No_Values){

    iarray list(No_Values);
    int precision = 0;
    if (No_Values>1){
        precision = (int) (end-start)/(No_Values-1);
    }
    else{
        precision  =0;
    }
    list[0] = start;
    for(int i=0;i<No_Values;i++){
        list[i] = start + i*precision;
    }
    return list;
}
void write_func_queue(darray &function,string file_name){
    ofstream out_p_rho_file;
    out_p_rho_file.open (OUTPUT_FILE_PATH+file_name);

    out_p_rho_file << "queue"<< " "<<"value"<<endl;
    for(int i=0;i<MAXQUEUE;i++) {
        out_p_rho_file << i *PRECISION_STATE<< " " << function[i]<<endl;
    }
    out_p_rho_file.close();
}

void write_func_bid(darray &function,string file_name){
    ofstream out_p_rho_file;
    out_p_rho_file.open (OUTPUT_FILE_PATH+file_name);
    out_p_rho_file << "bid"<< " "<<"value"<<endl;
    for(int i=0;i<MAXBID;i++) {
        out_p_rho_file << i *PRECISION_BID<< " " << function[i]<<endl;
    }
    out_p_rho_file.close();
}
