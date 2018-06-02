#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
using namespace std;

class Database
{
public:
Database(){
r = 0;
q = 0;
}
~Database(){}
// data
double r;
double q;
};

class TreeNode
{
public:
TreeNode(){ 
S = 0; 
V = 0; 
t = 0; 
}
~TreeNode(){}

// data
double S;
double V;
double t;
};

class Derivative
{
public:
virtual ~Derivative(){}

virtual double TerminalPayoff(double S) const{
return 0;
}

virtual int ValuationTests(TreeNode & node) const{
return 0;
}

// data
double T;

protected:
Derivative(){
T = 0;
}
};

class Option : public Derivative
{
public:
Option(){
K = 0;
isCall = false;
isAmerican = false;
}
virtual ~Option(){}

virtual double TerminalPayoff(double S) const;
virtual int ValuationTests(TreeNode & node) const;

// data
double K;
bool isCall;
bool isAmerican;
};

double Option::TerminalPayoff(double S) const{
if(isCall){
if(S > K){
return S - K;
}
else{
return 0.0;
}
}
else{
if(S < K){
return K - S;
}
else{
return 0.0;
}
}
}

int Option::ValuationTests(TreeNode & node) const{
if(isAmerican){
if(isCall){
node.V = max(node.V,max(node.S - K, 0.0));
}
else{
node.V = max(node.V, max(K - node.S, 0.0));
}
}
return 0;
}

class BinomialTree
{
public:
BinomialTree(int n);
~BinomialTree();

int FairValue(int n, const Derivative * p_derivative, const Database * p_db, double S, double sigma, double t0, double & FV);

private:
// methods
void Clear();
int Allocate(int n);

// data
int n_tree;
TreeNode **tree_nodes;
};

BinomialTree::BinomialTree(int n)
{
n_tree = 0;
tree_nodes = 0;
Allocate(n);
}

BinomialTree::~BinomialTree()
{
Clear();
}

void BinomialTree::Clear(){
if(n_tree > 0){
for(int i = 0; i <= n_tree; ++i){
delete[] tree_nodes[i];
}
delete[] tree_nodes;
}
}

int BinomialTree::Allocate(int n){
if (n <= n_tree){
return 0;
}

// deallocate old tree
Clear();

// allocate memory
try{
n_tree = n;
tree_nodes = new TreeNode*[n + 1];

for(int i = 0; i <= n; ++i){
tree_nodes[i] = new TreeNode[n + 1];
TreeNode * T_tmp = tree_nodes[i];
for(int j = 0; j <= n; ++j){
T_tmp[j].S = 0;
T_tmp[j].V = 0;
}
}
return 0;
}
catch(const bad_alloc& e){
return 1; 
}
}

int BinomialTree::FairValue(int n, const Derivative * p_derivative, const Database * p_db, double S, double sigma, double t0, double & FV){

FV = 0;
if(n < 1 || S <= 0 || p_derivative == NULL || p_db == NULL || p_derivative->T <= t0 || sigma <= 0.0){
return 1;
}

double dt = (p_derivative->T-t0)/double(n);
double df = exp(-p_db->r*dt);
double growth = exp((p_db->r - p_db->q)*dt);
double u = exp(sigma*sqrt(dt));
double d = 1.0/u;
double p_prob = (growth - d)/(u-d);
double q_prob = 1.0 - p_prob;

if(p_prob < 0.0 || p_prob > 1.0){
return 1;
}

Allocate(n);

TreeNode * node_tmp = tree_nodes[0];
node_tmp[0].S = S;
node_tmp[0].t = t0;

for(int i = 1; i <= n; ++i){

TreeNode * prev = tree_nodes[i-1];
node_tmp = tree_nodes[i];
node_tmp[0].S = prev[0].S * d;
double t = t0 + i*dt;
node_tmp[0].t = t;

for(int j = 1; j <= n; ++j){
node_tmp[j].S = node_tmp[j - 1].S * u * u;
node_tmp[j].t = t;
}
}

int i = n;
node_tmp = tree_nodes[i];
for(int j = 0; j <= n; ++j){ 
node_tmp[j].V = p_derivative->TerminalPayoff(node_tmp[j].S);
}

for(int i = n-1; i >= 0; --i){
node_tmp = tree_nodes[i];
TreeNode * node_next = tree_nodes[i+1];
for(int j = 0; j <= i; ++j){
node_tmp[j].V = df*(p_prob*node_next[j+1].V + q_prob*node_next[j].V);
p_derivative->ValuationTests(node_tmp[j]);
}
}

node_tmp = tree_nodes[0];
FV = node_tmp[0].V;
return 0;
}

class Straddle : public Derivative{
public:
Straddle(){
K = 0;
isAmerican = false;
}
virtual ~Straddle(){}

virtual double TerminalPayoff(double S) const;
virtual int ValuationTests(TreeNode & node) const;

//data
bool isAmerican;
double K;
};

double Straddle::TerminalPayoff(double S) const{
return abs(S - K);
}

int Straddle::ValuationTests(TreeNode & node) const{
if(isAmerican){
node.V = max(node.V,max(node.S - K,K - node.S)); 
}
return 0;
}

class BinaryOption : public Derivative{
public:
BinaryOption(){
K = 0;
isCall = false;
}
virtual ~BinaryOption(){}

virtual double TerminalPayoff(double S) const;
virtual int ValuationTests(TreeNode & node) const;

//data
double K;
bool isCall;
};

double BinaryOption::TerminalPayoff(double S) const{
if(isCall){
if(S >= K){
return 1.0;
}
else{
return 0.0;
}
}
else{
if(S < K){
return 1.0;
}
else{
return 0.0;
}
}
}

int BinaryOption::ValuationTests(TreeNode & node) const{
return 0;
}

class ConvertibleOption : public Derivative{
public:
ConvertibleOption(){
K = 0;
B = 0;
isAmerican = false;
}
virtual ~ConvertibleOption(){}

virtual double TerminalPayoff(double S) const;
virtual int ValuationTests(TreeNode & node) const;
 
//data
double K;
double B;
bool isAmerican;
};

double ConvertibleOption::TerminalPayoff(double S) const{
return max(S,K);
}
 
int ConvertibleOption::ValuationTests(TreeNode & node) const{
if(isAmerican){
node.V = max(node.V, node.S);
}
if(node.S >= B){
node.V = node.S;
}
return 0;
}

double cum_norm(double x){
const double root = sqrt(0.5);
return 0.5*(1.0 + erf(x*root));
}

double cbin(double S, double K, double q, double T, double t0, double sigma, double r){
double e = 2.71828182845904523536;
double x = pow(e,-r * (T - t0));
double top = log(S / K) + ((r - q)*(T - t0));
double bottom = sigma * (sqrt(T - t0));
double d1 = (top/bottom) + (.5*bottom);
double d2 = d1 - bottom;
double c = x * cum_norm(d2);
//double c = ((S*pow(e,-q*(T-t0)))*cum_norm(d1)) - ((K*pow(e,-r*(T-t0)))*cum_norm(d2));
return c;
}

double pbin(double S, double K, double q, double T, double t0, double sigma, double r){
double e = 2.71828182845904523536;
double x = pow(e,-r * (T - t0));
double top = log(S / K) + ((r - q)*(T - t0));
double bottom = sigma * (sqrt(T - t0));
double d1 = (top/bottom) + (.5*bottom);
double d2 = d1 - bottom;
double p = x * cum_norm(-d2);
//double p = ((K*pow(e,-r*(T-t0)))*cum_norm(-d2)) - ((S*pow(e,-q*(T-t0)))*cum_norm(-d1));
return p;
}

int main(){

/*

double r = 0.05;
double q = 0.00;
double T = 5.0;
double t0 = 0.00;
double B = 130;
Database db;
db.r = r;
db.q = q;
double S = 60.000;
double K = 100;
double sigma = 0.5;

int n = 100;
BinomialTree binom(n);

ConvertibleOption cOption;
cOption.K = K;
cOption.T = T;
cOption.B = B;
cOption.isAmerican = true;

double FV_cOption = 0;

binom.FairValue(n, &cOption, &db, S + 1, 0.3683, t0, FV_cOption);
cout << FV_cOption << endl;
binom.FairValue(n, &cOption, &db, S - 1, 0.3683, t0, FV_cOption);
cout << FV_cOption << endl;


for(double i = 0.3600; i <= 0.3700; i+=0.0001){
binom.FairValue(n, &cOption, &db, S, i, t0, FV_cOption);
cout << i << " " << FV_cOption << endl;
}

*/

}
