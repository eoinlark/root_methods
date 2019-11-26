#########################################
# Clearing the R Environment -------------
#########################################
rm(list=ls()) 

#########################################
# Loading Required Packages  ------------
#########################################
library(Rcpp)
library(inline)
library(rbenchmark)

#########################################
# Setting up some test data  ------------
# Data will be used to test root methods
#########################################
x <- c (12.262307 , 10.281078 , 10.287090 , 12.734039 ,
         11.731881 , 8.861998 , 12.246509 , 11.244818 ,
         9.696278 , 11.557572 , 11.112531 , 10.550190 ,
         9.018438 , 10.704774 , 9.515617 , 10.003247 ,
         10.278352 , 9.709630 , 10.963905 , 17.314814)


##################################################################################
# Defining 1st & 2nd Derivative Functions ------------
##################################################################################
incl <- '
double f1(NumericVector z, double root) //first derivative function
{
  double total = 0;                               //total is the summation of all z values
  for(int i=0; i < z.size(); i++ )                //loop over vectors
  {
   total += 2*(z[i]-root)/ (1+pow(z[i]-root,2));  //cumulatively adding as iterate through function
  }
  return total;
}

double f2(NumericVector z, double root) //second derivative function
{
  double total = 0;                   //total is the summation of all z values
  for(int i=0; i < z.size(); i++ )    //loop over vectors
  {
   total += 2*( pow(z[i]-root,2) - 1 )/  ( pow( 1 + pow(z[i]-root, 2),2) );
  }
  return total;
}
'

#########################################
# Bisect Function ------------
#########################################
body_bisect <- '
NumericVector xx(x);               //converting input x to a Numeric Vector
double r1 = as<double>(a);         //converting a to a double root 1 
double r2 = as<double>(b);         //converting b to a double root 2
double r3 = 0;                     //defining a third root, initiazlizing to zero
int num_its = as<int>(n);          //setting number of iterations, reading from input value n
double toler = as<double>(tol);    //setting tolerance, reading from input varialbe tolerance

for(int j = 0; j < num_its; j++)   //iterating through the loop to solve
{
  r3 = (r1 + r2)/2;
  double total_1 = f1(xx, r1);    //calculatng the sum of first derivative using function from incl
  double total_3 = f1(xx, r3);
 
  if (fabs(total_3) == 0) {break;}    //returning r3 if total is zero
  else if(total_3*total_1 >= 0.0f) {r1 = r3;}    // checking sign, setting r1 to m
  else {r2 = r3;}                                // else setting r2 to m  
  
  if ( fabs(r1-r2) <= toler) { break;} //if below threshold, returning the root
}

if (f2(xx, r3) < 0 ){return( wrap(r3) );} //Returning the root
else {return(wrap("No MLE found"));}
'

bisectRcpp <- cxxfunction(signature(x = "numeric", 
                                    a = "numeric", 
                                    b = "numeric", 
                                    n = "integer", 
                                    tol = "numeric"),
                         body = body_bisect,                 #linking to bisect function
                         includes = incl,                    #including differentiation functions
                         plugin = "Rcpp")

#########################################
# Newton Raphson ------------
#########################################
body_newton <- '

NumericVector xx(x);                //converting input x to a Numeric Vector
double x1_old = as<double>(b);      //setting "old" root as input value b
double x1_new = 0;                  //initializing new root to zero
int num_its = as<int>(n);           //setting number of iterations equal to input value n
double toler = as<double>(tol);     //setting tolerance toler to input value tol

for(int j = 0; j < num_its; j++)    //looping for the number of iterations
{
  double total_1 = f1(xx, x1_old);  //calculatng the sum of first derivative using function from incl
  double total_2 = f2(xx, x1_old);  //similarly for second derivative 
  
  x1_new = x1_old - (total_1 / total_2);  //calculating the new x1 root

  if (fabs(x1_new - x1_old) == 0) {break;}            //returning new root if solution is zero
  else if( fabs(x1_new - x1_old) <= toler) {break;}  //returning new root if difference is less than tolerance
  else {x1_old = x1_new;}                                             //setting old root equal to the new root and continuing to iterate
}

if (f2(xx, x1_new) < 0 ){return(wrap(x1_new));} //Returning the root
else {return(wrap("No MLE found"));}
'

newtonRcpp <- cxxfunction(signature(x = "numeric", b = "numeric", n = "integer", tol = "numeric"),
                           body = body_newton,
                           include = incl,
                           plugin = "Rcpp")

#########################################
# Secant Function ------------
#########################################
body_secant <- '

NumericVector xx(x);             //converting input x to a Numeric Vector
double x0 = as<double>(a);       //setting x0 as input value a
double x1 = as<double>(b);       //setting x1 as input value b
double x2 = 0;                   //initializing  x2 as zero  
int num_its = as<int>(n);        //setting number of iterations equal to n
double toler = as<double>(tol);  //setting tolerance equal to input value tol

for(int j = 0; j < num_its; j++)
{
  double total_0 = f1(xx, x0);      //calculatng the sum of first derivative using function from incl
  double total_1 = f1(xx, x1);
  x2 = x1 - total_1 * (x1 - x0)/(total_1 - total_0);  //calculating value of new root x2
  double total_2 = f1(xx, x2);                        //getting output of first derivative with x2
  
  if (total_2 == 0) {break;}              //returing x2 if derivative is zero
  else if( fabs(x2 - x1) <= toler) {break;}    //returing x2 if difference is less than threshold
  else {x0 = x1; x1 = x2; }                                 //updating values for next iteration      
}

if (f2(xx, x2) < 0 ){return( wrap(x2) );}   //Returning the root
else {return(wrap("No MLE found"));}
'

secantRcpp <- cxxfunction(signature(x = "numeric", a = "numeric", b = "numeric", n = "integer", tol = "numeric"),
                           body = body_secant,
                           include = incl,
                           plugin = "Rcpp")

#########################################
# Testing Functions ------------
#########################################
bisectRcpp(x, 9,11,1000,0.00000001)
bisectRcpp(x, 0,30,1000,0.00000001)
bisectRcpp(x, -10,50,1000,0.00000001)
bisectRcpp(x, 15,50,1000,0.00000001)

newtonRcpp(x, 11,100,0.00000001)
newtonRcpp(x, 9,100,0.00000001)
newtonRcpp(x, 13,100,0.00000001)

secantRcpp(x, 10,12, 1000, 0.00000001)
secantRcpp(x, 7, 12, 1000, 0.00000001)
secantRcpp(x, 8, 13, 1000, 0.00000001)
secantRcpp(x, 8, 14, 1000, 0.00000001)

#########################################
# Benchmarking Functions ------------
#########################################

benchmark(replications = 25000,
          bisectRcpp(x, a=9, b=11, n=1000, tol =0.00000001),
          secantRcpp(x, a=8, b=11, n=1000, tol =0.00000001),
          newtonRcpp(x, b=11, n=100, tol=0.00000001),
          order='relative')

#########################################
# Welfords Algorithim ------------  
#########################################

body_varCpp <- '
 NumericVector vec(x);
 NumericVector xx;

 int remNA = as<int>(na_rm);
 if(remNA == 1){ //if true
   xx = vec[!is_na(vec)]; //xx has missings removed
   }else{
   xx = vec; //xx does not have missings removed
 }


 int n = xx.size(); // length of vec
 double M_k = xx[ 0 ]; // mean of first k elements
 double S_k = 0;  // sum of squared differences of first k elements
 double M_prev;  // mean of first (k-1) elements

 //if there is only one element, return variance = 0
 if( n < 2){ return wrap(0); }  

 //iterate over 2nd to nth element 
 for( int k = 2; k <=n; k++ ){
  M_prev = M_k ; //mean of first (k-1) elements
  M_k += ( xx[ k-1 ] - M_k) / k ; //update M_k: mean of first k elements

  //update S_k sum of squared differences of first k elements
  S_k += ( xx[ k-1 ] - M_prev) * ( xx[ k-1 ] - M_k ); 
 }

 return wrap(  S_k / (n-1) );
'

#Compile Welfords Variance function
varCpp <- cxxfunction( signature( x = "numeric", na_rm = "logical"),
                       body = body_varCpp,
                       plugin = "Rcpp")


#########################################
# Testing Welfords Algorithim ------------  
#########################################


#generate some data
x <- rnorm(1000)

#call compiled code and validate code output
varCpp(x, na_rm = TRUE)
varCpp(x, na_rm = FALSE)
var(x)

#now set first 3 elements to missing
x[1:3] <- NA

varCpp(x, na_rm = TRUE)
varCpp(x, na_rm = FALSE)
var(x, na.rm = TRUE)

