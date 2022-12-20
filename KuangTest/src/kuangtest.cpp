#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace arma;
using namespace std;
#include <cmath>
double Fun1(double x)
{
  return 1/sqrt(2*datum::pi)*exp(-x*x/2)  ;
}


double Integral(double (*f)(double), double a, double b)
{
  double s, h;
  int n = 100, i;
  s = ((*f)(a) + (*f)(b)) / 2;
  h = (b - a) / n;
  
  for (i=1; i<n; i++)
  {
    s += (*f)(a + i * h);
  }
  return s * h ;
}



double combination(int n, int k){
  if(k==0)
  {
    return 1;
  }
  else{
  int temp=n-k;
  long long up=1;
  long long down=1;
  if(k>=temp){
    k=k+1;
      for(;k<=n;k++){
    up=up*k;
  }
      for(int i=1;i<=temp;i++){
    down=down*i;
  }
  }
  else{
    temp=temp+1;
    for(;temp<=n;temp++){
      up=up*temp;
    }
    for(int i=1;i<=k;i++){
      down=down*i;
    }
      }

  

  double result=double(up)/double(down);
  return result;
  }
}


//' Runs tests for 0-1 series and series with a given median.
//' 
//' @param run 0-1 series. If median is used, it can be any series to test its randomness.
//' @param median The defaulted value is 0.200164, it means nothing.
//' @param bigdata if the repetition of any element is over 16, it's better to set it as true. Here, central limit theorem will be applied to replace hyper-geometric distribution.
//' @return The p-value of a one-sided test.
//' @references 
//' Friedman J H, Rafsky L C. Multivariate generalizations of the Wald-Wolfowitz and Smirnov two-sample tests[J]. The Annals of Statistics, 1979: 697-717.
//' @examples
//' library(KuangTest)
//' a=c(1,0,0,0,1,1,0,1,0,0,1,0,0,1,1,1)
//' kuangtest(a)   
//' b=c(0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1)
//' kuangtest(b)
//' c<- c(5.36, 5.29, 5.58, 5.65, 5.57, 5.53, 5.62, 5.29, 5.44, 5.34, 5.79, 5.10)
//' kuangtest(c,median = median(c))
//' d<-c(1,0,0,0,1,1,0,1,0,0,1,0,0,1,1,1,0,0,0,1,1,0,1,0,0,1,0,0,1,1,1,0,0,0,1,1,0,1)
//' kuangtest(d,bigdata = TRUE)
// [[Rcpp::export]]
double kuangtest(arma::vec run, double median = 0.200164, bool bigdata=false){
  int n=run.n_elem;
  if(median!=0.200164)//mean is given
  {
    for(int i=0; i<n; i++)
    {
      if(run(i)<median)
      {
        run(i)=0;
      }
      else
      {
        run(i)=1;
      }
    }
    
  }
  
  
  
  
  int n1=sum(run);
  int n0=n-n1;
  int countall=1;
  int count0  =0;
  
  
  for(int i=0; i<n-1; i++)// Get Runs
  {
    if(run(i)!=run(i+1))
    {
      countall++;
    }
  }
  double R=countall;
  for(int i=0; i<n-1; i++)
  {
    if(run(i)!=run(i+1)&&run(i)==0)
    {
      count0++;
    }
  }
  if(run(n-1)==0){count0++;}
  
  
  
  if(bigdata==true||n1>=15||n0>=15)// Bigdata:Normal Distribution,N(0,1)
  {
    double E = 2*n1*n0/(n1+n0)+1;
    double sigma = sqrt(2*n1*n0*(2*n1*n0-n1-n0)/((n1+n0)*(n1+n0)*(n1+n0-1)));
    double Z = (R-E)/sigma;
    double result = Z;
    double bigdata_pvalue = Integral(Fun1, abs(Z), 50.0);
    cout<<"Runs Test(One-sided)"<<endl;
    cout<<" "<<endl;
    cout<<"The Statistics of standard normal distribution: "<< result<<endl;
    cout<<"The p value: "<< bigdata_pvalue<<endl;
    return bigdata_pvalue;
  }
  
  
  else{// Smalldata:Hypergeometric Distribution
    //  double test1=combination(16,8);
    //  double test2=16*15/8;
    
    //default series
    double accup=0;
    
    for(int r=2; r<=R; r++)
    {
      double P=0;
      if(r%2==0)
      {
        int k=r/2;
        P=2*combination(n1-1,k-1)*combination(n0-1,k-1)/combination(n,n1);
      }
      else{
        int k=(r-1)/2;
        P=(combination(n1-1,k)*combination(n0-1,k-1)+combination(n1-1,k-1)*combination(n0-1,k))/combination(n,n1);
      }
      accup=accup+P;
    }
    
    // vec result={R,accup,test1,test2,n1};
    double result = accup;
    cout<<"Runs Test(One-sided)"<<endl;
    cout<<" "<<endl;
    cout<<"The p value: "<< result<<endl;

    return result;
  }
  
}

//' Runs tests for 1-n series with repetition. This test is designed to judge the fairness of ranks.
//' 
//' @param rank 1-n series with repetition. 
//' @param max number of n.
//' @param method This test will generate many independent p-values according to each pair, Methods are used to adjust the level of significance.
//' @param alpha The defaulted value is 0.05.
//' @return The p-value of each one-sided test. Plus, A matrix contains each pair and whether we should reject H0 considering their p-values and the level of significance.
//' @references 
//' Friedman J H, Rafsky L C. Multivariate generalizations of the Wald-Wolfowitz and Smirnov two-sample tests[J]. The Annals of Statistics, 1979: 697-717.
//' 
//' Bonferroni, C. E. "Il calcolo delle assicurazioni su gruppi di teste." In Studi in Onore del Professore Salvatore Ortu Carboni. Rome: Italy, pp. 13-60, 1935.
//' 
//' Benjamini, Yoav, and Yosef Hochberg. "Controlling the false discovery rate: a practical and powerful approach to multiple testing." Journal of the Royal statistical society: series B (Methodological) 57.1 (1995): 289-300.
//' @examples
//' library(KuangTest)
//' test<-c(1,3,3,3,1,1,3,1,3,3,1,3,3,1,1,4,2,4,2,2,4,4,2,2,4,4,2,2,4,4,2,2,4,4,2)
//' Kfairnesstest(rank =test, max = 4,method = "BH")
//' Kfairnesstest(rank =test, max = 4,method = "B")
// [[Rcpp::export]]
arma::mat Kfairnesstest(arma::vec rank, int max, std::string method = "Bonferroni" , double alpha = 0.05)
  {
 
   
   int numberoftest = combination(max,2);
 
   //Create a mat to restore the results of tests
   mat R(3,numberoftest,fill::zeros);
  
  if(method=="Bonferroni"||method=="B")// The most direct way: Bonferroni
  {
    cout<<"Method: Bonferroni"<<endl;
    double Adjusted_alpha=alpha/numberoftest;
    //n(n-1)/2 extractions, from 1 to n 
    double testresult=0;
    int sequencefortest = 0;

    for(int i=1; i<=max-1; i++)
    {
      for(int j=i+1; j<=max; j++)
      { 
        
        // Store each pair into a memo;
        vec memo = rank.elem(find(rank==i||rank==j));
        memo.replace(i,0);
        memo.replace(j,1);
      
      cout<<"————————————————————————————————————"<<endl;
      cout<<"The total number of "<<i<<" and "<<j<<" : "<<memo.n_elem<<endl;

        //kuangtest for every single memo;
        testresult=kuangtest(memo);
        
        //Bonferroni Test
        
        //if(testresult<Adjusted_alpha){R.row(sequencefortest)<<vec({i,j,1});}
        if(testresult<Adjusted_alpha)
        {R(3*sequencefortest)   =i;
         R(3*sequencefortest+1) =j; 
         R(3*sequencefortest+2) =1; }
        else
        {R(3*sequencefortest)   =i;
         R(3*sequencefortest+1) =j; }
        cout<<endl;
        sequencefortest++;
        
      }
    }
 
   
    // Here, we've got all of p values for each pair
    cout<<"Each column presents a pair and a test result. 1 means reject H0."<<endl;
    cout<<"alternative hypothesis: nonrandomness"<<endl;
    
    
  }
  
  if(method=="Benjaminiand Hochberg"||method=="BH")// Ranked and Optimized
  { 
    cout<<"method: Benjaminiand Hochberg"<<endl;
    double Adjusted_alpha=alpha/numberoftest;
    //n(n-1)/2 extractions, from 1 to n 
    vec testresult(numberoftest,fill::zeros);
    vec ADJtestresult(numberoftest,fill::zeros);
    int sequencefortest = 0;
    
    for(int i=1; i<=max-1; i++)
    {
      for(int j=i+1; j<=max; j++)
      { 
        
        // Store each pair into a memo;
        vec memo = rank.elem(find(rank==i||rank==j));
        memo.replace(i,0);
        memo.replace(j,1);
        cout<<"————————————————————————————————————"<<endl;
        cout<<"The total number of "<<i<<" and "<<j<<" : "<<memo.n_elem<<endl;
       
        
        // kuangtest for every single memo;
        testresult(sequencefortest) = kuangtest(memo);
        //Benjaminiand Hochberg Test
        
        
        R(3*sequencefortest)   =i;
        R(3*sequencefortest+1) =j; 
        
        
        //if(testresult<Adjusted_alpha){R.row(sequencefortest)<<vec({i,j,1});}
  
        cout<<endl;
        sequencefortest++;
        
      }
    }
    
    for(int i=0; i<numberoftest; i++)
    {
      int rank=1;
      for(int j=0; j<numberoftest; j++)
      { 
        
        if(testresult(i)>testresult(j))
        {
          rank++;
        }
        
      }
      if(testresult(i)<Adjusted_alpha*rank)
      {
        ADJtestresult(i)=1;
      }
    }

    for(int i=0; i<numberoftest; i++)
    {
      
      R(3*i+2) =ADJtestresult(i); 
  
      }
    for(int i=numberoftest-1; i>=0; i--)
    {
      
      if(ADJtestresult(i)==1)
        {
        printf("The number of Signals is %d\n",i+1);
         i=-1;} 
    }

    // Here, we've got all of p values for each pair
    cout<<"Each column presents a pair and a test result. 1 means reject H0."<<endl;
    cout<<"alternative hypothesis: nonrandomness"<<endl;
   
    
  }
     return R;
 }
  
  
  
  

 
  


