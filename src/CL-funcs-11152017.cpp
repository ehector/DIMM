// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppParallel.h>
#include <math.h>
#include <RcppArmadillo.h>
using namespace Rcpp; 
using namespace arma;
using namespace RcppParallel;

struct logCLnormal_worker : public Worker
{
  // source
  const RVector<double> beta;
  const RMatrix<double> cov; 
  const RMatrix<double> block_y;
  const RMatrix<double> block_x;
  const int dim;
  
  // destination 
  double sum;
  
  // initialize with source and destination
  logCLnormal_worker(const NumericVector beta, const NumericMatrix cov, const NumericMatrix block_y, const NumericMatrix block_x, const int dim) 
    : beta(beta), cov(cov), block_y(block_y), block_x(block_x), dim(dim), sum(0) {}
  
  logCLnormal_worker(const logCLnormal_worker& little_worker, Split) 
    : beta(little_worker.beta), cov(little_worker.cov), block_y(little_worker.block_y), block_x(little_worker.block_x), 
      dim(little_worker.dim), sum(0) {}
  
  void operator()(std::size_t begin, std::size_t end) {
    double sigma_2r = 0;
    double sigma_2t = 0;
    double rho_rt = 0;
    double mean_r = 0;
    double mean_t = 0;
    
    for(unsigned int i = begin; i < end; i++) {
      int it = static_cast<int>(i);
      for(int r=0; r < (dim-1); r++) {
        for(int t=r+1; t < dim; t++){
          sigma_2r = cov(r,r);
          sigma_2t = cov(t,t);
          rho_rt = cov(r,t)/pow(sigma_2r*sigma_2t, 0.5);
          
          mean_r = 0;
          mean_t = 0;
          for (unsigned int k = 0; k < beta.length(); ++k) {
            mean_r += as_scalar(block_x(it*dim+r, k)*beta[k]);
            mean_t += as_scalar(block_x(it*dim+t, k)*beta[k]);
          }
          
          sum += (-log(pow(sigma_2r*sigma_2t*(1-pow(rho_rt,2.0)), 0.5)) - 
            (pow(block_y(r, it)-mean_r, 2.0)/sigma_2r + pow(block_y(t, it)-mean_t, 2.0)/sigma_2t - 
            2*rho_rt*(block_y(r, it)-mean_r)*(block_y(t, it)-mean_t)/pow(sigma_2r*sigma_2t,0.5))/
              (2*(1-pow(rho_rt, 2.0))));
        }
      }
    }
  }
  
  // join my sum with that of another logCLnormal_worker
  void join(const logCLnormal_worker& rhs) { 
    sum += rhs.sum;
  }
};


struct eenormalmean_worker : public Worker
{
  // source
  const RVector<double> beta;
  const RMatrix<double> cov; 
  const RMatrix<double> block_y;
  const RMatrix<double> block_x;
  const int dim;
  const int len;
  
  // destination 
  RMatrix<double> sum;
  
  // initialize with source and destination
  eenormalmean_worker(const NumericVector beta, const NumericMatrix cov, const NumericMatrix block_y, const NumericMatrix block_x, 
                      const int dim, const int len, NumericMatrix sum) 
    : beta(beta), cov(cov), block_y(block_y), block_x(block_x), dim(dim), len(len), sum(sum) {}
    
    void operator()(std::size_t begin, std::size_t end) {
      double sigma_2r = 0;
      double sigma_2t = 0;
      double rho_rt = 0;
      double mean_r = 0;
      double mean_t = 0;
      double sum_1 = 0;
      
      for(unsigned int i=begin; i < end; i++) {
        int it = static_cast<int>(i);
        for(int r=0; r < (dim-1); r++) {
          for(int t=(r+1); t < dim; t++){
            sigma_2r = cov(r,r);
            sigma_2t = cov(t,t);
            rho_rt = cov(r,t)/pow(sigma_2r*sigma_2t, 0.5);
            
            mean_r = 0;
            mean_t = 0;
            for (unsigned int k=0; k < beta.size(); k++){
              mean_r += block_x(it*dim+r, k)*beta[k];
              mean_t += block_x(it*dim+t, k)*beta[k];
            }
            
            for(unsigned int k=0; k < beta.size(); k++){
              sum_1 = 0;
              for (unsigned int q=0; q < beta.size(); q++){
                sum_1 += (block_x(it*dim+r, k)*block_x(it*dim+t, q) + block_x(it*dim+t, k)*block_x(it*dim+r, q))*beta[q];
              }
              sum(k,it) = sum(k,it) + 
                (block_x(it*dim+r,k)*(block_y(r, it)-mean_r)/(sigma_2r*(1-pow(rho_rt,2.0))) + 
                block_x(it*dim+t, k)*(block_y(t, it)-mean_t)/(sigma_2t*(1-pow(rho_rt,2.0))) - 
                rho_rt*(block_x(it*dim+t, k)*block_y(r, it)+block_x(it*dim+r,k)*block_y(t, it)-sum_1)/
                  ((1-pow(rho_rt,2.0))*pow(sigma_2r*sigma_2t,0.5))); 
            }
          }
        }
      }
    }
};

struct eenormalderivmean_worker : public Worker
{
  // source
  const int p;
  const RMatrix<double> cov; 
  const RMatrix<double> block_y;
  const RMatrix<double> block_x;
  const int dim;
  const int len;
  
  // destination 
  RMatrix<double> sum;
  
  // initialize with source and destination
  eenormalderivmean_worker(const int p, const NumericMatrix cov, const NumericMatrix block_y, const NumericMatrix block_x, 
                      const int dim, const int len, NumericMatrix sum) 
    : p(p), cov(cov), block_y(block_y), block_x(block_x), dim(dim), len(len), sum(sum) {}
  
  void operator()(std::size_t begin, std::size_t end) {

    double sigma_2r = 0;
    double sigma_2t = 0;
    double rho_rt = 0;
    int j = 0;
    
    for(unsigned int i=begin; i < end; i++) {
      int it = static_cast<int>(i);
      for(int r=0; r < (dim-1); r++) {
        for(int t=(r+1); t < dim; t++){
          sigma_2r = cov(r,r);
          sigma_2t = cov(t,t);
          rho_rt = cov(r,t)/pow(sigma_2r*sigma_2t, 0.5);
          
          for (int k=0; k < p; k++){
            for (int q=0; q < p; q++){
              j = it*p + q;
              sum(k,j) = sum(k,j) + (sigma_2t*block_x(it*dim+r,k)*block_x(it*dim+r,q) + sigma_2r*block_x(it*dim+t,k)*block_x(it*dim+t,q)-
                rho_rt*pow(sigma_2r*sigma_2t, 0.5)*(block_x(it*dim+r,k)*block_x(it*dim+t,q)+block_x(it*dim+t,k)*block_x(it*dim+r,q)))/
                (sigma_2r*sigma_2t*(1-pow(rho_rt,2.0)));
            }
          }
        }
      }
    }
  }
};

struct eenormalCSvar_worker : public Worker
{
  // source
  const RVector<double> beta;
  double sigma;
  double rho;
  const RMatrix<double> block_y;
  const RMatrix<double> block_x;
  const int dim;
  const int len;
  
  // destination 
  RMatrix<double> sum;
  
  // initialize with source and destination
  eenormalCSvar_worker(const NumericVector beta, const double sigma, const double rho, const NumericMatrix block_y, const NumericMatrix block_x, 
                      const int dim, const int len, NumericMatrix sum) 
    : beta(beta), sigma(sigma), rho(rho), block_y(block_y), block_x(block_x), dim(dim), len(len), sum(sum) {}
  
  void operator()(std::size_t begin, std::size_t end) {
    
    double mean_r = 0;
    double mean_t = 0;
    double sum_1 = 0;
    
    for(unsigned int i=begin; i < end; i++) {
      int it = static_cast<int>(i);
      for(int r=0; r < (dim-1); r++) {
        for(int t=(r+1); t < dim; t++){
          
          mean_r = 0;
          mean_t = 0;
          for (unsigned int k=0; k < beta.size(); k++){
            mean_r += block_x(it*dim+r, k)*beta[k];
            mean_t += block_x(it*dim+t, k)*beta[k];
          }
          
          for (unsigned int k=0; k < beta.size(); k++){
            sum_1 = 0;
            for (unsigned int q=0; q < beta.size(); q++){
              sum_1 += (block_x(it*dim+r, k)*block_x(it*dim+t, q) + block_x(it*dim+t, k)*block_x(it*dim+r, q))*beta[q];
            }
            
            sum(k,it) = sum(k, it) + (block_x(it*dim+r,k)*(block_y(r, i)-mean_r)/(pow(sigma,2.0)*(1-pow(rho,2.0))) + 
              block_x(it*dim+t,k)*(block_y(t, i)-mean_t)/(pow(sigma,2.0)*(1-pow(rho,2.0))) - 
              rho*(block_x(it*dim+t,k)*block_y(r, i)+block_x(it*dim+r,k)*block_y(t, i)-sum_1)/((1-pow(rho,2.0))*pow(sigma,2.0)));
          }
          
          sum(beta.size(),it) = sum(beta.size(),it) + 
            as_scalar(-2/sigma + (pow(block_y(r, i)-mean_r, 2.0)-2*rho*(block_y(r, i)-mean_r)*(block_y(t, i)-mean_t)+pow(block_y(t, i)-mean_t, 2.0))/
            ((1-pow(rho, 2.0))*pow(sigma,3.0)));
          
          sum(beta.size()+1, it) = sum(beta.size()+1, it) +
            as_scalar(rho/(1-pow(rho, 2.0))-(rho/pow(1-pow(rho,2.0),2.0))*(pow(block_y(r, i)-mean_r, 2.0)-2*rho*(block_y(r, i)-mean_r)*(block_y(t, i)-mean_t)+
            pow(block_y(t, i)-mean_t, 2.0))/pow(sigma, 2.0)+((block_y(r, i)-mean_r)*(block_y(t, i)-mean_t))/(pow(sigma,2.0)*(1-pow(rho,2.0))));
        
        }
      }
    }
  }
};

struct eenormalAR1var_worker : public Worker
{
  // source
  const RVector<double> beta;
  double sigma;
  double rho;
  const RMatrix<double> block_y;
  const RMatrix<double> block_x;
  const int dim;
  const int len;
  
  // destination 
  RMatrix<double> sum;
  
  // initialize with source and destination
  eenormalAR1var_worker(const NumericVector beta, const double sigma, const double rho, const NumericMatrix block_y, const NumericMatrix block_x, 
                       const int dim, const int len, NumericMatrix sum) 
    : beta(beta), sigma(sigma), rho(rho), block_y(block_y), block_x(block_x), dim(dim), len(len), sum(sum) {}
  
  void operator()(std::size_t begin, std::size_t end) {
    
    double mean_r = 0;
    double mean_t = 0;
    double sum_1 = 0;
    
    for(unsigned int i=begin; i < end; i++) {
      int it = static_cast<int>(i);
      for(int r=0; r < (dim-1); r++) {
        for(int t=(r+1); t < dim; t++){
          
          mean_r = 0;
          mean_t = 0;
          for (unsigned int k=0; k < beta.size(); k++){
            mean_r += block_x(it*dim+r, k)*beta[k];
            mean_t += block_x(it*dim+t, k)*beta[k];
          }
          
          for (unsigned int k=0; k < beta.size(); k++){
            sum_1 = 0;
            for (unsigned int q=0; q < beta.size(); q++){
              sum_1 += (block_x(it*dim+r, k)*block_x(it*dim+t, q) + block_x(it*dim+t, k)*block_x(it*dim+r, q))*beta[q];
            }
            sum(k,it) = sum(k,it) +
              (block_x(it*dim+r,k)*(block_y(r, i)-mean_r)/(pow(sigma,2.0)*(1-pow(rho,2.0*(t-r)))) + 
              block_x(it*dim+t,k)*(block_y(t, i)-mean_t)/(pow(sigma,2.0)*(1-pow(rho,2.0*(t-r)))) - 
              pow(rho,t-r)*(block_x(it*dim+t,k)*block_y(r, i)+block_x(it*dim+r,k)*block_y(t, i)-sum_1)/((1-pow(rho,2.0*(t-r)))*pow(sigma,2.0)));
          }
          
          sum(beta.size(),it) = sum(beta.size(),it) +
            as_scalar(-2/sigma + (pow(block_y(r, i)-mean_r, 2.0)-2*pow(rho,t-r)*(block_y(r, i)-mean_r)*(block_y(t, i)-mean_t)+pow(block_y(t, i)-mean_t, 2.0))/
            ((1-pow(rho, 2.0*(t-r)))*pow(sigma,3.0)));
          
          sum(beta.size()+1, it) = sum(beta.size()+1, it) +
            as_scalar((r-t)/(rho-pow(rho, 2.0*(r-t)+1)) - (t-r)*pow(rho, 2.0*(t-r)-1)*
            (pow(block_y(r, i)-mean_r, 2.0)-2*pow(rho, t-r)*(block_y(r, i)-mean_r)*(block_y(t, i)-mean_t)+pow(block_y(t, i)-mean_t,2.0))/
              (pow(1-pow(rho, 2.0*(t-r)), 2.0)*pow(sigma, 2.0)) +
                (t-r)*pow(rho, t-r-1)*(block_y(r, i)-mean_r)*(block_y(t, i)-mean_t)/
                  (pow(sigma, 2.0)*(1-pow(rho, 2.0*(t-r)))));
        }
      }
    }
  }
};

struct eenormalindvar_worker : public Worker
{
  // source
  const RVector<double> beta;
  double sigma;
  const RMatrix<double> block_y;
  const RMatrix<double> block_x;
  const int dim;
  const int len;
  
  // destination 
  RMatrix<double> sum;
  
  // initialize with source and destination
  eenormalindvar_worker(const NumericVector beta, const double sigma, const NumericMatrix block_y, const NumericMatrix block_x, 
                        const int dim, const int len, NumericMatrix sum) 
    : beta(beta), sigma(sigma), block_y(block_y), block_x(block_x), dim(dim), len(len), sum(sum) {}
  
  void operator()(std::size_t begin, std::size_t end) {

    double mean_r = 0;
    double mean_t = 0;
    
    for(unsigned int i=begin; i < end; i++) {
      int it = static_cast<int>(i);
      for(int r=0; r < (dim-1); r++) {
        for(int t=(r+1); t < dim; t++){
          
          mean_r = 0;
          mean_t = 0;
          for (unsigned int k=0; k < beta.size(); k++){
            mean_r += block_x(it*dim+r, k)*beta[k];
            mean_t += block_x(it*dim+t, k)*beta[k];
          }
          
          for (unsigned int k=0; k < beta.size(); k++){
            sum(k,it) = sum(k,it) +
              (block_x(it*dim+r,k)*(block_y(r, i)-mean_r) + block_x(it*dim+t,k)*(block_y(t, i)-mean_t))/pow(sigma, 2.0);
          }
          
          sum(beta.size(),it) = sum(beta.size(),it) - 2/sigma + (pow(block_y(r, i)-mean_r, 2.0) + pow(block_y(t, i)-mean_t, 2.0))/pow(sigma, 3.0);
        }
      }
    }
  }
};

// [[Rcpp::export]]
double logCLnormal(arma::vec beta, arma::mat cov, arma::mat block_y, arma::mat block_x, double m, double n)
{
  double sum = 0;
  int dim = m;
  int len = n;
  double sigma_2r = 0;
  double sigma_2t = 0;
  double rho_rt = 0;
  arma::vec xir_arma;
  arma::vec xit_arma;
  double mean_r = 0;
  double mean_t = 0;
    
  for(int i=0; i < len; i++) {
    for(int r=0; r < (dim-1); r++) {
      for(int t=r+1; t < dim; t++){
        sigma_2r = cov(r,r);
        sigma_2t = cov(t,t);
        rho_rt = cov(r,t)/pow(sigma_2r*sigma_2t, 0.5);
        
        xir_arma = block_x.row(i*dim+r).t();
        xit_arma = block_x.row(i*dim+t).t();
        mean_r = as_scalar(xir_arma.t()*beta);
        mean_t = as_scalar(xit_arma.t()*beta);
    
        sum += (-log(pow(sigma_2r*sigma_2t*(1-pow(rho_rt,2.0)), 0.5)) - 
        (pow(block_y(r, i)-mean_r, 2.0)/sigma_2r + pow(block_y(t, i)-mean_t, 2.0)/sigma_2t - 
        2*rho_rt*(block_y(r, i)-mean_r)*(block_y(t, i)-mean_t)/pow(sigma_2r*sigma_2t,0.5))/
        (2*(1-pow(rho_rt, 2.0))));
      }
    }
  }
  return sum;
}

// [[Rcpp::export]]
double logCLnormal_par(NumericVector beta, NumericMatrix cov, NumericMatrix block_y, NumericMatrix block_x, double m, double n)
{
  int dim = m;
  int len = n;
  
  // logCLnormal functor (pass input matrices)
  logCLnormal_worker little_worker(beta, cov, block_y, block_x, dim);
  
  // call parallelReduce to do the work
  parallelReduce(0, len, little_worker);
  
  // return the sum
  return little_worker.sum;
}

// [[Rcpp::export]]
List eenormalmean(arma::vec beta, arma::mat cov, arma::mat block_y, arma::mat block_x, double m, double n)
{
  int dim = m;
  int len = n;
  double sigma_2r = 0;
  double sigma_2t = 0;
  double rho_rt = 0;
  arma::vec xir_arma;
  arma::vec xit_arma;
  double mean_r = 0;
  double mean_t = 0;
  List sum(len);
  
  for(int i=0; i < len; i++) {
    arma::vec sum_1 = zeros<vec>(beta.size());
    for(int r=0; r < (dim-1); r++) {
      for(int t=(r+1); t < dim; t++){
        sigma_2r = cov(r,r);
        sigma_2t = cov(t,t);
        rho_rt = cov(r,t)/pow(sigma_2r*sigma_2t, 0.5);
        
        xir_arma = block_x.row(i*dim+r).t();
        xit_arma = block_x.row(i*dim+t).t();
        mean_r = as_scalar(xir_arma.t()*beta);
        mean_t = as_scalar(xit_arma.t()*beta);
        
        sum_1 += (xir_arma*(block_y(r, i)-mean_r)/(sigma_2r*(1-pow(rho_rt,2.0))) + 
          xit_arma*(block_y(t, i)-mean_t)/(sigma_2t*(1-pow(rho_rt,2.0))) - 
          rho_rt*(xit_arma*block_y(r, i)+xir_arma*block_y(t, i)-
          (xir_arma*xit_arma.t()+xit_arma*xir_arma.t())*beta)/((1-pow(rho_rt,2.0))*pow(sigma_2r*sigma_2t,0.5)));
      }
    }
    sum[i] = sum_1;
  }
  return sum;
}

// [[Rcpp::export]]
List eenormalmean_par(NumericVector beta, NumericMatrix cov, NumericMatrix block_y, NumericMatrix block_x, double m, double n)
{
  int dim = m;
  int len = n;
  NumericMatrix sum (beta.length(), len);
  
  
  // eenormalmean_par functor (pass input and output matrices)
  eenormalmean_worker little_worker(beta, cov, block_y, block_x, dim, len, sum);
  
  // call parallelFor to do the work
  parallelFor(0, len, little_worker); 
  
  // return the output matrix
  List output(len);
  for (int i=0; i<len; i++){
    arma::mat MAT(sum.column(i).begin(), beta.size(), 1, false);
    output[i] = MAT;
  }
  return output;
}

// [[Rcpp::export]]
List eenormalderivmean(arma::mat cov, arma::mat block_y, arma::mat block_x, double m, double n)
{
  int dim = m;
  int len = n;
  double sigma_2r = 0;
  double sigma_2t = 0;
  double rho_rt = 0;
  arma::vec xir_arma;
  arma::vec xit_arma;
  List sum(len);
  
  for(int i=0; i < len; i++) {
    arma::mat sum_1 = zeros<mat>(block_x.n_cols, block_x.n_cols);
    for(int r=0; r < (dim-1); r++) {
      for(int t=(r+1); t < dim; t++){
        sigma_2r = cov(r,r);
        sigma_2t = cov(t,t);
        rho_rt = cov(r,t)/pow(sigma_2r*sigma_2t, 0.5);
        
        xir_arma = block_x.row(i*dim+r).t();
        xit_arma = block_x.row(i*dim+t).t();
        
        sum_1 += (sigma_2r*xir_arma*xir_arma.t()+sigma_2t*xit_arma*xit_arma.t() - 
          rho_rt*pow(sigma_2r*sigma_2t, 0.5)*(xir_arma*xit_arma.t()+xit_arma*xir_arma.t()))/
            (sigma_2r*sigma_2t*(1-pow(rho_rt,2.0)));
      }
    }
    sum[i] = sum_1;
  }
  return sum;
}

// [[Rcpp::export]]
List eenormalderivmean_par(NumericMatrix cov, NumericMatrix block_y, NumericMatrix block_x, double m, double n)
{
  int dim = m;
  int len = n;
  int p = block_x.ncol();
  NumericMatrix sum (p, p*len);
  
  
  // eenormalmean_par functor (pass input and output matrices)
  eenormalderivmean_worker little_worker(p, cov, block_y, block_x, dim, len, sum);
  
  // call parallelFor to do the work
  parallelFor(0, len, little_worker); 
  
  // return the output matrix
  List output(len);
  for (int i=0; i<len; i++){
      arma::mat MAT(sum.column(i*p).begin(), p, p, false);
      output[i] = MAT;
  }
  return output;
}

// [[Rcpp::export]]
List eenormalCSvar(arma::vec beta, double sigma, double rho, arma::mat block_y, arma::mat block_x, double m, double n)
{
  int dim = m;
  int len = n;
  
  arma::vec xir_arma;
  arma::vec xit_arma;
  double mean_r = 0;
  double mean_t = 0;
  List sum(len);
  
  for(int i=0; i < len; i++) {
    arma::mat sum_1 = zeros<mat>(beta.size(),1);
    arma::mat sum_2 = zeros<mat>(1,1);
    arma::mat sum_3 = zeros<mat>(1,1);
    for(int r=0; r < (dim-1); r++) {
      for(int t=(r+1); t < dim; t++){
        
        xir_arma = block_x.row(i*dim+r).t();
        xit_arma = block_x.row(i*dim+t).t();
        mean_r = as_scalar(xir_arma.t()*beta);
        mean_t = as_scalar(xit_arma.t()*beta);
        
        sum_1 += (xir_arma*(block_y(r, i)-mean_r)/(pow(sigma,2.0)*(1-pow(rho,2.0))) + 
          xit_arma*(block_y(t, i)-mean_t)/(pow(sigma,2.0)*(1-pow(rho,2.0))) - 
          rho*(xit_arma*block_y(r, i)+xir_arma*block_y(t, i)-
          (xir_arma*xit_arma.t()+xit_arma*xir_arma.t())*beta)/((1-pow(rho,2.0))*pow(sigma,2.0)));
        
        sum_2 += -2/sigma + (pow(block_y(r, i)-mean_r, 2.0)-2*rho*(block_y(r, i)-mean_r)*(block_y(t, i)-mean_t)+pow(block_y(t, i)-mean_t, 2.0))/
          ((1-pow(rho, 2.0))*pow(sigma,3.0));
        
        sum_3 += rho/(1-pow(rho, 2.0))-(rho/pow(1-pow(rho,2.0),2.0))*(pow(block_y(r, i)-mean_r, 2.0)-2*rho*(block_y(r, i)-mean_r)*(block_y(t, i)-mean_t)+
          pow(block_y(t, i)-mean_t, 2.0))/pow(sigma, 2.0)+((block_y(r, i)-mean_r)*(block_y(t, i)-mean_t))/(pow(sigma,2.0)*(1-pow(rho,2.0)));
      }
    }
    sum[i] = join_cols(join_cols(sum_1, sum_2), sum_3);
  }
  return sum;
}

// [[Rcpp::export]]
List eenormalCSvar_par(NumericVector beta, double sigma, double rho, NumericMatrix block_y, NumericMatrix block_x, double m, double n)
{
  int dim = m;
  int len = n;
  int q = beta.size()+2;
  NumericMatrix sum (q, len);
  
  
  // eenormalmean_par functor (pass input and output matrices)
  eenormalCSvar_worker little_worker(beta, sigma, rho, block_y, block_x, dim, len, sum);
  
  // call parallelFor to do the work
  parallelFor(0, len, little_worker); 
  
  // return the output matrix
  List output(len);
  for (int i=0; i<len; i++){
    arma::mat MAT(sum.column(i).begin(), q, 1, false);
    output[i] = MAT;
  }
  return output;
}

// [[Rcpp::export]]
List eenormalAR1var(arma::vec beta, double sigma, double rho, arma::mat block_y, arma::mat block_x, double m, double n)
{
  int dim = m;
  int len = n;
  
  arma::vec xir_arma;
  arma::vec xit_arma;
  double mean_r = 0;
  double mean_t = 0;
  List sum(len);
  
  for(int i=0; i < len; i++) {
    arma::mat sum_1 = zeros<mat>(beta.size(),1);
    arma::mat sum_2 = zeros<mat>(1,1);
    arma::mat sum_3 = zeros<mat>(1,1);
    for(int r=0; r < (dim-1); r++) {
      for(int t=(r+1); t < dim; t++){
        
        xir_arma = block_x.row(i*dim+r).t();
        xit_arma = block_x.row(i*dim+t).t();
        mean_r = as_scalar(xir_arma.t()*beta);
        mean_t = as_scalar(xit_arma.t()*beta);
        
        sum_1 += (xir_arma*(block_y(r, i)-mean_r)/(pow(sigma,2.0)*(1-pow(rho,2.0*(t-r)))) + 
          xit_arma*(block_y(t, i)-mean_t)/(pow(sigma,2.0)*(1-pow(rho,2.0*(t-r)))) - 
          pow(rho,t-r)*(xit_arma*block_y(r, i)+xir_arma*block_y(t, i)-
          (xir_arma*xit_arma.t()+xit_arma*xir_arma.t())*beta)/((1-pow(rho,2.0*(t-r)))*pow(sigma,2.0)));
        
        sum_2 += -2/sigma + (pow(block_y(r, i)-mean_r, 2.0)-2*pow(rho,t-r)*(block_y(r, i)-mean_r)*(block_y(t, i)-mean_t)+pow(block_y(t, i)-mean_t, 2.0))/
          ((1-pow(rho, 2.0*(t-r)))*pow(sigma,3.0));
        
        sum_3 += (r-t)/(rho-pow(rho, 2.0*(r-t)+1)) - (t-r)*pow(rho, 2.0*(t-r)-1)*
          (pow(block_y(r, i)-mean_r, 2.0)-2*pow(rho, t-r)*(block_y(r, i)-mean_r)*(block_y(t, i)-mean_t)+pow(block_y(t, i)-mean_t,2.0))/
            (pow(1-pow(rho, 2.0*(t-r)), 2.0)*pow(sigma, 2.0)) +
              (t-r)*pow(rho, t-r-1)*(block_y(r, i)-mean_r)*(block_y(t, i)-mean_t)/
                (pow(sigma, 2.0)*(1-pow(rho, 2.0*(t-r))));
      }
    }
    sum[i] = join_cols(join_cols(sum_1, sum_2), sum_3);
  }
  return sum;
}

// [[Rcpp::export]]
List eenormalAR1var_par(NumericVector beta, double sigma, double rho, NumericMatrix block_y, NumericMatrix block_x, double m, double n)
{
  int dim = m;
  int len = n;
  int q = beta.size()+2;
  NumericMatrix sum (q, len);
  
  
  // eenormalmean_par functor (pass input and output matrices)
  eenormalAR1var_worker little_worker(beta, sigma, rho, block_y, block_x, dim, len, sum);
  
  // call parallelFor to do the work
  parallelFor(0, len, little_worker); 
  
  // return the output matrix
  List output(len);
  for (int i=0; i<len; i++){
    arma::mat MAT(sum.column(i).begin(), q, 1, false);
    output[i] = MAT;
  }
  return output;
}

// [[Rcpp::export]]
List eenormalindvar(arma::vec beta, double sigma, arma::mat block_y, arma::mat block_x, double m, double n)
{
  int dim = m;
  int len = n;
  
  arma::vec xir_arma;
  arma::vec xit_arma;
  double mean_r = 0;
  double mean_t = 0;
  List sum(len);
  
  for(int i=0; i < len; i++) {
    arma::mat sum_1 = zeros<mat>(beta.size(),1);
    arma::mat sum_2 = zeros<mat>(1,1);
    for(int r=0; r < (dim-1); r++) {
      for(int t=(r+1); t < dim; t++){
        
        xir_arma = block_x.row(i*dim+r).t();
        xit_arma = block_x.row(i*dim+t).t();
        mean_r = as_scalar(xir_arma.t()*beta);
        mean_t = as_scalar(xit_arma.t()*beta);
        
        sum_1 += (xir_arma*(block_y(r, i)-mean_r) + xit_arma*(block_y(t, i)-mean_t))/pow(sigma, 2.0);
        
        sum_2 += -2/sigma + (pow(block_y(r, i)-mean_r, 2.0) + pow(block_y(t, i)-mean_t, 2.0))/pow(sigma, 3.0);
        
      }
    }
    sum[i] = join_cols(sum_1, sum_2);
  }
  return sum;
}

// [[Rcpp::export]]
List eenormalindvar_par(NumericVector beta, double sigma, NumericMatrix block_y, NumericMatrix block_x, double m, double n)
{
  int dim = m;
  int len = n;
  int q = beta.size()+1;
  NumericMatrix sum (q, len);
  
  
  // eenormalmean_par functor (pass input and output matrices)
  eenormalindvar_worker little_worker(beta, sigma, block_y, block_x, dim, len, sum);
  
  // call parallelFor to do the work
  parallelFor(0, len, little_worker); 
  
  // return the output matrix
  List output(len);
  for (int i=0; i<len; i++){
    arma::mat MAT(sum.column(i).begin(), q, 1, false);
    output[i] = MAT;
  }
  return output;
}