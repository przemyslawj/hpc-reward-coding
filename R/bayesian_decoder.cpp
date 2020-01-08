#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <Rcpp.h>
#include <algorithm>
#include <vector>
#include <iostream>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::export]]
SEXP bayesmax(NumericVector& prior,
              NumericVector& likelihood,
              NumericVector& pv) {
  if (pv.size() == 0) {
    stop("Error, pv empty");
  }
  
  CharacterVector pvCellNames = pv.attr("names");
  List dimnames(likelihood.attr("dimnames"));
  CharacterVector modelCells = dimnames[1];
  
  IntegerVector likelihoodDim = likelihood.attr("dim");
  arma::cube likelihoodM(likelihood.begin(), likelihoodDim[0], likelihoodDim[1], likelihoodDim[2], false);
  
  std::unordered_map<int, int> pv2model_name(pvCellNames.size());
  for (int i = 0; i < pvCellNames.size(); ++i) {
    auto it = std::find(modelCells.begin(), modelCells.end(), pvCellNames[i]);
    if (it != modelCells.end()) {
      pv2model_name[i] = std::distance(modelCells.begin(), it);
      //std::cout << "Cell  " << pvCellNames[i] << " Found at " << pv2model_name[i]  << std::endl;    
    }
  }
  
  int nstim = likelihoodDim[0];
  
  double max_s = -1;
  double max_s_prob = -1.0;
  NumericVector probDensity(nstim);
  // Bayes rule (will assume P(s) uniform):
  // max_s [P(s|r)] ~ max_s P(r|s) * P(s)
  // Assuming conditionally independent activity of the cells:
  // P(r|s) = II_i P(r_i|s)
  for (int s = 0; s < nstim; ++s) {
    double s_prob = prior[s];;
    if (s_prob != NA_REAL) {
      for (int ci = 0; ci < pv2model_name.size(); ++ci) {
        int model_cell_i = pv2model_name[ci];
        double prob = likelihoodM(s, model_cell_i, pv[ci] - 1);
        //std::cout<<"Prob for s=" << s+1 << " cell=" << pvCellNames[ci] << " response=" << pv[ci] << " is=" << prob << std::endl;
        if (prob != NA_REAL) {
          s_prob = s_prob * prob;
        }
      }
      
      probDensity[s] = s_prob;
      if (max_s_prob <= s_prob) {
        max_s_prob = s_prob;
        max_s = s + 1;
      }
    } else{
      probDensity[s] = 0.0;
    }
  }
  
  List res = List::create(Named("s")=max_s, Named("prob")=max_s_prob);
  res["density"] = probDensity;
  return(res);
}

// [[Rcpp::export]]
double poisson_prob(double lambda, double k) {
  double prob = std::pow(lambda, k) * std::exp(-lambda) / std::tgamma(k + 1);
  return prob;
}

// [[Rcpp::export]]
SEXP bayesmax_mfr(NumericVector& prior,
                  NumericVector& mfr,
                  NumericVector& pv) {
  if (pv.size() == 0) {
    stop("Error, pv empty");
  }
  
  CharacterVector pvCellNames = pv.attr("names");
  List dimnames(mfr.attr("dimnames"));
  CharacterVector modelCells = dimnames[0];
  
  IntegerVector mfrDim = mfr.attr("dim");
  arma::mat frM(mfr.begin(), mfrDim[0], mfrDim[1], false);
  
  std::unordered_map<int, int> pv2model_name(pvCellNames.size());
  for (int i = 0; i < pvCellNames.size(); ++i) {
    auto it = std::find(modelCells.begin(), modelCells.end(), pvCellNames[i]);
    if (it != modelCells.end()) {
      pv2model_name[i] = std::distance(modelCells.begin(), it);
    }
  }
  
  int nstim = mfrDim[1];
  
  double max_s = -1;
  double max_s_prob = -1.0;
  NumericVector probDensity(nstim);
  // Bayes rule (will assume P(s) uniform):
  // max_s [P(s|r)] ~ max_s P(r|s) * P(s)
  // Assuming conditionally independent activity of the cells:
  // P(r|s) = II_i P(r_i|s)
  for (int s = 0; s < nstim; ++s) {
    double s_prob = prior[s];;
    if (s_prob != NA_REAL) {
      for (int ci = 0; ci < pv2model_name.size(); ++ci) {
        int model_cell_i = pv2model_name[ci];
        double prob = poisson_prob(frM(model_cell_i, s), pv[ci]);
        //std::cout<<"Prob for s=" << s+1 << " cell=" << pvCellNames[ci] << " mfr=" << frM(model_cell_i, s) <<" response=" << pv[ci] << " is=" << prob << std::endl;
        if (prob != NA_REAL) {
          s_prob = s_prob * prob;
        }
      }
      
      probDensity[s] = s_prob;
      if (max_s_prob <= s_prob) {
        max_s_prob = s_prob;
        max_s = s + 1;
      }
    } else{
      probDensity[s] = 0.0;
    }
  }
  
  List res = List::create(Named("s")=max_s, Named("prob")=max_s_prob);
  res["density"] = probDensity;
  return(res);
}
/*
// [[Rcpp::export]]
SEXP createBayesModel(int nstim, 
                      int nresponse,
                      IntegerMatrix& cellResponse,
                      IntegerVector& stimulus) {
  
  const int ncells = cellResponse.rows();
  const int N = stimulus.size();
  if(N != cellResponse.cols()) {
    stop("Expected equal lengths of cellResponse and stimulus");
  }
  
  IntegerVector stim_count(nstim, 0);
  IntegerMatrix cell_s_count(ncells, nstim);
  std::fill(cell_s_count.begin(), cell_s_count.end(), 0);
  //arma::Cube<int> cell_response_count(ncells, nstim, nresponse, arma::fill::zeros);
  // P(r_i | s)
  
  int totalStim = 0;
  for (int i = 0; i < N; ++i) {
    int stimBin = stimulus[i] - 1;
    if (stimBin >= nstim) {
      std::cout << "Error: Stimulus id greater than the passed count of stimuli (nstim)"  << std::endl;
      return List();
    }
    ++stim_count(stimBin);
    ++totalStim;
    
    for (int cell_id; cell_id < ncells; ++cell_id) {
      const int response = cellResponse(cell_id, i);
      //++cell_response_count(cell_id, stimBin, response);
      ++cell_s_count(cell_id, stimBin);
    }  
    std::cout << "Stim #" << i << std::endl;
  }
  std::cout << "Finished counts" << std::endl;
  
  
  NumericVector prior(nstim, 0.0);
  for (int s = 0; s < nstim; ++s) {
    prior[s] = (double) stim_count[s] / totalStim;
  }
  //arma::cube likelihood(ncells, nstim, nresponse, arma::fill::zeros);
  /*
  for (int r = 0; r < nresponse; ++r) {  
    for (int cell_id; cell_id < ncells; ++cell_id) {
      for (int s = 0; s < nstim; ++s) {
         // TODO add some non-zero likelihood for each response bin if visited.
         //likelihood(cell_id, s, r) = (double) cell_response_count(cell_id, s, r) / cell_s_count(cell_id, s);
      }
    }
  }

  return List::create(Named("prior") = prior);
                      //Named("likelihood") = wrap(likelihood));
}
*/
/*
df = data.frame(
  animal='A',
  date='2019-01-01',
  cell_id=rep(c(1,2),each=4),
  trial_id='trial1',
  time_bin=1:4,
  mean.trace=c(1:4, 9:6),
  mean.x=1:4,
  mean.y=1
)

responseM = df %>%
  select(cell_id, time_bin, mean.trace) %>%
  acast(cell_id ~ time_bin, value.var='mean.trace')
model = createBayesModel(4, 2, responseM, df$mean.x[1:4])

*/
/*** R
likelihoodM = array(0.1, dim=c(4,2,2))
colnames(likelihoodM) = c('1', '2')
likelihoodM[3,'1',1] = 0.9
likelihoodM[3,'2',2] = 0.9

prior=rep(0.25, 4)
pv = c(`1`=1,`2`=2)
res = bayesmax(prior, likelihoodM, pv)
expect_equal(3, res$s)



# bayesmax_mfr
mfr.matrix = matrix(0.01, nrow=2, ncol=4)
rownames(mfr.matrix) = c('1', '2')
mfr.matrix['1', 1] = 3
mfr.matrix['2', 2] = 5

prior=rep(0.25, 4)
pv = c(`1`=1,`2`=5)
res = bayesmax_mfr(prior, mfr.matrix, pv)
expect_equal(2, res$s)

*/

