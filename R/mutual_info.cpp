#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <Rcpp.h>
#include <algorithm>
#include <vector>
#include <iterator>
#include <cmath>
using namespace Rcpp;

//#define USEDEBUG

#ifdef USEDEBUG
#define Debug(x) std::cout << x
#else
#define Debug(x)
#endif


struct MI_Data {
  double MI;
  int totalResponseBins;
};


class BinnedResponseModel {
public:
  NumericMatrix prob_r_given_s;
  NumericVector prob_response;
  NumericVector prob_stim;
};


/** 
 * Mutual information calculated for a single stimulus bin.
 * 
 * \param stimCount count of times the stimulus was observed: P(s) * N
 * \param N total number of stimuli
 * \param r_count vector with the counts each response was observed across stimuli: P(r) * N
 * \param r_given_s_count vector with the counts each response was observed for the stimuli: P(r|s)
 */
MI_Data stimulus_mutual_info(double p_stim, 
                             NumericVector& p_response,
                             const NumericMatrix::Row& p_r_given_s) {
  int totalResponseBins = 0;
  double MI = 0.0;
  for (int r = 0; r < p_response.size(); ++r) {
    
    Debug("response=" << r);
    Debug(", p_stim=" << p_stim);
    Debug(", p_r_given_s=" << p_r_given_s[r]);
    Debug(", p_response=" << p_response[r]);
    if (p_r_given_s[r] > 0) {
      ++totalResponseBins;
      double r_MI = p_stim * p_r_given_s[r] * std::log2(p_r_given_s[r] / p_response[r]);
      Debug(", MI=" << r_MI);
      MI += r_MI;
    }
    Debug(std::endl);
  }
  MI_Data result;
  result.MI = MI;
  result.totalResponseBins = totalResponseBins;
  return result;
}

/**
 * Uses Bayesian formula to decode stimulus (s) given the response (r) from the model of responses given the stimulus.
 * 
 * \param r observed response for which the max probability stimulus will be returned
 * \param N total number of stimuli
 * \param r_counts vector with the counts each response was observed across stimuli: P(r) * N
 * \param s_counts vector with the counts each stimulus was observed: P(s) * N
 * \param r_given_s_count matrix indexed by [s,r] with the counts each response was observed for a stimulus: P(r|s)
 */
/*
int max_s_given_r(int r, int N,
                  const IntegerVector& r_counts,
                  const IntegerVector& s_counts,
                  const IntegerMatrix& r_given_s_count) {
  
  int max_s_i = -1;
  double max_s = -1.0;
  
  for (int s = 0; s < s_counts.size(); ++s) {
    double p_r_given_s = ((double) r_given_s_count(s,r)) / s_counts[s];
    
    Debug("response=" << r);
    Debug(", p_stim=" << p_stim);
    Debug(", p_r_given_s=" << p_r_given_s);
    Debug(", p_response=" << p_response);
    double p_s_given_r = p_r_given_s * p_stim / p_response;
    Debug(", p_s_given_r=" << p_s_given_r);
    if (p_s_given_r > max_s) {
      max_s = p_s_given_r;
      max_s_i = s;
    }
    Debug(std::endl);
  }
  return max_s_i;
}
*/


BinnedResponseModel createResponseModel(IntegerVector& response, 
                                        int nresponseBins,
                                        IntegerVector& stimulus, 
                                        int nstim) {

  arma::mat r_given_s(nstim, nresponseBins, arma::fill::zeros);
  std::vector<int> r_counts(nresponseBins, 0);
  std::vector<int> s_counts(nstim, 0);

  for (int i = 0; i < response.size(); ++i) {
    int responseBin = response[i] - 1;
    int stimBin = stimulus[i] - 1;
    if (stimBin > nstim) {
      std::cout << "Error: Stimulus id greater than the passed count of stimuli (nstim)"  << std::endl;
      continue;
    }
    
    ++s_counts[stimBin];
    ++r_given_s(stimBin, responseBin);
    ++r_counts[responseBin];
  }  

  int N = response.size(); 
  
  // Set probabilities
  BinnedResponseModel m = BinnedResponseModel();
  NumericVector p_stim(nstim, 0.0);
  NumericVector p_response(nresponseBins, 0.0);
  NumericMatrix p_r_given_s(nstim, nresponseBins);
  std::fill(p_r_given_s.begin(), p_r_given_s.end(), 0.0);
  
  for (int r = 0; r < nresponseBins; ++r) {
    p_response[r] = ((double) r_counts[r]) / N;
    for (int s = 0; s < nstim; ++s) {
      p_stim[s] = ((double) s_counts[s]) / N;
      // +1 every bin to have >0 probabilities of unssen bins
      p_r_given_s(s,r) = ((double) r_given_s(s,r) + 1) / (s_counts[s] + nresponseBins);
    }
  }
  
  m.prob_stim = p_stim;
  m.prob_response = p_response;
  m.prob_r_given_s = p_r_given_s;
  return m;
}


IntegerVector chunkShuffle(IntegerVector& trace,
                           IntegerVector& trialEnds, 
                           int shuffleChunkLength) {
  
  // Shuffle the trace keeping the order within small chunks
  int nchunks = std::ceil(((double) trace.size() / shuffleChunkLength));
  std::vector<int> chunkShuffle(nchunks);
  for (int j = 0; j < nchunks; ++j) {
    chunkShuffle[j] = j;
  }
  
  IntegerVector shuffledTrace(trace.size(), 0);
  
  // Shuffle chunks only within the same trial
  int trialStart = 0;
  for (int trial_i = 0; trial_i < trialEnds.size(); ++trial_i) {
    int trialEnd = trialEnds[trial_i] / shuffleChunkLength;
    std::random_shuffle(chunkShuffle.begin() + trialStart, chunkShuffle.begin() + trialEnd);
    trialStart = trialEnd;
  }
  
  // Offset from the start to randomize chunk boundaries
  int offset = std::rand() % ((int) shuffleChunkLength/2);
  for (int chunk = 0; chunk < nchunks; ++chunk) {
    int src_chunk_start = chunk * shuffleChunkLength + offset;
    int target_chunk_start = chunkShuffle[chunk] * shuffleChunkLength + offset;
    for (int j = 0; j < shuffleChunkLength; ++j) {
      int src_index = (src_chunk_start + j) % trace.size();
      int target_index = (target_chunk_start + j) % trace.size();
      shuffledTrace[target_index] = trace[src_index];
    }
  }
  
  return shuffledTrace;
}


// [[Rcpp::export]]
SEXP mutual_info(IntegerVector& response,
                 int nresponseBins,
                 IntegerVector& stimulus,
                 int nstim) {
  List result;
  result["mutual.info"] = 0.0;
  result["mutual.info.bias"] = 0.0;
  
  if (response.size() != stimulus.size()) {
    std::cout << "Error: response size need to equal stimulus size vector"  << std::endl;
    return result;
  }
  
  BinnedResponseModel m = createResponseModel(response, nresponseBins, stimulus, nstim);
  
  int totalResponseBins = 0;
  int totalStimuliBins = 0;
  double MI = 0.0;
  for (int s = 0; s < nstim; ++s) {
    double p_stim = m.prob_stim[s];
    if (p_stim > 0) {
      ++totalStimuliBins;
      
      MI_Data mi_Data = stimulus_mutual_info(p_stim, 
                                             m.prob_response, 
                                             m.prob_r_given_s.row(s));
      MI += mi_Data.MI;
      totalResponseBins += mi_Data.totalResponseBins;
    }
  }

  Debug("totalresponseBins=" << totalResponseBins << ", nresponseBins=" << nresponseBins << ", totalStimuliBins=" << totalStimuliBins << std::endl);
  double MI_bias = (totalResponseBins - nresponseBins - totalStimuliBins + 1 )  / (2 * response.size() * std::log(2));
  
  result["mutual.info"] = MI;
  result["mutual.info.bias"] = MI_bias;
  return result;
}


// [[Rcpp::export]]
SEXP mutual_info_with_shuffles(IntegerVector& response, 
                               int nresponseBins,
                               IntegerVector& stimulus,
                               int nstim,
                               IntegerVector& trialEnds,
                               int nshuffles,
                               int shuffleChunkLength ) {

  List result = mutual_info(response, nresponseBins, stimulus, nstim);
  
  NumericVector mis(nshuffles, 0.0);
  NumericVector mis_bias(nshuffles, 0.0);
  
  for (int i = 0; i < nshuffles; ++i) {
    IntegerVector shuffledResponse = chunkShuffle(response, trialEnds, shuffleChunkLength);
    Debug("Shuffled response: " << shuffledResponse << std::endl);
    List shuffledResult = mutual_info(shuffledResponse, nresponseBins, stimulus, nstim);
    mis[i] = shuffledResult["mutual.info"];
    mis_bias[i] = shuffledResult["mutual.info.bias"];
  }
  
  result["shuffle.mi"] = mis;
  result["shuffle.mi.bias"] = mis_bias;
  
  return result;
}

/*** R
# Perfect mutual info
response = c(c(1, 2, 1, 2, 2), rep(1, 10))
stimulus = c(c(1, 2, 1, 2, 2), rep(1, 10))
res = mutual_info_with_shuffles(response, 2, stimulus, 2, c(1, 10), 1, 2)
res$mutual.info
res$mutual.info.bias
*/

