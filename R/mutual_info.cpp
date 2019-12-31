#include <Rcpp.h>
#include <algorithm>
#include <vector>
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


class ResponseModel {
public:
  NumericVector field;
  IntegerMatrix r_given_s;
  IntegerVector r_counts;
  IntegerVector s_counts;
  
  ResponseModel(int nstim, int nresponseBins) {
    this->field = NumericVector(nstim, 0.0);
    this->r_given_s = IntegerMatrix(nstim, nresponseBins);
    this->r_counts = IntegerVector(nresponseBins, 0);
    this->s_counts = IntegerVector(nstim, 0);
  }
};


/** 
 * Mutual information calculated for a single stimulus bin.
 * 
 * \param stimCount count of times the stimulus was observed: P(s) * N
 * \param N total number of stimuli
 * \param r_count vector with the counts each response was observed across stimuli: P(r) * N
 * \param r_given_s_count vector with the counts each response was observed for the stimuli: P(r|s)
 */
MI_Data stimulus_mutual_info(int stimCount, int N, 
                             const IntegerVector& r_count,
                             const IntegerMatrix::Row& r_given_s_count) {
  int totalResponseBins = 0;
  double MI = 0.0;
  double p_stim = (double) stimCount / N;
  for (int r = 0; r < r_given_s_count.size(); ++r) {
    double p_response = ((double) r_count[r]) / N;
    double p_r_given_s = ((double) r_given_s_count[r]) / stimCount;
    
    Debug("response=" << r);
    Debug(", p_stim=" << p_stim);
    Debug(", p_r_given_s=" << p_r_given_s);
    Debug(", p_response=" << p_response);
    if (p_r_given_s > 0) {
      ++totalResponseBins;
      double r_MI = p_stim * p_r_given_s * std::log2(p_r_given_s / p_response);
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
int max_s_given_r(int r, int N,
                  const IntegerVector& r_counts,
                  const IntegerVector& s_counts,
                  const IntegerMatrix& r_given_s_count) {
  
  int max_s_i = -1;
  double max_s = -1.0;
  
  double p_response = ((double) r_counts[r]) / N;
  for (int s = 0; s < s_counts.size(); ++s) {
    double p_stim = (double) s_counts[s] / N;
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

ResponseModel createResponseModel(int nstim,
                                  NumericVector& response, 
                                  NumericVector& stimulus, 
                                  NumericVector& responseQuantiles) {
  int nresponseBins = responseQuantiles.size();
  ResponseModel m = ResponseModel(nstim, nresponseBins);

  for (int i = 0; i < response.size(); ++i) {
    auto responseBinIt = std::lower_bound(responseQuantiles.begin(), responseQuantiles.end(), response[i]);
    int responseBin = nresponseBins - 1;
    if (responseBinIt != responseQuantiles.end()) {
      responseBin = responseBinIt - responseQuantiles.begin();
    }
    
    int stimBin = stimulus[i];
    if (stimBin > nstim) {
      std::cout << "Error: Stimulus id greater than the passed count of stimuli (nstim)"  << std::endl;
      return m;
    }
    
    m.field(stimBin) += response[i];
    ++m.s_counts[stimBin];
    ++m.r_given_s(stimBin, responseBin);
    ++m.r_counts[responseBin];
  }  
  
  for (int s = 0; s < nstim; ++s) {
    m.field(s) = m.field(s) / m.s_counts[s];
  }
  
  return m;
}

// [[Rcpp::export]]
SEXP createBayesModel2(int nstim, 
                      int nresponse,
                      IntegerMatrix& cellResponse,
                      IntegerVector& stimulus) {
  
  return List();
}


NumericVector chunkShuffle(NumericVector& trace,
                           NumericVector& trialEnds, 
                           int shuffleChunkLength) {
  
  // Shuffle the trace keeping the order within small chunks
  int nchunks = std::ceil(((double) trace.size() / shuffleChunkLength));
  std::vector<int> chunkShuffle(nchunks);
  for (int j = 0; j < nchunks; ++j) {
    chunkShuffle[j] = j;
  }
  
  NumericVector shuffledTrace(trace.size(), 0);
  
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
SEXP mutual_info(NumericVector& response, 
                 NumericVector& stimulus, 
                 NumericVector& responseQuantiles,
                 int nstim) {
  List result;
  result["mutual.info"] = 0.0;
  result["mutual.info.bias"] = 0.0;
  
  int nresponseBins = responseQuantiles.size();
  
  if (response.size() != stimulus.size()) {
    std::cout << "Error: response size need to equal stimulus size vector"  << std::endl;
    return result;
  }
  
  ResponseModel m = createResponseModel(nstim, response, stimulus, responseQuantiles);
  
  int totalResponseBins = 0;
  int totalStimuliBins = 0;
  double MI = 0.0;
  for (int s = 0; s < nstim; ++s) {
    int stimCount = m.s_counts[s];
    if (stimCount > 0) {
      ++totalStimuliBins;
      
      MI_Data mi_Data = stimulus_mutual_info(stimCount, 
                                             response.size(), 
                                             m.r_counts, 
                                             m.r_given_s.row(s));
      MI += mi_Data.MI;
      totalResponseBins += mi_Data.totalResponseBins;
    }
  }

  Debug("totalresponseBins=" << totalResponseBins << ", nresponseBins=" << nresponseBins << ", totalStimuliBins=" << totalStimuliBins << std::endl);
  double MI_bias = (totalResponseBins - nresponseBins - totalStimuliBins + 1 )  / (2 * response.size() * std::log(2));
  
  result["mutual.info"] = MI;
  result["mutual.info.bias"] = MI_bias;
  result["field"] = m.field;
  return result;
}


// [[Rcpp::export]]
SEXP mutual_info_with_shuffles(NumericVector& response, 
                               NumericVector& stimulus, 
                               NumericVector& responseQuantiles,
                               NumericVector& trialEnds,
                               int nstim,
                               int nshuffles,
                               int shuffleChunkLength ) {

  List result = mutual_info(response, stimulus, responseQuantiles, nstim);
  
  NumericVector mis(nshuffles, 0.0);
  NumericVector mis_bias(nshuffles, 0.0);
  
  for (int i = 0; i < nshuffles; ++i) {
    NumericVector shuffledResponse = chunkShuffle(response, trialEnds, shuffleChunkLength);
    Debug("Shuffled response: " << shuffledResponse << std::endl);
    List shuffledResult = mutual_info(shuffledResponse, stimulus, responseQuantiles, nstim);
    mis[i] = shuffledResult["mutual.info"];
    mis_bias[i] = shuffledResult["mutual.info.bias"];
  }
  
  result["shuffle.mi"] = mis;
  result["shuffle.mi.bias"] = mis_bias;
  
  return result;
}



/*** R
# Perfect mutual info
response = c(c(0, 10, 0, 10, 10), rep(0, 10))
stimulus = c(c(0, 1, 0, 1, 1), rep(0, 10))
res = mutual_info_with_shuffles(response, stimulus, c(1, 10), c(length(response)), 2, 1, 2)
res$mutual.info
res$mutual.info.bias
*/

