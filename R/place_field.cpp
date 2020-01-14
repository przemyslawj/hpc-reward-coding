// Calculates place fields and two place cell related measures: Spatial information and Mutual information.
// The spatial information is tested for significance by shuffling the traces and comparing against the shuffled values
//
// The calculated mutual information is biased due to a finite and limited sampling of responses in the space.
// This bias is calculated using method from Panzeri and Treves, 1996.
//
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <Rcpp.h>
#include <cmath>
#include <float.h>
#include <vector>
#include <iostream>
#include <iterator>
#include <algorithm>
//#include "trace_algorithm.cpp" - can not be included using Rcpp


using namespace Rcpp;

//#define USEDEBUG

#ifdef USEDEBUG
#define Debug(x) std::cout << x
#else
#define Debug(x)
#endif

static const int NAN_FIELD = -1;
static const int binSizeX = 5;
static const int binSizeY = 5;

static const double MAX_X = 100;
static const double MAX_Y = 100;

static const int N_BINS_X = (int) std::floor(MAX_X / binSizeX);
static const int N_BINS_Y = (int) std::floor(MAX_Y / binSizeY);
static const int N_BINS = N_BINS_X * N_BINS_Y;

class SpatialInfoData {
public:
  NumericVector totalActivityMap;
  NumericVector occupancyMap;
  NumericVector fr;
  double mfr = 0.0;
  double SI = 0.0;
  double MI = 0.0;
  double MI_bias = 0.0;
  double space_sampling_factor = 0.0;
  double sparsity = 0.0;

  SpatialInfoData() {
    totalActivityMap = NumericVector(1, 0.0);
    occupancyMap = NumericVector(1, 0.0);
    fr = NumericVector(1, 0.0);
  }

  SpatialInfoData(NumericVector& totalActivityMap,
                  NumericVector& occupancyMap,
                  NumericVector& fr,
                  double mfr,
                  double SI,
                  double MI,
                  double MI_bias,
                  double space_sampling_factor,
                  double sparsity) {
    this->totalActivityMap = totalActivityMap;
    this->occupancyMap = occupancyMap;
    this->fr = fr;
    this->mfr = mfr;
    this->SI = SI;
    this->MI = MI;
    this->MI_bias = MI_bias;
    this->space_sampling_factor = space_sampling_factor;
    this->sparsity = sparsity;
  }
};

SpatialInfoData calculateSpatialInformation(NumericVector& bin_xy,
                                            NumericVector& trace,
                                            // Quantiles for binning trace values
                                            NumericVector& traceQuantiles,
                                            double minOccupancy) {

  NumericVector totalActivityMap = NumericVector(N_BINS, 0.0);
  NumericVector occupancyMap = NumericVector(N_BINS, 0.0);
  NumericVector fr = NumericVector(N_BINS, 0.0);

  int nresponse = traceQuantiles.size();
  arma::mat r_given_s(N_BINS, nresponse, arma::fill::zeros);

  std::vector<int> responseCount(traceQuantiles.size(), 0);

  // Calculate occupancy and total activity maps
  Debug("Trace Quantiles size=" << traceQuantiles.size() << std::endl);
  double mfr = 0;
  for (int i = 0; i < trace.size(); ++i) {
    
    int xy = bin_xy[i];

    occupancyMap[xy] += 1;
    totalActivityMap[xy] += trace[i];

    auto responseBinIt = std::lower_bound(traceQuantiles.begin(), traceQuantiles.end(), trace[i]);
    int responseBin = traceQuantiles.size() - 1;
    if (responseBinIt != traceQuantiles.end()) {
      responseBin = responseBinIt - traceQuantiles.begin();
    }

    ++r_given_s(xy,responseBin);
    ++responseCount[responseBin];

    mfr += trace[i] / trace.size();
  }
  Debug("Response 0 count: " << responseCount[0] << " Response 1 count:" << responseCount[1] << std::endl);

  // Calculate firing rates
  double fr_offset = 100.0;
  double sparsity = 0.0;
  for (int xy = 0; xy < N_BINS; ++xy) {
    if (occupancyMap[xy] >= minOccupancy) {
      fr[xy] = totalActivityMap[xy] / occupancyMap[xy];
      fr_offset = std::min(fr_offset, fr[xy]);
      double p_s = ((double) occupancyMap[xy]) / trace.size();
      sparsity += p_s * std::pow(fr[xy], 2) / std::pow(mfr, 2);
    } else {
      fr[xy] = NAN;
    }
  }

  //TODO: decide if want to keep the offset, if yes, uncomment below
  // fr_offset -= 10e-10;
  // Avoid FR==0, so the log's can be calculated everywhere
  fr_offset = -10e-10;
  // Update firing rates by the offset
  for (int xy = 0; xy < N_BINS; ++xy) {
    if (!std::isnan(fr[xy])) {
      fr[xy] -= fr_offset;
    }
  }

  // Calculate spatial information and the place field.
  mfr -= fr_offset;
  if (mfr == 0.0) {
    printf("MFR=0.0, calculation of spatial information aborted");
    return SpatialInfoData();
  }
  double MI = 0.0;
  double SI = 0.0;
  const int S = N_BINS;
  const int N = trace.size();
  int totalResponseBins = 0;
  int occupiedBins = 0;
  for (int xy = 0; xy < N_BINS; ++xy) {

    if (occupancyMap[xy] >= minOccupancy) {
      Debug(" xy=" << xy << std::endl);
      ++occupiedBins;
      double p_occupancy = (double) occupancyMap[xy] / N;
      if (fr[xy] > 0.0) {
        double r_SI = p_occupancy * fr[xy] / mfr * log2(fr[xy] / mfr);
        Debug("partial SI=" << r_SI << std::endl);
        SI += r_SI;
      }

      for (int r = 0; r < nresponse; ++r) {
        double p_response = ((double) responseCount[r]) / N;
        double p_r_given_s = ((double) r_given_s(xy, r)) / occupancyMap[xy];
        if (p_r_given_s > 0) {
          ++totalResponseBins;
          double r_MI = p_occupancy * p_r_given_s * std::log2(p_r_given_s / p_response);
          Debug("response=" << r);
          Debug(", p_r_given_s=" << p_r_given_s);
          Debug(", p_response=" << p_response);
          Debug(", p_occupancy=" << p_occupancy);
          Debug(", MI=" << r_MI <<std::endl);
          MI += r_MI;
          // This is bias just for xx,yy
          //C += (1 - p_r_given_s) / p_occupancy +
          //     (-p_r_given_s - 2 * p_r_given_s * p_r_given_s) / p_response - p_r_given_s;
        }
      }
    }
  }

  // Mutual information bias estimate from: Analytical estimates of limited sampling biases in different
  // information measures, Panzeri & Treves, 1996.
  double MI_bias = (totalResponseBins - nresponse - occupiedBins + 1 )  / (2 * N * std::log(2));
  Debug("totalResponseBins=" << totalResponseBins);
  Debug(",N=" << N);
  Debug(", MI_bias=" << MI_bias <<std::endl);
  double space_sampling_factor = ((double) occupiedBins) / S;
  return SpatialInfoData(totalActivityMap, occupancyMap, fr, mfr, SI, MI, MI_bias, space_sampling_factor, sparsity);
}


// [[Rcpp::export]]
int getNBinsXY() {
  return N_BINS_X;
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

NumericVector randomShift(NumericVector& trace,
                          NumericVector& trialEnds, 
                          int minShift) {
  NumericVector shiftedTrace(trace.size(), 0);
  
  int trialStart = 0;
  for (int trial_i = 0; trial_i < trialEnds.size(); ++trial_i) {
    int ntrial = trialEnds[trial_i] - trialStart;
    int shift = std::rand() % (ntrial - 2 * minShift) + minShift;
    for (int i = 0; i < ntrial; ++i) {
      int within_trial_i = (i + shift) % ntrial;
      shiftedTrace[trialStart + i] = trace[trialStart + within_trial_i];
    }
    
    trialStart = trialEnds[trial_i];
  }
  
  return shiftedTrace;
}

// [[Rcpp::export]]
SEXP getCppPlaceField(NumericVector& bin_xy,
                      NumericVector& trace,
                      NumericVector& traceQuantiles,
                      NumericVector& trialEnds,
                      int nshuffles,
                      int shuffleChunkLength,
                      double minOccupancy) {

  NumericVector shuffleSI = NumericVector(nshuffles, 0.0);
  NumericVector shuffleMI = NumericVector(nshuffles, 0.0);

  List result;
  result["field"] = NumericVector(N_BINS, 0.0);
  result["occupancy"] = NumericVector(N_BINS, 0.0);
  result["activity"] = NumericVector(N_BINS, 0.0);
  result["spatial.information"] = 0.0;
  result["spatial.information.perspike"] = 0.0;
  result["field.size.50"] = 0;
  result["field.size.25"] = 0;
  result["shuffle.si"]= shuffleSI;
  result["shuffle.mi"]= shuffleMI;
  result["space.sampling.factor"] = 0.0;
  result["sparsity"] = 0.0;

  if (trace.size() == 0) {
    return(result);
  }

  SpatialInfoData spatialInfoData = calculateSpatialInformation(bin_xy, trace, traceQuantiles, minOccupancy);
  NumericVector occupancyMap = spatialInfoData.occupancyMap;
  NumericVector totalActivityMap = spatialInfoData.totalActivityMap;
  NumericVector fr = spatialInfoData.fr;

  double maxField = 0.0;
  for (int xy = 0; xy < N_BINS; ++xy) {
    if (occupancyMap[xy] >= minOccupancy &&
        fr[xy] >= maxField) {
      maxField = fr[xy];
    }
  }
  // Populate NAs in the field with mean value and find field centre
  int nfield = 1; // avoid division by zero
  double fieldTotal = 0.0;
  const double FIELD_BINARY_THRESH = 0.5;
  double weightedFieldXy = 0.0;
  double totalWeights = 0.0;
  int nFieldBins50 = 0;
  int nFieldBins25 = 0;
  for (int xy = 0; xy < N_BINS; ++xy) {
    if (fr[xy] != NAN_FIELD) {
      fieldTotal += fr[xy];
      ++nfield;

      if (fr[xy] >= FIELD_BINARY_THRESH * maxField) {
        double weight = fr[xy];
        weightedFieldXy += xy * weight;
        totalWeights += weight;
        ++nFieldBins50;
      }
      if (fr[xy] >= 0.25 * maxField) {
        ++nFieldBins25;
      }
    }
  }

  for (int i = 0; i < nshuffles; ++i) {
    //NumericVector shuffledTrace = chunkShuffle(trace, trialEnds, shuffleChunkLength);
    NumericVector shuffledTrace = randomShift(trace, trialEnds, shuffleChunkLength * 2);
    SpatialInfoData shuffleData = calculateSpatialInformation(bin_xy, shuffledTrace, traceQuantiles, minOccupancy);
    shuffleSI[i] = shuffleData.SI;
    shuffleMI[i] = shuffleData.MI;
  }


  // Populate the result list
  result["field"] = fr;
  result["occupancy"] = occupancyMap;
  result["activity"] = totalActivityMap;
  result["spatial.information"] = (float) spatialInfoData.SI;
  result["spatial.information.perspike"] = spatialInfoData.SI / spatialInfoData.mfr;
  result["mfr"] = spatialInfoData.mfr;
  result["field.size.50"] = ((double) nFieldBins50) / N_BINS * 100.0;
  result["field.size.25"] = ((double) nFieldBins25) / N_BINS * 100.0;
  result["shuffle.si"] = shuffleSI;
  result["mutual.info"] = spatialInfoData.MI;
  result["mutual.info.bias"] = spatialInfoData.MI_bias;
  result["shuffle.mi"] = shuffleMI;
  result["space.sampling.factor"] = spatialInfoData.space_sampling_factor;
  result["sparsity"] = spatialInfoData.sparsity;
  return(result);
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//
/*** R
xy=0:9
trace=rep(0, 10)
trace[6:10] = 1

pf = getCppPlaceField(xy, trace, c(0.5, 0.9), c(7, length(trace)), 0, 2, 1)
#pf=with(cell.df, getCppPlaceField(smooth_trans_x, smooth_trans_y, deconv_trace, traceQuantiles, 10,2), 4)
pf$spatial.information
pf$mutual.info
pf$mfr

getNBinsXY()
*/

