// Calculates place fields and two place cell related measures: Spatial information and Mutual information.
// The spatial information is tested for significance by shuffling the traces and comparing against the shuffled values
//
// The calculated mutual information is biased due to a finite and limited sampling of responses in the space.
// This bias is calculated using method from Panzeri and Treves, 1996.
//

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

class SpatialInfoData {
public:
  NumericMatrix totalActivityMap;
  NumericMatrix occupancyMap;
  NumericMatrix fr;
  double mfr = 0.0;
  double SI = 0.0;
  double MI = 0.0;
  double MI_bias = 0.0;
  double space_sampling_factor = 0.0;
  double sparsity = 0.0;

  SpatialInfoData() {
    totalActivityMap = NumericMatrix(1,1);
    occupancyMap = NumericMatrix(1,1);
    fr = NumericMatrix(1,1);
  }

  SpatialInfoData(NumericMatrix& totalActivityMap,
                  NumericMatrix& occupancyMap,
                  NumericMatrix& fr,
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

SpatialInfoData calculateSpatialInformation(NumericVector& x,
                                            NumericVector& y,
                                            NumericVector& trace,
                                            // Quantiles for binning trace values
                                            NumericVector& traceQuantiles,
                                            int timebin_size) {

  NumericMatrix totalActivityMap = NumericMatrix(N_BINS_X,N_BINS_Y);
  NumericMatrix occupancyMap = NumericMatrix(N_BINS_X,N_BINS_Y);
  NumericMatrix fr = NumericMatrix(N_BINS_X,N_BINS_Y);
  std::fill(totalActivityMap.begin(), totalActivityMap.end(), 0);
  std::fill(occupancyMap.begin(), occupancyMap.end(), 0);
  std::fill(fr.begin(), fr.end(), 0);

  int nresponse = traceQuantiles.size();
  std::vector<std::vector<std::vector<int> > > binnedResponse(nresponse, std::vector<std::vector<int> > (N_BINS_X, std::vector<int>(N_BINS_Y,0)));

  std::vector<int> responseCount(traceQuantiles.size(), 0);

  // Calculate occupancy and total activity maps
  int binned_trace_size = trace.size() / timebin_size;
  Debug("Trace Quantiles size=" << traceQuantiles.size() << std::endl);
  double mfr = 0;
  std::vector<int> bin_x(timebin_size);
  std::vector<int> bin_y(timebin_size);
  for (int i = 0; i < trace.size(); i += timebin_size) {
    std::fill(bin_x.begin(), bin_x.end(), 0);
    std::fill(bin_y.begin(), bin_y.end(), 0);
    double trace_binval = 0.0;
    for (int j = 0; j < timebin_size; ++j) {
      bin_x[j] = std::min(N_BINS_X - 1, std::max(0, (int) std::floor(x[i+j] / binSizeX)));
      bin_y[j] = std::min(N_BINS_Y - 1, std::max(0, (int) std::floor(y[i+j] / binSizeY)));
      trace_binval += trace[i+j];
    }
    int xx = std::round(((double) std::accumulate(bin_x.begin(), bin_x.end(), 0)) / bin_x.size());
    int yy = std::round(((double) std::accumulate(bin_y.begin(), bin_y.end(), 0)) / bin_y.size());

    occupancyMap(xx,yy) += 1;
    trace_binval = trace_binval / timebin_size;
    totalActivityMap(xx,yy) += trace_binval;

    auto responseBinIt = std::lower_bound(traceQuantiles.begin(), traceQuantiles.end(), trace[i]);
    int responseBin = traceQuantiles.size() - 1;
    if (responseBinIt != traceQuantiles.end()) {
      responseBin = responseBinIt - traceQuantiles.begin();
    }

    ++binnedResponse[responseBin][xx][yy];
    ++responseCount[responseBin];

    mfr += trace_binval / binned_trace_size;
  }
  Debug("Response 0 count: " << responseCount[0] << " Response 1 count:" << responseCount[1] << std::endl);

  // Calculate firing rates
  double fr_offset = 100.0;
  double sparsity = 0.0;
  for (int yy = 0; yy < N_BINS_Y; ++yy) {
    for (int xx = 0; xx < N_BINS_X; ++xx) {
      if (occupancyMap(xx,yy) > 0) {
        fr(xx,yy) = totalActivityMap(xx,yy) / occupancyMap(xx,yy);
        fr_offset = std::min(fr_offset, fr(xx,yy));
        double p_s = (double) occupancyMap(xx, yy) / binned_trace_size;
        sparsity += p_s * std::pow(fr(xx,yy), 2) / std::pow(mfr, 2);
      } else {
        fr(xx,yy) = NAN;
      }
    }
  }

  //TODO: decide if want to keep the offset, if yes, uncomment below
  // fr_offset -= 10e-10;
  // Avoid FR==0, so the log's can be calculated everywhere
  fr_offset = -10e-10;
  // Update firing rates by the offset
  for (int yy = 0; yy < N_BINS_Y; ++yy) {
    for (int xx = 0; xx < N_BINS_X; ++xx) {
      if (!std::isnan(fr(xx,yy))) {
        fr(xx,yy) -= fr_offset;
      }
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
  const int S = N_BINS_X * N_BINS_Y;
  const int N = binned_trace_size;
  int totalResponseBins = 0;
  int occupiedBins = 0;
  for (int yy = 0; yy < N_BINS_Y; ++yy) {
    for (int xx = 0; xx < N_BINS_X; ++xx) {

      if (occupancyMap(xx,yy) > 0) {
        Debug(" xx=" << xx);
        Debug(" yy=" << yy << std::endl);
        ++occupiedBins;
        double p_occupancy = (double) occupancyMap(xx,yy) / N;
        if (fr(xx,yy) > 0.0) {
          double r_SI = p_occupancy * fr(xx,yy) / mfr * log2(fr(xx,yy) / mfr);
          Debug("partial SI=" << r_SI << std::endl);
          SI += r_SI;
        }

        for (int r = 0; r < nresponse; ++r) {
          double p_response = ((double) responseCount[r]) / N;
          double p_r_given_s = ((double) binnedResponse[r][xx][yy]) / occupancyMap(xx,yy);
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
SEXP getCppPlaceField(NumericVector& x,
                      NumericVector& y,
                      NumericVector& trace,
                      NumericVector& traceQuantiles,
                      NumericVector& trialEnds,
                      int nshuffles,
                      int shuffleChunkLength,
                      int timebin_size) {


  NumericMatrix field = NumericMatrix(N_BINS_X,N_BINS_Y);
  std::fill(field.begin(), field.end(), NAN_FIELD);

  NumericVector fieldCentre = NumericVector(2);
  NumericVector fieldMaxXY = NumericVector(2);
  NumericVector shuffleSI = NumericVector(nshuffles);
  NumericVector shuffleMI = NumericVector(nshuffles);
  std::fill(shuffleSI.begin(), shuffleSI.end(), 0);
  std::fill(shuffleMI.begin(), shuffleMI.end(), 0);

  List result;
  result["field"] = NumericMatrix(N_BINS_X,N_BINS_Y);
  result["occupancy"] = NumericMatrix(N_BINS_X,N_BINS_Y);
  result["activity"] = NumericMatrix(N_BINS_X,N_BINS_Y);
  result["spatial.information"] = 0.0;
  result["spatial.information.perspike"] = 0.0;
  result["field.centre"] = fieldCentre;
  result["field.size.50"] = 0;
  result["field.size.25"] = 0;
  result["field.max"] = 0;
  result["field.max.xy"] = fieldMaxXY;
  result["shuffle.si"]= shuffleSI;
  result["shuffle.mi"]= shuffleMI;
  result["space.sampling.factor"] = 0.0;
  result["sparsity"] = 0.0;

  if (trace.size() == 0) {
    return(result);
  }

  SpatialInfoData spatialInfoData = calculateSpatialInformation(x, y, trace, traceQuantiles,
                                                                timebin_size);
  NumericMatrix occupancyMap = spatialInfoData.occupancyMap;
  NumericMatrix totalActivityMap = spatialInfoData.totalActivityMap;
  NumericMatrix fr = spatialInfoData.fr;

  double maxField = 0.0;
  for (int yy = 0; yy < N_BINS_Y; ++yy) {
    for (int xx = 0; xx < N_BINS_X; ++xx) {
      if (occupancyMap(xx,yy) > 0) {
        field(xx,yy) = totalActivityMap(xx,yy) / occupancyMap(xx,yy);
        if (field(xx,yy) >= maxField) {
          maxField = field(xx,yy);
          fieldMaxXY[0] = xx * binSizeX;
          fieldMaxXY[1] = yy * binSizeY;
        }
      } else {
        field(xx,yy) = NAN_FIELD;
      }
    }
  }
  // Populate NAs in the field with mean value and find field centre
  int nfield = 1; // avoid division by zero
  double fieldTotal = 0.0;
  const double FIELD_BINARY_THRESH = 0.5;
  double weightedFieldX = 0.0;
  double weigthedFieldY = 0.0;
  double totalWeights = 0.0;
  int nFieldBins50 = 0;
  int nFieldBins25 = 0;
  for (int yy = 0; yy < N_BINS_Y; ++yy) {
    for (int xx = 0; xx < N_BINS_X; ++xx) {
      if (field(xx,yy) != NAN_FIELD) {
        fieldTotal += field(xx, yy);
        ++nfield;

        if (field(xx,yy) >= FIELD_BINARY_THRESH * maxField) {
          double weight = field(xx,yy);
          weightedFieldX += xx * weight;
          weigthedFieldY += yy * weight;
          totalWeights += weight;
          ++nFieldBins50;
        }
        if (field(xx,yy) >= 0.25 * maxField) {
          ++nFieldBins25;
        }
      }
    }
  }
  fieldCentre[0] = weightedFieldX / std::max(0.01, totalWeights) * binSizeX;
  fieldCentre[1] = weigthedFieldY / std::max(0.01, totalWeights) * binSizeY;

  // TODO: should remove the code below?
  //double mean_field = fieldTotal / nfield;
  //for (int yy = 0; yy < N_BINS_Y; ++yy) {
  //  for (int xx = 0; xx < N_BINS_X; ++xx) {
  //    if (field(xx,yy) == NAN_FIELD) {
  //      field(xx,yy) = mean_field;
  //    }
  //  }
  //}

  for (int i = 0; i < nshuffles; ++i) {

    NumericVector shuffledTrace = chunkShuffle(trace, trialEnds, shuffleChunkLength);
    //NumericVector shuffledTrace = randomShift(trace, trialEnds, shuffleChunkLength * 2);
    SpatialInfoData shuffleData = calculateSpatialInformation(x, y, shuffledTrace,
                                                              traceQuantiles, timebin_size);
    shuffleSI[i] = shuffleData.SI;
    shuffleMI[i] = shuffleData.MI;

  }


  // Populate the result list
  result["field"] = field;
  result["occupancy"] = occupancyMap;
  result["activity"] = totalActivityMap;
  result["spatial.information"] = (float) spatialInfoData.SI;
  result["spatial.information.perspike"] = spatialInfoData.SI / spatialInfoData.mfr;
  result["mfr"] = spatialInfoData.mfr;
  result["field.centre"] = fieldCentre;
  result["field.size.50"] = ((double) nFieldBins50) / (N_BINS_X * N_BINS_Y) * 100.0;
  result["field.size.25"] = ((double) nFieldBins25) / (N_BINS_X * N_BINS_Y) * 100.0;
  result["field.max"] = maxField;
  result["field.max.xy"] = fieldMaxXY;
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
x=0:9
y=0:9
trace=rep(0, 10)
trace[6:10] = 1

pf = getCppPlaceField(x,y,trace, c(0.5, 0.9), c(7, length(trace)), 0, 2, 1)
#pf=with(cell.df, getCppPlaceField(smooth_trans_x, smooth_trans_y, deconv_trace, traceQuantiles, 10,2), 4)
pf$spatial.information
pf$mutual.info
pf$mfr

getNBinsXY()
*/

