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
#include <random>
#include <iostream>
#include <iterator>
#include <algorithm>

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

static const int N_BINS_X = (int) std::floor(MAX_X / binSizeX) + 1;
static const int N_BINS_Y = (int) std::floor(MAX_Y / binSizeY) + 1;

class SpatialInfoData {
public:
  NumericMatrix totalActivityMap;
  NumericMatrix occupancyMap;
  NumericMatrix fr;
  double mfr;
  double SI;
  double MI;
  double MI_bias;

  SpatialInfoData() {
    totalActivityMap = NumericMatrix(1,1);
    occupancyMap = NumericMatrix(1,1);
  }
  
  SpatialInfoData(NumericMatrix& totalActivityMap,
                  NumericMatrix& occupancyMap,
                  NumericMatrix& fr,
                  double mfr,
                  double SI,
                  double MI,
                  double MI_bias) {
    this->totalActivityMap = totalActivityMap;
    this->occupancyMap = occupancyMap;
    this->fr = fr;
    this->mfr = mfr;
    this->SI = SI;
    this->MI = MI;
    this->MI_bias = MI_bias;
  }
};

SpatialInfoData calculateSpatialInformation(NumericVector& x, 
                                            NumericVector& y, 
                                            NumericVector trace,
                                            // Quantiles for binning trace values
                                            NumericVector& traceQuantiles,
                                            int timebin_size) {
  
  NumericMatrix totalActivityMap = NumericMatrix(N_BINS_X,N_BINS_Y);
  NumericMatrix occupancyMap = NumericMatrix(N_BINS_X,N_BINS_Y);
  NumericMatrix fr = NumericMatrix(N_BINS_X,N_BINS_Y);
  std::fill(totalActivityMap.begin(), totalActivityMap.end(), 0);
  std::fill(occupancyMap.begin(), occupancyMap.end(), 0);
  std::fill(fr.begin(), fr.end(), 0);
  
  std::vector<NumericMatrix> binnedResponse(traceQuantiles.size());
  for (int i = 0; i < traceQuantiles.size(); ++i) {
    binnedResponse[i] = NumericMatrix(N_BINS_X, N_BINS_Y);
    std::fill(binnedResponse[i].begin(), binnedResponse[i].end(), 0);
  }
  
  std::vector<int> responseCount(traceQuantiles.size());
  std::fill(responseCount.begin(), responseCount.end(), 0);
  
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
      bin_x[j] = std::min(N_BINS_X, std::max(0, (int) std::floor(x[i+j] / binSizeX)));
      bin_y[j] = std::min(N_BINS_Y, std::max(0, (int) std::floor(y[i+j] / binSizeY)));
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
    ++binnedResponse[responseBin](xx,yy);
    ++responseCount[responseBin];
    
    mfr += trace_binval / binned_trace_size;
  }
  Debug("Response 0 count: " << responseCount[0] << " Response 1 count:" << responseCount[1] << std::endl);
  
  // Calculate firing rates
  double fr_offset = 100.0;
  for (int yy = 0; yy < N_BINS_Y; ++yy) {
    for (int xx = 0; xx < N_BINS_X; ++xx) {
      if (occupancyMap(xx,yy) > 0) {
        fr(xx,yy) = totalActivityMap(xx,yy) / occupancyMap(xx,yy);
        fr_offset = std::min(fr_offset, fr(xx,yy));
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
  for (int yy = 0; yy < N_BINS_Y; ++yy) {
    for (int xx = 0; xx < N_BINS_X; ++xx) {
      
      if (occupancyMap(xx,yy) > 0) {
        Debug(" xx=" << xx);
        Debug(" yy=" << yy << std::endl);
        double p_occupancy = (double) occupancyMap(xx,yy) / N;
        if (fr(xx,yy) > 0.0) {
          double r_SI = p_occupancy * fr(xx,yy) / mfr * log2(fr(xx,yy) / mfr);
          Debug("partial SI=" << r_SI << std::endl);
          SI += r_SI;
        }
        
        for (int r = 0; r < traceQuantiles.size(); ++r) {
          double p_response = ((double) responseCount[r]) / N;
          double p_r_given_s = ((double) binnedResponse[r](xx,yy)) / occupancyMap(xx,yy);
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
  
  double MI_bias = (totalResponseBins - traceQuantiles.size() - S - 1 )  / (2 * N * std::log(2));
  Debug("totalResponseBins=" << totalResponseBins);
  Debug(",N=" << N);
  Debug(", MI_bias=" << MI_bias <<std::endl);
  return SpatialInfoData(totalActivityMap, occupancyMap, fr, mfr, SI, MI, MI_bias);
}


// [[Rcpp::export]]
int getNBinsXY() {
  return N_BINS_X;
}

// Exports the C++ function to R using the Rcpp::sourceCpp
// function.
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

  List result;
  result["field"] = NumericMatrix(N_BINS_X,N_BINS_Y);
  result["occupancy"] = NumericMatrix(N_BINS_X,N_BINS_Y);
  result["activity"] = NumericMatrix(N_BINS_X,N_BINS_Y);
  result["spatial.information"] = 0.0;
  result["spatial.information.perspike"] = 0.0;
  result["field.centre"] = fieldCentre;
  result["field.size"] = 0;
  result["field.max"] = 0;
  result["field.max.xy"] = fieldMaxXY;
  result["shuffle.si"]= shuffleSI;
  result["shuffle.mi"]= shuffleMI;
  
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
  int nFieldBins = 0;
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
          ++nFieldBins;
        }
      }
    }
  }
  fieldCentre[0] = weightedFieldX / std::max(0.01, totalWeights) * binSizeX;
  fieldCentre[1] = weigthedFieldY / std::max(0.01, totalWeights) * binSizeY;
  
  double mean_field = fieldTotal / nfield;
  for (int yy = 0; yy < N_BINS_Y; ++yy) {
    for (int xx = 0; xx < N_BINS_X; ++xx) {
      if (field(xx,yy) == NAN_FIELD) {
        field(xx,yy) = mean_field;
      }
    }
  }
  
  // Shuffle the trace keeping the order within small chunks 
  int nchunks = std::ceil(((double) trace.size() / shuffleChunkLength));
  std::vector<int> chunkShuffle(nchunks);
  for (int j = 0; j < nchunks; ++j) {
    chunkShuffle[j] = j;
  }
  for (int i = 0; i < nshuffles; ++i) {
    NumericVector shuffledTrace(trace.size());
    
    // Shuffle chunks only within the same trial
    int trialStart = 0;
    for (int trial_i = 0; trial_i < trialEnds.size(); ++trial_i) {
      int trialEnd = trialEnds[trial_i] / shuffleChunkLength + 1;
      std::random_shuffle(chunkShuffle.begin() + trialStart, chunkShuffle.begin() + trialEnd);  
      trialStart = trialEnd;
    }
    //for (int trial_i = 0; trial_i < trialEnds.size(); ++trial_i) {
    //  int trialLen = trialEnds[trial_i] / shuffleChunkLength - trialStart;
    //  for (int chunk_i = trialStart; chunk_i < trialEnds[trial_i] / shuffleChunkLength; ++i) {
    //    int j = chunk_i + rand() % (trialLen - chunk_i);
    //    std::swap(chunkShuffle[chunk_i], chunkShuffle[j]);
    //  }
    //  trialStart = (trialEnds[trial_i] - 1) / shuffleChunkLength;
    //}
    
    // Offset from the start to randomize chunk boundaries
    int offset = std::rand() % ((int) shuffleChunkLength/2);
    for (int chunk = 0; chunk < nchunks; ++chunk) {
      int trace_offset = chunk * shuffleChunkLength + offset;
      int start_i = chunkShuffle[chunk] * shuffleChunkLength + offset;
      for (int j = 0; j < shuffleChunkLength; ++j) {
        int trace_j = (trace_offset + j) % trace.size();
        int shuffled_j = (start_i + j) % trace.size();
        shuffledTrace[shuffled_j] = trace[trace_j];
      }
    }
    
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
  result["field.size"] = ((double) nFieldBins) / (N_BINS_X * N_BINS_Y) * 100.0;
  result["field.max"] = maxField;
  result["field.max.xy"] = fieldMaxXY;
  result["shuffle.si"] = shuffleSI;
  result["mutual.info"] = spatialInfoData.MI;
  result["mutual.info.bias"] = spatialInfoData.MI_bias;
  result["shuffle.mi"] = shuffleMI;
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

