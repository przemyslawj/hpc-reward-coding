#include <Rcpp.h>
#include <cmath>
#include <float.h>
#include <vector>
#include <random>
#include <iostream>

using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

static const int NAN_FIELD = -1;
static const int binSizeX = 3;
static const int binSizeY = 3;

class SpatialInfoData {
public:
  NumericMatrix totalActivityMap;
  NumericMatrix occupancyMap;
  NumericMatrix fr;
  double mfr;
  double SI;

  SpatialInfoData() {
    totalActivityMap = NumericMatrix(1,1);
    occupancyMap = NumericMatrix(1,1);
  }
  
  SpatialInfoData(NumericMatrix& totalActivityMap,
                  NumericMatrix& occupancyMap,
                  NumericMatrix& fr,
                  double mfr,
                  double SI) {
    this->totalActivityMap = totalActivityMap;
    this->occupancyMap = occupancyMap;
    this->fr = fr;
    this->mfr = mfr;
    this->SI = SI;
  }
};

SpatialInfoData calculateSpatialInformation(int nXbins, int nYbins, 
                                            NumericVector& x, 
                                            NumericVector& y, 
                                            NumericVector trace) {
  
  NumericMatrix totalActivityMap = NumericMatrix(nXbins,nYbins);
  NumericMatrix occupancyMap = NumericMatrix(nXbins,nYbins);
  NumericMatrix fr = NumericMatrix(nXbins,nYbins);
  std::fill(totalActivityMap.begin(), totalActivityMap.end(), 0);
  std::fill(occupancyMap.begin(), occupancyMap.end(), 0);
  std::fill(fr.begin(), fr.end(), 0);
  
  // Calculate occupancy and total activity maps
  double mfr = 0;
  for (int i = 0; i < trace.size(); ++i) {
    int xx = std::min(nXbins, std::max(0, (int) std::floor(x[i] / binSizeX)));
    int yy = std::min(nYbins, std::max(0, (int) std::floor(y[i] / binSizeY)));
    occupancyMap(xx,yy) += 1;
    totalActivityMap(xx,yy) += trace[i];
    
    mfr += trace[i] / trace.size();
  }
  
  // Calculate firing rates
  double fr_offset = 100.0;
  for (int yy = 0; yy < nYbins; ++yy) {
    for (int xx = 0; xx < nXbins; ++xx) {
      if (occupancyMap(xx,yy) > 0) {
        fr(xx,yy) = totalActivityMap(xx,yy) / occupancyMap(xx,yy);
        fr_offset = std::min(fr_offset, fr(xx,yy));
      } else {
        fr(xx,yy) = NAN;
      }
    }
  }
  
  // Avoid FR==0, so the log's can be calculated everywhere
  fr_offset -= 10e-10;
  
  // Update firing rates by the offset
  for (int yy = 0; yy < nYbins; ++yy) {
    for (int xx = 0; xx < nXbins; ++xx) {
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
  
  double SI = 0.0;
  for (int yy = 0; yy < nYbins; ++yy) {
    for (int xx = 0; xx < nXbins; ++xx) {
      if (occupancyMap(xx,yy) > 0 && !std::isnan(fr(xx,yy))) {
        double occupancyProb = (double) occupancyMap(xx,yy) / trace.size();
        SI += occupancyProb * fr(xx,yy) / mfr * log2(fr(xx,yy) / mfr);
      }
    }
  }
  
  return SpatialInfoData(totalActivityMap, occupancyMap, fr, mfr, SI);
}


// [[Rcpp::export]]
SEXP getCppPlaceField(NumericVector& x, 
                      NumericVector& y, 
                      NumericVector& trace, 
                      int nshuffles,
                      int shuffleChunkLength) {

  const double maxX = 100;
  const double maxY = 100;

  const int nXbins = (int) std::floor(maxX / binSizeX) + 1;
  const int nYbins = (int) std::floor(maxY / binSizeY) + 1;
  
  NumericMatrix field = NumericMatrix(nXbins,nYbins);
  std::fill(field.begin(), field.end(), NAN_FIELD);

  NumericVector fieldCentre = NumericVector(2);
  NumericVector fieldMaxXY = NumericVector(2);
  NumericVector shuffleSI = NumericVector(nshuffles);

  List result;
  result["field"] = NumericMatrix(nXbins,nYbins);
  result["occupancy"] = NumericMatrix(nXbins,nYbins);
  result["activity"] = NumericMatrix(nXbins,nYbins);
  result["spatial.information"] = 0.0;
  result["spatial.information.perspike"] = 0.0;
  result["field.centre"] = fieldCentre;
  result["field.size"] = 0;
  result["field.max"] = 0;
  result["field.max.xy"] = fieldMaxXY;
  result["shuffle.si"]= shuffleSI;
  
  if (trace.size() == 0) {
    return(result);
  }
  
  SpatialInfoData spatialInfoData = calculateSpatialInformation(nXbins, nYbins, x, y, trace);
  NumericMatrix occupancyMap = spatialInfoData.occupancyMap;
  NumericMatrix totalActivityMap = spatialInfoData.totalActivityMap;
  NumericMatrix fr = spatialInfoData.fr;
  
  double maxField = 0.0;
  for (int yy = 0; yy < nYbins; ++yy) {
    for (int xx = 0; xx < nXbins; ++xx) {
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
  for (int yy = 0; yy < nYbins; ++yy) {
    for (int xx = 0; xx < nXbins; ++xx) {
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
  for (int yy = 0; yy < nYbins; ++yy) {
    for (int xx = 0; xx < nXbins; ++xx) {
      if (field(xx,yy) == NAN_FIELD) {
        field(xx,yy) = mean_field;
      }
    }
  }
  
  // Shuffle the trace keeping the order within small chunks 
  auto rng = std::default_random_engine {};
  int nchunks = std::ceil(((double) trace.size() / shuffleChunkLength));
  std::vector<int> chunkShuffle(nchunks);
  for (int j = 0; j < nchunks; ++j) {
    chunkShuffle[j] = j;
  }
  for (int i = 0; i < nshuffles; ++i) {
    NumericVector shuffledTrace(trace.size());
    // TODO: shuffle only within the trial?
    std::shuffle(std::begin(chunkShuffle), std::end(chunkShuffle), rng);
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

    SpatialInfoData shuffleData = calculateSpatialInformation(nXbins, nYbins, x, y, shuffledTrace);
    //std::cout << shuffledTrace << std::endl;
    //std::cout << shuffleData.SI << std::endl;
    shuffleSI[i] = shuffleData.SI;
  }


  // Populate the result list
  result["field"] = field;
  result["occupancy"] = occupancyMap;
  result["activity"] = totalActivityMap;
  result["spatial.information"] = (float) spatialInfoData.SI;
  result["spatial.information.perspike"] = spatialInfoData.SI / spatialInfoData.mfr;
  result["mfr"] = spatialInfoData.mfr;
  result["field.centre"] = fieldCentre;
  result["field.size"] = ((double) nFieldBins) / (nXbins * nYbins) * 100.0;
  result["field.max"] = maxField;
  result["field.max.xy"] = fieldMaxXY;
  result["shuffle.si"] = shuffleSI;
  return(result);
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
x=0:10
y=0:10
trace=rep(0, 11)
trace[4:6] = 1
pf = getCppPlaceField(x,y,trace, 10, 2)
#pf=with(cell.df, getCppPlaceField(smooth_trans_x, smooth_trans_y, deconv_trace, 10))
pf$spatial.information
*/

