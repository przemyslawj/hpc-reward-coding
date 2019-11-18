#include <Rcpp.h>
#include <cmath>
#include <float.h>
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

// [[Rcpp::export]]
SEXP getCppPlaceField(NumericVector& x, NumericVector& y, NumericVector& trace) {
  const int NAN_FIELD = -1;
  const int binSizeX = 3;
  const int binSizeY = 3;

  double maxX = *std::max_element(x.begin(), x.end());
  double maxY = *std::max_element(y.begin(), y.end());

  const int nXbins = (int) std::floor(maxX / binSizeX) + 1;
  const int nYbins = (int) std::floor(maxY / binSizeY) + 1;

  NumericMatrix totalActivityMap = NumericMatrix(nXbins,nYbins);
  NumericMatrix occupancyMap = NumericMatrix(nXbins,nYbins);
  NumericMatrix field = NumericMatrix(nXbins,nYbins);
  std::fill(occupancyMap.begin(), occupancyMap.end(), 0);
  std::fill(totalActivityMap.begin(), totalActivityMap.end(), 0);
  std::fill(field.begin(), field.end(), NAN_FIELD);
  
  NumericVector fieldCentre = NumericVector(2);

  List result;
  result["field"] = field;
  result["occupancy"] = occupancyMap;
  result["activity"] = totalActivityMap;
  result["spatial.information"] = 0.0;
  result["spatial.information.perspike"] = 0.0;
  result["field.centre"] = fieldCentre;
  result["field.size"] = 0;
  if (trace.size() == 0) {
    return(result);
  }


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
  double fr[nXbins][nYbins];
  for (int yy = 0; yy < nYbins; ++yy) {
    for (int xx = 0; xx < nXbins; ++xx) {
      if (occupancyMap(xx,yy) > 0) {
        fr[xx][yy] = totalActivityMap(xx,yy) / occupancyMap(xx,yy);
        fr_offset = std::min(fr_offset, fr[xx][yy]);
      } else {
        fr[xx][yy] = NAN;
      }
    }
  }

  // Avoid FR==0, so the log's can be calculated everywhere
  fr_offset -= 10e-10;

  // Update firing rates by the offset
  for (int yy = 0; yy < nYbins; ++yy) {
    for (int xx = 0; xx < nXbins; ++xx) {
      if (!std::isnan(fr[xx][yy])) {
        fr[xx][yy] -= fr_offset;
      }
    }
  }

  // Calculate spatial information and the place field.
  mfr -= fr_offset;
  if (mfr == 0.0) {
    printf("MFR=0.0, calculation of spatial information aborted");
    return(result);
  }

  double SI = 0.0;
  double maxField = 0.0;
  for (int yy = 0; yy < nYbins; ++yy) {
    for (int xx = 0; xx < nXbins; ++xx) {
      if (occupancyMap(xx,yy) > 0 && !std::isnan(fr[xx][yy])) {
        double occupancyProb = (double) occupancyMap(xx,yy) / trace.size();
        SI += occupancyProb * fr[xx][yy] / mfr * log2(fr[xx][yy] / mfr);
        field(xx,yy) = totalActivityMap(xx,yy) / occupancyMap(xx,yy);
        maxField = std::max(field(xx,yy), maxField);
      } else {
        field(xx,yy) = NAN_FIELD;
      }
    }
  }

  // Populate NAs in the field with mean value and find field centre
  int nfield = 1; // avoid division by zero
  double fieldTotal = 0.0;
  const double FIELD_BINARY_THRESH = 0.8;
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


  // Populate the result list
  result["field"] = field;
  result["occupancy"] = occupancyMap;
  result["activity"] = totalActivityMap;
  result["spatial.information"] = (float) SI;
  result["spatial.information.perspike"] = SI/mfr;
  result["field.centre"] = fieldCentre;
  result["field.size"] = (int) std::round(nFieldBins / (nXbins * nYbins) * 100);
  return(result);
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
x=1:10
y=1:10
trace=rep(0, 10)
pf = getCppPlaceField(x,y,trace)
#pf=with(cell.df, getCppPlaceField(smooth_trans_x, smooth_trans_y, deconv_trace))
#pf$spatial.information
*/

