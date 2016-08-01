/****************************************************************************
 *
 * This file is part of the ViSP software.
 * Copyright (C) 2005 - 2015 by Inria. All rights reserved.
 *
 * This software is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * ("GPL") version 2 as published by the Free Software Foundation.
 * See the file LICENSE.txt at the root directory of this source
 * distribution for additional information about the GNU GPL.
 *
 * For using ViSP with software that can not be combined with the GNU
 * GPL, please contact Inria about acquiring a ViSP Professional
 * Edition License.
 *
 * See http://visp.inria.fr for more information.
 *
 * This software was developed at:
 * Inria Rennes - Bretagne Atlantique
 * Campus Universitaire de Beaulieu
 * 35042 Rennes Cedex
 * France
 *
 * If you have questions regarding the use of this file, please contact
 * Inria at visp@inria.fr
 *
 * This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
 * WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
 *
 * Description:
 * Convert image types.
 *
 * Authors:
 * Souriya Trinh
 *
 *****************************************************************************/

/*!
  \file vpCLAHE.cpp
  \brief Contrast Limited Adaptive Histogram Equalization (CLAHE).
*/

#include <visp3/imgproc/vpImgproc.h>
#include <visp3/core/vpImageConvert.h>
#include <visp3/core/vpMath.h>


static void clipHistogram(const std::vector<int> &hist, std::vector<int> &clippedHist, const unsigned int limit);

static std::vector<int> createHistogram(const int blockRadius, const int bins, const int blockXCenter, const int blockYCenter,
                                        const vpImage<unsigned char> &I);

static std::vector<double> createTransfer(const std::vector<int> &hist, const int limit);

static double transferValue(const int v, std::vector<int> &clippedHist);
static double transferValue(const int v, const std::vector<int> &hist, std::vector<int> &clippedHist, const int limit);


/*!
  \ingroup group_imgproc_brightness
  
   Adjust the contrast of a grayscale image locally using the Contrast Limited Adaptative Histogram Equalization method.
   The limit parameter allows to limit the slope of the transformation function to prevent the overamplification of noise.
   This method is a transcription of the CLAHE ImageJ plugin code by Stephan Saalfeld.

   \param I1 : The first grayscale image.
   \param I2 : The second grayscale image after application of the CLAHE method.
   \param blockRadius : The size of the local region around a pixel for which the histogram is equalized.
   This size should be larger than the size of features to be preserved.
   \param  bins : The number of histogram bins used for histogram equalization.
   The number of histogram bins should be smaller than the number of pixels in a block.
   \param slope : Limits the contrast stretch in the intensity transfer function.
   Very large values will let the histogram equalization do whatever it wants to do, that is result in maximal local contrast.
   The value 1 will result in the original image.
   \param fast : Use the fast but less accurate version of the filter. The fast version does not evaluate the intensity
   transfer function for each pixel independently but for a grid of adjacent boxes of the given block size only
   and interpolates for locations in between.
*/
void vp::clahe(const vpImage<unsigned char> &I1, vpImage<unsigned char> &I2, const unsigned int blockRadius,
    const unsigned int bins, const double slope, const bool fast) {
  I2 = vpImage<unsigned char>(I1.getHeight(), I1.getWidth());

  if(blockRadius > I1.getWidth() || blockRadius > I1.getHeight()) {
    return;
  }

  if(fast) {
    int blockSize = 2 * blockRadius + 1;
    int limit = ( int )( slope * blockSize * blockSize / bins + 0.5 );

    /* div */
    int nc = I1.getWidth() / blockSize;
    int nr = I1.getHeight() / blockSize;

    /* % */
    int cm = I1.getWidth() - nc * blockSize;
    std::vector<int> cs;

    switch (cm) {
    case 0:
      cs.resize(nc);
      for (int i = 0; i < nc; ++i) {
        cs[i] = i * blockSize + blockRadius + 1;
      }
      break;

    case 1:
      cs.resize(nc + 1);
      for (int i = 0; i < nc; ++i) {
        cs[i] = i * blockSize + blockRadius + 1;
      }
      cs[nc] = I1.getWidth() - blockRadius - 1;
      break;

    default:
      cs.resize(nc + 2);
      cs[0] = blockRadius + 1;
      for (int i = 0; i < nc; ++i) {
        cs[i + 1] = i * blockSize + blockRadius + 1 + cm / 2;
      }
      cs[nc + 1] = I1.getWidth() - blockRadius - 1;
    }

    int rm = I1.getHeight() - nr * blockSize;
    std::vector<int> rs;

    switch (rm) {
    case 0:
      rs.resize(nr);
      for (int i = 0; i < nr; ++i) {
        rs[i] = i * blockSize + blockRadius + 1;
      }
      break;

    case 1:
      rs.resize(nr + 1);
      for (int i = 0; i < nr; ++i) {
        rs[i] = i * blockSize + blockRadius + 1;
      }
      rs[nr] = I1.getHeight() - blockRadius - 1;
      break;

    default:
      rs.resize(nr + 2);
      rs[0] = blockRadius + 1;
      for (int i = 0; i < nr; ++i) {
        rs[i + 1] = i * blockSize + blockRadius + 1 + rm / 2;
      }
      rs[nr + 1] = I1.getHeight() - blockRadius - 1;
    }

    std::vector<int> hist;
    std::vector<double> tl;
    std::vector<double> tr;
    std::vector<double> br;
    std::vector<double> bl;
    for (int r = 0; r <= (int) rs.size(); ++r) {
      int r0 = std::max(0, r - 1);
      int r1 = std::min((int) rs.size() - 1, r);
      int dr = rs[r1] - rs[r0];


      hist = createHistogram(blockRadius, bins, cs[0], rs[r0], I1);
      tr = createTransfer(hist, limit);
      if (r0 == r1) {
        br = tr;
      }
      else {
        hist = createHistogram(blockRadius, bins, cs[0], rs[r1], I1);
        br = createTransfer(hist, limit);
      }

      int yMin = (r == 0 ? 0 : rs[r0]);
      int yMax = (r < (int) rs.size() ? rs[r1] : I1.getHeight() - 1);
      for (int c = 0; c <= (int) cs.size(); ++c) {
        int c0 = std::max(0, c - 1);
        int c1 = std::min((int) cs.size() - 1, c);
        int dc = cs[c1] - cs[c0];
        tl = tr;
        bl = br;
        if (c0 != c1) {
          hist = createHistogram(blockRadius, bins, cs[c1], rs[r0], I1);
          tr = createTransfer(hist, limit);
          if (r0 == r1) {
            br = tr;
          }
          else {
            hist = createHistogram(blockRadius, bins, cs[c1], rs[r1], I1);
            br = createTransfer(hist, limit);
          }
        }

        int xMin = (c == 0 ? 0 : cs[c0]);
        int xMax = (c < (int) cs.size() ? cs[c1] : I1.getWidth() - 1);
        for (int y = yMin; y < yMax; ++y) {
          float wy = (float) (rs[r1] - y) / dr;

          for (int x = xMin; x < xMax; ++x) {
            float wx = (float) (cs[c1] - x) / dc;
            int v = vpMath::round(I1[y][x] / 255.0 * bins);
            float t00 = tl[v];
            float t01 = tr[v];
            float t10 = bl[v];
            float t11 = br[v];
            float t0, t1;

            if (c0 == c1) {
              t0 = t00;
              t1 = t10;
            } else {
              t0 = wx * t00 + (1.0f - wx) * t01;
              t1 = wx * t10 + (1.0f - wx) * t11;
            }

            float t;
            if (r0 == r1) {
              t = t0;
            }
            else {
              t = wy * t0 + (1.0f - wy) * t1;
            }

            I2[y][x] = std::max( 0, std::min(255, vpMath::round(t * 255.0)) );
          }
        }
      }
    }

  } else {
    for(int y = 0; y < (int) I1.getHeight(); y++) {
      int yMin = std::max(0, y - (int) blockRadius);
      int yMax = std::min(I1.getHeight(), y + blockRadius + 1);
      int h = yMax-yMin;

      int xMin0 = 0;
      int xMax0 = blockRadius;

      std::vector<int> hist(bins+1);
      std::vector<int> clippedHist(bins+1);

      for(int yi = yMin; yi < yMax; yi++) {
        for(int xi = xMin0; xi < xMax0; xi++) {
          ++hist[vpMath::round((double) (I1[yi][xi] / 255.0 * bins))];
        }
      }

      for(int x = 0; x < (int) I1.getWidth(); x++) {
        int v = vpMath::round((double) (I1[y][x] / 255.0 * bins));
        int xMin = std::max(0, x - (int) blockRadius);
        int xMax = x+blockRadius+1;
        int w = std::min((int) I1.getWidth(), xMax) - xMin;
        int n = h*w;

        int limit = slope * n / bins + 0.5;

        if(xMin > 0) {
          int xMin1 = xMin-1;
          for(int yi = yMin; yi < yMax; yi++) {
            --hist[ vpMath::round((double) (I1[yi][xMin1] / 255.0 * bins ))];
          }
        }

        if(xMax <= (int) I1.getWidth()) {
          int xMax1 = xMax-1;
          for(int yi = yMin; yi < yMax; yi++) {
            ++hist[ vpMath::round((double) (I1[yi][xMax1] / 255.0 * bins ))];
          }
        }

        I2[y][x] = vpMath::round(transferValue(v, hist, clippedHist, limit) * 255.0);
      }
    }
  }
}

/*!
  \ingroup group_imgproc_brightness
  
  Adjust the contrast of a color image locally using the Contrast Limited Adaptative Histogram Equalization method.
  The limit parameter allows to limit the slope of the transformation function to prevent the overamplification of noise.
  This method is a transcription of the CLAHE ImageJ plugin code by Stephan Saalfeld.

  \param I1 : The first color image.
  \param I2 : The second color image after application of the CLAHE method.
  \param blockRadius : The size of the local region around a pixel for which the histogram is equalized.
  This size should be larger than the size of features to be preserved.
  \param  bins : The number of histogram bins used for histogram equalization.
  The number of histogram bins should be smaller than the number of pixels in a block.
  \param slope : Limits the contrast stretch in the intensity transfer function.
  Very large values will let the histogram equalization do whatever it wants to do, that is result in maximal local contrast.
  The value 1 will result in the original image.
  \param fast : Use the fast but less accurate version of the filter. The fast version does not evaluate the intensity
  transfer function for each pixel independently but for a grid of adjacent boxes of the given block size only
  and interpolates for locations in between.
*/
void vp::clahe(const vpImage<vpRGBa> &I1, vpImage<vpRGBa> &I2, const unsigned int blockRadius,
    const unsigned int bins, const double slope, const bool fast) {
  //Split
  vpImage<unsigned char> pR(I1.getHeight(), I1.getWidth());
  vpImage<unsigned char> pG(I1.getHeight(), I1.getWidth());
  vpImage<unsigned char> pB(I1.getHeight(), I1.getWidth());
  vpImage<unsigned char> pa(I1.getHeight(), I1.getWidth());

  vpImageConvert::split(I1, &pR, &pG, &pB, &pa);

  //Apply CLAHE independently on RGB channels
  vpImage<unsigned char> resR, resG, resB;
  clahe(pR, resR, blockRadius, bins, slope, fast);
  clahe(pG, resG, blockRadius, bins, slope, fast);
  clahe(pB, resB, blockRadius, bins, slope, fast);

  I2 = vpImage<vpRGBa>(I1.getHeight(), I1.getWidth());
  unsigned int size = I2.getWidth()*I2.getHeight();
  unsigned char *ptrStart = (unsigned char*) I2.bitmap;
  unsigned char *ptrEnd = ptrStart + size*4;
  unsigned char *ptrCurrent = ptrStart;

  unsigned int cpt = 0;
  while(ptrCurrent != ptrEnd) {
    *ptrCurrent = resR.bitmap[cpt];
    ++ptrCurrent;

    *ptrCurrent = resG.bitmap[cpt];
    ++ptrCurrent;

    *ptrCurrent = resB.bitmap[cpt];
    ++ptrCurrent;

    *ptrCurrent = pa.bitmap[cpt];
    ++ptrCurrent;

    cpt++;
  }
}

static void clipHistogram(const std::vector<int> &hist, std::vector<int> &clippedHist, const unsigned int limit) {
  int histlength = hist.size();

  for(int cpt = 0; cpt < histlength; cpt++) {
    clippedHist[cpt] = hist[cpt];
  }

  int clippedEntries = 0, clippedEntriesBefore;

  do {
    clippedEntriesBefore = clippedEntries;
    clippedEntries = 0;
    for (int i = 0; i < histlength; i++) {
      int d = (clippedHist)[i] - limit;
      if (d > 0) {
        clippedEntries += d;
        (clippedHist)[i] = limit;
      }
    }

    int d = clippedEntries / (histlength);
    int m = clippedEntries % (histlength);
    for (int i = 0; i < histlength; i++) {
      (clippedHist)[i] += d;
    }

    if (m != 0) {
      int s = (histlength - 1) / m;
      for (int i = s / 2; i < histlength; i += s) {
        ++((clippedHist)[i]);
      }
    }
  } while (clippedEntries != clippedEntriesBefore);
}

static std::vector<int> createHistogram(const int blockRadius, const int bins, const int blockXCenter,
    const int blockYCenter, const vpImage<unsigned char> &I) {
  std::vector<int> hist(bins + 1);

  int xMin = std::max( 0, blockXCenter - blockRadius );
  int yMin = std::max( 0, blockYCenter - blockRadius );
  int xMax = std::min( (int) I.getWidth(), blockXCenter + blockRadius + 1 );
  int yMax = std::min( (int) I.getHeight(), blockYCenter + blockRadius + 1 );

  for (int y = yMin; y < yMax; ++y) {
    for (int x = xMin; x < xMax; ++x) {
      ++hist[vpMath::round(I[y][x] / 255.0f * bins)];
    }
  }

  return hist;
}

static std::vector<double> createTransfer(const std::vector<int> &hist, const int limit) {
  std::vector<int> cdfs(hist.size());
  clipHistogram(hist, cdfs, limit);
  int hMin = hist.size() - 1;

  for (int i = 0; i < hMin; ++i) {
    if (cdfs[i] != 0) {
      hMin = i;
    }
  }
  int cdf = 0;
  for (int i = hMin; i < (int) hist.size(); ++i) {
    cdf += cdfs[i];
    cdfs[i] = cdf;
  }

  int cdfMin = cdfs[hMin];
  int cdfMax = cdfs[hist.size() - 1];

  std::vector<double> transfer(hist.size());
  for (int i = 0; i < (int) transfer.size(); ++i) {
    transfer[i] = (cdfs[i] - cdfMin) / (float) (cdfMax - cdfMin);
  }

  return transfer;
}

static double transferValue(const int v, std::vector<int> &clippedHist) {
  int clippedHistLength = clippedHist.size();
  int hMin = clippedHistLength - 1;
  for (int i = 0; i < hMin; i++) {
    if (clippedHist[i] != 0) {
      hMin = i;
    }
  }

  int cdf = 0;
  for (int i = hMin; i <= v; i++) {
    cdf += clippedHist[i];
  }

  int cdfMax = cdf;
  for (int i = v + 1; i < clippedHistLength; ++i) {
    cdfMax += clippedHist[i];
  }

  int cdfMin = clippedHist[ hMin ];
  return ( cdf - cdfMin ) / ( double )( cdfMax - cdfMin );
}

static double transferValue(const int v, const std::vector<int> &hist, std::vector<int> &clippedHist, const int limit) {
  clipHistogram(hist, clippedHist, limit);

  return transferValue(v, clippedHist);
}
