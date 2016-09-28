/****************************************************************************
 * Copyright (c) 2011, The University of Southampton and the individual contributors.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 *
 *   * Redistributions of source code must retain the above copyright notice,
 *  this list of conditions and the following disclaimer.
 *
 *   * Redistributions in binary form must reproduce the above copyright notice,
 *  this list of conditions and the following disclaimer in the documentation
 *  and/or other materials provided with the distribution.
 *
 *   * Neither the name of the University of Southampton nor the names of its
 *  contributors may be used to endorse or promote products derived from this
 *  software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
 * ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
 * ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
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
 * Contours extraction.
 *
 * Authors:
 * Souriya Trinh
 *
 *****************************************************************************/

/*!
  \file vpContour.cpp
  \brief Basic contours extraction.
*/

#include <map>
#include <visp3/imgproc/vpImgproc.h>


bool fromTo(const vpImagePoint &from, const vpImagePoint &to, vp::vpDirection &direction) {
  if (from == to) {
    return false;
  }

  if ( std::fabs(from.get_i() - to.get_i()) < std::numeric_limits<double>::epsilon() ) {
    if (from.get_j() < to.get_j()) {
      direction.m_direction = vp::EAST;
    } else {
      direction.m_direction = vp::WEST;
    }
  } else if (from.get_i() < to.get_i()) {
    if ( std::fabs(from.get_j() - to.get_j()) < std::numeric_limits<double>::epsilon() ) {
      direction.m_direction = vp::SOUTH;
    } else if (from.get_j() < to.get_j()) {
      direction.m_direction = vp::SOUTH_EAST;
    } else {
      direction.m_direction = vp::SOUTH_WEST;
    }
  } else {
    if ( std::fabs(from.get_j() - to.get_j()) < std::numeric_limits<double>::epsilon() ) {
      direction.m_direction = vp::NORTH;
    } else if (from.get_j() < to.get_j()) {
      direction.m_direction = vp::NORTH_EAST;
    } else {
      direction.m_direction = vp::NORTH_WEST;
    }
  }

  return true;
}

bool crossesEastBorder(const vpImage<int> &I, bool checked[8], const vpImagePoint &point) {
  vp::vpDirection direction;
  if ( !fromTo(point, vpImagePoint(point.get_i(), point.get_j() + 1), direction) ) {
    return false;
  }

  bool b = checked[(int) direction.m_direction];

  unsigned int i = (unsigned int) point.get_i();
  unsigned int j = (unsigned int) point.get_j();
  return I[i][j] != 0 && ( (unsigned int) point.get_j() == I.getWidth()-1 || b);
}

void addContourPoint(vpImage<int> &I, vp::vpContour *border, const vpImagePoint &point, bool checked[8], const int nbd) {
  border->m_points.push_back(point);

  unsigned int i = (unsigned int) point.get_i();
  unsigned int j = (unsigned int) point.get_j();

  if (crossesEastBorder(I, checked, point)) {
    I[i][j] = -nbd;
  } else if (I[i][j] == 1) {
    I[i][j] = nbd;
  }
}

bool followBorder(vpImage<int> &I, const vpImagePoint &ij, vpImagePoint &i2j2, vp::vpContour *border, const int nbd) {
  vp::vpDirection dir;
  if (!fromTo(ij, i2j2, dir)) {
    throw vpException(vpException::fatalError, "ij == i2j2");
  }

  vp::vpDirection trace = dir.clockwise();
  vpImagePoint i1j1(-1, -1);

  while (trace.m_direction != dir.m_direction) {
    vpImagePoint activePixel = trace.active(I, ij);

    if (activePixel.get_i() >= 0 && activePixel.get_j() >= 0) {
      i1j1 = activePixel;
      break;
    }

    trace = trace.clockwise();
  }

  if (i1j1.get_i() < 0 || i1j1.get_j() < 0) {
    return false;
  }

  i2j2 = i1j1;
  vpImagePoint i3j3 = ij;

  bool checked[8] = {
    false, false, false, false, false, false, false, false
  };

  while (true) {
    if (!fromTo(i3j3, i2j2, dir)) {
      throw vpException(vpException::fatalError, "i3j3 == i2j2");
    }

    trace = dir.counterClockwise();
    vpImagePoint i4j4(-1, -1);

    //resetChecked
    for (int cpt = 0; cpt < 8; cpt++) {
      checked[cpt] = false;
    }

    bool problem = false;
    while (!problem) {
      i4j4 = trace.active(I, i3j3);
      if (i4j4.get_i() >= 0 && i4j4.get_j() >= 0) {
        break;
      }

      checked[(int) trace.m_direction] = true;
      trace = trace.counterClockwise();

      problem = true;
      for (int cpt = 0; cpt < 8; cpt++) {
        if (!checked[cpt]) {
          problem = false;
          break;
        }
      }
    }

    if (problem) {
      return false;
    }

    addContourPoint(I, border, i3j3, checked, nbd);

    if (i4j4 == ij && i3j3 == i1j1) {
      break;
    }

    i2j2 = i3j3;
    i3j3 = i4j4;
  }

  return true;
}

bool isOuterBorderStart(const vpImage<int> &I, unsigned int i, unsigned int j) {
  return (I[i][j] == 1 && (j == 0 || I[i][j - 1] == 0));
}

bool isHoleBorderStart(const vpImage<int> &I, unsigned int i, unsigned int j) {
  return (I[i][j] >= 1 && (j == I.getWidth()-1 || I[i][j + 1] == 0));
}

/*!
  \ingroup group_imgproc_contours

  Extract contours from a binary image.

  \param I : Input image (0 means background, 1 means foreground, other values are not allowed).
  \param contour : Detected contours.
*/
void vp::extractContours(const vpImage<unsigned char> &I_original, vpContour &contour) {
  if (I_original.getSize() == 0) {
    return;
  }

  vpImage<int> I(I_original.getHeight(), I_original.getWidth());
  for (unsigned int cpt = 0; cpt < I_original.getSize(); cpt++) {
    I.bitmap[cpt] = I_original.bitmap[cpt];
  }

  int nbd = 1; //newest border
  int lnbd = 1; //last newest border

  //Background contour
  //By default the root contour is a hole contour
  vpContour *root = new vpContour(vp::CONTOUR_HOLE);

  std::map<int, vpContour*> borderMap;
  borderMap[lnbd] = root;

  for (unsigned int i = 0; i < I.getHeight(); i++) {
    lnbd = 1;

    for (unsigned int j = 0; j < I.getWidth(); j++) {
      int fji = I[i][j];

      bool isOuter = isOuterBorderStart(I, i, j);
      bool isHole = isHoleBorderStart(I, i, j);

      if (isOuter || isHole) { //else 1(c)
        vpContour *border = new vpContour;
        vpContour *borderPrime = NULL;
        vpImagePoint from(i, j);

        if (isOuter) {
          nbd++;
          from.set_j(from.get_j() - 1);

          border->m_contourType = vp::CONTOUR_OUTER;
          if (borderMap.find(lnbd) != borderMap.end()) {
            borderPrime = borderMap[lnbd];

            switch (borderPrime->m_contourType) {
              case vp::CONTOUR_OUTER:
                border->setParent(borderPrime->m_parent);
                break;

              case vp::CONTOUR_HOLE:
                border->setParent(borderPrime);
                break;

              default:
                break;
            }
          }
        } else {
          nbd++;

          if (fji > 1) {
            lnbd = fji;
          }

          if (borderMap.find(lnbd) != borderMap.end()) {
            borderPrime = borderMap[lnbd];
            from.set_j(from.get_j() + 1);
            border->m_contourType = vp::CONTOUR_HOLE;

            switch (borderPrime->m_contourType) {
              case vp::CONTOUR_OUTER:
                border->setParent(borderPrime);
                break;

              case vp::CONTOUR_HOLE:
                border->setParent(borderPrime->m_parent);
                break;

              default:
                break;
            }
          }
        }

        if (borderPrime != NULL)
        {
          vpImagePoint ij(i, j);
          if (followBorder(I, ij, from, border, nbd)) {
            if (border->m_points.empty()) {
              border->m_points.push_back(ij);
              I[i][j] = -nbd;
            }

            borderMap[nbd] = border;
          } else {
            I[i][j] = -nbd;
            border->m_points.clear();

            if (isOuter) {
              switch (borderPrime->m_contourType) {
                case vp::CONTOUR_OUTER:
                  border->removeParent(borderPrime->m_parent);
                  break;

                case vp::CONTOUR_HOLE:
                  border->removeParent(borderPrime);
                  break;

                default:
                  break;
              }
            }

            if (isHole) {
              switch (borderPrime->m_contourType) {
                case vp::CONTOUR_OUTER:
                  border->removeParent(borderPrime);
                  break;

                case vp::CONTOUR_HOLE:
                  border->removeParent(borderPrime->m_parent);
                  break;

                default:
                  break;
              }
            }

            delete border;
          }
        }
      }

      if (fji != 0 && fji != 1) {
       lnbd = std::abs(fji);
      }
    }
  }

  contour = *root;

  delete root;
  root = NULL;
}
