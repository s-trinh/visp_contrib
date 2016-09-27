/**
 * Copyright (c) 2011, The University of Southampton and the individual contributors.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 *
 *   * 	Redistributions of source code must retain the above copyright notice,
 * 	this list of conditions and the following disclaimer.
 *
 *   *	Redistributions in binary form must reproduce the above copyright notice,
 * 	this list of conditions and the following disclaimer in the documentation
 * 	and/or other materials provided with the distribution.
 *
 *   *	Neither the name of the University of Southampton nor the names of its
 * 	contributors may be used to endorse or promote products derived from this
 * 	software without specific prior written permission.
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
 */
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
 * Basic contours extraction.
 *
 * Authors:
 * Souriya Trinh
 *
 *****************************************************************************/



/*!
  \file vpContours.h
  \brief Basic contours extraction.

*/

#ifndef __vpContours_h__
#define __vpContours_h__

#include <visp3/core/vpImage.h>


namespace vp
{
  typedef enum {
    OUTER_CONTOUR,
    HOLE_CONTOUR,
    BACKGROUND_CONTOUR
  } vpContourType;

  typedef enum {
    NORTH, NORTH_EAST, EAST, SOUTH_EAST, SOUTH, SOUTH_WEST, WEST, NORTH_WEST, LAST_DIRECTION
  } vpDirectionType;

  struct vpContour {
    vpContourType m_contourType;
    std::vector<vpContour *> m_children;
    vpContour *m_parent;
    std::vector<vpImagePoint> m_points;

    vpContour() : m_contourType(vp::HOLE_CONTOUR), m_children(), m_parent(NULL), m_points() {
    }

    vpContour(const vpContourType &type) : m_contourType(type) {
    }

    void setParent(vpContour *parent) {
      m_parent = parent;
      parent->m_children.push_back(this);
    }
  };


  struct vpDirection {
    vpDirectionType m_direction;

    int dirx[8] = { 0, 1, 1, 1, 0, -1, -1, -1 };
    int diry[8] = { -1, -1, 0, 1, 1, 1, 0, -1 };

//    static Direction[] entry = new Direction[] {
//      WEST, WEST, NORTH, NORTH, EAST, EAST, SOUTH, SOUTH
//    };
//    static Direction[] ccentry = new Direction[] {
//      EAST, SOUTH, SOUTH, WEST, WEST, NORTH, NORTH, EAST
//    };

    vpDirection clockwise() {
      vpDirection direction;
      int directionSize = (int) LAST_DIRECTION;
      direction.m_direction = vpDirectionType ( ( (int) m_direction + 1) % directionSize );

      return direction;
    }

    vpDirection counterClockwise() {
      vpDirection direction;
      int directionSize = (int) LAST_DIRECTION;
      int idx = (int) m_direction - 1;
      idx = idx < 0 ? idx + directionSize : idx;
      direction.m_direction = vpDirectionType ( idx );

      return direction;
    }

    vpImagePoint active(const vpImage<int> &I, const vpImagePoint &point) {
      unsigned int yy = (unsigned int) (point.get_i() + diry[(int) m_direction]);
      unsigned int xx = (unsigned int) (point.get_j() + dirx[(int) m_direction]);

      if (xx < 0 || xx >= I.getWidth() || yy < 0 || yy >= I.getHeight()) {
        return vpImagePoint(-1, -1);
      }

      unsigned char pixel = I[yy][xx];
      return pixel != 0 ? vpImagePoint(yy, xx) : vpImagePoint(-1, -1);
    }
  };



  VISP_EXPORT void extractContours(const vpImage<unsigned char> &I_original, vpContour *contour);
}

#endif
