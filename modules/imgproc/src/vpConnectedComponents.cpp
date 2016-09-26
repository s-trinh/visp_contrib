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
 * Connected components.
 *
 * Authors:
 * Souriya Trinh
 *
 *****************************************************************************/

/*!
  \file vpConnectedComponents.cpp
  \brief Basic connected components.
*/

#include <map>
#include <set>
#include <queue>
#include <visp3/imgproc/vpImgproc.h>


void getNeighbors(const vpImage<unsigned char> &I, std::queue<vpImagePoint> &listOfNeighbors,
                  const unsigned int i, const unsigned int j,
                  const vp::vpConnectedConnexityType &connexity) {
  unsigned char currValue = I[i][j];

  if (connexity == vp::CONNECTED_CONNEXITY_4) {
    //Top
    if (I[i-1][j] == currValue) {
      listOfNeighbors.push(vpImagePoint(i-1, j));
    }

    //Left
    if (I[i][j-1] == currValue) {
      listOfNeighbors.push(vpImagePoint(i, j-1));
    }

    //Right
    if (I[i][j+1] == currValue) {
      listOfNeighbors.push(vpImagePoint(i, j+1));
    }

    //Bottom
    if (I[i+1][j] == currValue) {
      listOfNeighbors.push(vpImagePoint(i+1, j));
    }
  } else {
    for (int cpt1 = -1; cpt1 <= 1; cpt1++) {
      for (int cpt2 = -1; cpt2 <= 1; cpt2++) {
        //Everything except the current position
        if (cpt1 != 0 || cpt2 != 0) {
          if ( I[(int) i + cpt1][(int) j + cpt2] == currValue ) {
            listOfNeighbors.push(vpImagePoint( (int) i + cpt1, (int) j + cpt2));
          }
        }
      }
    }
  }
}

void visitNeighbors(vpImage<unsigned char> &I_copy, std::queue<vpImagePoint> &listOfNeighbors, vpImage<int> &labels,
                    const int current_label, const vp::vpConnectedConnexityType &connexity) {
  //Visit the neighbors
  while (!listOfNeighbors.empty()) {
    vpImagePoint imPt = listOfNeighbors.front();
    unsigned int i = (unsigned int) imPt.get_i();
    unsigned int j = (unsigned int) imPt.get_j();
    listOfNeighbors.pop();

    if (I_copy[i][j]) {
      getNeighbors(I_copy, listOfNeighbors, i, j, connexity);

      //Reset current position and set label
      I_copy[i][j] = 0;
      labels[i][j] = current_label;
    }
  }
}

/*!
  \ingroup group_imgproc_connected_components

  Perform connected components detection.

  \param I : Input image (0 means background).
  \param labels : Label image that contain for each position the component label.
  \param nbComponents : Number of connected components.
  \param connexity : Type of connexity.
*/
void vp::connectedComponents(const vpImage<unsigned char> &I, vpImage<int> &labels,
                             int &nbComponents, const vpConnectedConnexityType &connexity) {
  if (I.getSize() == 0) {
    return;
  }

  labels.resize(I.getHeight(), I.getWidth());

  vpImage<unsigned char> I_copy(I.getHeight()+2, I.getWidth()+2);
  // Copy and add border
  for (unsigned int i = 0; i < I_copy.getHeight(); i++) {
    if (i == 0 || i == I_copy.getHeight() - 1) {
      for (unsigned int j = 0; j < I_copy.getWidth(); j++) {
        I_copy[i][j] = 0;
      }
    } else {
      I_copy[i][0] = 0;
      memcpy(I_copy[i]+1, I[i-1], sizeof(unsigned char)*I.getWidth());
      I_copy[i][I_copy.getWidth() - 1] = 0;
    }
  }

  vpImage<int> labels_copy(I.getHeight()+2, I.getWidth()+2, 0);

  int current_label = 1;
  std::queue<vpImagePoint> listOfNeighbors;

  for (unsigned int cpt1 = 0; cpt1 < I.getHeight(); cpt1++) {
    unsigned int i = cpt1+1;

    for (unsigned int cpt2 = 0; cpt2 < I.getWidth(); cpt2++) {
      unsigned int j = cpt2+1;

      if (I_copy[i][j] && labels_copy[i][j] == 0) {
        //Get all the neighbors relative to the current position
        getNeighbors(I_copy, listOfNeighbors, i, j, connexity);

        //Reset current position and set label
        I_copy[i][j] = 0;
        labels_copy[i][j] = current_label;

        visitNeighbors(I_copy, listOfNeighbors, labels_copy, current_label, connexity);

        //Increment label
        current_label++;
      }
    }
  }

  for (unsigned int i = 0; i < labels.getHeight(); i++) {
    memcpy(labels[i], labels_copy[i+1]+1, sizeof(int)*labels.getWidth());
  }

  nbComponents = current_label - 1;
}

std::set<int> getNeighborLabels(const vpImage<unsigned char> &I, const vpImage<int> &labels,
                       const unsigned int i, const unsigned int j,
                       const vp::vpConnectedConnexityType &connexity) {
  unsigned char currValue = I[i][j];
  std::set<int> setOfNeighborLabels;

  if (connexity == vp::CONNECTED_CONNEXITY_4) {
    //Top
    if (I[i-1][j] == currValue && labels[i-1][j] != 0) {
      setOfNeighborLabels.insert( labels[i-1][j] );
    }

    //Left
    if (I[i][j-1] == currValue && labels[i][j-1] != 0) {
      setOfNeighborLabels.insert( labels[i][j-1] );
    }

  } else {
    //Top
    for (int cpt2 = -1; cpt2 <= 1; cpt2++) {
      if ( I[i-1][(int) j + cpt2] == currValue && labels[i-1][(int) j + cpt2] != 0 ) {
        setOfNeighborLabels.insert( labels[i-1][(int) j + cpt2] );
      }
    }

    //Left
    if (I[i][j-1] == currValue && labels[i][j-1] != 0) {
      setOfNeighborLabels.insert( labels[i][j-1] );
    }
  }

  return setOfNeighborLabels;
}

/*!
  \ingroup group_imgproc_connected_components

  Perform connected components detection.

  \param I : Input image (0 means background).
  \param labels : Label image that contain for each position the component label.
  \param nbComponents : Number of connected components.
  \param connexity : Type of connexity.
*/
void vp::connectedComponents2(const vpImage<unsigned char> &I, vpImage<int> &labels,
                             int &nbComponents, const vpConnectedConnexityType &connexity) {
  if (I.getSize() == 0) {
    return;
  }

  labels.resize(I.getHeight(), I.getWidth());

  vpImage<unsigned char> I_copy(I.getHeight()+2, I.getWidth()+2);
  // Copy and add border
  for (unsigned int i = 0; i < I_copy.getHeight(); i++) {
    if (i == 0 || i == I_copy.getHeight() - 1) {
      for (unsigned int j = 0; j < I_copy.getWidth(); j++) {
        I_copy[i][j] = 0;
      }
    } else {
      I_copy[i][0] = 0;
      memcpy(I_copy[i]+1, I[i-1], sizeof(unsigned char)*I.getWidth());
      I_copy[i][I_copy.getWidth() - 1] = 0;
    }
  }

  vpImage<int> labels_copy(I.getHeight()+2, I.getWidth()+2, 0);

  int current_label = 1;
  std::map<int, std::set<int> > equivalent_labels;

  //First pass
  for (unsigned int cpt1 = 0; cpt1 < I.getHeight(); cpt1++) {
    unsigned int i = cpt1+1;

    for (unsigned int cpt2 = 0; cpt2 < I.getWidth(); cpt2++) {
      unsigned int j = cpt2+1;

      if (I_copy[i][j]) {
        //Get neighbor labels
        std::set<int> neighbor_labels;
        neighbor_labels = getNeighborLabels(I_copy, labels_copy, i, j, connexity);

        if (neighbor_labels.empty()) {
          equivalent_labels[current_label].insert(current_label);
          labels_copy[i][j] = current_label;
          current_label++;
        } else {
          //Get the smallest label
          int smallest_label = *neighbor_labels.begin();
          labels_copy[i][j] = smallest_label;
          for (std::set<int>::const_iterator it = neighbor_labels.begin(); it != neighbor_labels.end(); ++it) {
            equivalent_labels[*it].insert(neighbor_labels.begin(), neighbor_labels.end());
          }

//          //Create same set?
//          for (std::set<int>::const_iterator it1 = neighbor_labels.begin(); it1 != neighbor_labels.end(); ++it1) {
//            for (std::set<int>::const_iterator it2 = neighbor_labels.begin(); it2 != neighbor_labels.end(); ++it2) {
//              if (*it1 != *it2) {
//                equivalent_labels[*it1].insert(equivalent_labels[*it2].begin(), equivalent_labels[*it2].end());
//              }
//            }
//          }
        }
      }
    }
  }

  for (std::map<int, std::set<int> >::iterator it1 = equivalent_labels.begin(); it1 != equivalent_labels.end(); ++it1) {
    std::set<int> current_set = it1->second;

    for (std::set<int>::iterator it2 = current_set.begin(); it2 != current_set.end(); ++it2) {
      if (it1->first != *it2) {
        std::set_union( it1->second.begin(), it1->second.end(),
                        equivalent_labels[*it2].begin(), equivalent_labels[*it2].end(),
                        std::inserter(it1->second, it1->second.begin()) );
      }
    }
  }

  std::cout << std::endl;
  for (unsigned int i = 1; i < labels_copy.getHeight()-1; i++) {
    for (unsigned int j = 1; j < labels_copy.getWidth()-1; j++) {
      std::cout << labels_copy[i][j] << " ";
    }

    std::cout << std::endl;
  }

  std::cout << std::endl;
  for (std::map<int, std::set<int> >::const_iterator it = equivalent_labels.begin(); it != equivalent_labels.end(); ++it) {
    std::cout << it->first << " is equivalent to ";
    for (std::set<int>::const_iterator it2 = it->second.begin(); it2 != it->second.end(); ++it2) {
      std::cout << *it2 << ", ";
    }

    std::cout << std::endl;
  }

  //Second pass
  for (unsigned int cpt1 = 0; cpt1 < I.getHeight(); cpt1++) {
    unsigned int i = cpt1+1;

    for (unsigned int cpt2 = 0; cpt2 < I.getWidth(); cpt2++) {
      unsigned int j = cpt2+1;

      if (I_copy[i][j]) {
        labels_copy[i][j] = *equivalent_labels[labels_copy[i][j]].begin();
      }
    }
  }

  for (unsigned int i = 0; i < labels.getHeight(); i++) {
    memcpy(labels[i], labels_copy[i+1]+1, sizeof(int)*labels.getWidth());
  }

  nbComponents = current_label - 1;
}
