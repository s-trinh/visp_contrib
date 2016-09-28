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
 * Test connected components.
 *
 * Authors:
 * Souriya Trinh
 *
 *****************************************************************************/

#include <iomanip>

#include <visp3/core/vpIoTools.h>
#include <visp3/core/vpImageTools.h>
#include <visp3/io/vpImageIo.h>
#include <visp3/io/vpParseArgv.h>
#include <visp3/imgproc/vpImgproc.h>


/*!
  \example testConnectedComponents.cpp

  \brief Test connected components.
*/

// List of allowed command line options
#define GETOPTARGS  "cdi:o:h"

void usage(const char *name, const char *badparam, std::string ipath, std::string opath, std::string user);
bool getOptions(int argc, const char **argv, std::string &ipath, std::string &opath, std::string user);

/*
  Print the program options.

  \param name : Program name.
  \param badparam : Bad parameter name.
  \param ipath: Input image path.
  \param opath : Output image path.
  \param user : Username.
 */
void usage(const char *name, const char *badparam, std::string ipath, std::string opath, std::string user)
{
  fprintf(stdout, "\n\
Test connected components.\n\
\n\
SYNOPSIS\n\
  %s [-i <input image path>] [-o <output image path>]\n\
     [-h]\n                 \
", name);

  fprintf(stdout, "\n\
OPTIONS:                                               Default\n\
  -i <input image path>                                %s\n\
     Set image input path.\n\
     From this path read \"ViSP-images/Klimt/Klimt.pgm\"\n\
     image.\n\
     Setting the VISP_INPUT_IMAGE_PATH environment\n\
     variable produces the same behaviour than using\n\
     this option.\n\
\n\
  -o <output image path>                               %s\n\
     Set image output path.\n\
     From this directory, creates the \"%s\"\n\
     subdirectory depending on the username, where \n\
     output result images are written.\n\
\n\
  -h\n\
     Print the help.\n\n",
    ipath.c_str(), opath.c_str(), user.c_str());

  if (badparam)
    fprintf(stdout, "\nERROR: Bad parameter [%s]\n", badparam);
}

/*!
  Set the program options.

  \param argc : Command line number of parameters.
  \param argv : Array of command line parameters.
  \param ipath: Input image path.
  \param opath : Output image path.
  \param user : Username.
  \return false if the program has to be stopped, true otherwise.

*/
bool getOptions(int argc, const char **argv, std::string &ipath, std::string &opath, std::string user)
{
  const char *optarg_;
  int c;
  while ((c = vpParseArgv::parse(argc, argv, GETOPTARGS, &optarg_)) > 1) {

    switch (c) {
    case 'i': ipath = optarg_; break;
      case 'o': opath = optarg_; break;
      case 'h': usage(argv[0], NULL, ipath, opath, user); return false; break;

    case 'c':
    case 'd':
      break;

    default:
      usage(argv[0], optarg_, ipath, opath, user); return false; break;
    }
  }

  if ((c == 1) || (c == -1)) {
    // standalone param or error
    usage(argv[0], NULL, ipath, opath, user);
    std::cerr << "ERROR: " << std::endl;
    std::cerr << "  Bad argument " << optarg_ << std::endl << std::endl;
    return false;
  }

  return true;
}

void printImage(const vpImage<unsigned char> &I, const std::string &name) {
  std::cout << "\n" << name << ":" << std::endl;

  std::cout << "   ";
  for (unsigned int j = 0; j < I.getWidth(); j++) {
    std::cout << std::setfill(' ') << std::setw(2) << j << " ";
  }
  std::cout << std::endl;

  for (unsigned int i = 0; i < I.getHeight(); i++) {
    std::cout << std::setfill(' ') << std::setw(2) << i << " ";

    for (unsigned int j = 0; j < I.getWidth(); j++) {
      std::cout << std::setfill(' ') << std::setw(2) << static_cast<unsigned int>(I[i][j]) << " ";
    }

    std::cout << std::endl;
  }
}

void displayContourInfo(const vp::vpContour &contour, const int level, std::vector<std::vector<vpImagePoint> > &contours) {
  std::cout << "\nContour:" << std::endl;
  std::cout << "\tlevel: " << level << std::endl;
  std::cout << "\tcontour type: " << (contour.m_contourType == vp::CONTOUR_OUTER ? "outer contour" : "hole contour") << std::endl;
  std::cout << "\tnb children: " << contour.m_children.size() << std::endl;

  contours.push_back(contour.m_points);

  for (std::vector<vp::vpContour *>::const_iterator it = contour.m_children.begin(); it != contour.m_children.end(); ++it) {
    displayContourInfo(**it, level+1, contours);
  }
}

int
main(int argc, const char ** argv)
{
  try {
    std::string env_ipath;
    std::string opt_ipath;
    std::string opt_opath;
    std::string ipath;
    std::string opath;
    std::string filename;
    std::string username;

    // Get the visp-images-data package path or VISP_INPUT_IMAGE_PATH environment variable value
    env_ipath = vpIoTools::getViSPImagesDataPath();

    // Set the default input path
    if (! env_ipath.empty())
      ipath = env_ipath;

    // Set the default output path
#if defined(_WIN32)
    opt_opath = "C:/temp";
#else
    opt_opath = "/tmp";
#endif

    // Get the user login name
    vpIoTools::getUserName(username);

    // Read the command line options
    if (getOptions(argc, argv, opt_ipath, opt_opath, username) == false) {
      exit (EXIT_FAILURE);
    }

    // Get the option values
    if (!opt_ipath.empty())
      ipath = opt_ipath;
    if (!opt_opath.empty())
      opath = opt_opath;

    // Append to the output path string, the login name of the user
    opath = vpIoTools::createFilePath(opath, username);

    // Test if the output path exist. If no try to create it
    if (vpIoTools::checkDirectory(opath) == false) {
      try {
        // Create the dirname
        vpIoTools::makeDirectory(opath);
      }
      catch (...) {
        usage(argv[0], NULL, ipath, opt_opath, username);
        std::cerr << std::endl
                  << "ERROR:" << std::endl;
        std::cerr << "  Cannot create " << opath << std::endl;
        std::cerr << "  Check your -o " << opt_opath << " option " << std::endl;
        exit(EXIT_FAILURE);
      }
    }

    // Compare ipath and env_ipath. If they differ, we take into account
    // the input path comming from the command line option
    if (!opt_ipath.empty() && !env_ipath.empty()) {
      if (ipath != env_ipath) {
        std::cout << std::endl
                  << "WARNING: " << std::endl;
        std::cout << "  Since -i <visp image path=" << ipath << "> "
                  << "  is different from VISP_IMAGE_PATH=" << env_ipath << std::endl
                  << "  we skip the environment variable." << std::endl;
      }
    }

    // Test if an input path is set
    if (opt_ipath.empty() && env_ipath.empty()){
      usage(argv[0], NULL, ipath, opt_opath, username);
      std::cerr << std::endl
                << "ERROR:" << std::endl;
      std::cerr << "  Use -i <visp image path> option or set VISP_INPUT_IMAGE_PATH "
                << std::endl
                << "  environment variable to specify the location of the " << std::endl
                << "  image path where test images are located." << std::endl << std::endl;
      exit(EXIT_FAILURE);
    }


    //
    // Here starts really the test
    //

    unsigned char image_data[14*10] = {
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 1, 1, 0, 0, 1, 0, 0, 0,
      0, 1, 1, 1, 1, 1, 1, 0, 0, 0,
      0, 0, 0, 1, 1, 1, 1, 1, 0, 0,
      0, 0, 1, 0, 0, 1, 0, 0, 1, 0,
      0, 0, 1, 0, 0, 0, 1, 0, 0, 0,
      0, 0, 0, 1, 1, 1, 0, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 1, 1, 1, 1, 0, 0, 0, 0,
      0, 0, 1, 0, 0, 0, 1, 0, 0, 0,
      0, 0, 1, 0, 1, 1, 0, 1, 0, 0,
      0, 0, 0, 1, 0, 0, 1, 0, 0, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    };

#define TEST_DATA 0
#if TEST_DATA
    vpImage<unsigned char> I(image_data, 14, 10, true);
    printImage(I, "I");
#else
    //Read Klimt.ppm
    filename = vpIoTools::createFilePath(ipath, "ViSP-images/Klimt/Klimt.pgm");
    vpImage<unsigned char> I;
    std::cout << "Read image: " << filename << std::endl;
    vpImageIo::read(I, filename);
    vpImageTools::binarise(I, (unsigned char) 127, (unsigned char) 255, (unsigned char) 0, (unsigned char) 1, (unsigned char) 1);

    //Remove 1-pixel border
    for (unsigned int i = 0; i < I.getHeight(); i++) {
      if (i == 0 || i == I.getHeight()-1) {
        for (unsigned int j = 0; j < I.getWidth(); j++) {
          I[i][j] = 0;
        }
      } else {
        I[i][0] = 0;
        I[i][I.getWidth()-1] = 0;
      }
    }

    vpImage<unsigned char> I2(I.getHeight(), I.getWidth());
    for (unsigned int cpt = 0; cpt < I2.getSize(); cpt++) {
      I2.bitmap[cpt] = 255*I.bitmap[cpt];
    }
    filename = vpIoTools::createFilePath(opath, "Klimt_contours_binarise.pgm");
    vpImageIo::write(I2, filename);
#endif

    cv::Mat matImg;
    vpImageConvert::convert(I, matImg);

    std::vector<std::vector<cv::Point> > contours_opencv;
    double t_opencv = vpTime::measureTimeMs();
    cv::findContours(matImg, contours_opencv, cv::RETR_TREE, cv::CHAIN_APPROX_NONE);
    t_opencv = vpTime::measureTimeMs() - t_opencv;
    std::cout << "OpenCV: nb contours=" << contours_opencv.size() << " ; t_opencv=" << t_opencv << " ms" << std::endl;

    vpImage<unsigned char> I_drawContours_opencv(I.getHeight(), I.getWidth(), 0);
    for (std::vector<std::vector<cv::Point> >::const_iterator it1 = contours_opencv.begin(); it1 != contours_opencv.end(); ++it1) {
      for (std::vector<cv::Point>::const_iterator it2 = it1->begin(); it2 != it1->end(); ++it2) {
        I_drawContours_opencv[it2->y][it2->x] = 255;
      }
    }

    vp::vpContour vp_contours;
    double t = vpTime::measureTimeMs();
    vp::extractContours(I, vp_contours);
    t = vpTime::measureTimeMs() - t;

    std::vector<std::vector<vpImagePoint> > contours;
    displayContourInfo(vp_contours, 0, contours);
    std::cout << "ViSP: nb contours=" << contours.size() << " ; t=" << t << " ms" << std::endl;

    vpImage<unsigned char> I_drawContours(I.getHeight(), I.getWidth(), 0);
    for (std::vector<std::vector<vpImagePoint> >::const_iterator it1 = contours.begin(); it1 != contours.end(); ++it1) {
      for (std::vector<vpImagePoint>::const_iterator it2 = it1->begin(); it2 != it1->end(); ++it2) {
        unsigned int i = (unsigned int) it2->get_i();
        unsigned int j = (unsigned int) it2->get_j();
        I_drawContours[i][j] = 255;
      }
    }

    std::cout << "(I_drawContours_opencv == I_drawContours)? " << (I_drawContours_opencv == I_drawContours) << std::endl;

    filename = vpIoTools::createFilePath(opath, "Klimt_contours_extracted.pgm");
    vpImageIo::write(I_drawContours, filename);

    filename = vpIoTools::createFilePath(opath, "Klimt_contours_extracted_opencv.pgm");
    vpImageIo::write(I_drawContours_opencv, filename);


    return EXIT_SUCCESS;
  }
  catch(vpException &e) {
    std::cerr << "Catch an exception: " << e.what() << std::endl;
    return EXIT_FAILURE;
  }
}
