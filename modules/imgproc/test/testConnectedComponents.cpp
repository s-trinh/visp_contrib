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

void printMatrix(const vpImage<int> &I, const std::string &name) {
  std::cout << "\n" << name << ":" << std::endl;
  for(unsigned int i = 0; i < I.getHeight(); i++) {
    for(unsigned int j = 0; j < I.getWidth(); j++) {
      std::cout << I[i][j] << " ";
    }
    std::cout << std::endl;
  }
}

void printMatrix(const vpImage<unsigned char> &I, const std::string &name) {
  std::cout << "\n" << name << ":" << std::endl;
  for(unsigned int i = 0; i < I.getHeight(); i++) {
    for(unsigned int j = 0; j < I.getWidth(); j++) {
      std::cout << static_cast<unsigned int>(I[i][j]) << " ";
    }
    std::cout << std::endl;
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

    //Read Klimt.ppm
    filename = vpIoTools::createFilePath(ipath, "ViSP-images/Klimt/Klimt.pgm");
    vpImage<unsigned char> I;
    std::cout << "Read image: " << filename << std::endl;
    vpImageIo::read(I, filename);
    vpImageTools::binarise(I, (unsigned char) 127, (unsigned char) 255, (unsigned char) 0, (unsigned char) 255, (unsigned char) 255);
    std::cout << "Image: " << I.getWidth() << "x" << I.getHeight() << std::endl;

    vpImage<int> labels_connex4;
    int nbComponents = 0;
    double t = vpTime::measureTimeMs();
    vp::connectedComponents(I, labels_connex4, nbComponents, vp::CONNECTED_CONNEXITY_4);
    t = vpTime::measureTimeMs() - t;
    std::cout << "\n4-connexity connected components:" << std::endl;
    std::cout << "Time: " << t << " ms" << std::endl;
    std::cout << "nbComponents=" << nbComponents << std::endl;

    vpImage<int> labels_connex8;
    t = vpTime::measureTimeMs();
    vp::connectedComponents(I, labels_connex8, nbComponents, vp::CONNECTED_CONNEXITY_8);
    t = vpTime::measureTimeMs() - t;
    std::cout << "\n8-connexity connected components:" << std::endl;
    std::cout << "Time: " << t << " ms" << std::endl;
    std::cout << "nbComponents=" << nbComponents << std::endl;


    //Save results
    vpImage<vpRGBa> labels_connex4_color(labels_connex4.getHeight(), labels_connex4.getWidth(), vpRGBa(0,0,0,0));
    for (unsigned int i = 0; i < labels_connex4.getHeight(); i++) {
      for (unsigned int j = 0; j < labels_connex4.getWidth(); j++) {
        if (labels_connex4[i][j] != 0) {
          labels_connex4_color[i][j] = vpRGBa(vpColor::getColor( (unsigned int) labels_connex4[i][j]).R,
                                              vpColor::getColor( (unsigned int) labels_connex4[i][j]).G,
                                              vpColor::getColor( (unsigned int) labels_connex4[i][j]).B);
        }
      }
    }

    filename = vpIoTools::createFilePath(opath, "Klimt_connected_components_4.ppm");
    vpImageIo::write(labels_connex4_color, filename);

    vpImage<vpRGBa> labels_connex8_color(labels_connex8.getHeight(), labels_connex8.getWidth(), vpRGBa(0,0,0,0));
    for (unsigned int i = 0; i < labels_connex8.getHeight(); i++) {
      for (unsigned int j = 0; j < labels_connex8.getWidth(); j++) {
        if (labels_connex8[i][j] != 0) {
          labels_connex8_color[i][j] = vpRGBa(vpColor::getColor( (unsigned int) labels_connex8[i][j]).R,
                                              vpColor::getColor( (unsigned int) labels_connex8[i][j]).G,
                                              vpColor::getColor( (unsigned int) labels_connex8[i][j]).B);
        }
      }
    }

    filename = vpIoTools::createFilePath(opath, "Klimt_connected_components_8.ppm");
    vpImageIo::write(labels_connex8_color, filename);


#if (VISP_HAVE_OPENCV_VERSION >= 0x030000)
    cv::Mat matImg;
    vpImageConvert::convert(I, matImg);

    cv::Mat matLabels_4;
    double t_opencv = vpTime::measureTimeMs();
    cv::connectedComponents(matImg, matLabels_4, 4);
    t_opencv = vpTime::measureTimeMs() - t_opencv;

    vpImage<int> labels_connex4_opencv((unsigned int) matLabels_4.rows, (unsigned int) matLabels_4.cols);
    for (int i = 0; i < matLabels_4.rows; i++) {
      for (int j = 0; j < matLabels_4.cols; j++) {
        labels_connex4_opencv[i][j] = matLabels_4.at<int>(i, j);
      }
    }

    std::cout << "\n4-connexity connected components (OpenCV):" << std::endl;
    std::cout << "Time: " << t_opencv << " ms" << std::endl;
    std::cout << "(labels_connex4_opencv == labels_connex4)? " << (labels_connex4_opencv == labels_connex4) << std::endl;
    if (labels_connex4_opencv != labels_connex4) {
      throw vpException(vpException::fatalError, "(labels_connex4_opencv != labels_connex4)");
    }


    cv::Mat matLabels_8;
    t_opencv = vpTime::measureTimeMs();
    cv::connectedComponents(matImg, matLabels_8, 8);
    t_opencv = vpTime::measureTimeMs() - t_opencv;

    vpImage<int> labels_connex8_opencv((unsigned int) matLabels_8.rows, (unsigned int) matLabels_8.cols);
    for (int i = 0; i < matLabels_8.rows; i++) {
      for (int j = 0; j < matLabels_8.cols; j++) {
        labels_connex8_opencv[i][j] = matLabels_8.at<int>(i, j);
      }
    }

    std::cout << "\n8-connexity connected components (OpenCV):" << std::endl;
    std::cout << "Time: " << t_opencv << " ms" << std::endl;
    std::cout << "(labels_connex8_opencv == labels_connex8)? " << (labels_connex8_opencv == labels_connex8) << std::endl;

    if (labels_connex8_opencv != labels_connex8) {
      throw vpException(vpException::fatalError, "(labels_connex8_opencv != labels_connex8)");
    }
#endif


#define TEST_DATA 0

#if TEST_DATA
    unsigned char image_data[9*17] = {
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
      0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0,
      0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0, 0,
      0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0,
      0, 0, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0,
      0, 1, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0,
      0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0,
      0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1, 0,
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    };

    vpImage<unsigned char> I_test(image_data, 9, 17, true);
#else
    vpImage<unsigned char> I_test = I;
#endif

    vpImage<int> labels_test_4;
    int nbComponents_4;
    double t2 = vpTime::measureTimeMs();
    vp::connectedComponents2(I_test, labels_test_4, nbComponents_4, vp::CONNECTED_CONNEXITY_4);
//    vp::connectedComponents2(I_test, labels_test_4, nbComponents_4, vp::CONNECTED_CONNEXITY_8);
    t2 = vpTime::measureTimeMs() - t2;
    std::cout << "t2=" << t2 << " ms" << std::endl;
    std::cout << "nbComponents_4=" << nbComponents_4 << std::endl;

#if TEST_DATA
    printMatrix(labels_test_4, "labels_test_4");
#else
    std::cout << "(labels_test_4 == labels_connex4)? " << (labels_test_4 == labels_connex4) << std::endl;
#endif


    return EXIT_SUCCESS;
  }
  catch(vpException &e) {
    std::cerr << "Catch an exception: " << e.what() << std::endl;
    return EXIT_FAILURE;
  }
}
