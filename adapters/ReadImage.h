/*=========================================================================

  Program:   C3D: Command-line companion tool to ITK-SNAP
  Module:    ReadImage.h
  Language:  C++
  Website:   itksnap.org/c3d
  Copyright (c) 2014 Paul A. Yushkevich
  
  This file is part of C3D, a command-line companion tool to ITK-SNAP

  C3D is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
 
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.

=========================================================================*/

#ifndef __ReadImage_h_
#define __ReadImage_h_

#include "ConvertAdapter.h"

template<class TPixel, unsigned int VDim>
class ReadImage : public ConvertAdapter<TPixel, VDim>
{
public:
  // Common typedefs
  CONVERTER_STANDARD_TYPEDEFS

  // Structure for additional information for reading images
  struct ImageInfo
  {
    char *dicom_series_id;
    ImageInfo() : dicom_series_id(NULL) {}
  };

  ReadImage(Converter *c) : c(c) {}

  /**
   * Second parameter is extra information for special IO needs, like DICOM
   */
  void operator() (const char *file, const ImageInfo &info = ImageInfo());

private:
  Converter *c;

};

#endif

