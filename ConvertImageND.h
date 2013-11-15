#ifndef __ConvertImageND_h_
#define __ConvertImageND_h_

#include "itkOrientedRASImage.h"
#include "itkImageRegionIteratorWithIndex.h"
#include "itkInterpolateImageFunction.h"
#include "ConvertException.h"
#include "itkMetaDataDictionary.h"

#include <iostream>
#include <cctype>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <cstdlib>

using namespace std;

template <class TPixel, unsigned int VDim> class ConvertAdapter;

template<class TPixel, unsigned int VDim>
class ImageConverter
{
public:

  // Image typedef
  typedef itk::OrientedRASImage<TPixel, VDim> ImageType;
  typedef itk::Image<TPixel, VDim> UnorientedImageType;
  typedef typename ImageType::Pointer ImagePointer;
  typedef typename ImageType::SizeType SizeType;
  typedef typename ImageType::IndexType IndexType;
  typedef typename ImageType::RegionType RegionType;
  typedef vnl_vector_fixed<double, VDim> RealVector;
  typedef vnl_vector_fixed<int, VDim> IntegerVector;
  typedef std::map<TPixel, vnl_vector_fixed<double, 4> > LabelToRGBAMap;

  // Complex stuff
  typedef std::complex<TPixel> ComplexPixel;
  typedef itk::OrientedRASImage<ComplexPixel, VDim> ComplexImageType;

  // Iterators
  typedef itk::ImageRegionIteratorWithIndex<ImageType> Iterator;
  typedef itk::ImageRegionConstIteratorWithIndex<ImageType> ConstIterator;

  // Interpolator
  typedef itk::InterpolateImageFunction<ImageType, double> Interpolator;
  
  ImageConverter();
  int ProcessCommandLine(int argc, char *argv[]);

  friend class ConvertAdapter<TPixel, VDim>;

  // Copy image on stack
  void CopyImage();

  // Get bounding box of an image
  void GetBoundingBox(ImageType *image, RealVector &bb0, RealVector &bb1);

  // Check if N images on the stack have the same dimensions as the top image
  // (pass zero to check all images on the stack)
  bool CheckStackSameDimensions(size_t n);

  // Read label to RGBA mapping from file (SNAP format)
  static LabelToRGBAMap ReadLabelToRGBAMap(const char *fname);

  // Print a matrix in a nice way
  void PrintMatrix(
    std::ostream &sout, vnl_matrix<double> mat, 
    const char *fmt = "%10.4f ", const char *prefix = "");

  // Set label set for split/merge operations
  typedef std::vector<double> LabelSet;
  LabelSet &GetSplitLabelSet() 
    { return m_SplitLabelSet; }

  Interpolator *GetInterpolator() const
    { return m_Interpolator; }

  void SetInterpolator(Interpolator *interp)
    { m_Interpolator = interp; }

  // Read vectors, etc from command line
  SizeType ReadSizeVector(const char *vec);
  IndexType ReadIndexVector(const char *vec);
  RealVector ReadRealVector(const char *vec, bool is_point);
  RealVector ReadRealSize(const char *vec);
  TPixel ReadIntensityValue(const char *vec);

  static void PrintCommandListing(std::ostream &out);

private:

  // Internal functions
  int ProcessCommand(int argc, char *argv[]);


  // Templated write function
  template<class TOutPixel> void TemplatedWriteImage(const char *file, double xRoundFactor);

  // Map of variable names
  typedef map<string, ImagePointer> ImageVariableMap;
  ImageVariableMap m_ImageVars;

  // Image interpolator for all interpolation commands
  itk::SmartPointer<Interpolator> m_Interpolator;

  // Implementation of the 'foreach' loop
  size_t ForEachLoop(int argc, char *argv[]);

  // Implementation of the 'accum' loop
  size_t AccumulateLoop(int argc, char *argv[]);

  // Type of loop we are in currently
  enum LoopType { LOOP_NONE = 0, LOOP_FOREACH, LOOP_ACCUM };
  LoopType m_LoopType;

  // Label set for split/merge
  LabelSet m_SplitLabelSet;

public:

  // Stack of images from the command line
  vector<ImagePointer> m_ImageStack;

  // Typeid of the image to be saved
  string m_TypeId;

  // Interpolation type
  string m_Interpolation;

  // Background intensity
  double m_Background;

  // Rounding factor (not used for float and double) is 0.5 unless -noround was specified
  double m_RoundFactor;

  // Whether SPM extensions are used
  bool m_FlagSPM;

  // Whether multicomponent images are split on read
  bool m_MultiComponentSplit;

  // Number of iterations for various algorithms
  size_t m_Iterations;

  // Root mean square error for anti-aliasing algorithm
  double m_AntiAliasRMS;

  // Level set algorithm parameters
  double m_LevSetCurvature, m_LevSetAdvection;

  // N3 and N4 filter parameters
  // distance (in mm) of the mesh resolution at the base level
  std::vector<double>  n4_spline_distance ; 
  // N3/N4 shrink factor
  int     n4_shrink_factor;
  // N3/N4 spline order
  int     n4_spline_order;
  // N3/N4 number of histogram bins
  int     n4_histogram_bins;
  // N3/N4 histogram blurring FWHM
  double  n4_fwhm;
  // N3/N4 convergence threshold
  double  n4_convergence_threshold;
  // N3/N4 weiner noise estimate for regularization
  double  n4_weiner_noise;
  // N3/N4 max number of itearation
  int     n4_max_iterations;
  // N4 use optimal sclaling
  bool    n4_optimal_scaling;
  // N3/N4 output bias field
  bool    n4_output_field;
  // N3/N4 use mask to estimate bias
  bool    n4_use_mask;
  // image metadata to be preserved
  itk::MetaDataDictionary m_Metadata;
  // command line history info
  std::string m_History;
  
  // How % is handled for intensity specs
  enum PercentIntensityMode { PIM_QUANTILE, PIM_FGQUANTILE, PIM_RANGE };
  PercentIntensityMode m_PercentIntensityMode;

  // Verbose output stream
  std::ostringstream devnull;
  std::ostream *verbose;
};

#endif

