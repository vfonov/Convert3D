#include "BiasFieldCorrectionN4.h"
#include "itkBSplineControlPointImageFilter.h"
//#include "itkBSplineControlPointImageFilter.h"
#include "itkConstantPadImageFilter.h"
#include "itkDivideImageFilter.h"
#include "itkExpImageFilter.h"
#include "itkExtractImageFilter.h"
//#include "itkN4MRIBiasFieldCorrectionImageFilter.h"

#include "itkN4BiasFieldCorrectionImageFilter.h"
#include "itkOtsuThresholdImageFilter.h"
#include "itkShrinkImageFilter.h"

template <class TPixel, unsigned int VDim>
void
BiasFieldCorrectionN4<TPixel, VDim>
::operator() ()
{
  // distance (in mm) of the mesh resolution at the base level
  std::vector<double>  n4_spline_distance = c->n4_spline_distance;
  int     n4_shrink_factor=c->n4_shrink_factor;
  int     n4_spline_order=c->n4_spline_order;
  int     n4_histogram_bins= c->n4_histogram_bins;
  double  n4_fwhm=c->n4_fwhm;
  double  n4_convergence_threshold=c->n4_convergence_threshold;
  double  n4_weiner_noise=c->n4_weiner_noise;
  int     n4_max_iterations=c->n4_max_iterations;
  bool    n4_optimal_scaling=c->n4_optimal_scaling;
  bool    n4_output_field=c->n4_output_field;
  bool    n4_use_mask=c->n4_use_mask;
  
  ImagePointer mri;
  ImagePointer mask;
  
  // Check input availability
  if(n4_use_mask)
  {
    if(c->m_ImageStack.size() < 2)
      throw ConvertException("No image and mask on stack");

    // Get image from stack
    mri = c->m_ImageStack.back();
    c->m_ImageStack.pop_back();

    mask = c->m_ImageStack.back();
    c->m_ImageStack.pop_back();

  } else {
    if(c->m_ImageStack.size() < 1)
      throw ConvertException("No images on stack");

    // Get image from stack
    mri = c->m_ImageStack.back();
    c->m_ImageStack.pop_back();
  }
  
  // Bias filter
  typedef itk::N4BiasFieldCorrectionImageFilter<ImageType, ImageType, ImageType> CorrecterType;
  typename CorrecterType::Pointer correcter = CorrecterType::New();

  *c->verbose << "N4 BiasFieldCorrection #" << c->m_ImageStack.size() << endl;
  *c->verbose << "  Shrink factor: " << n4_shrink_factor << endl;

  *c->verbose << "  Spline distance: ";

  for(int i=0;i<n4_spline_distance.size();i++) 
    *c->verbose << n4_spline_distance[i]<<" ";
  *c->verbose << std::endl;

  *c->verbose << "  Number of histogram bins: "<< n4_histogram_bins<<endl;
  *c->verbose << "  Weiner Filter noise: "<< n4_weiner_noise<<endl;
  *c->verbose << "  Bias FWHM: "<< n4_fwhm<<endl;
  *c->verbose << "  Max Number of Iterations: "<< n4_max_iterations<<endl;
  *c->verbose << "  Convergence Threshold: "<< n4_convergence_threshold<<endl;
  *c->verbose << "  Spline Order: "<< n4_spline_order<<endl;
  *c->verbose << "  use mask: "<< n4_use_mask<<endl;


  typename CorrecterType::ArrayType numberOfControlPoints;

  typename ImageType::IndexType inputImageIndex =
    mri->GetLargestPossibleRegion().GetIndex();
  typename ImageType::SizeType inputImageSize =
    mri->GetLargestPossibleRegion().GetSize();

  typename ImageType::PointType newOrigin = mri->GetOrigin();

  unsigned long lowerBound[VDim];
  unsigned long upperBound[VDim];

  // if spline distance doesn't have enough elements, just keep using the last one
  if(n4_spline_distance.size()<VDim)
    n4_spline_distance.resize(VDim,n4_spline_distance.back()); 

  for( unsigned int d = 0; d < VDim; d++ )
  {
    float domain = static_cast<float>( mri->
      GetLargestPossibleRegion().GetSize()[d] - 1 ) * mri->GetSpacing()[d];

    unsigned int numberOfSpans = static_cast<unsigned int>( vcl_ceil( domain / n4_spline_distance[d] ) );

    unsigned long extraPadding = static_cast<unsigned long>( ( numberOfSpans *
      n4_spline_distance[d] - domain ) / mri->GetSpacing()[d] + 0.5 );

    lowerBound[d] = static_cast<unsigned long>( 0.5 * extraPadding );
    upperBound[d] = extraPadding - lowerBound[d];
    newOrigin[d] -= ( static_cast<float>( lowerBound[d] ) * mri->GetSpacing()[d] );

    numberOfControlPoints[d] = numberOfSpans + correcter->GetSplineOrder();
  }
  correcter->SetNumberOfControlPoints( numberOfControlPoints );

  typedef itk::ConstantPadImageFilter<ImageType, ImageType> PadderType;
  typename PadderType::Pointer padder = PadderType::New();
  padder->SetInput( mri );
  padder->SetPadLowerBound( lowerBound );
  padder->SetPadUpperBound( upperBound );
  padder->SetConstant( 0 );
  padder->Update();

  // Set up a filter to shrink image by a factor
  typedef itk::ShrinkImageFilter<ImageType, ImageType> ShrinkerType;
  typename ShrinkerType::Pointer shrinker = ShrinkerType::New();
  shrinker->SetInput( padder->GetOutput() );
  shrinker->SetShrinkFactors( n4_shrink_factor );
  shrinker->Update();

  typedef itk::OtsuThresholdImageFilter<ImageType, ImageType> ThresholderType;
  typename ThresholderType::Pointer otsu = ThresholderType::New();

  if( !n4_use_mask )
  {
    // Compute mask using Otsu threshold
    otsu->SetInput( mri );
    otsu->SetNumberOfHistogramBins( 200 );
    otsu->SetInsideValue( 0 );
    otsu->SetOutsideValue( 1 );
    otsu->Update();
    mask = otsu->GetOutput();
    *c->verbose << "  Otsu threshold: "<<otsu->GetThreshold()<<endl;
  }

  typedef itk::ConstantPadImageFilter<ImageType, ImageType> MaskPadderType;
  typename MaskPadderType::Pointer maskPadder = MaskPadderType::New();

  maskPadder->SetInput( mask );
  maskPadder->SetPadLowerBound( lowerBound );
  maskPadder->SetPadUpperBound( upperBound );
  maskPadder->SetConstant( 0 );
  maskPadder->Update();

  // Shrink the mask
  typename ShrinkerType::Pointer maskshrinker = ShrinkerType::New();
  maskshrinker->SetInput( maskPadder->GetOutput() );
  maskshrinker->SetShrinkFactors( n4_shrink_factor );
  maskshrinker->Update();

  correcter->SetInput( shrinker->GetOutput() );
  correcter->SetMaskImage( maskshrinker->GetOutput() );

  // These parameters are pretty standard
  correcter->SetSplineOrder(  n4_spline_order );
  correcter->SetNumberOfHistogramBins( n4_histogram_bins );
  correcter->SetBiasFieldFullWidthAtHalfMaximum( n4_fwhm );
  correcter->SetConvergenceThreshold( n4_convergence_threshold );
  correcter->SetWienerFilterNoise( n4_weiner_noise );

  // You will probably want to have an option for the maximum number of
  //  iterations at each level, the shrink factor, and the spline distance.
  typename CorrecterType::ArrayType numberOfFittingLevels;
  numberOfFittingLevels.Fill( 3 );
  correcter->SetNumberOfFittingLevels( numberOfFittingLevels );
  typename CorrecterType::VariableSizeArrayType maximumNumberOfIterations;
  maximumNumberOfIterations.SetSize( 3 );
  maximumNumberOfIterations[0] = n4_max_iterations;
  maximumNumberOfIterations[1] = n4_max_iterations;
  maximumNumberOfIterations[2] = n4_max_iterations;

  correcter->SetMaximumNumberOfIterations( maximumNumberOfIterations );

  // Progress meter
/*  typedef  itk::CommandIterationUpdate<CorrecterType> CommandType;
  typename CommandType::Pointer observer = CommandType::New();
  correcter->AddObserver( itk::IterationEvent(), observer );*/
  correcter->Update();

  /**
   * Reconstruct the bias field at full image resolution.  Divide
   * the original input image by the bias field to get the final
   * corrected image.
   */
  typedef itk::BSplineControlPointImageFilter<typename
    CorrecterType::BiasFieldControlPointLatticeType, typename
    CorrecterType::ScalarImageType> BSplinerType;

  typename BSplinerType::Pointer bspliner = BSplinerType::New();
  bspliner->SetInput( correcter->GetLogBiasFieldControlPointLattice() );
  bspliner->SetSplineOrder( correcter->GetSplineOrder() );
  bspliner->SetSize( mri->GetLargestPossibleRegion().GetSize() );
  bspliner->SetOrigin( mri->GetOrigin() );
  bspliner->SetDirection( mri->GetDirection() );
  bspliner->SetSpacing( mri->GetSpacing() );
  bspliner->Update();

  typename ImageType::Pointer logField = ImageType::New();
  logField->SetOrigin(    bspliner->GetOutput()->GetOrigin() );
  logField->SetSpacing(   bspliner->GetOutput()->GetSpacing() );
  logField->SetRegions(   bspliner->GetOutput()->GetLargestPossibleRegion().GetSize() );
  logField->SetDirection( bspliner->GetOutput()->GetDirection() );
  logField->Allocate();

  itk::ImageRegionIterator<typename CorrecterType::ScalarImageType> ItB(
    bspliner->GetOutput(),
    bspliner->GetOutput()->GetLargestPossibleRegion() );

  itk::ImageRegionIterator<ImageType> ItF( logField,
    logField->GetLargestPossibleRegion() );

  for( ItB.GoToBegin(), ItF.GoToBegin(); !ItB.IsAtEnd(); ++ItB, ++ItF )
    {
    ItF.Set( ItB.Get()[0] );
    }

  typedef itk::ExpImageFilter<ImageType, ImageType> ExpFilterType;
  typename ExpFilterType::Pointer expFilter = ExpFilterType::New();
  expFilter->SetInput( logField );
  expFilter->Update();

  typedef itk::DivideImageFilter<ImageType, ImageType, ImageType> DividerType;
  typename DividerType::Pointer divider = DividerType::New();
  divider->SetInput1( mri );
  divider->SetInput2( expFilter->GetOutput() );
  divider->Update();

  typename ImageType::RegionType inputRegion;
  inputRegion.SetIndex( inputImageIndex );
  inputRegion.SetSize(  inputImageSize );

  typedef itk::ExtractImageFilter<ImageType, ImageType> CropperType;
  typename CropperType::Pointer cropper = CropperType::New();
  cropper->SetInput( divider->GetOutput() );
  cropper->SetExtractionRegion( inputRegion );
  cropper->Update();

  typename CropperType::Pointer biasFieldCropper = CropperType::New();
  biasFieldCropper->SetInput( expFilter->GetOutput() );
  biasFieldCropper->SetExtractionRegion( inputRegion );
  biasFieldCropper->Update();
  
  if( n4_output_field )
  {
    c->m_ImageStack.push_back( biasFieldCropper->GetOutput() );
  } else {
    c->m_ImageStack.push_back( cropper->GetOutput() );
  }
}

// Invocations
template class BiasFieldCorrectionN4<double, 2>;
template class BiasFieldCorrectionN4<double, 3>;
