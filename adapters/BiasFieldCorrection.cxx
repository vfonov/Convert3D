#include "BiasFieldCorrection.h"
#include "itkBSplineControlPointImageFilter.h"
#include "itkShrinkImageFilter.h"
#include "itkOtsuThresholdImageFilter.h"
#include "itkExpImageFilter.h"
#include "itkDivideImageFilter.h"
#include "itkN3MRIBiasFieldCorrectionImageFilter.h"

template <class TPixel, unsigned int VDim>
void
BiasFieldCorrection<TPixel, VDim>
::operator() ()
{
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

  // Report what the filter is doing
  *c->verbose << "N3 BiasFieldCorrection #" << c->m_ImageStack.size() << endl;
  *c->verbose << "  Spline distance: ";
  for(int i=0;i<n4_spline_distance.size();i++) 
    *c->verbose << n4_spline_distance[i]<<" ";
  *c->verbose << std::endl;
  
  // Set up a filter to shrink image by a factor
  typedef itk::ShrinkImageFilter<ImageType, ImageType> ShrinkerType;
  typename ShrinkerType::Pointer shrinker = ShrinkerType::New();
  shrinker->SetInput(mri);
  shrinker->SetShrinkFactors(n4_shrink_factor);
  shrinker->Update();

  typedef itk::OtsuThresholdImageFilter<ImageType, ImageType> ThresholderType;
  typename ThresholderType::Pointer otsu = ThresholderType::New();
  
  if(!n4_use_mask)
  {
    // Compute mask using Otsu threshold
    otsu->SetInput(mri);
    otsu->SetNumberOfHistogramBins( n4_histogram_bins );
    otsu->SetInsideValue( 0 );
    otsu->SetOutsideValue( 1 );
    otsu->Update();
    mask = otsu->GetOutput();
    *c->verbose << "   Otsu threshold: "<<otsu->GetThreshold()<<endl;
  }

  // Shrink the mask
  typename ShrinkerType::Pointer maskshrinker = ShrinkerType::New();
  maskshrinker->SetInput( mask);
  maskshrinker->SetShrinkFactors(n4_shrink_factor);
  maskshrinker->Update();

  // Bias filter
  typedef itk::N3MRIBiasFieldCorrectionImageFilter<ImageType, ImageType, ImageType> CorrecterType;
  typename CorrecterType::Pointer correcter = CorrecterType::New();
  
  ImagePointer _img=shrinker->GetOutput();
  ImagePointer _mask=maskshrinker->GetOutput();
  
  correcter->SetInput( shrinker->GetOutput() );
  correcter->SetMaskImage( maskshrinker->GetOutput() );

  // Progress meter
  // typedef CommandIterationUpdate<CorrecterType> CommandType;
  // typename CommandType::Pointer observer = CommandType::New();
  // correcter->AddObserver( itk::IterationEvent(), observer );
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
  logField->SetOrigin( bspliner->GetOutput()->GetOrigin() );
  logField->SetSpacing( bspliner->GetOutput()->GetSpacing() );
  logField->SetRegions( bspliner->GetOutput()->GetLargestPossibleRegion().GetSize() );
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

  if( n4_output_field )
  {
    c->m_ImageStack.push_back(expFilter->GetOutput());
  } else {
    typedef itk::DivideImageFilter<ImageType, ImageType, ImageType> DividerType;
    typename DividerType::Pointer divider = DividerType::New();
    divider->SetInput1( mri );
    divider->SetInput2( expFilter->GetOutput() );
    divider->Update();

    // Update
    c->m_ImageStack.push_back(divider->GetOutput());
  }
}

// Invocations
template class BiasFieldCorrection<double, 2>;
template class BiasFieldCorrection<double, 3>;
