import SimpleITK as sitk
from joblib import Parallel, delayed
import pandas as pd
import itertools
import numpy as np
import os
import pickle

import radiomics
from numba import uint64
from pandas.core.nanops import nanstd
from radiomics import generalinfo, getFeatureClasses, imageoperations
from radiomics.featureextractor import RadiomicsFeatureExtractor

from skimage.exposure import equalize_adapthist, rescale_intensity
from os import walk
from os.path import splitext
from os.path import join
import pathlib

import cv2
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#HELP FUNCTIONS
def bias_correction_sitk(image_sitk, otsu_threshold=False, shrink_factor=0):
    """Apply N4 Bias Correction."""
    if shrink_factor:
        # N4BiasFieldCorrectionImageFilter takes too long to run, shrink image
        mask_breast = sitk.OtsuThreshold(image_sitk, 0, 1)
        shrinked_image_sitk = sitk.Shrink(image_sitk, [shrink_factor] * image_sitk.GetDimension())
        shrinked_mask_breast = sitk.Shrink(mask_breast, [shrink_factor] * mask_breast.GetDimension())
        corrector = sitk.N4BiasFieldCorrectionImageFilter()
        tmp_image = corrector.Execute(shrinked_image_sitk, shrinked_mask_breast)
        log_bias_field = corrector.GetLogBiasFieldAsImage(image_sitk)
        corrected_image_sitk = image_sitk / sitk.Exp(log_bias_field)
    else:
        initial_img = image_sitk
        # Cast to float to enable bias correction to be used
        tmp_image = sitk.Cast(image_sitk, sitk.sitkFloat64)
        # Set zeroes to a small number to prevent division by zero
        tmp_image = sitk.GetArrayFromImage(tmp_image)
        tmp_image[tmp_image == 0] = np.finfo(float).eps
        tmp_image = sitk.GetImageFromArray(tmp_image)
        tmp_image.CopyInformation(initial_img)
        if otsu_threshold:
            maskImage = sitk.OtsuThreshold(tmp_image, 0, 1)
        # Apply image bias correction using N4 bias correction
        corrector = sitk.N4BiasFieldCorrectionImageFilter()
        if otsu_threshold:
            corrected_image_sitk = corrector.Execute(tmp_image, maskImage)
        else:
            corrected_image_sitk = corrector.Execute(tmp_image)

    return corrected_image_sitk
#-----------------------------------------------------------------------------------------------------------------------
def clip_image_sitk(image_sitk, percentiles=[1, 99]):
    #Clip intensity range of an image.
    image_array = sitk.GetArrayFromImage(image_sitk)
    image_array = image_array.ravel()
    # Drop all zeroes from array
    image_array = image_array[image_array != 0]
    lowerbound = np.percentile(image_array, percentiles[0])
    upperbound = np.percentile(image_array, percentiles[1])
    # Create clamping filter for clipping and set variables
    filter = sitk.ClampImageFilter()
    filter.SetLowerBound(float(lowerbound))
    filter.SetUpperBound(float(upperbound))

    # Execute
    clipped_image_sitk = filter.Execute(image_sitk)

    return clipped_image_sitk
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
def retrieve_data(ii, all_masks):

    file = pathlib.Path(ii)
    print(file)
    organized_data = []
    starts = [all_masks.index(l) for l in all_masks if l.startswith(str(file.parent))]
    # -------------------------------------------------------------------------------------------------------------------
    # -------------------------------------------------------------------------------------------------------------------
    # -------------------------------------------------------------------------------------------------------------------
    # -------------------------------------------------------------------------------------------------------------------
    #STEP1: load MRI
    image_ = sitk.ReadImage(ii, sitk.sitkFloat32)
    raw_image = image_
    # -------------------------------------------------------------------------------------------------------------------
    # -------------------------------------------------------------------------------------------------------------------
    # -------------------------------------------------------------------------------------------------------------------
    # -------------------------------------------------------------------------------------------------------------------
    # STEP2: N4 bias field correction
    # MR image
    image_ = bias_correction_sitk(image_, otsu_threshold=True, shrink_factor=0)
    # -------------------------------------------------------------------------------------------------------------------
    # -------------------------------------------------------------------------------------------------------------------
    # -------------------------------------------------------------------------------------------------------------------
    # -------------------------------------------------------------------------------------------------------------------
    # STEP3: Clipping
    # Intensity Clipping to Percentile Range
    image_ = clip_image_sitk(image_, percentiles=[0.1, 99.9])
    # -------------------------------------------------------------------------------------------------------------------
    # -------------------------------------------------------------------------------------------------------------------
    # -------------------------------------------------------------------------------------------------------------------
    # -------------------------------------------------------------------------------------------------------------------
    # STEP4: 0-1024 Normalization
    image_ = sitk.GetImageFromArray(rescale_intensity(sitk.GetArrayFromImage(image_), out_range=(0, 1024)))
    image_.CopyInformation(raw_image)
    # -------------------------------------------------------------------------------------------------------------------
    # -------------------------------------------------------------------------------------------------------------------
    # -------------------------------------------------------------------------------------------------------------------
    # -------------------------------------------------------------------------------------------------------------------
    # STEP5: Load Masks
    # -----load expert ROI

    # create unique mask SI>0 --> 1
    expert_mask = sitk.ReadImage(all_masks[np.array(starts).astype(int).item()])
    expert_mask_array = sitk.GetArrayFromImage(expert_mask)

    expert_mask_array_unique = np.copy(expert_mask_array)
    expert_mask_array_unique[expert_mask_array_unique > 0] = 1
    expert_mask_array_unique_image = sitk.GetImageFromArray(expert_mask_array_unique, isVector=False)
    # -------------------------------------------------------------------------------------------------------------------
    # -------------------------------------------------------------------------------------------------------------------
    # -------------------------------------------------------------------------------------------------------------------
    # -------------------------------------------------------------------------------------------------------------------
    # STEP6: get slices containing OVERALL ROI mask
    expert_mask_array_unique_image.CopyInformation(raw_image)

    bb, correctedMask = imageoperations.checkMask(image_, expert_mask_array_unique_image)
    inputImage_slices, maskImage_slices = imageoperations.cropToTumorMask(image_, expert_mask_array_unique_image, bb)

    tmp, maskImage_slices_raw = imageoperations.cropToTumorMask(image_, expert_mask, bb)

    if maskImage_slices.GetSize()[2] == 1:
        print('Image contains a single slice')
        organized_data = []
    else:
        # ----- split multi-mask into separate masks
        submasks = np.unique(sitk.GetArrayFromImage(maskImage_slices_raw))[np.unique(sitk.GetArrayFromImage(maskImage_slices_raw)) != 0]

        for iii in submasks:
            tmp = 1 * (sitk.GetArrayFromImage(maskImage_slices_raw) == iii)
            tmp = sitk.GetImageFromArray(tmp, isVector=False)

            tmp.CopyInformation(inputImage_slices)
            maskImage_slices_raw.CopyInformation(inputImage_slices)

            if tmp.GetSize()[2] > 1:
                organized_data.append(
                {'MRI': inputImage_slices, 'mask': tmp, 'image_case': ii, 'mask_case': all_masks[np.array(starts).astype(int).item()],
                 'mask_scenario': 'mask_' + str(iii), 'mask_entire': maskImage_slices_raw})

    return organized_data

#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
def return_radiomics(preprocess_list, input_bin_width):

    settings = {}
    settings['interpolator'] = None
    settings['level'] = 1
    settings['distances'] = [1]
    settings['binWidth'] = [0]
    settings['removeOutliers'] = None
    settings['normalize'] = False
    #settings['normalizeScale'] = 100
    #settings['voxelArrayShift'] = 300
    settings['resampledPixelSpacing'] = None

    extractor = radiomics.featureextractor.RadiomicsFeatureExtractor(**settings)
    extractor.settings.pop('binWidth')
    extractor.settings.update({'binWidth': input_bin_width})
    extractor.disableAllImageTypes()
    extractor.disableAllFeatures()
    extractor.enableImageTypeByName('Original')
    extractor.enableFeatureClassByName('shape')
    extractor.enableFeatureClassByName('firstorder')
    extractor.enableFeatureClassByName('glcm')
    extractor.enableFeatureClassByName('glrlm')
    extractor.enableFeatureClassByName('glszm')
    extractor.enableFeatureClassByName('gldm')
    extractor.enableFeatureClassByName('ngtdm')

    expert_mask_array = sitk.GetArrayFromImage(preprocess_list['mask_entire'])
    expert_mask_array[expert_mask_array > 0] = 1
    expert_mask_array = sitk.GetImageFromArray(expert_mask_array, isVector=False)

    #-----------------------------------------------------------------------------------
    expert_mask_array.SetOrigin(preprocess_list['mask_entire'].GetOrigin())
    expert_mask_array.SetSpacing(preprocess_list['mask_entire'].GetSpacing())
    expert_mask_array.SetDirection(preprocess_list['mask_entire'].GetDirection())
    # -----------------------------------------------------------------------------------

    d1_expert = extractor.execute( preprocess_list['MRI'] , preprocess_list['mask'], 1)
    d1_roi_all = extractor.execute(preprocess_list['MRI'], expert_mask_array, 1)

    # delete keys containing features of Minimum or Maximum
    for k in list(d1_expert.keys()):
        if 'firstorder_Maximum' in k:
            del d1_expert[k]

    for k in list(d1_roi_all.keys()):
        if 'firstorder_Maximum' in k:
            del d1_roi_all[k]

    for k in list(d1_expert.keys()):
        if 'firstorder_Minimum' in k:
            del d1_expert[k]

    for k in list(d1_roi_all.keys()):
        if 'firstorder_Minimum' in k:
            del d1_roi_all[k]

    for k in list(d1_expert.keys()):
        if 'diagnostics_' in k:
            del d1_expert[k]

    for k in list(d1_roi_all.keys()):
        if 'diagnostics_' in k:
            del d1_roi_all[k]
    ###########################################################################
    tmp = preprocess_list['image_case'].split('/')
    tmp = tmp[len(tmp)-3:len(tmp)]

    d1_expert_df = pd.DataFrame(d1_expert.items(), columns = ['rads' , 'value'])#, index = d1_expert.keys())
    d1_expert_df['patient_id'] = tmp[0]
    d1_expert_df['rads'] =  preprocess_list['mask_scenario'] + '__' +tmp[2].split('.')[0] + '__' + tmp[1].replace(' ','_') + '__' + d1_expert_df['rads'].values

    d1_roi_all_df = pd.DataFrame(d1_roi_all.items(), columns=['rads', 'value'])  # , index = d1_expert.keys())
    d1_roi_all_df['patient_id'] = tmp[0]
    d1_roi_all_df['rads'] = 'mask_all' + '__' + tmp[2].split('.')[0] + '__' + tmp[1].replace(' ','_') + '__' + d1_roi_all_df['rads'].values

    df_overall = pd.concat([d1_expert_df, d1_roi_all_df], ignore_index=True)
    df_overall['rads'] = df_overall['rads'].replace(replacers, regex=True)

    return df_overall
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------
#main_path = ....
#mri_folder = main_path + '/prospective/' OR '/retrospective/'
#manual_annotation_name = 'tumour_mask'


def my_main_funct(mri_folder):
    replacers = {'mask_1__': 'mask_necrosis__', 'mask_2__': 'mask_tumor__', 'mask_3__': 'mask_oedema__'}
    #-----------------------------------------------------------------------------------
    # STEP #1 Import Data
    barlist = list()
    for root, dirs, files in walk(mri_folder):
        for f in files:
            if splitext(f)[1].lower() == ".nii" or splitext(f)[1].lower() == ".gz":
                barlist.append(join(root, f))

    all_masks = [x.replace('\\', '/') for x in barlist if "mask" in x]    #--->for windows

    images_only = [x.replace('\\', '/') for x in barlist if "mask" not in x]   #---> for windows
    all_masks_list = []
    for i in range(len(images_only)):
        all_masks_list.append(all_masks)
    #-----------------------------------------------------------------------------------

    data__list = Parallel(n_jobs=-1)(delayed(retrieve_data)(i,j) for i,j in zip(images_only, all_masks_list))

    flat_list = []
    for xs in data__list:
        for x in xs:
            flat_list.append(x)

    #-----------------------------------------------------------------------------------
    input_bin_width = 5
    bin_width_list = []
    bin_width_list.extend(input_bin_width for i in range(len(flat_list)))

    rad_features = Parallel(n_jobs=-1)(delayed(return_radiomics)(i,j) for i,j in zip(flat_list,bin_width_list))
    rad_features = pd.concat(rad_features)
    rad_features_complete = rad_features.drop_duplicates(subset=['rads', 'patient_id'], keep='first')
    rad_features_complete = rad_features_complete.rename(columns={'rads': 'Radiomics Feature', 'value': 'Radiomics Value'})

    rad_features_complete = [v for k, v in rad_features_complete.groupby('patient_id')]

    return rad_features_complete