# Georeferencing_based_on_Ref_Img
Georeferencing_based_on_Ref_Img is a SW to determine pose of images based on refernce images through 1) Image Matching + 2) Bundle Adjustment. This repository provides codes of each module and distributed program.

## Elements of *Georeferencing_based_on_Ref_Img*
* Georeferencing_based_on_Ref_Img_exe(distributed)
* Image_Matching_based_on_Ref_Img
* BA_based_on_Ref_Img_Python
***

## Georeferencing_based_on_Ref_Img_exe
= 1) Image Matching + 2) Bundle Adjustment
* Prepare *INPUT*
  * Images
  * Img_List.txt
    * First line must have the name of the image want to determine
  * Params(input) for [BA ref](https://github.com/hwiyoung/BA_based_on_Ref_Img)
* Run <run.py>
  * 1 - IIM_IntegratedImageMatching.exe
  * 2 - Preprocess_IP_List.py
  * 3 - BBA_with_Reference_Images_v2.py
***

## Image_Matching_based_on_Ref_Img
* Prepare *INPUT*
 * Images
 * Img_List.txt
   * First line must have the name of the image want to determine
* Run <IIM_IntegratedImageMatching.exe>
***

## BA_based_on_Ref_Img_Python
* Prepare *INPUT*
  * Params(input) for [BA ref](https://github.com/hwiyoung/BA_based_on_Ref_Img)
* Run <Prerprocess_IP_List.py>
* Run <BBA_with_Reference_Images_v2.py>