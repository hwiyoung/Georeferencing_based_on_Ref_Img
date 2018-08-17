# Georeferencing_based_on_Ref_Img
Georeferencing_based_on_Ref_Img is a SW to determine pose of images based on refernce images through 1) Image Matching + 2) Bundle Adjustment. This repository provides codes of each module and distributed program.

## Elements of *Georeferencing_based_on_Ref_Img*
* Georeferencing_based_on_Ref_Img_exe(distributed ver.)
* Image_Matching_based_on_Ref_Img
* BA_based_on_Ref_Img_Python
***

## Georeferencing_based_on_Ref_Img_exe
= (1) Image Matching(IIM_IntegratedImageMatching.exe) + (2) Bundle Adjustment(1 - Preprocess_IP_List.py, 2 - BBA_with_Reference_Images_v2.py)
* Input (1)
  * ./Input/(images)
  * ./Input/Img_List.txt - First line must have the name of the image want to determine
* Output (1)
  * ./data_Pre/Log.txt
  * ./data_Pre/TP.txt
* Input (2-1)
  * ./data_Pre/TP.txt
  * ./data_Pre/ImgList.txt
  * ./data_Pre/indexRef.txt
  * ./data_Pre/EO_ref_true.txt
  * ./data_Pre/EO_smart_true.txt
* Output (2-1)
  * ./data_BBA/IP_m.txt
  * ./data_BBA/EO_c_ref.txt
  * ./data_BBA/EO_i_tar.txt
* Input (2-2)
  * ./data_BBA/IO_ref.txt
  * ./data_BBA/IO_tar.txt
  * ./data_BBA/EO_c_ref.txt
  * ./data_BBA/EO_i_tar.txt
  * ./data_BBA/IP_m.txt
* Output (2-2)
  * ./data_BBA/GP_i.txt
  * ./data_BBA/IA_summary.txt
  * ./data_BBA/IP_mf.txt
  * ./data_BBA/EO_cf.txt
  * ./data_BBA/EO_e.txt
  * ./data_BBA/GP_e.txt
* In [BA_based_on_Ref_Img](https://github.com/hwiyoung/BA_based_on_Ref_Img) , you can check information about input & output params
* Run <run.py>
  * 1 - IIM_IntegratedImageMatching.exe
  * 2 - Preprocess_IP_List.py
  * 3 - BBA_with_Reference_Images_v2.py
***

## Image_Matching_based_on_Ref_Img
* Prerequisite
  * OpenCV 3.x
* Input - Download files from [HERE(input)](http://bit.ly/2waSloD) !!!
  * ./Input/(images)
  * ./Input/Img_List.txt - First line must have the name of the image want to determine
* Output
  * ./data_Pre/Log.txt
  * ./data_Pre/TP.txt
  * ./ResultImage/(result_images) - Examples in [HERE(result images)](http://bit.ly/2BcrVJ6)
* Things to do before to execute
  * Paste dll files from [HERE(dll)](http://bit.ly/2Ms0HCr) to the path same with .sln if you just want to execute the sln
  * Settings in 'Visual Studio' - **Release/x64**
* Run <IIM_IntegratedImageMatching.exe>
***

## BA_based_on_Ref_Img_Python
= Bundle Adjustment(1 - Preprocess_IP_List.py, 2 - BBA_with_Reference_Images_v2.py)
* Prerequisite
  * numpy
* Input (1)
  * ./data_Pre/TP.txt
  * ./data_Pre/ImgList.txt
  * ./data_Pre/indexRef.txt
  * ./data_Pre/EO_ref_true.txt
  * ./data_Pre/EO_smart_true.txt
* Output (1)
  * ./data_BBA/IP_m.txt
  * ./data_BBA/EO_c_ref.txt
  * ./data_BBA/EO_i_tar.txt
* Input (2)
  * ./data_BBA/IO_ref.txt
  * ./data_BBA/IO_tar.txt
  * ./data_BBA/IP_m.txt
  * ./data_BBA/EO_c_ref.txt
  * ./data_BBA/EO_i_tar.txt  
* Output (2)
  * ./data_BBA/GP_i.txt
  * ./data_BBA/IA_summary.txt
  * ./data_BBA/IP_mf.txt
  * ./data_BBA/EO_cf.txt
  * ./data_BBA/EO_e.txt
  * ./data_BBA/GP_e.txt
* In [BA_based_on_Ref_Img](https://github.com/hwiyoung/BA_based_on_Ref_Img) , you can check information about input & output params
* Run <Prerprocess_IP_List.py>
* Run <BBA_with_Reference_Images_v2.py>