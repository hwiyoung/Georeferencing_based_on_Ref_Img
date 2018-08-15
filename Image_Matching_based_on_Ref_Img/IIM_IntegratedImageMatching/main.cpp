#include <opencv2\opencv.hpp>
#include <opencv2\highgui\highgui.hpp>
#include <opencv2\xfeatures2d\nonfree.hpp>
#include <opencv2\features2d\features2d.hpp>
#include <opencv2\calib3d\calib3d.hpp>

#include <iostream>
#include <fstream>
#include <vector>

#include <Windows.h>

#include "FeatureExtractor.h"
#include "IOcontroller.h"
#include "Timer.h"

using namespace cv;
using namespace std;
using namespace xfeatures2d;

void main()
{
	/*
		Set-up variables for Image Matching
	*/

	// input & output text file
	String path_img, path_config, path_result;
	String inImg, outTPF, outLogF;

	path_img = "./Input/";
	path_config = "./Input/";
	path_result = "./data_Pre/";

	inImg = path_config + "ImgList.txt";
	outTPF = path_result + "TP.txt";
	outLogF = path_result + "Log.txt";
	
	// file Input-Output controller
	IOcontroller IOhelper(inImg);
	
	// Img fie list to perform image matching
	vector<String> imgList;
	IOhelper.readConfigFile(imgList, path_img);

	// Feature detector, descriptor, matcher
	Ptr<FeatureDetector> detector;
	Ptr<DescriptorExtractor> descriptor;
	Ptr<DescriptorMatcher> matcher;
	
	detector = SIFT::create();
	descriptor = SIFT::create();
	matcher = DescriptorMatcher::create("FlannBased");

	FeatureExtractor mySTAR_ORB(detector, descriptor, matcher);
	mySTAR_ORB.setOutStream(outTPF, outLogF);

	vector<vector<DMatch>> matches;
	vector<int> gpid;
	vector<vector<KeyPoint>> keypoints;
	mySTAR_ORB.matchImgSeq(imgList, matches, gpid, keypoints);

	cout<< "The End" <<endl;
}
