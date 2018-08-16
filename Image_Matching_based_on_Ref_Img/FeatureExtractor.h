#include <opencv2\opencv.hpp>
#include <opencv2\highgui\highgui.hpp>
#include <opencv2\xfeatures2d\nonfree.hpp>
#include <opencv2\features2d\features2d.hpp>
#include <opencv2\calib3d\calib3d.hpp>

#include <iostream>
#include <fstream>            
#include <vector>

#include <Windows.h>

using namespace std;
using namespace cv;


class FeatureExtractor
{
private :

	/* Pointers to the feature point detector, descriptor, matcher */
	// feature detector 
	Ptr<FeatureDetector> detector;
	// feature descriptor
	Ptr<DescriptorExtractor> descriptor;
	// descriptor matcher
	Ptr<DescriptorMatcher> matcher;

	/* Parameters to the keypoint filtering setup */
	// NN ratio test
	float ratio_thresh;
	// RANSAC - refine F matrix
	bool refineF;
	// RANSAC - minimum distance to epipolar
	double epi_distance;
	// RANSAC - confidence level(probability)
	double confidence;

	/* Execution Time in milliseconds */
	vector<float> detectionTime;
	vector<float> descriptionTime;
	vector<float> knnMatchingTime1, knnMatchingTime2;
	vector<float> TPrefineTime;

	/* Print Matching Result*/
	ofstream outTP;
	ofstream outLog;


public :
	
	/* Constructor of object class */
	/*
		 1) supported detectors in OpenCV
		 "FAST"
		 "STAR"
		 "SIFT"
		 "SURF"
		 "ORB"
		 "BRISK"
		 "MSER"
		 "GFTT"
		 "HARRIS"
		 "Dense"
		 "SimpleBlob"

		 2) supported descriptors in OpenCV
		 "SIFT"
		 "SURF"
		 "BRIEF"
		 "BRISK"
		 "ORB"
		 "FREAK"

		 3) supported matcher in OpenCV
		 "BruteForce"
		 "BruteForce-L1"
		 "BruteForce-Hamming"
		 "BruteForce-Hamming(2)"
		 "FlannBased"
	*/
	FeatureExtractor(
		Ptr<FeatureDetector>& detect, 
		Ptr<DescriptorExtractor>& describe, 
		Ptr<DescriptorMatcher>& match);


	/* Setters of private class members */
	void setFeatureDetector(Ptr<FeatureDetector>& myDetector);

	void setFeatureDescriptor(Ptr<DescriptorExtractor>& myDescriptor);

	void setDescriptorMatcher(Ptr<DescriptorMatcher>& myMatcher);

	void setRatio(float r);

	void setRefineF(bool tf);

	void setEpiDistance(double d);

	void setConfidence(double c);

	void setOutStream(String& tp, String& log);

	
	/* Additional functions */
	int checkNNratio(vector<vector<DMatch>> &matches);

	void checkSymmetricity
	(vector<vector<DMatch>>& matches1, 
	vector<vector<DMatch>>& matches2, 
	vector<DMatch>& symMatches);

	Mat checkRANSAC
		(vector<DMatch>& matches,
		vector<KeyPoint>& keypoint1,
		vector<KeyPoint>& keypoint2,
		vector<DMatch>& validMatches,
		Mat& img1,
		Mat& img2,
		float& repeatability,
		int& corrCount);

	// Set train(2nd) image ID
	void setImgID(vector<DMatch>& matches, int imgID);

	void setGPID
		(vector<DMatch>& prevMatches,
		vector<DMatch>& curMatches,
		vector<int>& prevGPID,
		vector<int>& curGPID,
		vector<KeyPoint>& curKeypointA,
		vector<KeyPoint>& curKeypointB,
		vector<int>& addGPID);

	void setGPID
		(vector<DMatch>& initMatches,
		vector<int>& initGPID,
		vector<KeyPoint>& curKeypointA,
		vector<KeyPoint>& curKeypointB);

	void setCorrectMatchesMask
		(/*input*/ vector<vector<DMatch>>& matches1to2,
		vector<DMatch>& corrMatches,
		/*output*/ vector<vector<uchar>>& corrMatchesMask);

	/* Main function of object class */
	void matchImgSeq(/* input */ vector<String>& imgList,
		/* output */vector<vector<DMatch>>& matches,
		/*vector<vector<int>>& GPID,*/
		vector<int>& GPID,
		vector<vector<KeyPoint>>& keypoints);

};