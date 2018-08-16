#include <opencv2\opencv.hpp>
#include <opencv2\highgui\highgui.hpp>
#include <opencv2\xfeatures2d\nonfree.hpp>
#include <opencv2\features2d\features2d.hpp>
#include <opencv2\calib3d\calib3d.hpp>

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>

#include <Windows.h>

#include "FeatureExtractor.h"
#include "Timer.h"

using namespace std;
using namespace cv;


/* Constructor of object class */
FeatureExtractor::FeatureExtractor
	(Ptr<FeatureDetector>& detect, 
	Ptr<DescriptorExtractor>& describe, 
	Ptr<DescriptorMatcher>& match)
	: ratio_thresh(0.7f), refineF(true), epi_distance(2.0), confidence(0.99)
{
	detector = detect;
	descriptor = describe;
	matcher = match;

	detectionTime.clear();
	descriptionTime.clear();
	knnMatchingTime1.clear();
	knnMatchingTime2.clear();

	TPrefineTime.clear();
}

/* Setters of private class members */
void FeatureExtractor::setFeatureDetector(Ptr<FeatureDetector>& myDetector)
{
	detector = myDetector;
}

void FeatureExtractor::setFeatureDescriptor(Ptr<DescriptorExtractor>& myDescriptor)
{
	descriptor = myDescriptor;
}

void FeatureExtractor::setDescriptorMatcher(Ptr<DescriptorMatcher>& myMatcher)
{
	matcher = myMatcher;
}

void FeatureExtractor::setRatio(float r)
{
	ratio_thresh = r;
}

void FeatureExtractor::setRefineF(bool tf)
{
	refineF = tf;
}

void FeatureExtractor::setEpiDistance(double d)
{
	epi_distance = d;
}

void FeatureExtractor::setConfidence(double c)
{
	confidence = c;
}

void FeatureExtractor::setOutStream(String& tp, String& log)
{
	outTP = ofstream(tp);
	outLog = ofstream(log);

	if(!outTP)
	{
		cout<< "ERROR!! : Cannot open the Output TP File stream" << endl;
	}

	if(!outLog)
	{
		cout<< "ERROR!! : Cannot open the Output Log File stream" << endl;
	}

}

/* Additional functions */
int FeatureExtractor::checkNNratio
	(vector<vector<DMatch>> &matches)
{
	int removed_num = 0;

	float ratio_temp;
	vector<vector<DMatch>>::iterator matcheIterator;
	
	for(matcheIterator = matches.begin();
		matcheIterator != matches.end();
		++matcheIterator)
	{
		if(matcheIterator->size() > 1)
		{
			ratio_temp = (*matcheIterator)[0].distance / (*matcheIterator)[1].distance;
			
			if(ratio_temp > ratio_thresh)
			{
				matcheIterator->clear();
				removed_num++;
			}
		}

		else
		{
			matcheIterator->clear();
			removed_num++;
		}
	}

	return removed_num;
}

void FeatureExtractor::checkSymmetricity
(vector<vector<DMatch>>& matches1, 
vector<vector<DMatch>>& matches2, 
vector<DMatch>& symMatches)
{
	vector<vector<DMatch>>::iterator matcheIterator1;
	vector<vector<DMatch>>::iterator matcheIterator2;

	for(matcheIterator1 = matches1.begin(); 
		matcheIterator1 != matches1.end();
		++matcheIterator1)
	{
		if(matcheIterator1->size() < 2)
			continue;

		for(matcheIterator2 = matches2.begin();
			matcheIterator2 != matches2.end();
			++matcheIterator2)
		{
			if(matcheIterator2->size() < 2)
				continue;

			if((*matcheIterator1)[0].queryIdx == (*matcheIterator2)[0].trainIdx
				&& (*matcheIterator2)[0].queryIdx == (*matcheIterator1)[0].trainIdx)
			{
				symMatches.push_back(DMatch((*matcheIterator1)[0].queryIdx, (*matcheIterator1)[0].trainIdx, (*matcheIterator1)[0].distance));
				break;
			}
		}
	}
}

Mat FeatureExtractor::checkRANSAC
	(vector<DMatch>& matches,
	vector<KeyPoint>& keypoint1,
	vector<KeyPoint>& keypoint2,
	vector<DMatch>& validMatches,
	Mat& img1,
	Mat& img2,
	float& repeatability,
	int& corrCount)
{
	vector<Point2f> point1(matches.size());
	vector<Point2f> point2(matches.size());

	vector<DMatch>::const_iterator iter1;
	vector<DMatch>::const_iterator iter2;

	float x1, y1;
	float x2, y2;

	int idx1 = 0;
	
	for(iter1 = matches.begin();
		iter1 != matches.end();
		++iter1)
	{
		x1 = keypoint1[iter1->queryIdx].pt.x;
		y1 = keypoint1[iter1->queryIdx].pt.y;
		point1.at(idx1) = Point2f(x1, y1);

		x2 = keypoint2[iter1->trainIdx].pt.x;
		y2 = keypoint2[iter1->trainIdx].pt.y;
		point2.at(idx1) = Point2f(x2, y2);

		idx1++;
	}


	vector<uchar> inliers(point1.size(), 0);

	Mat fundamentalMat, homographMat;
	//fundamentalMat = findFundamentalMat(
	//	Mat(point1), Mat(point2),	// matching points
	//	inliers,		// match status (inlier or outlier)
	//	CV_FM_RANSAC,	// RANSAC method
	//	epi_distance,	// distance to epipolar line
	//	confidence);	// confidence probability
	
	homographMat = findHomography(point1, point2, inliers, CV_RANSAC, 5);

	// extract the surviving (inliers) matches
	vector<uchar>::const_iterator iterInliers = inliers.begin();
	vector<DMatch>::const_iterator iterMatches = matches.begin();

	// for all matches
	for( ; iterInliers != inliers.end(); ++iterInliers, ++iterMatches)
	{
		if(*iterInliers)	 // it is a valid match
		{
			validMatches.push_back(*iterMatches);
		}
	}

	if(refineF)
	{
		// The F matrix will be recomputed with
		// all accepted matches
		// Convert keypoints into Point2f
		// for final F computation

		vector<Point2f> point11(validMatches.size());
		vector<Point2f> point22(validMatches.size());

		int idx2 = 0;

		for(iter2 = validMatches.begin();
			iter2 != validMatches.end();
			++iter2)
		{
			x1 = keypoint1[iter2->queryIdx].pt.x;
			y1 = keypoint1[iter2->queryIdx].pt.y;
			point11.at(idx2) = Point2f(x1, y1);

			x2 = keypoint2[iter2->trainIdx].pt.x;
			y2 = keypoint2[iter2->trainIdx].pt.y;
			point22.at(idx2) = Point2f(x2, y2);

			idx2++;

		}

		/*fundamentalMat = findFundamentalMat(
			Mat(point11),
			Mat(point22),
			CV_FM_8POINT);*/
		if(point11.size() >= 4) {
			homographMat = findHomography(point11, point22, CV_RANSAC, 5);

			evaluateFeatureDetector(img1, img2, homographMat,
				&keypoint1, &keypoint2, repeatability, corrCount);

			cout << "" << endl;
		}
		else {
			repeatability = 0;
			corrCount = 0;

			cout << "" << endl;
		}
		

	}

	return homographMat;
}

void FeatureExtractor::setImgID(vector<DMatch>& matches, int imgID)
{
	for (int i = 0; i < matches.size(); i++)
	{
		matches[i].imgIdx = imgID;
	}
}

void FeatureExtractor::setGPID
(vector<DMatch>& prevMatches,
	vector<DMatch>& curMatches,
	vector<int>& prevGPID,
	vector<int>& curGPID,
	vector<KeyPoint>& curKeypointA,
	vector<KeyPoint>& curKeypointB,
	vector<int>& addGPID)
{
	int GPID;
	int maxGPID;
	int imgID = curMatches[0].imgIdx - 1;

	if (curMatches[0].imgIdx != 1)
	{
		vector<DMatch>::const_iterator prevMIter;
		vector<DMatch>::const_iterator curMIter;

		vector<int>::const_iterator prevGPIter, curGPIter;
		//vector<int> addGPID;

		prevGPIter = max_element(prevGPID.begin(), prevGPID.end());
		maxGPID = prevGPID.at(prevGPIter - prevGPID.begin());

		curGPID.clear();

		for (curMIter = curMatches.begin();
			curMIter != curMatches.end();
			++curMIter)
		{
			for (prevMIter = prevMatches.begin();
				prevMIter != prevMatches.end();
				++prevMIter)
			{
				if (prevMIter->queryIdx == curMIter->queryIdx)
					// Check previous and current matches by the acquried image
				{
					GPID = prevGPID.at(prevMIter - prevMatches.begin());
					break;
				}
			}

			if (prevMIter == prevMatches.end())
			{
				GPID = ++maxGPID;
			}
			// GP ID list of the current match by the acquried image
			curGPID.push_back(GPID);
		}


		bool gp_existence;

		for (int i = 0; i < curMatches.size(); i++)
		{
			gp_existence = false;

			for (int j = 0; j < prevGPID.size(); j++)
			{
				if (curGPID.at(i) == prevGPID.at(j))
				{
					// If the GP ID is already in GPID list,
					// we do not need to add it to the list
					gp_existence = true;
					break;
				}
			}

			if (!gp_existence)
			{
				// If the GP ID is not in GPID list,
				// It needs to add it to the list
				addGPID.push_back(curGPID.at(i));
				outTP <<
					0 << "\t" <<
					curGPID.at(i) << "\t" <<
					curKeypointA[curMatches.at(i).queryIdx].pt.x << "\t" <<
					curKeypointA[curMatches.at(i).queryIdx].pt.y << endl;
			}

		}
		
		//imgID++;

		for(int k = 0; k < curMatches.size(); k++)
		{
			outTP <<
				imgID << "\t" <<
				curGPID.at(k) << "\t" <<
				curKeypointB[curMatches.at(k).trainIdx].pt.x << "\t" <<
				curKeypointB[curMatches.at(k).trainIdx].pt.y << endl;
		}
	}
}

void FeatureExtractor::setGPID
		(vector<DMatch>& initMatches,
		vector<int>& initGPID,
		vector<KeyPoint>& curKeypointA,
		vector<KeyPoint>& curKeypointB)
{
	int gpid;
	int imgID = initMatches[0].imgIdx - 2;

	for(int i = 0; i < initMatches.size(); i++)
	{
		gpid = i + 1;
		initGPID.push_back(gpid);
	}


	for(int j = 0; j < initMatches.size(); j++)
	{
		/*outTP <<
			imgID << "\t" <<
			initGPID.at(j) << "\t" <<
			curKeypointA[initMatches.at(j).queryIdx].pt.x << "\t" <<
			curKeypointA[initMatches.at(j).queryIdx].pt.y << endl;*/
		outTP <<
			0 << "\t" <<
			initGPID.at(j) << "\t" <<
			curKeypointA[initMatches.at(j).queryIdx].pt.x << "\t" <<
			curKeypointA[initMatches.at(j).queryIdx].pt.y << endl;
	}

	imgID++;
	
	for(int j = 0; j < initMatches.size(); j++)
	{
		outTP <<
			imgID << "\t" <<
			initGPID.at(j) << "\t" <<
			curKeypointB[initMatches.at(j).trainIdx].pt.x << "\t" <<
			curKeypointB[initMatches.at(j).trainIdx].pt.y << endl;
	}
}

void FeatureExtractor::setCorrectMatchesMask
	(/*input*/ vector<vector<DMatch>>& matches1to2,
	vector<DMatch>& corrMatches,
	/*output*/ vector<vector<uchar>>& corrMatchesMask)
{
	int idxRow, idxCol;

	if(matches1to2.size() == corrMatchesMask.size())
	{
		for(size_t i = 0; i < corrMatches.size(); i++)
		{
			idxRow = corrMatches[i].queryIdx;

			for(size_t j = 0; j < matches1to2[idxRow].size(); j++)
			{
				if(matches1to2[idxRow][j].trainIdx == corrMatches[i].trainIdx)
				{
					corrMatchesMask[idxRow].push_back(1);
				}
				else
				{
					corrMatchesMask[idxRow].push_back(0);
				}
			}
			
		}
	}
}



/* Main function of object class */
void FeatureExtractor::matchImgSeq(/* input */ vector<String>& imgList,
		/* output */vector<vector<DMatch>>& matches,
		/*vector<vector<int>>& GPID,*/
		vector<int>& GPID,
		vector<vector<KeyPoint>>& keypoints)
{
	Mat imgA, imgB;
	int iter;

	vector<DMatch> curMatches;
	vector<int> curGPID;

	vector<DMatch> prevMatches;
	vector<int> prevGPID;

	vector<KeyPoint> keypointsA,keypointsB;
	Mat desc_vecA, desc_vecB;

	vector<vector<DMatch>> NNmatches1, NNmatches2;
	int NNremoved1, NNremoved2;

	vector<DMatch> symMatches;
	Mat homo_mat;

	Timer *myTimer = new Timer();

	vector<vector<DMatch>> allMatches1to2;
	vector<vector<uchar>> allCorectMatchesMask;
	vector<Point2f> recallPrecisionCurve;

	vector<int> addGPID;
	imgA = imread(imgList[0], CV_LOAD_IMAGE_GRAYSCALE);
	detector->detect(imgA, keypointsA);
	descriptor->compute(imgA, keypointsA, desc_vecA);
	keypoints.push_back(keypointsA);

	for (iter = 1; iter < imgList.size(); iter++)
	{

		cout << iter << endl;

		//outLog <<
		//	imgList[0] << "\t" << // Query Image
		//	imgList[iter] << "\t"; // Train Image

		outLog <<
			0 << "\t" << // Query Image
			iter << "\t"; // Train Image

		if (iter == 1)
		{
			//imgA = imread(imgList[iter], CV_LOAD_IMAGE_GRAYSCALE);
			imgB = imread(imgList[iter], CV_LOAD_IMAGE_GRAYSCALE);

			myTimer->start();
			//detector->detect(imgA, keypointsA);
			detector->detect(imgB, keypointsB);

			outLog <<
				myTimer->finish() << "\t"; // Detection Time

			myTimer->start();
			//descriptor->compute(imgA, keypointsA, desc_vecA);
			descriptor->compute(imgB, keypointsB, desc_vecB);

			outLog <<
				myTimer->finish() << "\t"; // Description Time


			//keypoints.push_back(keypointsA);
			keypoints.push_back(keypointsB);

			outLog <<
				keypointsA.size() << "\t" << // # of query keypoints
				keypointsB.size() << "\t"; // # of train keypoints

		}
		else
		{
			//imgA = *(&imgB);
			//keypointsA = *(&keypointsB);
			//desc_vecA = *(&desc_vecB);

			imgB.release();
			keypointsB.clear();
			desc_vecB.release();

			imgB = imread(imgList[iter], CV_LOAD_IMAGE_GRAYSCALE);

			myTimer->start();
			detector->detect(imgB, keypointsB);

			outLog <<
				myTimer->finish() << "\t"; // Detection Time

			myTimer->start();
			descriptor->compute(imgB, keypointsB, desc_vecB);

			outLog <<
				myTimer->finish() << "\t"; // Description Time

			keypoints.push_back(keypointsB);

			outLog <<
				keypointsA.size() << "\t" << // # of query keypoints
				keypointsB.size() << "\t"; // # of train keypoints

		}

		myTimer->start();
		matcher->knnMatch(desc_vecA, desc_vecB, NNmatches1, 2);

		outLog <<
			myTimer->finish() << "\t" << // Time for KNN1
			NNmatches1.size() << "\t"; // # of KNN1 matches

		myTimer->start();
		matcher->knnMatch(desc_vecB, desc_vecA, NNmatches2, 2);

		outLog <<
			myTimer->finish() << "\t" << // Time for KNN2
			NNmatches2.size() << "\t"; // # of KNN2 matches

		myTimer->start();
		NNremoved1 = checkNNratio(NNmatches1);

		outLog <<
			myTimer->finish() << "\t" << // Time for NN ratio test 1
			NNremoved1 << "\t"; // # of removed KNN1 matches


		myTimer->start();
		NNremoved2 = checkNNratio(NNmatches2);

		outLog <<
			myTimer->finish() << "\t" << // Time for NN ratio test 2
			NNremoved1 << "\t"; // # of removed KNN2 matches


		myTimer->start();
		checkSymmetricity(NNmatches1, NNmatches2, symMatches);

		outLog <<
			myTimer->finish() << "\t" << // Time for symmetricity check
			symMatches.size() << "\t"; // # of symmetric matches

		float repeatability;
		int corrCount;

		myTimer->start();
		if (symMatches.size() != 0) {
			homo_mat = checkRANSAC(symMatches, keypointsA, keypointsB, curMatches, imgA, imgB, repeatability, corrCount);
		}

		outLog <<
			myTimer->finish() << "\t" << // Time for RANSAC
			curMatches.size() << "\t" << // # of final matches
			repeatability << "\t" << // repeatability to evaluate detector
			corrCount << endl; // # of correpending pairs(point to point) for repeatability

		Mat matchingResult;
		drawMatches(imgA,
			keypointsA,
			imgB,
			keypointsB,
			curMatches,
			matchingResult,
			Scalar(255, 0, 0),
			Scalar(255, 255, 255),
			vector<char>(),
			DrawMatchesFlags::NOT_DRAW_SINGLE_POINTS);

		//imwrite("result_Query_"+imgList[iter], matchingResult);
		string output = "./ResultImage/result_Query_" + to_string(iter) + ".png";
		imwrite(output, matchingResult);



		/*vector<vector<DMatch>> matches1to2;
		vector<vector<uchar>> correctMatchesMask;
		vector<Point2f> rp_curve;

		cv::evaluateGenericDescriptorMatcher(imgA, imgB, homo_mat, keypointsA, keypointsB,
		&matches1to2, &correctMatchesMask, rp_curve, descMatch);

		allMatches1to2.insert(allMatches1to2.end(), matches1to2.begin(), matches1to2.end());
		allCorectMatchesMask.insert(allCorectMatchesMask.end(), correctMatchesMask.begin(), correctMatchesMask.end());*/


		setImgID(curMatches, iter + 1);

		matches.push_back(curMatches);

		// Current match exists
		if (curMatches.size() != 0) {
			// Not 1st match & 1st match exists
			if (iter != 1 && matches[0].size() != 0) {
				setGPID(matches[0], matches[iter - 1], GPID, curGPID, keypointsA, keypointsB, addGPID);
				for (int i = 0; i < addGPID.size(); i++) {
					GPID.push_back(addGPID.at(i));
				}
			}
			// Not 1st match & 1st match doesn't exist
			else if (iter != 1 && matches[0].size() == 0) {
				matches[0] = curMatches;
				setGPID(matches[0], curGPID, keypointsA, keypointsB);
				for (int i = 0; i < curGPID.size(); i++) {
					GPID.push_back(curGPID.at(i));
				}
			}
			// 1st match
			else {
				setGPID(matches[0], curGPID, keypointsA, keypointsB);
				for (int i = 0; i < curGPID.size(); i++) {
					GPID.push_back(curGPID.at(i));
				}
			}			
		}
		// Current match doesn't exist
		else {
			cout << "No matches!!!" << endl;
		}		


		//if (curMatches.size() == 0)
		//{
		//	cout << "No matches!!!" << endl;
		//	//exit(0);
		//}
		//if (iter == 1)
		//{
		//	if (curMatches.size() != 0) {
		//		setGPID(matches[0], curGPID, keypointsA, keypointsB);
		//		for (int i = 0; i < curGPID.size(); i++) {
		//			GPID.push_back(curGPID.at(i));
		//		}
		//	}
		//}
		//else
		//{
		//	// Current matches exist & 1st matches don't exist
		//	if (curMatches.size() != 0 && matches[0].size() == 0) {
		//		matches[0] = curMatches;
		//		//setGPID(matches[0], matches[iter-1], GPID[0], curGPID, keypointsA, keypointsB);
		//		//setGPID(matches[0], matches[iter - 1], GPID, curGPID, keypointsA, keypointsB, addGPID);
		//		setGPID(matches[0], curGPID, keypointsA, keypointsB);
		//		for (int i = 0; i < curGPID.size(); i++) {
		//			GPID.push_back(curGPID.at(i));
		//		}
		//	}
		//	// Current matches exist & 1st matches exist
		//	else if (curMatches.size() != 0 && matches[0].size() != 0) {
		//		setGPID(matches[0], matches[iter - 1], GPID, curGPID, keypointsA, keypointsB, addGPID);
		//		for (int i = 0; i < addGPID.size(); i++) {
		//			GPID.push_back(addGPID.at(i));
		//		}
		//	}
		//}

		//GPID.push_back(curGPID);

		
		

		NNmatches1.clear();
		NNmatches2.clear();
		symMatches.clear();
		curMatches.clear();
		curGPID.clear();

	}


//	OneWayDescriptorMatch *match = new OneWayDescriptorMatch();

//	Ptr<GenericDescriptorMatcher> descMatch = match;

	outTP.close();
	outLog.close();
}
	
	