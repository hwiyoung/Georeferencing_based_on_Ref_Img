#include <opencv2\opencv.hpp>
//#include <opencv2\highgui\highgui.hpp>
//#include <opencv2\nonfree\features2d.hpp>
//#include <opencv2\features2d\features2d.hpp>
//#include <opencv2\calib3d\calib3d.hpp>

#include <iostream>
#include <fstream>            
#include <vector>

#include <Windows.h>

using namespace std;
using namespace cv;

class IOcontroller
{
private:
	ifstream inImg;
	String inConfigFile;

public:
	IOcontroller(String& inF);

	void readConfigFile(vector<String>& imgList, String path);

	void writeTPFile(vector<vector<DMatch>> matches,
		vector<vector<int>> gpid,
		vector<vector<KeyPoint>> keypoints);

};