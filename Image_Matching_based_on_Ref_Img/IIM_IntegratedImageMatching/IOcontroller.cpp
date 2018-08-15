#include <opencv2\opencv.hpp>
//#include <opencv2\highgui\highgui.hpp>
//#include <opencv2\nonfree\features2d.hpp>
//#include <opencv2\features2d\features2d.hpp>
//#include <opencv2\calib3d\calib3d.hpp>

#include <iostream>
#include <fstream>            
#include <vector>

#include "IOcontroller.h"

using namespace std;
using namespace cv;


IOcontroller::IOcontroller(String& inF)
{
	inConfigFile = inF;

	inImg = ifstream(inConfigFile);

	if(!inImg)
	{
		cout<< "ERROR!! : Cannot open the Input Config File" << endl;
	}
}

void IOcontroller::readConfigFile(vector<String>& imgList, String path)
{
	String tmp;
	const int buffSize = 200;

	char buff[buffSize] = {0};

	while(!inImg.eof())
	{
		inImg.getline(buff, buffSize);
		tmp = path + String(buff);
		imgList.push_back(tmp);
	}
}

