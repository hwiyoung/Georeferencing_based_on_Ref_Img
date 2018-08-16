#include <opencv2\opencv.hpp>
#include <opencv2\highgui\highgui.hpp>
#include <opencv2\nonfree\features2d.hpp>
#include <opencv2\features2d\features2d.hpp>
#include <opencv2\calib3d\calib3d.hpp>

#include <iostream>
#include <fstream>
#include <vector>

#include <Windows.h>

using namespace cv;
using namespace std;

typedef struct IP *nodePointer;
typedef struct
{
	nodePointer leftLink;
	nodePointer rightLink;

	int imgID;
	int keypointID;
} IP;
