#include <Windows.h>
#include <iostream>
#include <stdio.h>

class Timer
{
private :
	INT64 startF, endF;
	INT64 frequency;

	float executionT;

public :

	void start();
	float finish();

};
