#include <Windows.h>
#include <iostream>
#include <stdio.h>

#include "Timer.h"

using namespace std;

void Timer::start()
{
	QueryPerformanceFrequency((LARGE_INTEGER*) (&frequency));

	QueryPerformanceCounter((LARGE_INTEGER*) (&startF));
}

float Timer::finish()
{
	QueryPerformanceCounter((LARGE_INTEGER*) (&endF));

	executionT = (float) (endF - startF) / frequency * (float) 1000.0;

	return executionT;
}