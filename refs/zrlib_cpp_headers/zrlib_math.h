#include<iostream>
#include<queue>
#include<stack>
#include<cstring>
#include<string>

/// switching on or off deps. on platforms (will be ava. during cmake process)
#define OPENCV_SUPPORT true
#define OPENGL3_SUPPORT false
#define OPENMESH_SUPPORT false

/// no spec. deps. needed
void zr_img_blur(unsigned char* data_ptr, int C, int H, int W) {
	for (int i = 0; i < W; i += 3) {
		for (int j = 0; j < 3 * H; j++) {
			// do some calcu. 
			float avg = data_ptr[i + W * j] + data_ptr[i + W * j + 1] + data_ptr[i + W * j + 2];
			avg = avg / 3;
			// set the val.
			// R,G,B
			data_ptr[i + W * j] = (char)avg;
			data_ptr[i + W * j + 1] = (char)avg;
			data_ptr[i + W * j + 2] = (char)avg;
		}
	}
}

void zr_img_sobel(unsigned char* data_ptr, int C, int H, int W) {
	C = 3; // currently, assume that the test image has 3 channels
	int stencil[3][3] =
	{
		{-1,0,1},
		{-2,0,2},
		{-1,0,1}
	};
	float val_r;
	float val_g;
	float val_b;
	for (int i = 0; i < W; i += C) {
		for (int j = 0; j < C * H; j++) {
			// do some calcu. 
			if (i >= 1 && j >= 1 && i <= W - C && j <= C * H - 1) {
				val_r = stencil[0][0] * data_ptr[i - 1 * C + W * (j - 1)]	+	stencil[0][2] * data_ptr[i + 1 * C + W * (j - 1)];
					   +stencil[1][0] * data_ptr[i - 1 * C + W * j] +		    stencil[1][2] * data_ptr[i + 1 * C + W * j];
					   +stencil[2][0] * data_ptr[i - 1 * C + W * (j - 1)] +     stencil[2][2] * data_ptr[i + 1 * C + W * (j + 1)];
				val_g = stencil[0][0] * data_ptr[i - 1 * C + W * (j - 1) + 1] +		stencil[0][2] * data_ptr[i + 1 * C + W * (j - 1) + 1];
						+stencil[1][0] * data_ptr[i - 1 * C + W * j + 1] +			stencil[1][2] * data_ptr[i + 1 * C + W * j + 1];
						+stencil[2][0] * data_ptr[i - 1 * C + W * (j - 1) + 1] +	stencil[2][2] * data_ptr[i + 1 * C + W * (j + 1) + 1];
				val_b = stencil[0][0] * data_ptr[i - 1 * C + W * (j - 1) + 2] +		stencil[0][2] * data_ptr[i + 1 * C + W * (j - 1) + 2];
						+stencil[1][0] * data_ptr[i - 1 * C + W * j + 2] +			stencil[1][2] * data_ptr[i + 1 * C + W * j + 2];
						+stencil[2][0] * data_ptr[i - 1 * C + W * (j - 1) + 2] +	stencil[2][2] * data_ptr[i + 1 * C + W * (j + 1) + 2];
				// set the val. // R,G,B
				if (val_r > 70)
					val_r = 70;
				if (val_g > 70)
					val_g = 70;
				if (val_b > 70)
					val_b = 70;
				data_ptr[i + W * j] =	  (char)val_r;
				data_ptr[i + W * j + 1] = (char)val_g;
				data_ptr[i + W * j + 2] = (char)val_b;
			}
		}
	}
}

/// code block that needs opencv support 
#if OPENCV_SUPPORT
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"
using namespace cv;
using namespace std;
Mat zr_img_sobel_opencv(Mat image) {
	char mask[3][3] = { {-1,-2,-1},
						{ 0, 0, 0},
						{ 1, 2, 1}};

	Mat temImage = image.clone();
	for (int i = 1; i < image.rows - 1; i++)
	{
		for (int j = 1; j < image.cols - 1; j++)
		{
			for (int k = 0; k < 3; k++)
			{
				int pixel1 = image.at<Vec3b>(i - 1, j - 1)[k] * mask[0][0];
				int pixel2 = image.at<Vec3b>(i,     j - 1)[k] * mask[1][0];
				int pixel3 = image.at<Vec3b>(i + 1, j - 1)[k] * mask[2][0];

				int pixel4 = image.at<Vec3b>(i - 1, j)[k] * mask[0][1];
				int pixel5 = image.at<Vec3b>(i,     j)[k] * mask[1][1];
				int pixel6 = image.at<Vec3b>(i + 1, j)[k] * mask[2][1];

				int pixel7 = image.at<Vec3b>(i - 1, j + 1)[k] * mask[0][2];
				int pixel8 = image.at<Vec3b>(i,     j + 1)[k] * mask[1][2];
				int pixel9 = image.at<Vec3b>(i + 1, j + 1)[k] * mask[2][2];

				int sum = pixel1 + pixel2 + pixel3 + pixel4 + pixel5 + pixel6 + pixel7 + pixel8 + pixel9;
				if (sum < 0)
					sum = 0;
				if (sum > 255)
					sum = 255;
				temImage.at<Vec3b>(i, j)[k] = sum;
			}
		}
	}
	return temImage;
}
#endif