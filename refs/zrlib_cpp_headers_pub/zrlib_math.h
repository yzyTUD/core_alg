#pragma once

#include<iostream>
#include<queue>
#include<stack>
#include<cstring>
#include<string>
#include<complex> // for using complex numbers
#include<cfloat>

#include "zrlib_util.h" // include util package in my lib.

/// switching on or off deps. on platforms (will be ava. during cmake process)
#define OPENCV_SUPPORT 1
#define OPENGL3_SUPPORT 0
#define OPENMESH_SUPPORT 0
#define OPENSSL_SUPPORT 0

// constants 
#if 0
GLfloat posi_array_triangle[] = {
	/*-1, -1, 1, 1,
	-1, 1, -1, 1,
	1, -1, -1, 1,*/
	-1, -1, 0, 1,
		1, -1, 0, 1,
		0,  1, 0, 1,

	/*-1, -1, 1, 1,
	1, -1, -1, 1,
	1 , 1 , 1 ,1,

	1, -1, -1, 1,
	-1, 1, -1, 1,
	1 , 1 , 1 ,1,

	-1, 1, -1, 1,
	-1, -1, 1, 1,
	1 , 1 , 1 ,1*/

	/*0, 1, -1, 1,
	-1, -1, 0, 1,
	1 , -1 , -1 ,1*/

};
#endif

// simple ones 
#if 1 
void zr_calcu_num_of_rabbits(int iter = 20) {
	long f1, f2;
	int i;
	f1 = 1, f2 = 2;
	for (i = 1; i <= iter; i++)
	{
		printf("%12ld %12ld", f1, f2);
		if (i % 2 == 0) std::cout << "\n";
		f1 = f1 + f2;
		f2 = f1 + f2;
	}
	std::cout << "\n";
}

int  helper_get_pixelvalue_mandelbrot(int x, int y, int width, int height) {
	using namespace std;
	complex<float> point((float)x / width - 1.5, (float)y / height - 0.5);
	// we divide by the image dimensions to get values smaller than 1
	// then apply a translation
	complex<float> z(0, 0);
	unsigned int nb_iter = 0;
	while (abs(z) < 2 && nb_iter <= 34) {
		z = z * z + point;
		nb_iter++;
	}
	if (nb_iter < 34) return 255;
	else return 0;
}

void zr_calcu_mandelbrot_output_as_ppm(float width = 600, float height = 600) {
	int number_of_pixels = width * height; 
	std::vector<int> mandelbrotImg;
	for (int i = 0; i < width; i++) {
		for (int j = 0; j < height; j++) {
			mandelbrotImg.push_back(helper_get_pixelvalue_mandelbrot(j, i, width, height));
		}
	}
	zr_output_ppm_image(mandelbrotImg, width, height, "mandelbrot.ppm");
}

void zr_calcu_series() {
	float n = 1;
	float s = 1;
	while (1 / n >= FLT_EPSILON * s) {
		n++;
		s += 1 / (n - 1);
	}
	printf("n: %f\n", n);
}
#endif

// linear_algebra
#if 1
int  zr_lu_decomposition_under_dev() {
	int n; double z; double** a; int i, j;
	printf("请输入矩阵A的大小:\n");
	scanf("%d", &n);
	a = (double**)malloc(sizeof(double*) * n);//
	for (i = 0; i < n; i++)a[i] = (double*)malloc(sizeof(double) * n);
	double** a_temp; a_temp = (double**)malloc(sizeof(double*) * n);//
	for (i = 0; i < n; i++)a_temp[i] = (double*)malloc(sizeof(double) * n);

	//LU
	int k; double** m; double max;
	m = (double**)malloc(sizeof(double*) * n);
	for (i = 0; i < n; i++) {
		m[i] = (double*)malloc(sizeof(double) * n);
	}
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++)
			m[i][j] = 0;
	}
	printf("please input matrix A\n");
	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++) {
			scanf("%lf ", &z);
			a_temp[i][j] = a[i][j] = z;
		}
	printf("A:\n");
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++)
			printf(" %lf ", a[i][j]);
		printf("\n");
	}
	double* b, * b_temp;
	double s;
	b = (double*)malloc(sizeof(double) * n);
	b_temp = (double*)malloc(sizeof(double) * n);
	printf("please give vector b\n");
	scanf("%lf", &s);
	for (i = 0; i < n; i++) {
		scanf("%lf", &s);
		b_temp[i] = b[i] = s;
	}
	printf("b is:\n");
	for (i = 0; i < n; i++)
		printf("%lf ", b[i]);
	printf("\n");

	//start LU
	for (k = 0; k < n - 1; k++) {
		if (a[k][k] == 0) {
			printf("zero zhuyuan！\n");
			return 0;
		}
		//
		for (i = k + 1; i < n; i++)
			m[i][k] = a[i][k] / a[k][k];

		for (j = k + 1; j < n; j++)
			for (i = k + 1; i < n; i++)
				a[i][j] = a[i][j] - m[i][k] * a[k][j];
	}

	printf("m is:\n");
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++)
			printf("%lf ", m[i][j]);
		printf("\n");
	}
	printf("L is:\n");
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++)
			printf("%lf ", m[i][j]);
		printf("\n");
	}
	printf("U is:\n");
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			if (i > j)
				a[i][j] = 0;
			printf("%lf ", a[i][j]);
		}
		printf("\n");
	}
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++)
			if (i == j)
				m[i][j] = 1;
	}

	printf("the solution of the equ. is:\n");
	//solve Ly=b;
	//下三角的前代法
	double* y;
	y = (double*)malloc(sizeof(double) * n);
	for (i = 0; i < n; i++)
		y[i] = 0;

	for (j = 0; j < n; j++) {
		if (m[j][j] == 0) {
			printf("err: sigular\n");
			return 0;
		}
		y[j] = b[j] / m[j][j];
		for (i = j + 1; i < n; i++) {
			b[i] = b[i] - m[i][j] * y[j];
		}
	}
	printf("y\n");
	for (i = 0; i < n; i++)
		printf("%lf ", y[i]);
	printf("\n");

	//solve Ux=y;
	//上三角回代法
	double* x;
	x = (double*)malloc(sizeof(double) * n);
	for (i = 0; i < n; i++)
		x[i] = 0;
	for (j = n - 1; j >= 0; j--) {
		if (a[j][j] == 0) {
			printf("矩阵奇异\n");
			return 0;
		}
		x[j] = y[j] / a[j][j];
		for (i = 0; i < j; i++)
			y[i] = y[i] - a[i][j] * x[j];
	}
	printf("输出x值\n");
	for (i = 0; i < n; i++)
		printf("%lf ", x[i]);
	printf("\n ");

	//计算残差
	printf("输出残差\n");
	double* r;
	r = (double*)malloc(sizeof(double) * n);
	for (i = 0; i < n; i++)
		r[i] = 0;
	for (i = 0; i < n; i++) {//i是行，遍历所有行
		for (j = 0; j < n; j++) {//一行中从左到右
			r[i] += a_temp[i][j] * x[j];
		}
	}

	double* can;
	can = (double*)malloc(sizeof(double) * n);
	for (i = 0; i < n; i++)
		can[i] = b_temp[i] - r[i];
	for (i = 0; i < n; i++)
		printf("%lf ", can[i]);

	free(r);
	free(can);
	free(x);
	free(y);
	free(m);
	free(b);
	free(a);
	free(b_temp);
	free(a_temp);
}
#endif

// convolution operations, image related 
#if 1
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
#endif

//
#if OPENSSL_SUPPORT
#include<openssl/md5.h>
void zr_calcu_md5() {
	unsigned char* data = "123";
	unsigned char md[16];
	int i;
	char tmp[3] = { '\0' }, buf[33] = { '\0' };
	MD5(data, strlen(data), md);
	for (i = 0; i < 16; i++) {
		sprintf(tmp, "%02x", md[i]);
		strcat(buf, tmp);
	}
	printf("%s\n", buf);
}
#endif

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
void zr_to_binary_image(Mat image, float illu_threshold) {
	//Mat temImage = image.clone();
	for (int i = 0; i < image.rows; i++)
	{
		for (int j = 0; j < image.cols; j++)
		{
			float sumtmp = 0;
			for (int k = 0; k < 3; k++)
			{
				sumtmp += image.at<Vec<uchar, 3>>(i, j)[k];
			}
			sumtmp = sumtmp / 3.0f;
			if (sumtmp > illu_threshold) {
				image.at<Vec<uchar, 3>>(i, j) = Vec<uchar, 3>(255, 255, 255);
			}
			else {
				image.at<Vec<uchar, 3>>(i, j) = Vec<uchar, 3>(0, 0, 0);
			}
		}
	}
}
#endif