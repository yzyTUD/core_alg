#pragma once

#include <iostream>
#include <ctime>
#include <string>
#include <fstream>
#include <stdio.h>

using namespace std;

/// switching on or off deps. on platforms (will be ava. during cmake process)
#define USE_LOCALTIME 1

#if USE_LOCALTIME
std::string zr_get_current_time() {
	time_t rawtime;
	struct tm* timeinfo;

	time(&rawtime);
	timeinfo = localtime(&rawtime);
	int Hour = timeinfo->tm_hour;
	int Min = timeinfo->tm_min;
	int Sec = timeinfo->tm_sec;

	return std::to_string(Hour) + "h" + std::to_string(Min) + "m" + std::to_string(Sec) + "s";
}
#endif

void zr_output_ppm_image(std::vector<int> img_data, int width, int height, string filename) {
	ofstream my_Image(filename);
	if (my_Image.is_open()) {
		my_Image << "P3\n" << width << " " << height << " 255\n";
		for (int i = 0; i < width; i++) {
			for (int j = 0; j < height; j++) {
				int val = img_data.at(i + width * j);
					//value(i, j);
				my_Image << val << ' ' << 0 << ' ' << 0 << "\n";
			}
		}
		my_Image.close();
	}
	else cout << "Could not open the file";
}

// file related 
const int BUFFER_SIZE = 1024;
char* buf = (char*)malloc(BUFFER_SIZE);
char* zr_str2chararray(string s) {
	int n = s.length();
	char* char_array = (char*)malloc(n + 1);
	strcpy(char_array, s.c_str());
	return char_array;
}
void zr_double_the_given_file(string filedir) {
	FILE* f = fopen(zr_str2chararray(filedir), "rb");
	FILE* new_f = fopen(zr_str2chararray(filedir + "_newdoubled.txt"), "wb+");
	int read = 0;
	while (read = fread(buf, 1, 1, f) != 0) {
		fwrite(buf, 1, read, new_f);
	}
	//rewind(f);
	while (read = fread(buf, 1, 1, f) != 0) {
		fwrite(buf, 1, read, new_f);
	}
	fclose(f);
	fclose(new_f);
}
void zr_bitwise_reverse(string filedir) {
	int BUFFER_SIZE = 1024;
	char *buffer, *inv_buffer;
	FILE *f1, *f2;
	int read;
	buffer = (char*)malloc(BUFFER_SIZE);
	inv_buffer = (char*)malloc(BUFFER_SIZE);
	long size;
	f1 = fopen(zr_str2chararray(filedir), "rb"); 
	f2 = fopen(zr_str2chararray(filedir + "_new.txt"), "wb+");

	//from the byte before last
	fseek(f1, 0, SEEK_END);
	size = ftell(f1);

	//printf("size is:%ld", SEEK_END - SEEK_SET);
	//printf("size: %ld\n",size);

	if (size > BUFFER_SIZE)//file is big enough 
	{
		fseek(f1, -1L * BUFFER_SIZE, SEEK_END);
		int aligned_itr = size / BUFFER_SIZE;

		//loop until first byte of f1 reached
		while (aligned_itr) {
			//returns number of bytes that have been read succ.
			read = fread(buffer, 1, BUFFER_SIZE, f1);

			//reverse inside buffer 
			int i;
			for (i = 0; i < BUFFER_SIZE; i++)
				inv_buffer[i] = buffer[BUFFER_SIZE - 1 - i];

			fwrite(inv_buffer, 1, read, f2);
			if (aligned_itr != 1)
				fseek(f1, -2L * BUFFER_SIZE, SEEK_CUR); //has to be a long type of number 
			aligned_itr--;
		}

		//non-aligned part
		//move back...
		fseek(f1, -1L * BUFFER_SIZE, SEEK_CUR);
		int unaligned_itr = size - size / BUFFER_SIZE * BUFFER_SIZE;

		//start bytewise reverse...
		fseek(f1, -1L, SEEK_CUR);
		//loop until first byte of f1 reached
		while (unaligned_itr) {
			//returns number of bytes that have been read succ.
			read = fread(buffer, 1, 1, f1);
			fwrite(buffer, 1, read, f2);
			fseek(f1, -2L, SEEK_CUR); //has to be a long type of number 
			unaligned_itr--;
		}
	}
	else {// bytewise
		fseek(f1, -1L, SEEK_END);

		//loop until first byte of f1 reached
		//size lower than buffersize, no need to assign itr num...
		while (size) {
			//returns number of bytes that have been read succ.
			read = fread(buffer, 1, 1, f1);
			fwrite(buffer, 1, read, f2);
			fseek(f1, -2L, SEEK_CUR); //has to be a long type of number 
			size--;
		}
	}

	//close file 
	fclose(f1); fclose(f2);
	free(buffer);
	std::cout << "succ reverse!" << std::endl;
}