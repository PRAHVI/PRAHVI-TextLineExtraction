#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <opencv2/opencv.hpp>
#include <random>

#define PI 3.14

using namespace cv;
using namespace std;

const int DIAMETER_VALUES[] = {64, 128, 256};

void mserExtractor(Mat&, vector<cv::Rect>&);
void getLinePointsAtDeg(int, cv::Point, int, cv::Point&, cv::Point&); 
cv::Point getRectCenterPoint(cv::Rect); 
void getLinePointsAtMidPoint(int, int, int, int, cv::Point&, cv::Point&, cv::Point&); 
void projectionProfile(cv::Rect, vector<cv::Rect>, int, int, Mat&, Mat&); 
bool rectInLine(cv::Rect, cv::Point, cv::Point); 
bool rectInCircle(cv::Rect, cv::Point, int); 




int main(int argc, char** argv) {
	Mat in_grey_img = imread(argv[1], IMREAD_GRAYSCALE);
	vector<cv::Rect> mser_bbox;	

	// Original Photo Display
	mserExtractor(in_grey_img, mser_bbox);
	
	//namedWindow("Projection Profile");
	namedWindow("Image");
	imshow("Image", in_grey_img);
	waitKey(0);

	// Main Loop
	srand(time(NULL));
	while(1) {	
		Mat img_region, projection_img;
		in_grey_img.copyTo(img_region);

		int rand_rect_ind = rand() % mser_bbox.size();
		int rand_deg = rand() % 360;
		int rand_diameter = DIAMETER_VALUES[rand() % 3];
		cv::Rect rand_cc = mser_bbox[rand_rect_ind];


		cout << "Degree: " << rand_deg << endl;
		cout << "Diameter: " << rand_diameter << endl;

		projectionProfile(rand_cc, mser_bbox, rand_diameter, rand_deg, img_region, projection_img);

		cout << endl;
		imshow("Image", img_region);
		imshow("Projection Profile", projection_img);

		int key_code = waitKey(0);	
		//cout << "The KeyCode pressed was " << key_code << endl;
		
		if (key_code == 1048689) {
			break;
		}
	}
		

}

void projectionProfile(cv::Rect cc, vector<cv::Rect> mser_bbox, int diameter, int degree, Mat& img_region, Mat& projection_img) {
	cv::Point center = getRectCenterPoint(cc);
	cv::Point p1, p2;
	int radius = diameter / 2;
	vector<cv::Rect> ccs_in_region;
		

	for (cv::Rect r : mser_bbox) {
		if (rectInCircle(r, center, diameter)) {
			ccs_in_region.push_back(r);
		}
	}		

	rectangle(img_region, cc, Scalar(0), -1);
	circle(img_region, center, radius, Scalar(0), 5);
	getLinePointsAtDeg(degree, center, radius, p1, p2);
	line(img_region, p1, p2, Scalar(0), 5);
	
	LineIterator proj_axis(img_region, p1, p2);
	vector<int> histogram(proj_axis.count, 0);

	cout << "Iterator length: " << proj_axis.count << endl;
	for (int i = 0; i < proj_axis.count; i++, proj_axis++) {
		cv::Point poi = proj_axis.pos();
		cv::Point start_point, end_point;
		getLinePointsAtMidPoint(i, proj_axis.count, radius, degree, poi, start_point, end_point);

		for (cv::Rect r : ccs_in_region) {
			if (rectInLine(r, start_point, end_point)) {
				histogram[i] += 1;
			}
		}

		// Testing Purposes
		//if (i % 40 == 0) {
		//	line(img_region, start_point, end_point, Scalar(0), 5);
		//	circle(img_region, poi, 5, Scalar(255), -1);
		//	imshow("Image", img_region);
		//	waitKey(0);
		//}	
	}
	
	projection_img = Mat(diameter*2, diameter, CV_8UC1, Scalar(0));	
	float proj_scale = static_cast<float>(histogram.size()) / diameter;
	for(int i = 0; i < diameter; i++) {
		int hist_ind = i * proj_scale;
		cv::Point bottomPoint = cv::Point(i, diameter*2);
		cv::Point topPoint = cv::Point(i, (diameter*2) - histogram[hist_ind]);
		
		line(projection_img, bottomPoint, topPoint, Scalar(255));	
	}
	//img_region.at<uchar>(p1) = 255;
}


bool rectInLine(cv::Rect r, cv::Point start, cv::Point end) {
	bool cond1, cond2;
	float slope = static_cast<float>(end.y - start.y) / (end.x - start.x);
	float intercept = static_cast<float>(end.y) - slope * end.x;

	cv::Point top_left, top_right, bottom_left, bottom_right;

	top_left = r.tl();
	top_right = r.tl();
	top_right.x += r.width;
	bottom_right = r.br();
	bottom_left = r.br();
	bottom_left.x -= r.width;

	cond1 = (static_cast<float>(bottom_left.y) > slope * bottom_left.x + intercept) != \
			(static_cast<float>(top_right.y) > slope * top_right.x + intercept);	
	cond2 = (static_cast<float>(bottom_right.y) > slope * bottom_right.x + intercept) != \
			(static_cast<float>(top_left.y) > slope * top_left.x + intercept);	

	return (cond1 || cond2); 
}

bool rectInCircle(cv::Rect r, cv::Point center, int radius) {
	cv::Point r_center = getRectCenterPoint(r); 
	return pow(r_center.x - center.x, 2) + pow(r_center.y - center.y, 2) < pow(radius, 2); 
}

void getLinePointsAtMidPoint(int ind_pos, int line_it_count, int radius, int deg, cv::Point& poi, cv::Point& start, cv::Point& end) {
	int x = radius * 2 * (static_cast<float>(ind_pos) / line_it_count) - radius; 	
	int y = sqrt(pow(radius, 2) - pow(x, 2));
	int del_x = y * cos(static_cast<float>((90 + deg) % 360) * PI / 180.0);
	int del_y = -y * sin(static_cast<float>((90 + deg) % 360) * PI / 180.0);
	
	//if (ind_pos % 40 == 0) {
	//	cout << "\tdelta x = " << del_x << endl;
	//	cout << "\tdelta y = " << del_y << endl;
	//	cout << "\tx = " << x << endl;
	//	cout << "\ty = " << y << endl;
	//	cout << endl;
	//}

	start.x = poi.x + del_x;
	start.y = poi.y + del_y;

	end.x = poi.x - del_x;
	end.y = poi.y - del_y;
}

void getLinePointsAtDeg(int degree, cv::Point center, int radius, cv::Point& p1, cv::Point& p2) {
	int del_x, del_y;

	del_x = radius * cos(static_cast<float> (degree) * PI / 180.0);	
	del_y = - radius * sin(static_cast<float> (degree) * PI / 180.0);	

	//cout << "del_x = " << del_x / static_cast<float>(N) << endl;
	//cout << "del_y = " << del_y / static_cast<float>(N) << endl;

	
	p1.x = center.x + del_x;
	p1.y = center.y + del_y;

	p2.x = center.x - del_x;
	p2.y = center.y - del_y;
}

cv::Point getRectCenterPoint(cv::Rect r) {
	return (r.br() + r.tl())*0.5; 
}

void mserExtractor(Mat& image, vector<cv::Rect>& mser_bbox) {
	Ptr<MSER> mserExtractor = MSER::create();
	vector<vector<cv::Point> > mserContours;

	mserExtractor->detectRegions(image, mserContours, mser_bbox);

	// Draw Contours and Bounding Boxes
	for (cv::Rect r : mser_bbox) {
		rectangle(image, r, Scalar(255));
	}

	for (vector<cv::Point> v : mserContours) {
		for (cv::Point p : v) {
			image.at<uchar>(p.y, p.x) = 255;
		}
	}

}
