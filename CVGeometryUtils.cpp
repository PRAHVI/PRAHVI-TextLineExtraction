#include <opencv2/opencv.hpp>
#include "CVGeometryUtils.h"
#include <math.h>

using namespace cv;

Point getRectCenter(Rect r) {
	return (r.br() + r.tl()) * 0.5; 
}


void getLinePointsThroughRegionCenterAtDegree(double deg, const Point &center, int radius, Point &p1, Point &p2) {
	int del_x, del_y;

	del_x = radius * cos(deg * M_PI / 180.0);	
	del_y = - radius * sin(deg * M_PI / 180.0);	

	p1.x = center.x + del_x;
	p1.y = center.y + del_y;

	p2.x = center.x - del_x;
	p2.y = center.y - del_y;

}


void getPerpendicularLinePoints(int pos, int radius, double deg, LineIterator &l, Point &p1, Point &p2) {	
	Point cur_point = l.pos();
	double percent_traveled = static_cast<double>(pos) / l.count;
	int x = radius * (2.0 * percent_traveled - 1);
	int y = sqrt(pow(radius, 2) - pow(x, 2)); 
	int del_x = y * cos( (90.0 + deg) * (M_PI / 180.0)); 
	int del_y = y * sin( (90.0 + deg) * (M_PI / 180.0));

	p1.x = cur_point.x + del_x;
	p1.y = cur_point.y + del_y;	

	p2.x = cur_point.x - del_x;
	p2.y = cur_point.y - del_y;
}	


bool isRectInLine(const Rect &r, const Point &start, const Point &end) {
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

	cond1 = (static_cast<float>(bottom_left.y) > slope * bottom_left.x + intercept) != 
			(static_cast<float>(top_right.y) > slope * top_right.x + intercept);	

	cond2 = (static_cast<float>(bottom_right.y) > slope * bottom_right.x + intercept) != 
			(static_cast<float>(top_left.y) > slope * top_left.x + intercept);	

	return (cond1 || cond2); 

}

bool isRectInCircle(const cv::Rect &r, const cv::Point &center, int radius) {
	Point rect_center = getRectCenter(r);
	return pow(rect_center.x - center.x, 2) + pow(rect_center.y - center.y, 2) < pow(radius, 2);
}

