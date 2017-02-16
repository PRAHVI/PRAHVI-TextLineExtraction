#ifndef CVGEOMETRYUTILS_H
#define CVGEOMETRYUTILS_H

extern cv::Point getRectCenter(cv::Rect);
extern void getLinePointsThroughRegionCenterAtDegree(double, cv::Point, int, cv::Point, cv::Point);
extern void getPerpendicularLinePoints(int, double, cv::Point, int, cv::Point, cv::Point);
extern bool isRectInLine(cv::Rect, cv::Point, cv::Point);

#endif 
