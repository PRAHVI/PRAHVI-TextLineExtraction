#ifndef CVGEOMETRYUTILS_H
#define CVGEOMETRYUTILS_H

extern cv::Point getRectCenter(cv::Rect);
extern void getLinePointsThroughRegionCenterAtDegree(double, const cv::Point&, int, cv::Point&, cv::Point&);
extern void getPerpendicularLinePoints(int pos, int radius, double deg, cv::LineIterator &l, cv::Point&, cv::Point&);
extern bool isRectInLine(const cv::Rect &, const cv::Point&, const cv::Point&);
extern bool isRectInCircle(const cv::Rect &, const cv::Point&, int);

#endif 
