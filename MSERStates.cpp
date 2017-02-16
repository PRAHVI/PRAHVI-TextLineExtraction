#include "MSERStates.h"
#include <opencv2/opencv.hpp>
#include "CVGeometryUtils.h"
#include <math.h>
#include <fftw3.h>

using namespace std;
using namespace cv;


void MSERStates::generateHistogram(int siteID, const MSERState &siteState, vector<int> &histogram) {
	Rect mser = msers[siteID];
	Point center = getRectCenter(mser);
	int radius = siteState.scale.first / 2;
	double orientation = siteState.orientation;

	Point p1, p2;
	LineIterator center_axis;

	getLinePointsThroughRegionCenterAtDegree(orientation, center, radius, p1, p2);	
	
	center_axis = LineIterator(img, p1, p2);

	/* Init histogram to zeros */
	histogram.resize(center_axis.count, 0);

	/* Counts # cc perpendicular to the axis line */
	for (int i = 0; i < center_axis.count; i++, center_axis++) {
		Point poi = center_axis.pos();
		Point perp_p1, perp_p2;

		getPerpendicularLinePoints(radius, orientaion, poi, center_axis.count, perp_p1, perp_p2);

		for (Rect r : msers) {
			if (isrRctInLine(r, perp_p1, perp_p2)) 
				histogram[i]++;
		}
	}


}


void MSERStates::encodeLabelToState(int label, MSERState &s) {
	int k_orientation = label / qunatized_scale_factor;
	int scale_index = label % quantized_scale_factor;

	s.orientation = static_cast<double>(k_orientation) * M_PI / static_cast<double>(quantized_orientation_factor);
	s.scale = SCALES[scale_index];
}


double MSERStates::dataCost(int siteID, int label) {
	MSERState siteState;
	vector<int> histogram;

	encodeLabelToState(label, stateState);
	
	generateHistogram(siteID, siteState, histogram);

	/* FFT */

	/* Since original diameter of projection region shrinks at
	 * different angles, we must change the k phase to a value that
	 * makes N/k == histogram.size()/k_new more or less invariant 
	 */
	fftw_complex *x, *X;
	fftw_plan p;

	double k_new = static_cast<double>(histogram.size()) * static_cast<double>(siteState.scale.second) / static_cast<double>(siteState.scale.first);
	
	int k;
	if ( (k_new - static_cast<double>(static_cast<int>(k_new)) ) >  0.5)
		k = ciel(k_new);
	else
		k = floor(k_new);

	x = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * histogram.size());
	X = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * histogram.size());
	p = fftw_plan_dft(histogram.size(), in, out, FFTW_FORWARD, FFTW_ESTIMATE);
	
	for (int i = 0; i < histogram.size(); i++) {
		x[i][0] = static_cast<double>(histogram[i]);
		x[i][1] = 0.0;
	}

	fftw_execute(p);

	double mag_X0_squared = pow(X[0][0], 2) + pow(X[0][1], 2);
	double mag_Xk_squared = pow(X[0][0], 2) + pow(X[0][1], 2);
	double v_p1 = -log10(mag_Xk_squared / mag_X0_squared);

	/* V_2 Calculation: number of non-zero vals in the projection profile */
	int len_set_x_neq_0 = 0;
	
	for (int i = 0; i < histogram.size(); i++) {
		if (histogram[i] != 0) 
			len_set_x_neq_0++;
		
	double v_p2 = log10(static_cast<double>(len_set_x_neq_0) / static_cast<double>(histogram.size()));


	return (LAMBDA * v_p1) + ((1 - LAMBDA) * v_p2);

}

double MSERStates::smoothCost(int siteID1, int sitdID2, int l1, int l2) {
	Rect mser1 = msers[siteID1], mser2 = msers[siteID2];
	MSERtate s1, s2;
	Point site1_center = getRectCenter(mser1), site2_center = getRectCenter(mser2);
	
	encodeLabelToState(l1, s1);
	encodeLabelToState(l2, s2);

	double distance_squared_site1_site2 =  pow( site1_center.x - site2_center.x, 2) + 
										   pow( site1_center.y - site2_center.y, 2) ;
	double u_1_2 = stateDistance(l1, l2);

	double s1_squared = pow(static_cast<double>(s1.scale.first) / static_cast<double>(s1.scale.second), 2);
	double s2_squared = pow(static_cast<double>(s1.scale.first) / static_cast<double>(s1.scale.second), 2);

	return u_1_2 * exp(- (BETA * distance_squared_site1_site2) / ( s1_squared + s2_squared ) ); 
}

double MSERtates::stateDistance(int label1, int label2) {
	int s1_index, s2_index, theta1_index, theta2_index;
	s1_index = label_1 % quantized_scale_factor;
	s2_index = label_2 % quantized_scale_factor;
	theta1_index = (theta1_index / quantized_scale_factor) % ((quantized_orientation_factor / 2) + 1);
	theta2_index = (theta2_index / quantized_scale_factor) % ((quantized_orientation_factor/ 2) + 1);
 
	return abs(s1_index - s2_index) + abs(theta1_index - theta2_index);

}
