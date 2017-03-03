#include "MSERStates.h"
#include "CVGeometryUtils.h"
#include <math.h>
#include <fftw3.h>

using namespace std;
using namespace cv;


/* 
 * Function: MSERStates::generateHistogram
 *
 * Description: calculates the histogram of the projection profile of 
 *				a mser connected component with a state of siteState
 *
 */	
void MSERStates::generateHistogram(int siteID, const MSERState &siteState, vector<int> &histogram) {
	Rect mser = this->msers[siteID];
	Point center = getRectCenter(mser);
	int radius = siteState.scale.first / 2;
	double orientation = siteState.orientation;
	vector<Rect*> msers_in_region;

	/* First iterate and find the msers in the circular region so we dont waste time later */
	for (Rect r : msers) {
		if (isRectInCircle(r, center, radius))
			msers_in_region.push_back(&r);
	}

	Point p1, p2;

	getLinePointsThroughRegionCenterAtDegree(orientation, center, radius, p1, p2);	
	
	LineIterator center_axis(this->img, p1, p2);

	/* Init histogram to zeros */
	histogram.resize(center_axis.count, 0);

	/* Counts # cc perpendicular to the axis line */
	for (int i = 0; i < center_axis.count; i++, center_axis++) {
		Point perp_p1, perp_p2;

		getPerpendicularLinePoints(i, radius, orientation, center_axis, perp_p1, perp_p2);

		for (Rect r : msers) {
			if (isRectInLine(r, perp_p1, perp_p2)) 
				histogram[i]++;
		}
	}


}

/*
 * Function: MSERStates::encodeLabelToState
 *
 * Description: given a labeling, it converts it into an 
 *					- orientation angle 

 *					- scale(int N, int k)
 */
void MSERStates::encodeLabelToState(int label, MSERState &s) {
	int k_orientation = label / quantized_scale_factor;
	int scale_index = label % quantized_scale_factor;

	s.orientation = (static_cast<double>(k_orientation) * M_PI / 
					static_cast<double>(quantized_orientation_factor)) * (180.0 / M_PI);
	s.scale = SCALES[scale_index];
}


/*
 * Function: MSERStates::dataCost
 *
 * Description: Assigns a energy to a site given a label.
 *				This function is passed to the GCOptimization Library
 *
 *				The energy calculation is defined in the following paper:
 *					Text-Line Detection in Camera-Captured Document Images
 *					Using the State Estimation of Connected Components
 *
 *					url: http://ieeexplore.ieee.org/document/7563454/
 */
int MSERStates::dataCost(int siteID, int label) {
	MSERState siteState;
	vector<int> histogram;

	encodeLabelToState(label, siteState);
	
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
		k = ceil(k_new);
	else
		k = floor(k_new);

	x = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * histogram.size());
	X = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * histogram.size());
	p = fftw_plan_dft_1d(histogram.size(), x, X, FFTW_FORWARD, FFTW_ESTIMATE);
	
	for (int i = 0; i < histogram.size(); i++) {
		x[i][0] = static_cast<double>(histogram[i]);
		x[i][1] = 0.0;
	}

	fftw_execute(p);

	double mag_X0_squared = pow(X[0][0], 2) + pow(X[0][1], 2);
	double mag_Xk_squared = pow(X[k][0], 2) + pow(X[k][1], 2);
	double v_p1 = -log10(mag_Xk_squared / mag_X0_squared);

	/* V_2 Calculation: number of non-zero vals in the projection profile */
	int len_set_x_neq_0 = 0;
	
	for (int i = 0; i < histogram.size(); i++) {
		if (histogram[i] != 0) 
			len_set_x_neq_0++;
	}
		
	double v_p2 = log10(static_cast<double>(len_set_x_neq_0) / static_cast<double>(histogram.size()));


	return static_cast<int>((LAMBDA * v_p1) + ((1.0 - LAMBDA) * v_p2));

}

/*
 * Function: MSER::smoothCost
 *
 * Description: Assigns an enery to two sites with their respective labels.
 *				This function is passed to the GCOptimization Library
 *
 *				The energy calculation is defined in the following paper:
 *					Text-Line Detection in Camera-Captured Document Images
 *					Using the State Estimation of Connected Components
 *
 *					url: http://ieeexplore.ieee.org/document/7563454/
 *				
 */
int MSERStates::smoothCost(int siteID1, int siteID2, int l1, int l2) {
	Rect mser1 = this->msers[siteID1]; 
	Rect mser2 = this->msers[siteID2];
	MSERState s1, s2;
	Point site1_center = getRectCenter(mser1), site2_center = getRectCenter(mser2);
	
	encodeLabelToState(l1, s1);
	encodeLabelToState(l2, s2);

	double distance_squared_site1_site2 =  pow( site1_center.x - site2_center.x, 2) + 
										   pow( site1_center.y - site2_center.y, 2) ;
	double u_1_2 = stateDistance(l1, l2);

	double s1_squared = pow(static_cast<double>(s1.scale.first) / static_cast<double>(s1.scale.second), 2);
	double s2_squared = pow(static_cast<double>(s1.scale.first) / static_cast<double>(s1.scale.second), 2);

	return static_cast<int>(u_1_2 * exp(- (BETA * distance_squared_site1_site2) / ( s1_squared + s2_squared ) )); 
}

/*
 * Function MSERStates::stateDistance
 *
 * Description: Computes the distances between two states, which in this case will be 
 *				interpreted by int labels
 *
 *				The distance between states is calculated as follows:
 *				| f_p - f_q | = | s_p - s_q | + | theta_p - theta_q |
 *
 *				This comes from the following paper:
 *					Text-Line Detection in Camera-Captured Document Images
 *					Using the State Estimation of Connected Components
 *
 *					url: http://ieeexplore.ieee.org/document/7563454/
 */

double MSERStates::stateDistance(int label1, int label2) {
	int orientation_index1 = (label1 / quantized_scale_factor) % ((quantized_orientation_factor / 2) + 1);
	int scale_index1 = label1 % quantized_scale_factor;
	int orientation_index2 = (label2 / quantized_scale_factor) % ((quantized_orientation_factor / 2) + 1);
	int scale_index2 = label2 % quantized_scale_factor;


	return abs(scale_index1 - scale_index2) + abs(orientation_index1 - orientation_index2);

}

string pointToString(const Point &p) 
{
	return to_string(p.x) + "," + to_string(p.y);
}

/* 
 * Function MSERStates::generateDelaunayNeighbors
 * Description: Computes the neighbors for each point and allocates and stores them 
 *				in neighborIndexes 2D array. 
 *				IMPORTANT: This function is not responsible for deallocation of memory
 *
 */
void MSERStates::generateDelaunayNeighbors(int *numNeighbors, int **neighborsIndexes, int **neighborWeights) 
{
	Size size = this->img.size();
	Rect rect(0, 0, size.width, size.height);
	MSERIDMap mser_id_map;

	Subdiv2D subdiv(rect);

	vector<Point2f> points;
	vector<Vec6f> triangleList;
		
	for (int i = 0; i < this->msers.size(); i++) {
		Rect r = this->msers[i];
		Point rect_center = getRectCenter(r);
		points.push_back(rect_center);
		mser_id_map.insert(MSERIDMap::value_type(pointToString(rect_center), i));
	}

	for (int i = 0; i < points.size(); i++) {
		subdiv.insert(points[i]);
	}

	/////////////////////////////////////////////////////
	//				Neighbor Hashing                  //	
	////////////////////////////////////////////////////
	subdiv.getTriangleList(triangleList);
	vector<Point> pts(3);
	NeighborMap nm;
	
	for (int i = 0; i < triangleList.size(); i++) {
		Vec6f t = triangleList[i];
		pts[0] = Point(cvRound(t[0]), cvRound(t[1]));
		pts[1] = Point(cvRound(t[2]), cvRound(t[3]));
		pts[2] = Point(cvRound(t[4]), cvRound(t[5]));
		
		for (int j = 0; j < pts.size(); j++) {
			int n1_id = mser_id_map[pointToString(pts[ (j + 1) % pts.size()] )];
			int n2_id = mser_id_map[pointToString(pts[(j + 2) % pts.size()])];
			int id = mser_id_map[pointToString(pts[j])];

			NeighborMap::const_iterator got = nm.find(id);

			if (got == nm.end()) {
				Neighbors ns( {n1_id, n2_id} );
				nm.insert(NeighborMap::value_type(id, ns));
			} else {
				nm[id].insert(n1_id);
				nm[id].insert(n2_id);
			}
		}
	}

	//////////////////////////////////////////////////////
	//				ALLOCATION						   //
	////////////////////////////////////////////////////
	int num_sites = this->msers.size();	
	numNeighbors = new int[num_sites];
	neighborsIndexes = new int*[num_sites];
	neighborWeights = new int*[num_sites];


	for (int i = 0; i < num_sites; i++) {
		//TODO: Could be that nm[i] does not exist
		NeighborMap::const_iterator got = nm.find(i);

		if (got == nm.end()) {
			numNeighbors[i] = 0;
		} else {
			int arr_size = nm[i].size();
			numNeighbors[i] = arr_size;
			neighborsIndexes[i] = new int[arr_size];
			neighborWeights[i] = new int[arr_size];

			Neighbors::const_iterator it = nm[i].begin();
			for (int j = 0; it != nm[i].end(); it++, j++) {
				neighborsIndexes[i][j] = *it;
				neighborWeights[i][j] = 1;

			}
		}

		

	}

	for (int i = 0; i < num_sites; i++) {
		cout << "numNeighbors[" << i << "] = " << numNeighbors[i] << endl;
	}

}

int MSERStates::static_dataCost(int siteID, int label, void *object)
{
	return static_cast<MSERStates*>(object)->dataCost(siteID, label);
}

int MSERStates::static_smoothCost(int siteID1, int siteID2, int l1, int l2, void *object)
{
	return static_cast<MSERStates*>(object)->smoothCost(siteID1, siteID2, l1, l2);
}
