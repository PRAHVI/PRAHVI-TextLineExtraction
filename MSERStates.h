/* Name: MSERStates.h
 * Auther: Abe Millan
 * Description: This class enclouses all the functions that the connected components will need 
 *				in order to run the label graph cut optimizer.
 */

#ifndef CCSTATE_H
#define CCSTATE_H

#include <opencv2/opencv.hpp>
#include <utility> /* Pair */
#include <vector>
#include <unordered_map>
#include <set>

#define Pair std::make_pair
typedef std::vector<std::pair <int, int> > QuantScales;
typedef std::pair<int, int> Scale;
typedef std::unordered_map<int, set<int> > NeighborMap;
typedef std::set<int> Neighbors;
typedef std::unordered_map<string, int> MSERIDMap;


typedef struct MSERState {
	double orientation;
	Scale scale;
} MSERState;


class MSERStates {

public:

	MSERStates(cv::Mat &image, 
	           std::vector<cv::Rect> &msers, 
	           int qo = 32, 
	           int qs = 10, 
	           double LAMBDA = 0.5,
	           double BETA = 0.125) : img(image), 
	                                  msers(msers),
	                                  quantized_orientation_factor(qo),
	                                  quantized_scale_factor(qs), 
	                                  LAMBDA(LAMBDA),
	                                  BETA(BETA) {}
	
	const QuantScales SCALES = { Pair(64, 5),  /* 12.8 */
                                 Pair(64, 4),  /* 16.0 */
                                 Pair(64, 3),  /* 21.3 */
                                 Pair(128, 5), /* 25.6 */ 
                                 Pair(128, 4), /* 32.0 */
                                 Pair(128, 3), /* 42.7 */
                                 Pair(256, 5), /* 51.2 */
                                 Pair(256, 4), /* 64.0 */
                                 Pair(256, 3), /* 85.3 */ 
                                 Pair(256, 2)  /* 128.0 */
							   };

	double dataCost(int siteID, int label);	
	double smoothCost(int siteID1, int siteID2, int l1, int l2); 
	cv::Mat& draw();
	void updateState(int siteID, int label);
	void generateDelaunayNeighbors(int *numNeighors, int **neighborsIndexes);


private:

	cv::Mat &img;
	std::vector<cv::Rect> &msers;
	std::vector<MSERState> msers_states;
	int quantized_orientation_factor;
	int quantized_scale_factor;
	double LAMBDA;
	double BETA;


	void encodeLabelToState(int label, MSERState& s);
	double stateDistance(int label1, int label2);
	double getOrientationFromQuantize(int k);
	void generateHistogram(int siteID, 
	                       const MSERState &siteState, 
						   std::vector<int> &histogram);
};

#endif
