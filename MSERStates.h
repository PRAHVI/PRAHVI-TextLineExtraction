/* Name: MSERStates.h
 * Auther: Abe Millan
 * Description: This class enclouses all the functions that the connected components will need 
 *				in onder to run the label graph cut optimizer.
 */

#ifndef CCSTATE_H
#define CCSTATE_H

#include <utility> /* Pair */
#include <vector>

#define Pair std::make_pair
typedef std::vector<std::pair <int, int> > QuantScales;
typedef std::pair<int, int> Scale;


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
			   double BETA = 0.125) : img(image) 
									  msers(msers) 
									  quantized_orientation_factor(qo)
									  quantized_scale_factor(qs) {};
	
	static QuantScales SCALES[] = { Pair(64, 5),  /* 12.8 */
									Pair(64, 4),  /* 16.0 */
									Pair(64, 3),  /* 21.3 */
									Pair(128, 5), /* 25.6 */ 
									Pair(128, 4), /* 32.0 */
									Pair(128, 3), /* 42.7 */
									Pair(256, 5), /* 51.2 */
									Pair(256, 4), /* 64.0 */
									Pair(256, 3), /* 85.3 */ 
									Pair(256, 2)  /* 128.0 */
							      }

	double dataCost(int siteID, int label);	
	double smoothCost(int siteID1, int sitdID2, int l1, int l2); 


private:

	cv::Mat &img;
	std::vector<cv::Rect> &msers;
	int quantized_orientation_factor;
	int quantized_scales_factor;
	double LAMBDA;
	double BETA;


	void encodeLabelToState(int label, MSERState& s);
	double stateDistance(int label1, int label2);
	double getOrientationFromQuantize(int k);
	void generateHistogram(int siteID, const MSERState &siteState, std::vector<int> &histogram);
};

#endif



// Flow
GCoptimizationGeneralGraph gc = GCoptimizationGeneralGraph(num_ccs, num_label);

// Set the data costs funtions
...
// May want to set the labels randomly
gc.setLabelOrder(true) // This visits the ccs in random

gc.expansion();

for (r in cv::Rect) {
	gc.whatLabel(getSite(r));


// Cluster
	
