#include <opencv2/opencv.hpp>
#include "MSERStates.h"
#include <GCoptimization.h>
#include "CVGeometryUtils.h"
#include <iostream>

#define TESTING 1

using namespace std;
using namespace cv;

//void testDataCost();
void mserExtractor(Mat& image, vector<cv::Rect>& msers_bbox);

int main(int argc, char** argv) {

	//if (TESTING) 
//		testDataCost();

	Mat img = imread("test_images/0000.jpg", IMREAD_GRAYSCALE);
	vector<Rect> msers_bbox;
	mserExtractor(img, msers_bbox);

	int num_labels = 10 * 32;
	int *numNeighbors;
	int **neighborIndexes;
	int **neighborWeights;

	MSERStates image_states = MSERStates(img, msers_bbox);
	image_states.generateDelaunayNeighbors(numNeighbors, neighborIndexes, neighborWeights);

	for (int i = 0; i < msers_bbox.size(); i++) {
		cout << "numNeighbors[" << i << "] = " << numNeighbors[i] << endl;
	}

	cout << msers_bbox.size() << endl;
	try {
		GCoptimizationGeneralGraph *gc = new GCoptimizationGeneralGraph(msers_bbox.size(), num_labels);
		gc->setDataCost(MSERStates::static_dataCost, &image_states);
		gc->setSmoothCost(MSERStates::static_smoothCost, &image_states);

		gc->setAllNeighbors(numNeighbors, neighborIndexes, neighborWeights);
		//gc->setLabelOrder(true);
		//gc->expansion();


		delete gc;
	} catch (GCException e) {
		e.Report();
	}

	//gc.setDataCost(MSERStates::static_dataCost, &image_states);
	//gc.setSmoothCost(MSERStates::static_smoothCost, &image_states);
	//gc.setLabelOrder(true);
	//gc.expansion();
	
	//for (int i = 0; i < msers_bbox.size(); i++) {
	//	image_states.updateState(i, gc.whatLabel(i));
	//}

	//Mat states_image = image_states.draw();
	return 0;
}	

//void inputState(int event, int x, int y, int flags, void* userdata) {
//	if (event == EVENT_LBUTTONDOWN) {
//		void** realData = (void**) userdata;
//		int rand_siteID = *(int *) realData[0];
//		MSERStates *ms = (MSERStates *) realData[1]; 
//		int k_orientation, scale;
//
//		cout << "Input orientation (0 - 31): ";
//		cin >> k_orientation;
//		cout << endl;
//		cout << "Input Scale (0 - 9): ";
//		cin >> scale;
//		cout << endl;
//
//		double result = ms->dataCost(rand_siteID, 10*k_orientation + scale);
//
//		cout << "Energy (Remember the lower the better state estimate): " << result << endl;
//		cout << endl;
//		int key_code = waitKey(0);
//	}
//}
//
//void testDataCost() {
//	Mat image = imread("test_images/0063.jpg", IMREAD_GRAYSCALE);
//	vector<Rect> msers_bbox;	
//
//	mserExtractor(image, msers_bbox);
//
//	MSERStates ms = MSERStates(image, msers_bbox);
//
//	srand(time(NULL));
//	while(1) {
//		Mat image_one_cc, proj_profile;
//		image.copyTo(image_one_cc);
//
//		/* Show the selected CC */
//		int rand_siTeID = rand() % msers_bbox.size();
//		Point center = getRectCenter(msers_bbox[rand_siteID]);
//		circle(image_one_cc, center, 100, Scalar(0), 5);
//
//		namedWindow("Image", 1);
//		void *data[2] = { &rand_siteID, &ms };
//		imshow("Image", image_one_cc); 
//		setMouseCallback("Image", inputState, data);
//
//		//int key_code = waitKey(0);
//		//if (key_code == 1048689)
//	//		break;
//		
//		
//		
//	}
//}

void mserExtractor(Mat& image, vector<cv::Rect>& msers_bbox) {
	Ptr<MSER> mserExtractor = MSER::create();
	vector<vector<Point> > mser_contours;
	mserExtractor->detectRegions(image, mser_contours, msers_bbox);
}
