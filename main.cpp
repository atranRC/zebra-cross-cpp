#include <iostream>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
//#include "highgui.h"
#include <vector>
#include <cmath>
#include <iterator>
#include <opencv/cv.hpp>

using namespace std;
using namespace cv;

template<typename T>
std::ostream & operator<<(std::ostream & os, std::vector<T> vec)
{
    os<<"{ ";
    std::copy(vec.begin(), vec.end(), std::ostream_iterator<T>(os, " "));
    os<<"}";
    return os;
}

struct _Point {
    double x;
    double y;
};

struct Line {
    double m;
    double b;
};
Point median(vector<Point> arr) {
    vector<int> xs;
    vector<int> ys;
    for ( auto &i : arr ) {
        xs.push_back(i.x);
        ys.push_back(i.y);
    }
    auto first = xs.begin();
    auto last = xs.end();
    auto middle = first + (last - first) / 2;
    nth_element(first, middle, last);
    int x_median = *middle;

    first = ys.begin();
    last = ys.end();
    middle = first + (last - first) / 2;
    nth_element(first, middle, last);
    int y_median = *middle;

    return Point(x_median, y_median);
}

//vector<vector<double>> transpose(vector<)

Line lineCalc(double vx, double vy, double x0, double y0) {
    int scale = 10;
    double x1 = x0 + scale + vx;
    double y1 = y0 + scale + vy;
    Line l;
    l.m = (y1 - y0) / (x1 - x0);
    l.b = y1 - l.m * x1;
    return l;
}

double angle(_Point pt1, _Point pt2) {
    double inner_product = pt1.x * pt2.x + pt1.y + pt2.y;
    double len1 = hypot(pt1.x, pt1.y);
    double len2 = hypot(pt2.x, pt2.y);
    double a = acos(inner_product / (len1 * len2));
    return a * 180 / M_PI;
}

_Point lineIntersect(Line l1, Line l2) {
    double a1 = -l1.m;
    double b1 = 1;
    double c1 = b1;

    double a2 = -l2.m;
    double b2 = 1;
    double c2 = b2;

    double d = a1 * b2 - a2 * b1;
    double dx = c1 * b2 - c2 * b1;
    double dy = a1 * c2 - a2 * c1;
    _Point p;
    p.x = dx / d;
    p.y = dy / d;
    return p;
}

_Point process(Mat &im) {
    double y = im.size[0];
    double x = im.size[1];

    // TODO: Tweak these values and see changes
    double radius = 250;
    double thresh = 170;
    double bw_width = 170;

    vector<double> bxLeft;
    vector<double> byLeft;
    vector<double> bxRight;
    vector<double> byRight;
    vector<Point> boundedLeft;
    vector<Point> boundedRight;
    vector<Point> bxbyLeftArray;
    vector<Point> bxbyRightArray;

    Scalar lower(170, 170, 170);
    Scalar higher(255, 255, 255);
    Mat mask;
    inRange(im, lower, higher, mask);

    int erodeSize = (int) y / 30;
    Mat erodeStructure = getStructuringElement(MORPH_RECT, Point(erodeSize, 1));
    Mat erodeMat;
    erode(mask, erodeMat, erodeStructure);
    vector<vector<Point> > contours;
    vector<Vec4i> hierarchy;
    findContours(erodeMat, contours, hierarchy, RETR_TREE, CHAIN_APPROX_NONE);

    for ( size_t i = 0; i < contours.size(); i++ ) {
        Rect re = boundingRect(contours[i]);
        int bx = re.x, by = re.y, bw = re.width, bh = re.height;
        if (re.width > bw_width) {
            line(im, Point(bx, by), Point(bx + bw, by + bh), Scalar(0, 255, 0), 2);
            bxRight.push_back(bx + bw);
            byRight.push_back(by + bh);
            bxLeft.push_back(bx);
            byLeft.push_back(by);
            bxbyRightArray.push_back(Point(bx + bw, by + bh));
            bxbyLeftArray.push_back(Point(bx, by));
            circle(im, Point(bx, by), 5, Scalar(0,250,250), 2);
            circle(im, Point(bx + bw, by + bh), 5, Scalar(250,250,0), 2);
        }
    }
    Point medianL = median(bxbyLeftArray);
    Point medianR = median(bxbyRightArray);

    for (auto &i : bxbyLeftArray) {
        if (pow(medianL.x - i.x, 2) + pow(medianL.y - i.y, 2) < pow(radius, 2)) {
            boundedLeft.push_back(i);
        }
    }

    for (auto &i : bxbyRightArray) {
        if (pow(medianR.x - i.x, 2) + pow(medianR.y - i.y, 2) < pow(radius, 2)) {
            boundedRight.push_back(i);
        }
    }

    bxLeft.clear();
    byLeft.clear();
    for ( auto &i : boundedLeft ) {
        bxLeft.push_back(i.x);
        byLeft.push_back(i.y);
    }
    bxRight.clear();
    byRight.clear();
    for ( auto &i : boundedRight ) {
        bxRight.push_back(i.x);
        byRight.push_back(i.y);
    }

    cout<<bxLeft;
    transpose(bxLeft, bxLeft);
    cout<<bxLeft;

//    findHomography()

//    transpose(bxRight, bxRight);

    _Point D;
    D.x = 2;
    D.y = 3;
    return D;
}

//    bxLeftT = np.array([bxLeft]).transpose()
//    bxRightT = np.array([bxRight]).transpose()
//
//    # run ransac for LEFT
//    model_ransac = linear_model.RANSACRegressor(linear_model.LinearRegression())
//    ransacX = model_ransac.fit(bxLeftT, byLeft)
//    inlier_maskL = model_ransac.inlier_mask_  # right mask
//
//    # run ransac for RIGHT
//    ransacY = model_ransac.fit(bxRightT, byRight)
//    inlier_maskR = model_ransac.inlier_mask_  # left mask
//
//    # draw RANSAC selected circles
//    for i, element in enumerate(boundedRight[inlier_maskR]):
//        # print(i,element[0])
//        cv2.circle(im, (element[0], element[1]), 10, (250, 250, 250), 2)  # circles -> right line
//
//    for i, element in enumerate(boundedLeft[inlier_maskL]):
//        # print(i,element[0])
//        cv2.circle(im, (element[0], element[1]), 10, (100, 100, 250), 2)  # circles -> right line
//
//    # 6. Calcuate the intersection point of the bounding lines
//    # unit vector + a point on each line
//    vx, vy, x0, y0 = cv2.fitLine(boundedLeft[inlier_maskL], cv2.DIST_L2, 0, 0.01, 0.01)
//    vx_R, vy_R, x0_R, y0_R = cv2.fitLine(boundedRight[inlier_maskR], cv2.DIST_L2, 0, 0.01, 0.01)
//
//    # get m*x+b
//    m_L, b_L = lineCalc(vx, vy, x0, y0)
//    m_R, b_R = lineCalc(vx_R, vy_R, x0_R, y0_R)
//
//    # calculate intersention
//    intersectionX, intersectionY = lineIntersect(m_R, b_R, m_L, b_L)
//
//    # 7. draw the bounding lines and the intersection point
//    m = radius * 10
//    if (intersectionY < H / 2):
//        cv2.circle(im, (int(intersectionX), int(intersectionY)), 10, (0, 0, 255), 15)
//        cv2.line(im, (x0 - m * vx, y0 - m * vy), (x0 + m * vx, y0 + m * vy), (255, 0, 0), 3)
//        cv2.line(im, (x0_R - m * vx_R, y0_R - m * vy_R), (x0_R + m * vx_R, y0_R + m * vy_R), (255, 0, 0), 3)
//
//    # 8. calculating the direction vector
//    POVx = W / 2  # camera POV - center of the screen
//    POVy = H / 2  # camera POV - center of the screen
//
//    Dx = -int(intersectionX - POVx)  # regular x,y axis coordinates
//    Dy = -int(intersectionY - POVy)  # regular x,y axis coordinates
//
//    # focal length in pixels = (image width in pixels) * (focal length in mm) / (CCD width in mm)
//    focalpx = int(W * 4.26 / 6.604)  # all in mm
//
//    end = timeit.timeit()  # STOP TIMER
//    time_ = end - start
//
//    print('DELTA (x,y from POV):' + str(Dx) + ',' + str(Dy))
//    return im, Dx, Dy


int main() {
    Mat img = imread("/home/nathnel/PycharmProjects/innovatefpgacv/zebra_crossing/images/main.jpeg");

    _Point x = process(img);

//    cout<<x;
    imshow("Image", img);
//
    waitKey(0);

    return 0;
}
