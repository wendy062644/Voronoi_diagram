#include <bits/stdc++.h>
using namespace std;

struct Point {
    double x, y;
};

vector<pair<Point, Point>> E;

bool compare(Point a, Point b)
{
	return (a.x < b.x) || (a.x == b.x && a.y < b.y);
}

Point calculateMidPoint(const Point& A, const Point& B) {
    return {(A.x + B.x) / 2.0, (A.y + B.y) / 2.0};
}

Point calculateIntersection(double slope1, Point mid1, double slope2, Point mid2) {
    // 方程式：y - mid1.y = slope1 * (x - mid1.x)
    // y - mid2.y = slope2 * (x - mid2.x)
    double x = (slope2 * mid2.x - slope1 * mid1.x + mid1.y - mid2.y) / (slope2 - slope1);
    double y = slope1 * (x - mid1.x) + mid1.y;
    return {x, y};
}

pair<Point, Point> calculatePerpendicularBisector(const Point& A, const Point& B) {
    double x = (A.x + B.x) / 2;
    double y = (A.y + B.y) / 2;

    Point startPoint, endPoint;

    // 水平或垂直
    if(B.y - A.y == 0) {
        startPoint.x = 0;
        startPoint.x = A.y;
        endPoint.x = 600;
        endPoint.y = A.y;
        return {startPoint, endPoint};
    }
    if(B.x - A.x == 0) {
        startPoint.x = A.x;
        startPoint.y = 0;
        endPoint.x = A.x;
        endPoint.y = 600;
        return {startPoint, endPoint};
    }

    double m = -(1.0/((B.y - A.y) / (B.x - A.x))); // 垂直平分線斜率

    if(m > 0) { // 正斜率
        if((600-y)/m > 600-x) { // 中點延伸會碰到x邊界
            startPoint.x = 600;
            startPoint.y = y + (600-x)*m;
        }
        else { // 碰到y邊界
            startPoint.x = x + (600-y)/m;
            startPoint.y = 600;
        }

        if((y-0)/m > x) { // 中點延伸會碰到x邊界
            endPoint.x = 0;
            endPoint.y = y + -x*m;
        }
        else { // 碰到y邊界
            endPoint.x = x + -y/m;
            endPoint.y = 0;
        }
    }
    else {
        if((-y)/m > 600-x) { // 中點延伸會碰到x邊界
            startPoint.x = 600;
            startPoint.y = y + (600-x)*m;
        }
        else { // 碰到y邊界
            startPoint.x = x + (-y)/m;
            startPoint.y = 0;
        }

        if((y-600)/m > x) { // 中點延伸會碰到x邊界
            endPoint.x = 0;
            endPoint.y = y - x*m;
        }
        else { // 碰到y邊界
            endPoint.x = x + (600-y)/m;
            endPoint.y = 600;
        }
    }
    //cout << x << '\t' << y << endl;
    //cout << startPoint.x << '\t' << startPoint.y << '\t' << endPoint.x << '\t' << endPoint.y << '\t' << m << endl;
    return {startPoint, endPoint};
}

void Voronoi_Diagram (vector<Point> &v, int L, int R) {
    if(R-L == 2) {
        E.push_back(calculatePerpendicularBisector(v[L], v[L+1]));
        return;
    }
    else if(R-L == 3) {

        // 計算兩條垂直平分線的交點
        double slope1 = -(v[L+1].x - v[L].x) / (v[L+1].y - v[L].y);
        double slope2 = -(v[L+2].x - v[L+1].x) / (v[L+2].y - v[L+1].y);

        Point mid1 = calculateMidPoint(v[L], v[L+1]);
        Point mid2 = calculateMidPoint(v[L+1], v[L+2]);

        // 對於三個點，計算三條垂直平分線並求交點
        pair<Point, Point> bisector1 = calculatePerpendicularBisector(v[L], v[L+1]);
        pair<Point, Point> bisector2 = calculatePerpendicularBisector(v[L+1], v[L+2]);
        pair<Point, Point> bisector3 = calculatePerpendicularBisector(v[L], v[L+2]);

        Point intersection = calculateIntersection(slope1, mid1, slope2, mid2);

        if(intersection.x > (v[L].x+v[L+1].x)/2) {
            bisector1.first.x = intersection.x;
            bisector1.first.y = intersection.y;
        }
        else {
            bisector1.second.x = intersection.x;
            bisector1.second.y = intersection.y;
        }

        if(intersection.x > (v[L+1].x+v[L+2].x)/2) {
            bisector2.first.x = intersection.x;
            bisector2.first.y = intersection.y;
        }
        else {
            bisector2.second.x = intersection.x;
            bisector2.second.y = intersection.y;
        }

        if(intersection.x > (v[L].x+v[L+2].x)/2) {
            bisector3.first.x = intersection.x;
            bisector3.first.y = intersection.y;
        }
        else {
            bisector3.second.x = intersection.x;
            bisector3.second.y = intersection.y;
        }
        // 輸出交點，這裡可以選擇使用交點進行更多處理
        //cout << "Intersection at: (" << intersection.x << ", " << intersection.y << ")" << endl;

        // 也可以選擇將這些垂直平分線加到結果中
        E.push_back(bisector1);
        E.push_back(bisector2);
        E.push_back(bisector3);

        return;
    }
    Voronoi_Diagram(v, L, (L+R)/2); // left
    Voronoi_Diagram(v, (L+R)/2, R); // right
    // merge
}

int main() {
    int n;
    while(cin >> n && n) {
        vector<Point> v(n);
        for(int i = 0; i < n; ++i) {
            cin >> v[i].x >> v[i].y;
        }
        sort(v.begin(), v.end(), compare);
        for(Point p : v)
            cout << p.x << ' ' << p.y << endl;
        Voronoi_Diagram(v, 0, n);
        cout << E.size() << endl;
        for(pair<Point, Point> p : E) {
            printf("E %lf %lf %lf %lf\n", p.first.x, p.first.y, p.second.x, p.second.y);
        }
        E.clear();
    }
    return 0;
}
