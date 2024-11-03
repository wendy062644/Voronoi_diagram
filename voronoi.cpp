#include <bits/stdc++.h>

using namespace std;

struct Point {
    double x, y;

    bool operator==(const Point& other) const {
        return x == other.x && y == other.y;
    }

    bool operator<(const Point& other) const {
        return x < other.x || x == other.x && y < other.y;
    }
};

struct Edge {
    Point start, end;
};

vector<Point> points;
vector<Edge> edges;

pair<Point, Point> calculatePerpendicularBisector(const Point& A, const Point& B) {
    double x = (A.x + B.x) / 2;
    double y = (A.y + B.y) / 2;

    Point startPoint, endPoint;

    if (B.y - A.y == 0) { // 水平
        startPoint.x = x;
        startPoint.y = 0;
        endPoint.x = x;
        endPoint.y = 600;
        return {startPoint, endPoint};
    }
    if (B.x - A.x == 0) { // 垂直
        startPoint.x = 0;
        startPoint.y = y;
        endPoint.x = 600;
        endPoint.y = y;
        return {startPoint, endPoint};
    }

    double m = -(1.0 / ((B.y - A.y) / (B.x - A.x)));

    if (m > 0) {
        if ((600 - y) / m > 600 - x) {
            endPoint.x = 600;
            endPoint.y = y + (600 - x) * m;
        } else {
            endPoint.x = x + (600 - y) / m;
            endPoint.y = 600;
        }

        if (y / m > x) {
            startPoint.x = 0;
            startPoint.y = y - x * m;
        } else {
            startPoint.x = x - y / m;
            startPoint.y = 0;
        }
    } else {
        if (-y / m > 600 - x) {
            endPoint.x = 600;
            endPoint.y = y + (600 - x) * m;
        } else {
            endPoint.x = x - y / m;
            endPoint.y = 0;
        }

        if ((y - 600) / m > x) {
            startPoint.x = 0;
            startPoint.y = y - x * m;
        } else {
            startPoint.x = x + (600 - y) / m;
            startPoint.y = 600;
        }
    }
    return {startPoint, endPoint};
}

Point calculateMidPoint(const Point& A, const Point& B) {
    return {(A.x + B.x) / 2.0, (A.y + B.y) / 2.0};
}

Point calculateIntersection(double slope1, Point mid1, double slope2, Point mid2) {
    if (isinf(slope1)) {
        return {mid1.x, slope2 * (mid1.x - mid2.x) + mid2.y};
    }
    if (isinf(slope2)) {
        return {mid2.x, slope1 * (mid2.x - mid1.x) + mid1.y};
    }

    // y - mid1.y = slope1 * (x - mid1.x)
    // y - mid2.y = slope2 * (x - mid2.x)
    double x = (slope2 * mid2.x - slope1 * mid1.x + mid1.y - mid2.y) / (slope2 - slope1);
    double y = slope1 * (x - mid1.x) + mid1.y;
    return {x, y};
}

void recursiveVoronoi(int L, int R) {
    if (R - L == 2) {
        pair<Point, Point> t = calculatePerpendicularBisector(points[L], points[L + 1]);
        edges.push_back({t.first, t.second});
        return;
    } else if (R - L == 3) {

        if(points[L+1].x == points[L].x && points[L+2].x == points[L].x
           || points[L+1].y == points[L].y && points[L+2].y == points[L].y
           || (points[L+1].y - points[L].y)*(points[L+2].x - points[L+1].x) == (points[L+2].y - points[L+1].y)*(points[L+1].x - points[L].x)) { // 垂直 or 水平
            recursiveVoronoi(L, L+2);
            recursiveVoronoi(L+1, R);
            return;
        }

        double slope1 = -(points[L+1].x - points[L].x) / (points[L+1].y - points[L].y);
        double slope2 = -(points[L+2].x - points[L+1].x) / (points[L+2].y - points[L+1].y);

        Point mid1 = calculateMidPoint(points[L], points[L+1]);
        Point mid2 = calculateMidPoint(points[L+1], points[L+2]);

        pair<Point, Point> bisector1 = calculatePerpendicularBisector(points[L], points[L+1]);
        pair<Point, Point> bisector2 = calculatePerpendicularBisector(points[L+1], points[L+2]);
        pair<Point, Point> bisector3 = calculatePerpendicularBisector(points[L], points[L+2]);
        pair<Point, Point> temp;

        Point intersection = calculateIntersection(slope1, mid1, slope2, mid2);

        if(intersection.x >= 0 && intersection.x <= 600 && intersection.y >= 0 && intersection.y <= 600) {
            double AB2 = pow(points[L].x - points[L+1].x, 2) + pow(points[L].y - points[L+1].y, 2);
            double BC2 = pow(points[L+1].x - points[L+2].x, 2) + pow(points[L+1].y - points[L+2].y, 2);
            double CA2 = pow(points[L].x - points[L+2].x, 2) + pow(points[L].y - points[L+2].y, 2);
            if(AB2 >= BC2 && AB2 >= CA2) {
                swap(bisector1, bisector3);
                swap(points[L+1], points[L+2]);
            }
            else if(BC2 >= AB2 && BC2 >= CA2) {
                swap(bisector2, bisector3);
                swap(points[L], points[L+1]);
            }
            if (calculateMidPoint(points[L], points[L+1]).x > intersection.x) {
                    bisector1.first = intersection;
            } else {
                bisector1.second = intersection;
            }
            if (calculateMidPoint(points[L+1], points[L+2]).x > intersection.x) {
                bisector2.first = intersection;
            } else {
                bisector2.second = intersection;
            }
            if(AB2 + BC2 <= CA2 || AB2 + CA2 <= BC2 || BC2 + CA2 <= AB2) { // 直角 or 鈍角
                if (calculateMidPoint(points[L], points[L+2]).x > intersection.x) {
                    bisector3.second = intersection;
                } else {
                    bisector3.first = intersection;
                }
            }
            else { // 銳角
                if (calculateMidPoint(points[L], points[L+2]).x > intersection.x) {
                    bisector3.first = intersection;
                } else {
                    bisector3.second = intersection;
                }
            }
            if(AB2 >= BC2 && AB2 >= CA2) {
                swap(points[L+1], points[L+2]);
            }
            else if(BC2 >= AB2 && BC2 >= CA2) {
                swap(points[L], points[L+2]);
            }
        }
        else {
            recursiveVoronoi(L, L+2);
            recursiveVoronoi(L+1, R);
            return;
        }

        edges.push_back({bisector1.first, bisector1.second});
        edges.push_back({bisector2.first, bisector2.second});
        edges.push_back({bisector3.first, bisector3.second});
        return;
    }

    int mid = (L + R) / 2;
    recursiveVoronoi(L, mid);
    recursiveVoronoi(mid, R);
}

void read_points(const string& filename) {
    ifstream infile(filename);
    double x, y;
    while (infile >> x >> y) {
        Point point = {x, y};
        if (find(points.begin(), points.end(), point) == points.end()) {
            points.push_back(point);
        }
    }
}

void write_edges(const string& filename) {
    ofstream outfile(filename);
    for (const auto& edge : edges) {
        outfile << edge.start.x << " " << edge.start.y << " " << edge.end.x << " " << edge.end.y << "\n";
    }
}

int main() {
    read_points("points.txt");
    sort(points.begin(), points.end());
    recursiveVoronoi(0, points.size());
    cout << points.size() << endl;
    write_edges("lines.txt");
    return 0;
}
