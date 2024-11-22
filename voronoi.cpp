#include <bits/stdc++.h>
using namespace std;

struct Point {
    double x, y;

    bool operator==(const Point& other) const {
        return x == other.x && y == other.y || fabs(x - other.x) < 0.01 && fabs(y - other.y) < 0.01;
    }

    Point& operator=(const Point& other) {
        if (this != &other) {
            x = other.x;
            y = other.y;
        }
        return *this;
    }

    bool operator<(const Point& other) const {
        return x < other.x || x == other.x && y < other.y;
    }

    bool operator!=(const Point& other) const {
        return fabs(x - other.x) >= 0.01 && fabs(y - other.y) >= 0.01;
    }
};

struct Edge {
    Point start, end, A, B;
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

// 判斷點是否在基準線的左側
double crossProduct(const Point& p, const Point& a, const Point& b) {
    return (b.x - a.x) * (p.y - a.y) - (b.y - a.y) * (p.x - a.x);
}

// 判斷點 P 是否在線段 AB 上
bool isOnSegment(const Point& A, const Point& B, const Point& P) {

    if(fabs(A.x - B.x) < 0.01 && abs(P.x - A.x) < 0.01) return true;
    if(fabs(A.y - B.y) < 0.01 && abs(P.y - A.y) < 0.01) return true;

    double crossProduct = (B.x - A.x) * (P.y - A.y) - (B.y - A.y) * (P.x - A.x);

    if (fabs(crossProduct) > 0.01) {
        return false;
    }

    return min(A.x, B.x) <= P.x && P.x <= max(A.x, B.x) &&
           min(A.y, B.y) <= P.y && P.y <= max(A.y, B.y);
}

pair<Point, Point> adjustPointToBoundary(const Point& p1, const Point& p2) {
    Point startPoint = p1, endPoint = p2;
    if (p2.y - p1.y == 0) { // 水平
        startPoint.x = 0;
        startPoint.y = p2.y;
        endPoint.x = 600;
        endPoint.y = p2.y;
        return {startPoint, endPoint};
    }
    if (p2.x - p1.x == 0) { // 垂直
        startPoint.x = p2.x;
        startPoint.y = 0;
        endPoint.x = p2.x;
        endPoint.y = 600;
        return {startPoint, endPoint};
    }
    double m = (p1.y-p2.y) / (p1.x-p2.x);

    if (p2.x > 600) {  // 右邊界
        endPoint.x = 600;
        endPoint.y = p1.y + (600 - p1.x) * m;
    } else if (p2.x < 0) {  // 左邊界
        endPoint.x = 0;
        endPoint.y = p1.y + (0 - p1.x) * m;
    }

    if (endPoint.y > 600) {  // 上邊界
        endPoint.y = 600;
        endPoint.x = p1.x + (600 - p1.y) / m;
    } else if (endPoint.y < 0) {  // 下邊界
        endPoint.y = 0;
        endPoint.x = p1.x - (p1.y / m);
    }
    if (p1.x > 600) {
        startPoint.x = 600;
        startPoint.y = p2.y + (600 - p2.x) * m;
    } else if (p1.x < 0) {
        startPoint.x = 0;
        startPoint.y = p2.y + (0 - p2.x) * m;
    }

    if (startPoint.y > 600) {
        startPoint.y = 600;
        startPoint.x = p2.x + (600 - p2.y) / m;
    } else if (startPoint.y < 0) {
        startPoint.y = 0;
        startPoint.x = p2.x - (p2.y / m);
    }
    return {startPoint, endPoint};
}

vector<Point> convexHull(Point* begin, Point* end) {
    vector<Point> points(begin, end + 1);

    if (points.size() < 4) return points;

    sort(points.begin(), points.end());

    vector<Point> hull;

    for (const auto& point : points) {
        while (hull.size() >= 2 && crossProduct(hull[hull.size() - 2], hull[hull.size() - 1], point) <= 0) {
            hull.pop_back();
        }
        hull.push_back(point);
    }

    size_t t = hull.size();
    for (int i = points.size() - 2; i >= 0; --i) {
        while (hull.size() > t && crossProduct(hull[hull.size() - 2], hull[hull.size() - 1], points[i]) <= 0) {
            hull.pop_back();
        }
        hull.push_back(points[i]);
    }

    hull.pop_back();

    return hull;
}

double distanceToSegment(const Point& A, const Point& B, const Point& P) {
    double ABx = B.x - A.x;
    double ABy = B.y - A.y;
    double APx = P.x - A.x;
    double APy = P.y - A.y;

    double dotProduct = ABx * APx + ABy * APy;
    double lengthSquared = ABx * ABx + ABy * ABy;

    double t = max(0.0, min(1.0, dotProduct / lengthSquared));
    double projX = A.x + t * ABx;
    double projY = A.y + t * ABy;

    return sqrt((P.x - projX) * (P.x - projX) + (P.y - projY) * (P.y - projY));
}

void Swap(Point &a, Point &b) {
    if(a.y > b.y) swap(a, b);
}

vector<Edge> two_point(const Point A, const Point B) {
    pair<Point, Point> t = calculatePerpendicularBisector(A, B);
    Swap(t.first, t.second);
    Point a = A, b = B;
    Swap(a, b);
    return vector<Edge>{Edge{t.first, t.second, a, b}};
}

vector<Edge> trangle(const Point A, const Point B, const Point C, Point &t) {
    if(B.x == A.x && C.x == A.x
        || B.y == A.y && C.y == A.y
        || (B.y - A.y)*(C.x - B.x) == (C.y - B.y)*(B.x - A.x)) { // 垂直 or 水平
        return vector<Edge>{
            Edge{two_point(A, B)[0].start, two_point(A, B)[0].end, A, B},
            Edge{two_point(B, C)[0].start, two_point(B, C)[0].end, B, C}
        };
    }
    double slope1 = -(B.x - A.x) / (B.y - A.y);
    double slope2 = -(C.x - B.x) / (C.y - B.y);

    Point mid1 = calculateMidPoint(A, B);
    Point mid2 = calculateMidPoint(B, C);

    pair<Point, Point> bisector1 = calculatePerpendicularBisector(A, B);

    t = calculateIntersection(slope1, mid1, slope2, mid2);
    if(t.x >= 0 && t.x <= 600 && t.y >= 0 && t.y <= 600) {
        if(distanceToSegment(bisector1.first, C, t) > distanceToSegment(bisector1.second, C, t)) {
            bisector1.first = t;
        }else{
            bisector1.second = t;
        }
        return vector<Edge>{
            Edge{bisector1.first, bisector1.second, A, B}
            //Edge{bisector2.first, bisector2.second, b, c},
            //Edge{bisector3.first, bisector3.second, a, c}
        };
    }
    return vector<Edge>{};
}

Point IntersectionPoint(const Edge& A, const Edge& B) {
    double slope1 = (A.end.y - A.start.y) / (A.end.x - A.start.x);
    double slope2 = (B.end.y - B.start.y) / (B.end.x - B.start.x);

    double x = ((slope1 * A.start.x - A.start.y) - (slope2 * B.start.x - B.start.y)) / (slope1 - slope2);
    double y = slope1 * (x - A.start.x) + A.start.y;

    return Point{x, y};
}

pair<Point, Point> HyperPlane(const vector<Point>& leftConvexHull, const vector<Point>& rightConvexHull) {
    for(int i = 0; i < leftConvexHull.size(); i++) {
        for(int j = 0; j < rightConvexHull.size(); j++) {
            bool b = 1;
            for(int k = 0; k < leftConvexHull.size() && b; k++) {
                if(crossProduct(leftConvexHull[i], rightConvexHull[j], leftConvexHull[k]) < 0) {
                    b = 0;
                }
            }
            for(int k = 0; k < rightConvexHull.size() && b; k++) {
                if(crossProduct(leftConvexHull[i], rightConvexHull[j], rightConvexHull[k]) < 0) {
                    b = 0;
                }
            }
            if(b) return {leftConvexHull[i], rightConvexHull[j]};
        }
    }
}

optional<Point> doLinesIntersect(const Edge& A, const Edge& B) {
    // 計算每個點與線段的相對位置
    double d1 = crossProduct(B.start, B.end, A.start);
    double d2 = crossProduct(B.start, B.end, A.end);
    double d3 = crossProduct(A.start, A.end, B.start);
    double d4 = crossProduct(A.start, A.end, B.end);

    // 判斷是否相交
    if (d1 * d2 < 0 && d3 * d4 < 0) {
        // 計算交點
        double a1 = A.end.y - A.start.y;
        double b1 = A.start.x - A.end.x;
        double c1 = a1 * A.start.x + b1 * A.start.y;

        double a2 = B.end.y - B.start.y;
        double b2 = B.start.x - B.end.x;
        double c2 = a2 * B.start.x + b2 * B.start.y;

        double determinant = a1 * b2 - a2 * b1;

        if (determinant != 0) { // 兩條線段相交
            double x = (b2 * c1 - b1 * c2) / determinant;
            double y = (a1 * c2 - a2 * c1) / determinant;
            Point intersection = {x, y};

            // 確保交點在兩條線段的範圍內
            if (isOnSegment(A.start, A.end, intersection) && isOnSegment(B.start, B.end, intersection)) {
                return intersection;
            }
        }
    }

    // 判斷是否在端點上
    if (d1 == 0 && isOnSegment(B.start, B.end, A.start)) return A.start;
    if (d2 == 0 && isOnSegment(B.start, B.end, A.end)) return A.end;
    if (d3 == 0 && isOnSegment(A.start, A.end, B.start)) return B.start;
    if (d4 == 0 && isOnSegment(A.start, A.end, B.end)) return B.end;

    return nullopt;
}

vector<Edge> constructConvexHull(vector<Edge>& voronoiEdges, bool Left) {
    vector<Edge> Convexhull;

    for(int i = 0; i < voronoiEdges.size(); i++) {
        if(voronoiEdges[i].start.x == voronoiEdges[i].end.x || voronoiEdges[i].start.y == voronoiEdges[i].end.y) {
            Convexhull.push_back(voronoiEdges[i]);
            continue;
        }
        if(voronoiEdges[i].start.y > voronoiEdges[i].end.y) {
            swap(voronoiEdges[i].start, voronoiEdges[i].end);
        }
        double m = (voronoiEdges[i].end.y - voronoiEdges[i].start.y)/(voronoiEdges[i].end.x - voronoiEdges[i].start.x);
        if(voronoiEdges[i].start.y == 0 || voronoiEdges[i].start.x == 0 || voronoiEdges[i].start.y == 600 || voronoiEdges[i].start.x == 600 ||
           voronoiEdges[i].end.y == 0 || voronoiEdges[i].end.x == 0 || voronoiEdges[i].end.y == 600 || voronoiEdges[i].end.x == 600) {
            Convexhull.push_back(voronoiEdges[i]);
        }
    }
    return Convexhull;
}

vector<Point> ConvexHull_Point(const vector<Edge>& voronoiEdges) {
    auto compare = [](const Point& p1, const Point& p2) {
        return p1.y > p2.y || p1.y == p2.y && p1.x > p2.x;
    };
    set<Point, decltype(compare)> ConvexHull(compare);
    for(Edge i : voronoiEdges) {
        ConvexHull.insert(i.A);
        ConvexHull.insert(i.B);
    }
    vector<Point> result;
    for(Point i : ConvexHull) {
        result.push_back(i);
    }
    return result;
}

double calculateslope(Edge E) {
    return (E.start.y - E.end.y)/(E.start.x - E.end.x);
}

vector<Edge> recursiveVoronoi(int L, int R) {
    if(L+1 == R) return vector<Edge>{};
    vector<Edge> voronoi;
    if (R - L == 2) {
        pair<Point, Point> t = calculatePerpendicularBisector(points[L], points[L+1]);
        Swap(t.first, t.second);
        Point A = points[L], B = points[L+1];
        Swap(A, B);
        return vector<Edge>{Edge{t.first, t.second, A, B}};
    } else if (R - L == 3) {
        if(points[L+1].x == points[L].x && points[L+2].x == points[L].x
           || points[L+1].y == points[L].y && points[L+2].y == points[L].y
           || (points[L+1].y - points[L].y)*(points[L+2].x - points[L+1].x) == (points[L+2].y - points[L+1].y)*(points[L+1].x - points[L].x)) { // 垂直 or 水平
            return vector<Edge>{
                Edge{recursiveVoronoi(L, L+2)[0].start, recursiveVoronoi(L, L+2)[0].end, points[L], points[L+1]},
                Edge{recursiveVoronoi(L+1, R)[0].start, recursiveVoronoi(L+1, R)[0].end, points[L+1], points[L+2]}
            };
        }
        double slope1 = -(points[L+1].x - points[L].x) / (points[L+1].y - points[L].y);
        double slope2 = -(points[L+2].x - points[L+1].x) / (points[L+2].y - points[L+1].y);

        Point mid1 = calculateMidPoint(points[L], points[L+1]);
        Point mid2 = calculateMidPoint(points[L+1], points[L+2]);

        pair<Point, Point> bisector1 = calculatePerpendicularBisector(points[L], points[L+1]);
        pair<Point, Point> bisector2 = calculatePerpendicularBisector(points[L+1], points[L+2]);
        pair<Point, Point> bisector3 = calculatePerpendicularBisector(points[L], points[L+2]);

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
            Point t = calculateMidPoint(points[L], points[L+1]);
            if(isOnSegment(intersection, bisector1.first, t)) {
                bisector1.second = intersection;
            }
            else {
                bisector1.first = intersection;
            }

            t = calculateMidPoint(points[L+1], points[L+2]);
            if(isOnSegment(intersection, bisector2.first, t)) {
                bisector2.second = intersection;
            }
            else {
                bisector2.first = intersection;
            }

            Swap(bisector3.second, bisector3.first);
            if(AB2 + BC2 == CA2 || AB2 + CA2 == BC2 || BC2 + CA2 == AB2) { // 直角
                if(calculateslope(Edge{points[L+1], intersection}) > 0) {
                    if(points[L+1].x > intersection.x) {
                        bisector3.first = intersection;
                    } else {
                        bisector3.second = intersection;
                    }
                }
                else {
                    if(points[L+1].x > intersection.x) {
                        bisector3.second = intersection;
                    } else {
                        bisector3.first = intersection;
                    }
                }
            }
            else if(AB2 + BC2 < CA2 || AB2 + CA2 < BC2 || BC2 + CA2 < AB2) { // 鈍角
                if(calculateslope(Edge{calculateMidPoint(points[L], points[L+2]), intersection}) > 0) {
                    if(calculateMidPoint(points[L], points[L+2]).x > intersection.x) {
                        bisector3.first = intersection;
                    } else {
                        bisector3.second = intersection;
                    }
                }
                else {
                    if(calculateMidPoint(points[L], points[L+2]).x > intersection.x) {
                        bisector3.second = intersection;
                    } else {
                        bisector3.first = intersection;
                    }
                }
            }
            else { // 銳角
                if(calculateslope(Edge{calculateMidPoint(points[L], points[L+2]), intersection}) > 0) {
                    if(calculateMidPoint(points[L], points[L+2]).x > intersection.x) {
                        bisector3.second = intersection;
                    } else {
                        bisector3.first = intersection;
                    }
                }
                else {
                    if(calculateMidPoint(points[L], points[L+2]).x > intersection.x) {
                        bisector3.first = intersection;
                    } else {
                        bisector3.second = intersection;
                    }
                }

            }
        }
        else {
            sort(points.begin()+L, points.begin()+L+3, [](Point a, Point b) {
                    return a.y > b.y;
                 });
            return vector<Edge>{
                Edge{recursiveVoronoi(L, L+2)[0].start, recursiveVoronoi(L, L+2)[0].end, points[L+1], points[L]},
                Edge{recursiveVoronoi(L+1, R)[0].start, recursiveVoronoi(L+1, R)[0].end, points[L+2], points[L+1]}
            };
        }

        Swap(bisector1.first, bisector1.second);
        Swap(bisector2.first, bisector2.second);
        Swap(bisector3.first, bisector3.second);
        Point A = points[L], B = points[L+1], C = points[L+2];
        return vector<Edge>{
            Edge{bisector1.first, bisector1.second, (A.y < B.y ? A:B), (A.y > B.y ? A:B)},
            Edge{bisector2.first, bisector2.second, (B.y < C.y ? B:C), (B.y > C.y ? B:C)},
            Edge{bisector3.first, bisector3.second, (A.y < C.y ? A:C), (A.y > C.y ? A:C)}
        };
    }
    int mid = (L + R) / 2;
    vector<Edge> RightEdge = recursiveVoronoi(L, mid);
    vector<Edge> LeftEdge = recursiveVoronoi(mid, R);

    vector<Edge> RightConvexhull_Edge = constructConvexHull(RightEdge, 0);
    vector<Edge> LeftConvexhull_Edge = constructConvexHull(LeftEdge, 1);

    vector<Point> RightConvexhull = ConvexHull_Point(RightConvexhull_Edge);
    vector<Point> LeftConvexhull = ConvexHull_Point(LeftConvexhull_Edge);

    ofstream outfile("convexhull.txt");
    for (const auto& edge : RightConvexhull) {
        outfile << "R " << edge.x << " " << edge.y << "\n";
    }
    outfile << "\n";
    for (const auto& edge : LeftConvexhull) {
        outfile << "L " << edge.x << " " << edge.y << "\n";
    }
    outfile << "\n";
    outfile.close();

    vector<Edge> mid_edge, temp;


    pair<Point, Point> hyperplane = HyperPlane(LeftConvexhull, RightConvexhull);

    LeftConvexhull.erase(find(LeftConvexhull.begin(), LeftConvexhull.end(), hyperplane.first));
    RightConvexhull.erase(find(RightConvexhull.begin(), RightConvexhull.end(), hyperplane.second));
    double r_highest = 601, l_highest = 601;
    while(LeftConvexhull.size() && RightConvexhull.size()) {
        Point last_p = {0, 0};
        bool R_Edge = 0;
        int index = -1;
        pair<Point, Point> hp = calculatePerpendicularBisector(hyperplane.first, hyperplane.second);

        Swap(hp.first, hp.second);

        if(mid_edge.size()) {
            hp.second = prev(mid_edge.end()) -> start;
        }

        for(int i = 0; i < RightConvexhull_Edge.size(); i++) {
            if(RightConvexhull_Edge[i].B == hyperplane.second || RightConvexhull_Edge[i].A == hyperplane.second ||
               RightConvexhull_Edge[i].B == hyperplane.first || RightConvexhull_Edge[i].A == hyperplane.first) {
                Point intersect = IntersectionPoint(RightConvexhull_Edge[i], Edge{hp.first, hp.second});
                if(intersect.x < 0 || intersect.x > 600 || intersect.y < 0 || intersect.y > 600) {
                    //temp.push_back(RightConvexhull_Edge[i]);
                    //RightConvexhull_Edge.erase(RightConvexhull_Edge.begin() + i);
                    continue;
                }
                else {
                    if(intersect.y > last_p.y && isOnSegment(RightConvexhull_Edge[i].start, RightConvexhull_Edge[i].end, intersect)) {
                        last_p = intersect;
                        index = i;
                        R_Edge = 1;
                    }
                }
            }
        }
        for(int i = 0; i < LeftConvexhull_Edge.size(); i++) {
            if(LeftConvexhull_Edge[i].B == hyperplane.first || LeftConvexhull_Edge[i].A == hyperplane.first ||
               LeftConvexhull_Edge[i].B == hyperplane.second || LeftConvexhull_Edge[i].A == hyperplane.second) {
                Point intersect = IntersectionPoint(LeftConvexhull_Edge[i], Edge{hp.first, hp.second});
                if(intersect.x < 0 || intersect.x > 600 || intersect.y < 0 || intersect.y > 600) {
                    //temp.push_back(LeftConvexhull_Edge[i]);
                    //LeftConvexhull_Edge.erase(LeftConvexhull_Edge.begin() + i);
                    continue;
                }
                else {
                    if(intersect.y > last_p.y && isOnSegment(LeftConvexhull_Edge[i].start, LeftConvexhull_Edge[i].end, intersect)) {
                        last_p = intersect;
                        index = i;
                        R_Edge = 0;
                    }
                }
            }
        }
        if(index != -1 && R_Edge) {
            hp.first = last_p;
            Point t;
            if(crossProduct(RightConvexhull_Edge[index].end, hp.first, hp.second) < 0) {
                t = RightConvexhull_Edge[index].end;
                RightConvexhull_Edge[index].end = last_p;
            }
            else {
                t = RightConvexhull_Edge[index].start;
                RightConvexhull_Edge[index].start = last_p;
            }

            temp.push_back(RightConvexhull_Edge[index]);
            RightConvexhull_Edge.erase(RightConvexhull_Edge.begin() + index);

            for(int i = 0; i < RightConvexhull_Edge.size(); i++) {
                if((RightConvexhull_Edge[i].start == t || RightConvexhull_Edge[i].end == t) &&
                   (RightConvexhull_Edge[i].start.x > t.x || RightConvexhull_Edge[i].end.x > t.x)) {
                    RightConvexhull_Edge.erase(RightConvexhull_Edge.begin() + i);
                    break;
                }
            }
            for(int i = 0; i < RightConvexhull.size(); i++) {
                pair<Point, Point> p = calculatePerpendicularBisector(hyperplane.first, RightConvexhull[i]);
                Swap(p.first, p.second);
                if(isOnSegment(p.first, p.second, last_p) && RightConvexhull[i].y < r_highest) {
                    hyperplane.second = RightConvexhull[i];
                    r_highest = RightConvexhull[i].y;
                    RightConvexhull.erase(RightConvexhull.begin()+i);
                    break;
                }
            }
        }
        else if(index != -1) {
            hp.first = last_p;
            Point t;
            if(crossProduct(LeftConvexhull_Edge[index].start, hp.first, hp.second) > 0) {
                t = LeftConvexhull_Edge[index].start;
                LeftConvexhull_Edge[index].start = last_p;
            }
            else {
                t = LeftConvexhull_Edge[index].end;
                LeftConvexhull_Edge[index].end = last_p;
            }
            temp.push_back(LeftConvexhull_Edge[index]);
            LeftConvexhull_Edge.erase(LeftConvexhull_Edge.begin() + index);

            for(int i = 0; i < LeftConvexhull_Edge.size(); i++) {
                if((LeftConvexhull_Edge[i].start == t || LeftConvexhull_Edge[i].end == t) &&
                   (LeftConvexhull_Edge[i].start.x < t.x || LeftConvexhull_Edge[i].end.x < t.x)) {
                    LeftConvexhull_Edge.erase(LeftConvexhull_Edge.begin() + i);
                    break;
                }
            }
            for(int i = 0; i < LeftConvexhull.size(); i++) {
                pair<Point, Point> p = calculatePerpendicularBisector(LeftConvexhull[i], hyperplane.second);
                Swap(p.first, p.second);
                if(isOnSegment(p.first, p.second, last_p) && LeftConvexhull[i].y < l_highest) {
                    hyperplane.first = LeftConvexhull[i];
                    l_highest = LeftConvexhull[i].y;
                    LeftConvexhull.erase(LeftConvexhull.begin()+i);
                    break;
                }
            }
        }
        else {
            double L_dis = 600*600, R_dis = 600*600, L_index = -1, R_index = -1;
            for(int i = 0; i < RightConvexhull.size(); i++) {
                if(pow(hyperplane.first.x - RightConvexhull[i].x, 2) + pow(hyperplane.first.y - RightConvexhull[i].y, 2) < R_dis) {
                    R_dis = pow(hyperplane.first.x - RightConvexhull[i].x, 2) + pow(hyperplane.first.y - RightConvexhull[i].y, 2);
                    R_index = i;
                }
            }
            for(int i = 0; i < LeftConvexhull.size(); i++) {
                if(pow(hyperplane.second.x - LeftConvexhull[i].x, 2) + pow(hyperplane.second.y - LeftConvexhull[i].y, 2) < L_dis) {
                    L_dis = pow(hyperplane.second.x - LeftConvexhull[i].x, 2) + pow(hyperplane.second.y - LeftConvexhull[i].y, 2);
                    L_index = i;
                }
            }
            if(R_dis < L_dis) {
                hyperplane.second = RightConvexhull[R_index];
                RightConvexhull.erase(RightConvexhull.begin() + R_index);
            }
            else {
                hyperplane.first = LeftConvexhull[L_index];
                LeftConvexhull.erase(LeftConvexhull.begin() + L_index);
            }
        }
        mid_edge.push_back(Edge{hp.first, hp.second});
    }

    while(RightConvexhull.size()) {
        Point last_p = {0, 0};
        int index = -1;
        pair<Point, Point> hp = calculatePerpendicularBisector(hyperplane.first, hyperplane.second);
        Swap(hp.first, hp.second);
        if(mid_edge.size()) hp.second = prev(mid_edge.end()) -> start;

        for(int i = 0; i < RightConvexhull_Edge.size(); i++) {
            if(RightConvexhull_Edge[i].B == hyperplane.second || RightConvexhull_Edge[i].A == hyperplane.second ||
               RightConvexhull_Edge[i].B == hyperplane.first || RightConvexhull_Edge[i].A == hyperplane.first) {
                Point intersect = IntersectionPoint(RightConvexhull_Edge[i], Edge{hp.first, hp.second});
                if(intersect.x < 0 || intersect.x > 600 || intersect.y < 0 || intersect.y > 600) {
                    temp.push_back(RightConvexhull_Edge[i]);
                    RightConvexhull_Edge.erase(RightConvexhull_Edge.begin() + i);
                    continue;
                }
                else {
                    Point t = IntersectionPoint(RightConvexhull_Edge[i], Edge{hp.first, hp.second});
                    if(t.y > last_p.y && t.y < hp.second.y) {
                        last_p = t;
                        index = i;
                    }
                }
            }
        }
        if(index != -1) {
            hp.first = last_p;
            Point t;
            if(crossProduct(RightConvexhull_Edge[index].end, hp.first, hp.second) < 0) {
                t = RightConvexhull_Edge[index].end;
                RightConvexhull_Edge[index].end = last_p;
            }
            else {
                t = RightConvexhull_Edge[index].start;
                RightConvexhull_Edge[index].start = last_p;
            }

            temp.push_back(RightConvexhull_Edge[index]);
            RightConvexhull_Edge.erase(RightConvexhull_Edge.begin() + index);

            for(int i = 0; i < RightConvexhull_Edge.size(); i++) {
                if((RightConvexhull_Edge[i].start == t || RightConvexhull_Edge[i].end == t) &&
                   (RightConvexhull_Edge[i].start.x > t.x || RightConvexhull_Edge[i].end.x > t.x)) {
                    RightConvexhull_Edge.erase(RightConvexhull_Edge.begin() + i);
                    break;
                }
            }

            for(int i = 0; i < RightConvexhull.size(); i++) {
                pair<Point, Point> p = calculatePerpendicularBisector(hyperplane.first, RightConvexhull[i]);
                Swap(p.first, p.second);
                if(isOnSegment(p.first, p.second, last_p) && RightConvexhull[i].y < r_highest) {
                    hyperplane.second = RightConvexhull[i];
                    RightConvexhull.erase(RightConvexhull.begin()+i);
                    break;
                }
            }
        }
        else if(index != -1) {
            double R_dis = 600*600, R_index = -1;
            for(int i = 0; i < RightConvexhull.size(); i++) {
                if(pow(hyperplane.first.x - RightConvexhull[i].x, 2) + pow(hyperplane.first.y - RightConvexhull[i].y, 2) < R_dis &&
                    RightConvexhull[i].y < r_highest) {
                    R_dis = pow(hyperplane.first.x - RightConvexhull[i].x, 2) + pow(hyperplane.first.y - RightConvexhull[i].y, 2);
                    R_index = i;
                }
            }
            if(R_index != -1) {
                r_highest = RightConvexhull[R_index].y;
                hyperplane.second = RightConvexhull[R_index];
                RightConvexhull.erase(RightConvexhull.begin() + R_index);
            }
        }
        mid_edge.push_back(Edge{hp.first, hp.second});
        if(index == -1) break;
    }

    while(LeftConvexhull.size()) {
        Point last_p = {0, 0};
        int index = -1;
        pair<Point, Point> hp = calculatePerpendicularBisector(hyperplane.first, hyperplane.second);
        Swap(hp.first, hp.second);
        if(mid_edge.size()) hp.second = prev(mid_edge.end()) -> start;
        for(int i = 0; i < LeftConvexhull_Edge.size(); i++) {
            if(LeftConvexhull_Edge[i].B == hyperplane.first || LeftConvexhull_Edge[i].A == hyperplane.first) {
                Point intersect = IntersectionPoint(LeftConvexhull_Edge[i], Edge{hp.first, hp.second});
                if(intersect.x < 0 || intersect.x > 600 || intersect.y < 0 || intersect.y > 600) {
                    temp.push_back(LeftConvexhull_Edge[i]);
                    LeftConvexhull_Edge.erase(LeftConvexhull_Edge.begin() + i);
                    continue;
                }
                else {
                    Point t = IntersectionPoint(LeftConvexhull_Edge[i], Edge{hp.first, hp.second});
                    if(t.y > last_p.y && t.y < hp.second.y) {
                        last_p = t;
                        index = i;
                    }
                }
            }
        }
        if(index != -1) {
            hp.first = last_p;
            Point t;
            if(crossProduct(LeftConvexhull_Edge[index].start, hp.first, hp.second) > 0) {
                t = LeftConvexhull_Edge[index].start;
                LeftConvexhull_Edge[index].start = last_p;
            }
            else {
                t = LeftConvexhull_Edge[index].end;
                LeftConvexhull_Edge[index].end = last_p;
            }
            temp.push_back(LeftConvexhull_Edge[index]);
            LeftConvexhull_Edge.erase(LeftConvexhull_Edge.begin() + index);

            for(int i = 0; i < LeftConvexhull_Edge.size(); i++) {
                if((LeftConvexhull_Edge[i].start == t || LeftConvexhull_Edge[i].end == t) &&
                   (LeftConvexhull_Edge[i].start.x < t.x || LeftConvexhull_Edge[i].end.x < t.x)) {
                    LeftConvexhull_Edge.erase(LeftConvexhull_Edge.begin() + i);
                    break;
                }
            }
            for(int i = 0; i < LeftConvexhull.size(); i++) {
                pair<Point, Point> p = calculatePerpendicularBisector(LeftConvexhull[i], hyperplane.second);
                Swap(p.first, p.second);
                if(isOnSegment(p.first, p.second, last_p) && LeftConvexhull[i].y < l_highest) {
                    hyperplane.first = LeftConvexhull[i];
                    LeftConvexhull.erase(LeftConvexhull.begin()+i);
                    break;
                }
            }
        }
        else if(index != -1) {
            double L_dis = 600*600, L_index = -1;
            for(int i = 0; i < LeftConvexhull.size(); i++) {
                if(pow(hyperplane.second.x - LeftConvexhull[i].x, 2) + pow(hyperplane.second.y - LeftConvexhull[i].y, 2) < L_dis &&
                   LeftConvexhull[i].y < l_highest) {
                    L_dis = pow(hyperplane.second.x - LeftConvexhull[i].x, 2) + pow(hyperplane.second.y - LeftConvexhull[i].y, 2);
                    L_index = i;
                }
            }
            if(L_index != -1) {
                l_highest = LeftConvexhull[L_index].y;
                hyperplane.first = LeftConvexhull[L_index];
                LeftConvexhull.erase(LeftConvexhull.begin() + L_index);
            }
        }
        mid_edge.push_back(Edge{hp.first, hp.second});
        if(index == -1) break;
    }

    //cout << LeftConvexhull.size() << ' ' << RightConvexhull.size() << endl;

    if(!LeftConvexhull.size() && !RightConvexhull.size()) {
        Point last_p = {0, 0};
        pair<Point, Point> hp = calculatePerpendicularBisector(hyperplane.first, hyperplane.second);
        Swap(hp.first, hp.second);
        if(mid_edge.size()) hp.second = prev(mid_edge.end()) -> start;
        mid_edge.push_back(Edge{hp.first, hp.second});
    }

    outfile.open("hyperplane.txt");
    for (const auto& edge : mid_edge) {
        outfile << edge.end.x << " " << edge.end.y << "\n";
    }
    outfile << mid_edge[mid_edge.size()-1].start.x << " " << mid_edge[mid_edge.size()-1].start.y << "\n";
    outfile.close();

    voronoi.insert(voronoi.end(), RightConvexhull_Edge.begin(), RightConvexhull_Edge.end());
    voronoi.insert(voronoi.end(), LeftConvexhull_Edge.begin(), LeftConvexhull_Edge.end());
    voronoi.insert(voronoi.end(), temp.begin(), temp.end());
    voronoi.insert(voronoi.end(), mid_edge.begin(), mid_edge.end());
    return voronoi;
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
        //cout << edge.start.x << ' ' << edge.start.y << ' ' << edge.end.x << ' ' << edge.end.y << endl;
    }
    outfile << "\n";
}

int main() {
    read_points("points.txt");
    sort(points.begin(), points.end());
    edges = recursiveVoronoi(0, points.size());
    write_edges("lines.txt");
    return 0;
}
