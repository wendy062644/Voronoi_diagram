#include <bits/stdc++.h>
using namespace std;

struct Point {
    double x, y;

    bool operator==(const Point& other) const {
        return x == other.x && y == other.y;
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
        return x != other.x && y != other.y;
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
        cout << t.x << ' ' << t.y << endl;
        return vector<Edge>{
            Edge{bisector1.first, bisector1.second, A, B}
            //Edge{bisector2.first, bisector2.second, b, c},
            //Edge{bisector3.first, bisector3.second, a, c}
        };
    }
    return vector<Edge>{};
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
            voronoiEdges.erase(voronoiEdges.begin()+i--);
            continue;
        }
        double m = (voronoiEdges[i].start.y - voronoiEdges[i].end.y)/(voronoiEdges[i].start.x - voronoiEdges[i].end.x);
        if(Left) {
            if(m > 0 && (voronoiEdges[i].start.y == 0 || voronoiEdges[i].start.x == 0)) {
                Convexhull.push_back(voronoiEdges[i]);
                voronoiEdges.erase(voronoiEdges.begin()+i--);
            }
            else if(voronoiEdges[i].end.y == 0 || voronoiEdges[i].end.x == 0) {
                Convexhull.push_back(voronoiEdges[i]);
                voronoiEdges.erase(voronoiEdges.begin()+i--);
            }
        }
        else if(!Left) {
            if(m > 0 && voronoiEdges[i].end.y == 600 || voronoiEdges[i].end.x == 600) {
                Convexhull.push_back(voronoiEdges[i]);
                voronoiEdges.erase(voronoiEdges.begin()+i--);
            }
            else if(voronoiEdges[i].start.y == 600 || voronoiEdges[i].start.x == 600) {
                Convexhull.push_back(voronoiEdges[i]);
                voronoiEdges.erase(voronoiEdges.begin()+i--);
            }
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

vector<Edge> recursiveVoronoi(int L, int R) {
    cout << L << ' ' << R << endl;
    vector<Edge> voronoi;
    vector<Point> Convex_Hull = convexHull(&points[L], &points[R - 1]), Points;
    if(Convex_Hull.size() == R - L - 1) {
        Point p;
        int len = Convex_Hull.size();
        for(int i = L; i < R; i++) {
            for (int j = 0; j < len; j++) {
                if (points[i] == Convex_Hull[j])
                    break;
                if (j == len - 1) {
                    p = points[i];
                }
            }
        }

        for(int i = 0; i < len; i++) {
            Point A = Convex_Hull[i], B = Convex_Hull[(i + 1) % len], t;
            vector<Edge> edge = trangle(A, B, p, t);
            //cout << edge[0].start.x << ' ' << edge[0].start.y << ' ' << edge[0].end.x << ' ' << edge[0].end.y << endl;
            Points.push_back(t);
            voronoi.insert(voronoi.end(), edge.begin(), edge.end());
        }
        for(int i = 0; i < Points.size(); i++) {
            pair<Point, Point> temp = adjustPointToBoundary(Points[i], Points[(i+1)%int(Points.size())]);
            //cout << temp.first.x << ' ' << temp.first.y << ' ' << temp.second.x << ' ' << temp.second.y  << endl;
            voronoi.push_back({temp.first, temp.second , p, Convex_Hull[i]});
        }

        return voronoi;
    }
    else if (R - L == 2) {
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
            return vector<Edge>{
                Edge{recursiveVoronoi(L, L+2)[0].start, recursiveVoronoi(L, L+2)[0].end, points[L], points[L+1]},
                Edge{recursiveVoronoi(L+1, R)[0].start, recursiveVoronoi(L+1, R)[0].end, points[L+1], points[L+2]}
            };
        }

        Swap(bisector1.first, bisector1.second);
        Swap(bisector2.first, bisector2.second);
        Swap(bisector3.first, bisector3.second);
        Point A = points[L], B = points[L+1], C = points[L+2];
        Swap(B, C);
        Swap(A, C);
        Swap(A, B);
        return vector<Edge>{
            Edge{bisector1.first, bisector1.second, A, B},
            Edge{bisector2.first, bisector2.second, B, C},
            Edge{bisector3.first, bisector3.second, A, C}
        };
    }
    int mid = (L + R) / 2, l_index = 0, r_index = 0;
    vector<Edge> RightEdge = recursiveVoronoi(L, mid);
    vector<Edge> LeftEdge = recursiveVoronoi(mid, R);

    vector<Edge> RightConvexhull_Edge = constructConvexHull(RightEdge, 0);
    vector<Edge> LeftConvexhull_Edge = constructConvexHull(LeftEdge, 1);

    vector<Point> RightConvexhull = ConvexHull_Point(RightConvexhull_Edge);
    vector<Point> LeftConvexhull = ConvexHull_Point(LeftConvexhull_Edge);

    voronoi.insert(voronoi.end(), RightEdge.begin(), RightEdge.end());
    voronoi.insert(voronoi.end(), LeftEdge.begin(), LeftEdge.end());

    //cout << RightConvexhull_Edge[0].B.x << ' ' << RightConvexhull_Edge[0].B.y << endl;
    //cout << LeftConvexhull_Edge[0].B.x << ' ' << LeftConvexhull_Edge[0].B.y << endl;
    //cout << RightConvexhull_Edge.size() << ' ' << LeftConvexhull_Edge.size() << endl;

    while(l_index < LeftConvexhull_Edge.size() && r_index < RightConvexhull_Edge.size()) {
        pair<Point, Point> t = calculatePerpendicularBisector(LeftConvexhull_Edge[l_index].B, RightConvexhull_Edge[r_index].B);
        Swap(t.first, t.second);
        Point A = LeftConvexhull_Edge[l_index].B, B = RightConvexhull_Edge[r_index].B;
        Swap(A, B);
        Edge edge = Edge{t.first, t.second, A, B};
        if(LeftConvexhull_Edge[l_index].B.x == RightConvexhull_Edge[r_index].B.x
            || LeftConvexhull_Edge[l_index].B.y == RightConvexhull_Edge[r_index].B.y) {
            voronoi.push_back(LeftConvexhull_Edge[l_index]);
            voronoi.push_back(RightConvexhull_Edge[r_index]);
            voronoi.push_back(edge);
            l_index++;
            r_index++;
            continue;
        }
        auto L = doLinesIntersect(edge, LeftConvexhull_Edge[l_index]);
        auto R = doLinesIntersect(edge, RightConvexhull_Edge[r_index]);
        Point M = calculateMidPoint(LeftConvexhull_Edge[l_index].B, RightConvexhull_Edge[r_index].B);
        if(pow(M.x - L -> x, 2) + pow(M.y - L -> y, 2) <= pow(M.x - R -> x, 2) + pow(M.y - R -> y, 2)) {
            edge.start = *L;
            double m = (LeftConvexhull_Edge[l_index].start.y - LeftConvexhull_Edge[l_index].end.y)/(LeftConvexhull_Edge[l_index].start.x - LeftConvexhull_Edge[l_index].end.x);
            if(m > 0)
                LeftConvexhull_Edge[l_index].start = *L;
            else if(m < 0)
                LeftConvexhull_Edge[l_index].end = *L;

            voronoi.push_back(LeftConvexhull_Edge[l_index]);
            voronoi.push_back(edge);
            l_index++;
        }
        else {
            edge.start = *R;
            double m = (RightConvexhull_Edge[r_index].start.y - RightConvexhull_Edge[r_index].end.y)/(RightConvexhull_Edge[r_index].start.x - RightConvexhull_Edge[r_index].end.x);
            if(m > 0)
                RightConvexhull_Edge[r_index].end = *R;
            else if(m < 0)
                RightConvexhull_Edge[r_index].start = *R;
            voronoi.push_back(RightConvexhull_Edge[r_index]);
            voronoi.push_back(edge);
            r_index++;
        }
        /*if(l_index+1 == LeftConvexhull_Edge.size() && r_index+1 == RightConvexhull_Edge.size()) {
            if(LeftConvexhull_Edge[l_index].A.y > RightConvexhull_Edge[r_index].A.y) {
                l_index++;
            }
            else {
                r_index++;
            }
        }
        else if(l_index+1 == LeftConvexhull_Edge.size()) {
            if(LeftConvexhull_Edge[l_index].A.y > RightConvexhull_Edge[r_index+1].B.y) {
                l_index++;
            }
            else {
                r_index++;
            }
        }
        else if(r_index+1 == RightConvexhull_Edge.size()){
            if(LeftConvexhull_Edge[l_index+1].B.y > RightConvexhull_Edge[r_index].A.y) {
                l_index++;
            }
            else {
                r_index++;
            }
        }
        else {
            if(LeftConvexhull_Edge[l_index+1].B.y > RightConvexhull_Edge[r_index+1].B.y) {
                l_index++;
            }
            else {
                r_index++;
            }
        }*/
    }
    //cout << l_index << ' ' << r_index << endl;
    while(l_index < LeftConvexhull_Edge.size()) {
        pair<Point, Point> t = calculatePerpendicularBisector(LeftConvexhull_Edge[l_index].B, RightConvexhull_Edge[r_index-1].A);
        Swap(t.first, t.second);
        t.second = prev(voronoi.end()) -> start;
        Point A = LeftConvexhull_Edge[l_index].B, B = RightConvexhull_Edge[r_index-1].A;
        Swap(A, B);
        Edge edge = Edge{t.first, t.second, A, B};
        auto L = doLinesIntersect(edge, LeftConvexhull_Edge[l_index]);
        Point M = calculateMidPoint(LeftConvexhull_Edge[l_index].B, RightConvexhull_Edge[r_index-1].A);
        edge.start = *L;
        double m = (LeftConvexhull_Edge[l_index].start.y - LeftConvexhull_Edge[l_index].end.y)/(LeftConvexhull_Edge[l_index].start.x - LeftConvexhull_Edge[l_index].end.x);
        if(m > 0)
            LeftConvexhull_Edge[l_index].start = *L;
        else if(m < 0)
            LeftConvexhull_Edge[l_index].end = *L;
        voronoi.push_back(LeftConvexhull_Edge[l_index]);
        voronoi.push_back(edge);
        l_index++;
    }
    while(r_index < RightConvexhull_Edge.size()) {
        pair<Point, Point> t = calculatePerpendicularBisector(LeftConvexhull_Edge[l_index-1].A, RightConvexhull_Edge[r_index].B);
        Swap(t.first, t.second);
        t.second = prev(voronoi.end()) -> start;
        Point A = LeftConvexhull_Edge[l_index-1].A, B = RightConvexhull_Edge[r_index].B;
        Swap(A, B);
        Edge edge = Edge{t.first, t.second, A, B};
        auto R = doLinesIntersect(edge, RightConvexhull_Edge[r_index]);
        Point M = calculateMidPoint(LeftConvexhull_Edge[l_index-1].A, RightConvexhull_Edge[r_index].B);

        edge.start = *R;
        double m = (RightConvexhull_Edge[r_index].start.y - RightConvexhull_Edge[r_index].end.y)/(RightConvexhull_Edge[r_index].start.x - RightConvexhull_Edge[r_index].end.x);
        if(m > 0)
            RightConvexhull_Edge[r_index].end = *R;
        else if(m < 0)
            RightConvexhull_Edge[r_index].start = *R;
        voronoi.push_back(RightConvexhull_Edge[r_index]);
        voronoi.push_back(edge);
        r_index++;
    }
    pair<Point, Point> t = calculatePerpendicularBisector(LeftConvexhull_Edge[l_index-1].A, RightConvexhull_Edge[r_index-1].A);
    Swap(t.first, t.second);
    t.second = prev(voronoi.end()) -> start;
    Point A = LeftConvexhull_Edge[l_index-1].A, B = RightConvexhull_Edge[r_index-1].A;
    Edge edge = Edge{t.first, t.second, A, B};
    voronoi.push_back(edge);
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
}

int main() {
    read_points("points.txt");
    sort(points.begin(), points.end());
    edges = recursiveVoronoi(0, points.size());

    write_edges("lines.txt");
    return 0;
}
