#include <stdio.h>
#include <iostream>
#include <math.h>
#include <algorithm>
#include <iomanip>
#include <vector>

const double EPS = 1E-10;

class Point;
class Line;
class Segment;
class PointInt;

class Point {
public:
    Point() : _x(0), _y(0) {};
    Point(const double & x, const double & y) : _x(x), _y(y) {};
    Point(const Point & s) : _x(s.x()), _y(s.y()) {};

    double x() const
    {
        return _x;
    }

    double y() const
    {
        return _y;
    }

    void setX(double newX)
    {
        _x = newX;
    }

    void setY(double newY)
    {
        _y = newY;
    }

    double len()
    {
        return sqrt(_x * _x + _y * _y);
    }

    Point operator+(const Point & p)
    {
        return Point(_x + p.x(), _y + p.y());
    }

    Point operator-(const Point & p)
    {
        return Point(_x - p.x(), _y - p.y());
    }

    double operator^(const Point & p)
    {
        return _x * p.y() - _y * p.x();
    }

    double operator*(const Point & p)
    {
        return _x * p.x() + _y * p.y();
    }
    
    Point operator*(double k)
    {
        return Point(k * this->_x, k * this->_y);
    }

    bool operator==(const Point & p)
    {
        return (fabs(_x - p.x()) <= EPS) && (fabs(_y - p.y()) <= EPS);
    }

    bool operator!=(const Point & p)
    {
        return !((*this) == p);
    }

    friend std::istream& operator>>(std::istream& is, Point& p)
    {
        double x, y;
        is >> x >> y;
        p.setX(x);
        p.setY(y);
        return is;
    }

    bool cmp(const Point & first, const Point & second)
    {

    }
    
private:
    double _x;
    double _y;
};

class PointInt {
public:
    PointInt() : _x(0), _y(0) {};
    PointInt(const int & x, const int & y) : _x(x), _y(y) {};
    PointInt(const PointInt & s) : _x(s.x()), _y(s.y()) {};

    int x() const
    {
        return _x;
    }

    int y() const
    {
        return _y;
    }

    void setX(int newX)
    {
        _x = newX;
    }

    void setY(int newY)
    {
        _y = newY;
    }

    long long len2()
    {
        return 1LL * _x * _x + _y * _y;
    }

    PointInt operator-(PointInt & s)
    {
        return PointInt(_x - s.x(), _y - s.y());
    }

    long long operator^(const PointInt & p)
    {
        return 1LL * _x * p.y() - 1LL * _y * p.x();
    }

    bool operator==(const PointInt & p)
    {
        return (_x == p.x() && _y == p.y());
    }

    friend std::istream& operator>>(std::istream& is, PointInt& p)
    {
        int x, y;
        is >> x >> y;
        p.setX(x);
        p.setY(y);
        return is;
    }

    bool operator<(PointInt& s)
    {
        if (_x < s.x())
            return true;
        else
            if (_x == s.x())
                if (_y < s.y())
                    return true;
                else
                    return false;
            else
                return false;
    }

private:
    int _x;
    int _y;
};

class Line{
public:
    struct LineParams {
        LineParams(double a, double b, double c) : a(a), b(b), c(c) {};
        double a;
        double b;
        double c;
    };

    Line(const Point & first, const Point & second)
    {
        a = first.y() - second.y();
        b = - first.x() + second.x();
        c = second.y() * first.x() - second.x() * first.y();
    }
    
    Line(double a, double b, double c) : a(a), b(b), c(c) {};

    double distanceToPoint(const Point & s)
    {
        return (a * s.x() + b * s.y() + c) / sqrt(a * a + b * b);
    }

    LineParams getParams() const
    {
        return LineParams(a, b, c);
    }
    
    static bool parallel(const Line & first, const Line & second)
    {
        Line::LineParams param1 = first.getParams();
        Line::LineParams param2 = second.getParams();
        //solving using Kramor's method
        double d = param1.a * param2.b - param1.b * param2.a;
        if (fabs(d) <= EPS)
            return true;
        return false;
    }

    static bool same(const Line & first, const Line & second)
    {
        if (!parallel(first, second))
            return false;
        Line::LineParams param1 = first.getParams();
        Line::LineParams param2 = second.getParams();
        //solving using Kramor's method
        double c = param1.a * param2.c - param1.c * param2.a;
        if (fabs(c) <= EPS)
            return true;
        return false;
    }

    static Point intersect(const Line & first, const Line & second)
    {
        Line::LineParams param1 = first.getParams();
        Line::LineParams param2 = second.getParams();
        //solving using Kramor's method
        double d = param1.a * param2.b - param1.b * param2.a;
        if (fabs(d) <= EPS) {
            std::cout << "NO" << std::endl;
            return Point();
        } else {
            double x = (-param1.c * param2.b + param2.c * param1.b) / d;
            double y = (param1.a * (-param2.c) + param1.c * param2.a) / d;
            return Point(x, y);
        }
    }
private:
    double a, b, c;
};

class Segment{
public:
    Segment(const Point & first, const Point & second)
    {
        begin = first;
        end = second;
    }

    Point getBegin() const
    {
        return begin;
    }

    Point getEnd() const
    {
        return end;
    }

    double len()
    {
        return (end - begin).len();
    }

    bool hasPoint(const Point & s)
    {
        if (begin == s || end == s)
            return true;
        if (begin == end)
            return false;
        Line line(begin, end);
        if (fabs(line.distanceToPoint(s)) > EPS)
            return false;
        if ((begin - s) * (end - s) <= EPS)
            return true;
        else
            return false;
    }

    double distanceToPoint(const Point & s)
    {
        if (begin == end)
            return (begin - s).len();
        Line l(begin, end);
        double dist = l.distanceToPoint(s);
        Point normal = Point(l.getParams().a, l.getParams().b);
        Point n = Point(s) - (normal * (dist / normal.len()));
        if (hasPoint(n))
            return fabs(dist);
        else
            return std::min((begin - s).len(), (end - s).len());
    }

    static bool intersect(const Segment & first, const Segment & second)
    {
        Segment ab = first;
        Segment cd = second;
        if (ab.len() <= EPS) 
            return cd.hasPoint(first.getBegin());
        if (cd.len() <= EPS)
            return ab.hasPoint(second.getBegin());

        Line a(first.getBegin(), first.getEnd());
        Line b(second.getBegin(), second.getEnd());
        if (Line::parallel(a, b)) {
            if (Line::same(a, b))
                if (ab.hasPoint(second.getBegin()) || 
                     ab.hasPoint(second.getEnd()) || 
                     cd.hasPoint(first.getBegin()) || 
                     cd.hasPoint(first.getEnd()))
                    return true;
                else
                    return false;
            else 
                return false;
        }
        Point intersect = Line::intersect(a, b);
        if (ab.hasPoint(intersect) && cd.hasPoint(intersect))
            return true;
        return false;
    }

private:
    Point begin;
    Point end;
};

void pointToSegment()
{
    int x, y;
    std::cin >> x >> y;
    Point a(x, y);
    std::cin >> x >> y;
    Point b(x, y);
    std::cin >> x >> y;
    Point c(x, y);
    Segment bc = Segment(b, c);
    std::cout << std::fixed << std::setprecision(10);
    std::cout << bc.distanceToPoint(a) << std::endl;
}

void twoLines()
{
    double x, y, z;
    std::cin >> x >> y >> z;
    Line a(x, y, z);
    std::cin >> x >> y >> z;
    Line b(x, y, z);
    Point intersect = Line::intersect(a, b);
    std::cout << std::fixed << std::setprecision(10);
    std::cout << intersect.x() << ' ' << intersect.y() << std::endl;
}

void twoSegments()
{
    int koord[8];
    for (int i = 0; i < 8; ++i)
        std::cin >> koord[i];
    Segment first(Point(koord[0], koord[1]), Point(koord[2], koord[3]));
    Segment second(Point(koord[4], koord[5]), Point(koord[6], koord[7]));
    if (Segment::intersect(first, second))
        std::cout << "YES" << std::endl;
    else
        std::cout << "NO" << std::endl;
}

long long sq(int a)
{
    return 1LL * a * a;
}

void circleIntersect(PointInt c, int r, Line l, int k)
{
    Point centre(c.x(), c.y());
    double dist = l.distanceToPoint(centre);
    Point normal(l.getParams().a, l.getParams().b);
    Point h = centre - (normal * (dist / normal.len()));
    std::cout << std::fixed << std::setprecision(10);
    std::cout << h.x() << ' ' << h.y() << std::endl;
    if (k != 1) {
        double move = sqrt(double(r) * r - dist * dist);
        Point dir(-l.getParams().b, l.getParams().a);
        Point a = h + dir * (move / dir.len());
        Point b = h - dir * (move / dir.len());

        std::cout << fabs(dist) << ' ' << move << std::endl;
        std::cout << a.x() << ' ' << a.y() << std::endl;
        std::cout << b.x() << ' ' << b.y() << std::endl;
    }
}

void twocircles()
{
    int n;
    std::cin >> n;
    for (int i = 0; i < n; ++i) {
        PointInt c1;
        int r1, r2;
        PointInt c2;
        std::cin >> c1 >> r1 >> c2 >> r2;
        if (c1 == c2)
            if (r1 == r2) {
                std::cout << 3 << std::endl;
                continue;
            } else {
                std::cout << 0 << std::endl;
                continue;
        }
        long long distC = (c2 - c1).len2();
        if (distC > sq(r1 + r2) || distC < sq(r1 - r2)) {
            std::cout << 0 << std::endl;
            continue;
        }
        int k = 0;
        if (distC == sq(r1 + r2) || distC == sq(r1 - r2))
            k = 1;
        else
            k = 2;
        std::cout << k << std::endl;
        long long a, b, c;
        a = 2 * (c2.x() - c1.x());
        b = 2 * (c2.y() - c1.y());
        c = (c1.x() - c2.x()) * (c1.x() + c2.x()) +
             (c1.y() - c2.y()) * (c1.y() + c2.y()) +
             (r2 - r1) * (r2 + r1);
        Line l(a, b, c);
        circleIntersect(c1, r1, l, k);
    }
}

void circleAndPoint()
{
    PointInt c, s;
    int r;
    std::cout << std::fixed << std::setprecision(10);
    std::cin >> c >> r >> s;
    if ((s - c).len2() < sq(r)) {
        std::cout << 0 << std::endl;
    } else {
        if ((s - c).len2() == sq(r)) {
            std::cout << 1 << std::endl;
            std::cout << s.x() << ' ' << s.y() << std::endl;
        } else { 
            std::cout << 2 << std::endl;

            Point o(s.x(), s.y());
            Point centre(c.x(), c.y());
            double l = (o - centre).len();
            double x = 1.0 * sq(r) / l;
            Point op = (centre - o) * ((l - x) / l);
            Point normal(-op.y(), op.x());
            double h = sqrt(double(r) * r - x * x);
            normal = normal * (h / normal.len());

            std::cout << (op + o).x() << ' ' << (op + o).y() << std::endl;
            std::cout << l - x << ' ' << h << std::endl;
            Point a = op + o + normal;
            Point b = op + o - normal;
            std::cout << a.x() << ' ' << a.y() << std::endl;
            std::cout << b.x() << ' ' << b.y() << std::endl;
        }
    }
}

void polygonArea()
{
    int n;
    std::cin >> n;
    double square = 0;
    int x, y;
    std::cin >> x >> y;
    Point p(x, y), q;
    Point first(x, y);
    for (int i = 1; i < n; ++i) {
        std::cin >> x >> y;
        q = Point(x, y);
        square += p ^ q;
        p = q;
    }
    square += p ^ first;
    square = 0.5 * fabs(square);
    std::cout << std::fixed << std::setprecision(6);
    std::cout << square << std::endl;
}

int signf(double f)
{
    if (f > 0)
        return 1;
    else
        if (f == 0)
            return 0;
        else
            return -1;
}

void polygonConvex()
{
    int n;
    std::cin >> n;
    double square = 0;
    Point p, q, r, first, second;
    std::cin >> first >> second >> r;
    p = first;
    q = second;
    int sign = 0;
    int k = 3;
    while (((q - p) ^ (r - q)) == 0) {
        p = q;
        q = r;
        std::cin >> r;
        ++k;
    }
    sign = signf((q - p) ^ (r - q));
    for (int i = k; i < n; ++i) {
        std::cin >> r;
        if (signf((q - p) ^ (r - q)) != sign) {
            std::cout << "NO" << std::endl;
            return;
        }
        p = q;
        q = r;
    }
    if (k == n) {
        p = q;
        q = r;
    }
    r = first;
    if (signf((q - p) ^ (r - q)) != sign) {
        std::cout << "NO" << std::endl;
        return;
    }
    p = q;
    q = r;
    r = second;
    if (signf((q - p) ^ (r - q)) != sign) {
        std::cout << "NO" << std::endl;
        return;
    }
    std::cout << "YES" << std::endl;
}

void hull()
{
    int n;
    std::cin >> n;
    std::vector<PointInt> dots(n);
    int mini = 0;
    for (int i = 0; i < n; ++i) {
        std::cin >> dots[i];
        if (dots[i] < dots[mini])
            mini = i;
    }
    PointInt minpoint = dots[mini];
    std::swap(dots[mini], dots[0]);
    std::sort(dots.begin() + 1, dots.end(), [&](PointInt a, PointInt b){
        long long x = (a - minpoint) ^ (b - minpoint);
        if (x > 0)
            return true;
        else
            if (x == 0)
                return (a - minpoint).len2() < (b - minpoint).len2();
            else
                return false;
    });
    std::vector<PointInt> ans;
    int k = 0;
    for (int i = 0; i < n; ++i) {
        if (i < n - 1 && dots[i] == dots[i + 1])
            continue;
        while (k > 1 && ((ans[k - 1] - ans[k - 2]) ^ (dots[i] - ans[k - 1])) <= 0) {
            ans.pop_back();
            --k;
        }
        ans.push_back(dots[i]);
        ++k;
    }
    std::cout << ans.size() << std::endl;
    for (int i = 0; i < ans.size(); ++i)
        std::cout << ans[i].x() << ' ' << ans[i].y() << std::endl;
    //std::cin >> minpoint;
}

int main() {
    freopen("hull.in", "r", stdin);
    freopen("hull.out", "w", stdout);
    //pointToSegment();
    //twoLines();
    //twoSegments();
    //twocircles();
    //circleAndPoint();
    //polygonArea();
    //polygonConvex();
    hull();
    return 0;
}