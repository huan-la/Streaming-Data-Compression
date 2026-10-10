#ifndef ALGEBRAIC_CONVEX_HPP
#define ALGEBRAIC_CONVEX_HPP

#include <vector>
#include "algebraic/function.hpp"

class UpperHull {
    private:
        std::vector<Point2D> points;

    public:
        void concat(UpperHull& l_cvx, int offset=0) {
            for (int i=0; i<l_cvx.size(); i++) {
                this->points.push_back(
                    Point2D(l_cvx.at(i).x+offset, l_cvx.at(i).y)
                );
            }
        }

        void append(Point2D point, bool triangle_check = true) {
            this->points.push_back(point);
            if (!triangle_check) return;

            while (this->points.size() > 2) {
                int size = this->points.size();
                Point2D& p1 = this->points[size-1];
                Point2D& p2 = this->points[size-2];
                Point2D& p3 = this->points[size-3];

                double cross = (p2.x - p3.x)*(p1.y - p3.y) - (p2.y - p3.y)*(p1.x - p3.x);
                if (cross >= -EPS) { 
                    this->points[size-2] = points[size-1];
                    this->points.pop_back();       
                } else {
                    break;
                }
            }
        }

        void append_backward(Point2D point, bool triangle_check = true) {
            this->points.insert(this->points.begin(), point);
            if (!triangle_check) return;

            while (this->points.size() > 2) {
                Point2D& p1 = this->points[0];
                Point2D& p2 = this->points[1];
                Point2D& p3 = this->points[2];

                double cross = (p2.x - p3.x)*(p1.y - p3.y) - (p2.y - p3.y)*(p1.x - p3.x);
                if (cross >= -EPS) {
                    this->points.erase(this->points.begin()+1);
                }
                else {
                    break;
                }
            }
        }

        Point2D& at(int i) {
            return this->points[i];
        }

        int size() {
            return this->points.size();
        }

        void clear() {
            this->points.clear();
        }

        Point2D& pop() {
            Point2D& p = this->points.back();
            this->points.pop_back();

            return p;
        }

        void erase_from_begin(int length) {
            this->points.erase(this->points.begin(), this->points.begin() + length);
        }

        void erase_from_end(int index) {
            this->points.erase(this->points.begin() + index, this->points.end());
        }
};

class LowerHull {
    private:
        std::vector<Point2D> points;

    public:
        void concat(LowerHull& u_cvx, int offset=0) {
            for (int i=0; i<u_cvx.size(); i++) {
                this->points.push_back(
                    Point2D(u_cvx.at(i).x+offset, u_cvx.at(i).y)
                );
            }
        }
            
        void append(Point2D point, bool triangle_check = true) {
            this->points.push_back(point);
            if (!triangle_check) return;

            while (this->points.size() > 2) {
                int size = this->points.size();
                Point2D& p1 = this->points[size-1];
                Point2D& p2 = this->points[size-2];
                Point2D& p3 = this->points[size-3];

                double cross = (p2.x - p3.x)*(p1.y - p3.y) - (p2.y - p3.y)*(p1.x - p3.x);
                if (cross <= EPS) { 
                    this->points[size-2] = points[size-1]; 
                    this->points.pop_back();         
                } else {
                    break;
                }
            }
        }

        void append_backward(Point2D point, bool triangle_check = true) {
            this->points.insert(this->points.begin(), point);
            if (!triangle_check) return;

            while (this->points.size() > 2) {
                Point2D& p1 = this->points[0];
                Point2D& p2 = this->points[1];
                Point2D& p3 = this->points[2];

                double cross = (p2.x - p3.x)*(p1.y - p3.y) - (p2.y - p3.y)*(p1.x - p3.x);
                if (cross <= EPS) {
                    this->points.erase(this->points.begin()+1);
                }
                else {
                    break;
                }
            }
        }

        Point2D& at(int i) {
            return this->points[i];
        }

        int size() {
            return this->points.size();
        }

        void clear() {
            this->points.clear();
        }

        Point2D& pop() {
            Point2D& p = this->points.back();
            this->points.pop_back();

            return p;
        }

        void erase_from_begin(int length) {
            this->points.erase(this->points.begin(), this->points.begin() + length);
        }

        void erase_from_end(int index) {
            this->points.erase(this->points.begin() + index, this->points.end());
        }
};

class ConvexHull {
    public:
        LowerHull lower;
        UpperHull upper;

    public:
        void append(Point2D point) {
            this->lower.append(point);
            this->upper.append(point);
        }

        Point2D& at(int i) {
            if (i < this->lower.size()) {
                return this->lower.at(i);
            }
            else {
                return this->upper.at(i-this->lower.size()+1);
            }
        }

        Line side(int i) {
            Point2D p1 = this->at_ccw(i);
            Point2D p2 = this->at_ccw(i + 1);

            return Line::line(p1, p2);
        }

        Point2D& at_cw(int i) {
            i = i % this->size();
            if (i < this->upper.size()) {
                return this->upper.at(i);
            }
            else {
                return this->lower.at(this->size()-i);
            }
        }

        Point2D& at_ccw(int i) {
            i = i % this->size();
            if (i < this->lower.size()) {
                return this->lower.at(i);
            }
            else {
                return this->upper.at(this->size()-i);
            }
        }

        int rightmost_index_cw() {
            return this->upper.size() - 1;
        }

        int rightmost_index_ccw() {
            return this->lower.size() - 1;
        }

        int size() {
            if (this->lower.size() + this->upper.size() == 0) return 0;
            else if (this->lower.size() + this->upper.size() == 2) return 1;
            else return this->lower.size() + this->upper.size() - 2;
        }

        void clear() {
            this->upper.clear();
            this->lower.clear();
        }
};

#endif