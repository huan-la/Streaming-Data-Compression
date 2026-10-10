#ifndef ALGEBRAIC_FUNCTION_HPP
#define ALGEBRAIC_FUNCTION_HPP

#include <string>
#include <cmath>

struct Point2D {
    double x;
    double y;

    Point2D() {
        this->x = 0;
        this->y = 0;
    }

    Point2D(double x, double y) {
        this->x = x;
        this->y = y;
    }
};

// Geometry function support timeseries analysis
class Line {
    private:
        double slope;
        double intercept;

    public:
        Line() {
            this->slope = 0;
            this->intercept = 0;
        }

        Line(double slope, double intercept) {
            this->slope = slope;
            this->intercept = intercept;
        }

        double subs(double x) {
            return this->slope * x + this->intercept;
        }

        void set(double slope, double intercept) {
            this->slope = slope;
            this->intercept = intercept;
        }

        double get_slope() {
            return this->slope;
        }

        double get_intercept() {
            return this->intercept;
        }

        double get_root(double shift = 0, double scale = 1) {
            if (this->slope == 0) return INFINITY;
            else return (shift - scale * this->intercept) / this->slope;
        }

        static Line line(Point2D p1, Point2D p2) {
            double slope = (p1.y - p2.y) / (p1.x - p2.x);
            double intercept = p1.y - slope*p1.x;

            return Line(slope, intercept);
        }

        static Line line(double slope, Point2D p) {
            return Line(slope, p.y - slope*p.x);
        }

        static Point2D intersection(Line l1, Line l2) {
            double x = (l1.intercept - l2.intercept) / (l2.slope - l1.slope);
            double y = l1.intercept + l1.slope*x;

            return Point2D(x, y); 
        }
};

// Source for manipulating polynomial function
class Polynomial {
    public:
        int degree;
        std::vector<double> coefficients;    // coefficient degree starts from 0 

    public:
        template <typename T>
        Polynomial(int k, const T* coeffs, bool reverse = false) {
            degree = k;
            this->coefficients.resize(k + 1);

            for (int i = 0; i <= k; i++) {
                int idx = reverse ? k - i : i;
                this->coefficients[idx] = static_cast<double>(coeffs[i]);
            }
        }

        template <typename T>
        Polynomial(int k, const std::vector<T>& coeffs, bool reverse = false) {
            degree = k;
            this->coefficients.resize(k + 1);

            for (int i = 0; i <= k; i++) {
                int idx = reverse ? k - i : i;
                this->coefficients[idx] = static_cast<double>(coeffs[i]);
            }
        }

        Polynomial(float coeff) {
            this->degree = 0;
            this->coefficients.push_back((double) coeff);
        }

        Polynomial(double coeff) {
            this->degree = 0;
            this->coefficients.push_back(coeff);
        }

        double subs(double indeterminate) const {
            double result = this->coefficients[0];
            for (int i=1; i<this->degree+1; i++) {
                result += this->coefficients[i]*pow(indeterminate, i);
            }

            return result;
        }

        template <typename T>
        void set_coefficient(T* coeffs, bool reverse = false) {
            for (int i = 0; i <= this->degree; i++) {
                int idx = reverse ? this->degree - i : i;
                this->coefficients[idx] = static_cast<double>(coeffs[i]);
            }
        }

        double get_coefficient(int degree) {
            return this->coefficients[degree];
        }

        std::string str() const {
            std::string s = "";
            for (int i=this->degree; i>0; i--) {
                s += std::to_string(this->coefficients[i]) + "x^" + std::to_string(i) + " + ";
            }

            return s + std::to_string(this->coefficients[0]);
        }
};

#endif