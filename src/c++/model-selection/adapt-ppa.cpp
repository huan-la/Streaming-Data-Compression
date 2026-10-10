#include "model-selection/model-selection.hpp"

namespace AdaptPPA {
    // Begin: material
    LinearSegment::~LinearSegment() {
        if (this->p_start != nullptr) delete this->p_start;
        if (this->p_end != nullptr) delete this->p_end;
        if (this->pivot != nullptr) delete this->pivot;         
        if (this->u_line != nullptr) delete this->u_line;
        if (this->l_line != nullptr) delete this->l_line;
        this->u_cvx.clear(); this->l_cvx.clear();
    }

    int LinearSegment::getLength() {
        return this->p_end->x - this->p_start->x + 1;
    }

    Line LinearSegment::getLine() {
        return Line(
            (this->u_line->get_slope() + this->l_line->get_slope()) / 2,
            (this->u_line->get_intercept() + this->l_line->get_intercept()) / 2
        );
    }

    void LinearSegment::translation(Point2D* p) {
        Line* n_u_line = new Line(this->u_line->get_slope(), this->u_line->subs(p->x));
        Line* n_l_line = new Line(this->l_line->get_slope(), this->l_line->subs(p->x));
        
        this->p_end->x = this->getLength() - 1;
        this->p_start->x = 0;
        
        delete this->u_line; this->u_line = n_u_line;
        delete this->l_line; this->l_line = n_l_line;
    }

    void LinearSegment::approximate(Point2D& p, double error) {        
        if (this->pivot == nullptr) {
            this->p_start = new Point2D(p.x, p.y);
            this->pivot = new Point2D(p.x, p.y);
            this->u_cvx.append(Point2D(p.x, p.y - error));
            this->l_cvx.append(Point2D(p.x, p.y + error));
        }
        else if (this->u_line == nullptr) {
            Line u_line = Line::line(
                Point2D(this->pivot->x, this->pivot->y - error), 
                Point2D(p.x, p.y + error)
            );
            Line l_line = Line::line(
                Point2D(this->pivot->x, this->pivot->y + error), 
                Point2D(p.x, p.y - error)
            );

            this->u_line = new Line(u_line.get_slope(), u_line.get_intercept());
            this->l_line = new Line(l_line.get_slope(), l_line.get_intercept());
            this->u_cvx.append(Point2D(p.x, p.y - error));
            this->l_cvx.append(Point2D(p.x, p.y + error));
        }
        else {
            if (this->l_line->subs(p.x) - p.y - error > EPS 
                || p.y - error - this->u_line->subs(p.x) > EPS) 
            {
                this->p_end = new Point2D(p.x - 1, p.y);
                this->is_complete = true;
            }
            else {
                bool update_u = p.y + error < this->u_line->subs(p.x);
                bool update_l = p.y - error > this->l_line->subs(p.x);
                
                if (update_u) {
                    int index = -1;
                    double min_slp = INFINITY;

                    for (int i=0; i<this->u_cvx.size(); i++) {
                        Line line = Line::line(this->u_cvx.at(i), Point2D(p.x, p.y + error));
                        if (line.get_slope() < min_slp) {
                            index = i;
                            min_slp = line.get_slope();
                            this->u_line->set(line.get_slope(), line.get_intercept());
                        }
                    }

                    if (index >= 0) this->u_cvx.erase_from_begin(index);
                }
                if (update_l) {
                    int index = -1;
                    double max_slp = -INFINITY;

                    for (int i=0; i<this->l_cvx.size(); i++) {
                        Line line = Line::line(this->l_cvx.at(i), Point2D(p.x, p.y - error));
                        if (line.get_slope() > max_slp) {
                            index = i;
                            max_slp = line.get_slope();
                            this->l_line->set(line.get_slope(), line.get_intercept());
                        }
                    }

                    if (index >= 0) this->l_cvx.erase_from_begin(index);
                }

                if (update_u) this->l_cvx.append(Point2D(p.x, p.y + error));
                if (update_l) this->u_cvx.append(Point2D(p.x, p.y - error));
            }
        }
    }
    // End: material 

    // Begin: compression
    bool Compression::__merge_check() {
        bool flag = true;
        Line line_1 = this->seg_1->getLine();
        Line line_2 = this->seg_2->getLine();        

        if (this->p2->x <= this->p1->x || this->p2->x >= this->p3->x) flag = false;
        else {
            Eigen::MatrixXd A(3, 3);
            A << this->p1->x*this->p1->x, this->p1->x, 1,
                 this->p2->x*this->p2->x, this->p2->x, 1,
                 this->p3->x*this->p3->x, this->p3->x, 1;
        
            Eigen::VectorXd b(3);
            b << this->p1->y, this->p2->y, this->p3->y;
        
            Eigen::VectorXd x = A.colPivHouseholderQr().solve(b);

            double extreme_x = (-x(1)) / (2*x(0));
            double extreme_y = x(0)*extreme_x*extreme_x + x(1)*extreme_x + x(2);
            
            if (line_1.get_slope() * line_2.get_slope() > 0) {
                if (extreme_x > this->p1->x && extreme_x < this->p3->x) flag = false;
            }
            else {  
                if (std::abs(extreme_y - this->p2->y) > this->error) flag = false;
            }
        }

        if (flag) {
            int n_direction = this->seg_1->getLine().get_slope() > this->seg_2->getLine().get_slope() ? -1 : 1;
            if (this->direction != n_direction) {
                this->direction = n_direction;
                this->degree += 1;
                if (this->degree > this->max_degree) {
                    flag = false;
                    this->degree = this->max_degree;
                }
            }
        }

        return flag;
    }

    double* Compression::__approximate(int pivot) {
        int n = this->window.size() - pivot;
        Eigen::MatrixXd A(n, this->degree + 1);
        Eigen::VectorXd b(n);

        for (int i = 0; i < n; ++i) {                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
            b(i) = this->window[i].y;
            for (int j = 0; j <= this->degree; ++j) {
                A(i, j) = std::pow(this->window[i].x, j);
            }
        }

        double* coeffs = new double[this->degree+1];
        Eigen::VectorXd coefficients = (((A.transpose() * A).inverse()) * A.transpose()) * b;
        for (int i=0; i<=this->degree; i++) coeffs[i] = coefficients(i);

        return coeffs;
    }

    void Compression::initialize(int count, char** params) {
        this->error = atof(params[0]);
        this->max_degree = atoi(params[1]);

        this->p1 = new Point2D(-1, -1);
        this->p2 = new Point2D(-1, -1);
        this->p3 = new Point2D(-1, -1);
        this->seg_1 = new LinearSegment();
    }

    void Compression::finalize() {
        this->seg_2->is_complete = true;
        this->seg_2->p_end = new Point2D(this->window.back().x, this->window.back().y);
        
        if (this->seg_2->getLength() < 2) {
            this->yield();
        }
        else {
            this->p3->x = this->seg_2->p_end->x;
            this->p3->y = this->seg_2->getLine().subs(this->p3->x);

            bool flag = this->__merge_check();
            int n_direction = this->seg_1->getLine().get_slope() > this->seg_2->getLine().get_slope() ? -1 : 1;

            if (flag) {
                if (direction != n_direction) {
                    direction = n_direction;
                    this->degree += 1;
                }
                if (this->degree > this->max_degree) {
                    flag = true;
                    this->degree = this->max_degree;
                }
            }

            if (flag) {
                this->yield();
            }
            else {
                this->yield();
                this->degree = 1;
                
                delete this->seg_1;
                this->seg_1 = seg_2;
                this->seg_1->translation(this->p2);

                this->seg_2 = nullptr;
                this->yield();
            }
        }
        
        if (this->p1 != nullptr) delete this->p1;
        if (this->p2 != nullptr) delete this->p2;
        if (this->p3 != nullptr) delete this->p3;
        if (this->seg_1 != nullptr) delete this->seg_1;
        if (this->seg_2 != nullptr) delete this->seg_2;
    }

    BinObj* Compression::serialize() {
        BinObj* obj = new BinObj;
        
        if (this->degree == 1) {
            Line line = this->seg_1->getLine();
            int embedded = this->seg_1->getLength() << 2 | 1;

            obj->put(VariableByteEncoding::encode(embedded));
            obj->put((float) line.get_intercept());
            obj->put((float) line.get_slope());
        }
        else {
            int pivot = this->seg_2->is_complete ? this->seg_2->getLength() : 0;
            int embedded = (this->window.size() - pivot) << 2 | this->degree;
            double* coeffs = this->__approximate(pivot);
            
            obj->put(VariableByteEncoding::encode(embedded));
            for (int i = 0; i <= this->degree; i++) {
                obj->put((float) coeffs[i]);
            }
            delete[] coeffs;
        }

        return obj;
    }
    
    void Compression::compress(Univariate* data) {
        Point2D p(this->window.size(), data->get_value());

        if (!this->seg_1->is_complete) {
            this->seg_1->approximate(p, this->error);
            if (this->seg_1->is_complete) {
                this->p1->x = 0; this->p1->y = this->seg_1->getLine().subs(0);

                this->seg_2 = new LinearSegment();
                this->seg_2->approximate(p, this->error);
            }
        }
        else if (!this->seg_2->is_complete) {
            this->seg_2->approximate(p, this->error);
            if (this->seg_2->is_complete) {
                this->p3->x = this->seg_2->p_end->x;
                this->p3->y = this->seg_2->getLine().subs(this->seg_2->p_end->x);
                Point2D intersection = Line::intersection(this->seg_1->getLine(), this->seg_2->getLine());
                this->p2->x = intersection.x;
                this->p2->y = intersection.y;

                if (this->__merge_check()) {
                    delete this->seg_1;
                    this->seg_1 = seg_2;
                    this->seg_2 = new LinearSegment();
                    this->seg_2->approximate(p, this->error);

                    this->p1->x = this->seg_1->p_start->x; this->p1->y = this->seg_1->getLine().subs(this->p1->x);
                }
                else {
                    this->yield();
                    this->window = { this->window.end() - this->seg_2->getLength(), this->window.end() };
                    for (int i=0; i<this->window.size(); i++) this->window[i].x = i; 
                    
                    delete this->seg_1; 
                    this->seg_1 = this->seg_2; 
                    this->seg_1->translation(this->seg_2->p_start);

                    this->p1->x = 0; this->p1->y = this->seg_1->getLine().subs(0);

                    p = Point2D(this->window.size(), p.y);
                    this->seg_2 = new LinearSegment();
                    this->seg_2->approximate(p, this->error);

                    this->degree = 1; 
                    this->direction = 0;
                }
            }
        }
        this->window.push_back(p);
    }

    // Begin: decompression
    void Decompression::initialize(int count, char** params) {
        // Do nothing
    }

    void Decompression::finalize() {
        // Do nothing
    }

    CSVObj* Decompression::decompress(BinObj* compress_data) {
        CSVObj* base_obj = nullptr;
        CSVObj* prev_obj = nullptr;

        long embedded = VariableByteEncoding::decode(compress_data);
        int degree = (embedded & 3);
        int length = embedded >> 2;

        float* coefficients = new float[degree+1];
        for (int i = 0; i <= degree; i++) {
            coefficients[i] = compress_data->getFloat();
        }

        Polynomial polynomial(degree, coefficients);
        for (int i = 0; i < length; i++) {
            if (base_obj == nullptr) {
                base_obj = new CSVObj;
                base_obj->pushData(std::to_string(this->basetime + i * this->interval));
                base_obj->pushData(std::to_string(polynomial.subs(i)));

                prev_obj = base_obj;
            }
            else {
                CSVObj* obj = new CSVObj;
                obj->pushData(std::to_string(this->basetime + i * this->interval));
                obj->pushData(std::to_string(polynomial.subs(i)));

                prev_obj->setNext(obj);
                prev_obj = obj;
            }
        }

        this->basetime += length * this->interval;

        delete[] coefficients;
        return base_obj;
    }
    // End: decompression
};