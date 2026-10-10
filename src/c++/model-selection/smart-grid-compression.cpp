#include "model-selection/model-selection.hpp"

namespace SmartGridCompression {
    // Begin: material
    ConstantModel::ConstantModel() {
        this->degree = 0;
        this->length = 0;
        
        this->value = 0;
        this->min = INFINITY;
        this->max = -INFINITY;
    }

    double ConstantModel::getValue() {
        return this->value;
    }

    double ConstantModel::getCompressionRatio() {
        if (this->length == 0) return -1;       // Ignore this model
        else {
            return (double) this->length / 7;
        } 
    }

    void ConstantModel::clear() {
        this->degree = 0;
        this->length = 0;
        
        this->min = INFINITY;
        this->max = -INFINITY;
    }

    bool ConstantModel::approximate(double bound, std::vector<Point2D>& segment) {
        if (this->length == 0) {
            for (Point2D& p : segment) {
                this->min = this->min < p.y ? this->min : p.y;
                this->max = this->max > p.y ? this->max : p.y;
                
                if (this->max - this->min <= 2 * bound) {
                    this->value = (this->max + this->min) / 2;
                    this->length++;
                }
                else {
                    return false;
                }
            }
            
            return true;
        }
        else {
            Point2D& p = segment.back();
            this->min = this->min < p.y ? min : p.y;
            this->max = this->max > p.y ? max : p.y;
            
            if (max - min <= 2 * bound) {
                this->value = (max + min) / 2;
                this->length++;
                return true;
            }
            else {
                return false;
            }
        }
    }

    bool LinearModel::__verify(double bound, Line& line) {
        if (this->cvx.size() <= 2) {        
            for (int i=0; i<this->cvx.size(); i++) {
                Point2D p = this->cvx.at(i);
                if (std::abs(line.subs(p.x) - p.y) > bound) return false;
            }

            return true;
        }
        else {
            double A_dis = std::abs(line.subs(this->A.x) - this->A.y);
            double B_dis = std::abs(line.subs(this->B.x) - this->B.y);
            double C_dis = std::abs(line.subs(this->C.x) - this->C.y);

            return std::max({ A_dis, B_dis, C_dis} ) <= bound;
        }
    }

    double LinearModel::__distance(Line line, Point2D p) {
        return std::abs(p.y - line.subs(p.x)) 
            / sqrt(line.get_slope() * line.get_slope() + 1);
    }

    int LinearModel::__x_external(int x, int x1, int x2) {
        // 0:  x_external
        // 1:  x_external to the right
        // -1: x_external to the left

        if (x >= x1 && x >= x2) return 1;
        else if (x <= x1 && x <= x2) return -1;
        else return 0;
    }

    int LinearModel::__search(int side_index, int prev_v_l) {
        Line line = this->cvx.side(side_index);
        int v = -1; double prev_dis = -INFINITY;
        
        // search in ccw order
        // if pj(i+1) < pj(i) return i
        for (int i=prev_v_l; i<=this->cvx.size(); i++) {
            double dis = this->__distance(line, this->cvx.at_ccw(i));            
            if (prev_dis < dis) {
                prev_dis = dis;
                v = i;
            }
            else break;
        }

        return v;
    }

    Line LinearModel::__approx(const std::vector<Point2D>& segment) {
        if (this->cvx.size() <= 2) return Line::line(this->cvx.at(0), this->cvx.at(1));
        else {
            // find the A, B and C
            Point2D n_A(0, 0), n_B(0,0), n_C(0,0);
            int prev_v_l = this->cvx.rightmost_index_ccw();

            // find the first v(l1)
            for (int i=0; i<this->cvx.size(); i++) {
                int v_l = this->__search(i, prev_v_l);                        
                int external = this->__x_external(this->cvx.at_ccw(v_l).x, this->cvx.at_ccw(i).x, this->cvx.at_ccw(i+1).x);

                if (external == 1) {
                    prev_v_l = v_l;
                }
                else if (external == 0) {
                    n_A = this->cvx.at_ccw(i);
                    n_B = this->cvx.at_ccw(i+1);
                    n_C = this->cvx.at_ccw(v_l);
                    break;
                }
                else {
                    n_C = this->cvx.at_ccw(i);
                    for (int j = prev_v_l; j < v_l; j++) {
                        if (this->__x_external(n_C.x, this->cvx.at_ccw(j).x, this->cvx.at_ccw(j+1).x) == 0) {
                           n_A = this->cvx.at_ccw(j);
                           n_B = this->cvx.at_ccw(j+1);
                           break; 
                        }                      
                    }
                    break;
                }
            }
            
            // find the optimal approximation line
            Line line = Line::line(n_A, n_B);
            double a = line.get_slope();
            double b = line.get_intercept();
            double c = (n_C.y - line.get_slope() * n_C.x + line.get_intercept()) / 2;
            
            this->A.x = n_A.x; this->A.y = n_A.y;
            this->B.x = n_B.x; this->B.y = n_B.y;
            this->C.x = n_C.x; this->C.y = n_C.y;

            return Line(a, c);
        }
    }

    LinearModel::LinearModel() {
        this->degree = 1;
        this->length = 0;
        this->line = nullptr;
    }

    Line* LinearModel::getLine() {
        return this->line;
    }

    double LinearModel::getCompressionRatio() {
        if (this->length == 0) return -1;       // Ignore this model
        else {
            return (double) this->length / 11;
        } 
    }

    void LinearModel::clear() {
        this->length = 0;
        this->cvx.clear();
        
        this->A.x = 0; this->A.y = 0;
        this->B.x = 0; this->B.y = 0;
        this->C.x = 0; this->C.y = 0;
        
        
        if (this->line != nullptr) {
            delete this->line;
            this->line = nullptr;
        }
    }

    bool LinearModel::approximate(double bound, std::vector<Point2D>& segment) {
        if (this->length == 0) {
            for (Point2D& p : segment) {
                this->cvx.append(p);
                if (this->cvx.size() > 1) {
                    Line line = this->__approx(segment);
                    if (this->__verify(bound, line)) {
                        if (this->line != nullptr) delete this->line;
                        this->line = new Line(line.get_slope(), line.get_intercept());
                    }
                    else return false; 
                }

                this->length++;
            }
            
            return true;
        }
        else {
            Point2D& p = segment.back();
            this->cvx.append(p);

            if (this->cvx.size() > 1) {            
                Line line = this->__approx(segment);
                if (this->__verify(bound, line)) {
                    if (this->line != nullptr) delete this->line;
                    this->line = new Line(line.get_slope(), line.get_intercept());
                }
                else return false;  
            }
            this->length++;

            return true;
        }
    }

    PolynomialModel::PolynomialModel(int degree) {
        this->degree = degree;
        this->length = 0;   
        this->polynomial = nullptr;
    }

    Polynomial* PolynomialModel::getPolynomial() {
        return this->polynomial;
    }

    double PolynomialModel::getCompressionRatio() {
        if (this->length == 0) return -1;       // Ignore this model
        else {
            return (double) this->length / (3 + 4 * (this->degree + 1));
        } 
    }

    void PolynomialModel::clear() {
        this->length = 0;
        if (this->polynomial != nullptr) {
            delete this->polynomial;
            this->polynomial = nullptr;
        }
    }

    bool PolynomialModel::approximate(double bound, std::vector<Point2D>& segment) {
        if (this->polynomial != nullptr) {
            Eigen::VectorXd x(this->degree + 2);                    // decision variables
            Eigen::VectorXd c(this->degree + 2);                    // objective coefficients
            Eigen::MatrixXd A(2, this->degree + 2);  // constraint matrix
            Eigen::VectorXd b(2);                    // constraint bound

            for (int i=0; i<this->degree+2; i++) {
                if (i != this->degree + 1) c(i) = 0.0;
                else c(i) = 1.0;
            }

            for (int j=this->degree; j>=0; j--) {
                A(0, this->degree-j) = -pow(segment.back().x, j);
                A(1, this->degree-j) = pow(segment.back().x, j);
            }
            A(0, this->degree+1) = -1.0;
            A(1, this->degree+1) = -1.0;

            b(0) = -segment.back().y;
            b(1) = segment.back().y;
        

            double minobj = this->sdlp.warm_linprog(c, A, b, x);
            if (minobj == INFINITY || minobj == -INFINITY || minobj > bound) {
                return false;
            }
            else {
                float* coefficients = new float[this->degree+1];
                for (int i = 0; i <= this->degree; i++) {
                    coefficients[this->degree-i] = x(i);
                }

                this->polynomial->set_coefficient<float>(coefficients);
                this->length = segment.size();

                delete[] coefficients;
                return true;
            }
        }
        else {
            Eigen::VectorXd x(this->degree + 2);                    // decision variables
            Eigen::VectorXd c(this->degree + 2);                    // objective coefficients
            Eigen::MatrixXd A(2*segment.size(), this->degree + 2);  // constraint matrix
            Eigen::VectorXd b(2*segment.size());                    // constraint bound

            for (int i=0; i<this->degree+2; i++) {
                if (i != this->degree + 1) c(i) = 0.0;
                else c(i) = 1.0;
            }

            for (int i=0; i<segment.size(); i++) {
                for (int j=this->degree; j>=0; j--) {
                    A(2*i, this->degree-j) = -pow(segment[i].x, j);
                    A(2*i+1, this->degree-j) = pow(segment[i].x, j);
                }
                A(2*i, this->degree+1) = -1.0;
                A(2*i+1, this->degree+1) = -1.0;

                b(2*i) = -segment[i].y;
                b(2*i+1) = segment[i].y;
            }

            double minobj = this->sdlp.linprog(c, A, b, x);
            if (minobj == INFINITY || minobj == -INFINITY || minobj > bound) {
                return false;
            }
            else {
                float* coefficients = new float[this->degree+1];
                for (int i = 0; i <= this->degree; i++) {
                    coefficients[this->degree-i] = x(i);
                }

                this->polynomial = new Polynomial(this->degree, coefficients);
                this->length = segment.size();
                
                delete[] coefficients;
                return true;
            }
            
        }
    }
    // End: material

    // Begin: compression
    void Compression::initialize(int count, char** params) {
        this->error = atof(params[0]);
        this->max_degree = atoi(params[1]);

        for (int i=0; i<=this->max_degree; i++) {
            if (i == 0) this->models.push_back(new ConstantModel());
            else if (i == 1) this->models.push_back(new LinearModel());
            else this->models.push_back(new PolynomialModel(i));
        }
    }

    void Compression::finalize() {
        while (this->window.size() > 0) {
            this->__choose_best_model();
            this->yield();
            
            this->window.erase(this->window.begin(), this->window.begin() + this->models[this->chosen_model]->length);
            for (int i=0; i<this->window.size(); i++) this->window[i].x = i;
            for (int i=0; i<=this->max_degree; i++) this->models[i]->clear();
            this->curr_degree = 0; this->chosen_model = -1;

            if (this->window.size() > 0) {
                for (int k=this->curr_degree; k<=this->max_degree; k++) {
                    if (this->__approximate(k)) {
                        this->curr_degree = k;
                        break;
                    }
                }
            }
        }

        for (int i=0; i<=this->max_degree; i++) {
            this->models[i]->clear();
            delete this->models[i];
        }        
    }

    bool Compression::__approximate(int degree) {
        bool success = true;
        if (degree == 0) {
            ConstantModel* constantModel = (ConstantModel*) this->models[degree];
            success = constantModel->approximate(this->error, this->window);
        }
        else if (degree == 1) {
            LinearModel* linearModel = (LinearModel*) this->models[degree];
            success = linearModel->approximate(this->error, this->window);
        }
        else {
            if (this->window.size() > degree) {
                PolynomialModel* polynomialModel = (PolynomialModel*) this->models[degree];
                success = polynomialModel->approximate(this->error, this->window);
            }
        }
        return success;
    }

    void Compression::__choose_best_model() {
        int index = 0;
        double max_gain = -INFINITY;

        for (int i=0; i<this->models.size(); i++) {
            double gain = this->models[i]->getCompressionRatio();
            if (gain > max_gain) {
                index = i;
                max_gain = gain;
            }
        }
        this->chosen_model = index;
    }

    BinObj* Compression::serialize() {
        BinObj* obj = new BinObj;
        Model* model = this->models[this->chosen_model];

        int embedded = model->length << 2 | model->degree;

        if (model->degree == 0) {
            obj->put(VariableByteEncoding::encode(embedded));
            obj->put((float) ((ConstantModel*) model)->getValue());
        }
        else if (model->degree == 1) {
            Line* line = ((LinearModel*) model)->getLine();
            obj->put(VariableByteEncoding::encode(embedded));
            obj->put((float) line->get_intercept());
            obj->put((float) line->get_slope());
        }
        else {
            Polynomial* polynomial = ((PolynomialModel*) model)->getPolynomial();
            obj->put(VariableByteEncoding::encode(embedded));
            for (int i = 0; i <= polynomial->degree; i++) {
                obj->put((float) polynomial->coefficients[i]);
            }
        }

        return obj;
    }

    void Compression::compress(Univariate* data) {
        Point2D p(this->window.size(), data->get_value());
        this->window.push_back(p);

        bool success = false;
        for (int k=this->curr_degree; k<=this->max_degree; k++) {
            if (this->__approximate(k)) {
                this->curr_degree = k;
                success = true;
                break;
            }
        }

        if (!success) {
            this->__choose_best_model();
            this->yield();

            this->window.erase(this->window.begin(), this->window.begin() + this->models[this->chosen_model]->length);
            for (int i=0; i<this->window.size(); i++) this->window[i].x = i;

            for (int i=0; i<=this->max_degree; i++) {
                this->models[i]->clear();
            }

            this->curr_degree = 0;
            this->chosen_model = -1;
        }
    }
    // End: compression

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
}