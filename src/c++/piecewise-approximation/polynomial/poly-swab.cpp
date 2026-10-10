#include "piecewise-approximation/polynomial.hpp"

namespace PolySwab {
    std::vector<double> Approximator::__interpolate(int degree, std::vector<double>& data) {
        std::vector<int> indices;
        // Get evenly spaced indices between starting and ending data points
        for (int i = 0; i <= degree; i++) {
            int idx = (long)i * (data.size() - 1) / degree;
            indices.push_back(idx);
        }

        Eigen::MatrixXd A(degree+1, degree+1);
        Eigen::VectorXd b(degree+1);
        for (int i = 0; i < indices.size(); i++) {
            b(i) = data[indices[i]];
            for (int j = 0; j <= degree; j++) {
                A(i, j) = std::pow(indices[i], j);
            }
        }
        Eigen::VectorXd c = A.colPivHouseholderQr().solve(b);

        return std::vector<double>(c.data(), c.data() + c.size());
    }

    // For runtime optimization: Use cached (previous) model if available.
    // Only approximate or compute a new model when the previous one cannot be reused.
    std::vector<double> Approximator::__regression(int degree, std::vector<double>& data) {
        Eigen::MatrixXd A(data.size(), degree + 1);
        Eigen::VectorXd b(data.size());
        for (int i = 0; i < data.size(); ++i) {                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
            b(i) = data[i];
            for (int j = 0; j <= degree; ++j) {
                A(i, j) = std::pow(i, j);
            }
        }
        Eigen::VectorXd c = (((A.transpose() * A).inverse()) * A.transpose()) * b;
        
        return std::vector<double>(c.data(), c.data() + c.size());
    }

    Polynomial Approximator::approximate(std::string mode, int degree, std::vector<double>& data) {
        std::vector<double> coeffs;

        if (mode == "regression") 
            coeffs = Approximator::__regression(degree, data);
        else if (mode == "interpolate") 
            coeffs = Approximator::__interpolate(degree, data);
    
        return Polynomial(degree, coeffs);
    }

    // Mean square error calculate
    double Approximator::cal_error(std::vector<double>& segment, Polynomial& model) {
        double error = 0;
        for (int i=0; i<segment.size(); i++) {
            error += std::pow(segment[i]-model.subs(i), 2);
        }

        return error / segment.size();
    }

    // Verify new line satisfies infinity bound or not
    bool Grouper::bound_check(std::vector<double>& segment, Polynomial& model, double error, bool skip) {
        if (skip) {
            Point2D p(segment.size()-1, segment.back());
            return std::abs(p.y-model.subs(p.x)) <= error;
        }
        else {
            for (int i=0; i<segment.size(); i++) {
                // Immediately terminate when individual error exceed the allowable threshold
                if (std::abs(segment[i]-model.subs(i)) > error) return false;
            }
            return true;
        }
    }

    void Grouper::merge(std::vector<double>& s1, std::vector<double>& s2, std::string mode) {
        int offset = mode == "interpolate" ? 1 : 0;
        s1.insert(s1.end(), s2.begin() + offset, s2.end());
    }

    double Grouper::merge_cost(std::vector<double>& s1, std::vector<double>& s2, std::string mode, int degree, double error) {
        int offset = mode == "interpolate" ? 1 : 0;
        std::vector<double> s = s1;
        s.insert(s.end(), s2.begin() + offset, s2.end());
        Polynomial model = Approximator::approximate(mode, degree, s);

        if (!Grouper::bound_check(s, model, error)) return INFINITY;
        else return Approximator::cal_error(s, model);
    }

    // Begin: compression
    void Compression::__bottom_up() {
        for (int i=this->m_err.size(); i<this->segments.size()-1; i++) {
            this->m_err.push_back(Grouper::merge_cost(
                this->segments[i], this->segments[i+1], 
                this->mode, this->degree, this->error
            ));
        }

        auto it = std::min_element(this->m_err.begin(), this->m_err.end());
        while (*it != INFINITY && this->segments.size() > 1) {
            int index = std::distance(this->m_err.begin(), it);
            Grouper::merge(this->segments[index], this->segments[index+1], this->mode);

            this->m_err.erase(this->m_err.begin() + index);
            this->segments.erase(this->segments.begin() + index + 1);

            this->models.erase(this->models.begin() + index + 1);
            Polynomial n_model = Approximator::approximate(this->mode, this->degree, this->segments[index]);
            this->models[index] = n_model;

            if (index != 0) {
                this->m_err[index-1] = Grouper::merge_cost(
                    this->segments[index-1], this->segments[index],
                    this->mode, this->degree, this->error
                );
            }
            if (index < this->m_err.size()) {
                this->m_err[index] = Grouper::merge_cost(
                    this->segments[index], this->segments[index+1],
                    this->mode, this->degree, this->error
                );
            }
            
            it = std::min_element(this->m_err.begin(), this->m_err.end());  
        }
    }

    bool Compression::__sliding_window() {
        if (this->models.size() == this->segments.size()) {
            Polynomial model = Approximator::approximate(this->mode, this->degree, this->window);
            this->models.push_back(model);
        }
        else if (this->mode == "regression") {
            Polynomial& model = this->models.back();
            if (Grouper::bound_check(this->window, model, this->error, true)) return true;
        }

        Polynomial n_model = Approximator::approximate(this->mode, this->degree, this->window);
        if (!Grouper::bound_check(this->window, n_model, this->error)) return false;
        
        this->models[this->models.size()-1] = n_model;
        return true;
    }

    void Compression::initialize(int count, char** params) {
        this->error = atof(params[0]);
        this->mode = params[1];
        this->degree = atoi(params[2]);
        this->n_segment = atoi(params[3]);
    }

    void Compression::finalize() {
        if (this->window.size() >= this->degree + 1) {
            this->segments.push_back(this->window);
        }

        this->__bottom_up();
        while (this->segments.size() != 0) {
            this->yield();
            this->models.erase(this->models.begin());
            this->segments.erase(this->segments.begin());
        }
        this->m_err.clear();
    }

    void Compression::compress(Univariate* data) {
        this->window.push_back(data->get_value());

        if (this->window.size() > this->degree) {
            if (!this->__sliding_window()) {
                this->window.pop_back();
                this->segments.push_back(this->window);
                
                if (this->mode == "interpolate") {
                    // Add last data point of previous segment to ensure connectivity
                    double tail = this->window.back();
                    this->window.clear();
                    this->window.push_back(tail);
                    this->window.push_back(data->get_value());

                    if (this->window.size() == this->degree + 1) {
                        this->__sliding_window();
                    }
                }
                else if (this->mode == "regression") {
                    this->window.clear();
                    this->window.push_back(data->get_value());
                }
            }
        }

        // Perform bottom up once n_segment sliding window segments have accumulated
        if (this->segments.size() == this->n_segment) {
            this->__bottom_up();
            this->yield();
            this->models.erase(this->models.begin());
            this->segments.erase(this->segments.begin());
            if (!this->m_err.empty()) this->m_err.erase(this->m_err.begin());
        }
    }

    BinObj* Compression::serialize() {
        BinObj* obj = new BinObj;
        std::vector<double>& segment = this->segments[0];
        Polynomial& model = this->models[0];

        if (this->mode == "interpolate") {
            if (this->first) {
                obj->put((float) model.coefficients[0]);
                this->first = false;
            }
            
            obj->put(VariableByteEncoding::encode(segment.size()));
            for (int i = 1; i <= this->degree; i++) {
                obj->put((float) model.coefficients[i]);
            }
        }
        else if (this->mode == "regression") {
            obj->put(VariableByteEncoding::encode(segment.size()));
            for (int i = 0; i <= this->degree; i++) {
                obj->put((float) model.coefficients[i]);
            }
        }
        
        return obj;
    }
    // End: compression

    // Begin: decompression
    void Decompression::initialize(int count, char** params) {
        this->degree = atoi(params[2]);   
        this->mode = params[1];
    }

    void Decompression::finalize() {
        // Do nothing
    }

    CSVObj* Decompression::decompress(BinObj* compress_data) {
        CSVObj* base_obj = nullptr;
        CSVObj* prev_obj = nullptr;

        int start = 0;
        unsigned long length = 0;
        Polynomial* model = nullptr;
        float* coefficients = new float[this->degree+1];
        int offset = this->mode == "interpolate" ? 1 : 0;

        if (this->mode == "interpolate") {
            start = 1;
            if (this->pivot == INFINITY) {
                start = 0;
                this->pivot = compress_data->getFloat();
            }

            length = VariableByteEncoding::decode(compress_data);
            for (int i = 1; i <= this->degree; i++) {
                coefficients[i] = compress_data->getFloat();
            }
            coefficients[0] = this->pivot;
            model = new Polynomial(this->degree, coefficients);            
        }
        else if (this->mode == "regression") {
            length = VariableByteEncoding::decode(compress_data);
            for (int i = 0; i <= this->degree; i++) {
                coefficients[i] = compress_data->getFloat();
            }
            model = new Polynomial(this->degree, coefficients);
        }

        for (int i=start; i<length; i++) {
            if (base_obj == nullptr) {
                base_obj = new CSVObj;
                base_obj->pushData(std::to_string(this->basetime + i * this->interval));
                base_obj->pushData(std::to_string(model->subs(i)));

                prev_obj = base_obj;
            }
            else {
                CSVObj* obj = new CSVObj;
                obj->pushData(std::to_string(this->basetime + i * this->interval));
                obj->pushData(std::to_string(model->subs(i)));

                prev_obj->setNext(obj);
                prev_obj = obj;
            }
        }
        
        this->pivot = model->subs(length-1);
        this->basetime += (length - offset) * this->interval;

        delete model;
        delete[] coefficients;

        return base_obj;
    }
    // End: decompression
};