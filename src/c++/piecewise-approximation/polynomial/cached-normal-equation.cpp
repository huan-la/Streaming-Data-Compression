#include "piecewise-approximation/polynomial.hpp"

namespace CachedNormalEquation {

    Model::Model(Polynomial* function) {
        this->function = function;
    }

    Model::~Model() {
        delete this->function;
        this->error = -1;
        this->length = 0;
    }

    // Begin: compression
    // Using eigen library
    Polynomial* Compression::__calPolynomial() {
        Eigen::VectorXd theta;
        if (this->cache.find(this->window.size()) == this->cache.end()) {
            Eigen::MatrixXd X(this->window.size(), this->degree + 1);
            Eigen::VectorXd y(this->window.size());

            for (int i = 0; i < this->window.size(); i++) {
                for (int k=0; k<=this->degree; k++) {
                    X(i, k) = std::pow(i, k);
                } 
                y(i) = this->window[i].y;
            }

            Eigen::MatrixXd X_T_X_inv_X_T = (X.transpose() * X).inverse() * X.transpose();
            this->cache.insert({this->window.size(), X_T_X_inv_X_T});

            theta = X_T_X_inv_X_T * y;
        }
        else {
            Eigen::VectorXd y(this->window.size());
            for (int i=0; i<this->window.size(); i++) {
                y(i) = this->window[i].y;
            }

            Eigen::MatrixXd X_T_X_inv_X_T = this->cache.find(this->window.size())->second;
            theta = X_T_X_inv_X_T * y;
        }
        
        float* coeffs = new float[this->degree+1];
        for (int i = 0; i <= this->degree; ++i) {
            coeffs[i] = theta[i];
        }

        Polynomial* model = new Polynomial(degree, coeffs);
        
        delete[] coeffs;
        return model;
    }

    bool Compression::__approxSuccess(Model* model) {
        if (this->mode == "individual") {
            if (model->error == -1) {
                model->error = -INFINITY;
                for (int i = 0; i < this->window.size(); i++) {
                    float error = std::abs(model->function->subs(i) - this->window[i].y);
                    model->error = model->error < error ? error : model->error;
                }
            }
            else {
                float error = std::abs(model->function->subs(this->window.size()-1) - this->window[this->window.size()-1].y);
                model->error = model->error < error ? error : model->error;
            }
            
            return model->error < this->error;
        }
        else if (mode == "accumulate") {
            if (model->error == -1) {
                model->error = 0;
                for (int i = 0; i < this->window.size(); i++) {
                    model->error += std::abs(model->function->subs(i) - this->window[i].y);
                }
            }
            else {
                model->error += std::abs(model->function->subs(this->window.size()-1) - this->window[this->window.size()-1].y);
            }
            
            return model->error < this->error * this->window.size();
        }

        return false;
    }

    void Compression::initialize(int count, char** params) {
        this->error = atof(params[0]);
        this->mode = params[1];
        this->degree = atoi(params[2]);
    }

    void Compression::finalize() {
        if (this->window.size() >= this->degree + 1) {
            this->yield();
            delete this->model;
        }
    }

    void Compression::compress(Univariate* data) {
        this->window.push_back(Point2D(this->window.size(), data->get_value()));

        if (this->window.size() == this->degree + 1) {
            model = new Model(__calPolynomial());
            this->model->length = this->window.size();
        }
        if (this->window.size() > this->degree + 1) {
            if (!__approxSuccess(this->model)) {
                Model* n_model = new Model(__calPolynomial());
                if (!__approxSuccess(n_model)) {
                    this->yield();
                    window = {Point2D(0, data->get_value())};

                    delete this->model;
                    delete n_model;
                }
                else {
                    delete this->model;
                    this->model = n_model;
                    this->model->length = this->window.size();
                }
            }
            else {
                this->model->length = this->window.size();
            }
        }
    }

    BinObj* Compression::serialize() {
        BinObj* obj = new BinObj;
        
        obj->put(VariableByteEncoding::encode(this->model->length));
        for (int i = 0; i <= this->degree; i++) {
            obj->put((float) this->model->function->coefficients[i]);
        }

        return obj;
    }
    // End: compression

    // Begin: decompression
    void Decompression::initialize(int count, char** params) {
        this->degree = atoi(params[2]);   
    }

    void Decompression::finalize() {
        // Do nothing
    }

    CSVObj* Decompression::decompress(BinObj* compress_data) {
        CSVObj* base_obj = nullptr;
        CSVObj* prev_obj = nullptr;

        unsigned long length = VariableByteEncoding::decode(compress_data);
        float* coefficients = new float[this->degree+1];

        for (int i = 0; i <= this->degree; i++) {
            coefficients[i] = compress_data->getFloat();
        }
        Polynomial model(this->degree, coefficients);

        for (int i=0; i<length; i++) {
            if (base_obj == nullptr) {
                base_obj = new CSVObj;
                base_obj->pushData(std::to_string(this->basetime + i * this->interval));
                base_obj->pushData(std::to_string(model.subs(i)));

                prev_obj = base_obj;
            }
            else {
                CSVObj* obj = new CSVObj;
                obj->pushData(std::to_string(this->basetime + i * this->interval));
                obj->pushData(std::to_string(model.subs(i)));

                prev_obj->setNext(obj);
                prev_obj = obj;
            }
        }

        delete[] coefficients;
        this->basetime += length * this->interval;
        return base_obj;
    }
    // End: decompression
};