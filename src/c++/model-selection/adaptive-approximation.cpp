#include "model-selection/model-selection.hpp"

namespace AdaptiveApproximation {
    // Begin: material
    Model::~Model() {
        if (this->coeffs != nullptr) delete this->coeffs;
        if (this->pivot != nullptr) delete this->pivot;
    }

    int Model::getDim() {
        if (this->type == "exponent") return 1;
        else if (this->type == "linear") return 1;
        else if (this->type == "quadratic") return 2;
        else if (this->type == "cubic") return 3;

        return 0;
    }

    int Model::getLen() {
        return this->length;
    }
    
    double* Model::getCoeffs() {
        return this->coeffs;
    }

    std::string Model::getType() {
        return this->type;
    }

    ExpoFunction::ExpoFunction() {
        this->type = "exponent";
        this->coeffs = new double[2];
    }

    double ExpoFunction::subs(double x) {
        return this->coeffs[1] * exp(this->coeffs[0] * x);
    }
    
    bool ExpoFunction::approximate(double bound, double data) {
        Point2D p(this->length, data);

        if (this->pivot == nullptr) {
            this->coeffs[1] = p.y;
            this->pivot = new Point2D(p.x, p.y);
        }
        else {
            double n_upper = (log(p.y + bound) - log(this->pivot->y)) / (p.x - this->pivot->x);
            double n_lower = (log(p.y - bound) - log(this->pivot->y)) / (p.x - this->pivot->x);
            
            if (n_upper < this->lower || n_lower > this->upper) { 
                this->complete = true; 
                return false; 
            }
            else {
                this->upper = n_upper < this->upper ? n_upper : this->upper;
                this->lower = n_lower > this->lower ? n_lower : this->lower;
                this->coeffs[0] = (this->upper + this->lower) / 2;
            }
        }

        this->length++;
        return true;
    }

    LinearFunction::LinearFunction() {
        this->type = "linear";
        this->coeffs = new double[2];
    }

    double LinearFunction::subs(double x) {
        return this->coeffs[0] * x + this->coeffs[1];
    }
    
    bool LinearFunction::approximate(double bound, double data) {
        Point2D p(this->length, data);

        if (this->pivot == nullptr) {
            this->coeffs[1] = p.y;
            this->pivot = new Point2D(p.x, p.y);
        }
        else {
            double n_upper = (p.y - this->pivot->y + bound) / (p.x - this->pivot->x);
            double n_lower = (p.y - this->pivot->y - bound) / (p.x - this->pivot->x);

            if (n_upper < this->lower || n_lower > this->upper) { 
                this->complete = true; 
                return false; 
            }
            else {
                this->upper = n_upper < this->upper ? n_upper : this->upper;
                this->lower = n_lower > this->lower ? n_lower : this->lower;
                this->coeffs[0] = (this->upper + this->lower) / 2;
            }
        }

        this->length++;
        return true;
    }

    PolyFunction::PolyFunction(int degree, std::string type) {
        this->type = type;
        this->degree = degree;
        this->coeffs = new double[degree+1];
    }

    double PolyFunction::subs(double x) {
        Polynomial func(this->degree, this->coeffs, true);
        return func.subs(x);
    }

    bool PolyFunction::approximate(double bound, double data) {
        Point2D p(this->length, data);

        if (this->pivot == nullptr) {
            this->coeffs[this->degree] = p.y;
            this->pivot = new Point2D(p.x, p.y);
        }
        else {
            Eigen::VectorXd x(this->degree);                    
            Eigen::VectorXd c(this->degree);                    
            Eigen::MatrixXd A(2, this->degree);  
            Eigen::VectorXd b(2);                   

            for (int i=0; i<this->degree; i++) c(i) = 0.0;
            for (int j=this->degree; j>0; j--) {
                A(0, this->degree-j) = pow(p.x, j);
                A(1, this->degree-j) = -pow(p.x, j);
            }
            
            b(0) = p.y-this->pivot->y+bound;
            b(1) = this->pivot->y-p.y+bound;

            double minobj = INFINITY;
            if (this->length == 1) minobj = this->sdlp.linprog(c, A, b, x);
            else minobj = this->sdlp.warm_linprog(c, A, b, x);

            if (minobj == INFINITY) { this->complete = true; return false; }
            else for (int i=0; i<this->degree; i++) this->coeffs[i] = x(i);
        }

        this->length++;
        return true;
    }    
    // End: material

    // Begin: compression
    void Compression::__init_candidates() {
        for (std::string& type : this->types) {
            if (type == "exponent") {
                std::vector<Model*> candidate;
                candidate.push_back(new ExpoFunction());
                this->candidates.push_back(candidate);
            }
            else if (type == "linear") {
                std::vector<Model*> candidate;
                candidate.push_back(new LinearFunction());
                this->candidates.push_back(candidate);
            }
            else if (type == "quadratic") {
                std::vector<Model*> candidate;
                candidate.push_back(new PolyFunction(2, "quadratic"));
                this->candidates.push_back(candidate);
            }
            else if (type == "cubic") {
                std::vector<Model*> candidate;
                candidate.push_back(new PolyFunction(3, "cubic"));
                this->candidates.push_back(candidate);
            }
        }
    }

    void Compression::__choose_best_model() {
        double min_param = INFINITY;

        for (int i=0; i<this->candidates.size(); i++) {
            double num_param = 0;
            for (Model* model : this->candidates[i]) {
                if (model->complete) num_param += model->getDim();
                else if (model->getLen() > 2) num_param += model->getDim() * this->penalty; 
            }

            if (num_param < min_param) { 
                min_param = num_param; 
                this->best_model = i; 
            }
        }
    }

    void Compression::initialize(int count, char** params) {
        this->error = atof(params[0]);
        this->penalty = atof(params[1]);
        for (int i=2; i<count; i++) 
            this->types.push_back(std::string(params[i]));

        this->__init_candidates();
    }

    BinObj* Compression::serialize() {
        BinObj* obj = new BinObj;
        std::vector<Model*>& candidate = this->candidates[this->best_model];

        if (this->first) {
            this->first = false;
            obj->put((float) candidate[0]->subs(0));
        }
        
        for (Model* model : candidate) {
            if (model->complete) {
                double* coeffs = model->getCoeffs();
                int embedded = model->getLen() << 2;

                if (model->getType() == "exponent") embedded |= 0; 
                else if (model->getType() == "linear") embedded |= 1;
                else if (model->getType() == "quadratic") embedded |= 2;
                else if (model->getType() == "cubic") embedded |= 3;
                
                obj->put(VariableByteEncoding::encode(embedded));
                for (int i = 0; i < model->getDim(); i++) {
                    obj->put((float) coeffs[i]);
                }
            }
        }

        return obj;
    }

    void Compression::finalize() {
        for (std::vector<Model*>& candidate : this->candidates) {
            for (Model* model : candidate) model->complete = true;
        }

        this->__choose_best_model();
        this->yield();

        for (std::vector<Model*>& candidate : this->candidates) {
            for (Model* model : candidate) delete model;
        }
    }

    void Compression::compress(Univariate* data) {
        if (this->Fcf == false) {
            bool flag = true;
            for (std::vector<Model*>& candidate : this->candidates) {
                Model* model = candidate.back();
                if (!model->approximate(this->error, data->get_value())) {
                    if (model->getType() == "exponent") candidate.push_back(new ExpoFunction());
                    else if (model->getType() == "linear") candidate.push_back(new LinearFunction());
                    else if (model->getType() == "quadratic") candidate.push_back(new PolyFunction(2, "quadratic"));
                    else if (model->getType() == "cubic") candidate.push_back(new PolyFunction(3, "cubic"));

                    candidate.back()->approximate(this->error, model->subs(model->getLen() - 1));       
                    candidate.back()->approximate(this->error, data->get_value());
                }

                if (candidate.size() == 1) flag = false;
            }
            
            if (flag) {
                this->__choose_best_model();
                Model* model = this->candidates[this->best_model].back();
                if (model->getLen() == 2) this->Frv = true;
                else this->Fcf = true; 
            }
        }
        else {
            Model* model = this->candidates[this->best_model].back();
            if (!model->approximate(this->error, data->get_value())) {
                this->Fcf = false;
                this->Frv = true;
            }
        }

        if (this->Frv) {
            this->yield();
            Model* model = this->candidates[this->best_model].back();
            double last = model->complete
                ? model->subs(model->getLen() - 1) : model->subs(0);
            
            for (std::vector<Model*>& candidate : this->candidates) {
                for (Model* model : candidate) delete model;
            }
            this->candidates.clear();
            this->__init_candidates();

            for (std::vector<Model*>& candidate : this->candidates) {
                candidate.back()->approximate(this->error, last);
                candidate.back()->approximate(this->error, data->get_value());
            }

            this->best_model = -1;
            this->Fcf = false;
            this->Frv = false;
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

        int start = 1;
        if (this->pivot == INFINITY) {
            start = 0;
            this->pivot = compress_data->getFloat();
        }

        long embedded = VariableByteEncoding::decode(compress_data);
        int type = (embedded & 3);
        int length = embedded >> 2;

        if (type == 0) {
            double growth = compress_data->getFloat();
            double scale = this->pivot;

            for (int i = start; i < length; i++) {
                if (base_obj == nullptr) {
                    base_obj = new CSVObj;
                    base_obj->pushData(std::to_string(this->basetime + i * this->interval));
                    base_obj->pushData(std::to_string(scale * exp(growth * i)));

                    prev_obj = base_obj;
                }
                else {
                    CSVObj* obj = new CSVObj;
                    obj->pushData(std::to_string(this->basetime + i * this->interval));
                    obj->pushData(std::to_string(scale * exp(growth * i)));

                    prev_obj->setNext(obj);
                    prev_obj = obj;
                }
            }
            this->pivot = scale * exp(growth * (length - 1));
        }
        else if (type == 1) {
            double slope = compress_data->getFloat();
            double intercept = this->pivot;
            
            for (int i = start; i < length; i++) {
                if (base_obj == nullptr) {
                    base_obj = new CSVObj;
                    base_obj->pushData(std::to_string(this->basetime + i * this->interval));
                    base_obj->pushData(std::to_string(slope * i + intercept));

                    prev_obj = base_obj;
                }
                else {
                    CSVObj* obj = new CSVObj;
                    obj->pushData(std::to_string(this->basetime + i * this->interval));
                    obj->pushData(std::to_string(slope * i + intercept));

                    prev_obj->setNext(obj);
                    prev_obj = obj;
                }
            }
            this->pivot = slope * (length - 1) + intercept;
        }
        else if (type == 2) {
            float* coefficients = new float[3];
            for (int i = 0; i < 2; i++) {
                coefficients[i] = compress_data->getFloat();
            }
            coefficients[2] = this->pivot;
            Polynomial polynomial(2, coefficients, true);
            delete[] coefficients;

            for (int i = start; i < length; i++) {
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
            this->pivot = polynomial.subs(length - 1);
        }
        else if (type == 3) {
            float* coefficients = new float[4];
            for (int i = 0; i < 3; i++) {
                coefficients[i] = compress_data->getFloat();
            }
            coefficients[3] = this->pivot;
            Polynomial polynomial(3, coefficients, true);
            delete[] coefficients;

            for (int i = start; i < length; i++) {
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
            this->pivot = polynomial.subs(length - 1);
        }
        
        this->basetime += (length - 1) * this->interval;
        return base_obj;
    }
    // End: decompression
}