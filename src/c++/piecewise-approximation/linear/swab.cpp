#include "piecewise-approximation/linear.hpp"

namespace Swab {
    Segment::Segment(std::vector<double>& data, UpperHull& l_cvx, LowerHull& u_cvx) {
        this->data = data;
        this->u_cvx = u_cvx;
        this->l_cvx = l_cvx;
    }

    Line Approximator::interpolate(std::vector<double>& data) {
        Point2D first(0, data.front());
        Point2D last(data.size()-1, data.back());

        return Line::line(first, last);
    }

    Line Approximator::regression(int size, double acc, double avg_x, double avg_y, long square) {
        double variance = (double) square / size - avg_x * avg_x;
        double covariance = (double) acc / size - avg_x * avg_y;
        double slope = covariance / variance;
        double intercept = avg_y - slope * avg_x;

        return Line(slope, intercept);
    }

    Line Approximator::regression(std::vector<double>& data) {
        double avg_x = 0;
        double avg_y = 0;
        double acc = 0;
        unsigned long square = 0;

        for (unsigned long i = 0; i < data.size(); i++) {
            avg_x += (double)i;
            avg_y += (double)data[i];
            acc += (double)i * data[i];
            square += i * i;
        }
        avg_x /= data.size(); avg_y /= data.size();

        double variance = (double) square / data.size() - avg_x * avg_x;
        double covariance = (double) acc / data.size() - avg_x * avg_y;
        double slope = covariance / variance;
        double intercept = avg_y - slope * avg_x;

        return Line(slope, intercept);
    }

    // Mean square error calculate
    double Approximator::cal_error(std::vector<double>& data, Line& line) {
        double error = 0;
        for (int i=0; i<data.size(); i++) {
            error += std::pow(data[i]-line.subs(i), 2);
        }

        return error / data.size();
    }

    // Verify new line satisfies infinity bound or not
    bool Grouper::bound_check(UpperHull& l_cvx, LowerHull& u_cvx, Line& line, int offset) {
        // Immediately terminate when individual error exceed the allowable threshold
        for (int i=0; i<l_cvx.size(); i++) {
            Point2D p = l_cvx.at(i);
            if (line.subs(p.x + offset) - p.y < 0) return false;
        }

        for (int i=0; i<u_cvx.size(); i++) {
            Point2D p = u_cvx.at(i);
            if (p.y - line.subs(p.x + offset) < 0) return false;
        }

        return true;
    }

    void Grouper::merge(Segment& s1, Segment& s2, std::string mode) {
        int offset = mode == "interpolate" ? 1 : 0;
        s1.u_cvx.concat(s2.u_cvx, s1.data.size() - offset);
        s1.l_cvx.concat(s2.l_cvx, s1.data.size() - offset);
        s1.data.insert(s1.data.end(), s2.data.begin() + offset, s2.data.end());
    }

    double Grouper::merge_cost(Segment& s1, Segment& s2, std::string mode) {
        int offset = mode == "interpolate" ? 1 : 0;
        std::vector<double> s = s1.data;
        s.insert(s.end(), s2.data.begin() + offset, s2.data.end());
        Line line = mode == "interpolate" ?
            Approximator::interpolate(s) : Approximator::regression(s);
    
        if (!Grouper::bound_check(s1.l_cvx, s1.u_cvx, line)) return INFINITY;
        else if (!Grouper::bound_check(s2.l_cvx, s2.u_cvx, line, s1.data.size()-offset)) return INFINITY;
        else return Approximator::cal_error(s, line);
    }

    // Begin: compression
    bool Compression::__sliding_window() {
        Point2D p(this->window.size()-1, this->window.back());

        if (this->mode == "interpolate") {
            if (this->window.size() > 2) {
                Line line = Approximator::interpolate(this->window);
                
                return Grouper::bound_check(this->l_cvx, this->u_cvx, line) && 
                    std::abs(line.subs(p.x) - p.y) <= this->error;
            }
        }
        else if (this->mode == "regression") {
            this->accumulate += p.x*p.y;
            this->accumulate_square += p.x*p.x;
            this->average_x = (this->average_x * (this->window.size() - 1) + p.x) / this->window.size();
            this->average_y = (this->average_y * (this->window.size() - 1) + p.y) / this->window.size();
            
            if (this->window.size() > 2) {
                Line line = Approximator::regression(
                    this->window.size(), this->accumulate, this->average_x, 
                    this->average_y, this->accumulate_square
                );
                
                return Grouper::bound_check(this->l_cvx, this->u_cvx, line) && 
                    std::abs(line.subs(p.x) - p.y) <= this->error;
            }
        }
        
        return true;
    }

    void Compression::__bottom_up() {
        for (int i=this->m_err.size(); i<this->segments.size()-1; i++) {
            this->m_err.push_back(Grouper::merge_cost(
                this->segments[i], this->segments[i+1], this->mode
            ));
        }

        auto it = std::min_element(this->m_err.begin(), this->m_err.end());
        while (*it != INFINITY && this->segments.size() > 1) {
            int index = std::distance(this->m_err.begin(), it);
            Grouper::merge(this->segments[index], this->segments[index+1], this->mode);

            this->m_err.erase(this->m_err.begin() + index);
            this->segments.erase(this->segments.begin() + index + 1);
            
            if (index != 0) {
                this->m_err[index-1] = Grouper::merge_cost(
                    this->segments[index-1], this->segments[index], this->mode
                );
            }
            if (index < this->m_err.size()) {
                this->m_err[index] = Grouper::merge_cost(
                    this->segments[index], this->segments[index+1], this->mode
                );
            }
            it = std::min_element(this->m_err.begin(), this->m_err.end());  
        }
    }

    void Compression::initialize(int count, char** params) {
        this->error = atof(params[0]);
        this->mode = params[1];
        this->n_segment = atoi(params[2]);
    }

    void Compression::finalize() {
        if (this->window.size() >= 2) {
            this->segments.push_back(Segment(
                this->window, this->l_cvx, this->u_cvx
            ));
        }

        this->__bottom_up();
        while (this->segments.size() != 0) {
            this->yield();
            this->segments.erase(this->segments.begin());
        }
        this->m_err.clear();
    }

    void Compression::compress(Univariate* data) {
        this->window.push_back(data->get_value());
        
        if (!this->__sliding_window()) {
            this->window.pop_back();
            this->segments.push_back(Segment(
                this->window, this->l_cvx, this->u_cvx
            ));
            
            if (this->mode == "interpolate") {
                // Add last data point of previous segment to ensure connectivity
                double tail = this->window.back();
                this->window.clear();
                this->window.push_back(tail);
                this->window.push_back(data->get_value());

                this->l_cvx.clear();
                this->l_cvx.append(Point2D(0, tail - this->error));
                this->u_cvx.clear();
                this->u_cvx.append(Point2D(0, tail + this->error));
            }
            else if (this->mode == "regression") {
                this->window.clear();
                this->window.push_back(data->get_value());
                this->l_cvx.clear();
                this->u_cvx.clear();

                this->average_x = 0;
                this->average_y = data->get_value();
                this->accumulate = 0;
                this->accumulate_square = 0;
            }
        }

        // Perform bottom up once n_segment sliding window segments have accumulated
        if (this->segments.size() == this->n_segment) {
            this->__bottom_up();
            this->yield();
            this->segments.erase(this->segments.begin());
            if (!this->m_err.empty()) this->m_err.erase(this->m_err.begin());
        }

        this->l_cvx.append(Point2D(this->window.size()-1, data->get_value() - this->error));
        this->u_cvx.append(Point2D(this->window.size()-1, data->get_value() + this->error));
    }

    BinObj* Compression::serialize() {
        BinObj* obj = new BinObj;
        Segment& segment = this->segments[0];

        if (this->mode == "interpolate") {
            if (this->first) {
                obj->put((float) segment.data.front());
                this->first = false;
            }
            
            obj->put(VariableByteEncoding::encode(segment.data.size()));
            obj->put((float) segment.data.back());
        }
        else if (this->mode == "regression") {
            Line line = Approximator::regression(segment.data);
            obj->put(VariableByteEncoding::encode(segment.data.size()));
            obj->put((float) line.get_intercept());
            obj->put((float) line.get_slope());
        }
        
        return obj;
    }
    // End: compression

    // Begin: decompression
    void Decompression::initialize(int count, char** params) {
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
        Line line(0, 0);
        int offset = this->mode == "interpolate" ? 1 : 0;

        if (this->mode == "interpolate") {
            start = 1;
            if (this->pivot == INFINITY) {
                start = 0;
                this->pivot = compress_data->getFloat();
            }

            length = VariableByteEncoding::decode(compress_data);
            float last = compress_data->getFloat();
            line = Line::line(Point2D(0, this->pivot), Point2D(length-1, last)); 
            this->pivot = last;
        }
        else if (this->mode == "regression") {
            length = VariableByteEncoding::decode(compress_data);
            float intercept = compress_data->getFloat();
            float slope = compress_data->getFloat();
            line = Line(slope, intercept);
        }

        for (int i=start; i<length; i++) {
            if (base_obj == nullptr) {
                base_obj = new CSVObj;
                base_obj->pushData(std::to_string(this->basetime + i * this->interval));
                base_obj->pushData(std::to_string(line.subs(i)));

                prev_obj = base_obj;
            }
            else {
                CSVObj* obj = new CSVObj;
                obj->pushData(std::to_string(this->basetime + i * this->interval));
                obj->pushData(std::to_string(line.subs(i)));

                prev_obj->setNext(obj);
                prev_obj = obj;
            }
        }
        
        this->basetime += (length - offset) * this->interval;
        return base_obj;
    }
    // End: decompression
};