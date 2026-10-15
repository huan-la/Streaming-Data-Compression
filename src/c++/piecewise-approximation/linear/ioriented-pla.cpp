#include "piecewise-approximation/linear.hpp"

namespace IOrientedPLA {
    SDLP sdlp;

    bool Compression::__feasible_cone(BinObj* obj) {
        if (std::abs(this->u_line->get_slope()-this->l_line->get_slope()) < EPS) return false;

        Point2D p = Line::intersection(*this->u_line, *this->l_line);
        long u_left = static_cast<long>(std::ceil(this->u_line->get_intercept()));
        long l_left = static_cast<long>(std::floor(this->l_line->get_intercept()));
        long u_right = static_cast<long>(std::floor(this->u_line->subs(this->length)));
        long l_right = static_cast<long>(std::ceil(this->l_line->subs(this->length)));
        
        if (u_right >= l_right) {
            unsigned long value = ZigZagEncoding::encode(u_right);
            Line line = Line::line(p, Point2D(this->length, u_right));

            long embedded = this->length << 3 | 1;
            obj->put(VariableByteEncoding::encode(embedded));
            obj->put((float) line.get_slope());
            obj->put(VariableByteEncoding::encode(value));
        
            return true;
        }
        else if (l_left >= u_left) {
            unsigned long value = ZigZagEncoding::encode(u_left);
            Line line = Line::line(p, Point2D(0, u_left));

            long embedded = this->length << 3 | 2;
            obj->put(VariableByteEncoding::encode(embedded));
            obj->put((float) line.get_slope());
            obj->put(VariableByteEncoding::encode(value));

            return true;
        }   

        long root = this->length;
        if (std::abs(this->u_line->get_slope()) < EPS) {
            if (this->l_line->get_root() > p.x) {
                root = static_cast<long>(std::ceil(this->l_line->get_root()));
                if (std::abs(root - p.x) < EPS) root++;
            }
            else {
                root = static_cast<long>(std::floor(this->l_line->get_root()));
                if (std::abs(root - p.x) < EPS) root--;
            }
        }
        else if (std::abs(this->l_line->get_slope()) < EPS) {
            if (this->u_line->get_root() > p.x) {
                root = static_cast<long>(std::ceil(this->u_line->get_root()));
                if (std::abs(root - p.x) < EPS) root++;
            }
            else {
                root = static_cast<long>(std::floor(this->u_line->get_root()));
                if (std::abs(root - p.x) < EPS) root--;
            }
        }
        else {
            // u_root: root that higher than l_root
            // Not the root of u_line
            double u_root = this->u_line->get_root() > this->l_line->get_root() ? this->u_line->get_root() : this->l_line->get_root();
            double l_root = this->u_line->get_root() <= this->l_line->get_root() ? this->u_line->get_root() : this->l_line->get_root();

            if (this->u_line->get_slope() * this->l_line->get_slope() < 0) {
                root = static_cast<long>(std::floor(l_root));
                if (std::abs(root - p.x) < EPS) root--;
            }
            else {
                u_root = static_cast<long>(std::floor(u_root));
                l_root = static_cast<long>(std::ceil(l_root));
                
                if (u_root >= l_root) {
                    root = l_root;
                    if (std::abs(root - p.x) < EPS) {
                        if (root + 1 <= u_root) root++;
                        else root = this->length;
                    }
                }
            }
        }

        if (root != this->length) {
            Point2D p2(root, 0);
            Line line = Line::line(p2, p);
            unsigned long value = ZigZagEncoding::encode(root);
            long embedded = this->length << 3 | 3;

            obj->put(VariableByteEncoding::encode(embedded));
            obj->put(VariableByteEncoding::encode(value));
            obj->put((float) line.subs(this->length));

            return true;
        }

        return false;
    }

    bool Compression::__feasible_rectangle(BinObj* obj) {
        double g_dis = -INFINITY;
        Line u_bound(INFINITY, INFINITY);
        Line l_bound(INFINITY, INFINITY);

        // Using proposed algorithm
        auto solve_min_dis = [&](auto& A, auto& B, bool upper) {
            double l_dis = -INFINITY;

            for (int i = 0, end = B.size() - 1; i < A.size() - 1; ++i) {
                Line line = Line::line(A.at(i), A.at(i + 1));

                double m = line.get_slope();
                double b = line.get_intercept();
                double best = upper ? INFINITY : -INFINITY;

                for (int j = end; j >= 0; --j) {
                    double x = B.at(j).y - m * B.at(j).x;

                    if ((upper && x < best) || (!upper && x > best))
                        best = x;
                    else {
                        end = j + 1;
                        break;
                    }
                }

                double d = upper ? best - b : b - best;
                if (d <= l_dis) break;

                l_dis = d;
                if (l_dis > g_dis) {
                    g_dis = l_dis;
                    (upper ? l_bound : u_bound).set(m, b);
                    (upper ? u_bound : l_bound).set(m, best);
                }
            }
        };

        solve_min_dis(this->u_cvx, this->l_cvx, true);
        solve_min_dis(this->l_cvx, this->u_cvx, false);

        // // Using Seidel LP
        // Eigen::VectorXd x(3); 
        // Eigen::VectorXd c(3);                    
        // Eigen::MatrixXd A(this->u_cvx.size() + this->l_cvx.size(), 3); 
        // Eigen::VectorXd b(this->u_cvx.size() + this->l_cvx.size());

        // c(0) = 0; c(1) = -1.0; c(2) = 1.0;
        // for (int i=0; i<this->u_cvx.size(); i++) {
        //     Point2D& p = this->u_cvx.at(i);
        //     A(i, 0) = -p.x; A(i, 1) = 0; A(i, 2) = -1.0; b(i) = -p.y;
        // }
        // for (int i=0; i<this->l_cvx.size(); i++) {
        //     Point2D& p = this->l_cvx.at(i);
        //     int size = this->u_cvx.size() + i;
        //     A(size, 0) = p.x; A(size, 1) = 1.0; A(size, 2) = 0; b(size) = p.y;
        // }

        // double minobj = sdlp.linprog(c, A, b, x);
        // if (minobj == INFINITY || minobj == -INFINITY) {
        //     return false;
        // }
        // else {
        //     u_bound.set(x(0), x(1));
        //     l_bound.set(x(0), x(2));
        // }

        long u_right = static_cast<long>(std::floor(u_bound.subs(this->length-1)));
        long l_right = static_cast<long>(std::ceil(l_bound.subs(this->length-1)));
        long u_left = static_cast<long>(std::floor(u_bound.get_intercept()));
        long l_left = static_cast<long>(std::ceil(l_bound.get_intercept()));
        
        if (u_left >= l_left && u_right >= l_right) {
            unsigned long value_1 = ZigZagEncoding::encode(u_left);
            unsigned long value_2 = ZigZagEncoding::encode(u_right);
            
            if (value_1 < 15 && value_2 < 15) {
                // Elias Gamma Encoding is more efficient
                obj->put(VariableByteEncoding::encode(this->length << 3 | 7));
                EliasGammaEncoding::encode(value_1 + 1, obj);
                EliasGammaEncoding::encode(value_2 + 1, obj);
            }
            else {
                // Variable Byte Encoding is more efficient
                obj->put(VariableByteEncoding::encode(this->length << 3 | 4));
                obj->put(VariableByteEncoding::encode(value_1));
                obj->put(VariableByteEncoding::encode(value_2));
            }

            return true;
        }
        else if (u_right >= l_right) {
            double u_root = u_bound.get_root() > l_bound.get_root() ? u_bound.get_root() : l_bound.get_root();
            double l_root = u_bound.get_root() <= l_bound.get_root() ? u_bound.get_root() : l_bound.get_root();
            long i_u_root = static_cast<long>(std::floor(u_root));
            long i_l_root = static_cast<long>(std::ceil(l_root));

            if (i_u_root >= i_l_root && i_l_root < 0) {
                unsigned long value = ZigZagEncoding::encode(u_right);
                unsigned long z_root = -i_l_root;
                long embedded = this->length << 3 | 5;

                obj->put(VariableByteEncoding::encode(embedded));
                obj->put(VariableByteEncoding::encode(z_root));
                obj->put(VariableByteEncoding::encode(value));
                
                return true;
            }       
        }
        else if (u_left >= l_left) {
            double u_root = u_bound.get_root() > l_bound.get_root() ? u_bound.get_root() : l_bound.get_root();
            double l_root = u_bound.get_root() <= l_bound.get_root() ? u_bound.get_root() : l_bound.get_root();
            long i_u_root = static_cast<long>(std::floor(u_root));
            long i_l_root = static_cast<long>(std::ceil(l_root));

            if (i_u_root >= i_l_root && i_u_root >= this->length) {
                unsigned long value = ZigZagEncoding::encode(u_left);
                unsigned long z_root = i_u_root;
                long embedded = this->length << 3 | 6;

                obj->put(VariableByteEncoding::encode(embedded));
                obj->put(VariableByteEncoding::encode(z_root));
                obj->put(VariableByteEncoding::encode(value));
                
                return true;
            }
        }

        return false;
    }

    // Begin: compression
    void Compression::initialize(int count, char** params) {
        this->error = atof(params[0]);
    }

    void Compression::finalize() {
        if (this->length >= 2) this->yield();
        
        if (this->pivot != nullptr) delete this->pivot;
        if (this->u_line != nullptr) delete this->u_line;
        if (this->l_line != nullptr) delete this->l_line;
        this->u_cvx.clear();
        this->l_cvx.clear();
    }

    BinObj* Compression::serialize() {
        BinObj* obj = new BinObj;

        // if (this->__feasible_rectangle(obj)) { return obj; }
        if (this->__feasible_cone(obj)) { return obj; }
        else {
            float slope = (this->u_line->get_slope() + this->l_line->get_slope()) / 2;
            float intercept = (this->u_line->get_intercept() + this->l_line->get_intercept()) / 2;

            long embedded = this->length << 3 | 0;
            obj->put(VariableByteEncoding::encode(embedded));
            obj->put((float) slope);
            obj->put((float) intercept);
        }

        return obj;
    }
    
    void Compression::compress(Univariate* data) {
        Point2D p(this->length++, data->get_value());

        if (this->pivot == nullptr) {
            this->pivot = new Point2D(p.x, p.y);
            this->u_cvx.append(Point2D(p.x, p.y - this->error));
            this->l_cvx.append(Point2D(p.x, p.y + this->error));
        }
        else if (this->u_line == nullptr) {
            Line u_line = Line::line(
                Point2D(this->pivot->x, this->pivot->y-this->error), 
                Point2D(p.x, p.y+this->error)
            );
            Line l_line = Line::line(
                Point2D(this->pivot->x, this->pivot->y+this->error), 
                Point2D(p.x, p.y-this->error)
            );

            this->u_line = new Line(u_line.get_slope(), u_line.get_intercept());
            this->l_line = new Line(l_line.get_slope(), l_line.get_intercept());
            this->u_cvx.append(Point2D(p.x, p.y - this->error));
            this->l_cvx.append(Point2D(p.x, p.y + this->error));
        }
        else {
            if (this->l_line->subs(p.x) - p.y - this->error > EPS || p.y - this->error - this->u_line->subs(p.x) > EPS) {
                this->length--;
                this->yield();
                
                delete this->u_line; this->u_line = nullptr;
                delete this->l_line; this->l_line = nullptr;
                this->u_cvx.clear(); this->l_cvx.clear();

                this->pivot->x = 0; this->pivot->y = p.y;
                u_cvx.append(Point2D(0, p.y - this->error));
                l_cvx.append(Point2D(0, p.y + this->error));
                this->length = 1;
            }
            else {
                bool update_u = p.y + this->error < this->u_line->subs(p.x);
                bool update_l = p.y - this->error > this->l_line->subs(p.x);
                
                if (update_u) {
                    int index = -1;
                    double min_slp = INFINITY;

                    for (int i=0; i<this->u_cvx.size(); i++) {
                        Line line = Line::line(this->u_cvx.at(i), Point2D(p.x, p.y + this->error));
                        if (line.get_slope() < min_slp) {
                            index = i;
                            min_slp = line.get_slope();
                            this->u_line->set(line.get_slope(), line.get_intercept());
                        }
                    }
                    this->u_cvx.erase_from_begin(index);
                }
                if (update_l) {
                    int index = -1;
                    double max_slp = -INFINITY;

                    for (int i=0; i<this->l_cvx.size(); i++) {
                        Line line = Line::line(this->l_cvx.at(i), Point2D(p.x, p.y - this->error));
                        if (line.get_slope() > max_slp) {
                            index = i;
                            max_slp = line.get_slope();
                            this->l_line->set(line.get_slope(), line.get_intercept());
                        }
                    }
                    this->l_cvx.erase_from_begin(index);
                }

                if (update_u) this->l_cvx.append(Point2D(p.x, p.y + this->error));
                if (update_l) this->u_cvx.append(Point2D(p.x, p.y - this->error));
            }
        }
    }
    // End: compression

    // Begin: decompression
    void Decompression::initialize(int count, char** params) {

    }

    void Decompression::finalize() {
        // Do nothing
    }
    
    CSVObj* Decompression::decompress(BinObj* compress_data) {
        CSVObj* base_obj = nullptr;
        CSVObj* prev_obj = nullptr;

        unsigned long embedded = VariableByteEncoding::decode(compress_data);
        int flag = embedded & 0b111L;
        int length = embedded >> 3;
        float slp = 0;
        float intercept = 0;

        if (flag == 0) {
            slp = compress_data->getFloat();
            intercept = compress_data->getFloat();
        }
        else if (flag == 1) {
            slp = compress_data->getFloat();
            long value = ZigZagEncoding::decode(VariableByteEncoding::decode(compress_data));
            Line line = Line::line(slp, Point2D(length, value));
            intercept = line.get_intercept();
        }
        else if (flag == 2) {
            slp = compress_data->getFloat();
            long value = ZigZagEncoding::decode(VariableByteEncoding::decode(compress_data));
            Line line = Line::line(slp, Point2D(0, value));
            intercept = line.get_intercept();
        }
        else if (flag == 3) {
            long root = ZigZagEncoding::decode(VariableByteEncoding::decode(compress_data));
            double value = compress_data->getFloat();
            Line line = Line::line(Point2D(root, 0), Point2D(length, value));
            
            slp = line.get_slope();
            intercept = line.get_intercept();
        }
        else if (flag == 4) {
            long value_1 = ZigZagEncoding::decode(VariableByteEncoding::decode(compress_data));
            long value_2 = ZigZagEncoding::decode(VariableByteEncoding::decode(compress_data));
            Line line = Line::line(Point2D(0, (double) value_1), Point2D(length-1, (double) value_2));

            slp = line.get_slope();
            intercept = line.get_intercept();
        }
        else if (flag == 5) {
            long root = VariableByteEncoding::decode(compress_data);
            long value = ZigZagEncoding::decode(VariableByteEncoding::decode(compress_data));

            Line line = Line::line(Point2D(-root, 0), Point2D(length-1, (double) value));
            slp = line.get_slope();
            intercept = line.get_intercept();
        }
        else if (flag == 6) {
            long root = VariableByteEncoding::decode(compress_data);
            long value = ZigZagEncoding::decode(VariableByteEncoding::decode(compress_data));
            
            Line line = Line::line(Point2D(root, 0), Point2D(0, (double) value));
            slp = line.get_slope();
            intercept = line.get_intercept();
        }
        else if (flag == 7) {
            long value_1 = ZigZagEncoding::decode(EliasGammaEncoding::decode(compress_data) - 1);
            long value_2 = ZigZagEncoding::decode(EliasGammaEncoding::decode(compress_data) - 1);
            Line line = Line::line(Point2D(0, (double) value_1), Point2D(length-1, (double) value_2));
        
            slp = line.get_slope();
            intercept = line.get_intercept();
            compress_data->flushBits();
        }

        for (int i=0; i<length; i++) {
            if (base_obj == nullptr) {
                base_obj = new CSVObj;
                base_obj->pushData(std::to_string(this->basetime + i * this->interval));
                base_obj->pushData(std::to_string(slp * i + intercept));

                prev_obj = base_obj;
            }
            else {
                CSVObj* obj = new CSVObj;
                obj->pushData(std::to_string(this->basetime + i * this->interval));
                obj->pushData(std::to_string(slp * i + intercept));

                prev_obj->setNext(obj);
                prev_obj = obj;
            }
        }

        this->basetime += length * this->interval;
        return base_obj;
    }
    // End: decompression
};
