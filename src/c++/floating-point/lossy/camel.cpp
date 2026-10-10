#include "floating-point/lossy.hpp"

namespace Camel {
    // Begin: compression
    void Compression::initialize(int count, char** params) {
        this->buffer_size = atoi(params[0]);
        this->decimal_count = atoi(params[1]);

        this->bitstream = new BinObj;
    }

    void Compression::finalize() {
        this->__compress(Camel::END_SIGNAL);
        this->yield();
        
        delete this->bitstream;
        this->bitstream = nullptr;
        std::cout << "size: " << this->size << "\n";
    }

    void Compression::__compress_integer(int int_value, int intSignal) {
        int diff = int_value - this->storedVal;
        this->bitstream->put(intSignal, 1);
        this->size += 1;

        if (diff >= -1 && diff <= 1) {
            this->bitstream->put((diff + 1), 2); 
            this->size += 2;
        } 
        else{
            this->bitstream->put(3, 2);
            this->bitstream->put(diff >= 0, 1);
            diff = diff < 0 ? -diff : diff;
            diff = diff > Camel::END_DIFF ? Camel::END_DIFF : diff;
            
            if (diff >=2 && diff < 8) { 
                this->bitstream->put(0, 1);
                this->bitstream->put(diff, 3);
                this->size += 7;
            } 
            else {
                this->bitstream->put(1, 1); 
                this->bitstream->put(diff, 16);
                this->size += 20;
            }
        }

        this->storedVal = int_value;
    }

    void Compression::__compress_decimal(int decimal_value, int decimal_count) {
        if (decimal_count != 0) {
            this->bitstream->put(decimal_count-1, 2);
            this->size += 2;

            int threshold = Camel::THRESHOLDS[decimal_count-1];
            int m = decimal_value;
    
            if (decimal_value - threshold >= 0) {  
                this->bitstream->put(true, 1);
                m = decimal_value % threshold; 
                uint32_t xored = std::bit_cast<uint32_t>(((float) m / Camel::POWERS[decimal_count]+1)) ^ 
                    std::bit_cast<uint32_t>((float) decimal_value / Camel::POWERS[decimal_count]+1);

                this->bitstream->put(xored >> 23 - decimal_count, decimal_count);
                this->size += 1 + decimal_count;
            } 
            else {  
                this->bitstream->put(false, 1);
                this->size += 1;
            }
            
            if (decimal_count == 1) { 
                this->bitstream->put(m, 3);
                this->size += 3;
            } 
            else if (decimal_count == 2) {
                if (m < 8) {
                    this->bitstream->put(0, 1);
                    this->bitstream->put(m, 3);
                    this->size += 4;
                }  
                else {
                    this->bitstream->put(1, 1);
                    this->bitstream->put(m, 5);
                    this->size += 6;
                }
            } 
            else if (decimal_count == 3) {
                if (m < 2) {
                    this->bitstream->put(0, 2);
                    this->bitstream->put(m, 1);
                    this->size += 3;
                }
                else if (m < 8){
                    this->bitstream->put(1, 2);
                    this->bitstream->put(m, 3);
                    this->size += 5;
                }
                else if (m < 32) {
                    this->bitstream->put(2, 2);
                    this->bitstream->put(m, 5);
                    this->size += 7;
                }
                else {
                    this->bitstream->put(3, 2);
                    this->bitstream->put(m, Camel::MVALUEBITS[decimal_count-1]);
                    this->size += 2;
                    this->size += Camel::MVALUEBITS[decimal_count-1];
                }
            } 
            else {
                if (m < 16) {
                    this->bitstream->put(0, 2);
                    this->bitstream->put(m, 4);
                    this->size += 6;
                }
                else if (m < 64){
                    this->bitstream->put(1, 2);
                    this->bitstream->put(m, 6);
                    this->size += 8;
                }
                else if (m < 256) {
                    this->bitstream->put(2, 2);
                    this->bitstream->put(m, 8);
                    this->size += 10;
                }
                else {
                    this->bitstream->put(3, 2);
                    this->bitstream->put(m, Camel::MVALUEBITS[decimal_count-1]);
                    this->size += 2;
                    this->size += Camel::MVALUEBITS[decimal_count-1];
                }
            }
        }
    }

    std::pair<int, int> Compression::__decimal_count(float value) {
        float factor = 1;
        int decimal_count = 0;
        int decimal_value = 0;

        value = std::abs(value);
        while (std::abs(value * factor - std::round(value * factor)) > EPS) {
            factor *= 10.0;
            decimal_count++;
        }
        
        if (decimal_count == 0) decimal_count = 1;
        if (decimal_count > 0 && decimal_count <= this->decimal_count) 
            decimal_value = (int) std::round(value * Camel::POWERS[decimal_count]) % Camel::POWERS[decimal_count];
        else {
            decimal_value = (int) std::round(value * Camel::POWERS[this->decimal_count]) % Camel::POWERS[this->decimal_count];
            decimal_count = this->decimal_count;
        }
         
        return std::make_pair(std::abs(decimal_value), decimal_count);
    }

    void Compression::__compress(float value) {
        // Dirty trick to deal with rounding problem
        value = std::round(value * Camel::POWERS[this->decimal_count]) / Camel::POWERS[this->decimal_count];
        
        if (this->first) {
            this->size += 32;
            this->first = false;
            this->bitstream->put(std::bit_cast<uint32_t>(value), 32);
            this->storedVal = (int) value;
        }
        else if (value == Camel::END_SIGNAL) {
            // Dirty trick: end signal always have the maximum integer difference.
            // And it shoule be divided by 2 to prevent numerical overflow
            this->__compress_integer(INT_MAX / 2, 1);
        }
        else {
            int intSignal = value < 0 ? 0 : 1;
            std::pair<int, int> decimal = this->__decimal_count(value);

            this->__compress_integer((int) value, intSignal);
            this->__compress_decimal(decimal.first, decimal.second);
        }
    }

    void Compression::compress(Univariate* data) {
        float value = data->get_value();
        this->__compress(value);

        if (++this->count == this->buffer_size) {
            this->yield();
            this->count = 0;
        }
    }

    BinObj* Compression::serialize() {
        BinObj* obj = this->bitstream;
        this->bitstream = new BinObj;
        
        return obj;
    }
    // End: compression

    // Begin: decompression
    void Decompression::initialize(int count, char** params) {
        this->buffer_size = atoi(params[0]);
    }

    void Decompression::finalize() {
        // Do nothing
    }

    std::pair<int, int> Decompression::__decompress_integer(BinObj* compress_data){
        int intSignal = compress_data->getBits(1);
        int integerNum = compress_data->getBits(2);

        int diffVal;
        if (integerNum < 3) {
            diffVal = integerNum - 1;
        } 
        else {
            int diffSymbol = compress_data->getBits(1);
            if (compress_data->getBits(1) == 0) diffVal = compress_data->getBits(3);
            else diffVal = compress_data->getBits(16);

            diffVal = (diffSymbol == 0 ? -1: 1) * diffVal;
        }
        if (std::abs(diffVal) != Camel::END_DIFF) {
            this->storedVal += diffVal;

            return std::make_pair(
                (this->storedVal == 0 && intSignal == 0) ? 0 : (this->storedVal >= 0 ? 1 : 0), std::abs(storedVal)
            );
        }

        return std::make_pair(1, Camel::END_DIFF);
    }

    std::pair<int, float> Decompression::__decompress_decimal(BinObj* compress_data){
        int decimal_count = compress_data->getBits(2) + 1;
        int isM = compress_data->getBits(1);
        uint32_t xored = 0, m_int = 0;

        if (isM == 1) xored = compress_data->getBits(decimal_count) << (23 - decimal_count);
        if (decimal_count <= 1) m_int = compress_data->getBits(3);
        else if (decimal_count ==2) {
            if (compress_data->getBits(1)) m_int = compress_data->getBits(5);
            else m_int = compress_data->getBits(3);
        } 
        else if (decimal_count == 3) {
            int flag = compress_data->getBits(2);
            if (flag == 0) m_int = compress_data->getBits(1);
            else if (flag == 1) m_int = compress_data->getBits(3);
            else if (flag == 2) m_int = compress_data->getBits(5);
            else m_int = compress_data->getBits(Camel::MVALUEBITS[decimal_count-1]);
        } 
        else {
            int flag = compress_data->getBits(2);
            if (flag == 0) m_int = compress_data->getBits(4);
            else if (flag == 1) m_int = compress_data->getBits(6);
            else if (flag == 2) m_int = compress_data->getBits(8);
            else m_int = compress_data->getBits(Camel::MVALUEBITS[decimal_count - 1]);
        }

        float decimalVal, m;
        float scale = (float) Camel::POWERS[decimal_count];
        if (isM == 1){
            m = m_int / scale + 1;
            uint32_t m_prime = std::bit_cast<uint32_t>(m);
            uint32_t decimalLong = xored ^ m_prime;
            decimalVal = std::bit_cast<float>(decimalLong) - 1;
        } 
        else {
            m = (float) m_int / Camel::POWERS[decimal_count];
            decimalVal = m;
        }
    
        return std::make_pair(decimal_count, std::round(decimalVal * scale) / scale);
    }

    float Decompression::__decompress(BinObj* compress_data) {
        if (first) {
            first = false;
            float firstVal = std::bit_cast<float>(compress_data->getBits(32));
            this->storedVal = (int) firstVal;

            return firstVal;
        } 
        else {
            std::pair<int, int> intVal = this->__decompress_integer(compress_data);
            if (intVal.second != Camel::END_DIFF) {
                std::pair<float, float> decimal = this->__decompress_decimal(compress_data);

                return (intVal.first == 0 ? -1: 1) * ( (double) intVal.second + (double) decimal.second);
            }

            return Camel::END_SIGNAL;
        }
    }

    CSVObj* Decompression::decompress(BinObj* compress_data) {
        CSVObj* base_obj = nullptr;
        CSVObj* prev_obj = nullptr;

        int count = 0;
        while (count++ != this->buffer_size) {
            float value = this->__decompress(compress_data);
            if (value == Camel::END_SIGNAL) break;
            else {
                if (base_obj == nullptr) {
                    base_obj = new CSVObj;
                    base_obj->pushData(std::to_string(this->basetime));
                    base_obj->pushData(std::to_string(value));

                    prev_obj = base_obj;
                }
                else {
                    CSVObj* obj = new CSVObj;
                    obj->pushData(std::to_string(this->basetime));
                    obj->pushData(std::to_string(value));

                    prev_obj->setNext(obj);
                    prev_obj = obj;
                }
                this->basetime += this->interval;
            }            
        }
        
        compress_data->flushBits();
        return base_obj;    
    }
    // End: decompression
};