#include "floating-point/lossless.hpp"

namespace Gorilla {
    // Begin: compression
    void Compression::initialize(int count, char** params) {
        this->buffer_size = atoi(params[0]);
        this->bitstream = new BinObj;
    }

    void Compression::finalize() {
        this->__compress(std::bit_cast<uint32_t>(Gorilla::END_SIGNAL));
        this->yield();
        delete this->bitstream;
        this->bitstream = nullptr;
        std::cout << "size: " << this->size << "\n";
    }

    void Compression::__compress(uint32_t value) {
        if (this->first) {
            this->bitstream->put(value, 32);
            this->size += 32;
            this->first = false;
        }
        else {
            uint32_t xored = this->storedVal ^ value;

            if (xored == 0) {
                // Write 0
                this->bitstream->put(0, 1);
                this->size += 1;
            } 
            else {
                int leadingZeros = std::countl_zero(xored);
                int trailingZeros = std::countr_zero(xored);

                // Check overflow of leading? Cannot be 32!
                if (leadingZeros >= 16) leadingZeros = 15;

                // Store bit '1'
                this->bitstream->put(1, 1);
                size += 1;

                if (leadingZeros >= this->storedLeadingZeros && trailingZeros >= this->storedTrailingZeros) {
                    this->bitstream->put(0, 1);
                    int significantBits = 32 - this->storedLeadingZeros - this->storedTrailingZeros;
                    this->bitstream->put(xored >> this->storedTrailingZeros, significantBits);
                    this->size += 1 + significantBits;
                } 
                else {
                    this->bitstream->put(1, 1);
                    this->bitstream->put(leadingZeros, 4); // Number of leading zeros in the next 5 bits

                    int significantBits = 32 - leadingZeros - trailingZeros;
                    if (significantBits == 32) this->bitstream->put(0, 5); // Length of meaningful bits in the next 6 bits
                    else this->bitstream->put(significantBits, 5); // Length of meaningful bits in the next 6 bits

                    this->bitstream->put(xored >> trailingZeros, significantBits); // Store the meaningful bits of xored
                    
                    this->storedLeadingZeros = leadingZeros;
                    this->storedTrailingZeros = trailingZeros;

                    this->size += 1 + 4 + 5 + significantBits;
                }
            }
        }

        this->storedVal = value;
    }

    void Compression::compress(Univariate* data) {
        float value = data->get_value();
        this->__compress(std::bit_cast<uint32_t>(value));

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

    float Decompression::__decompress(BinObj* compress_data) {
        if (this->first) {
            this->first = false;
            this->storedVal = compress_data->getBits(32);
        }
        else {
            if (compress_data->getBits(1) == 1) {
                if (compress_data->getBits(1) == 1) {
                    // New leading and trailing zeros
                    this->storedLeadingZeros = compress_data->getBits(4);
                    int significantBits = compress_data->getBits(5);

                    if(significantBits == 0) significantBits = 32;
                    this->storedTrailingZeros = 32 - significantBits - this->storedLeadingZeros;
                }

                uint32_t value = compress_data->getBits(32 - this->storedLeadingZeros - this->storedTrailingZeros);
                value = value << this->storedTrailingZeros;
                value = this->storedVal ^ value;
                this->storedVal = value;
            }
        }

        return std::bit_cast<float>(this->storedVal);
    }

    CSVObj* Decompression::decompress(BinObj* compress_data) {
        CSVObj* base_obj = nullptr;
        CSVObj* prev_obj = nullptr;

        int count = 0;
        while (count++ != this->buffer_size) {
            float value = this->__decompress(compress_data);
            if (value == Gorilla::END_SIGNAL) break;
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