#include "floating-point/lossless.hpp"

namespace Chimp {
    // Begin: compression
    void Compression::initialize(int count, char** params) {
        this->buffer_size = atoi(params[0]);
        this->mode = params[1];
        this->bitstream = new BinObj;

        if (this->mode == "basic") this->__compress = &Compression::__basic_compress;
        else if (this->mode == "128") this->__compress = &Compression::__128_compress;

        if (this->mode == "128") {
            this->previousValues = 64;
            this->previousValuesLog2 = (int) log2(this->previousValues);
            this->threshold = 5 + previousValuesLog2;
            this->setLsb = (int) pow(2, threshold + 1) - 1;
            this->indices = new int[(int) pow(2, threshold + 1)];
            this->storedValues = new uint32_t[previousValues];
        }
    }

    void Compression::finalize() {
        (this->*__compress)(Chimp::END_SIGNAL);
        this->yield();

        delete this->bitstream;
        this->bitstream = nullptr;
        std::cout << "size: " << this->size << "\n";

        if (this->mode == "128") {
            delete this->indices;
            delete this->storedValues;
        }
    }

    void Compression::__basic_compress(uint32_t value) {
        if (this->first) {
            this->first = false;
            this->storedVal = value;
            this->bitstream->put(this->storedVal, 32);
            this->size += 32;
        }
        else {
            uint32_t xored = this->storedVal ^ value;
            if (xored == 0) {
                // Write 0
                this->bitstream->put(0, 2);
                this->storedLeadingZeros = 33;
                this->size += 2;
            } 
            else {
                int leadingZeros = leadingRound[std::countl_zero(xored)];
                int trailingZeros = std::countr_zero(xored);

                if (trailingZeros > 5) {
                    int significantBits = 32 - leadingZeros - trailingZeros;
                    this->bitstream->put(1, 2);
                    this->bitstream->put(leadingRepresentation[leadingZeros], 3);
                    this->bitstream->put(significantBits, 5);
                    this->bitstream->put(xored >> trailingZeros, significantBits); // Store the meaningful bits of XOR
                    this->size += 10 + significantBits;
                    this->storedLeadingZeros = 33;
                } 
                else if (leadingZeros == this->storedLeadingZeros) {
                    this->bitstream->put(2, 2);
                    int significantBits = 32 - leadingZeros;
                    this->bitstream->put(xored, significantBits);
                    this->size += 2 + significantBits;
                } 
                else {
                    this->storedLeadingZeros = leadingZeros;
                    int significantBits = 32 - leadingZeros;
                    this->bitstream->put(24 + leadingRepresentation[leadingZeros], 5);
                    this->bitstream->put(xored, significantBits);
                    this->size += 5 + significantBits;
                }
            }
            this->storedVal = value;
        }
    }

    void Compression::__128_compress(uint32_t value) {
        if (this->first) {
            this->first = false;
            this->storedValues[this->current] = value;
            this->bitstream->put(value, 32);
            indices[value & setLsb] = this->index;
            this->size += 32;
        }
        else {
            uint32_t xored;
            int previousIndex;
            int trailingZeros;
            int key = value & setLsb;
            int currIndex = indices[key];

            if ((this->index - currIndex) < this->previousValues) {
                uint32_t tempXor = value ^ this->storedValues[currIndex % this->previousValues];
                trailingZeros = std::countr_zero(tempXor);

                if (trailingZeros > this->threshold) {
                    previousIndex = currIndex % this->previousValues;
                    xored = tempXor;
                } 
                else {
                    previousIndex = this->index % this->previousValues;
                    xored = this->storedValues[previousIndex] ^ value;
                    trailingZeros = std::countr_zero(xored);
                }
            } 
            else {
                previousIndex = this->index % this->previousValues;
                xored = this->storedValues[previousIndex] ^ value;
                trailingZeros = std::countr_zero(xored);
            }

            if (xored == 0) {
                // Write 0
                this->bitstream->put(0, 2);
                this->bitstream->put(previousIndex, this->previousValuesLog2);
                
                this->size += 2 + this->previousValuesLog2;
                this->storedLeadingZeros = 33;
            } 
            else {
                int leadingZeros = std::countl_zero(xored);

                if (trailingZeros > this->threshold) {
                    int significantBits = 32 - leadingRound[leadingZeros] - trailingZeros;
                    this->bitstream->put(256 * (this->previousValues + previousIndex) + 32 * leadingRepresentation[leadingZeros] + significantBits, this->previousValuesLog2 + 10);
                    this->bitstream->put(xored >> trailingZeros, significantBits); // Store the meaningful bits of XOR
                    this->size += 10 + significantBits + this->previousValuesLog2;
                    this->storedLeadingZeros = 33;
                } 
                else if (leadingRound[leadingZeros] == storedLeadingZeros) {
                    this->bitstream->put(1, 1);
                    this->bitstream->put(0, 1);
                    int significantBits = 32 - leadingRound[leadingZeros];
                    this->bitstream->put(xored, significantBits);
                    this->size += 2 + significantBits;
                } 
                else {
                    this->storedLeadingZeros = leadingRound[leadingZeros];
                    int significantBits = 32 - leadingRound[leadingZeros];
                    this->bitstream->put(16 + 8 + leadingRepresentation[leadingZeros], 5);
                    this->bitstream->put(xored, significantBits);
                    this->size += 5 + significantBits;
                }
            }
            this->current = ((this->current + 1) % this->previousValues);
            this->storedValues[this->current] = std::bit_cast<uint32_t>(value);
            this->indices[key] = ++this->index;
        }
    }

    void Compression::compress(Univariate* data) {
        float value = data->get_value();
        (this->*__compress)(std::bit_cast<uint32_t>(value));

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
        this->previousValues = 64;
        this->buffer_size = atoi(params[0]);
        this->mode = params[1];

        if (this->mode == "basic") this->__decompress = &Decompression::__basic_decompress;
        else if (this->mode == "128") this->__decompress = &Decompression::__128_decompress;

        if (this->mode == "128") {
            this->previousValuesLog2 = (int) log2(previousValues);
            this->storedValues = new uint32_t[previousValues];
        }
    }

    void Decompression::finalize() {
        if (this->mode == "128") delete this->storedValues;
    }

    float Decompression::__basic_decompress(BinObj* compress_data) {
        if (this->first) {
            this->first = false;
            this->storedVal = compress_data->getBits(32);
        }
        else {
            if (compress_data->getBits(1) == 1) {
                if (compress_data->getBits(1) == 1) {
                    // New leading zeros
                    this->storedLeadingZeros = leadingRepresentation[
                        compress_data->getBits(3)
                    ];
                }
                int significantBits = 32 - this->storedLeadingZeros;
                if (significantBits == 0) significantBits = 32;
                
                uint32_t value = compress_data->getBits(32 - this->storedLeadingZeros);
                value = this->storedVal ^ value;
                this->storedVal = value;
            } 
            else if (compress_data->getBits(1) == 1) {
                this->storedLeadingZeros = leadingRepresentation[compress_data->getBits(3)];
                int significantBits = compress_data->getBits(5);
                if (significantBits == 0) significantBits = 32;
                
                this->storedTrailingZeros = 32 - significantBits - this->storedLeadingZeros;
                uint32_t value = compress_data->getBits(32 - this->storedLeadingZeros - this->storedTrailingZeros);
                value <<= this->storedTrailingZeros;
                value = this->storedVal ^ value;
                this->storedVal = value;
            }
        }

        if (this->storedVal == Chimp::END_SIGNAL) return Chimp::END_SIGNAL;
        else return std::bit_cast<float>(this->storedVal);
    }

    float Decompression::__128_decompress(BinObj* compress_data) {
        uint32_t value = 0;

        if (this->first) {
            this->first = false;
            value = compress_data->getBits(32);
        } 
        else {
            if (compress_data->getBits(1) == 1) {
                if (compress_data->getBits(1) == 1) {
                    // New leading zeros
                    this->storedLeadingZeros = leadingRepresentation[
                        compress_data->getBits(3)
                    ];
                } 

                int significantBits = 32 - this->storedLeadingZeros;
                if (significantBits == 0) significantBits = 32;
                
                value = compress_data->getBits(32 - this->storedLeadingZeros);
                value = this->storedVal ^ value;
            } 
            else if (compress_data->getBits(1) == 1) {
                int fill = this->previousValuesLog2 + 8;
                uint32_t temp = compress_data->getBits(fill);
                int index = (temp >> (fill -= this->previousValuesLog2)) & ((1 << this->previousValuesLog2) - 1);
                this->storedLeadingZeros = leadingRepresentation[(temp >> (fill -= 3)) & ((1 << 3) - 1)];
                int significantBits = (temp >> (fill -= 5)) & ((1 << 5) - 1);
                this->storedVal = this->storedValues[index];
                if (significantBits == 0) significantBits = 32;

                this->storedTrailingZeros = 32 - significantBits - this->storedLeadingZeros;
                value = compress_data->getBits(32 - this->storedLeadingZeros - this->storedTrailingZeros);
                value <<= this->storedTrailingZeros;
                value = this->storedVal ^ value;    
            } 
            else {
                // else -> same value as before
                uint32_t index = compress_data->getBits(this->previousValuesLog2);
                value = this->storedValues[index];
            }
        }

        if (value != Chimp::END_SIGNAL) {
            this->storedVal = value;
            this->current = (this->current + 1) % this->previousValues;
            this->storedValues[this->current] = this->storedVal;

            return std::bit_cast<float>(value);
        }

        return Chimp::END_SIGNAL;
    }

    CSVObj* Decompression::decompress(BinObj* compress_data) {
        CSVObj* base_obj = nullptr;
        CSVObj* prev_obj = nullptr;

        int count = 0;
        while (count++ != this->buffer_size) {
            float value = (this->*__decompress)(compress_data);
            if (value == Chimp::END_SIGNAL) break;
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