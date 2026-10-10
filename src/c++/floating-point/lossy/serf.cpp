#include "floating-point/lossy.hpp"

namespace Serf {
    // Begin: Serf material
    uint32_t SerfUtils32::FindAppInt(float min, float max, float v, uint32_t last_int, float max_diff) {
        if (min >= 0) {
            // both positive
            return FindAppInt(min, max, 0, v, last_int, max_diff);
        } 
        else if (max <= 0) {
            // both negative
            return FindAppInt(-max, -min, 0x80000000, v, last_int, max_diff);
        } 
        else if (last_int >> 31 == 0) {
            // consider positive part only, to make more leading zeros
            return FindAppInt(0, max, 0, v, last_int, max_diff);
        } 
        else {
            // consider negative part only, to make more leading zeros
            return FindAppInt(0, -min, 0x80000000, v, last_int, max_diff);
        }
    }

    uint32_t SerfUtils32::FindAppInt(float min_float, float max_float, uint32_t sign, float original, uint32_t last_int, float max_diff) {
        // may be negative zero
        uint32_t min = std::bit_cast<uint32_t>(min_float) & 0x7fffffff;
        uint32_t max = std::bit_cast<uint32_t>(max_float);
        int leading_zeros = std::countl_zero(min ^ max);
        int32_t front_mask = 0xffffffff << (32 - leading_zeros);
        int shift = 32 - leading_zeros;
        uint32_t result_int;
        float diff;
        uint32_t append;

        while (shift >= 0) {
            uint32_t front = front_mask & min;
            uint32_t rear = (~front_mask) & last_int;

            append = rear | front;
            if (append >= min && append <= max) {
                result_int = append ^ sign;
                diff = std::bit_cast<float>(result_int) - original;
                if (diff >= -max_diff && diff <= max_diff) {
                    return result_int;
                }
            }

            append = (append + kBitWeight[shift]) & 0x7fffffff;  // may be overflow
            if (append <= max) {
                // append must be greater than min
                result_int = append ^ sign;
                diff = std::bit_cast<float>(result_int) - original;
                if (diff >= -max_diff && diff <= max_diff) {
                    return result_int;
                }
            }

            front_mask = front_mask >> 1;

            --shift;
        }

        // we do not find a satisfied value, so we return the original value
        return std::bit_cast<uint32_t>(original);
    }
    // End: Serf material

    // Begin: compression
    void Compression::initialize(int count, char** params) {
        this->error = atof(params[0]);
        this->buffer_size = atoi(params[1]);
        this->mode = params[2];
        this->bitstream = new BinObj;

        if (this->mode == "xor") {
            this->__compress = &Compression::__xor_compress;
        }
        else if (this->mode == "qt") {
            this->error = this->error * 0.99f;
            this->__compress = &Compression::__qt_compress;
        }
    }

    void Compression::finalize() {
        (this->*__compress)(Serf::END_SIGNAL);
        this->yield();

        delete this->bitstream;
        this->bitstream = nullptr;
        std::cout << "size: " << this->size << "\n";
    }

    
    void Compression::__qt_compress(float value) {
        if (this->first) {
            this->first = false;
            this->bitstream->put(std::bit_cast<uint32_t>(value), 32);
            this->storedVal = value;
            this->size += 32;
        }
        else if (value == Serf::END_SIGNAL) {
            // Quite inefficient here but it only occurs at the end of stream
            this->size += EliasGammaEncoding::encode(INT_MAX, this->bitstream);
        }
        else {
            long q = static_cast<long>(std::round((value - this->storedVal) / (2 * this->error)));
            float recover_value = this->storedVal + 2 * this->error * q;
        
            this->size += EliasGammaEncoding::encode(ZigZagEncoding::encode(q) + 1, this->bitstream);
            this->storedVal = recover_value;
        }
    }

    void Compression::__xor_compress(float value) {
        float this_val;
        // note we cannot let > maxDiff, because kNan - v > maxDiff is always false
        if (std::abs(this->storedVal - value) > this->error) {
            // in our implementation, we do not consider special cases and overflow case
            this_val = std::bit_cast<float>(
                SerfUtils32::FindAppInt(value - this->error, value + this->error, value, 
                    std::bit_cast<uint32_t>(this->storedVal), this->error)
            );
        } 
        else {
            // let current value be the last value, making an XORed value of 0.
            this_val = this->storedVal;
        }

        // this_val -> value
        uint32_t xor_result = std::bit_cast<uint32_t>(this->storedVal) ^ std::bit_cast<uint32_t>(this_val);

        if (xor_result == 0) {
            this->bitstream->put(1, 2);
            this->size += 2;
        } 
        else {
            int leading_count = std::countl_zero(xor_result);
            int trailing_count = std::countr_zero(xor_result);
            int leading_zeros = leading_round_[leading_count];
            int trailing_zeros = trailing_round_[trailing_count];

            if (leading_zeros >= this->stored_leading_zeros_ && trailing_zeros >= this->stored_trailing_zeros_ &&
                (leading_zeros - this->stored_leading_zeros_) + (trailing_zeros - this->stored_trailing_zeros_)
                < 1 + Serf::LEADING_BITS_PER_VALUE + Serf::TRAILING_BITS_PER_VALUE) 
            {
                // case 1
                int center_bits = 32 - this->stored_leading_zeros_ - this->stored_trailing_zeros_;
                int len = 1 + center_bits;
                if (len > 32) {
                    this->bitstream->put(1, 1);
                    this->bitstream->put(xor_result >> this->stored_trailing_zeros_, center_bits);
                } 
                else {
                    this->bitstream->put((1 << center_bits) | (xor_result >> this->stored_trailing_zeros_), 1 + center_bits);
                }
                this->size += len;
            } 
            else {
                this->stored_leading_zeros_ = leading_zeros;
                this->stored_trailing_zeros_ = trailing_zeros;
                int center_bits = 32 - this->stored_leading_zeros_ - this->stored_trailing_zeros_;

                // case 00
                int len = 2 + Serf::LEADING_BITS_PER_VALUE + Serf::TRAILING_BITS_PER_VALUE + center_bits;
                if (len > 32) {
                    this->bitstream->put((leading_representation_[this->stored_leading_zeros_] 
                        << Serf::TRAILING_BITS_PER_VALUE) | this->trailing_representation_[this->stored_trailing_zeros_], 
                        2 + Serf::LEADING_BITS_PER_VALUE + Serf::TRAILING_BITS_PER_VALUE);

                    this->bitstream->put(xor_result >> this->stored_trailing_zeros_, center_bits);
                } 
                else {
                    this->bitstream->put((((leading_representation_[this->stored_leading_zeros_]
                        << Serf::TRAILING_BITS_PER_VALUE) | this->trailing_representation_[this->stored_trailing_zeros_])
                        << center_bits) | (xor_result >> this->stored_trailing_zeros_), len);
                }
                this->size += len;
            }
        }

        this->storedVal = this_val;
    }

    void Compression::compress(Univariate* data) {
        float value = data->get_value();
        (this->*__compress)(value);

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
        this->error = atof(params[0]);
        this->buffer_size = atoi(params[1]);
        this->mode = params[2];

        if (this->mode == "xor") {
            this->__decompress = &Decompression::__xor_decompress;
        }
        else if (this->mode == "qt") {
            this->error = this->error * 0.99f;
            this->__decompress = &Decompression::__qt_decompress;
        }
    }

    void Decompression::finalize() {
        // Do nothing
    }

    float Decompression::__qt_decompress(BinObj* compress_data) {
        if (this->first) {
            this->first = false;
            this->storedVal = std::bit_cast<float>(compress_data->getBits(32));
        }
        else {
            long elias_res = EliasGammaEncoding::decode(compress_data);
            if (elias_res == INT_MAX) {
                return Serf::END_SIGNAL;
            }
            else {
                long decodeValue = ZigZagEncoding::decode(elias_res - 1);
                float recoverValue = this->storedVal + 2 * this->error * decodeValue;
                this->storedVal = recoverValue;
            }
        }

        return this->storedVal;
    }

    float Decompression::__xor_decompress(BinObj* compress_data) {
        int center_bits;
        uint32_t value = std::bit_cast<uint32_t>(this->storedVal);
        
        if (compress_data->getBits(1) == 1) {
            // case 1
            center_bits = 32 - this->stored_leading_zeros_ - this->stored_trailing_zeros_;
            value = compress_data->getBits(center_bits) << this->stored_trailing_zeros_;
            value = std::bit_cast<uint32_t>(this->storedVal) ^ value;
        } 
        else if (compress_data->getBits(1) == 0) {
            // case 00
            int lead_and_trail = compress_data->getBits(Serf::LEADING_BITS_PER_VALUE + Serf::TRAILING_BITS_PER_VALUE);
            int lead = lead_and_trail >> Serf::TRAILING_BITS_PER_VALUE;
            int trail = ~(0xffff << Serf::TRAILING_BITS_PER_VALUE) & lead_and_trail;
            this->stored_leading_zeros_ = leading_representation_[lead];
            this->stored_trailing_zeros_ = trailing_representation_[trail];
            center_bits = 32 - this->stored_leading_zeros_ - this->stored_trailing_zeros_;

            value = compress_data->getBits(center_bits) << this->stored_trailing_zeros_;
            value = std::bit_cast<uint32_t>(this->storedVal) ^ value;
        }

        this->storedVal = std::bit_cast<float>(value);
        return this->storedVal;
    }

    CSVObj* Decompression::decompress(BinObj* compress_data) {
        CSVObj* base_obj = nullptr;
        CSVObj* prev_obj = nullptr;

        int count = 0;
        while (count++ != this->buffer_size) {
            float value = (this->*__decompress)(compress_data);
            if (value == Serf::END_SIGNAL) break;
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