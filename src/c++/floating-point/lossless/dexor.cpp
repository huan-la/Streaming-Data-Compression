#include "floating-point/lossless.hpp"
#include <charconv>
#include <limits>

namespace DeXOR {
    // Float32 port of the authors' DoubleDeXOREncoder.Method.Decimal_XOR
    // (Native) and DoubleDeXORDecoder.Native. Keep their forward prefix search,
    // four control codes and independent adaptive exception state.
    // Begin: DeXORTools materials
    const int DeXORTools::cost[7] = {0, 4, 7, 10, 14, 17, 20};

    // Double intermediates avoid rounding decimal powers to binary32 before use.
    // The five-bit q field covers [-20, 11]; q + delta is at most 17.
    const double DeXORTools::P10[38] = {
        1e-20, 1e-19, 1e-18, 1e-17, 1e-16, 1e-15, 1e-14,
        1e-13, 1e-12, 1e-11, 1e-10, 1e-9, 1e-8, 1e-7,
        1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1.0,
        1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9,
        1e10, 1e11, 1e12, 1e13, 1e14, 1e15, 1e16, 1e17
    };

    // Reciprocals of the same double literals, rounded once at compile time.
    const double DeXORTools::I10[38] = {
        1 / 1e-20, 1 / 1e-19, 1 / 1e-18, 1 / 1e-17, 1 / 1e-16, 1 / 1e-15, 1 / 1e-14,
        1 / 1e-13, 1 / 1e-12, 1 / 1e-11, 1 / 1e-10, 1 / 1e-9, 1 / 1e-8, 1 / 1e-7,
        1 / 1e-6, 1 / 1e-5, 1 / 1e-4, 1 / 1e-3, 1 / 1e-2, 1 / 1e-1, 1.0,
        1 / 1e1, 1 / 1e2, 1 / 1e3, 1 / 1e4, 1 / 1e5, 1 / 1e6, 1 / 1e7, 1 / 1e8, 1 / 1e9,
        1 / 1e10, 1 / 1e11, 1 / 1e12, 1 / 1e13, 1 / 1e14, 1 / 1e15, 1 / 1e16, 1 / 1e17
    };

    double DeXORTools::getP10(int pow) {
        return DeXORTools::P10[pow - DeXORTools::MIN_Q];
    }

    int DeXORTools::getP2(int pow) {
        return 1 << pow;
    }

    bool DeXORTools::__isEnd(float value, int q) {
        int index = q - DeXORTools::MIN_Q;
        double pow = DeXORTools::P10[index];
        double scaled = std::round(static_cast<double>(value) * DeXORTools::I10[index]);
        if (scaled == 0 && value != 0) return false;

        // A decimal coordinate is useful only if its integer grid preserves the
        // input float. An absolute epsilon would collapse small values to zero.
        double decimal = scaled * pow;
        static constexpr double max_float = std::numeric_limits<float>::max();
        if (std::fabs(decimal) > max_float) return false;
        float decoded = static_cast<float>(decimal);
        // getEnd handles both zeros separately; finite nonzero floats compare
        // equal exactly when their binary32 representations match.
        return decoded == value;
    }

    int DeXORTools::getEnd(float value, int last_end) {
        if (value == 0) return 0;

        int q = last_end;
        if (q < DeXORTools::MIN_Q) q = DeXORTools::MIN_Q;
        else if (q > DeXORTools::MAX_Q) q = DeXORTools::MAX_Q;
        if (DeXORTools::__isEnd(value, q)) {
            while (q < DeXORTools::MAX_Q && DeXORTools::__isEnd(value, q + 1)) ++q;
        }
        else {
            // The current coordinate has already failed; probe each lower one once.
            do {
                if (q == DeXORTools::MIN_Q) return DeXORTools::MIN_Q - 1;
                --q;
            } while (!DeXORTools::__isEnd(value, q));
        }

        return q;
    }

    int DeXORTools::decimalBits(int dp) {
        return DeXORTools::cost[dp];
    }

    double DeXORTools::truncate(double value) {
        // Preserve the float32 port's near-integer tolerance and signed zeros.
        double nearest = std::round(value);
        static constexpr double float_epsilon = std::numeric_limits<float>::epsilon();
        double eps = std::fabs(value) * float_epsilon;
        if (eps < 1e-8) eps = 1e-8;
        if (std::fabs(value - nearest) <= eps) return nearest;
        return std::trunc(value);
    }

    uint32_t DeXORTools::segment(uint32_t v, int st, int ed) {
        int len = ed - st + 1;
        uint32_t mask = UINT32_MAX >> (32 - len);
        return (v >> (32 - ed)) & mask;
    }
    // End: DeXORTools materials

    // Begin: compression
    void Compression::initialize(int count, char** params) {
        this->buffer_size = atoi(params[0]);
        this->bitstream = new BinObj;
        // Reserve bytes only; BinObj still writes the authors' fields in order.
        this->bitstream->byte_vector.reserve(6 * std::clamp(this->buffer_size, 1, 4096));
    }

    void Compression::finalize() {
        this->__compress(DeXOR::END_SIGNAL, true);
        this->yield();

        delete this->bitstream;
        this->bitstream = nullptr;
        std::cout << "size: " << this->size << "\n";
    }

    void Compression::__compress(float value, bool end) {
        uint32_t raw = std::bit_cast<uint32_t>(value);
        if (!std::isfinite(value)) {
            this->__exceptionHandle(raw, end);
            return;
        }

        int q = DeXORTools::getEnd(value, this->previous_q);
        if (q < DeXORTools::MIN_Q) {
            this->__exceptionHandle(raw, end);
            return;
        }

        const double* powers = DeXORTools::P10 + (q - DeXORTools::MIN_Q);
        const double* scales = DeXORTools::I10 + (q - DeXORTools::MIN_Q);
        double current_value = value;
        double reference_value = this->previous_value;
        int delta = 0;
        double alpha = 0;
        while (delta <= DeXORTools::MAX_DELTA) {
            double pow = powers[delta];
            double scale = scales[delta];
            double a = DeXORTools::truncate(current_value * scale);
            double b = DeXORTools::truncate(reference_value * scale);
            if (a == b) {
                alpha = b * pow;
                break;
            }
            ++delta;
        }
        // A failed float32 prefix search already determines the exception path.
        if (delta > DeXORTools::MAX_DELTA) {
            this->__exceptionHandle(raw, end);
            return;
        }

        double pow = powers[0];
        double beta = std::round((current_value - alpha) / pow);
        if (this->__cannotEncode(raw, alpha, beta, pow, delta)) {
            this->__exceptionHandle(raw, end);
            return;
        }

        uint32_t magnitude = static_cast<uint32_t>(std::fabs(beta));
        int bits = DeXORTools::cost[delta];
        // Author's controls: 00 new q/delta, 01 new delta, 10 reuse alpha;
        // 11 above selects the independent exponent-coded exception path.
        bool same_q = q == this->previous_q;
        uint32_t record;
        int record_bits = 2;
        // 10 reuses the decoder's alpha; preserve its sign as well as its value.
        if (same_q && delta == this->previous_delta &&
            alpha == this->previous_alpha &&
            (alpha != 0 || std::signbit(alpha) == std::signbit(this->previous_alpha))) {
            record = 0b10;
        }
        else {
            record = same_q ? 0b01 : 0b00;
            if (!same_q) {
                record = (record << 5) | static_cast<uint32_t>(q - DeXORTools::MIN_Q);
                record_bits += 5;
                this->previous_q = q;
            }
            record = (record << DeXORTools::DELTA_BITS) | static_cast<uint32_t>(delta);
            record_bits += DeXORTools::DELTA_BITS;
            this->previous_delta = delta;
            this->previous_alpha = alpha;
        }

        if (alpha == 0) {
            record = (record << 1) | ((raw >> 31) ^ 1u);
            ++record_bits;
        }
        // The author's control | q? | delta? | sign? | suffix layout is at most
        // 2 + 5 + 3 + 1 + 20 = 31 bits for float32, so one BinObj write suffices.
        record = (record << bits) | magnitude;
        record_bits += bits;
        this->bitstream->put(record, record_bits);
        this->size += record_bits;
        this->previous_value = value;
    }

    bool Compression::__cannotEncode(uint32_t raw, double alpha, double beta, double pow, int delta) {
        // Float32 counterpart of the authors' cannotDecodeBitwise helper.
        if (delta > DeXORTools::MAX_DELTA) return true;
        double magnitude = std::fabs(beta);
        if (!(magnitude < (1u << DeXORTools::cost[delta]))) return true;
        double sign = alpha == 0 ? (raw >> 31 ? -1 : 1) : (alpha > 0 ? 1 : -1);
        float decoded = static_cast<float>(alpha + sign * static_cast<uint32_t>(magnitude) * pow);
        return std::bit_cast<uint32_t>(decoded) != raw;
    }

    void Compression::__exceptionHandle(uint32_t raw, bool end) {
        int exp = (raw >> 23) & 0xffu;
        int delta = exp - this->previous_exp;
        int bias = (1 << (this->EL - 1)) - 1;

        if (delta >= -bias && delta <= bias) {
            uint32_t header = (0b11u << this->EL) | static_cast<uint32_t>(delta + bias);
            uint32_t payload = ((raw >> 31) << 23) | (raw & 0x7fffffu);
            // Keep 11 | exponent delta | sign | fraction in the author's order.
            // Wider records retain two writes to stay within BinObj's 32-bit API.
            if (this->EL <= 6) this->bitstream->put((header << 24) | payload, this->EL + 26);
            else {
                this->bitstream->put(header, this->EL + 2);
                this->bitstream->put(payload, 24);
            }
            this->size += this->EL + 26;

            if (this->EL > 1) {
                int smaller_bias = (1 << (this->EL - 2)) - 1;
                if (delta >= -smaller_bias && delta <= smaller_bias) ++this->contract_step;
                else this->contract_step = 0;
                if (this->contract_step == this->rho) {
                    --this->EL;
                    this->contract_step = 0;
                }
            }
        }
        else {
            // 11 followed by the all-ones exponent escape is one all-ones field.
            this->bitstream->put((1u << (this->EL + 2)) - 1, this->EL + 2);
            this->bitstream->put(raw, 32);
            this->size += this->EL + 34;
            this->contract_step = 0;
            // Eight exponent bits require at most nine signed delta bits.
            if (this->EL < 9) ++this->EL;
        }
        this->previous_exp = exp;
        // Exceptions do not replace the decimal path's reference or metadata.
        // Distinguish an infinity data point from the harness's completion marker.
        if (raw == 0x7f800000u) {
            this->bitstream->put(!end, 1);
            ++this->size;
        }
    }

    void Compression::compress(Univariate* data) {
        this->__compress(static_cast<float>(data->get_value()));
        if (++this->count == this->buffer_size) {
            this->yield();
            this->count = 0;
        }
    }

    BinObj* Compression::serialize() {
        BinObj* obj = this->bitstream;
        this->bitstream = new BinObj;
        this->bitstream->byte_vector.reserve(6 * std::clamp(this->buffer_size, 1, 4096));
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
        int con = compress_data->getBits(2);
        if (con == 3) {
            float value = this->__exceptionDecode(compress_data);
            if (std::bit_cast<uint32_t>(value) == 0x7f800000u) {
                this->finished = compress_data->getBits(1) == 0;
            }
            return value;
        }

        if (con == 0 || con == 1) {
            if (con == 0) this->previous_q = static_cast<int>(compress_data->getBits(5)) + DeXORTools::MIN_Q;
            this->previous_delta = compress_data->getBits(DeXORTools::DELTA_BITS);
            double pow = DeXORTools::getP10(this->previous_q + this->previous_delta);
            double scale = DeXORTools::I10[this->previous_q + this->previous_delta - DeXORTools::MIN_Q];
            this->previous_alpha = DeXORTools::truncate(static_cast<double>(this->previous_value) * scale) * pow;
        }

        double sign = this->previous_alpha > 0 ? 1 : -1;
        if (this->previous_alpha == 0) sign = compress_data->getBits(1) ? 1 : -1;
        uint32_t magnitude = compress_data->getBits(DeXORTools::decimalBits(this->previous_delta));
        double beta = sign * magnitude * DeXORTools::getP10(this->previous_q);
        this->previous_value = static_cast<float>(this->previous_alpha + beta);
        return this->previous_value;
    }

    float Decompression::__exceptionDecode(BinObj* compress_data) {
        int bias = DeXORTools::getP2(this->EL - 1) - 1;
        int delta = static_cast<int>(compress_data->getBits(this->EL)) - bias;
        uint32_t raw;

        if (delta >= -bias && delta <= bias) {
            this->previous_exp += delta;
            uint32_t sign = compress_data->getBits(1);
            uint32_t fraction = compress_data->getBits(23);
            raw = (sign << 31) | (static_cast<uint32_t>(this->previous_exp) << 23) | fraction;

            if (this->EL > 1) {
                int smaller_bias = DeXORTools::getP2(this->EL - 2) - 1;
                if (delta >= -smaller_bias && delta <= smaller_bias) ++this->contract_step;
                else this->contract_step = 0;
                if (this->contract_step == this->rho) {
                    --this->EL;
                    this->contract_step = 0;
                }
            }
        }
        else {
            raw = compress_data->getBits(32);
            this->previous_exp = DeXORTools::segment(raw, 2, 9);
            this->contract_step = 0;
            if (this->EL < 9) ++this->EL;
        }

        return std::bit_cast<float>(raw);
    }

    CSVObj* Decompression::decompress(BinObj* compress_data) {
        CSVObj* base_obj = nullptr;
        CSVObj* prev_obj = nullptr;

        int count = 0;
        while (!this->finished && count++ < this->buffer_size) {
            float value = this->__decompress(compress_data);
            if (value == DeXOR::END_SIGNAL) break;
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
