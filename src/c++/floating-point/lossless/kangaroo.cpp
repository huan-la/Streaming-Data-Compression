#include "floating-point/lossless.hpp"
#include <cerrno>
#include <charconv>
#include <cstdlib>
#include <limits>
#include <stdexcept>

namespace Kangaroo {
    namespace {
        // Construct the small decimal-coordinate tables at compile time,
        // removing ceil/floor calls from the per-sample path.
        struct DecimalTables {
            int bits[46];
            int precision[288];

            constexpr DecimalTables() : bits{}, precision{} {
                for (int alpha = -1; alpha < 45; ++alpha) {
                    double value = alpha * 3.32192809488736234787;
                    int integer = static_cast<int>(value);
                    this->bits[alpha + 1] = integer + (value > integer);
                }
                for (int shift = -137; shift <= 150; ++shift) {
                    double value = (shift - 1) * 0.30102999566398119521;
                    int integer = static_cast<int>(value);
                    this->precision[shift + 137] = integer - (value < integer) + 1;
                }
            }
        };
        constexpr DecimalTables DECIMAL_TABLES;

        int positiveParameter(const char* text) {
            char* end = nullptr;
            errno = 0;
            long value = std::strtol(text, &end, 10);
            if (errno == ERANGE || end == text || *end != '\0' || value < 1 || value > INT_MAX) {
                throw std::invalid_argument("Kangaroo size and window must be positive integers");
            }
            return static_cast<int>(value);
        }

        void remember(std::vector<uint32_t>& history, uint64_t& index, int& slot, uint32_t value, uint32_t repeats) {
            // Repeated samples occupy the window too, without repeating the
            // reference search or decimal erasure for every sample in a run.
            int window = history.size();
            if (repeats == 1) {
                history[slot] = value;
                if (++slot == window) slot = 0;
            }
            else if (repeats >= history.size()) {
                std::fill(history.begin(), history.end(), value);
                slot = (static_cast<uint64_t>(slot) + repeats) % history.size();
            }
            else {
                for (uint32_t i = 0; i < repeats; ++i) {
                    history[slot] = value;
                    if (++slot == window) slot = 0;
                }
            }
            index += repeats;
        }
    }

    // Begin: Kangaroo32Utils materials
    const double Kangaroo32Utils::P10[45] = {
        1.0, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9,
        1e10, 1e11, 1e12, 1e13, 1e14, 1e15, 1e16, 1e17, 1e18, 1e19,
        1e20, 1e21, 1e22, 1e23, 1e24, 1e25, 1e26, 1e27, 1e28, 1e29,
        1e30, 1e31, 1e32, 1e33, 1e34, 1e35, 1e36, 1e37, 1e38, 1e39,
        1e40, 1e41, 1e42, 1e43, 1e44
    };

    // The same three-bit leading-zero buckets as the framework's Chimp32.
    const int Kangaroo32Utils::leading[8] = {0, 6, 10, 12, 14, 16, 18, 20};

    int Kangaroo32Utils::__eraseBits(int exponent, int alpha) {
        return 23 - exponent - DECIMAL_TABLES.bits[alpha + 1];
    }

    int Kangaroo32Utils::__decimalPlaces(uint32_t value, int exponent) {
        int trailing = std::countr_zero(value);
        // A zero mantissa cannot gain trailing zeros through erasure.
        if (trailing >= 23) return -1;
        // Only alpha <= floor((22-E-trailing)*log10(2)) can erase more
        // bits than are already zero. Decimal grids are nested, so if the
        // finest eligible grid fails to round-trip, every coarser one fails.
        int alpha = DECIMAL_TABLES.precision[23 - exponent - trailing + 137] - 1;
        if (alpha < 0) return -1;
        double magnitude = std::bit_cast<float>(value & 0x7fffffffu);
        double scaled = magnitude * Kangaroo32Utils::P10[alpha];
        // The eligible-grid bound also guarantees scaled < 2^23. Round
        // positive coordinates without a library call or half-way addition.
        uint32_t coordinate = static_cast<uint32_t>(scaled);
        coordinate += scaled - coordinate >= 0.5;
        float decimal = static_cast<float>(coordinate / Kangaroo32Utils::P10[alpha]);
        if (std::bit_cast<uint32_t>(decimal) != (value & 0x7fffffffu)) return -1;
        // The grid spacing exceeds one float ULP, making its matching
        // coordinate unique. Removing decimal zeros finds the same smallest
        // alpha as the former linear search.
        while (alpha > 0 && coordinate % 10 == 0) {
            coordinate /= 10;
            --alpha;
        }
        return alpha;
    }

    int Kangaroo32Utils::__recoverPrecision(uint32_t value) {
        int exponent = static_cast<int>((value >> 23) & 0xffu) - 127;
        int trailing = std::countr_zero(value);
        // Invert L(alpha) <= trailing < L(alpha - 1), including the
        // ceiling in L(alpha). Using floor((23-E-trailing)*log10(2))
        // alone misidentifies the precision of some bit-flipped values.
        return DECIMAL_TABLES.precision[23 - exponent - trailing + 137];
    }

    uint32_t Kangaroo32Utils::erase(uint32_t value, int& tag) {
        tag = 0b111;
        int exponent_field = (value >> 23) & 0xffu;
        // Retain zeros, subnormals, infinities and every NaN payload verbatim.
        if (exponent_field == 0 || exponent_field == 255) return value;

        int exponent = exponent_field - 127;
        if (exponent >= 23) return value;
        int alpha = Kangaroo32Utils::__decimalPlaces(value, exponent);
        if (alpha < 0) return value;

        int length = std::min(23, Kangaroo32Utils::__eraseBits(exponent, alpha));
        uint32_t erased = value & (UINT32_MAX << length);
        int trailing = std::countr_zero(erased);
        if (trailing <= std::countr_zero(value)) return value;

        int upper = Kangaroo32Utils::__eraseBits(exponent, alpha - 1);
        int candidate_tag = 0;
        if (trailing >= upper) {
            // The precision marker must stay inside the 23-bit mantissa.
            if (upper < 1 || upper > 23) return value;
            erased |= 1u << (upper - 1);
            candidate_tag = 0b10;
        }

        if (Kangaroo32Utils::__recoverPrecision(erased) != alpha ||
            Kangaroo32Utils::recover(erased, candidate_tag) != value) return value;
        tag = candidate_tag;
        return erased;
    }

    uint32_t Kangaroo32Utils::recover(uint32_t value, int tag) {
        if (tag == 0b111) return value;

        int alpha = Kangaroo32Utils::__recoverPrecision(value);
        if (tag == 0b10) {
            int trailing = std::countr_zero(value);
            value &= ~(1u << trailing);
        }

        float erased = std::bit_cast<float>(value);
        double scale = Kangaroo32Utils::P10[alpha];
        double magnitude = std::ceil(std::abs(static_cast<double>(erased)) * scale) / scale;
        float recovered = std::copysign(static_cast<float>(magnitude), erased);
        return std::bit_cast<uint32_t>(recovered);
    }

    int Kangaroo32Utils::leadingCode(int zeros) {
        static const unsigned char codes[32] = {
            0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 3, 3, 4, 4,
            5, 5, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7
        };
        return codes[zeros];
    }

    int Kangaroo32Utils::leadingZeros(int code) {
        return Kangaroo32Utils::leading[code];
    }
    // End: Kangaroo32Utils materials

    // Begin: compression
    void Compression::initialize(int count, char** params) {
        if (count != 2) throw std::invalid_argument("Kangaroo requires size and window parameters");
        this->buffer_size = positiveParameter(params[0]);
        this->window_size = positiveParameter(params[1]);
        this->power_of_two = std::has_single_bit(static_cast<uint32_t>(this->window_size));
        this->window_mask = this->window_size - 1;
        this->index_bits = std::bit_width(static_cast<uint32_t>(this->window_size - 1));
        this->stored_values.assign(this->window_size, 0);
        this->history_index = 0;
        this->history_slot = 0;
        this->run_count = 0;
        this->count = 0;
        this->size = 0;
        std::fill(this->previous_indices, this->previous_indices + 33, UINT64_MAX);
        std::fill(this->lookup_mask, this->lookup_mask + 33, 0);
        std::fill(this->cache_valid, this->cache_valid + 256, false);
        this->cache_enabled = true;
        this->cache_queries = 0;
        this->cache_hits = 0;
        this->bitstream = new BinObj;
    }

    void Compression::finalize() {
        this->__flushRun();
        // 110 with run length zero terminates the stream independently of
        // float values, keeping both infinities available as ordinary data.
        this->__write(0b110, 3);
        this->__writeRun(0);
        this->yield();

        delete this->bitstream;
        this->bitstream = nullptr;
        std::vector<uint32_t>().swap(this->stored_values);
        std::cout << "size: " << this->size << "\n";
    }

    void Compression::__write(uint64_t value, int bits) {
        this->size += bits;
        if (bits == 0) return;
        // Write the packed record directly to BinObj's byte buffer. Casting
        // each chunk to a byte discards its consumed high bits, avoiding the
        // masks and second put() call needed to split a 64-bit record.
        int position = this->bitstream->wBitPos;
        if (position != 0) {
            int available = 8 - position;
            if (bits <= available) {
                this->bitstream->byte_vector.back() |= static_cast<unsigned char>(value << (available - bits));
                this->bitstream->wBitPos = (position + bits) % 8;
                return;
            }
            bits -= available;
            this->bitstream->byte_vector.back() |= static_cast<unsigned char>(value >> bits);
        }
        // Short records need at most two pushes. Grow longer records once,
        // then fill their bytes without a vector capacity check per byte.
        if (bits <= 16) {
            if (bits >= 8) {
                bits -= 8;
                this->bitstream->byte_vector.push_back(static_cast<unsigned char>(value >> bits));
            }
            if (bits != 0) this->bitstream->byte_vector.push_back(static_cast<unsigned char>(value << (8 - bits)));
            this->bitstream->wBitPos = bits % 8;
            return;
        }
        size_t start = this->bitstream->byte_vector.size();
        this->bitstream->byte_vector.resize(start + (bits + 7) / 8);
        unsigned char* output = this->bitstream->byte_vector.data() + start;
        while (bits >= 8) {
            bits -= 8;
            *output++ = static_cast<unsigned char>(value >> bits);
        }
        if (bits != 0) *output = static_cast<unsigned char>(value << (8 - bits));
        this->bitstream->wBitPos = bits;
    }

    void Compression::__writeRun(uint32_t length) {
        // Base-128 groups, low group first; MSB=1 means another group follows.
        do {
            uint32_t byte = length & 0x7fu;
            length >>= 7;
            this->__write(byte | (length != 0 ? 0x80u : 0), 8);
        } while (length != 0);
    }

    uint64_t Compression::__findReference(uint32_t value, int trailing) {
        uint64_t previous = this->previous_indices[trailing];
        if (previous == UINT64_MAX || this->history_index - previous > static_cast<uint64_t>(this->window_size)) {
            return UINT64_MAX;
        }

        int leading = std::countl_zero(value ^ this->previous_values[trailing]);
        if (leading == 32) return previous;
        if ((this->lookup_mask[trailing] & (1u << leading)) != 0) {
            uint64_t reference = this->lookup[trailing][leading];
            if (this->history_index - reference <= static_cast<uint64_t>(this->window_size)) return reference;
        }
        return previous;
    }

    void Compression::__remember(uint32_t value, uint32_t repeats, int trailing) {
        uint64_t previous = this->previous_indices[trailing];
        if (previous == UINT64_MAX || this->history_index - previous > static_cast<uint64_t>(this->window_size)) {
            this->lookup_mask[trailing] = 0;
        }
        else {
            int leading = std::countl_zero(value ^ this->previous_values[trailing]);
            if (leading < 32) {
                // A new adjacent comparison below q invalidates its old hop;
                // equality records the boundary that improves a future XOR
                // with q leading zeros. Lower q entries remain valid.
                this->lookup[trailing][leading] = previous;
                this->lookup_mask[trailing] = (this->lookup_mask[trailing] & (UINT32_MAX >> (31 - leading))) | (1u << leading);
            }
        }
        remember(this->stored_values, this->history_index, this->history_slot, value, repeats);
        this->previous_indices[trailing] = this->history_index - 1;
        this->previous_values[trailing] = value;
        this->previous_value = value;
    }

    void Compression::__xorCompress(uint32_t value, int tag, int trailing) {
        uint64_t packet = tag;
        int bits = tag == 0 ? 1 : (tag == 0b10 ? 2 : 3);
        if (this->history_index == 0) {
            this->__write((packet << 32) | value, bits + 32);
            return;
        }

        uint64_t reference = this->__findReference(value, trailing);
        uint32_t previous;
        if (reference == UINT64_MAX) {
            previous = this->previous_value;
            packet = (packet << 2) | 0b01;
            bits += 2;
            // Six bits also represent 32 for the all-zero binary32 value.
            packet = (packet << 6) | trailing;
            bits += 6;
        }
        else {
            uint32_t slot = this->power_of_two ? (reference & this->window_mask) : (reference % this->window_size);
            previous = this->stored_values[slot];
            if (previous == value) {
                packet <<= 2; // 00: identical transformed value
                bits += 2;
                this->__write((packet << this->index_bits) | slot, bits + this->index_bits);
                return;
            }
            packet = (packet << 1) | 0b1;
            bits += 1;
            packet = (packet << this->index_bits) | slot;
            bits += this->index_bits;
        }

        int code = Kangaroo32Utils::leadingCode(std::countl_zero(value ^ previous));
        int leading = Kangaroo32Utils::leadingZeros(code);
        packet = (packet << 3) | code;
        bits += 3;
        int center = std::max(0, 31 - leading - trailing);
        uint32_t difference = center > 0 ? (value ^ previous) >> (trailing + 1) : 0;
        // Concatenate the existing fields before writing; their layout is
        // unchanged. Very large reference indices may exceed one 64-bit packet.
        if (bits + center > 64) {
            this->__write(packet, bits);
            this->__write(difference, center);
        }
        else this->__write((packet << center) | difference, bits + center);
        // The current value's low trailing bits are zero and its next bit is
        // one. The decoder restores these directly, including for fallback
        // references whose XOR has fewer trailing zeros than the value.
    }

    void Compression::__flushRun() {
        if (this->run_count == 0) return;
        if (this->run_count > 1) {
            this->__write(0b110, 3);
            this->__writeRun(this->run_count);
        }

        int tag;
        uint32_t value;
        if (!this->cache_enabled) {
            value = Kangaroo32Utils::erase(this->run_value, tag);
        }
        else {
            uint32_t key = ((this->run_value ^ (this->run_value >> 16)) * 0x9e3779b1u) >> 24;
            if (this->cache_valid[key] && this->cache_values[key] == this->run_value) {
                value = this->cache_erased[key];
                tag = this->cache_tags[key];
                if (this->cache_queries < 256) ++this->cache_hits;
            }
            else {
                value = Kangaroo32Utils::erase(this->run_value, tag);
                this->cache_values[key] = this->run_value;
                this->cache_erased[key] = value;
                this->cache_tags[key] = tag;
                this->cache_valid[key] = true;
            }
            // Skip memoization on streams that do not reuse values. This
            // controls computation only; the erasure result is identical.
            if (this->cache_queries < 256 && ++this->cache_queries == 256 && this->cache_hits < 8) {
                this->cache_enabled = false;
            }
        }
        int trailing = std::countr_zero(value);
        this->__xorCompress(value, tag, trailing);
        this->__remember(value, this->run_count, trailing);
        this->run_count = 0;
    }

    void Compression::__compress(float value) {
        uint32_t raw = std::bit_cast<uint32_t>(value);
        if (this->run_count != 0 && raw != this->run_value) this->__flushRun();
        this->run_value = raw;
        ++this->run_count;
    }

    void Compression::compress(Univariate* data) {
        this->__compress(data->get_value());
        if (++this->count == this->buffer_size) {
            this->yield();
            this->count = 0;
        }
    }

    BinObj* Compression::serialize() {
        // Runs end at emission boundaries so each ordinary block reconstructs
        // exactly buffer_size samples while reference state survives the yield.
        this->__flushRun();
        BinObj* obj = this->bitstream;
        this->bitstream = new BinObj;
        return obj;
    }
    // End: compression

    // Begin: decompression
    void Decompression::initialize(int count, char** params) {
        if (count != 2) throw std::invalid_argument("Kangaroo requires size and window parameters");
        this->buffer_size = positiveParameter(params[0]);
        this->window_size = positiveParameter(params[1]);
        this->index_bits = std::bit_width(static_cast<uint32_t>(this->window_size - 1));
        this->stored_values.assign(this->window_size, 0);
        this->history_index = 0;
        this->history_slot = 0;
        this->repeat_count = 0;
        this->finished = false;
    }

    void Decompression::finalize() {
        std::vector<uint32_t>().swap(this->stored_values);
    }

    uint32_t Decompression::__readRun(BinObj* compress_data) {
        uint32_t length = 0;
        for (int shift = 0; shift <= 28; shift += 7) {
            uint32_t byte = compress_data->getBits(8);
            length |= (byte & 0x7fu) << shift;
            if ((byte & 0x80u) == 0) break;
        }

        return length;
    }

    uint32_t Decompression::__xorDecompress(BinObj* compress_data) {
        if (this->history_index == 0) return compress_data->getBits(32);

        int control = compress_data->getBits(1);
        int trailing;
        uint32_t previous;
        if (control != 0) {
            uint32_t reference = compress_data->getBits(this->index_bits);
            previous = this->stored_values[reference];
            trailing = std::countr_zero(previous);
        }
        else if (compress_data->getBits(1) == 0) {
            uint32_t reference = compress_data->getBits(this->index_bits);
            return this->stored_values[reference];
        }
        else {
            previous = this->stored_values[(this->history_index - 1) % this->window_size];
            trailing = compress_data->getBits(6);
        }

        int leading = Kangaroo32Utils::leadingZeros(compress_data->getBits(3));
        int center = std::max(0, 31 - leading - trailing);
        uint32_t difference = compress_data->getBits(center);
        if (trailing == 32) return 0;
        if (trailing == 31) return 0x80000000u;
        uint32_t upper = previous ^ (difference << (trailing + 1));
        return (upper & (UINT32_MAX << (trailing + 1))) | (1u << trailing);
    }

    float Decompression::__decompress(BinObj* compress_data) {
        if (this->finished) return Kangaroo::END_SIGNAL;
        if (this->repeat_count != 0) {
            --this->repeat_count;
            return std::bit_cast<float>(this->repeat_value);
        }

        int tag = compress_data->getBits(1);
        if (tag != 0) {
            tag = (tag << 1) | compress_data->getBits(1);;
            if (tag == 0b11) tag = (tag << 1) | compress_data->getBits(1);;
        }

        uint32_t repeats = 1;
        if (tag == 0b110) {
            repeats = this->__readRun(compress_data);
            if (repeats == 0) {
                this->finished = true;
                return Kangaroo::END_SIGNAL;
            }
            tag = compress_data->getBits(1);
            if (tag != 0) {
                tag = (tag << 1) | compress_data->getBits(1);;
                if (tag == 0b11) tag = (tag << 1) | compress_data->getBits(1);;
            }
        }

        uint32_t value = this->__xorDecompress(compress_data);
        this->repeat_value = Kangaroo32Utils::recover(value, tag);
        this->repeat_count = repeats - 1;
        remember(this->stored_values, this->history_index, this->history_slot, value, repeats);
        return std::bit_cast<float>(this->repeat_value);
    }

    CSVObj* Decompression::decompress(BinObj* compress_data) {
        CSVObj* base_obj = nullptr;
        CSVObj* prev_obj = nullptr;

        int count = 0;
        while (!this->finished && count++ < this->buffer_size) {
            float value = this->__decompress(compress_data);
            if (value == Kangaroo::END_SIGNAL) break;
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
