#include "floating-point/lossless.hpp"

namespace Elf {
    // Begin: Utils material
    double Elf32Utils::LOG_2_10 = log2(10);
    int Elf32Utils::mapSPGreater1[10] = {1, 10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000, 1000000000};
    float Elf32Utils::mapSPLess1[11] = {1, 0.1f, 0.01f, 0.001f, 0.0001f, 0.00001f, 0.000001f, 0.0000001f, 0.00000001f, 0.000000001f, 0.0000000001f};
    int Elf32Utils::f[21] = {0, 4, 7, 10, 14, 17, 20, 24, 27, 30, 34, 37, 40, 44, 47, 50, 54, 57, 60, 64, 67};

    float Elf32Utils::map10iP[21] = {1.0f, 1.0E1f, 1.0E2f, 1.0E3f, 1.0E4f, 1.0E5f, 1.0E6f, 1.0E7f, 1.0E8f, 1.0E9f, 
                        1.0E10f, 1.0E11f, 1.0E12f, 1.0E13f, 1.0E14f, 1.0E15f, 1.0E16f, 1.0E17f, 1.0E18f, 1.0E19f, 1.0E20f};

    float Elf32Utils::map10iN[21] = {1.0f, 1.0E-1f, 1.0E-2f, 1.0E-3f, 1.0E-4f, 1.0E-5f, 1.0E-6f, 1.0E-7f, 1.0E-8f, 1.0E-9f, 1.0E-10f, 
                        1.0E-11f, 1.0E-12f, 1.0E-13f, 1.0E-14f, 1.0E-15f, 1.0E-16f, 1.0E-17f, 1.0E-18f, 1.0E-19f, 1.0E-20f};


    int Elf32Utils::getFAlpha(int alpha) {
        if (alpha >= 21) return (int) std::ceil(alpha * LOG_2_10);
        else return f[alpha];
    }

    std::pair<int, int> Elf32Utils::getAlphaAndBetaStar(float v, int lastBetaStar) {
        if (v < 0) v = -v;

        std::pair<int, int> alphaAndBetaStar;
        std::pair<int, int> spAnd10iNFlag = getSPAnd10iNFlag(v);
        int beta = getSignificantCount(v, spAnd10iNFlag.first, lastBetaStar);
        alphaAndBetaStar.first = beta - spAnd10iNFlag.first - 1;
        alphaAndBetaStar.second = spAnd10iNFlag.second == 1 ? 0 : beta;
        
        return alphaAndBetaStar;
    }

    float Elf32Utils::roundUp(float v, int alpha) {
        float scale = get10iP(alpha);

        if (v < 0) return (float) (std::floor(v * scale) / scale); 
        else return (float) (std::ceil(v * scale) / scale);
    }

    int Elf32Utils::getSignificantCount(float v, int sp, int lastBetaStar) {
        int i;

        if (lastBetaStar != INT_MAX && lastBetaStar != 0) i = lastBetaStar - sp - 1 > 1 ? lastBetaStar - sp - 1 : 1;
        else if (lastBetaStar == INT_MAX) i = 8 - sp - 1;
        else if (sp >= 0) i = 1;
        else i = -sp;

        float temp = v * get10iP(i);
        int tempInt = (int) temp;
        while (tempInt != temp) {
            i++;
            temp = v * get10iP(i);
            tempInt = (int) temp;
        }

        // There are some bugs for those with high significand, e.g., 995455.44
        // So we should further check
        if (temp / get10iP(i) != v) return 8;
        else {
            while (i > 0 && tempInt % 10 == 0) {
                i--;
                tempInt = tempInt / 10;
            }
            return sp + i + 1;
        }
    }

    float Elf32Utils::get10iP(int i) {
        if (i >= 21) return stof("1.0E" + std::to_string(i));
        else return map10iP[i];
    }

    float Elf32Utils::get10iN(int i) {
        if (i >= 21) return stof("1.0E-" + std::to_string(i));
        else return map10iN[i];
    }

    int Elf32Utils::getSP(double v) {
        if (v >= 1) {
            int i = 0;
            while (i < 10 - 1) {
                if (v < mapSPGreater1[i + 1]) return i;
                i++;
            }
        } 
        else {
            int i = 1;
            while (i < 11) {
                if (v >= mapSPLess1[i]) return -i;
                i++;
            }
        }

        return (int) std::floor(log10(v));
    }

    std::pair<int, int> Elf32Utils::getSPAnd10iNFlag(double v) {
        std::pair<int, int> spAnd10iNFlag;
        if (v >= 1) {
            int i = 0;
            while (i < 10 - 1) {
                if (v < mapSPGreater1[i + 1]) {
                    spAnd10iNFlag.first = i;
                    return spAnd10iNFlag;
                }
                i++;
            }
        } 
        else {
            int i = 1;
            while (i < 11) {
                if (v >= mapSPLess1[i]) {
                    spAnd10iNFlag.first = -i;
                    spAnd10iNFlag.second = v == mapSPLess1[i] ? 1 : 0;
                    return spAnd10iNFlag;
                }
                i++;
            }
        }

        double log10v = log10(v);
        spAnd10iNFlag.first = (int) std::floor(log10v);
        spAnd10iNFlag.second = log10v == (long) log10v ? 1 : 0;

        return spAnd10iNFlag;
    }
    // End: Utils material

    // Begin: compression
    void Compression::initialize(int count, char** params) {
        this->buffer_size = atoi(params[0]);
        this->bitstream = new BinObj;
    }

    void Compression::finalize() {
        this->__compress(std::bit_cast<float>(Elf::END_SIGNAL));
        this->yield();

        delete this->bitstream;
        this->bitstream = nullptr;
        std::cout << "size: " << this->size << "\n";
    }

    void Compression::__xorCompression(uint32_t value) {
        if (this->first) {
            this->first = false;
            this->storedVal = value;
            int trailingZeros = std::countr_zero(value);
            this->bitstream->put(trailingZeros, 6);
            this->bitstream->put(this->storedVal >> trailingZeros, 32 - trailingZeros);
            this->size += 38 - trailingZeros;
        }
        else {
            uint32_t xored = this->storedVal ^ value;

            if (xored == 0) {
                // case 01
                this->bitstream->put(1, 2);
                this->size += 2;
            } 
            else {
                int leadingZeros = leadingRound[std::countl_zero(xored)];
                int trailingZeros = std::countr_zero(xored);

                if (leadingZeros == this->storedLeadingZeros && trailingZeros >= this->storedTrailingZeros) {
                    // case 00
                    int centerBits = 32 - this->storedLeadingZeros - this->storedTrailingZeros;
                    int len = 2 + centerBits;
                    if(len > 32) {
                        this->bitstream->put(0, 2);
                        this->bitstream->put(xored >> this->storedTrailingZeros, centerBits);
                    } 
                    else this->bitstream->put(xored >> this->storedTrailingZeros, len);
                    this->size += len;
                } 
                else {
                    this->storedLeadingZeros = leadingZeros;
                    this->storedTrailingZeros = trailingZeros;
                    int centerBits = 32 - this->storedLeadingZeros - this->storedTrailingZeros;
                    if (centerBits <= 8) {
                        // case 10
                        this->bitstream->put((((0x2 << 3) | leadingRepresentation[storedLeadingZeros]) << 3) | (centerBits & 0x7), 8);
                        this->bitstream->put(xored >> storedTrailingZeros, centerBits);
                        this->size += 8 + centerBits;
                    } 
                    else {
                        // case 11
                        this->bitstream->put((((0x3 << 3) | leadingRepresentation[this->storedLeadingZeros]) << 5) | (centerBits & 0x1f), 10);
                        this->bitstream->put(xored >> this->storedTrailingZeros, centerBits);
                        this->size += 10 + centerBits;
                    }
                }
                this->storedVal = value;
            }
        }
    }

    void Compression::__compress(float value) {
        uint32_t vInt = std::bit_cast<uint32_t>(value);
        uint32_t vPrimeInt;

        if (value == 0.0 || value == INFINITY) {
            this->bitstream->put(2, 2);
            this->size += 2;
            vPrimeInt = vInt;
        } else {
            // C1: v is a normal or subnormal
            std::pair<int, int> alphaAndBetaStar = Elf32Utils::getAlphaAndBetaStar(value, this->lastBetaStar);
            int e = (vInt >> 23) & 0xff;
            int gAlpha = Elf32Utils::getFAlpha(alphaAndBetaStar.first) + e - 127;
            int eraseBits = 23 - gAlpha;
            int mask = 0xffffffff << eraseBits;
            int delta = (~mask) & vInt;
            if (delta != 0 && eraseBits > 3) {
                if(alphaAndBetaStar.second == this->lastBetaStar) {
                    this->bitstream->put(0, 1);
                    this->size += 1;
                } 
                else {
                    this->bitstream->put(alphaAndBetaStar.second | 0x18, 5);
                    this->lastBetaStar = alphaAndBetaStar.second;
                    this->size += 5;
                }
                vPrimeInt = mask & vInt;
            } 
            else {
                this->bitstream->put(2, 2);
                this->size += 2;
                vPrimeInt = vInt;
            }
        }
        
        this->__xorCompression(vPrimeInt);
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

    float Decompression::__recoverVByBetaStar(BinObj* compress_data) {
        float v;
        float vPrime = this->__xorDecompress(compress_data);
        int sp = Elf32Utils::getSP(abs(vPrime));
        
        if (this->lastBetaStar == 0) {
            v = Elf32Utils::get10iN(-sp - 1);
            if (vPrime < 0) v = -v;
        } 
        else {
            int alpha = this->lastBetaStar - sp - 1;
            v = Elf32Utils::roundUp(vPrime, alpha);
        }
        
        return v;
    }

    float Decompression::__xorDecompress(BinObj* compress_data) {
        uint32_t value;

        if (this->first) {
            this->first = false;
            int trailingZeros = compress_data->getBits(6);
            value = compress_data->getBits(32 - trailingZeros) << trailingZeros;
        } 
        else {
            int centerBits, leadAndCenter;
            int flag = compress_data->getBits(2);
            switch (flag) {
                case 3:
                    // case 11
                    leadAndCenter = compress_data->getBits(8);
                    this->storedLeadingZeros = leadingRepresentation[leadAndCenter >> 5];
                    centerBits = leadAndCenter & 0x1f;
                    if(centerBits == 0) centerBits = 32;
                    
                    this->storedTrailingZeros = 32 - this->storedLeadingZeros - centerBits;
                    value = compress_data->getBits(centerBits) << this->storedTrailingZeros;
                    value = this->storedVal ^ value;
                    break;
                case 2:
                    // case 10
                    leadAndCenter = compress_data->getBits(6); 
                    this->storedLeadingZeros = leadingRepresentation[leadAndCenter >> 3];
                    centerBits = leadAndCenter & 0x7;
                    if(centerBits == 0) centerBits = 8;

                    this->storedTrailingZeros = 32 - this->storedLeadingZeros - centerBits;
                    value = compress_data->getBits(centerBits) << this->storedTrailingZeros;
                    value = this->storedVal ^ value;
                    break;
                case 1:
                    // case 01, we do nothing, the same value as before
                    value = this->storedVal;
                    break;
                default:
                    // case 00
                    centerBits = 32 - this->storedLeadingZeros - this->storedTrailingZeros;
                    value = compress_data->getBits(centerBits) << this->storedTrailingZeros;
                    value = this->storedVal ^ value;
                    break;
            }
        }

        if (value != Elf::END_SIGNAL) {
            this->storedVal = value;
            return std::bit_cast<float>(value);
        }
        
        return Elf::END_SIGNAL;
    }

    float Decompression::__decompress(BinObj* compress_data) {
        if (compress_data->getBits(1) == 0) {
            return this->__recoverVByBetaStar(compress_data);               // case 0
        } 
        else if (compress_data->getBits(1) == 0) {
            return this->__xorDecompress(compress_data);                    // case 10
        } 
        else {
            this->lastBetaStar = compress_data->getBits(3);          // case 11
            return this->__recoverVByBetaStar(compress_data);
        }
    }

    CSVObj* Decompression::decompress(BinObj* compress_data) {
        CSVObj* base_obj = nullptr;
        CSVObj* prev_obj = nullptr;

        int count = 0;
        while (count++ != this->buffer_size) {
            float value = this->__decompress(compress_data);
            if (value == Elf::END_SIGNAL) break;
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