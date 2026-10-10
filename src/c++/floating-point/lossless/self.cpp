#include "floating-point/lossless.hpp"

namespace SElf {
    // Begin: Huffman encoding material
    Code::Code() {
        this->code = 0;
        this->length = 0;
    }

    Code::Code(int code, int length) {
        this->code = code;
        this->length = length;
    }

    Node::~Node() {
        if (children != nullptr) {
            delete children[0];
            delete children[1];
            delete[] children;
        }
    }

    Node::Node() : data(-1), frequency(0), height(0) {
        this->children = new Node*[2];
        this->children[0] = nullptr;
        this->children[1] = nullptr;
    }

    Node::Node(int data) : data(data), frequency(0), height(0) {
        this->children = new Node*[2];
        this->children[0] = nullptr;
        this->children[1] = nullptr;
    }

    Node::Node(int data, int frequency, int height)
         : data(data), frequency(frequency), height(height) 
    {   
        this->children = new Node*[2];
        this->children[0] = nullptr;
        this->children[1] = nullptr;
    }

    std::vector<Code> HuffmanEncode::getHuffmanCodes(int* frequencies) {
        std::vector<Code> huffmanCodes(9);

        std::priority_queue<Node*, std::vector<Node*>, CompareNode> pq;
        // Construct priority queue
        for (int i = 0; i < 9; i++) {
            Node* node = new Node(i, frequencies[i], 0);
            pq.push(node);
        }

        // Construct Huffman tree
        while (pq.size() > 1) {
            Node* left = pq.top(); pq.pop();
            Node* right = pq.top(); pq.pop();
            int height = (left->height > right->height) ? left->height + 1 : right->height + 1;

            Node* newNode = new Node(-1, left->frequency + right->frequency, height);
            newNode->children[0] = left;
            newNode->children[1] = right;

            pq.push(newNode);
        }
        Node* root = pq.top();
        
        generateHuffmanCodes(huffmanCodes, root, 0, 0);

        delete root;
        return huffmanCodes;
    }

    void HuffmanEncode::freeTree(Node* root) {
        if (root == nullptr) return;
        freeTree(root->children[0]);
        freeTree(root->children[1]);
        delete root;
    }

    void HuffmanEncode::generateHuffmanCodes(std::vector<Code>& huffmanCodes, Node* root, int code, int length) {
        if (root != nullptr) {
            if (root->data >= 0) {
                huffmanCodes[root->data] = Code(code, length);
            }

            generateHuffmanCodes(huffmanCodes, root->children[0], code << 1, length + 1);
            generateHuffmanCodes(huffmanCodes, root->children[1], (code << 1) | 1, length + 1);
        }
    }

    bool HuffmanEncode::nodeCompare(Node& a, Node& b) {
        if (a.frequency != b.frequency)
            return a.frequency > b.frequency;  // min-heap
        return a.height > b.height;
    }

    Node* HuffmanEncode::buildHuffmanTree(std::vector<Code> huffmanCodes) {
        Node* root = new Node(-1);
        Node* curNode = root;

        for (int value = 0; value < 9; value++) {
            int code = huffmanCodes[value].code;
            int length = huffmanCodes[value].length - 1;

            while (length >= 0) {
                int signal = (code >> length) & 1;

                if (curNode->children[signal] == nullptr) {
                    curNode->children[signal] = new Node(-1);
                }

                curNode = curNode->children[signal];
                --length;
            }

            curNode->data = value;
            curNode = root;
        }

        return root;
    }
    // End: Huffman encoding material

    // Begin: Utils material
    double SElf32Utils::LOG_2_10 = log2(10);
    int SElf32Utils::mapSPGreater1[10] = {1, 10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000, 1000000000};
    float SElf32Utils::mapSPLess1[11] = {1, 0.1f, 0.01f, 0.001f, 0.0001f, 0.00001f, 0.000001f, 0.0000001f, 0.00000001f, 0.000000001f, 0.0000000001f};
    int SElf32Utils::f[21] = {0, 4, 7, 10, 14, 17, 20, 24, 27, 30, 34, 37, 40, 44, 47, 50, 54, 57, 60, 64, 67};

    float SElf32Utils::map10iP[21] = {1.0f, 1.0E1f, 1.0E2f, 1.0E3f, 1.0E4f, 1.0E5f, 1.0E6f, 1.0E7f, 1.0E8f, 1.0E9f, 
                        1.0E10f, 1.0E11f, 1.0E12f, 1.0E13f, 1.0E14f, 1.0E15f, 1.0E16f, 1.0E17f, 1.0E18f, 1.0E19f, 1.0E20f};

    float SElf32Utils::map10iN[21] = {1.0f, 1.0E-1f, 1.0E-2f, 1.0E-3f, 1.0E-4f, 1.0E-5f, 1.0E-6f, 1.0E-7f, 1.0E-8f, 1.0E-9f, 1.0E-10f, 
                        1.0E-11f, 1.0E-12f, 1.0E-13f, 1.0E-14f, 1.0E-15f, 1.0E-16f, 1.0E-17f, 1.0E-18f, 1.0E-19f, 1.0E-20f};

    int SElf32Utils::getFAlpha(int alpha) {
        if (alpha >= 21) return (int) std::ceil(alpha * LOG_2_10);
        else return f[alpha];
    }

    std::pair<int, int> SElf32Utils::getAlphaAndBetaStar(float v, int lastBetaStar) {
        if (v < 0) v = -v;

        std::pair<int, int> alphaAndBetaStar;
        std::pair<int, int> spAnd10iNFlag = SElf32Utils::getSPAnd10iNFlag(v);
        int beta = SElf32Utils::getSignificantCount(v, spAnd10iNFlag.first, lastBetaStar);
        alphaAndBetaStar.first = beta - spAnd10iNFlag.first - 1;
        alphaAndBetaStar.second = spAnd10iNFlag.second == 1 ? 0 : beta;
        
        return alphaAndBetaStar;
    }

    float SElf32Utils::roundUp(float v, int alpha) {
        float scale = SElf32Utils::get10iP(alpha);

        if (v < 0) return (float) (std::floor(v * scale) / scale); 
        else return (float) (std::ceil(v * scale) / scale);
    }

    int SElf32Utils::getSignificantCount(float v, int sp, int lastBetaStar) {
        int i;
        if (lastBetaStar != INT_MAX && lastBetaStar != 0) i = lastBetaStar - sp - 1 > 1 ? lastBetaStar - sp - 1 : 1;
        else if (sp >= 0) i = 1;
        else i = -sp;

        float temp = v * SElf32Utils::get10iP(i);
        int tempInt = (int) temp;
        while (tempInt != temp) {
            i++;
            temp = v * SElf32Utils::get10iP(i);
            tempInt = (int) temp;
        }

        // There are some bugs for those with high significand, e.g., 995455.44
        // So we should further check
        if (temp / SElf32Utils::get10iP(i) != v) return 8;
        else {
            while (i > 0 && tempInt % 10 == 0) {
                i--;
                tempInt = tempInt / 10;
            }
            return sp + i + 1;
        }
    }

    float SElf32Utils::get10iP(int i) {
        if (i >= 21) return stof("1.0E" + std::to_string(i));
        else return map10iP[i];
    }

    float SElf32Utils::get10iN(int i) {
        if (i >= 21) return stof("1.0E-" + std::to_string(i));
        else return map10iN[i];
    }

    int SElf32Utils::getSP(float v) {
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

    std::pair<int, int> SElf32Utils::getSPAnd10iNFlag(float v) {
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
        this->__compress(std::bit_cast<float>(SElf::END_SIGNAL));
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
            this->size += 6;
            if (trailingZeros < 32) {
                this->bitstream->put(this->storedVal >> (trailingZeros + 1), 31 - trailingZeros);
                this->size += 31 - trailingZeros;
            }
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
                int trailingZeros = trailingRound[std::countr_zero(xored)];

                if (leadingZeros >= this->storedLeadingZeros && trailingZeros >= this->storedTrailingZeros &&
                    (leadingZeros - this->storedLeadingZeros) + (trailingZeros - this->storedTrailingZeros) < 1 + leadingBitsPerValue + trailingBitsPerValue) 
                {
                    // case 1
                    int centerBits = 32 - this->storedLeadingZeros - this->storedTrailingZeros;
                    int len = 1 + centerBits;

                    if (len > 32) {
                        this->bitstream->put(1, 1);
                        this->bitstream->put(xored >> this->storedTrailingZeros, centerBits);
                    } 
                    else {
                        uint32_t packed = (1 << centerBits) | (xored >> this->storedTrailingZeros);
                        this->bitstream->put(packed, len);
                    }
                    this->size += len;
                } 
                else {
                    this->storedLeadingZeros = leadingZeros;
                    this->storedTrailingZeros = trailingZeros;
                    int centerBits = 32 - this->storedLeadingZeros - this->storedTrailingZeros;

                    // case 00
                    int len = 2 + leadingBitsPerValue + trailingBitsPerValue + centerBits;

                    uint32_t header =
                        (leadingRepresentation[this->storedLeadingZeros] << trailingBitsPerValue) |
                        trailingRepresentation[this->storedTrailingZeros];

                    if (len > 32) {
                        this->bitstream->put(header, 2 + leadingBitsPerValue + trailingBitsPerValue);
                        this->bitstream->put(xored >> this->storedTrailingZeros, centerBits);
                    } 
                    else {
                        uint32_t packed = (header << centerBits) | (xored >> this->storedTrailingZeros);
                        this->bitstream->put(packed, len);
                    }

                    this->size += len;
                }

                this->storedVal = value;
            }
        }
    }

    void Compression::__addValueFirst(float value) {
        uint32_t vInt = std::bit_cast<uint32_t>(value);
        uint32_t vPrimeInt;

        if (value == 0.0 || value == INFINITY) {
            this->bitstream->put(2, 2); // case 10
            this->size += 2;
            vPrimeInt = vInt;
            this->frequency[8]++;
        } 
        else {
            // C1: v is a normal or subnormal
            std::pair<int, int> alphaAndBetaStar = SElf32Utils::getAlphaAndBetaStar(value, this->lastBetaStar);
            int e = ((int) (vInt >> 23)) & 0xff;
            int gAlpha = SElf32Utils::getFAlpha(alphaAndBetaStar.first) + e - 127;
            int eraseBits = 23 - gAlpha;
            long mask = 0xffffffff << eraseBits;
            long delta = (~mask) & vInt;
            if (delta != 0 && eraseBits > 3) {  // C2
                if (alphaAndBetaStar.second == this->lastBetaStar) {
                    this->bitstream->put(false, 1);    // case 0
                    this->size += 1;
                } 
                else {
                    this->bitstream->put(alphaAndBetaStar.second | 0x18, 5);  // case 11, 2 + 4 = 6
                    this->size += 5;
                    this->lastBetaStar = alphaAndBetaStar.second;
                }
                vPrimeInt = mask & vInt;
                this->frequency[alphaAndBetaStar.second]++;
            } 
            else {
                this->bitstream->put(2, 2); // case 10
                this->size += 2;
                vPrimeInt = vInt;
                this->frequency[8]++;
            }
        }
        
        this->__xorCompression(vPrimeInt);
    }

    void Compression::__addValueHuffman(float value) {
        uint32_t vInt = std::bit_cast<uint32_t>(value);
        uint32_t vPrimeInt;

        if (value == 0.0 || value == INFINITY) {
            this->bitstream->put(this->huffmanCode[8].code, this->huffmanCode[8].length); // not erase
            this->size += this->huffmanCode[8].length;
            vPrimeInt = vInt;
            this->frequency[8]++;
        } 
        else {
            // C1: v is a normal or subnormal
            std::pair<int, int> alphaAndBetaStar = SElf32Utils::getAlphaAndBetaStar(value, this->lastBetaStar);
            int e = ((int) (vInt >> 23)) & 0xff;
            int gAlpha = SElf32Utils::getFAlpha(alphaAndBetaStar.first) + e - 127;
            int eraseBits = 23 - gAlpha;
            long mask = 0xffffffff << eraseBits;
            long delta = (~mask) & vInt;
            if (delta != 0 && eraseBits > 3) {  // C2
                this->bitstream->put(this->huffmanCode[alphaAndBetaStar.second].code, this->huffmanCode[alphaAndBetaStar.second].length);  // case 11, 2 + 4 = 6
                this->size += this->huffmanCode[alphaAndBetaStar.second].length;
                this->lastBetaStar = alphaAndBetaStar.second;
                vPrimeInt = mask & vInt;
                this->frequency[alphaAndBetaStar.second]++;
            } 
            else {
                this->bitstream->put(this->huffmanCode[8].code, this->huffmanCode[8].length); // not erase
                this->size += this->huffmanCode[8].length;
                vPrimeInt = vInt;
                this->frequency[8]++;
            }
        }

        this->__xorCompression(vPrimeInt);
    }

    void Compression::__compress(float value) {
        if (!this->isFirstBlock) {
            this->__addValueHuffman(value);
        } else {
            this->__addValueFirst(value);
        }
    }

    void Compression::compress(Univariate* data) {
        float value = data->get_value();
        this->__compress(value);
        
        if (++this->count == this->buffer_size) {
            this->yield();
            this->isFirstBlock = false;

            this->huffmanCode = HuffmanEncode::getHuffmanCodes(this->frequency);
            
            for (int i=0; i<9; i++) this->frequency[i] = 0;
            this->lastBetaStar = INT_MAX;
            this->count = 0;
            this->first = true;
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
        
        int sp = SElf32Utils::getSP(std::abs(vPrime));
        if (this->lastBetaStar == 0) {
            v = SElf32Utils::get10iN(-sp - 1);
            if (vPrime < 0) {
                v = -v;
            }
        } 
        else {
            int alpha = this->lastBetaStar - sp - 1;
            v = SElf32Utils::roundUp(vPrime, alpha);
        }
        return v;
    }

    float Decompression::__xorDecompress(BinObj* compress_data) {
        uint32_t value;

        if (this->first) {
            this->first = false;
            int trailingZeros = compress_data->getBits(6);
            value = trailingZeros >= 32 ? 0 :
                ((compress_data->getBits(31 - trailingZeros) << 1) + 1) << trailingZeros;        
        }
        else {
            int centerBits;

            if (compress_data->getBits(1) == 1) {
                // case 1
                centerBits = 32 - this->storedLeadingZeros - this->storedTrailingZeros;
                value = compress_data->getBits(centerBits) << this->storedTrailingZeros;
                value = this->storedVal ^ value;
            } 
            else if (compress_data->getBits(1) == 0) {
                // case 00
                int leadAndTrail = compress_data->getBits(leadingBitsPerValue + trailingBitsPerValue);
                int lead = leadAndTrail >> trailingBitsPerValue;
                int trail = ~(0xffff << trailingBitsPerValue) & leadAndTrail;
                this->storedLeadingZeros = leadingRepresentation[lead];
                this->storedTrailingZeros = trailingRepresentation[trail];
                centerBits = 32 - this->storedLeadingZeros - this->storedTrailingZeros;

                value = compress_data->getBits(centerBits) << storedTrailingZeros;
                value = this->storedVal ^ value;
            }
            else value = this->storedVal;
        }
        
        if (value != SElf::END_SIGNAL) {
            this->storedVal = value;
            return std::bit_cast<float>(value);
        }
        
        return SElf::END_SIGNAL;
    }

    float Decompression::__decompress(BinObj* compress_data) {
        if (this->isFirstBlock) {
            float v;
            if (compress_data->getBits(1) == 0) {
                v = this->__recoverVByBetaStar(compress_data);               // case 0
                this->frequency[this->lastBetaStar]++;
            } else if (compress_data->getBits(1) == 0) {
                v = this->__xorDecompress(compress_data);        // case 10
                this->frequency[8]++;
            } else {
                this->lastBetaStar = compress_data->getBits(3);          // case 11
                v = this->__recoverVByBetaStar(compress_data);
                this->frequency[lastBetaStar]++;
            }
            return v;
        }
        else {
            float v;
            Node* current = this->root;

            while (true) {
                current = current->children[compress_data->getBits(1)];
                if (current->data >= 0) {
                    if (current->data != 8) {
                        this->lastBetaStar = current->data;
                        v = this->__recoverVByBetaStar(compress_data);
                        this->frequency[this->lastBetaStar]++;
                    } 
                    else {
                        v = this->__xorDecompress(compress_data);
                        this->frequency[8]++;
                    }
                    break;
                }
            }
            return v;
        }
    }

    CSVObj* Decompression::decompress(BinObj* compress_data) {
        CSVObj* base_obj = nullptr;
        CSVObj* prev_obj = nullptr;

        int count = 0;
        while (count++ != this->buffer_size) {
            float value = this->__decompress(compress_data);
            if (value == SElf::END_SIGNAL) break;
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

        // this->frequency[8]--;
        delete this->root;

        std::vector<Code> huffmanCode = HuffmanEncode::getHuffmanCodes(this->frequency);
        this->root = HuffmanEncode::buildHuffmanTree(huffmanCode);
        for (int i=0 ;i<9; i++) this->frequency[i] = 0;

        this->lastBetaStar = INT_MAX;
        this->isFirstBlock = false;
        this->first = true;
        
        compress_data->flushBits();
        return base_obj;    
    }
    // End: decompression
};