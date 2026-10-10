#include "base-c.hpp"
#include "base-d.hpp"

namespace Gorilla {
    // Source paper: Gorilla: A Fast, Scalable, In-Memory Time Series Database
    // Source path: src/floating-point/lossless/gorilla.cpp

    inline float END_SIGNAL = INFINITY;

    class Compression : public BaseCompression {
        private:
            long size = 0;
            int count = 0;
            int buffer_size = 0;
            BinObj* bitstream = nullptr;

            bool first = true;
            uint32_t storedVal;
            uint32_t storedLeadingZeros = INT_MAX;
            uint32_t storedTrailingZeros = INT_MAX;
            
            void __compress(uint32_t value);

        protected:
            void compress(Univariate* data) override;
            BinObj* serialize() override;

        public:
            Compression(std::string output) : BaseCompression(output) {}            
            void initialize(int count, char** params) override;
            void finalize() override;
    };

    class Decompression : public BaseDecompression {
        private:
            int buffer_size = 0;
            bool first = true;
            uint32_t storedVal;
            uint32_t storedLeadingZeros = INT_MAX;
            uint32_t storedTrailingZeros = INT_MAX;

            float __decompress(BinObj* compress_data);

        protected:
            CSVObj* decompress(BinObj* compress_data) override;

        public:
            Decompression(std::string output, int interval, std::time_t basetime) : BaseDecompression(output, interval, basetime) {} 
            void initialize(int count, char** params) override;
            void finalize() override;
    };
};

namespace Chimp {
    // Source paper: Chimp: Efficient Lossless Floating Point Compression for Time Series Databases
    // Source path: src/floating-point/lossless/chimp.cpp

    inline uint32_t END_SIGNAL = std::bit_cast<uint32_t>(INFINITY);

    class Compression : public BaseCompression {
        private:
            long size = 0;    
            int count = 0;
            int buffer_size = 0;
            std::string mode = "";
            BinObj* bitstream = nullptr;
            void (Compression::*__compress) (uint32_t);

            bool first = true;
            int setLsb;
            int* indices;
            int index = 0;
            int current = 0;
            int threshold;
            int previousValues;
            int previousValuesLog2;
            int storedLeadingZeros = INT_MAX;
            uint32_t storedVal;
            uint32_t* storedValues;
            
            short leadingRepresentation[32] = {
                0, 0, 0, 0, 0, 0, 1, 1,
                1, 1, 2, 2, 3, 3, 4, 4,
                5, 5, 6, 6, 7, 7, 7, 7,
                7, 7, 7, 7, 7, 7, 7, 7
            };

            short leadingRound[32] = {
                0, 0, 0, 0, 0, 0, 6, 6,
                6, 6, 10, 10, 12, 12, 14, 14,
                16, 16, 18, 18, 20, 20, 20, 20,
                20, 20, 20, 20, 20, 20, 20, 20
            };

            void __basic_compress(uint32_t value);
            void __128_compress(uint32_t value);

        protected:
            void compress(Univariate* data) override;
            BinObj* serialize() override;

        public:
            Compression(std::string output) : BaseCompression(output) {}            
            void initialize(int count, char** params) override;
            void finalize() override;
    };

    class Decompression : public BaseDecompression {
        private:
            int buffer_size = 0;
            std::string mode = "";
            float (Decompression::*__decompress) (BinObj*);

            bool first = true;
            int current = -1;
            int previousValues;
            int previousValuesLog2;
            int storedLeadingZeros = INT_MAX;
            int storedTrailingZeros = 0;
            uint32_t storedVal;
            uint32_t* storedValues;
            
            short leadingRepresentation[8] = {0, 6, 10, 12, 14, 16, 18, 20};

            float __basic_decompress(BinObj* compress_data);
            float __128_decompress(BinObj* compress_data);

        protected:
            CSVObj* decompress(BinObj* compress_data) override;

        public:
            Decompression(std::string output, int interval, std::time_t basetime) : BaseDecompression(output, interval, basetime) {} 
            void initialize(int count, char** params) override;
            void finalize() override;
    };
};

namespace Elf {
    // Source paper: Elf: Erasing-based Lossless Floating-Point Compression
    // Source path: src/floating-point/lossless/elf.cpp

    inline uint32_t END_SIGNAL = std::bit_cast<uint32_t>(INFINITY);

    class Elf32Utils {
        private:
            static int f[21];
            static float map10iP[21];
            static float map10iN[21];
            static int mapSPGreater1[10];
            static float mapSPLess1[11];
            static double LOG_2_10;

        public:
            static int getFAlpha(int alpha);
            static std::pair<int, int> getAlphaAndBetaStar(float v, int lastBetaStar);
            static float roundUp(float v, int alpha);
            static int getSignificantCount(float v, int sp, int lastBetaStar);
            static float get10iP(int i);
            static float get10iN(int i);
            static int getSP(double v);
            static std::pair<int, int> getSPAnd10iNFlag(double v);
    };

    class Compression : public BaseCompression {
        private:
            long size = 0;
            int count = 0;
            int buffer_size = 0;
            BinObj* bitstream = nullptr;

            bool first = true;
            uint32_t storedVal = 0;
            int lastBetaStar = INT_MAX;
            int storedLeadingZeros = INT_MAX;
            int storedTrailingZeros = INT_MAX;

            short leadingRepresentation[32] = {
                0, 0, 0, 0, 0, 0, 1, 1,
                1, 1, 2, 2, 3, 3, 4, 4,
                5, 5, 6, 6, 7, 7, 7, 7,
                7, 7, 7, 7, 7, 7, 7, 7
            };

            short leadingRound[32] = {
                0, 0, 0, 0, 0, 0, 6, 6,
                6, 6, 10, 10, 12, 12, 14, 14,
                16, 16, 18, 18, 20, 20, 20, 20,
                20, 20, 20, 20, 20, 20, 20, 20
            };

            void __xorCompression(uint32_t value);
            void __compress(float value);

        protected:
            void compress(Univariate* data) override;
            BinObj* serialize() override;

        public:
            Compression(std::string output) : BaseCompression(output) {}            
            void initialize(int count, char** params) override;
            void finalize() override;
    };

    class Decompression : public BaseDecompression {
        private:
            bool first = true;
            int buffer_size = 0;
            uint32_t storedVal = 0;
            int lastBetaStar = INT_MAX;
            int storedLeadingZeros = INT_MAX;
            int storedTrailingZeros = INT_MAX;
            short leadingRepresentation[8] = {0, 6, 10, 12, 14, 16, 18, 20};

            float __recoverVByBetaStar(BinObj* compress_data);
            float __xorDecompress(BinObj* compress_data);
            float __decompress(BinObj* compress_data);

        protected:
            CSVObj* decompress(BinObj* compress_data) override;

        public:
            Decompression(std::string output, int interval, std::time_t basetime) : BaseDecompression(output, interval, basetime) {} 
            void initialize(int count, char** params) override;
            void finalize() override;
    };
};

namespace SElf {
    // Source paper: Adaptive Encoding Strategies for Lossless Floating-Point Compression
    // Source path: src/floating-point/lossless/self.cpp

    inline uint32_t END_SIGNAL = std::bit_cast<uint32_t>(INFINITY);

    struct Code {
        int code;
        int length;

        Code();
        Code(int code, int length);
    };

    struct Node {
        int data;
        int frequency;
        Node** children;
        int height; // leaf node is 0
        
        Node();
        ~Node();
        Node(int data);
        Node(int data, int frequency, int height);
    };

    struct CompareNode {
        bool operator()(const Node* a, const Node* b) const {
            if (a->frequency != b->frequency)
                return a->frequency > b->frequency;

            return a->height > b->height;
        }
    };

    class HuffmanEncode {
        private:
            static void freeTree(Node* root);
            static void generateHuffmanCodes(std::vector<Code>& huffmanCodes, Node* root, int code, int length);
            static bool nodeCompare(Node& a, Node& b);

        public:
            static std::vector<Code> getHuffmanCodes(int* frequencies);
            static Node* buildHuffmanTree(std::vector<Code> huffmanCodes);
    };

    class SElf32Utils {
        private:
            static int f[21];
            static float map10iP[21];
            static float map10iN[21];
            static int mapSPGreater1[10];
            static float mapSPLess1[11];
            static double LOG_2_10;

        public:
            static int getFAlpha(int alpha);
            static std::pair<int, int> getAlphaAndBetaStar(float v, int lastBetaStar);
            static float roundUp(float v, int alpha);
            static int getSignificantCount(float v, int sp, int lastBetaStar);
            static float get10iP(int i);
            static float get10iN(int i);
            static int getSP(float v);
            static std::pair<int, int> getSPAnd10iNFlag(float v);
    };

    class Compression : public BaseCompression {
        private:
            long size = 0;
            int count = 0;
            int buffer_size = 0;
            BinObj* bitstream = nullptr;

            short leadingRepresentation[32] = {
                    0, 0, 0, 0, 0, 0, 0, 0,
                    1, 1, 1, 1, 2, 2, 2, 2,
                    3, 3, 3, 3, 3, 3, 3, 3,
                    3, 3, 3, 3, 3, 3, 3, 3,
            };
            short leadingRound[32] = {
                    0, 0, 0, 0, 0, 0, 0, 0,
                    8, 8, 8, 8, 12, 12, 12, 12,
                    16, 16, 16, 16, 16, 16, 16, 16,
                    16, 16, 16, 16, 16, 16, 16, 16
            };
            short trailingRepresentation[32] = {
                    0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0,
                    1, 1, 1, 1, 1, 1, 1, 1,
                    1, 1, 1, 1, 1, 1, 1, 1
            };
            short trailingRound[32] = {
                    0, 0, 0, 0, 0, 0, 0, 0,
                    0, 0, 0, 0, 0, 0, 0, 0,
                    16, 16, 16, 16, 16, 16, 16, 16,
                    16, 16, 16, 16, 16, 16, 16, 16
            };

            int storedLeadingZeros = INT_MAX;
            int storedTrailingZeros = INT_MAX;
            
            bool first = true;
            uint32_t storedVal = 0;
            int leadingBitsPerValue = 2;
            int trailingBitsPerValue = 1;
                
            bool isFirstBlock = true;
            int lastBetaStar = INT_MAX;
            int frequency[9] = { 0, 0, 0, 0, 0, 0, 0, 0, 0}; 
            std::vector<Code> huffmanCode;

            void __xorCompression(uint32_t value);
            void __addValueFirst(float value);
            void __addValueHuffman(float value);
            void __compress(float value);

        protected:
            void compress(Univariate* data) override;
            BinObj* serialize() override;

        public:
            Compression(std::string output) : BaseCompression(output) {}            
            void initialize(int count, char** params) override;
            void finalize() override;
    };

    class Decompression : public BaseDecompression {
        private:
            int buffer_size = 0;
            int lastBetaStar = INT_MAX;
            bool isFirstBlock = true;
            Node* root;
            int frequency[9] = { 0, 0, 0, 0, 0, 0, 0, 0, 0};
            
            bool first = true;
            uint32_t storedVal = 0;
            int storedLeadingZeros = INT_MAX;
            int storedTrailingZeros = INT_MAX;
            short leadingRepresentation[8] = {0, 8, 12, 16};
            short trailingRepresentation[8] = {0, 16};
            int leadingBitsPerValue = 2;
            int trailingBitsPerValue = 1;

            float __recoverVByBetaStar(BinObj* compress_data);
            float __xorDecompress(BinObj* compress_data);
            float __decompress(BinObj* compress_data);

        protected:
            CSVObj* decompress(BinObj* compress_data) override;

        public:
            Decompression(std::string output, int interval, std::time_t basetime) : BaseDecompression(output, interval, basetime) {} 
            void initialize(int count, char** params) override;
            void finalize() override;
    };
};

namespace DeXOR {
    // Source paper: DeXOR: Enabling XOR in Decimal Space for Streaming Lossless Compression of Floating-point Data
    // Source path: src/c++/floating-point/lossless/dexor.cpp

    inline float END_SIGNAL = INFINITY;

    class DeXORTools {
        private:
            // The encoder reuses consecutive powers during the author's search.
            friend class Compression;
            friend class Decompression;
            static const int cost[7];
            static const double P10[38];
            static const double I10[38];
            static bool __isEnd(float value, int q);

        public:
            static constexpr int MIN_Q = -20;
            static constexpr int MAX_Q = 11;
            // Float32: suffix lengths 0..6 use at most 20 bits; delta needs 3 bits.
            static constexpr int MAX_DELTA = 6;
            static constexpr int DELTA_BITS = 3;

            static double getP10(int pow);
            static int getP2(int pow);
            static int getEnd(float value, int last_end);
            static int decimalBits(int dp);
            static double truncate(double value);
            static uint32_t segment(uint32_t v, int st, int ed);
    };

    class Compression : public BaseCompression {
        private:
            long size = 0;
            int count = 0;
            int buffer_size = 0;
            BinObj* bitstream = nullptr;

            float previous_value = 0;
            int previous_q = 0;
            int previous_delta = 0;
            double previous_alpha = 0;

            int previous_exp = 127;
            int EL = 1;
            int contract_step = 0;
            int rho = 8;
            bool __cannotEncode(uint32_t raw, double alpha, double beta, double pow, int delta);
            void __exceptionHandle(uint32_t raw, bool end);

        protected:
            void __compress(float value, bool end = false);
            void compress(Univariate* data) override;
            BinObj* serialize() override;

        public:
            Compression(std::string output) : BaseCompression(output) {}
            void initialize(int count, char** params) override;
            void finalize() override;
    };

    class Decompression : public BaseDecompression {
        private:
            int buffer_size = 0;
            bool finished = false;

            float previous_value = 0;
            int previous_q = 0;
            int previous_delta = 0;

            int previous_exp = 127;
            int EL = 1;
            int contract_step = 0;

            double previous_alpha = 0;

            int rho = 8;

            float __exceptionDecode(BinObj* compress_data);

        protected:
            float __decompress(BinObj* compress_data);
            CSVObj* decompress(BinObj* compress_data) override;

        public:
            Decompression(std::string output, int interval, std::time_t basetime) : BaseDecompression(output, interval, basetime) {}
            void initialize(int count, char** params) override;
            void finalize() override;
    };
};

namespace Kangaroo {
    // Source paper: Kangaroo: Efficient Lossless Floating-Point Compression via Dynamic Reference Selection
    // Source path: src/floating-point/lossless/kangaroo.cpp

    inline float END_SIGNAL = INFINITY;

    class Kangaroo32Utils {
        private:
            static const double P10[45];
            static const int leading[8];

            static int __eraseBits(int exponent, int alpha);
            static int __decimalPlaces(uint32_t value, int exponent);
            static int __recoverPrecision(uint32_t value);

        public:
            static uint32_t erase(uint32_t value, int& tag);
            static uint32_t recover(uint32_t value, int tag);
            static int leadingCode(int zeros);
            static int leadingZeros(int code);
    };

    class Compression : public BaseCompression {
        private:
            long size = 0;
            int count = 0;
            int buffer_size = 0;
            int window_size = 0;
            BinObj* bitstream = nullptr;

            int index_bits = 0;
            uint64_t history_index = 0;
            std::vector<uint32_t> stored_values;
            uint64_t previous_indices[33];
            uint64_t lookup[33][32];
            uint32_t lookup_mask[33];
            uint32_t previous_values[33];
            uint32_t previous_value = 0;
            int history_slot = 0;
            bool power_of_two = false;
            uint32_t window_mask = 0;

            // Memoize deterministic erasure, including nonconsecutive repeats.
            uint32_t cache_values[256];
            uint32_t cache_erased[256];
            unsigned char cache_tags[256];
            bool cache_valid[256];
            bool cache_enabled = true;
            int cache_queries = 0;
            int cache_hits = 0;

            uint32_t run_value = 0;
            uint32_t run_count = 0;

            void __write(uint64_t value, int bits);
            void __writeRun(uint32_t length);
            void __flushRun();
            void __xorCompress(uint32_t value, int tag, int trailing);
            uint64_t __findReference(uint32_t value, int trailing);
            void __remember(uint32_t value, uint32_t repeats, int trailing);

        protected:
            void __compress(float value);
            void compress(Univariate* data) override;
            BinObj* serialize() override;

        public:
            Compression(std::string output) : BaseCompression(output) {}
            void initialize(int count, char** params) override;
            void finalize() override;
    };

    class Decompression : public BaseDecompression {
        private:
            int buffer_size = 0;
            int window_size = 0;
            int index_bits = 0;
            uint64_t history_index = 0;
            std::vector<uint32_t> stored_values;

            uint32_t repeat_value = 0;
            uint32_t repeat_count = 0;
            bool finished = false;
            int history_slot = 0;

            uint32_t __readRun(BinObj* compress_data);
            uint32_t __xorDecompress(BinObj* compress_data);

        protected:
            float __decompress(BinObj* compress_data);
            CSVObj* decompress(BinObj* compress_data) override;

        public:
            Decompression(std::string output, int interval, std::time_t basetime) : BaseDecompression(output, interval, basetime) {}
            void initialize(int count, char** params) override;
            void finalize() override;
    };
};
