#include "base-c.hpp"
#include "base-d.hpp"


namespace CachedNormalEquation {
    // Source paper: Fast Piecewise Polynomial Fitting of Time-Series Data for Streaming Computing
    // Source path: src/piecewise-approximation/non-linear/normal-equation.cpp
    
    struct Model {
        int length = 0;
        float error = -1;
        Polynomial* function = nullptr;

        Model(Polynomial* function);
        ~Model();
    };

    class Compression : public BaseCompression {
        private:
            int degree = -1;
            double error = 0;
            std::string mode = "";

            Model* model = nullptr;
            std::vector<Point2D> window;

            std::map<int, Eigen::MatrixXd> cache;       // Use eigen library
            
            Polynomial* __calPolynomial();
            bool __approxSuccess(Model* model);

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
            int degree = 0;

        protected:
            CSVObj* decompress(BinObj* compress_data) override;
        
        public:
            Decompression(std::string output, int interval, std::time_t basetime) : BaseDecompression(output, interval, basetime) {} 
            void initialize(int count, char** params) override;
            void finalize() override;
    };
};

namespace PolySwab {
    // Source paper: An Online Algorithm for Segmenting Time Series
    // Source path: src/piecewise-approximation/non-linear/poly-swab.cpp

    class Approximator {
        private:
            static std::vector<double> __interpolate(int degree, std::vector<double>& data);
            static std::vector<double> __regression(int degree, std::vector<double>& data);
        
        public:
            static Polynomial approximate(std::string mode, int degree, std::vector<double>& data);
            static double cal_error(std::vector<double>& segment, Polynomial& model);
    };

    class Grouper {
        public:
            static void merge(std::vector<double>& s1, std::vector<double>& s2, std::string mode);
            static double merge_cost(std::vector<double>& s1, std::vector<double>& s2, std::string mode, int degree, double error);
            static bool bound_check(std::vector<double>& segment, Polynomial& model, double error, bool skip=false);
    };

    class Compression : public BaseCompression {
        private:
            int degree = -1;
            int n_segment = 0;
            double error = 0;
            std::string mode = "";

            bool first = true;
            std::vector<double> window;
            std::vector<Polynomial> models; 
            std::vector<double> m_err;
            std::vector<std::vector<double>> segments;

            void __bottom_up();
            bool __sliding_window();

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
            int degree = -1;
            std::string mode = "";
            double pivot = INFINITY;
            
        protected:
            CSVObj* decompress(BinObj* compress_data) override;
        
        public:
            Decompression(std::string output, int interval, std::time_t basetime) : BaseDecompression(output, interval, basetime) {} 
            void initialize(int count, char** params) override;
            void finalize() override;
    };
};