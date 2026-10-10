#include "base-c.hpp"
#include "base-d.hpp"

namespace AdaptiveApproximation {
    // Source paper: An Adaptive Algorithm for Online Time Series Segmentation with Error Bound Guarantee
    // Source path: src/model-selection/adaptive-approximation.cpp
    class Model {
        protected:
            int length = 0;
            std::string type = "";
            Point2D* pivot = nullptr;
            double* coeffs = nullptr;

        public:
            bool complete = false;

            ~Model();
            int getDim();
            int getLen();
            double* getCoeffs();
            std::string getType();

            virtual double subs(double x) = 0;
            virtual bool approximate(double bound, double data) = 0;
    };

    class ExpoFunction : public Model {
        private:
            double upper = INFINITY;
            double lower = -INFINITY;

        public:
            ExpoFunction();
            double subs(double x) override;
            bool approximate(double bound, double data) override;
    };

    class LinearFunction : public Model {
        private:
            double upper = INFINITY;
            double lower = -INFINITY;

        public:
            LinearFunction();
            double subs(double x) override;
            bool approximate(double bound, double data) override;
    };

    class PolyFunction : public Model {
        private:
            int degree;
            SDLP sdlp;

        public:
            PolyFunction(int degree, std::string type);
            double subs(double x) override;
            bool approximate(double bound, double data) override;
    };

    class Compression : public BaseCompression {
        private:
            double penalty = 0;
            double error = 0;
            std::vector<std::string> types;

            bool first = true;
            bool Fcf = false;
            bool Frv = false;
            int best_model = -1;
            std::vector<std::vector<Model*>> candidates;

            void __init_candidates();
            void __choose_best_model();
            
        protected:
            void compress(Univariate* data) override;
            BinObj* serialize() override;

        public:
            Compression(std::string output) : BaseCompression(output) {}      
            void initialize(int count, char** params) override;
            void finalize() override;
    };

    class Decompression : public BaseDecompression {
        protected:
            double pivot = INFINITY;
            CSVObj* decompress(BinObj* compress_data) override;

        public:
            Decompression(std::string output, int interval, std::time_t basetime) : BaseDecompression(output, interval, basetime) {}
            void initialize(int count, char** params) override;
            void finalize() override;
    };
};

namespace SmartGridCompression {
    // Source paper: A time-series compression technique and its application to the smart grid
    // Source path: src/model-selection/smart-grid-compression.cpp

    class Model {
        public:
            int degree;
            int length;

            virtual void clear() = 0;
            virtual double getCompressionRatio() = 0;
            virtual bool approximate(double bound, std::vector<Point2D>& segment) = 0;
    };

    // Constant approximate
    class ConstantModel : public Model {
        private:
            double value;
            double min;
            double max;

        public:
            ConstantModel();
            double getValue();

            void clear() override;
            double getCompressionRatio() override;
            bool approximate(double bound, std::vector<Point2D>& segment) override;
    };

    // Linear approximate    
    class LinearModel : public Model {
        private:
            Line* line;
            Point2D A; 
            Point2D B; 
            Point2D C;
            ConvexHull cvx;

            bool __verify(double bound, Line& line);
            double __distance(Line line, Point2D p);
            int __x_external(int x, int x1, int x2);
            int __search(int side_index, int prev_v_l);
            Line __approx(const std::vector<Point2D>& segment);

        public:
            LinearModel();
            Line* getLine();

            void clear() override;
            double getCompressionRatio() override;
            bool approximate(double bound, std::vector<Point2D>& segment) override;
    };

    // Polynomial approximate
    class PolynomialModel : public Model {
        private:
            SDLP sdlp;
            Polynomial* polynomial;

        public:
            PolynomialModel(int degree);
            Polynomial* getPolynomial();

            void clear() override;
            double getCompressionRatio() override;
            bool approximate(double bound, std::vector<Point2D>& segment) override;
    };

    class Compression : public BaseCompression {
        private:
            int max_degree = 0;
            double error = 0;

            int curr_degree = 0;
            int chosen_model = -1;
            std::vector<Point2D> window;
            std::vector<Model*> models;

            bool __approximate(int degree);
            void __choose_best_model();

        protected:
            void compress(Univariate* data) override;
            BinObj* serialize() override;

        public:
            Compression(std::string output) : BaseCompression(output) {}      
            void initialize(int count, char** params) override;
            void finalize() override;
    };

    class Decompression : public BaseDecompression {
        protected:
            CSVObj* decompress(BinObj* compress_data) override;

        public:
            Decompression(std::string output, int interval, std::time_t basetime) : BaseDecompression(output, interval, basetime) {}
            void initialize(int count, char** params) override;
            void finalize() override;
    };
};

namespace AdaptPPA {
    // Source paper: Shape-Aware Compression of Streaming Data: An Online Dynamic Polynomial Selection Approach
    // Source path: src/model-selection/adapt-ppa.cpp
    class LinearSegment {
        private:
            UpperHull u_cvx;
            LowerHull l_cvx;
            Point2D* pivot = nullptr;
            
        public:
            bool is_complete = false;
            Line* u_line = nullptr;
            Line* l_line = nullptr;
            Point2D* p_start = nullptr;
            Point2D* p_end = nullptr;
            
            ~LinearSegment();

            Line getLine();
            int getLength();
            void translation(Point2D* p);
            void approximate(Point2D& p, double error);
    };

    class Compression : public BaseCompression {
        private:
            int max_degree = 0;
            double error = 0;

            int degree = 1;
            int direction = 0;
            int slope = 0;

            Point2D* p1 = nullptr;
            Point2D* p2 = nullptr;
            Point2D* p3 = nullptr;

            LinearSegment* seg_1 = nullptr;
            LinearSegment* seg_2 = nullptr;

            SDLP sdlp;
            std::vector<Point2D> window;

            bool __merge_check();
            double* __approximate(int pivot);

        protected:
            void compress(Univariate* data) override;
            BinObj* serialize() override;

        public:
            Compression(std::string output) : BaseCompression(output) {}      
            void initialize(int count, char** params) override;
            void finalize() override;
    };

    class Decompression : public BaseDecompression {
        protected:
            CSVObj* decompress(BinObj* compress_data) override;

        public:
            Decompression(std::string output, int interval, std::time_t basetime) : BaseDecompression(output, interval, basetime) {}
            void initialize(int count, char** params) override;
            void finalize() override;
    };
};