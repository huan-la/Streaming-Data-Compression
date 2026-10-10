#ifndef BASE_COMPRESSION_HPP
#define BASE_COMPRESSION_HPP
#define EPS 0.0000001

#include <map>
#include <bit>
#include <cmath>
#include <vector>
#include <queue>
#include <deque>
#include <utility>
#include <cstdint>
#include <climits>
#include <iostream>
#include <algorithm>
#include <eigen/Eigen>

#include "system/io.hpp"
#include "system/monitor.hpp"
#include "algebraic/function.hpp"
#include "algebraic/convex.hpp"
#include "algebraic/sdlp.hpp"
#include "timeseries.hpp"

using namespace std::chrono;

class BaseCompression {
    private:
        bool y_trigger = false;
        BinObj* headObj = nullptr;
        BinObj* currObj = nullptr;
        IterIO* comFile = nullptr;

    protected:
        virtual void compress(Univariate* data) = 0;
        virtual BinObj* serialize() = 0;

        void yield() {
            this->y_trigger = true;

            if (this->headObj == nullptr) {
                this->headObj = this->serialize();
                this->currObj = this->headObj;
            }
            else {
                this->currObj->setNext(this->serialize());
                this->currObj = (BinObj*) this->currObj->getNext();
            }
        }

    public:
        virtual void initialize(int count, char** params) = 0;
        virtual void finalize() = 0;

        BaseCompression(std::string output) {
            this->comFile = new IterIO(output, false);
        }

        ~BaseCompression() {
            // Write compression data
            this->comFile->close();
            IOObj::clear(this->headObj);
        }

        void clear_buffer() {
            IOObj::clear(this->headObj);
            this->headObj = nullptr;
            this->currObj = nullptr;
        }

        BinObj* process(Univariate* data) {
            this->compress(data);
            
            if (this->y_trigger) {
                this->y_trigger = false;
                this->comFile->writeBin(this->headObj);

                return this->headObj;
            }

            return nullptr;
        }

        BinObj* complete() {
            this->finalize();

            if (this->y_trigger) {
                this->y_trigger = false;
                this->comFile->writeBin(this->headObj);

                return this->headObj;
            }

            return nullptr;
        }
};

#endif