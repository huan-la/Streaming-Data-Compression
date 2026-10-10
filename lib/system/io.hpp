#ifndef SYSTEM_IO_HPP
#define SYSTEM_IO_HPP

#include <string>
#include <fstream>
#include <vector>
#include <sstream>

const char HEX[16] = {'0', '1', '2', '3', '4', '5', '6', '7',
                      '8', '9', 'A', 'B', 'C', 'D', 'E', 'F'};

class IOObj {
    public:
        virtual IOObj* getNext() = 0;
        virtual void setNext(IOObj* obj) = 0;
        virtual std::string toStr() = 0;

        static void clear(IOObj* head) {
            IOObj* curr = head;
            IOObj* next = nullptr;

            while (curr != nullptr) {
                next = curr->getNext();
                delete curr;
                curr = next;
            }
        }
};

class BinObj : public IOObj {
    public:
        int wBitPos = 0;
        int rBitPos = 0;
        int rIndex = 0;
        std::vector<unsigned char> byte_vector;
        BinObj* next = nullptr;

    public:
        void setNext(IOObj* obj) override {
            this->next = (BinObj*) obj;
        }

        IOObj* getNext() override {
            return this->next;
        }

        std::string toStr() override {
            BinObj* curr = this;
            std::string str = "";
            
            while (curr != nullptr) {
                for (unsigned char byte : this->byte_vector) {
                    str += HEX[(byte & 0xf0) >> 4];
                    str += HEX[(byte & 0x0f) >> 0];
                }

                curr = (BinObj*) curr->getNext();
            }

            return str;
        }

        void put(unsigned char data) {
            this->byte_vector.push_back(data);
        }

        void put(std::vector<unsigned char> data) {
            for (unsigned char& byte : data) {
                this->put(byte);
            }
        }

        unsigned char getByte() {
            unsigned char data = this->byte_vector[0];
            this->byte_vector.erase(this->byte_vector.begin());
            return data;
        }

        std::vector<unsigned char> getBytes(int n, bool last) {
            if (last) {
                std::vector<unsigned char> bytes(this->byte_vector.end() - n, this->byte_vector.end());
                this->byte_vector.resize(this->byte_vector.size() - n);
                
                return bytes;
            }
            else {
                std::vector<unsigned char> bytes(this->byte_vector.begin(), this->byte_vector.begin() + n);
                this->byte_vector.erase(this->byte_vector.begin(), this->byte_vector.begin() + n);

                return bytes;
            }
        }

        void put(uint32_t value, int bitCount) {
            if (bitCount > 0) {
                if (bitCount <= 8 - this->wBitPos) {
                    if (this->wBitPos == 0) this->byte_vector.push_back(value << (8 - bitCount));
                    else this->byte_vector.back() |= value << (8 - bitCount - this->wBitPos);
                    
                    this->wBitPos = (bitCount + this->wBitPos) % 8;
                }
                else {
                    int i = bitCount - (8 - this->wBitPos);
                    if (this->wBitPos != 0) this->byte_vector.back() |= value >> i;
                    else this->byte_vector.push_back(value >> i);

                    value = value & (0xFFFFFFFF >> (32 - i));
                    while (i > 0) {
                        if (i > 7) {
                            this->wBitPos = 0;
                            this->byte_vector.push_back(value >> (i - 8));
                        }   
                        else {
                            this->wBitPos = i;
                            this->byte_vector.push_back(value << (8 - i));
                        } 
                        
                        i -= 8;
                        if (i > 0) value = value & ((1u << i) - 1);
                    }
                }
            }
        }

        uint32_t getBits(int bitCount) {
            uint32_t value = 0;

            if (bitCount > 0) {
                if (bitCount <= 8 - this->rBitPos) {
                    value = ((0xFF >> this->rBitPos) & this->byte_vector[this->rIndex]) >> (8 - this->rBitPos - bitCount);
                    this->rBitPos = (this->rBitPos + bitCount) % 8;
                    if (this->rBitPos == 0) this->rIndex++;
                }
                else {
                    int i = bitCount - (8 - this->rBitPos);
                    value = (0xFF >> this->rBitPos) & this->byte_vector[this->rIndex++];
                    
                    while (i > 0) {
                        if (i > 7) {
                            value = (value << 8) | this->byte_vector[this->rIndex++];
                            this->rBitPos = 0;
                        }
                        else {
                            this->rBitPos = i;
                            value = (value << i) | ((this->byte_vector[this->rIndex] & (0xFF << (8 - i))) >> (8 - i)); 
                        }

                        i -= 8;
                    }
                }
            }

            return value;
        }

        void flushBits() {
            // Flush remaining bit in byte
            if (this->getSize() != 0) {
                rBitPos = 0;
                this->byte_vector.clear();
            }
        }

        void put(short data) {
            this->byte_vector.push_back((data >> 8) & 0xff);
            this->byte_vector.push_back((data >> 0) & 0xff);
        }

        short getShort() {
            short data = 0;
            data += (short) this->byte_vector[0] << 8;
            data += (short) this->byte_vector[1] << 0;

            this->byte_vector.erase(this->byte_vector.begin(), this->byte_vector.begin() + 2);
            return data;
        }

        void put(int data) {
            this->byte_vector.push_back((data >> 24) & 0xff);
            this->byte_vector.push_back((data >> 16) & 0xff);
            this->byte_vector.push_back((data >> 8) & 0xff);
            this->byte_vector.push_back((data >> 0) & 0xff);
        }

        int getInt() {
            int data = 0;
            data += (int) this->byte_vector[0] << 24;
            data += (int) this->byte_vector[1] << 16;
            data += (int) this->byte_vector[2] << 8;
            data += (int) this->byte_vector[3] << 0;

            this->byte_vector.erase(this->byte_vector.begin(), this->byte_vector.begin() + 4);
            return data;
        }

        void put(long data) {
            this->byte_vector.push_back((data >> 56) & 0xff);
            this->byte_vector.push_back((data >> 48) & 0xff);
            this->byte_vector.push_back((data >> 40) & 0xff);
            this->byte_vector.push_back((data >> 32) & 0xff);
            this->byte_vector.push_back((data >> 24) & 0xff);
            this->byte_vector.push_back((data >> 16) & 0xff);
            this->byte_vector.push_back((data >> 8) & 0xff);
            this->byte_vector.push_back((data >> 0) & 0xff);
        }

        long getLong() {
            long data = 0;
            data += (long) this->byte_vector[0] << 56;
            data += (long) this->byte_vector[1] << 48;
            data += (long) this->byte_vector[2] << 40;
            data += (long) this->byte_vector[3] << 32;
            data += (long) this->byte_vector[4] << 24;
            data += (long) this->byte_vector[5] << 16;
            data += (long) this->byte_vector[6] << 8;
            data += (long) this->byte_vector[7] << 0;

            this->byte_vector.erase(this->byte_vector.begin(), this->byte_vector.begin() + 8);
            return data;
        }

        void put(float data) {
            unsigned char* ptr = (unsigned char*) &data;
            for (int i=0; i< sizeof(float); i++) {
                this->byte_vector.push_back(ptr[i]);
            }
        }

        float getFloat() {
            float data = 0;
            *((unsigned char*)(&data) + 3) = this->byte_vector[3];
            *((unsigned char*)(&data) + 2) = this->byte_vector[2];
            *((unsigned char*)(&data) + 1) = this->byte_vector[1];
            *((unsigned char*)(&data) + 0) = this->byte_vector[0];

            this->byte_vector.erase(this->byte_vector.begin(), this->byte_vector.begin() + 4);
            return data;
        }

        void put(double data) {
            unsigned char* ptr = (unsigned char*) &data;
            for (int i = 0; i < sizeof(double); i++) {
                this->byte_vector.push_back(ptr[i]);
            }
        }

        double getDouble() {
            double data = 0;
            *((unsigned char*)(&data) + 7) = this->byte_vector[7];
            *((unsigned char*)(&data) + 6) = this->byte_vector[6];
            *((unsigned char*)(&data) + 5) = this->byte_vector[5];
            *((unsigned char*)(&data) + 4) = this->byte_vector[4];
            *((unsigned char*)(&data) + 3) = this->byte_vector[3];
            *((unsigned char*)(&data) + 2) = this->byte_vector[2];
            *((unsigned char*)(&data) + 1) = this->byte_vector[1];
            *((unsigned char*)(&data) + 0) = this->byte_vector[0];

            this->byte_vector.erase(this->byte_vector.begin(), this->byte_vector.begin() + 8);
            return data;
        }

        void put(std::string data) {
            for (unsigned char c : data) {
                this->byte_vector.push_back(c);
            }
        }

        std::string getString(int length) {
            std::string data = "";
            for (int i = 0; i < length; i++) {
                data += (char) this->byte_vector[i];
            }

            this->byte_vector.erase(this->byte_vector.begin(), this->byte_vector.begin() + length);
            return data;
        }

        int getSize() {
            return this->byte_vector.size();
        }

        unsigned char* serialize() {
            return &this->byte_vector[0];
        }
};

class CSVObj : public IOObj {
    private:
        std::vector<std::string> data;
        CSVObj* next = nullptr;

    public:
        void setNext(IOObj* obj) override {
            this->next = (CSVObj*) obj;
        }

        IOObj* getNext() override {
            return this->next;
        }

        std::string toStr() override {
            std::string str = "";
            for (std::string ele : this->data) {
                str += "," + ele;
            }

            return str.substr(1);
        }

        std::vector<std::string> getData() {
            return this->data;
        }

        std::string getData(int i) {
            return this->data.at(i);
        }

        void pushData(std::string data) {
            this->data.push_back(data);
        }
};

class VariableByteEncoding {
    public:
        static std::vector<unsigned char> encode(unsigned long number) {
            std::vector<unsigned char> bytes;
            do {
                unsigned char byte = number % 128; // get last 7 bits
                bytes.insert(bytes.begin(), byte);
                number /= 128;
            } while (number > 0);

            // Set the continuation bit (MSB = 1) on the last byte
            bytes.back() |= 0x80;

            return bytes;
        }

        static unsigned long decode(BinObj* obj) {
            long number = 0;

            while (true) {
                unsigned char byte = obj->getByte();
                if (byte & 0x80) { // this is the last byte (MSB = 1)
                    number = number * 128 + (byte & 0x7F);
                    break;
                }
                number = number * 128 + byte;
            }

            return number;
        }
};

class EliasGammaEncoding {
    private:
        constexpr static double kLog2Table[17] = {
            INFINITY, 0, 1, 1.584962500721156, 2, 2.321928094887362, 2.584962500721156, 2.807354922057604, 3, 3.169925001442312, 
            3.321928094887362, 3.459431618637297, 3.584962500721156, 3.700439718141092, 3.807354922057604, 3.906890595608519, 4
        };

    public:
        static int encode(unsigned long number, BinObj* bitstream) {
            if (number == 1) {
                bitstream->put(1, 1);
                return 1;
            }
            else {
                int n = std::floor(log2(number));
                bitstream->put(0, n);
                bitstream->put(number, n + 1);

                return n + n + 1;
            }
        }

        static unsigned long decode(BinObj* bitstream) {
            int n = 0;
            while (!bitstream->getBits(1)) n++;
            return n == 0 ? 1 :  (1 << n) | bitstream->getBits(n);
        }
};

class ZigZagEncoding {
    public:
        static unsigned long encode(long number) {
            return (number << 1) ^ (number >> 63);
        }

        static long decode(unsigned long number) {
            return (number >> 1) ^ -(number & 1);
        }
};

class DigitEncoding {
    // Based on Huffman encoding only 0->9 digits
    public:
        static int encode(int digit, BinObj* bitstream) {
            switch (digit) {
                case 0: { bitstream->put(0b000, 3); return 3; }
                case 1: { bitstream->put(0b001, 3); return 3; }
                case 2: { bitstream->put(0b010, 3); return 3; }
                case 3: { bitstream->put(0b011, 3); return 3; }
                case 4: { bitstream->put(0b100, 3); return 3; }
                case 5: { bitstream->put(0b101, 3); return 3; }
                case 6: { bitstream->put(0b1100, 4); return 4; }
                case 7: { bitstream->put(0b1101, 4); return 4; }
                case 8: { bitstream->put(0b1110, 4); return 4; }
                case 9: { bitstream->put(0b1111, 4); return 4; }
                default: return 0;
            }
        }

        static unsigned char decode(BinObj* bitstream) {
            if (bitstream->getBits(1) == 0) return bitstream->getBits(2); // 0xx
            else if (bitstream->getBits(1) == 0) return 4 + bitstream->getBits(1); // 10x
            else return 6 + bitstream->getBits(2); // 11x
        }
};

class BatchIO {
    public:
        static long write(std::string filename, IOObj* obj, bool append = false) {
            long count = 0;
            std::ofstream file;

            if (!append) file.open(filename);
            else file.open(filename, std::ios::app);

            while (obj != nullptr) {
                file << obj->toStr();
                obj = obj->getNext();
                if (obj != nullptr) {
                    file << "\n";
                }
                count++;
            }

            file.close();

            return count;
        }

        static CSVObj* readCSV(std::string filename) {
            std::ifstream file(filename);

            std::string line;
            CSVObj* head = nullptr;
            CSVObj* current = nullptr;

            while(std::getline(file, line)) {
                CSVObj* data = new CSVObj;
                std::stringstream ss(line);
                std::string ele;

                while(std::getline(ss, ele, ',')) {
                    data->pushData(ele);
                }

                if (head == nullptr) {
                    head = data;
                    current = head;
                }
                else {
                    current->setNext(data);
                    current = data;
                }
            }

            file.close();
            return head;
        }
};

class IterIO {
    private:
        std::fstream file;
        bool mode;

    public:
        IterIO(std::string filename, bool read, bool binary = false, bool append = false) {
            auto flag = read ? std::ios::in : std::ios::out;
            if (binary) flag = flag | std::ios::binary;
            if (append) flag = flag | std::ios::app;
            
            this->file.open(filename, flag);
        }

        void close() {
            if (this->file.is_open()) {
                this->file.close();
            }
        }

        long write(std::string data, bool endline = true) {
            this->file << data;
            if (endline) {
                this->file << "\n";
            }

            return 1;
        }

        long write(IOObj* obj, bool endline = true) {
            long count = 0;

            while (obj != nullptr) {
                this->file << obj->toStr();
                if (endline) {
                    this->file << "\n";
                }

                obj = obj->getNext();
                count++;
            }

            return count;
        }

        long writeBin(BinObj* obj) {
            long count = 0;
            
            while (obj != nullptr) {
                this->file.write((char*) obj->serialize(), obj->getSize());
                obj = (BinObj*) obj->getNext();
                count++;
            }

            return count;
        }

        CSVObj* readCSV() {
            std::string ele;
            std::string line;
            CSVObj* data = nullptr;
            std::stringstream ss(line);

            if (std::getline(this->file, line)) {
                data = new CSVObj;
                while(std::getline(ss, ele, ',')) {
                    data->pushData(ele);
                }
            }

            return data;
        }

        BinObj* readBin() {
            unsigned char* byte = new unsigned char;
            BinObj* obj = new BinObj;

            while (true) {
                this->file.read((char*) byte, sizeof(unsigned char));
                if (this->file.eof()) break;
                else obj->put(*byte);
            }

            return obj;
        }
};

#endif
