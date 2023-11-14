// MyClass.h
#ifndef MYCLASS_H
#define MYCLASS_H

#include <string>


class MyClass {
public:
    MyClass(const std::string& out_filename, const std::string& in_filename) : out_filename(out_filename), in_filename(in_filename) {}
    void process_data();

private:
    std::string out_filename;
    std::string in_filename;
};

#endif // MYCLASS_H
